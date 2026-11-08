"""Input data: index of the raw files, choice of the vertical component, station names, reading one processing day.

The user gives only data.raw_dir and data.station_file.  Everything else is read from the files:
  * files: any layout, any names, any format ObsPy reads (format detected per file); files that are not seismic
    data are listed and skipped;
  * station, network, location, channel, sampling rate, times: from the record headers;
  * vertical component (data.channel = auto): a channel code ending in Z, or a SAC header with cmpinc = 0; a station
    with a single channel and nothing that says otherwise uses that channel (reported); one location / channel /
    rate per station (the one with most data);
  * names: the station table is read in whatever form it has and matched to the header codes (lant/stations.py);
    from then on every station is NET.STA as written in 01_clean/stations.csv.
One processing day runs from data.day_start_utc_hour UTC to the same hour on the next day and is labelled YYYYDDD
by the UTC date of its start.
"""
import csv
import fnmatch
import json
import os
from collections import defaultdict

import numpy as np

from . import stations as ST


class DataError(ValueError):
    pass


# --------------------------------------------------------------------------------------------- station tables
def read_stations(path):
    """Any station table -> {NET.STA: dict(lon, lat, elev)} without matching to data (names as in the table)."""
    if not path:
        raise DataError("data.station_file is not set")
    if not os.path.exists(path):
        raise DataError("station file %s not found" % path)
    try:
        rows, _ = ST.read_table(path)
    except ST.StationError as exc:
        raise DataError(str(exc))
    out = {}
    for r in rows:
        net, sta = r["network"], r["station"]
        if net is None:
            net, sta = ST.split_name(sta)
        out["%s.%s" % (net or "XX", sta)] = dict(lon=r["lon"], lat=r["lat"], elev=r["elev"])
    return out


def canonical_table(out_dir):
    return os.path.join(out_dir, "01_clean", "stations.csv")


def write_table(path, coords):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["network", "station", "latitude", "longitude", "elevation"])
        for name, c in sorted(coords.items()):
            n, s = name.split(".", 1)
            w.writerow([n, s, "%.7f" % c["lat"], "%.7f" % c["lon"], "" if not np.isfinite(c["elev"]) else "%.1f" % c["elev"]])


def read_canonical(path):
    out = {}
    for r in csv.DictReader(open(path, newline="")):
        out["%s.%s" % (r["network"], r["station"])] = dict(lon=float(r["longitude"]), lat=float(r["latitude"]),
                                                           elev=float(r["elevation"]) if r["elevation"] else np.nan)
    return out


def resolve_names(names, known, what):
    """User-given station names (any of 'YY.1001', 'YY1001', '1001') -> canonical names."""
    out = []
    for n in names or []:
        n = str(n)
        if n in known:
            out.append(n)
            continue
        net, sta = ST.split_name(n)
        cand = [k for k in known if k.split(".", 1)[1] == (sta if net is None else sta) or k.replace(".", "") == n]
        if net:
            cand = [k for k in cand if k.split(".")[0] == net] or cand
        if len(cand) != 1:
            raise DataError("%s: station %r %s" % (what, n, "matches %s" % cand if cand else "is not in the station table"))
        out.append(cand[0])
    return out


def used_stations(cfg):
    """The stations used: the matched station table of the clean step minus data.exclude_stations."""
    p = canonical_table(cfg["project"]["output_dir"])
    if not os.path.exists(p):
        raise DataError("%s not found: run `lant run clean` first" % p)
    sta = read_canonical(p)
    ex = set(resolve_names(cfg["data"]["exclude_stations"], sta, "data.exclude_stations"))
    return {k: v for k, v in sta.items() if k not in ex}


# --------------------------------------------------------------------------------------------- raw file index
def find_files(raw_dir, patterns, ignore):
    if not raw_dir or not os.path.isdir(raw_dir):
        raise DataError("data.raw_dir %s is not a directory" % raw_dir)
    out = []
    for root, dirs, files in os.walk(raw_dir, followlinks=True):
        dirs[:] = sorted(d for d in dirs if not d.startswith("."))
        for f in sorted(files):
            if any(fnmatch.fnmatch(f, p) for p in patterns) and not any(fnmatch.fnmatch(f, p) for p in ignore):
                out.append(os.path.join(root, f))
    return out


def _scan_one(args):
    path, fmt = args
    import obspy
    try:
        st = obspy.read(path, format=None if fmt == "auto" else fmt, headonly=True)
    except Exception as exc:          # not a seismic data file, or corrupt
        return path, None, "%s: %s" % (type(exc).__name__, str(exc)[:150])
    rows = []
    for tr in st:
        s = tr.stats
        coord, cmpinc = None, None
        if "sac" in s:
            sac = s.sac
            la, lo = sac.get("stla"), sac.get("stlo")
            if la is not None and lo is not None and abs(la) <= 90 and la != -12345.0 and lo != -12345.0:
                el = sac.get("stel", -12345.0)
                coord = (float(lo), float(la), float(el) if el != -12345.0 else np.nan)
            ci = sac.get("cmpinc")
            if ci is not None and ci != -12345.0:
                cmpinc = float(ci)
        rows.append(dict(file=path, net=s.network, sta=s.station, channel=s.channel, location=s.location,
                         fs=float(s.sampling_rate), start=s.starttime.timestamp,
                         end=s.endtime.timestamp + 1.0 / s.sampling_rate, npts=int(s.npts),
                         format=str(s.get("_format", "")), coord=coord, cmpinc=cmpinc))
    return path, rows, None


def build_index(cfg, nproc=4):
    """Scan the headers of every candidate file -> index dict (all channels)."""
    from multiprocessing import Pool
    d = cfg["data"]
    files = find_files(d["raw_dir"], d["file_patterns"], d["ignore_patterns"])
    if not files:
        raise DataError("no candidate files under %s" % d["raw_dir"])
    with Pool(nproc) as pool:
        res = pool.map(_scan_one, [(f, d["format"]) for f in files], chunksize=16)
    traces, unreadable = [], []
    for path, rows, err in res:
        if err:
            unreadable.append((path, err))
        else:
            traces.extend(rows)
    return dict(raw_dir=d["raw_dir"], n_files=len(files), unreadable=unreadable, traces=traces)


def is_vertical(tr):
    if tr.get("cmpinc") is not None:
        return abs(tr["cmpinc"]) < 1.0
    return tr["channel"].upper().endswith("Z")


def choose_channel(traces, pattern):
    """Per header station: the records of its vertical component.  Returns (rows by (net, sta), notes)."""
    by = defaultdict(list)
    for t in traces:
        by[(t["net"], t["sta"])].append(t)
    out, notes = {}, []
    for code, rows in by.items():
        if pattern != "auto":
            sel = [r for r in rows if fnmatch.fnmatch(r["channel"], pattern)]
        else:
            sel = [r for r in rows if is_vertical(r)]
            chans = sorted({r["channel"] for r in rows})
            if not sel and len(chans) == 1:
                sel = rows
                notes.append("%s.%s: only channel %s, not recognisable as vertical: used" % (code[0], code[1], chans[0]))
            elif not sel:
                notes.append("%s.%s: channels %s, none recognisable as vertical: station not used (set data.channel)"
                             % (code[0], code[1], ",".join(chans)))
        if sel:
            out[code] = sel
    return out, notes


def best_sensor(rows):
    """One location / channel / sampling rate per station: the one with most data."""
    dur = defaultdict(float)
    for r in rows:
        dur[(r["location"], r["channel"], r["fs"])] += r["end"] - r["start"]
    best = max(dur, key=dur.get)
    other = {"%s.%s %g Hz" % k: round(v / 3600, 2) for k, v in dur.items() if k != best}
    return [r for r in rows if (r["location"], r["channel"], r["fs"]) == best], ("%s.%s %g Hz" % best, other)


def day_label(t0):
    import obspy
    u = obspy.UTCDateTime(t0)
    return "%04d%03d" % (u.year, u.julday)


def day_windows(traces, cfg):
    """Processing days that hold data: {label: (t0, [trace rows])}."""
    import obspy
    h = cfg["data"]["day_start_utc_hour"]
    lo = obspy.UTCDateTime(cfg["data"]["start_date"]).timestamp if cfg["data"]["start_date"] else -np.inf
    hi = obspy.UTCDateTime(cfg["data"]["end_date"]).timestamp + 86400 if cfg["data"]["end_date"] else np.inf
    days = defaultdict(list)
    for tr in traces:
        u = obspy.UTCDateTime(tr["start"])
        t0 = obspy.UTCDateTime(u.year, u.month, u.day, h).timestamp
        if t0 > tr["start"]:
            t0 -= 86400
        while t0 < tr["end"]:
            if t0 + 86400 > lo and t0 < hi:
                days[t0].append(tr)
            t0 += 86400
    return {day_label(t0): (t0, rows) for t0, rows in sorted(days.items())}


def processing_rows(cfg):
    """Trace rows of the cleaned archive (or of the raw files when clean.write_archive is false)."""
    p = os.path.join(cfg["project"]["output_dir"], "01_clean", "archive_index.json")
    if not os.path.exists(p):
        raise DataError("%s not found: run `lant run clean` first" % p)
    return json.load(open(p))["traces"]


def station_rows(traces, stations):
    """{station: [rows]} for the used stations (rows carry the canonical name in 'station')."""
    out = defaultdict(list)
    for t in traces:
        if t["station"] in stations:
            out[t["station"]].append(t)
    return dict(out), {}


def read_day(rows, t0, fs):
    """Traces of one station for the day starting at t0 -> (data at fs, valid mask) or (None, None).

    Samples placed on the UTC grid of their rate (whole-sample shift; Lanczos interpolation only when the offset
    exceeds 0.05 sample; the cleaned archive is already on the grid), median removed, zero-phase 8th-order
    Butterworth low-pass at 0.8 x the new Nyquist frequency, decimation.  A decimated sample is valid only if the raw
    data within +-1 s of it were all present."""
    import obspy
    from scipy.signal import butter, sosfiltfilt
    raw_fs = rows[0]["fs"]
    n_raw = int(round(86400 * raw_fs))
    x = np.zeros(n_raw, np.float64)
    ok = np.zeros(n_raw, bool)
    want = (rows[0]["net"], rows[0]["sta"], rows[0]["location"], rows[0]["channel"], raw_fs)
    for path in sorted({r["file"] for r in rows}):
        st = obspy.read(path, starttime=obspy.UTCDateTime(t0 - 1), endtime=obspy.UTCDateTime(t0 + 86401))
        for tr in st:
            s = tr.stats
            if (s.network, s.station, s.location, s.channel, s.sampling_rate) != want:
                continue
            off = (s.starttime.timestamp - t0) * raw_fs
            k0 = int(round(off))
            if abs(off - k0) > 0.05:
                tr.interpolate(raw_fs, starttime=obspy.UTCDateTime(t0 + k0 / raw_fs), method="lanczos", a=20)
            d = np.asarray(tr.data, dtype=np.float64)
            a, b = max(k0, 0), min(k0 + d.size, n_raw)
            if b > a:
                x[a:b] = d[a - k0:b - k0]
                ok[a:b] = True
    if not ok.any():
        return None, None
    x[ok] -= np.median(x[ok])
    x[~ok] = 0.0
    q = int(round(raw_fs / fs))
    sos = butter(8, 0.8 * fs / 2, btype="low", fs=raw_fs, output="sos")
    y = sosfiltfilt(sos, x)[::q].astype(np.float32)
    bad = ~ok
    if bad.any():
        kern = int(raw_fs)
        cs = np.concatenate([[0], np.cumsum(bad)])
        idx = np.arange(n_raw)
        lo = np.clip(idx - kern, 0, n_raw)
        hi = np.clip(idx + kern + 1, 0, n_raw)
        bad = (cs[hi] - cs[lo]) > 0
    return y, ~bad[::q]
