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
One processing day runs from data.day_start_utc_hour UTC to the same hour on the next day (a fixed 24-h UTC window,
kept as its UTC start t0).  Its label YYYYDDD (file names, start_date / end_date, stack.exclude_station_days) is the
date of t0 by data.day_label: "utc_start" (the UTC date of the start; the default) or a fixed offset "UTC+8",
"UTC-05:30" (the date of the start at that offset, e.g. the local date when the day starts at local midnight).
Daylight-saving time is not followed: the windows stay fixed 24-h UTC windows.  01_clean/processing_days.txt lists
every label with its UTC window.
"""
import csv
import fnmatch
import json
import os
import re
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


def read_quiet(path, **kw):
    """obspy.read without the warnings of the reader -> (stream, number of warnings).  Recorder files often carry
    blocks that are not data (padding of pre-allocated files): the miniSEED reader skips them with one warning
    each, hundreds of thousands per file; they are counted, not printed."""
    import warnings
    import obspy
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        st = obspy.read(path, **kw)
    return st, len(w)


def _scan_one(args):
    path, fmt = args
    try:
        st, nwarn = read_quiet(path, format=None if fmt == "auto" else fmt, headonly=True)
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
    return path, rows, nwarn


def build_index(cfg, nproc=4):
    """Scan the headers of every candidate file -> index dict (all channels)."""
    from multiprocessing import Pool
    d = cfg["data"]
    files = find_files(d["raw_dir"], d["file_patterns"], d["ignore_patterns"])
    if not files:
        raise DataError("no candidate files under %s" % d["raw_dir"])
    # headers read by an earlier run are used again for the files that have not changed (path, size, time)
    cache_path = os.path.join(cfg["project"]["output_dir"], "01_clean", "header_cache.json")
    cache, sig = {}, {}
    try:
        c = json.load(open(cache_path))
        if c.get("format") == d["format"]:
            cache = c["files"]
    except (OSError, ValueError, KeyError):
        pass
    for f in files:
        st = os.stat(f)
        sig[f] = [st.st_size, int(st.st_mtime)]
    todo = [f for f in files if f not in cache or cache[f][0] != sig[f]]
    with Pool(nproc) as pool:
        new = {r[0]: r for r in pool.map(_scan_one, [(f, d["format"]) for f in todo], chunksize=16)}
    res = [new[f] if f in new else (f, cache[f][1], cache[f][2]) for f in files]
    try:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        with open(cache_path + ".tmp", "w") as fh:
            json.dump(dict(format=d["format"], files={r[0]: [sig[r[0]], r[1], r[2]] for r in res}), fh)
        os.replace(cache_path + ".tmp", cache_path)
    except OSError:
        pass
    traces, unreadable, skipped = [], [], []
    for path, rows, err in res:
        if rows is None:
            unreadable.append((path, err))
        else:
            traces.extend(rows)
            if err:
                skipped.append((path, int(err)))
    return dict(raw_dir=d["raw_dir"], n_files=len(files), unreadable=unreadable, traces=traces, reader_warnings=skipped)


def orientation(tr):
    """'vertical', 'horizontal' or 'unknown' from the SAC inclination (cmpinc: 0 or 180 vertical, 90 horizontal) or
    the last letter of the SEED channel code (Z vertical; N, E, R, T, 1, 2 horizontal)."""
    ci = tr.get("cmpinc")
    if ci is not None:
        if abs(ci) < 1.0 or abs(abs(ci) - 180.0) < 1.0:
            return "vertical"
        if abs(abs(ci) - 90.0) < 1.0:
            return "horizontal"
    last = tr["channel"].upper()[-1:] if tr["channel"] else ""
    if last == "Z":
        return "vertical"
    if last in ("N", "E", "R", "T", "1", "2"):
        return "horizontal"
    return "unknown"


def is_down(tr):
    """Vertical sensor pointing down (SAC cmpinc = 180): its polarity is reversed."""
    ci = tr.get("cmpinc")
    return ci is not None and abs(abs(ci) - 180.0) < 1.0


def is_vertical(tr):
    return orientation(tr) == "vertical"


def choose_channel(traces, pattern):
    """Per header station: the records of its vertical component.  Returns (rows by (net, sta), notes).

    data.channel = auto: the channels recognised as vertical (lant.data.orientation).  A station whose channels are
    all recognised as horizontal is not used.  A station with a single channel of unknown orientation uses it and
    is listed as "orientation unknown: confirm" (docs/ai_workflow.md: the user must confirm).  A pattern selects
    the channels by name instead."""
    by = defaultdict(list)
    for t in traces:
        by[(t["net"], t["sta"])].append(t)
    out, notes = {}, []
    for code, rows in by.items():
        chans = sorted({r["channel"] for r in rows})
        if pattern != "auto":
            sel = [r for r in rows if fnmatch.fnmatch(r["channel"], pattern)]
        else:
            sel = [r for r in rows if orientation(r) == "vertical"]
            unknown = sorted({r["channel"] for r in rows if orientation(r) == "unknown"})
            if not sel and len(chans) == 1 and unknown:
                sel = rows
                notes.append("%s.%s: only channel %s, orientation unknown: used - confirm that it is vertical (or set "
                             "data.channel)" % (code[0], code[1], chans[0]))
            elif not sel and not unknown:
                notes.append("%s.%s: channels %s are horizontal: station not used (only vertical-component "
                             "correlations are supported)" % (code[0], code[1], ",".join(chans)))
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


def label_offset_s(spec):
    """data.day_label -> offset (s) added to the UTC start before taking its date."""
    if spec in (None, "utc_start"):
        return 0.0
    if isinstance(spec, (int, float)) and not isinstance(spec, bool):
        h = float(spec)
    else:
        m = re.fullmatch(r"\s*UTC\s*([+-])\s*(\d{1,2})(?::(\d{2}))?\s*", str(spec), re.I)
        if not m:
            raise ValueError("data.day_label must be 'utc_start', a fixed offset such as 'UTC+8' or 'UTC-05:30', or "
                             "a number of hours; got %r" % (spec,))
        h = (1 if m.group(1) == "+" else -1) * (int(m.group(2)) + int(m.group(3) or 0) / 60.0)
    if not -14 <= h <= 14:
        raise ValueError("data.day_label: offset %g h outside -14..+14 h" % h)
    return h * 3600.0


def day_label(t0, cfg=None, offset_s=None):
    """Label YYYYDDD of the processing day that starts at t0 (UTC timestamp) under data.day_label."""
    import obspy
    off = offset_s if offset_s is not None else label_offset_s((cfg or {}).get("data", {}).get("day_label"))
    u = obspy.UTCDateTime(t0 + off)
    return "%04d%03d" % (u.year, u.julday)


def label_date(t0, cfg):
    """YYYY-MM-DD of the label of the processing day starting at t0 (what start_date / end_date compare with)."""
    import obspy
    return (obspy.UTCDateTime(t0 + label_offset_s(cfg["data"].get("day_label")))).strftime("%Y-%m-%d")


def utc_date(t):
    """YYYY-MM-DD of a UTC timestamp (record spans and the like; not a processing-day label)."""
    import obspy
    return obspy.UTCDateTime(t).strftime("%Y-%m-%d")


def utc_window(t0):
    import obspy
    return "%s - %s UTC" % (obspy.UTCDateTime(t0).strftime("%Y-%m-%d %H:%M"),
                            obspy.UTCDateTime(t0 + 86400).strftime("%Y-%m-%d %H:%M"))


def day_windows(traces, cfg):
    """Processing days that hold data: {label: (t0, [trace rows])}.  data.start_date / end_date select by the date
    of the label (data.day_label), so one date is one processing day."""
    import obspy
    h = cfg["data"]["day_start_utc_hour"]
    off = label_offset_s(cfg["data"].get("day_label"))
    lo = str(cfg["data"]["start_date"]) if cfg["data"]["start_date"] else None
    hi = str(cfg["data"]["end_date"]) if cfg["data"]["end_date"] else None
    days = defaultdict(list)
    for tr in traces:
        u = obspy.UTCDateTime(tr["start"])
        t0 = obspy.UTCDateTime(u.year, u.month, u.day, h).timestamp
        if t0 > tr["start"]:
            t0 -= 86400
        while t0 < tr["end"]:
            d = obspy.UTCDateTime(t0 + off).strftime("%Y-%m-%d")
            if (lo is None or d >= lo) and (hi is None or d <= hi):
                days[t0].append(tr)
            t0 += 86400
    return {day_label(t0, offset_s=off): (t0, rows) for t0, rows in sorted(days.items())}


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


def _decimate_piece(d, q, sos, pad, r0):
    """Zero-phase low-pass + decimation of one contiguous piece, in blocks (memory independent of its length):
    the filtered samples at the piece indices r0, r0 + q, r0 + 2q, ..."""
    from scipy.signal import sosfiltfilt
    n = d.size
    idx = np.arange(r0, n, q)
    out = np.empty(idx.size, np.float32)
    block = max(q * 4096, (1 << 22) // q * q)
    for a in range(0, n, block):
        b = min(n, a + block)
        lo, hi = max(0, a - pad), min(n, b + pad)
        y = sosfiltfilt(sos, d[lo:hi])
        k0 = int(np.searchsorted(idx, a))
        k1 = int(np.searchsorted(idx, b))
        out[k0:k1] = y[idx[k0:k1] - lo]
    return out


def snap_to_grid(tr, t0, fs):
    """Put a trace on the UTC sample grid of its rate (t0 + k / fs), in place, the same way for the cleaned archive
    and for direct reading: within 1e-3 sample it is on the grid; within 0.05 sample its start time is moved to the
    nearest grid point (no change of the samples); beyond that it is resampled by Lanczos interpolation (a = 20)
    starting at the first grid point inside the record.  -> "on_grid" | "shifted" | "interpolated"."""
    import obspy
    off = (tr.stats.starttime.timestamp - t0) * fs
    k0 = int(round(off))
    if abs(off - k0) <= 1e-3:
        return "on_grid"
    if abs(off - k0) <= 0.05:
        tr.stats.starttime = obspy.UTCDateTime(t0 + k0 / fs)
        return "shifted"
    tr.data = np.asarray(tr.data, np.float64)
    tr.interpolate(fs, starttime=obspy.UTCDateTime(t0 + int(np.ceil(off)) / fs), method="lanczos", a=20)
    return "interpolated"


def read_day(rows, t0, fs):
    """Traces of one station for the day starting at t0 -> (data at fs, valid mask) or (None, None).

    Samples placed on the UTC grid of their rate (whole-sample shift; Lanczos interpolation only when the offset
    exceeds 0.05 sample; the cleaned archive is already on the grid), non-finite samples treated as gaps, the sign
    reversed for a vertical sensor pointing down (cmpinc 180; the cleaned archive is already corrected), the raw rate
    must be a whole multiple of fs, the median of the day removed, zero-phase
    8th-order Butterworth low-pass at 0.8 x the new Nyquist frequency and decimation, piece by piece and in blocks
    (so that the memory does not grow with the raw sampling rate).  A decimated sample is valid only if the raw data
    within +-max(1 s, 50 output samples) of it were all present."""
    import obspy
    from scipy.signal import butter
    raw_fs = rows[0]["fs"]
    q = int(round(raw_fs / fs))
    if q < 1 or abs(raw_fs / fs - q) > 1e-6:
        raise ValueError("sampling rate %g Hz is not a whole multiple of correlation.fs = %g Hz" % (raw_fs, fs))
    nout = int(round(86400 * fs))
    polarity = -1.0 if is_down(rows[0]) else 1.0
    pad = 200 * q                                   # filter run-in (raw samples) on each side of a block
    guard = int(np.ceil(max(1.0, 50.0 / fs) * raw_fs))
    extra = (pad + guard) / raw_fs + 1.0            # read this much beyond the day so its edges are not data ends
    want = (rows[0]["net"], rows[0].get("sta_header", rows[0]["sta"]), rows[0]["location"], rows[0]["channel"], raw_fs)
    st = obspy.Stream()
    for path in sorted({r["file"] for r in rows}):
        part, _ = read_quiet(path, format=rows[0].get("format") or None, starttime=obspy.UTCDateTime(t0 - extra),
                             endtime=obspy.UTCDateTime(t0 + 86400 + extra))
        for tr in part:
            s = tr.stats
            if (s.network, s.station, s.location, s.channel, s.sampling_rate) == want:
                st.append(tr)
    pieces = []
    if len(st) > 1:                                 # pieces in different files that join up are one piece
        try:
            st.merge(method=-1)
            st.merge(method=0, fill_value=None)
        except Exception:
            pass
    for tr in st.split():
        if np.issubdtype(tr.data.dtype, np.floating) and not np.isfinite(tr.data).all():
            tr.data = np.ma.masked_invalid(tr.data)      # NaN / inf samples are gaps
            parts = tr.split()
        else:
            parts = [tr]
        for tr2 in parts:
            if tr2.stats.npts < 2:
                continue
            snap_to_grid(tr2, t0, raw_fs)
            k0 = int(round((tr2.stats.starttime.timestamp - t0) * raw_fs))
            pieces.append((k0, np.asarray(tr2.data) * polarity))
    # a piece shorter than the guard band on both sides gives no valid sample (and is too short for the zero-phase
    # filter): it is left out, the other pieces of the day are kept
    pieces = [(a, d) for a, d in pieces if d.size > 2 * guard + 1 and a + d.size > 0 and a < nout * q]
    if not pieces:
        return None, None
    tot = sum(d.size for _, d in pieces)
    stride = max(1, tot // 1000000)
    med = float(np.median(np.concatenate([d[::stride] for _, d in pieces]).astype(np.float64)))
    sos = butter(8, 0.8 * fs / 2, btype="low", fs=raw_fs, output="sos")
    y = np.zeros(nout, np.float32)
    ok = np.zeros(nout, bool)
    for a, d in sorted(pieces, key=lambda t: t[0]):
        j0 = max(0, -(-a // q))                     # first output sample whose raw index j*q lies in the piece
        j1 = min(nout, (a + d.size - 1) // q + 1)
        if j1 <= j0:
            continue
        r0 = j0 * q - a
        lo = max(0, r0 - pad - guard)
        hi = min(d.size, (j1 - 1) * q - a + 1 + pad + guard)
        z = _decimate_piece(d[lo:hi].astype(np.float64) - med, q, sos, pad, r0 - lo)
        y[j0:j0 + min(z.size, j1 - j0)] = z[:j1 - j0]
        g0 = max(j0, -(-(a + guard) // q))
        g1 = min(j1, (a + d.size - 1 - guard) // q + 1)
        if g1 > g0:
            ok[g0:g1] = True
    return y, ok
