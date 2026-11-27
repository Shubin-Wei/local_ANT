"""Raw records: the index of the files, which records belong to which station, and reading one station-day.

Index: every file under data.raw_dir is read once with ObsPy (headers only; any format ObsPy recognises, any folder
layout).  The index is the table 1_cc/record_index.csv (file, network, station, location, channel, sampling rate,
start, end); a file whose size and modification time did not change is not read again.  Files ObsPy cannot read are
listed in the log and skipped.

Stations: a record belongs to a station of the table when its header network and station codes are the table's
(NET.STA), or, for a table without a network column, when its station code is.  Only channels matching
data.channel (default "*Z", the vertical component) are used; a station with records on several location / channel
/ rate combinations uses the one with most data (the log lists the others).  Names that do not match are reported,
never guessed.

Station-day (read_day): a UTC day 00:00-24:00, on the exact sample grid of the raw rate; gaps stay gaps (no
interpolation, no zero filling); the median is removed; the response is removed when data.response is given (to
ground velocity, ObsPy's remove_response, pre-filter flat in the correlation band; Bensen et al. 2007 sect. 2.1);
zero-phase 8th-order Butterworth low-pass at 0.8 x the new Nyquist frequency and decimation by a whole factor.
"""

import csv
import fnmatch
import os
from collections import defaultdict

import numpy as np

DAY = 86400.0


class DataError(ValueError):
    pass


# ----------------------------------------------------------------------------------------------------- the index
INDEX_COLUMNS = ["file", "size", "mtime", "network", "station", "location", "channel", "fs", "start", "end", "cmpinc"]


def find_files(raw_dir):
    if not raw_dir or not os.path.isdir(raw_dir):
        raise DataError("data.raw_dir %s is not a directory" % raw_dir)
    out = []
    for root, dirs, files in os.walk(raw_dir, followlinks=True):
        dirs[:] = sorted(d for d in dirs if not d.startswith("."))
        out += [os.path.join(root, f) for f in sorted(files) if not f.startswith(".")]
    return out


def read_quiet(path, **kw):
    """obspy.read without printing its warnings (recorder files often carry padding blocks, one warning each)."""
    import warnings

    import obspy

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return obspy.read(path, **kw)


def _scan(path):
    st = os.stat(path)
    try:
        traces = read_quiet(path, headonly=True)
    except Exception as exc:
        return path, None, "%s: %s" % (type(exc).__name__, str(exc)[:120])
    rows = []
    for tr in traces:
        s = tr.stats
        ci = s.sac.get("cmpinc") if "sac" in s else None
        rows.append(dict(file=path, size=st.st_size, mtime=st.st_mtime_ns, network=s.network, station=s.station,
                         location=s.location, channel=s.channel, fs=float(s.sampling_rate),
                         start=s.starttime.timestamp, end=s.endtime.timestamp + 1.0 / s.sampling_rate,
                         cmpinc="" if ci is None or ci == -12345.0 else float(ci)))  # fmt: skip
    return path, rows, None


def build_index(raw_dir, index_csv, n_cores=1, log=print):
    """Scan the headers (reusing index_csv for unchanged files) -> (rows, unreadable [(file, reason)])."""
    files = find_files(raw_dir)
    if not files:
        raise DataError("no files under %s" % raw_dir)
    old = defaultdict(list)
    if os.path.exists(index_csv):
        with open(index_csv, newline="") as fh:
            for r in csv.DictReader(fh):
                old[r["file"]].append(r)
    rows, todo = [], []
    for f in files:
        st = os.stat(f)
        prev = old.get(f)
        if prev and int(prev[0]["size"]) == st.st_size and int(prev[0]["mtime"]) == st.st_mtime_ns:
            rows += [_typed(r) for r in prev if r["network"] != "" or r["station"] != ""]
        else:
            todo.append(f)
    log("  %d files: %d indexed before, %d to read" % (len(files), len(files) - len(todo), len(todo)))
    unreadable = []
    if todo:
        from multiprocessing import Pool

        with Pool(max(1, min(n_cores, len(todo)))) as pool:
            for path, got, err in pool.imap(_scan, todo, chunksize=8):
                if got is None:
                    unreadable.append((path, err))
                else:
                    rows += got
    rows.sort(key=lambda r: (r["file"], r["start"]))
    os.makedirs(os.path.dirname(index_csv), exist_ok=True)
    with open(index_csv + ".tmp", "w", newline="") as fh:
        w = csv.DictWriter(fh, INDEX_COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in INDEX_COLUMNS})
    os.replace(index_csv + ".tmp", index_csv)
    return rows, unreadable


def _typed(r):
    out = dict(r)
    for k in ("fs", "start", "end"):
        out[k] = float(r[k])
    out["size"], out["mtime"] = int(r["size"]), int(r["mtime"])
    out["cmpinc"] = float(r["cmpinc"]) if r["cmpinc"] not in ("", None) else ""
    return out


# ------------------------------------------------------------------------------------------------ the stations
def match(rows, table, channel="*Z", exclude=()):
    """Records of each table station -> (records {name: [rows]}, notes [str]).

    The name of a station is the table's: NET.STA, or STA for a table without networks (then the header network
    must be unique for that code)."""
    notes = []
    by_code = defaultdict(list)
    for r in rows:
        if fnmatch.fnmatch(r["channel"], channel):
            by_code[(r["network"], r["station"])].append(r)
    found = {}
    for name, s in table.items():
        if name in exclude:
            continue
        if s["network"]:
            got = by_code.get((s["network"], s["station"]), [])
        else:
            nets = sorted({n for (n, c) in by_code if c == s["station"]})
            if len(nets) > 1:
                notes.append("%s: records of networks %s carry this station code; give the network in the table"
                             % (name, ", ".join(nets)))  # fmt: skip
                continue
            got = by_code.get((nets[0], s["station"]), []) if nets else []
        if not got:
            continue
        found[name] = _one_sensor(name, got, notes)
    table_codes = {(s["network"], s["station"]) if s["network"] else s["station"] for s in table.values()}
    extra = sorted({"%s.%s" % k for k in by_code if k not in table_codes and k[1] not in table_codes})
    if extra:
        notes.append("records of %d stations not in the table (not used): %s%s"
                     % (len(extra), ", ".join(extra[:20]), " ..." if len(extra) > 20 else ""))  # fmt: skip
    missing = sorted(set(table) - set(found) - set(exclude))
    if missing:
        notes.append("%d table stations without %s records: %s%s"
                     % (len(missing), channel, ", ".join(missing[:20]), " ..." if len(missing) > 20 else ""))  # fmt: skip
    return found, notes


def _one_sensor(name, rows, notes):
    """One location / channel / rate per station: the one with most data."""
    dur = defaultdict(float)
    for r in rows:
        dur[(r["location"], r["channel"], r["fs"])] += r["end"] - r["start"]
    best = max(dur, key=dur.get)
    if len(dur) > 1:
        notes.append("%s: used %s.%s %g Hz; also present: %s" % (name, best[0], best[1], best[2], ", ".join(
            "%s.%s %g Hz (%.1f h)" % (k[0], k[1], k[2], v / 3600) for k, v in dur.items() if k != best)))  # fmt: skip
    return [r for r in rows if (r["location"], r["channel"], r["fs"]) == best]


def days(records, start_date=None, end_date=None):
    """UTC days (timestamps of 00:00) that hold records, within data.start_date / end_date."""
    import obspy

    lo = obspy.UTCDateTime(str(start_date)).timestamp if start_date else -np.inf
    hi = obspy.UTCDateTime(str(end_date)).timestamp if end_date else np.inf
    out = set()
    for rows in records.values():
        for r in rows:
            t = np.floor(r["start"] / DAY) * DAY
            while t < r["end"]:
                if lo <= t <= hi:
                    out.add(t)
                t += DAY
    return sorted(out)


def day_name(t0):
    import obspy

    return obspy.UTCDateTime(t0).strftime("%Y-%m-%d")


# -------------------------------------------------------------------------------------------------- the responses
def load_responses(path):
    """Inventory of a response file or of every readable file of a folder (StationXML, RESP, dataless SEED ...)."""
    import obspy

    files = [path] if os.path.isfile(path) else sorted(
        os.path.join(r, f) for r, _, fs in os.walk(path) for f in fs)  # fmt: skip
    inv, bad = None, []
    for f in files:
        try:
            part = obspy.read_inventory(f)
        except Exception as exc:
            bad.append("%s (%s)" % (f, str(exc)[:60]))
            continue
        inv = part if inv is None else inv + part
    if inv is None:
        raise DataError(
            "data.response: no file under %s is a response file ObsPy reads: %s" % (path, "; ".join(bad[:3]))
        )
    return inv


def response_id(inv, row, t0):
    """SEED id of the response of a record, exact NET.STA.LOC.CHA; None when the inventory has none."""
    import obspy

    sid = "%s.%s.%s.%s" % (row["network"], row["station"], row["location"], row["channel"])
    try:
        inv.get_response(sid, obspy.UTCDateTime(t0))
        return sid
    except Exception:
        return None


def remove_response(tr, inv, sid, band):
    """Response of `sid` removed from a gap-free trace, to ground velocity, in place.  Pre-filter: cosine edges from
    half the lower band edge to it and from the upper edge to 1.25 x it (flat inside the band, no water level)."""
    f2, f3 = band
    nyq = 0.5 * tr.stats.sampling_rate
    pf = [0.5 * f2, f2, min(f3, 0.8 * nyq), min(1.25 * f3, nyq)]
    tr.stats.response = inv.get_response(sid, tr.stats.starttime)
    frac = min(0.05, (2.0 / pf[0]) / (tr.stats.npts * tr.stats.delta))
    tr.data = tr.data.astype(np.float64)
    tr.remove_response(output="VEL", pre_filt=pf, water_level=None, zero_mean=True, taper=True, taper_fraction=frac)
    del tr.stats.response
    return tr


# ----------------------------------------------------------------------------------------------- one station-day
def _decimate(d, q, sos, pad, r0):
    """Zero-phase low-pass of one gap-free piece in blocks (memory independent of its length) and its samples
    r0, r0 + q, r0 + 2q, ..."""
    from scipy.signal import sosfiltfilt

    idx = np.arange(r0, d.size, q)
    out = np.empty(idx.size, np.float32)
    block = max(q * 4096, (1 << 22) // q * q)
    for a in range(0, d.size, block):
        b = min(d.size, a + block)
        lo, hi = max(0, a - pad), min(d.size, b + pad)
        y = sosfiltfilt(sos, d[lo:hi])
        k0, k1 = int(np.searchsorted(idx, a)), int(np.searchsorted(idx, b))
        out[k0:k1] = y[idx[k0:k1] - lo]
    return out


def read_day(rows, t0, fs, inv=None, band=None):
    """Records of one station for the UTC day starting at t0 -> (data at fs [86400 fs], valid mask) or (None, None).

    rows: the index rows of the station (one location / channel / rate).  A decimated sample is valid only when the
    raw data within +-max(1 s, 50 output samples) of it were all present (the low-pass needs them).  A record with
    SAC cmpinc = 180 (sensor pointing down) is multiplied by -1."""
    import obspy
    from scipy.signal import butter

    raw_fs = rows[0]["fs"]
    q = int(round(raw_fs / fs))
    if q < 1 or abs(raw_fs / fs - q) > 1e-6:
        raise DataError("raw sampling rate %g Hz is not a whole multiple of cc.fs_hz = %g Hz" % (raw_fs, fs))
    nout = int(round(DAY * fs))
    pad = 200 * q
    guard = int(np.ceil(max(1.0, 50.0 / fs) * raw_fs))
    extra = (pad + guard) / raw_fs + 1.0
    want = (rows[0]["network"], rows[0]["station"], rows[0]["location"], rows[0]["channel"])
    st = obspy.Stream()
    for path in sorted({r["file"] for r in rows if r["end"] > t0 - extra and r["start"] < t0 + DAY + extra}):
        part = read_quiet(path, starttime=obspy.UTCDateTime(t0 - extra), endtime=obspy.UTCDateTime(t0 + DAY + extra))
        for tr in part:
            s = tr.stats
            if (s.network, s.station, s.location, s.channel) == want and abs(s.sampling_rate - raw_fs) < 1e-6:
                st.append(tr)
    if len(st) > 1:
        st.merge(method=-1)  # identical overlaps joined; differing overlaps become gaps below
        st.merge(method=0, fill_value=None)
    polarity = -1.0 if rows[0].get("cmpinc") not in ("", None) and abs(abs(float(rows[0]["cmpinc"])) - 180) < 1 else 1.0
    pieces = []
    for tr in st.split():
        if np.issubdtype(tr.data.dtype, np.floating) and not np.isfinite(tr.data).all():
            tr.data = np.ma.masked_invalid(tr.data)
            parts = tr.split()
        else:
            parts = [tr]
        for p in parts:
            if p.stats.npts <= 2 * guard + 1:
                continue
            off = (p.stats.starttime.timestamp - t0) * raw_fs
            if abs(off - round(off)) > 0.05:  # not on the sample grid: resample onto it (Lanczos)
                p.data = np.asarray(p.data, np.float64)
                p.interpolate(raw_fs, starttime=obspy.UTCDateTime(t0 + np.ceil(off) / raw_fs), method="lanczos", a=20)
            if inv is not None:
                sid = response_id(inv, rows[0], t0)
                if sid is None:
                    raise DataError("data.response has no response for %s.%s.%s.%s" % want)
                remove_response(p, inv, sid, band)
            k0 = int(round((p.stats.starttime.timestamp - t0) * raw_fs))
            pieces.append((k0, np.asarray(p.data, np.float64) * polarity))
    pieces = [(a, d) for a, d in pieces if a + d.size > 0 and a < nout * q]
    if not pieces:
        return None, None
    med = float(np.median(np.concatenate([d[:: max(1, d.size // 100000)] for _, d in pieces])))
    sos = butter(8, 0.8 * fs / 2, btype="low", fs=raw_fs, output="sos")
    y = np.zeros(nout, np.float32)
    ok = np.zeros(nout, bool)
    for a, d in sorted(pieces, key=lambda t: t[0]):
        j0 = max(0, -(-a // q))
        j1 = min(nout, (a + d.size - 1) // q + 1)
        if j1 <= j0:
            continue
        r0 = j0 * q - a
        lo = max(0, r0 - pad - guard)
        hi = min(d.size, (j1 - 1) * q - a + 1 + pad + guard)
        z = _decimate(d[lo:hi] - med, q, sos, pad, r0 - lo)
        y[j0 : j0 + min(z.size, j1 - j0)] = z[: j1 - j0]
        g0 = max(j0, -(-(a + guard) // q))
        g1 = min(j1, (a + d.size - 1 - guard) // q + 1)
        if g1 > g0:
            ok[g0:g1] = True
    return y, ok
