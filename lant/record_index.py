"""Raw-record discovery and cached header indexing; waveform processing stays in lant.data."""

import fnmatch
import json
import os

import numpy as np


class DataError(ValueError):
    pass


def file_signature(path):
    """Cheap metadata identity for local cache invalidation; not a content checksum."""
    st = os.stat(path)
    return [st.st_size, st.st_mtime_ns, st.st_ctime_ns, st.st_dev, st.st_ino]


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


READ_ERROR = "read error"


def _scan_one(args):
    path, fmt = args
    try:
        st, nwarn = read_quiet(path, format=None if fmt == "auto" else fmt, headonly=True)
    except OSError as exc:  # the file could not be read now (disk, network share): not kept in the header cache
        return path, None, "%s %s: %s" % (READ_ERROR, type(exc).__name__, str(exc)[:150])
    except Exception as exc:  # not a seismic data file, or corrupt
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
        rows.append(
            dict(
                file=path,
                net=s.network,
                sta=s.station,
                channel=s.channel,
                location=s.location,
                fs=float(s.sampling_rate),
                start=s.starttime.timestamp,
                end=s.endtime.timestamp + 1.0 / s.sampling_rate,
                npts=int(s.npts),
                format=str(s.get("_format", "")),
                coord=coord,
                cmpinc=cmpinc,
            )
        )
    return path, rows, nwarn


def _peek_one(path):
    """The first and the last record of a miniSEED file (two small reads) -> (path, rows) with one row per channel
    seen, marked header_only; rows None when the file is not plain miniSEED of fixed record length or the two
    records belong to different stations (the file is then read in full)."""
    try:
        from obspy.io.mseed.util import get_record_information

        first = get_record_information(path)
        size, reclen = first["filesize"], first["record_length"]
        if not reclen or size % reclen:
            return path, None
        last = get_record_information(path, offset=size - reclen)
    except Exception:
        return path, None
    if (first["network"], first["station"]) != (last["network"], last["station"]):
        return path, None
    rows = []
    for rec in {(first["channel"], first["location"]): first, (last["channel"], last["location"]): last}.values():
        fs = float(rec["samp_rate"])
        rows.append(
            dict(
                file=path,
                net=first["network"],
                sta=first["station"],
                channel=rec["channel"],
                location=rec["location"],
                fs=fs,
                start=first["starttime"].timestamp,
                end=last["endtime"].timestamp + 1.0 / fs,
                npts=None,
                format="MSEED",
                coord=None,
                cmpinc=None,
                header_only="first and last record",
            )
        )
    return path, rows


def build_index(cfg):
    """Scan the headers of every candidate file -> index dict (all channels).

    A miniSEED file whose first and last records belong to a station that cannot be one of the station table
    (station_match.could_match) is not read in full: its index rows come from those two records (header_only), so
    that the header codes seen stay the same for the matching, and the clean step never uses its records."""
    from multiprocessing import Pool

    from . import resources
    from .processes import pool_results

    d = cfg["data"]
    files = find_files(d["raw_dir"], d["file_patterns"], d["ignore_patterns"])
    if not files:
        raise DataError("no candidate files under %s" % d["raw_dir"])
    # Metadata reuse is an I/O optimization, not a content checksum. Nanosecond mtime and ctime detect finer edits
    # than the former integer seconds, subject to the filesystem's metadata resolution. Legacy signatures miss once.
    cache_path = os.path.join(cfg["project"]["output_dir"], "work", "01_clean", "header_cache.json")
    cache, sig = {}, {}
    try:
        c = json.load(open(cache_path))
        if c.get("format") == d["format"] and isinstance(c.get("files"), dict):
            cache = c["files"]
    except (OSError, ValueError, KeyError):
        pass
    for f in files:
        sig[f] = file_signature(f)
    from .station_match import could_match

    relevant = could_match(cfg) if d["format"] in ("auto", "MSEED") else None

    def kept_peek(entry):  # a cached first-and-last-record entry still stands while its station stays irrelevant
        return relevant is not None and entry[1] and not any(relevant(r["net"], r["sta"]) for r in entry[1])

    todo = [
        f for f in files if f not in cache or cache[f][0] != sig[f] or (len(cache[f]) > 3 and not kept_peek(cache[f]))
    ]
    new = {}
    if todo and relevant is not None:
        from multiprocessing import Pool

        from . import resources
        from .processes import pool_results

        plan = resources.plan(cfg, len(todo), None, "files (first and last records)")
        with Pool(plan["parallel"], initializer=resources.limit_threads, initargs=(1,)) as pool:
            for path, rows in pool_results(pool, _peek_one, todo, chunksize=64, ordered=True):
                if rows and not any(relevant(r["net"], r["sta"]) for r in rows):
                    new[path] = (path, rows, 0, "peek")
        todo = [f for f in todo if f not in new]
    peeked = sum(1 for f in files if f in new or (f not in todo and len(cache.get(f, ())) > 3))
    print(
        "  headers: %d cached, %d to scan%s"
        % (len(files) - len(todo) - len(new), len(todo), "; %d files of stations not in the table: first and last "
           "record only" % peeked if peeked else ""), flush=True,
    )  # fmt: skip
    if todo:
        plan = resources.plan(cfg, len(todo), None, "files (header scan)")
        with Pool(plan["parallel"], initializer=resources.limit_threads, initargs=(1,)) as pool:
            tasks = [(f, d["format"]) for f in todo]
            new.update({r[0]: r for r in pool_results(pool, _scan_one, tasks, chunksize=16, ordered=True)})
            resources.note_workers(w.pid for w in pool._pool)
    res = [new[f] if f in new else (f, *cache[f][1:]) for f in files]
    try:
        if new or set(cache) != set(files):
            os.makedirs(os.path.dirname(cache_path), exist_ok=True)
            with open(cache_path + ".tmp", "w") as fh:
                kept = (r for r in res if not (r[1] is None and str(r[2]).startswith(READ_ERROR)))  # retried
                json.dump(dict(format=d["format"], files={r[0]: [sig[r[0]], *r[1:]] for r in kept}), fh)
            os.replace(cache_path + ".tmp", cache_path)
    except OSError:
        pass
    traces, unreadable, skipped = [], [], []
    for path, rows, err, *_ in res:
        if rows is None:
            unreadable.append((path, err))
        else:
            traces.extend(rows)
            if err:
                skipped.append((path, int(err)))
    return dict(raw_dir=d["raw_dir"], n_files=len(files), unreadable=unreadable, traces=traces, reader_warnings=skipped,
                header_only=peeked)  # fmt: skip
