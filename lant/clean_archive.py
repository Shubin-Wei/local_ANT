"""Clean and reuse one station-day archive; whole-step orchestration stays in lant.clean."""

import json
import os

import numpy as np

from . import data as D
from . import pulse as PU
from .record_index import file_signature

ARCHIVE_VERSION = 2  # raise when _clean_station_day changes what it writes: older station-days are redone


def _task_key(name, rows, t0, opt):
    """Fingerprint of everything a cleaned station-day depends on: its raw files (path, size, modification time,
    the time span used), the day and the settings.  The same key = the archived file can be used again."""
    import hashlib

    files = []
    for r in sorted(rows, key=lambda r: (r["file"], r["start"])):
        files.append([r["file"], file_signature(r["file"]), r["start"], r["end"], r.get("clip"), r.get("cmpinc")])
    resp = opt.get("response")
    from .runtime_identity import code_digest

    what = [
        ARCHIVE_VERSION,
        code_digest("clean"),
        name,
        t0,
        files,
        opt["format"],
        opt["min_len"],
        opt["encoding"],
        opt["codes"][name],
        opt.get("label_offset", 0.0),
        [resp["path"], resp["ids"][name], resp["pre_filt"], resp["identity"]] if resp else None,
    ]
    return hashlib.sha1(json.dumps(what, sort_keys=True, default=str).encode()).hexdigest()


def _clean_station_day(task):
    name, rows, t0, opt = task
    import obspy

    fs = rows[0]["fs"]
    want = (rows[0]["net"], rows[0].get("sta_header", rows[0]["sta"]), rows[0]["location"], rows[0]["channel"], fs)
    day = D.day_label(t0, offset_s=opt.get("label_offset", 0.0))
    d = os.path.join(opt["archive"], name)
    fn = os.path.join(d, "%s.%s.%s.%s.mseed" % (name, want[2], want[3], day))
    key = _task_key(name, rows, t0, opt)
    try:  # done by an earlier run with the same inputs and settings?
        old = json.load(open(fn + ".done.json"))
        if old["key"] == key and (not old["rows"] or os.path.exists(fn)):
            for row in old["rows"]:
                row["file"] = fn
            return name, old["rows"], dict(old["stat"], reused=1)
    except (OSError, ValueError, KeyError):
        pass

    def finish(rows_out, stat):
        os.makedirs(d, exist_ok=True)
        with open(fn + ".done.json.tmp", "w") as fh:
            json.dump(dict(key=key, rows=rows_out, stat=stat), fh, default=float)
        os.replace(fn + ".done.json.tmp", fn + ".done.json")
        return name, rows_out, stat

    st = obspy.Stream()
    for path in sorted({r["file"] for r in rows}):
        try:
            part, _ = D.read_quiet(
                path,
                format=None if opt["format"] == "auto" else opt["format"],
                starttime=obspy.UTCDateTime(t0),
                endtime=obspy.UTCDateTime(t0 + 86400 - 0.5 / fs),
            )
        except Exception as exc:
            return name, None, dict(error="%s: %s" % (path, exc))
        for tr in part:
            s = tr.stats
            if (s.network, s.station, s.location, s.channel, s.sampling_rate) == want:
                st.append(tr)
    clips = [r["clip"] for r in rows if r.get("clip")]
    if clips:  # a moved node: only its time at this point
        st.trim(
            obspy.UTCDateTime(max(t0, min(c[0] for c in clips))),
            obspy.UTCDateTime(min(t0 + 86400, max(c[1] for c in clips))),
            nearest_sample=False,
        )
        st = obspy.Stream([tr for tr in st if tr.stats.npts > 0])
    stat = dict(
        station=name,
        day=day,
        traces_in=len(st),
        shifted=0,
        interpolated=0,
        duplicates_or_identical_overlaps=0,
        conflicting_overlap_s=0.0,
        dropped_short=0,
        hours=0.0,
        gaps=0,
        nonfinite_samples=0,
        polarity_reversed=int(D.is_down(rows[0])),
        integer=all(np.issubdtype(tr.data.dtype, np.integer) for tr in st),
    )
    if not len(st):
        return finish(None, stat)
    for tr in st:  # 0. non-finite samples are gaps; a downward sensor is turned over
        if np.issubdtype(tr.data.dtype, np.floating):
            bad = ~np.isfinite(tr.data)
            if bad.any():
                stat["nonfinite_samples"] += int(bad.sum())
                tr.data = np.ma.masked_invalid(tr.data)
        if stat["polarity_reversed"]:
            tr.data = -tr.data
    if stat["nonfinite_samples"]:
        st = st.split()
    for tr in st:  # 1. onto the UTC sample grid of the day (data.snap_to_grid)
        kind = D.snap_to_grid(tr, t0, fs)
        if kind == "shifted":
            stat["shifted"] += 1
        elif kind == "interpolated":
            stat["interpolated"] += 1
            stat["integer"] = False
    if not stat["integer"]:
        for tr in st:
            tr.data = tr.data.astype(np.float32)
    n0 = len(st)
    st.merge(method=-1)  # 2. exact duplicates / identical overlaps
    stat["duplicates_or_identical_overlaps"] = n0 - len(st)
    for g in st.get_gaps():
        if g[6] < 0:
            stat["conflicting_overlap_s"] += float(-g[6])
    st.merge(method=0, fill_value=None)  # remaining overlaps with different values -> masked
    out = obspy.Stream()
    for tr in st.split():
        if tr.stats.npts / fs < opt["min_len"]:
            stat["dropped_short"] += 1
            continue
        out.append(tr)
    if not len(out):
        return finish(None, stat)
    stat["hours"] = round(sum(tr.stats.npts for tr in out) / fs / 3600, 4)
    stat["gaps"] = len(out) - 1
    fold = PU.fold_1s([(int(round((tr.stats.starttime.timestamp - t0) * fs)), tr.data) for tr in out], fs)
    if fold:
        stat.update(pulse_ratio=round(fold["ratio"], 3), pulse_power=fold["power"], pulse_limit=round(fold["limit"], 3))
    if opt.get("response"):  # 3. instrument response (gap-free pieces)
        from . import response as R

        try:
            inv = R.load(opt["response"]["path"])
            for tr in out:
                R.remove(tr, inv, opt["response"]["ids"][name], opt["response"]["pre_filt"])
        except Exception as exc:
            return name, None, dict(error="%s %s: response removal failed (%s)" % (name, stat["day"], exc))
        stat["integer"] = False
    net, sta, loc, cha = opt["codes"][name]  # miniSEED header codes (limited lengths), see storage_codes
    for tr in out:
        tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel = net, sta, loc, cha
    os.makedirs(d, exist_ok=True)
    enc = opt["encoding"]
    if enc == "auto":
        enc = "STEIM2" if stat["integer"] else "FLOAT32"
    for tr in out:
        if enc in ("STEIM1", "STEIM2", "INT32"):
            if not np.issubdtype(tr.data.dtype, np.integer):
                enc = "FLOAT32"
        if enc == "FLOAT32":
            tr.data = tr.data.astype(np.float32)
    out.write(fn + ".tmp", format="MSEED", encoding=enc, reclen=4096)
    back = obspy.read(fn + ".tmp", headonly=True)  # what was written is what the index says
    if sorted((t.id, t.stats.npts) for t in back) != sorted((t.id, t.stats.npts) for t in out):
        return (
            name,
            None,
            dict(
                error="%s: the written file does not read back as written (%s / %s)"
                % (fn, sorted(t.id for t in back), sorted(t.id for t in out))
            ),
        )
    os.replace(fn + ".tmp", fn)
    rows_out = [
        dict(
            file=fn,
            station=name,
            net=net,
            sta=sta,
            channel=cha,
            location=loc,
            fs=fs,
            start=tr.stats.starttime.timestamp,
            end=tr.stats.endtime.timestamp + 1 / fs,
            npts=int(tr.stats.npts),
        )
        for tr in out
    ]
    return finish(rows_out, stat)
