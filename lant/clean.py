"""Step 01: raw-data cleaning.  Any directory layout, any file names, any format ObsPy reads.

1. index: the headers of every candidate file under data.raw_dir (data.format = auto: MSEED, SAC, GSE2, SEG2,
   SEISAN, ...); records of other channels (data.channel) are ignored; files that are not seismic data are listed;
2. the vertical component of every station (data.channel = auto: channel ending in Z, or SAC cmpinc = 0);
3. names: the station table (any common form) is matched to the header codes (lant/stations.py) and written as
   01_clean/stations.csv; coordinates found in SAC headers go to stations_from_headers.csv and are compared;
   per station one sensor: one location/channel/sampling rate (the one with most data); the others are listed;
   the "auto" parameters are derived from the data (lant/auto.py, records/auto_parameters.json);
4. per station and processing day (data.day_start_utc_hour), the records are
     - placed on the UTC sample grid of their rate (start shifted when the offset is <= 0.05 sample, Lanczos
       interpolation when larger; both counted),
     - merged: exact duplicates and identical overlaps removed; overlaps with different values are NOT
       resolved by choosing one of them - the overlapping samples are left out (a gap), and counted;
     - traces shorter than data.min_trace_length_s dropped;
   gaps stay gaps (nothing is filled or interpolated across a gap); the sample values are not changed otherwise
   (no filtering, no instrument correction, no resampling);
5. writes 01_clean/archive/NET.STA/NET.STA.LOC.CHA.YYYYDDD.mseed (one file per station-day; clean.encoding) and
   archive_index.json, which every later step reads; report.txt / station_days.csv / coverage.png list what was
   found and done.  With clean.write_archive = false nothing is written and the later steps read the raw files
   through the index (steps 2-4 are then applied on the fly, except the overlap rule).
"""
import csv
import json
import os
import shutil
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np

from . import data as D
from . import paths


def _clean_station_day(task):
    name, rows, t0, opt = task
    import obspy
    fs = rows[0]["fs"]
    want = (rows[0]["net"], rows[0]["sta"], rows[0]["location"], rows[0]["channel"], fs)
    st = obspy.Stream()
    for path in sorted({r["file"] for r in rows}):
        try:
            part = obspy.read(path, format=None if opt["format"] == "auto" else opt["format"],
                              starttime=obspy.UTCDateTime(t0), endtime=obspy.UTCDateTime(t0 + 86400 - 0.5 / fs))
        except Exception as exc:
            return name, None, dict(error="%s: %s" % (path, exc))
        for tr in part:
            s = tr.stats
            if (s.network, s.station, s.location, s.channel, s.sampling_rate) == want:
                st.append(tr)
    stat = dict(station=name, day=D.day_label(t0), traces_in=len(st), shifted=0, interpolated=0,
                duplicates_or_identical_overlaps=0, conflicting_overlap_s=0.0, dropped_short=0, hours=0.0, gaps=0,
                integer=all(np.issubdtype(tr.data.dtype, np.integer) for tr in st))
    if not len(st):
        return name, None, stat
    for tr in st:                                    # 1. onto the UTC sample grid of the day
        off = (tr.stats.starttime.timestamp - t0) * fs
        k0 = int(round(off))
        if abs(off - k0) > 1e-3:
            if abs(off - k0) <= 0.05:
                tr.stats.starttime = obspy.UTCDateTime(t0 + k0 / fs)
                stat["shifted"] += 1
            else:
                tr.data = tr.data.astype(np.float64)
                tr.interpolate(fs, starttime=obspy.UTCDateTime(t0 + int(np.ceil(off)) / fs), method="lanczos", a=20)
                stat["interpolated"] += 1
                stat["integer"] = False
    if not stat["integer"]:
        for tr in st:
            tr.data = tr.data.astype(np.float32)
    n0 = len(st)
    st.merge(method=-1)                              # 2. exact duplicates / identical overlaps
    stat["duplicates_or_identical_overlaps"] = n0 - len(st)
    for g in st.get_gaps():
        if g[6] < 0:
            stat["conflicting_overlap_s"] += float(-g[6])
    st.merge(method=0, fill_value=None)              # remaining overlaps with different values -> masked
    out = obspy.Stream()
    for tr in st.split():
        if tr.stats.npts / fs < opt["min_len"]:
            stat["dropped_short"] += 1
            continue
        out.append(tr)
    if not len(out):
        return name, None, stat
    stat["hours"] = round(sum(tr.stats.npts for tr in out) / fs / 3600, 4)
    stat["gaps"] = len(out) - 1
    net, sta = name.split(".", 1)
    for tr in out:
        tr.stats.network, tr.stats.station = net, sta
    d = os.path.join(opt["archive"], name)
    os.makedirs(d, exist_ok=True)
    fn = os.path.join(d, "%s.%s.%s.%s.mseed" % (name, want[2], want[3], stat["day"]))
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
    os.replace(fn + ".tmp", fn)
    rows_out = [dict(file=fn, station=name, net=net, sta=sta, channel=want[3], location=want[2], fs=fs,
                     start=tr.stats.starttime.timestamp, end=tr.stats.endtime.timestamp + 1 / fs, npts=int(tr.stats.npts))
                for tr in out]
    return name, rows_out, stat


def header_coordinates(traces):
    c = defaultdict(list)
    for t in traces:
        if t.get("coord"):
            c[(t["net"], t["sta"])].append(t["coord"])
    return {k: tuple(np.median(np.array(v, float), axis=0)) for k, v in c.items()}


def match(cfg, traces, out, rep):
    """Header codes -> canonical NET.STA of the station table; writes 01_clean/stations.csv."""
    from . import stations as ST
    d = cfg["data"]
    hc = header_coordinates(traces)
    if hc:
        coords = {"%s.%s" % ((k[0] or d["network_default"] or "XX"), k[1]): dict(lon=v[0], lat=v[1], elev=v[2]) for k, v in hc.items()}
        D.write_table(os.path.join(out, "stations_from_headers.csv"), coords)
        rep.append("  coordinates found in the headers of %d stations -> stations_from_headers.csv" % len(hc))
    if not d["station_file"]:
        raise D.DataError("data.station_file is not set%s" % (
            "; the headers hold coordinates: check 01_clean/stations_from_headers.csv and use it as the station file" if hc else ""))
    try:
        rows, how = ST.read_table(d["station_file"])
    except ST.StationError as exc:
        raise D.DataError(str(exc))
    rep.append("  station table %s: %d rows; %s" % (d["station_file"], len(rows), how))
    codes = {(t["net"], t["sta"]) for t in traces}
    forced = {}
    for k, v in (d["rename"] or {}).items():
        n, s_ = k.split(".", 1)
        forced[(n, s_)] = v
    if d["network_default"]:
        for t in traces:
            if not t["net"]:
                forced.setdefault(("", t["sta"]), "%s.%s" % (d["network_default"], t["sta"]))
    mapping, coords, mrep = ST.match_names(rows, codes - set(forced))
    for code, name in forced.items():
        if code in codes:
            mapping[code] = name
    table_names = set(coords)
    lost = [n for n in set(forced.values()) if n not in table_names]
    if lost:
        raise D.DataError("data.rename / network_default give names that are not in the station table: %s" % ", ".join(sorted(lost)))
    if not mapping:
        ex_t = next(iter(table_names))
        ex_d = "%s.%s" % next(iter(codes))
        raise D.DataError("no station of the table could be matched to the record headers (table e.g. %s, headers e.g. %s); "
                          "give data.rename" % (ex_t, ex_d))
    used = {name: coords[name] for name in set(mapping.values())}
    D.write_table(D.canonical_table(cfg["project"]["output_dir"]), used)
    for line in mrep["renamed"][:50]:
        rep.append("    name matched: " + line)
    if len(mrep["renamed"]) > 50:
        rep.append("    ... %d more matched names" % (len(mrep["renamed"]) - 50))
    for key, txt in (("ambiguous", "WARNING ambiguous names (not used): "), ("table_without_data", "WARNING table stations without data: "),
                     ("data_without_table", "WARNING stations in the data but not in the table (not used): ")):
        if mrep[key]:
            rep.append("  %s%d: %s" % (txt, len(mrep[key]), ", ".join(mrep[key][:30]) + (" ..." if len(mrep[key]) > 30 else "")))
    mism = []
    from .signal import geodesic_km
    for code, (lo, la, el) in hc.items():
        if code in mapping:
            c = used[mapping[code]]
            dkm = float(geodesic_km(lo, la, c["lon"], c["lat"])[0])
            if dkm > 0.05:
                mism.append((mapping[code], round(dkm, 3)))
    if mism:
        rep.append("  WARNING %d stations differ by more than 50 m between the headers and the table (the table is used): %s"
                   % (len(mism), ", ".join("%s %.2f km" % m for m in sorted(mism, key=lambda m: -m[1])[:20])))
    for t in traces:
        t["station"] = mapping.get((t["net"], t["sta"]))
    return used, dict(matched=len(mapping), renamed=len(mrep["renamed"]), ambiguous=mrep["ambiguous"],
                      table_without_data=mrep["table_without_data"], data_without_table=mrep["data_without_table"],
                      coordinate_mismatch_km=dict(mism))


def run(cfg):
    from . import auto
    out = paths.mk(paths.P(cfg, "01_clean"))
    nproc = int(cfg["clean"]["nproc"])
    print("  scanning headers under %s ..." % cfg["data"]["raw_dir"], flush=True)
    index = D.build_index(cfg, nproc=nproc)
    json.dump(index, open(os.path.join(out, "raw_index.json"), "w"))
    fmts = Counter(t["format"] for t in index["traces"])
    rep = ["raw-data cleaning (lant clean)", "  raw directory %s" % index["raw_dir"],
           "  candidate files %d; not seismic data or unreadable %d; records %d; formats: %s"
           % (index["n_files"], len(index["unreadable"]), len(index["traces"]),
              ", ".join("%s %d" % kv for kv in fmts.most_common()))]
    for pth, err in index["unreadable"][:10]:
        rep.append("    skipped %s (%s)" % (pth, err[:80]))
    if not index["traces"]:
        raise D.DataError("no seismic record found under %s" % index["raw_dir"])
    chans, notes = D.choose_channel(index["traces"], cfg["data"]["channel"])
    rep.append("  vertical component (data.channel = %s): %d stations; channels used: %s"
               % (cfg["data"]["channel"], len(chans), ", ".join(sorted({r["channel"] for v in chans.values() for r in v}))))
    rep += ["    " + n for n in notes]
    traces = [r for v in chans.values() for r in v]
    used, msum = match(cfg, traces, out, rep)
    rows, several = {}, {}
    by = defaultdict(list)
    for t in traces:
        if t["station"]:
            by[t["station"]].append(t)
    for s, rr in by.items():
        rows[s], (kept, other) = D.best_sensor(rr)
        if other:
            several[s] = dict(kept=kept, other_hours=other)
            rep.append("  %s: several sensors/rates, kept %s; not used: %s" % (s, kept, other))
    decisions, geom = auto.infer(cfg, {k: used[k] for k in rows}, [r for v in rows.values() for r in v])
    auto.write(cfg, decisions, geom)
    cfg = auto.apply(cfg)
    rep.append("  parameters from the data (records/auto_parameters.json):")
    for k, v in decisions.items():
        if v["reason"] != "given in the configuration":
            rep.append("    %s = %s  (%s)" % (k, json.dumps(v["value"]) if not isinstance(v["value"], str) else v["value"], v["reason"]))
    ex = set(D.resolve_names(cfg["data"]["exclude_stations"], used, "data.exclude_stations"))
    rows = {k: v for k, v in rows.items() if k not in ex}
    days = D.day_windows([r for v in rows.values() for r in v], cfg)
    nsd = sum(len({t["station"] for t in trs}) for _, trs in days.values())
    est_bytes = nsd * max(r[0]["fs"] for r in rows.values()) * 86400 * 4
    rep.append("  %d stations used (%d excluded by data.exclude_stations), %d processing days, %d station-days"
               % (len(rows), len(ex), len(days), nsd))
    archive = os.path.join(out, "archive")
    summary = dict(files=index["n_files"], unreadable=len(index["unreadable"]), formats=dict(fmts), stations=len(rows),
                   days=len(days), station_days=nsd, channel_notes=notes, names=msum, several_sensors=several,
                   sampling_rates=sorted({r[0]["fs"] for r in rows.values()}),
                   record_span_utc=[D.day_label(min(t["start"] for t in traces)), D.day_label(max(t["end"] for t in traces))],
                   auto_parameters={k: v["value"] for k, v in decisions.items()})
    if not cfg["clean"]["write_archive"]:
        idx = dict(source="raw", traces=[r for v in rows.values() for r in v
                                         if r["end"] - r["start"] >= cfg["data"]["min_trace_length_s"]])
        json.dump(idx, open(os.path.join(out, "archive_index.json"), "w"))
        rep.append("  clean.write_archive is false: later steps read the raw files")
        open(os.path.join(out, "report.txt"), "w").write("\n".join(rep) + "\n")
        print("\n".join(rep))
        return summary
    free = shutil.disk_usage(out).free
    if free < cfg["clean"]["min_free_disk_factor"] * est_bytes:
        raise D.DataError("the cleaned archive needs up to %.1f GB, %.1f GB are free on %s: free space, or set "
                          "clean.write_archive: false" % (est_bytes / 1e9, free / 1e9, out))
    opt = dict(format=cfg["data"]["format"], min_len=float(cfg["data"]["min_trace_length_s"]), archive=archive,
               encoding=cfg["clean"]["encoding"])
    tasks = []
    for lab, (t0, trs) in days.items():
        bys = defaultdict(list)
        for r in trs:
            bys[r["station"]].append(r)
        for s, rr in bys.items():
            tasks.append((s, rr, t0, opt))
    stats, new_rows = [], []
    with Pool(nproc) as pool:
        for k, (name, rr, stat) in enumerate(pool.imap_unordered(_clean_station_day, tasks, chunksize=4)):
            stats.append(stat)
            if rr:
                new_rows.extend(rr)
            if (k + 1) % 200 == 0:
                print("  cleaned %d / %d station-days" % (k + 1, len(tasks)), flush=True)
    errors = [s for s in stats if "error" in s]
    for e in errors:
        rep.append("  ERROR %s" % e["error"])
    stats = [s for s in stats if "error" not in s]
    keys = ["station", "day", "hours", "gaps", "traces_in", "shifted", "interpolated", "duplicates_or_identical_overlaps",
            "conflicting_overlap_s", "dropped_short"]
    with open(os.path.join(out, "station_days.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for s in sorted(stats, key=lambda s: (s["station"], s["day"])):
            w.writerow([s[k] for k in keys])
    tot = {k: sum(s[k] for s in stats) for k in keys[3:]}
    rep.append("  archive %s: %d station-day files, %.1f station-hours; gaps %d; records shifted onto the grid %d, "
               "interpolated %d; duplicates/identical overlaps removed %d; conflicting overlaps left out %.1f s; "
               "pieces shorter than %.0f s dropped %d"
               % (archive, len({(s['station'], s['day']) for s in stats if s['hours'] > 0}), sum(s["hours"] for s in stats),
                  tot["gaps"], tot["shifted"], tot["interpolated"], tot["duplicates_or_identical_overlaps"],
                  tot["conflicting_overlap_s"], opt["min_len"], tot["dropped_short"]))
    json.dump(dict(source="archive", traces=new_rows), open(os.path.join(out, "archive_index.json"), "w"))
    open(os.path.join(out, "report.txt"), "w").write("\n".join(rep) + "\n")
    coverage_figure(stats, sorted(days), os.path.join(out, "coverage.png"))
    summary.update(archive_station_days=len(stats), station_hours=round(sum(s["hours"] for s in stats), 1),
                   errors=len(errors), totals=tot)
    print("\n".join(rep))
    return summary


def coverage_figure(stats, days, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    stas = sorted({s["station"] for s in stats})
    M = np.zeros((len(stas), len(days)))
    for s in stats:
        M[stas.index(s["station"]), days.index(s["day"])] = s["hours"]
    fig, ax = plt.subplots(figsize=(max(6, 0.25 * len(days) + 3), max(4, 0.12 * len(stas) + 1.5)))
    im = ax.imshow(M, aspect="auto", cmap="viridis", vmin=0, vmax=24, interpolation="nearest")
    ax.set_xticks(range(len(days)))
    ax.set_xticklabels([d[-3:] for d in days], rotation=90, fontsize=6)
    step = max(1, len(stas) // 60)
    ax.set_yticks(range(0, len(stas), step))
    ax.set_yticklabels(stas[::step], fontsize=5)
    ax.set_xlabel("Processing day (day of year)")
    ax.set_ylabel("Station")
    fig.colorbar(im, ax=ax, label="Hours of data")
    ax.set_title("Data coverage after cleaning", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
