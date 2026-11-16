"""Step 01: raw-data cleaning.  Any directory layout, any file names, any format ObsPy reads.

1. index: the headers of every candidate file under data.raw_dir (data.format = auto: MSEED, SAC, GSE2, SEG2,
   SEISAN, ...); records of other channels (data.channel) are ignored; files that are not seismic data are listed;
2. the vertical component of every station (data.channel = auto: channel ending in Z, or SAC cmpinc = 0);
3. names: the station table (any common form) is matched to the header codes (lant/stations.py) and written as
   work/01_clean/stations.csv; coordinates found in SAC headers go to stations_from_headers.csv and are compared
   (lant/station_match.py, with the geometry checks);
   per station one sensor: one location/channel/sampling rate (the one with most data); the others are listed;
   instruments: with data.response the instrument responses are removed (lant/response.py; stations of different
   models can be used together); without it the records are used as they are, which is exact when every station
   has the same model - an instrument column of the station table that names several models is refused;
   the "auto" parameters follow from the station geometry, the records and data.velocity_range_kms (lant/auto.py,
   .lant/auto_parameters.json);
4. per station and processing day (data.day_start_utc_hour), the records are
     - placed on the UTC sample grid of their rate (start shifted when the offset is <= 0.05 sample, Lanczos
       interpolation when larger; both counted),
     - merged: exact duplicates and identical overlaps removed; overlaps with different values are NOT
       resolved by choosing one of them - the overlapping samples are left out (a gap), and counted;
     - traces shorter than data.min_trace_length_s dropped;
   gaps stay gaps (nothing is filled or interpolated across a gap); the sample values are not changed otherwise
   (no filtering, no resampling), except by the response removal when data.response is given;
5. writes work/01_clean/archive/NET.STA/NET.STA.LOC.CHA.YYYYDDD.mseed (one file per station-day; clean.encoding) and
   archive_index.json, which every later step reads; report.txt / station_days.csv / coverage.png list what was
   found and done.  With clean.write_archive = false nothing is written and the later steps read the raw files
   through the index (steps 2-4 are then applied on the fly, except the overlap rule).
"""

import csv
import datetime
import json
import os
import shutil
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np

from . import data as D
from . import paths
from . import pulse as PU
from . import station_match as SM
from .clean_archive import _clean_station_day
from .processes import pool_results
from .progress import track


def storage_codes(names, rows):
    """miniSEED header codes for the archive: network <= 2, station <= 5, location <= 2, channel <= 3 characters.
    The canonical NET.STA (and the original location / channel) are used when they fit; otherwise network "XX" and
    a unique station code "Z" + base-36 index, recorded in the archive index next to the canonical name."""
    out, used = {}, set()
    for i, n in enumerate(sorted(names)):
        net, sta = n.split(".", 1)
        loc, cha = rows[n][0]["location"], rows[n][0]["channel"]
        if len(net) > 2 or len(sta) > 5 or not sta.isascii() or not net.isascii():
            net, sta = "XX", "Z" + np.base_repr(i, 36).rjust(4, "0")
        loc = loc if len(loc) <= 2 else ""
        cha = cha if len(cha) <= 3 else "HHZ"
        key = (net, sta)
        if key in used:
            raise RuntimeError("two stations would get the archive code %s.%s" % key)
        used.add(key)
        out[n] = (net, sta, loc, cha)
    return out


def write_input_list(index, path):
    """work/01_clean/input_files.csv: every candidate file with size, modification time and what was read from it (the
    identity of the input of this run; a changed file shows up as a changed size or time)."""
    recs = defaultdict(list)
    for t in index["traces"]:
        recs[t["file"]].append("%s.%s.%s.%s" % (t["net"], t["sta"], t["location"], t["channel"]))
    bad = dict(index["unreadable"])
    for t in index["traces"]:
        if t.get("header_only"):
            bad[t["file"]] = "station not in the table: only the first and last records read"
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["path", "bytes", "modified_utc", "records", "not_read_because"])
        for f in sorted(set(recs) | set(bad)):
            try:
                st = os.stat(f)
                size, mt = st.st_size, datetime.datetime.utcfromtimestamp(st.st_mtime).isoformat(timespec="seconds")
            except OSError:
                size, mt = "", ""
            w.writerow([f, size, mt, " ".join(sorted(set(recs.get(f, [])))), bad.get(f, "")])


def instruments(cfg, used, rows):
    """How the instrument responses are treated -> dict(removed, ...).
    data.response given: every station needs a response (ids: station -> SEED id in the response files).
    Not given: the records are used as they are; this is exact when all stations carry the same instrument model
    (a common response cancels in the cross-coherence).  Different models do not cancel - a 5-Hz geophone against a
    broadband sensor differs in phase by 40-170 degrees between 0.5 and 10 Hz - so an instrument column of the
    station table that names several models is refused."""
    models = Counter(v["instrument"] for v in used.values() if v.get("instrument"))
    path = cfg["data"]["response"]
    if path:
        from . import response as R

        try:
            inv = R.load(path)
        except R.ResponseError as exc:
            raise D.DataError(str(exc)) from exc
        ids, miss = {}, []
        for name, rr in rows.items():
            r = rr[0]
            sid = R.lookup(
                inv, r["net"], r["sta"], r["location"], r["channel"], min(x["start"] for x in rr) + 1
            ) or R.lookup(
                inv,
                name.split(".", 1)[0],
                name.split(".", 1)[1],
                r["location"],
                r["channel"],
                min(x["start"] for x in rr) + 1,
            )
            if sid:
                ids[name] = sid
            else:
                miss.append(name)
        if miss:
            raise D.DataError(
                "data.response %s holds no (unique) response for %d of %d stations at the time of their "
                "records: %s.  Add their responses, or leave the stations out (data.exclude_stations)"
                % (path, len(miss), len(rows), ", ".join(sorted(miss)[:15]))
            )
        return dict(
            removed=True,
            path=path,
            identity=R.identity(path),
            ids=ids,
            pre_filt=[round(x, 6) for x in R.pre_filter(cfg)],
            models=dict(models),
        )
    if len(models) > 1:
        raise D.DataError(
            "the station table names different instruments (%s) and no response files are given: the "
            "responses of different models do not cancel in the cross-coherence and the pairs between "
            "them would carry phase errors of several percent of the velocity or more.  Give the "
            "response files in data.response, or process one model at a time (data.exclude_stations)"
            % ", ".join("%s: %d" % kv for kv in models.most_common())
        )
    return dict(removed=False, models=dict(models), named=sum(models.values()), stations=len(used))


def run(cfg):
    from . import auto

    out = paths.mk(paths.P(cfg, "01_clean"))
    print("  scanning headers under %s ..." % cfg["data"]["raw_dir"], flush=True)
    index = D.build_index(cfg)
    json.dump(index, open(os.path.join(out, "raw_index.json"), "w"))
    fmts = Counter(t["format"] for t in index["traces"])
    write_input_list(index, os.path.join(out, "input_files.csv"))
    rep = [
        "raw-data cleaning (lant clean)",
        "  raw directory %s" % index["raw_dir"],
        "  candidate files %d; not seismic data or unreadable %d; records %d; formats: %s"
        % (
            index["n_files"],
            len(index["unreadable"]),
            len(index["traces"]),
            ", ".join("%s %d" % kv for kv in fmts.most_common()),
        ),
    ]
    if index.get("header_only"):
        rep.append(
            "  %d of these files are of stations that cannot be in the station table (station_match.could_match): "
            "only their first and last records were read" % index["header_only"]
        )
    for pth, err in index["unreadable"][:10]:
        rep.append("    skipped %s (%s)" % (pth, err[:80]))
    rw = index.get("reader_warnings") or []
    if rw:
        rep.append(
            "  %d files gave reader warnings (%d in all; e.g. blocks that are not data, skipped by the reader); "
            "most: %s (%d)" % (len(rw), sum(n for _, n in rw), *max(rw, key=lambda it: it[1]))
        )
    if not index["traces"]:
        raise D.DataError("no seismic record found under %s" % index["raw_dir"])
    rows_tab = SM.table_rows(cfg, index["traces"], out, rep)
    chans, notes = D.choose_channel([t for t in index["traces"] if t["sta"]], cfg["data"]["channel"])
    rep.append(
        "  vertical component (data.channel = %s): %d stations; channels used: %s"
        % (cfg["data"]["channel"], len(chans), ", ".join(sorted({r["channel"] for v in chans.values() for r in v})))
    )
    rep += ["    " + n for n in notes]
    traces = [r for v in chans.values() for r in v]
    used, msum = SM.match(cfg, traces, out, rep, rows_tab)
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
    ex = set(D.resolve_names(cfg["data"]["exclude_stations"], used, "data.exclude_stations"))
    keep = {k: v for k, v in rows.items() if k not in ex}  # the geometry of the stations that will be used
    decisions, geom = auto.infer(cfg, {k: used[k] for k in keep}, [r for v in keep.values() for r in v])
    raw_rates = sorted({v[0]["fs"] for v in keep.values()})
    decisions.update(auto.from_velocity(cfg, geom, raw_rates))
    auto.reset(cfg, auto.REFERENCE_KEYS)  # measured again by the reference step on the new data
    auto.write(cfg, decisions, geom)
    cfg = auto.apply(cfg)
    rep.append("  parameters from the data (.lant/auto_parameters.json):")
    for k, v in decisions.items():
        if v["reason"] != "given in the configuration":
            rep.append(
                "    %s = %s  (%s)"
                % (k, json.dumps(v["value"]) if not isinstance(v["value"], str) else v["value"], v["reason"])
            )
    rows = keep
    SM.geometry_check({k: used[k] for k in keep})
    inst = instruments(cfg, {k: used[k] for k in keep}, keep)
    if inst["removed"]:
        rep.append(
            "  instrument responses removed (%s): %d stations, to ground velocity, pre-filter %s Hz"
            % (inst["path"], len(inst["ids"]), inst["pre_filt"])
        )
    else:
        rep.append(
            "  instrument responses NOT removed (no data.response): exact only if every station carries the "
            "same instrument model; the station table %s"
            % (
                "names one model for all stations (%s)" % next(iter(inst["models"]))
                if len(inst["models"]) == 1 and inst["named"] == inst["stations"]
                else "does not name the models"
            )
        )
    at = cfg.get("_array_type") or {}
    rep.append(
        "  array type %s%s"
        % (
            cfg["project"]["array_type"],
            " (given in the configuration)" if at.get("given") else ": " + at.get("text", ""),
        )
    )
    days = D.day_windows([r for v in rows.values() for r in v], cfg)
    nsd = sum(len({t["station"] for t in trs}) for _, trs in days.values())
    est_bytes = nsd * max(r[0]["fs"] for r in rows.values()) * 86400 * 4
    rep.append(
        "  %d stations used (%d excluded by data.exclude_stations), %d processing days, %d station-days"
        % (len(rows), len(ex), len(days), nsd)
    )
    from . import cc as CC

    coords = {n: (v["lon"], v["lat"]) for n, v in D.used_stations(cfg).items()}
    npair = CC.pair_indices(sorted(rows), coords, cfg["correlation"]["max_distance_km"])[0].size
    res = CC.resources(cfg, len(rows), npair, len(days), max(r[0]["fs"] for r in rows.values()))
    from .resources import shortages

    # Advance notice only: clean itself needs neither, memory may be freed before cc, and the clean report is the
    # evidence for deciding.  cc checks again and refuses; the archive's own disk check follows below.
    for problem in shortages(res, out, est_bytes if cfg["clean"]["write_archive"] else 0):
        import sys

        print("[lant] before cc: %s" % problem, file=sys.stderr)
        rep.append("  WARNING before cc: %s" % problem)
    from .config import DEFAULTS

    if cfg["project"]["array_type"] == "line" and list(cfg["data"]["velocity_range_kms"]) == list(
        DEFAULTS["data"]["velocity_range_kms"]
    ):
        import sys

        print(
            "[lant] line array uses the broad default velocity prior; review its scientific suitability "
            "and resource estimate before cc. No prior was changed.",
            file=sys.stderr,
        )
    rep.append("  " + CC.resources_text(cfg, res))
    with open(os.path.join(out, "processing_days.txt"), "w") as fh:
        fh.write("# label (data.day_label = %s)  UTC window\n" % cfg["data"].get("day_label", "utc_start"))
        for lab, (t0, _) in days.items():
            fh.write("%s  %s\n" % (lab, D.utc_window(t0)))
    if days:
        lab0, (t00, _) = next(iter(days.items()))
        rep.append(
            "  labels: data.day_label = %s, e.g. %s = %s (all in work/01_clean/processing_days.txt)"
            % (cfg["data"].get("day_label", "utc_start"), lab0, D.utc_window(t00))
        )
    archive = os.path.join(out, "archive")
    summary = dict(
        day_label=cfg["data"].get("day_label", "utc_start"),
        day_start_utc_hour=cfg["data"]["day_start_utc_hour"],
        files=index["n_files"],
        unreadable=len(index["unreadable"]),
        formats=dict(fmts),
        stations=len(rows),
        days=len(days),
        station_days=nsd,
        channel_notes=notes,
        names=msum,
        several_sensors=several,
        sampling_rates=sorted({r[0]["fs"] for r in rows.values()}),
        record_span_utc=[D.utc_date(min(t["start"] for t in traces)), D.utc_date(max(t["end"] for t in traces))],
        response={k: v for k, v in inst.items() if k != "ids"},
        array_type=dict(at, type=cfg["project"]["array_type"]),
        auto_parameters={k: v["value"] for k, v in decisions.items()},
    )
    if not cfg["clean"]["write_archive"]:
        idx = dict(
            source="raw",
            traces=[r for v in rows.values() for r in v if r["end"] - r["start"] >= cfg["data"]["min_trace_length_s"]],
        )
        json.dump(idx, open(os.path.join(out, "archive_index.json"), "w"))
        rep.append("  clean.write_archive is false: later steps read the raw files")
        rep.append(
            "  1-s periodic signal: NOT CHECKED (the records are read only when the archive is written); "
            "correlation.remove_1s_comb = auto removes nothing - set it to true if the loggers are known to "
            "produce it, or write the archive to have it tested"
        )
        summary.update(pulse_1s=dict(tested=0, stations=[], note="not checked: clean.write_archive is false"))
        open(os.path.join(out, "report.txt"), "w").write("\n".join(rep) + "\n")
        print("\n".join(rep))
        return summary
    free = shutil.disk_usage(out).free
    if free < cfg["clean"]["min_free_disk_factor"] * est_bytes:
        raise D.DataError(
            "the cleaned archive needs up to %.1f GB, %.1f GB are free on %s: free space, or set "
            "clean.write_archive: false" % (est_bytes / 1e9, free / 1e9, out)
        )
    opt = dict(
        format=cfg["data"]["format"],
        min_len=float(cfg["data"]["min_trace_length_s"]),
        archive=archive,
        encoding=cfg["clean"]["encoding"],
        codes=storage_codes(list(rows), rows),
        label_offset=D.label_offset_s(cfg["data"].get("day_label")),
        response=inst if inst["removed"] else None,
    )
    renamed = {n: c for n, c in opt["codes"].items() if "%s.%s" % c[:2] != n}
    if renamed:
        rep.append(
            "  %d names too long for miniSEED headers, archived under codes: %s"
            % (len(renamed), ", ".join("%s -> %s.%s" % (n, c[0], c[1]) for n, c in sorted(renamed.items())[:10]))
        )
    tasks = []
    for lab, (t0, trs) in days.items():
        bys = defaultdict(list)
        for r in trs:
            bys[r["station"]].append(r)
        for s, rr in bys.items():
            tasks.append((s, rr, t0, opt))
    stats, new_rows = [], []
    archive_files = {}
    from . import resources

    # One station-day at the highest raw rate, about 6 copies of 8 bytes per sample (the record, its float64 copy,
    # the filtered copy and the complex spectrum of the response removal): an engineering estimate.
    plan = resources.plan(cfg, len(tasks), 6 * 8 * 86400 * max(r[0]["fs"] for r in rows.values()), "station-days")
    with Pool(plan["parallel"], initializer=resources.limit_threads, initargs=(plan["threads_each"],)) as pool:
        for k, (name, rr, stat) in enumerate(
            track(pool_results(pool, _clean_station_day, tasks, chunksize=4), len(tasks), "station-days")
        ):
            stats.append(stat)
            if rr:
                new_rows.extend(rr)
                archive_files[(name, stat["day"])] = sorted({r["file"] for r in rr})
        resources.note_workers(w.pid for w in pool._pool)
    # the workers finish in any order: a fixed order makes two runs of the same data write the same files
    new_rows.sort(key=lambda r: (r["station"], r["start"], r["file"]))
    stats.sort(key=lambda s: (str(s.get("station")), str(s.get("day"))))
    reused = sum(1 for s in stats if s.pop("reused", 0))
    if reused:
        rep.append(
            "  %d of %d station-days taken from the archive of an earlier run (same raw files and settings); "
            "%d cleaned now" % (reused, len(tasks), len(tasks) - reused)
        )
    errors = [s for s in stats if "error" in s]
    source_map = os.path.join(out, "archive_sources.csv")
    with open(source_map, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["station", "day", "source_file", "archive_file", "operation"])
        for name, rr, t0, task_opt in tasks:
            day = D.day_label(t0, offset_s=task_opt.get("label_offset", 0.0))
            source_files = sorted({r["file"] for r in rr})
            outputs = archive_files.get((name, day), [])
            for source in source_files:
                w.writerow(
                    [
                        name,
                        day,
                        source,
                        ";".join(os.path.relpath(p, cfg["project"]["output_dir"]) for p in outputs),
                        ("merge_and_clean" if len(source_files) > 1 else "clean")
                        if outputs
                        else "no_output_see_report",
                    ]
                )
    rep.append("  raw files were not modified or deleted; raw-to-archive mapping: work/01_clean/archive_sources.csv")
    for e in errors:
        rep.append("  ERROR %s" % e["error"])
    stats = [s for s in stats if "error" not in s]
    keys = [
        "station",
        "day",
        "hours",
        "gaps",
        "traces_in",
        "shifted",
        "interpolated",
        "duplicates_or_identical_overlaps",
        "conflicting_overlap_s",
        "dropped_short",
        "nonfinite_samples",
        "polarity_reversed",
    ]
    with open(os.path.join(out, "station_days.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for s in sorted(stats, key=lambda s: (s["station"], s["day"])):
            w.writerow([s[k] for k in keys])
    tot = {k: sum(s[k] for s in stats) for k in keys[3:]}
    rep.append(
        "  archive %s: %d station-day files, %.1f station-hours; gaps %d; records shifted onto the grid %d, "
        "interpolated %d; duplicates/identical overlaps removed %d; conflicting overlaps left out %.1f s; "
        "pieces shorter than %.0f s dropped %d; non-finite samples treated as gaps %d; station-days of downward "
        "(cmpinc 180) sensors turned over %d"
        % (
            archive,
            len({(s["station"], s["day"]) for s in stats if s["hours"] > 0}),
            sum(s["hours"] for s in stats),
            tot["gaps"],
            tot["shifted"],
            tot["interpolated"],
            tot["duplicates_or_identical_overlaps"],
            tot["conflicting_overlap_s"],
            opt["min_len"],
            tot["dropped_short"],
            tot["nonfinite_samples"],
            tot["polarity_reversed"],
        )
    )
    pulse = PU.pulse_summary(stats, os.path.join(out, "pulse_1s.csv"))
    rep.append(
        "  1-s periodic signal (records folded modulo 1 s, F test at %g per station-day): %s"
        % (
            PU.PULSE_P,
            "not testable (no whole-number sampling rate)"
            if not pulse["tested"]
            else "%d of %d stations show it in more than half of their days%s; it is subtracted from those stations "
            "before the correlation (correlation.remove_1s_comb = auto)"
            % (
                len(pulse["stations"]),
                pulse["tested"],
                " (%s; up to %.2g %% of the variance)"
                % (
                    ", ".join(pulse["stations"][:10]) + (" ..." if len(pulse["stations"]) > 10 else ""),
                    100 * pulse["max_power"],
                )
                if pulse["stations"]
                else "",
            ),
        )
    )
    for row in new_rows:
        row["file"] = os.path.relpath(row["file"], cfg["project"]["output_dir"])
    json.dump(
        dict(source="archive", traces=new_rows, pulse_1s=pulse["stations"]),
        open(os.path.join(out, "archive_index.json"), "w"),
    )
    open(os.path.join(out, "report.txt"), "w").write("\n".join(rep) + "\n")
    coverage_figure(stats, sorted(days), os.path.join(out, "coverage.png"))
    summary.update(
        pulse_1s=pulse,
        archive_station_days=len(stats),
        station_hours=round(sum(s["hours"] for s in stats), 1),
        errors=len(errors),
        totals=tot,
        failed=[e["error"] for e in errors],
    )
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
