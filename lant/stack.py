"""Step 03: stacking, the station diagnosis (lant.station_qc), the stack that the later steps read.

stack.h5 holds, for every pair of the station list:
  stack        [npair, nt]          segment-weighted mean over all hours
  daily        [nblock, npair, nt]  sub-stacks of stack.block_hours hours (jackknife), labelled YYYYDDD[_hh]
  daily_count  [nblock, npair]      segments in each sub-stack
  nseg         [npair]              segments used
  dist_km, az_deg, lon1, lat1, lon2, lat2, pairs

Diagnosis: lant.station_qc (nothing is removed automatically: the candidates go to work/03_stack/candidates.txt
and are judged by the user, `lant confirm stack [--exclude-candidates]`).
Removed: data.exclude_stations (clean step), stack.exclude_station_days / exclude_station_hours written by the
user (a day, a sub-stack YYYYDDD_hh, or a range "YYYYDDD-" / "YYYYDDD-YYYYDDD"), and the candidates when
stack.exclude_candidates is true (`lant confirm stack --exclude-candidates`).
"""

import glob
import hashlib
import json
import os

import h5py
import numpy as np
from numba import njit, prange

from . import data as D
from . import paths
from .signal import geodesic_km
from .station_qc import EXPECTED_BY_CHANCE, station_qc

ALL_HOURS_VERSION = 1  # raise when build_stack changes what it writes: an all-hours stack kept from before is redone


@njit(cache=True, parallel=True)
def _add_hour(hourly, cnt, gi, total, daily_block):
    """One hour into the sums, pair by pair in parallel; per sample the operations of the former NumPy form:
    x = float64(hourly) * count; total += x; daily_block += float32(x)."""
    for p in prange(gi.size):
        g = gi[p]
        c = cnt[p]
        for t in range(hourly.shape[1]):
            x = np.float64(hourly[p, t]) * c
            total[g, t] += x
            daily_block[g, t] += np.float32(x)


def day_files(cfg):
    """The day files of the processing days of this configuration (not every day_*.h5 that happens to be in the
    directory: files of days outside start_date / end_date, or left from an earlier run, are not stacked)."""
    stations = D.used_stations(cfg)
    rows, _ = D.station_rows(D.processing_rows(cfg), stations)
    days = D.day_windows([r for v in rows.values() for r in v], cfg)
    d = paths.P(cfg, paths.CCF)
    files = [os.path.join(d, "day_%s.h5" % lab) for lab in sorted(days)]
    have = [f for f in files if os.path.exists(f)]
    t0_of = {os.path.join(d, "day_%s.h5" % lab): t0 for lab, (t0, _) in days.items()}
    wrong = []
    for f in have:
        with h5py.File(f, "r") as h:
            if float(h.attrs.get("t0_utc", np.nan)) != t0_of[f]:
                wrong.append(os.path.basename(f))
    if wrong:
        raise RuntimeError(
            "%d day files in %s cover another UTC window than the processing day of their label (were "
            "data.day_label or data.day_start_utc_hour changed?): %s; run `lant run cc` again after "
            "deleting them" % (len(wrong), d, " ".join(wrong[:10]))
        )
    if not have:
        raise RuntimeError("no day files in %s: run `lant run cc` first" % d)
    missing = [lab for lab in sorted(days) if not os.path.exists(os.path.join(d, "day_%s.h5" % lab))]
    if missing:
        from . import records

        cc = records.read_summary(cfg, "cc") or {}
        few = {x.split(":")[0] for x in cc.get("skipped_days", [])}
        unfinished = [lab for lab in missing if lab not in few]
        if unfinished:
            raise RuntimeError(
                "%d processing days have no day file although cc did not report them as days with fewer "
                "than two stations (cc not finished or failed): %s; run `lant run cc` again"
                % (len(unfinished), " ".join(unfinished[:10]))
            )
        print(
            "  note: %d processing days have fewer than two stations and no day file: %s"
            % (len(missing), " ".join(missing[:10])),
            flush=True,
        )
    others = sorted(set(glob.glob(os.path.join(d, "day_*.h5"))) - set(files))
    if others:
        print(
            "  note: %d day files in %s are not processing days of this configuration and are left out"
            % (len(others), d),
            flush=True,
        )
    return have


def build_stack(cfg, out, exclude_hours=None, block_hours=None):
    """Stack the hourly files -> out.  exclude_hours: {(day, hour): set(stations)} removed before stacking."""
    stations = D.used_stations(cfg)
    names = sorted(stations)
    from .cc import pair_indices

    iu, ju = pair_indices(
        names, {k: (v["lon"], v["lat"]) for k, v in stations.items()}, cfg["correlation"]["max_distance_km"]
    )
    pair_names = ["%s_%s" % (names[i], names[j]) for i, j in zip(iu, ju)]
    index = {p: k for k, p in enumerate(pair_names)}
    npair = len(pair_names)
    files = day_files(cfg)
    with h5py.File(files[0], "r") as h:
        nt = h["hourly"].shape[2]
        attrs = dict(h.attrs)
    bh = block_hours or cfg["stack"]["block_hours"]
    nblk = 24 // bh
    nday = len(files) * nblk
    drop = exclude_hours or {}
    import numba

    from . import resources

    threads = resources.plan(cfg, 1, None, "stack: pairs of an hour, in parallel")["threads_each"]
    previous = numba.get_num_threads()
    numba.set_num_threads(max(1, min(threads, numba.config.NUMBA_NUM_THREADS)))
    total = np.zeros((npair, nt), np.float64)
    nseg = np.zeros(npair, np.int64)
    dcount = np.zeros((nday, npair), np.int32)
    days = []
    tmp = out + ".tmp"
    ho = h5py.File(tmp, "w")
    dset = ho.create_dataset("daily", (nday, npair, nt), np.float32, chunks=(1, npair, nt))
    for fi, f in enumerate(files):
        daily = np.zeros((nblk, npair, nt), np.float32)  # the sub-stacks of this day only: written, then freed
        day = os.path.basename(f)[4:11]
        days += [day] if nblk == 1 else ["%s_%02d" % (day, b * bh) for b in range(nblk)]
        with h5py.File(f, "r") as h:
            from .cc import SETTINGS

            for k in SETTINGS:
                if h.attrs.get(k) != attrs.get(k):
                    raise RuntimeError(
                        "%s was computed with %s = %s, the first day file with %s" % (f, k, h.attrs[k], attrs[k])
                    )
            pn = [p.decode() for p in h["pairs"][:]]
            gi = np.array([index.get(p, -1) for p in pn])
            keep = gi >= 0
            gi = gi[keep]
            cnt = h["count"][:][:, keep].astype(np.float64)
            kept_pairs = [q for q, k_ in zip(pn, keep) if k_]
            for hr in range(24):
                bad = drop.get((day, hr))
                if bad:
                    cnt[hr, [any(s_ in bad for s_ in q.split("_")) for q in kept_pairs]] = 0
            for hr in range(24):
                if not cnt[hr].any():
                    continue
                hour = np.ascontiguousarray(h["hourly"][hr][keep])
                _add_hour(hour, np.ascontiguousarray(cnt[hr]), gi.astype(np.int64), total, daily[hr // bh])
                dcount[fi * nblk + hr // bh, gi] += cnt[hr].astype(np.int32)
            nseg[gi] += cnt.sum(axis=0).astype(np.int64)
        dset[fi * nblk : (fi + 1) * nblk] = daily / np.maximum(dcount[fi * nblk : (fi + 1) * nblk], 1)[:, :, None]
        print("  stacked", os.path.basename(f), flush=True)
    numba.set_num_threads(previous)
    total /= np.maximum(nseg, 1)[:, None]
    lon1 = np.array([stations[names[i]]["lon"] for i in iu])
    lat1 = np.array([stations[names[i]]["lat"] for i in iu])
    lon2 = np.array([stations[names[j]]["lon"] for j in ju])
    lat2 = np.array([stations[names[j]]["lat"] for j in ju])
    dist, az = geodesic_km(lon1, lat1, lon2, lat2)
    with ho as h:
        h.create_dataset("stack", data=total.astype(np.float32))
        h.create_dataset("daily_count", data=dcount)
        h.create_dataset("nseg", data=nseg)
        h.create_dataset("pairs", data=paths.ascii_names(pair_names))
        h.create_dataset("days", data=paths.ascii_names(days))
        for k, v in dict(lon1=lon1, lat1=lat1, lon2=lon2, lat2=lat2, dist_km=dist, az_deg=az).items():
            h.create_dataset(k, data=v)
        h.attrs.update({k: v for k, v in attrs.items() if k != "t0_utc"})
        h.attrs["block_hours"] = bh
    os.replace(tmp, out)
    return dict(pairs=npair, pairs_with_data=int((nseg > 0).sum()), blocks=nday, segments=int(nseg.sum()))


def exclude_station_days(stack_path, bad_list, out):
    """Rebuild a stack without the listed sub-stacks: an entry (station, label) removes the sub-stack with that
    label (YYYYDDD or YYYYDDD_hh, as written by the station-day QC), and an entry with a whole day YYYYDDD
    removes every sub-stack of that day."""
    bad = {}
    for s, d in bad_list:
        bad.setdefault(d, set()).add(s)
    with h5py.File(stack_path, "r") as h, h5py.File(out + ".tmp", "w") as o:
        pairs = [p.decode() for p in h["pairs"][:]]
        days = [d.decode() for d in h["days"][:]]
        cnt = h["daily_count"][:].astype(np.int64)
        members = [set(p.split("_")) for p in pairs]
        dropped = 0
        for k, dy in enumerate(days):
            sel = bad.get(dy, set()) | (
                bad.get(dy[:7], set()) if dy != dy[:7] else set()
            )  # the block, or its whole day
            if not sel:
                continue
            for i, m in enumerate(members):
                if cnt[k, i] and (m & sel):
                    dropped += cnt[k, i]
                    cnt[k, i] = 0
        npair, nlag = h["stack"].shape
        acc = np.zeros((npair, nlag), np.float64)
        dd = o.create_dataset("daily", (len(days), npair, nlag), np.float32, chunks=(1, npair, nlag))
        for k in range(len(days)):
            x = h["daily"][k]
            x[cnt[k] == 0] = 0
            dd[k] = x
            acc += x.astype(np.float64) * cnt[k][:, None]
        tot = cnt.sum(0)
        o["stack"] = (acc / np.maximum(tot, 1)[:, None]).astype(np.float32)
        o["daily_count"] = cnt.astype(np.int32)
        o["nseg"] = tot
        for k in ("days", "az_deg", "dist_km", "lat1", "lat2", "lon1", "lon2", "pairs"):
            o[k] = h[k][:]
        for k, v in h.attrs.items():
            o.attrs[k] = v
        before = int(h["daily_count"][:].astype(np.int64).sum())
    os.replace(out + ".tmp", out)
    return dict(
        dropped_segments=int(dropped),
        segments_before=before,
        dropped_percent=round(100.0 * dropped / max(before, 1), 3),
        pairs_with_data=int((tot > 0).sum()),
    )


def expand_days(entries, days):
    """[station, label] entries with ranges "YYYYDDD-" / "YYYYDDD-YYYYDDD" -> [(station, label)] over `days`."""
    out = []
    for e in entries:
        s, d = str(e[0]), str(e[1])
        if "-" in d:
            a, b = d.split("-", 1)
            out += [(s, x) for x in days if x[:7] >= a and (not b or x[:7] <= b)]
        else:
            out.append((s, d))
    return out


def all_hours_key(cfg):
    """Fingerprint of everything the stack of every hour depends on: the day files (path, size, times), the stations
    and their coordinates, the pair distance limit, the sub-stack length and ALL_HOURS_VERSION."""
    from .record_index import file_signature

    stations = D.used_stations(cfg)
    what = [
        ALL_HOURS_VERSION,
        [[f, file_signature(f)] for f in day_files(cfg)],
        [[n, stations[n]["lon"], stations[n]["lat"]] for n in sorted(stations)],
        cfg["correlation"]["max_distance_km"],
        cfg["stack"]["block_hours"],
    ]
    return hashlib.sha1(json.dumps(what, default=str).encode()).hexdigest()


def run(cfg):
    """The stack of every hour, the station diagnosis on it, then the product without what is removed.  When the
    product of the last run is the stack of every hour and its inputs are unchanged (all_hours.json: the key and the
    file's signature), it is not built again: the usual case is the run after `lant confirm stack
    --exclude-candidates`, which only removes sub-stacks."""
    from .record_index import file_signature

    d = paths.mk(paths.P(cfg, paths.STACK))
    raw = os.path.join(d, "stack_all.h5")  # every hour; kept only until the product is written
    s = cfg["stack"]
    clean = paths.stack_product(cfg)
    note = os.path.join(d, "all_hours.json")
    key = all_hours_key(cfg)
    reused = None
    try:
        kept = json.load(open(note))
        if kept["key"] == key and os.path.exists(clean) and kept["signature"] == file_signature(clean):
            reused = kept["summary"]
    except (OSError, ValueError, KeyError):
        pass
    if reused is not None:
        raw = clean  # the product of the last run is the stack of every hour
        summary = {"all_hours": dict(reused, reused_from_last_run=True)}
        print("  the stack of every hour of the last run is reused (inputs unchanged)", flush=True)
    else:
        summary = {"all_hours": build_stack(cfg, raw)}
    if os.path.exists(note):
        os.remove(note)
    with h5py.File(raw, "r") as h:
        days = [x.decode() for x in h["days"][:]]
    cand_days = {}
    if s["station_qc"]["enabled"]:
        cand_days, summary["station_qc"] = station_qc(cfg, raw, os.path.join(d, "station_qc"))
    L = [
        "# candidates of the stack step (lant stack): nothing below is removed unless confirmed with",
        "#   lant confirm stack --exclude-candidates",
        "[station sub-stacks]  station label reason",
    ]
    L += ["%s %s  %s" % (k[0], k[1], v) for k, v in sorted(cand_days.items())]
    open(os.path.join(d, "candidates.txt"), "w").write("\n".join(L) + "\n")
    summary["candidates"] = dict(station_blocks=len(cand_days), expected_by_chance=EXPECTED_BY_CHANCE)
    # what is removed: the user's lists, and the candidates when confirmed so
    ex_days = expand_days(s["exclude_station_days"] or [], days)
    ex_hours = {}
    for e in s["exclude_station_hours"] or []:
        ex_hours.setdefault((str(e[1]), int(e[2])), set()).add(str(e[0]))
    if s["exclude_candidates"]:
        ex_days += list(cand_days)
    with open(os.path.join(d, "removed.txt"), "w") as fh:
        fh.write(
            "# removed from the stack: the user's lists%s\n"
            % (" and the confirmed candidates" if s["exclude_candidates"] else "")
        )
        for k in sorted(set(ex_days)):
            fh.write("%s %s\n" % k)
        for (day, hr), sts in sorted(ex_hours.items()):
            for st in sorted(sts):
                fh.write("%s %s %02d\n" % (st, day, hr))
    if ex_hours:
        summary["product"] = build_stack(cfg, clean, exclude_hours=ex_hours)
        if ex_days:
            summary["removed_station_days"] = exclude_station_days(clean, sorted(set(ex_days)), clean)
        if raw != clean:
            os.remove(raw)
    elif ex_days:
        summary["product"] = exclude_station_days(raw, sorted(set(ex_days)), clean)
        if raw != clean:
            os.remove(raw)
    else:
        if raw != clean:
            os.replace(raw, clean)  # nothing removed: the stack of every hour is the product
        summary["product"] = "the stack of every hour (nothing removed)"
        all_hours = {k: v for k, v in summary["all_hours"].items() if k != "reused_from_last_run"}
        with open(note + ".tmp", "w") as fh:
            json.dump(dict(key=key, signature=file_signature(clean), summary=all_hours), fh)
        os.replace(note + ".tmp", note)
    summary["removed"] = dict(station_blocks=len(set(ex_days)), station_hours=sum(len(v) for v in ex_hours.values()))
    return summary
