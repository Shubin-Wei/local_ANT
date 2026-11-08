"""Step 03: stacking, station-day / station-hour screening, the cleaned stack.

stack.h5 holds, for every pair of the station list:
  stack        [npair, nt]          segment-weighted mean over all hours
  daily        [nblock, npair, nt]  sub-stacks of stack.block_hours hours (jackknife), labelled YYYYDDD[_hh]
  daily_count  [nblock, npair]      segments in each sub-stack
  nseg         [npair]              segments used
  dist_km, az_deg, lon1, lat1, lon2, lat2, pairs

Screening (rules of the method; the thresholds are project values, see docs/parameters.md):
  areal arrays, station-days (station_day_qc): 1-4 Hz daily sub-stacks of the pairs 2-10 km long against the
      pair's total stack in the surface-wave window (1.8-3.5 km/s, +-0.5 s), best lag within +-0.5 s; per
      station and day the median over >= 10 partners of the correlation and of the common shift.  Listed when the
      correlation is more than corr_drop below that day's array median or |shift| >= shift_s.  Stations with a
      field note (moved, dug out) are cut from the first day after which every evaluable day is more than
      field_drop below the day median.  (Idea: Stehly et al. 2007; Sens-Schonfelder 2008.)
  lines, station-hours (station_hour_qc): zero-lag correlation of each hourly CCF (1-4 Hz, +-10 s, pairs up to
      6 km) with the pair's total stack, median over >= 3 partners; hours below the threshold, or with too few
      partners to be checked, are removed.
"""
import glob
import os

import h5py
import numpy as np
from scipy.signal import butter, sosfiltfilt

from . import data as D
from . import paths
from .signal import best_lag, geodesic_km


def day_files(cfg):
    files = sorted(glob.glob(os.path.join(paths.P(cfg, paths.CCF), "day_*.h5")))
    if not files:
        raise RuntimeError("no day files in %s: run `lant run cc` first" % paths.P(cfg, paths.CCF))
    return files


def build_stack(cfg, out, exclude_hours=None, block_hours=None):
    """Stack the hourly files -> out.  exclude_hours: {(day, hour): set(stations)} removed before stacking."""
    stations = D.used_stations(cfg)
    names = sorted(stations)
    iu, ju = np.triu_indices(len(names), 1)
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
    total = np.zeros((npair, nt), np.float64)
    nseg = np.zeros(npair, np.int64)
    daily = np.zeros((nday, npair, nt), np.float32)
    dcount = np.zeros((nday, npair), np.int32)
    days = []
    for fi, f in enumerate(files):
        day = os.path.basename(f)[4:11]
        days += [day] if nblk == 1 else ["%s_%02d" % (day, b * bh) for b in range(nblk)]
        with h5py.File(f, "r") as h:
            for k in ("fs", "maxlag", "seg_len", "seg_step", "fmin", "fmax"):
                if h.attrs[k] != attrs[k]:
                    raise RuntimeError("%s was computed with %s = %s, the first day file with %s" % (f, k, h.attrs[k], attrs[k]))
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
                x = h["hourly"][hr][keep].astype(np.float64) * cnt[hr][:, None]
                total[gi] += x
                daily[fi * nblk + hr // bh, gi] += x.astype(np.float32)
                dcount[fi * nblk + hr // bh, gi] += cnt[hr].astype(np.int32)
            nseg[gi] += cnt.sum(axis=0).astype(np.int64)
        print("  stacked", os.path.basename(f), flush=True)
    total /= np.maximum(nseg, 1)[:, None]
    daily /= np.maximum(dcount, 1)[:, :, None]
    lon1 = np.array([stations[names[i]]["lon"] for i in iu])
    lat1 = np.array([stations[names[i]]["lat"] for i in iu])
    lon2 = np.array([stations[names[j]]["lon"] for j in ju])
    lat2 = np.array([stations[names[j]]["lat"] for j in ju])
    dist, az = geodesic_km(lon1, lat1, lon2, lat2)
    tmp = out + ".tmp"
    with h5py.File(tmp, "w") as h:
        h.create_dataset("stack", data=total.astype(np.float32))
        h.create_dataset("daily", data=daily, chunks=(1, npair, nt))
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


def station_day_qc(cfg, stack_path, out_dir):
    """Station-day list from the preliminary stack -> (flags {(station, day): reason}, report text)."""
    q = cfg["stack"]["station_day_qc"]
    paths.mk(out_dir)
    h = h5py.File(stack_path, "r")
    pairs = [p.decode() for p in h["pairs"][:]]
    dist = h["dist_km"][:]
    nseg = h["nseg"][:]
    fs = float(h.attrs["fs"])
    nl = h["stack"].shape[1]
    lag = (np.arange(nl) - nl // 2) / fs
    days = [x.decode() for x in h["days"][:]]
    nday = len(days)
    stas = sorted({s for p in pairs for s in p.split("_")})
    lo, hi = q["band_hz"]
    sos = butter(4, [lo / (fs / 2), hi / (fs / 2)], "band", output="sos")
    dc = h["daily_count"][:]
    use_p = np.where((nseg > 0) & (dist >= q["dist_km"][0]) & (dist <= q["dist_km"][1]))[0]
    percorr = np.full((nday, len(use_p)), np.nan)
    pershift = np.full((nday, len(use_p)), np.nan)
    maxs = int(0.5 * fs)
    for c0 in range(0, len(use_p), 2000):
        ids = use_p[c0:c0 + 2000]
        S = sosfiltfilt(sos, h["stack"][ids].astype(np.float64), axis=-1)
        Dd = sosfiltfilt(sos, h["daily"][:, ids].astype(np.float64), axis=-1)
        for qq, k in enumerate(ids):
            w = (np.abs(lag) >= dist[k] / 3.5 - 0.5) & (np.abs(lag) <= dist[k] / 1.8 + 0.5)
            for dd in range(nday):
                if dc[dd, k] < q["min_segments"]:
                    continue
                l_, c_ = best_lag(S[qq][w], Dd[dd, qq][w], maxs)
                percorr[dd, c0 + qq] = c_
                pershift[dd, c0 + qq] = l_ / fs
        print("  station-day QC: pairs %d / %d" % (min(c0 + 2000, len(use_p)), len(use_p)), flush=True)
    sd_corr = np.full((len(stas), nday), np.nan)
    sd_shift = np.full((len(stas), nday), np.nan)
    names = [pairs[k].split("_") for k in use_p]
    for i, s in enumerate(stas):
        mine = np.array([k for k in range(len(use_p)) if s in names[k]], int)
        if not mine.size:
            continue
        sgn = np.array([1 if names[k][0] == s else -1 for k in mine])
        for dd in range(nday):
            c = percorr[dd, mine]
            m = np.isfinite(c)
            if m.sum() >= q["min_partners"]:
                sd_corr[i, dd] = np.median(c[m])
                sd_shift[i, dd] = np.median(pershift[dd, mine][m] * sgn[m])
    np.savez(os.path.join(out_dir, "station_day_ccf.npz"), stations=np.array(stas), days=np.array(days),
             corr=sd_corr, shift=sd_shift)
    with np.errstate(all="ignore"):
        day_med = np.nanmedian(sd_corr, axis=0)
    low_day = day_med < np.nanmedian(day_med) - q["field_drop"]
    flags = {}
    for i, s in enumerate(stas):
        for dd in range(nday):
            c, t_ = sd_corr[i, dd], sd_shift[i, dd]
            if not np.isfinite(c):
                continue
            why = []
            if c < day_med[dd] - q["corr_drop"]:
                why.append("corr %.2f (day %.2f)" % (c, day_med[dd]))
            if abs(t_) >= q["shift_s"]:
                why.append("shift %+.3f s" % t_)
            if why:
                flags[(s, days[dd])] = "CCF: " + "; ".join(why)
    notes = []
    for s in cfg["data"]["field_note_stations"] or []:
        if s not in stas:
            notes.append("    %s: not in the stack" % s)
            continue
        i = stas.index(s)
        low = np.isfinite(sd_corr[i]) & (sd_corr[i] < day_med - q["field_drop"])
        ev = np.isfinite(sd_corr[i])
        mine_pairs = [k for k, p in enumerate(pairs) if s in p.split("_")]
        has = np.array([dc[dd][mine_pairs].any() for dd in range(nday)])
        start = None
        for dd in range(nday):
            if low[dd] and not low_day[dd] and np.all(low[dd:][ev[dd:]]):
                start = dd
                break
        if start is not None:
            for dd in range(start, nday):
                if ev[dd]:
                    flags.setdefault((s, days[dd]), "field note: persistently low from day %s on, corr %.2f (day %.2f)"
                                     % (days[start], sd_corr[i, dd], day_med[dd]))
                elif has[dd]:
                    flags.setdefault((s, days[dd]), "field note: not evaluable, inside the cut from day %s" % days[start])
        notes.append("    %s: %s" % (s, "cut from day %s" % days[start] if start is not None else "no persistent anomaly, nothing cut"))
    h.close()
    L = ["station-day QC (lant stack) from %s" % stack_path,
         "  band %.1f-%.1f Hz, pairs %.1f-%.1f km, partners >= %d, sub-stacks with >= %d segments"
         % (lo, hi, q["dist_km"][0], q["dist_km"][1], q["min_partners"], q["min_segments"]),
         "  rule: correlation < day median - %.2f, or |common shift| >= %.2f s" % (q["corr_drop"], q["shift_s"]),
         "  field-note stations: %s" % (", ".join(cfg["data"]["field_note_stations"] or []) or "none"), *notes,
         "  %d stations, %d days, %d station-days evaluated" % (len(stas), nday, np.isfinite(sd_corr).sum()),
         "  day medians: " + " ".join("%s:%.2f" % (d, v) for d, v in zip(days, day_med)),
         "", "  listed station-days: %d" % len(flags)]
    for (s, d), why in sorted(flags.items()):
        L.append("    %s %s  %s" % (s, d, why))
    txt = "\n".join(L) + "\n"
    open(os.path.join(out_dir, "station_day_qc.txt"), "w").write(txt)
    _station_day_figure(out_dir, stas, days, sd_corr, sd_shift, day_med, flags, q)
    return flags, dict(evaluated=int(np.isfinite(sd_corr).sum()), listed=len(flags), day_median=dict(zip(days, day_med)))


def _station_day_figure(out_dir, stas, days, sd_corr, sd_shift, day_med, flags, q):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nday = len(days)
    fig, axs = plt.subplots(1, 2, figsize=(15, 6))
    dev = sd_corr - day_med[None, :]
    axs[0].scatter(np.tile(np.arange(nday), len(stas)), dev.ravel(), s=5, c="0.6")
    axs[1].scatter(np.tile(np.arange(nday), len(stas)), sd_shift.ravel(), s=5, c="0.6")
    for (s, d), why in flags.items():
        i, dd = stas.index(s), days.index(d)
        col = "C0" if why.startswith("field") else "r"
        axs[0].plot(dd, dev[i, dd], "o", color=col, ms=5)
        axs[1].plot(dd, sd_shift[i, dd], "o", color=col, ms=5)
    axs[0].axhline(-q["corr_drop"], color="r", lw=0.8)
    for y in (-q["shift_s"], q["shift_s"]):
        axs[1].axhline(y, color="r", lw=0.8)
    for ax, yl in zip(axs, ("Correlation minus the day's array median", "Common time shift (s)")):
        ax.set_xticks(range(nday))
        ax.set_xticklabels([d[-3:] for d in days], rotation=90, fontsize=7)
        ax.set_xlabel("Day of year")
        ax.set_ylabel(yl)
    axs[0].set_title("Station-days against the station's total stack (red: listed; blue: field-note rule)", fontsize=10)
    axs[1].set_title("Common time shift of each station-day (red lines: threshold)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "station_day_qc.png"), dpi=140)
    plt.close(fig)


def exclude_station_days(stack_path, bad_list, out):
    """Rebuild a stack without the listed (station, day) sub-stacks (block labels starting with the day)."""
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
            sel = bad.get(dy[:7])
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
    return dict(dropped_segments=int(dropped), segments_before=before,
                dropped_percent=round(100.0 * dropped / max(before, 1), 3), pairs_with_data=int((tot > 0).sum()))


def station_hour_qc(cfg, stack_path, out_dir):
    """Lines: score of every station-hour -> {(day, hour): set(stations)} to drop, table written to out_dir."""
    q = cfg["stack"]["station_hour_qc"]
    paths.mk(out_dir)
    with h5py.File(stack_path, "r") as h:
        fs = float(h.attrs["fs"])
        L = int(round(float(h.attrs["maxlag"]) * fs))
        k = int(q["maxlag_s"] * fs)
        pairs = [p.decode() for p in h["pairs"][:]]
        dist = h["dist_km"][:]
        ref = h["stack"][:, L - k:L + k + 1].astype(np.float64)
        nseg = h["nseg"][:]
    sos = butter(4, list(q["band_hz"]), btype="band", fs=fs, output="sos")
    ref = sosfiltfilt(sos, ref, axis=1)
    idx = {p: i for i, p in enumerate(pairs)}
    rows = []
    for f in day_files(cfg):
        day = os.path.basename(f)[4:11]
        with h5py.File(f, "r") as h:
            hp = [p.decode() for p in h["pairs"][:]]
            cnt = h["count"][:]
            gi = np.array([idx.get(p, -1) for p in hp])
            for hr in range(24):
                live = np.where((cnt[hr] > 0) & (gi >= 0))[0]
                if live.size == 0:
                    continue
                live = live[(dist[gi[live]] <= q["max_dist_km"]) & (nseg[gi[live]] > 0)]
                if live.size == 0:
                    continue
                x = sosfiltfilt(sos, h["hourly"][hr][live][:, L - k:L + k + 1].astype(np.float64), axis=1)
                r = ref[gi[live]]
                cc = (x * r).sum(1) / np.sqrt((x * x).sum(1) * (r * r).sum(1) + 1e-30)
                sta = {}
                for c, j in zip(cc, live):
                    s1, s2 = hp[j].split("_")
                    sta.setdefault(s1, []).append(c)
                    sta.setdefault(s2, []).append(c)
                for s, v in sta.items():
                    sc = float(np.median(v)) if len(v) >= q["min_partners"] else np.nan
                    rows.append((s, day, hr, len(v), sc))
    with open(os.path.join(out_dir, "station_hour_qc.txt"), "w") as fh:
        fh.write("# station day hour partners score  (score nan: fewer than %d partners, not checkable)\n" % q["min_partners"])
        for r in rows:
            fh.write("%s %s %02d %3d %+.3f\n" % r)
    drop = {}
    bad = [r for r in rows if not r[4] >= q["threshold"]]
    for s, day, hr, _, _ in bad:
        drop.setdefault((day, hr), set()).add(s)
    with open(os.path.join(out_dir, "exclude_station_hours.txt"), "w") as fh:
        fh.write("# station day hour  (score < %.2f or not checkable)\n" % q["threshold"])
        for r in bad:
            fh.write("%s %s %02d\n" % r[:3])
    sc = np.array([r[4] for r in rows])
    return drop, dict(station_hours=len(rows), excluded=len(bad), not_checkable=int((~np.isfinite(sc)).sum()),
                      score_percentiles={p: float(np.nanpercentile(sc, p)) for p in (5, 25, 50)} if np.isfinite(sc).any() else {})


def run(cfg):
    d = paths.mk(paths.P(cfg, paths.STACK))
    raw = os.path.join(d, "stack_raw.h5")
    s = cfg["stack"]
    summary = {"raw": build_stack(cfg, raw)}
    clean = os.path.join(d, "stack_clean.h5")
    extra = [tuple(x) for x in (s["exclude_station_days"] or [])]
    if s["station_hour_qc"]["enabled"]:
        drop, summary["station_hour_qc"] = station_hour_qc(cfg, raw, os.path.join(d, "station_hour_qc"))
        summary["clean"] = build_stack(cfg, clean, exclude_hours=drop)
        if extra:
            summary["extra_station_days"] = exclude_station_days(clean, extra, clean)
    elif s["station_day_qc"]["enabled"]:
        flags, summary["station_day_qc"] = station_day_qc(cfg, raw, os.path.join(d, "station_day_qc"))
        lst = sorted(set(flags) | set(extra))
        with open(os.path.join(d, "exclude_station_days.txt"), "w") as fh:
            fh.write("# station  day  reason (generated by lant from the rule in stack.station_day_qc; "
                     "entries of stack.exclude_station_days marked 'config')\n")
            for k in lst:
                fh.write("%s %s  %s\n" % (k[0], k[1], flags.get(k, "config")))
        summary["clean"] = exclude_station_days(raw, lst, clean)
    else:
        summary["clean"] = exclude_station_days(raw, extra, clean)
    return summary
