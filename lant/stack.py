"""Step 03: stacking, the station diagnosis, the stack that the later steps read.

stack.h5 holds, for every pair of the station list:
  stack        [npair, nt]          segment-weighted mean over all hours
  daily        [nblock, npair, nt]  sub-stacks of stack.block_hours hours (jackknife), labelled YYYYDDD[_hh]
  daily_count  [nblock, npair]      segments in each sub-stack
  nseg         [npair]              segments used
  dist_km, az_deg, lon1, lat1, lon2, lat2, pairs

Diagnosis (station_qc; areal arrays and lines alike; nothing is removed automatically: the candidates go to
03_stack/candidates.txt and are judged by the user, `lant confirm stack [--exclude-candidates]`):
  every sub-stack of a pair is compared with the stack of all the OTHER sub-stacks of that pair (leave-one-out, so
  that a sub-stack is not compared with itself): the correlation coefficient and the time shift (peak of their
  cross-correlation, parabolic interpolation between samples), over the whole correlation band, in the lag window
  |t| <= distance / (0.6 x slowest velocity).  Per station and sub-stack: the median over its pairs with the
  `partners` nearest stations (at least `min_partners` of them with data) of the correlation and of the shift
  counted towards that station, with the standard error of that median over the pairs.  A clock error or a moved /
  disturbed sensor shows as a shift or a low correlation common to all pairs of one station (idea: Stehly et al.
  2007; Sens-Schonfelder 2008).
  Candidates: outliers of the two distributions over all station sub-stacks, by the median and the scaled median
  absolute deviation (sigma = 1.4826 MAD; Leys et al. 2013) with Chauvenet's criterion for the number N of station
  sub-stacks: a value is a candidate when fewer than half a value that far out is expected among N normal values.
  The correlation enters as Fisher's z = atanh(r) (near-normal for correlation coefficients), relative to the
  array median of the same sub-stack (the noise level changes with time), low side only; the shift both sides, at
  least SHIFT_FLOOR of a sample, and beyond the same number of its own standard errors (a shift shared by the
  pairs of the station; sub-stacks too noisy to measure a shift scatter widely over the pairs and are not listed
  for their shift).  With fewer than MIN_SAMPLES station
  sub-stacks no threshold can be set and nothing is listed.
Removed: data.exclude_stations (clean step), stack.exclude_station_days / exclude_station_hours written by the
user (a day, a sub-stack YYYYDDD_hh, or a range "YYYYDDD-" / "YYYYDDD-YYYYDDD"), and the candidates when
stack.exclude_candidates is true (`lant confirm stack --exclude-candidates`).
"""
import glob
import os

import h5py
import numpy as np

from . import data as D
from . import paths
from .signal import geodesic_km


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
        raise RuntimeError("%d day files in %s cover another UTC window than the processing day of their label (were "
                           "data.day_label or data.day_start_utc_hour changed?): %s; run `lant run cc` again after "
                           "deleting them" % (len(wrong), d, " ".join(wrong[:10])))
    if not have:
        raise RuntimeError("no day files in %s: run `lant run cc` first" % d)
    missing = [lab for lab in sorted(days) if not os.path.exists(os.path.join(d, "day_%s.h5" % lab))]
    if missing:
        from . import records
        cc = records.read_summary(cfg, "cc") or {}
        few = {x.split(":")[0] for x in cc.get("skipped_days", [])}
        unfinished = [lab for lab in missing if lab not in few]
        if unfinished:
            raise RuntimeError("%d processing days have no day file although cc did not report them as days with fewer "
                               "than two stations (cc not finished or failed): %s; run `lant run cc` again"
                               % (len(unfinished), " ".join(unfinished[:10])))
        print("  note: %d processing days have fewer than two stations and no day file: %s"
              % (len(missing), " ".join(missing[:10])), flush=True)
    others = sorted(set(glob.glob(os.path.join(d, "day_*.h5"))) - set(files))
    if others:
        print("  note: %d day files in %s are not processing days of this configuration and are left out"
              % (len(others), d), flush=True)
    return have


def build_stack(cfg, out, exclude_hours=None, block_hours=None):
    """Stack the hourly files -> out.  exclude_hours: {(day, hour): set(stations)} removed before stacking."""
    stations = D.used_stations(cfg)
    names = sorted(stations)
    from .cc import pair_indices
    iu, ju = pair_indices(names, {k: (v["lon"], v["lat"]) for k, v in stations.items()},
                          cfg["correlation"]["max_distance_km"])
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
    dcount = np.zeros((nday, npair), np.int32)
    days = []
    tmp = out + ".tmp"
    ho = h5py.File(tmp, "w")
    dset = ho.create_dataset("daily", (nday, npair, nt), np.float32, chunks=(1, npair, nt))
    for fi, f in enumerate(files):
        daily = np.zeros((nblk, npair, nt), np.float32)      # the sub-stacks of this day only: written, then freed
        day = os.path.basename(f)[4:11]
        days += [day] if nblk == 1 else ["%s_%02d" % (day, b * bh) for b in range(nblk)]
        with h5py.File(f, "r") as h:
            from .cc import SETTINGS
            for k in SETTINGS:
                if h.attrs.get(k) != attrs.get(k):
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
                daily[hr // bh, gi] += x.astype(np.float32)
                dcount[fi * nblk + hr // bh, gi] += cnt[hr].astype(np.int32)
            nseg[gi] += cnt.sum(axis=0).astype(np.int64)
        dset[fi * nblk:(fi + 1) * nblk] = daily / np.maximum(dcount[fi * nblk:(fi + 1) * nblk], 1)[:, :, None]
        print("  stacked", os.path.basename(f), flush=True)
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


MIN_SAMPLES = 10          # station sub-stacks needed before outliers of their distribution can be named
SHIFT_FLOOR = 0.25        # samples: smallest time shift listed (empirical: the precision of the interpolated peak)


def _shift_corr(a, b, w0, w1, maxs):
    """Correlation coefficient and lag (samples, sub-sample by a parabola through the peak) of b against a: the
    lag s maximising corr(a[t], b[t + s]) for t in [w0, w1), |s| <= maxs.  A positive lag: b is late."""
    x = a[w0:w1]
    y = b[w0 - maxs:w1 + maxs]
    num = np.correlate(y, x, "valid")                         # num[s + maxs] = sum_t x[t] y[t + s]
    e = np.concatenate([[0.0], np.cumsum(y * y)])
    ey = e[len(x):] - e[:-len(x)]
    c = num / np.sqrt(np.maximum((x * x).sum() * ey, 1e-300))
    k = int(np.argmax(c))
    lag = float(k - maxs)
    if 0 < k < len(c) - 1:
        den = c[k - 1] - 2 * c[k] + c[k + 1]
        if den < 0:
            lag += 0.5 * (c[k - 1] - c[k + 1]) / den
    return float(c[k]), lag


EXPECTED_BY_CHANCE = 1      # candidates expected among normal values (0.5 per test); no more: listed, not gated


def outlier_limits(dev, shift, fs):
    """Thresholds of the station diagnosis from the two distributions -> dict, or None with too few samples."""
    from scipy.stats import norm
    dv, sh = dev[np.isfinite(dev)], shift[np.isfinite(shift)]
    n = int(dv.size)
    if n < MIN_SAMPLES:
        return None
    mad = lambda v: 1.4826 * float(np.median(np.abs(v - np.median(v))))
    k1, k2 = float(norm.isf(0.5 / n)), float(norm.isf(0.25 / n))
    return dict(n=n, k_corr=round(k1, 2), k_shift=round(k2, 2), corr_median=float(np.median(dv)), corr_sigma=mad(dv),
                corr_min=float(np.median(dv) - k1 * mad(dv)), shift_median=float(np.median(sh)), shift_sigma=mad(sh),
                shift_max=float(max(k2 * mad(sh), SHIFT_FLOOR / fs)))


def station_qc(cfg, stack_path, out_dir):
    """Station diagnosis from the stack of every hour -> (flags {(station, sub-stack label): reason}, summary)."""
    q = cfg["stack"]["station_qc"]
    paths.mk(out_dir)
    v_lo = float(cfg["data"]["velocity_range_kms"][0])
    h = h5py.File(stack_path, "r")
    pairs = [p.decode() for p in h["pairs"][:]]
    dist = h["dist_km"][:]
    nseg = h["nseg"][:].astype(np.float64)
    fs = float(h.attrs["fs"])
    nl = h["stack"].shape[1]
    mid = nl // 2
    days = [x.decode() for x in h["days"][:]]
    nday = len(days)
    stas = sorted({s for p in pairs for s in p.split("_")})
    si = {s: i for i, s in enumerate(stas)}
    dc = h["daily_count"][:].astype(np.float64)
    names = [p.split("_") for p in pairs]
    near = {}
    for k in np.where(nseg > 0)[0]:
        for s in names[k]:
            near.setdefault(s, []).append((dist[k], k))
    use_p = np.array(sorted({k for v in near.values() for _, k in sorted(v)[:int(q["partners"])]}), int)
    tsearch = 1.0 / np.sqrt(float(h.attrs["fmin"]) * float(h.attrs["fmax"]))
    maxs = max(2, int(round(tsearch * fs)))
    percorr = np.full((nday, len(use_p)), np.nan)
    pershift = np.full((nday, len(use_p)), np.nan)
    step = max(1, int(2e8 // (8 * max(nday, 1) * nl)))         # pairs read at a time (about 200 MB)
    for c0 in range(0, len(use_p), step):
        ids = use_p[c0:c0 + step]
        S = h["stack"][ids].astype(np.float64)
        Dd = h["daily"][:, ids].astype(np.float64)
        for qq, k in enumerate(ids):
            half = min(mid - maxs, int(np.ceil(dist[k] / (0.6 * v_lo) * fs)) + maxs)
            w0, w1 = mid - half, mid + half + 1
            cnt = dc[:, k]
            if not (cnt > 0).any():
                continue
            enough = cnt >= 0.5 * np.median(cnt[cnt > 0])
            for dd in np.where(enough & (cnt < nseg[k]))[0]:
                loo = (S[qq] * nseg[k] - Dd[dd, qq] * cnt[dd]) / (nseg[k] - cnt[dd])
                percorr[dd, c0 + qq], lag = _shift_corr(loo, Dd[dd, qq], w0, w1, maxs)
                pershift[dd, c0 + qq] = lag / fs
        if c0 // step % 20 == 0:
            print("  station diagnosis: pairs %d / %d" % (min(c0 + step, len(use_p)), len(use_p)), flush=True)
    h.close()
    sd_corr = np.full((len(stas), nday), np.nan)
    sd_shift = np.full((len(stas), nday), np.nan)
    sd_shift_se = np.full((len(stas), nday), np.nan)      # standard error of that median over the pairs
    present = np.zeros((len(stas), nday), bool)
    mine = {s: [] for s in stas}
    for col, k in enumerate(use_p):
        for pos, s in enumerate(names[k]):
            mine[s].append((col, 1 if pos == 0 else -1))
    for s, lst in mine.items():
        if not lst:
            continue
        cols, sgn = np.array([c for c, _ in lst]), np.array([g for _, g in lst], float)
        for dd in range(nday):
            c = percorr[dd, cols]
            m = np.isfinite(c)
            present[si[s], dd] = m.any()
            if m.sum() >= q["min_partners"]:
                sd_corr[si[s], dd] = np.median(c[m])
                v = pershift[dd, cols][m] * sgn[m]
                sd_shift[si[s], dd] = np.median(v)
                # a clock error or a moved sensor shifts every pair of the station alike; noise scatters them.
                # Standard error of the median: 1.2533 x (1.4826 MAD) / sqrt(n)
                sd_shift_se[si[s], dd] = 1.2533 * 1.4826 * np.median(np.abs(v - np.median(v))) / np.sqrt(m.sum())
    zc = np.arctanh(np.clip(sd_corr, -0.999999, 0.999999))      # Fisher's z: near-normal for correlation coefficients
    with np.errstate(all="ignore"):
        day_med = np.nanmedian(sd_corr, axis=0)
        dev = zc - np.nanmedian(zc, axis=0)[None, :]
    lim = outlier_limits(dev, sd_shift, fs)
    np.savez(os.path.join(out_dir, "station_qc.npz"), stations=np.array(stas), blocks=np.array(days), corr=sd_corr,
             shift=sd_shift, shift_se=sd_shift_se, block_median=day_med)
    flags = {}
    if lim:
        for i, s in enumerate(stas):
            for dd in range(nday):
                if not np.isfinite(dev[i, dd]):
                    continue
                why = []
                if dev[i, dd] < lim["corr_min"]:
                    why.append("corr %.2f (array %.2f)" % (sd_corr[i, dd], day_med[dd]))
                if abs(sd_shift[i, dd] - lim["shift_median"]) > lim["shift_max"] and \
                        abs(sd_shift[i, dd] - lim["shift_median"]) > lim["k_shift"] * sd_shift_se[i, dd]:
                    why.append("shift %+.4f s (+-%.4f s over its pairs)" % (sd_shift[i, dd], sd_shift_se[i, dd]))
                if why:
                    flags[(s, days[dd])] = "; ".join(why)
    judged = int(np.isfinite(sd_corr).sum())
    L = ["station diagnosis (lant stack) from %s" % stack_path,
         "  each sub-stack against the stack of the other sub-stacks of its pair; band %.4g-%.4g Hz; lags within "
         "distance / %.3g km/s; shift searched +-%.3g s" % (float(cfg["correlation"]["fmin_hz"]),
                                                             float(cfg["correlation"]["fmax_hz"]), 0.6 * v_lo, maxs / fs),
         "  per station: its pairs with the %d nearest stations (%d pairs in all), at least %d with data"
         % (q["partners"], len(use_p), q["min_partners"]),
         "  %d stations, %d sub-stacks: %d station sub-stacks judged, %d with data but fewer than %d pairs (not judged)"
         % (len(stas), nday, judged, int(present.sum()) - judged, q["min_partners"])]
    if lim:
        L += ["  correlation as Fisher's z = atanh(r), minus the array median of the sub-stack: median %+.3f, sigma "
              "%.3f (1.4826 MAD); candidate "
              "below %+.3f (%.2f sigma: Chauvenet's criterion for N = %d)"
              % (lim["corr_median"], lim["corr_sigma"], lim["corr_min"], lim["k_corr"], lim["n"]),
              "  time shift: median %+.5f s, sigma %.5f s; candidate beyond +-%.5f s of the median (%.2f sigma, at "
              "least %.2f sample) and beyond %.2f standard errors of the station's own median over its pairs (a "
              "shift common to the pairs of a station, not a scatter of poorly measured ones)"
              % (lim["shift_median"], lim["shift_sigma"], lim["shift_max"], lim["k_shift"], SHIFT_FLOOR, lim["k_shift"])]
    else:
        L.append("  fewer than %d station sub-stacks could be judged: no thresholds, nothing listed (NOT CHECKED)" % MIN_SAMPLES)
    L += ["", "  candidates: %d   (about %d expected by chance alone among normal values: Chauvenet's criterion lets "
          "0.5 through per test, correlation and time shift)" % (len(flags), EXPECTED_BY_CHANCE)]
    for (s, d), why in sorted(flags.items()):
        L.append("    %s %s  %s" % (s, d, why))
    open(os.path.join(out_dir, "station_qc.txt"), "w").write("\n".join(L) + "\n")
    _station_figure(out_dir, stas, days, dev, sd_shift, flags, lim)
    return flags, dict(evaluated=judged, not_judged=int(present.sum()) - judged, candidates=len(flags), limits=lim,
                       checked=bool(lim), pairs=int(len(use_p)))


def _station_figure(out_dir, stas, days, dev, sd_shift, flags, lim):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nday = len(days)
    fig, axs = plt.subplots(1, 2, figsize=(15, 6))
    x = np.tile(np.arange(nday), len(stas))
    axs[0].scatter(x, dev.ravel(), s=5, c="0.6")
    axs[1].scatter(x, sd_shift.ravel(), s=5, c="0.6")
    for (s, d) in flags:
        i, dd = stas.index(s), days.index(d)
        axs[0].plot(dd, dev[i, dd], "o", color="r", ms=5)
        axs[1].plot(dd, sd_shift[i, dd], "o", color="r", ms=5)
    if lim:
        axs[0].axhline(lim["corr_min"], color="r", lw=0.8)
        for sg in (-1, 1):
            axs[1].axhline(lim["shift_median"] + sg * lim["shift_max"], color="r", lw=0.8)
    step = max(1, nday // 40)
    for ax, yl in zip(axs, ("atanh(correlation) minus the array median of the sub-stack", "Time shift (s)")):
        ax.set_xticks(range(0, nday, step))
        ax.set_xticklabels([days[k][4:] for k in range(0, nday, step)], rotation=90, fontsize=7)
        ax.set_xlabel("Sub-stack (day of year_hour)")
        ax.set_ylabel(yl)
    axs[0].set_title("Station sub-stacks against the stack of the others (red: candidates)", fontsize=10)
    axs[1].set_title("Time shift of each station sub-stack (red lines: limits)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "station_qc.png"), dpi=140)
    plt.close(fig)


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
            sel = bad.get(dy, set()) | (bad.get(dy[:7], set()) if dy != dy[:7] else set())   # the block, or its whole day
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


def run(cfg):
    d = paths.mk(paths.P(cfg, paths.STACK))
    raw = os.path.join(d, "stack_all.h5")          # every hour; kept only until the product is written
    s = cfg["stack"]
    summary = {"all_hours": build_stack(cfg, raw)}
    clean = paths.stack_product(cfg)
    with h5py.File(raw, "r") as h:
        days = [x.decode() for x in h["days"][:]]
    cand_days = {}
    if s["station_qc"]["enabled"]:
        cand_days, summary["station_qc"] = station_qc(cfg, raw, os.path.join(d, "station_qc"))
    L = ["# candidates of the stack step (lant stack): nothing below is removed unless confirmed with",
         "#   lant confirm stack --exclude-candidates",
         "[station sub-stacks]  station label reason"]
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
        fh.write("# removed from the stack: the user's lists%s\n" % (" and the confirmed candidates" if s["exclude_candidates"] else ""))
        for k in sorted(set(ex_days)):
            fh.write("%s %s\n" % k)
        for (day, hr), sts in sorted(ex_hours.items()):
            for st in sorted(sts):
                fh.write("%s %s %02d\n" % (st, day, hr))
    if ex_hours:
        summary["product"] = build_stack(cfg, clean, exclude_hours=ex_hours)
        if ex_days:
            summary["removed_station_days"] = exclude_station_days(clean, sorted(set(ex_days)), clean)
        os.remove(raw)
    elif ex_days:
        summary["product"] = exclude_station_days(raw, sorted(set(ex_days)), clean)
        os.remove(raw)
    else:
        os.replace(raw, clean)                     # nothing removed: the stack of every hour is the product
        summary["product"] = "the stack of every hour (nothing removed)"
    summary["removed"] = dict(station_blocks=len(set(ex_days)), station_hours=sum(len(v) for v in ex_hours.values()))
    return summary
