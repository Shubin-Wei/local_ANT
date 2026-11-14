"""Station diagnosis of the stack step (lant.stack): which station sub-stacks differ from the rest.

Areal arrays and lines alike; nothing is removed automatically: the candidates go to
work/03_stack/candidates.txt and are judged by the user, `lant confirm stack [--exclude-candidates]`.
  Every sub-stack of a pair is compared with the stack of all the OTHER sub-stacks of that pair (leave-one-out, so
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
"""

import os

import h5py
import numpy as np

from . import paths

MIN_SAMPLES = 10  # station sub-stacks needed before outliers of their distribution can be named
SHIFT_FLOOR = 0.25  # samples: smallest time shift listed (empirical: the precision of the interpolated peak)


def _shift_corr(a, b, w0, w1, maxs):
    """Correlation coefficient and lag (samples, sub-sample by a parabola through the peak) of b against a: the
    lag s maximising corr(a[t], b[t + s]) for t in [w0, w1), |s| <= maxs.  A positive lag: b is late."""
    x = a[w0:w1]
    y = b[w0 - maxs : w1 + maxs]
    num = np.correlate(y, x, "valid")  # num[s + maxs] = sum_t x[t] y[t + s]
    e = np.concatenate([[0.0], np.cumsum(y * y)])
    ey = e[len(x) :] - e[: -len(x)]
    c = num / np.sqrt(np.maximum((x * x).sum() * ey, 1e-300))
    k = int(np.argmax(c))
    lag = float(k - maxs)
    if 0 < k < len(c) - 1:
        den = c[k - 1] - 2 * c[k] + c[k + 1]
        if den < 0:
            lag += 0.5 * (c[k - 1] - c[k + 1]) / den
    return float(c[k]), lag


EXPECTED_BY_CHANCE = 1  # candidates expected among normal values (0.5 per test); no more: listed, not gated


def outlier_limits(dev, shift, fs):
    """Thresholds of the station diagnosis from the two distributions -> dict, or None with too few samples."""
    from scipy.stats import norm

    dv, sh = dev[np.isfinite(dev)], shift[np.isfinite(shift)]
    n = int(dv.size)
    if n < MIN_SAMPLES:
        return None
    mad = lambda v: 1.4826 * float(np.median(np.abs(v - np.median(v))))
    k1, k2 = float(norm.isf(0.5 / n)), float(norm.isf(0.25 / n))
    return dict(
        n=n,
        k_corr=round(k1, 2),
        k_shift=round(k2, 2),
        corr_median=float(np.median(dv)),
        corr_sigma=mad(dv),
        corr_min=float(np.median(dv) - k1 * mad(dv)),
        shift_median=float(np.median(sh)),
        shift_sigma=mad(sh),
        shift_max=float(max(k2 * mad(sh), SHIFT_FLOOR / fs)),
    )


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
    use_p = np.array(sorted({k for v in near.values() for _, k in sorted(v)[: int(q["partners"])]}), int)
    tsearch = 1.0 / np.sqrt(float(h.attrs["fmin"]) * float(h.attrs["fmax"]))
    maxs = max(2, int(round(tsearch * fs)))
    percorr = np.full((nday, len(use_p)), np.nan)
    pershift = np.full((nday, len(use_p)), np.nan)
    step = max(1, int(2e8 // (8 * max(nday, 1) * nl)))  # pairs read at a time (about 200 MB)
    for c0 in range(0, len(use_p), step):
        ids = use_p[c0 : c0 + step]
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
    sd_shift_se = np.full((len(stas), nday), np.nan)  # standard error of that median over the pairs
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
    zc = np.arctanh(np.clip(sd_corr, -0.999999, 0.999999))  # Fisher's z: near-normal for correlation coefficients
    with np.errstate(all="ignore"):
        day_med = np.nanmedian(sd_corr, axis=0)
        dev = zc - np.nanmedian(zc, axis=0)[None, :]
    lim = outlier_limits(dev, sd_shift, fs)
    np.savez(
        os.path.join(out_dir, "station_qc.npz"),
        stations=np.array(stas),
        blocks=np.array(days),
        corr=sd_corr,
        shift=sd_shift,
        shift_se=sd_shift_se,
        block_median=day_med,
    )
    flags = {}
    if lim:
        for i, s in enumerate(stas):
            for dd in range(nday):
                if not np.isfinite(dev[i, dd]):
                    continue
                why = []
                if dev[i, dd] < lim["corr_min"]:
                    why.append("corr %.2f (array %.2f)" % (sd_corr[i, dd], day_med[dd]))
                if (
                    abs(sd_shift[i, dd] - lim["shift_median"]) > lim["shift_max"]
                    and abs(sd_shift[i, dd] - lim["shift_median"]) > lim["k_shift"] * sd_shift_se[i, dd]
                ):
                    why.append("shift %+.4f s (+-%.4f s over its pairs)" % (sd_shift[i, dd], sd_shift_se[i, dd]))
                if why:
                    flags[(s, days[dd])] = "; ".join(why)
    judged = int(np.isfinite(sd_corr).sum())
    L = [
        "station diagnosis (lant stack) from %s" % stack_path,
        "  each sub-stack against the stack of the other sub-stacks of its pair; band %.4g-%.4g Hz; lags within "
        "distance / %.3g km/s; shift searched +-%.3g s"
        % (float(cfg["correlation"]["fmin_hz"]), float(cfg["correlation"]["fmax_hz"]), 0.6 * v_lo, maxs / fs),
        "  per station: its pairs with the %d nearest stations (%d pairs in all), at least %d with data"
        % (q["partners"], len(use_p), q["min_partners"]),
        "  %d stations, %d sub-stacks: %d station sub-stacks judged, %d with data but fewer than %d pairs (not judged)"
        % (len(stas), nday, judged, int(present.sum()) - judged, q["min_partners"]),
    ]
    if lim:
        L += [
            "  correlation as Fisher's z = atanh(r), minus the array median of the sub-stack: median %+.3f, sigma "
            "%.3f (1.4826 MAD); candidate "
            "below %+.3f (%.2f sigma: Chauvenet's criterion for N = %d)"
            % (lim["corr_median"], lim["corr_sigma"], lim["corr_min"], lim["k_corr"], lim["n"]),
            "  time shift: median %+.5f s, sigma %.5f s; candidate beyond +-%.5f s of the median (%.2f sigma, at "
            "least %.2f sample) and beyond %.2f standard errors of the station's own median over its pairs (a "
            "shift common to the pairs of a station, not a scatter of poorly measured ones)"
            % (lim["shift_median"], lim["shift_sigma"], lim["shift_max"], lim["k_shift"], SHIFT_FLOOR, lim["k_shift"]),
        ]
    else:
        L.append(
            "  fewer than %d station sub-stacks could be judged: no thresholds, nothing listed (NOT CHECKED)"
            % MIN_SAMPLES
        )
    L += [
        "",
        "  candidates: %d   (about %d expected by chance alone among normal values: Chauvenet's criterion lets "
        "0.5 through per test, correlation and time shift)" % (len(flags), EXPECTED_BY_CHANCE),
    ]
    for (s, d), why in sorted(flags.items()):
        L.append("    %s %s  %s" % (s, d, why))
    open(os.path.join(out_dir, "station_qc.txt"), "w").write("\n".join(L) + "\n")
    _station_figure(out_dir, stas, days, dev, sd_shift, flags, lim)
    return flags, dict(
        evaluated=judged,
        not_judged=int(present.sum()) - judged,
        candidates=len(flags),
        limits=lim,
        checked=bool(lim),
        pairs=int(len(use_p)),
    )


def _station_figure(out_dir, stas, days, dev, sd_shift, flags, lim):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    nday = len(days)
    fig, axs = plt.subplots(1, 2, figsize=(15, 6))
    x = np.tile(np.arange(nday), len(stas))
    axs[0].scatter(x, dev.ravel(), s=5, c="0.6")
    axs[1].scatter(x, sd_shift.ravel(), s=5, c="0.6")
    for s, d in flags:
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
