"""Cross-correlation figures.

moveout     all pairs of the final product against distance in a band (default: the screening band of the
            stack step): top, the traces averaged in distance bins of 1/250 of the longest distance as a grey
            image; bottom, wiggles averaged in bins of max(1/50 of it, 10 x it / pairs) (>= 3-10 pairs); two reference moveouts (the fast and
            slow edges of the velocity axis, figure_band) on both lag sides
stack_time  one pair: a diagnostic built from the raw hourly correlations of step 02 (not the final product:
            station-hour screening is not applied; hours of listed station-days left out): the
            stack, the hourly correlations in time order (each row scaled to its peak), and the SNR of the symmetric component for the stack of the first n
            hours against n (signal |t| = d/v_fast ... d/v_slow, noise |t| > d/v_slow + 2 periods of the band
            bottom; the same windows are shaded on the stack)
"""

import glob
import os

import h5py
import numpy as np
from scipy.signal import butter, detrend, sosfiltfilt

from .. import paths
from . import style
from .style import plt


def figure_band(cfg):
    """Band and the fast / slow reference velocities of the correlation figures: the measurement band and velocity
    axis once the reference step has set them, before that the correlation band and data.velocity_range_kms."""
    d, c = cfg["dispersion"], cfg["correlation"]
    auto = lambda v: isinstance(v, str)
    band = (
        (d["fmin_hz"], d["fmax_hz"]) if not (auto(d["fmin_hz"]) or auto(d["fmax_hz"])) else (c["fmin_hz"], c["fmax_hz"])
    )
    v0, v1 = (d["vmin"], d["vmax"]) if not (auto(d["vmin"]) or auto(d["vmax"])) else cfg["data"]["velocity_range_kms"]
    return (float(band[0]), float(band[1])), (float(v1), float(v0))


def hz(f):
    """A frequency for a figure label: two significant digits, no trailing zeros."""
    return "%g" % float("%.2g" % f)


def nice_hz(f):
    """A band edge for a figure: whole Hz from 2 Hz, 0.5 Hz from 0.5 Hz, 0.1 Hz below."""
    return float(round(f)) if f >= 2 else round(2 * f) / 2 if f >= 0.5 else max(0.1, round(10 * f) / 10)


def figure_bands(lo, hi, nmax=4):
    """Bands of the moveout figures: the measurement band with its ends rounded (to whole Hz from 2 Hz, to 0.5 Hz
    from 0.5 Hz, to 0.1 Hz below), cut at the round frequencies 1-2-5 x 10^k inside it; pieces narrower than
    2/3 octave joined to a neighbour, at most nmax bands (the narrowest joined first)."""
    lo, hi = nice_hz(lo), nice_hz(hi)
    edges = [lo] + [m * 10.0**k for k in range(-3, 4) for m in (1, 2, 5) if lo * 1.05 < m * 10.0**k < hi / 1.05] + [hi]
    edges = sorted(edges)
    while len(edges) > 2:
        w = np.log2(np.array(edges[1:]) / np.array(edges[:-1]))
        i = int(np.argmin(w))
        if w[i] >= 2 / 3 and len(edges) - 1 <= nmax:
            break
        if i == 0:
            del edges[1]
        elif i == len(w) - 1:
            del edges[-2]
        else:  # join with the narrower neighbour
            del edges[i if w[i - 1] < w[i + 1] else i + 1]
    return [(edges[k], edges[k + 1]) for k in range(len(edges) - 1)]


def band_window(disp_path, band, min_samples=20):
    """From the measurements of the disp step, for a moveout figure of `band`: the distance up to which pairs were
    measured in it (99th percentile of their distances) and the 2nd and 98th percentiles of the group velocities
    measured in it: slow edge the 2nd percentile of the group velocities, fast edge the 98th percentile of the phase
    velocities (the wave train lies between them under normal dispersion); None when the band has fewer than
    min_samples accepted samples."""
    with h5py.File(disp_path, "r") as h:
        f, ok, vg, vc, dist = h["freqs"][:], h["ok_c"][:], h["vg"][:], h["vc"][:], h["dist_km"][:]
    m = (f >= band[0]) & (f <= band[1])
    acc = ok[:, m] & np.isfinite(vg[:, m]) & np.isfinite(vc[:, m])
    if acc.sum() < min_samples:
        return None
    d = dist[acc.any(axis=1)]
    vs, vf = np.percentile(vg[:, m][acc], 2), np.percentile(vc[:, m][acc], 98)
    return dict(
        dmax=float(np.percentile(d, 99)),
        vslow=max(0.1, float(np.floor(vs * 10) / 10)),
        vfast=float(np.ceil(vf * 10) / 10),
    )


def moveout(cfg, stack_path, out, band=None, maxlag=None, vref=None, dmax=None):
    style.apply(cfg["plots"]["font_scale"])
    b, (vf, vs) = figure_band(cfg)
    band = band or b
    vref = vref or ((vf, "#ff1a1a"), (vs, "#00b4ff"))
    with h5py.File(stack_path, "r") as h:
        fs, full = float(h.attrs["fs"]), float(h.attrs["maxlag"])
        use = np.where(h["nseg"][:] > 0)[0]
        d, S = h["dist_km"][:][use], h["stack"][:][use].astype(float)
    o = np.argsort(d, kind="stable")
    d, S = d[o], S[o]
    maxlag_auto = not maxlag
    maxlag = maxlag or min(full, d.max() / vref[-1][0] + 2.0 / band[0])
    tax = np.linspace(-full, full, S.shape[1])
    keep = np.abs(tax) <= maxlag + 0.5 / fs
    t = tax[keep]
    mat = sosfiltfilt(
        butter(4, list(band), btype="bandpass", fs=fs, output="sos"), detrend(S, type="constant", axis=1), axis=1
    )[:, keep]
    mat /= np.abs(mat).max(axis=1, keepdims=True) + 1e-20
    BIN, IMG = max(d.max() / 50, 10 * d.max() / len(d)), d.max() / 250  # about 10 pairs per wiggle bin at least
    min_pairs = max(3, min(10, len(d) // 50))
    edges = np.arange(0.0, d.max() + BIN, BIN)
    ib0 = np.clip(np.digitize(d, edges) - 1, 0, len(edges) - 2)
    full_bins = [i for i in range(len(edges) - 1) if (ib0 == i).sum() >= min_pairs]
    dtop = edges[full_bins[-1] + 1] if full_bins else d.max()  # the distances with enough pairs: no
    if dmax:  # and only as far as pairs were measured
        dtop = min(dtop, dmax * 1.05)
    sel = d <= dtop  # empty or single-pair rows at the top
    d, mat = d[sel], mat[sel]
    if maxlag_auto:
        nkeep = np.abs(t) <= min(full, dtop / vref[-1][0] + 2.0 / band[0]) + 0.5 / fs
        t, mat = t[nkeep], mat[:, nkeep]
    ib = np.clip(np.digitize(d, edges) - 1, 0, len(edges) - 2)
    ctr, B = [], []
    for i in range(len(edges) - 1):
        m = ib == i
        if m.sum() >= min_pairs:
            tr = mat[m].mean(axis=0)
            B.append(tr / (np.abs(tr).max() + 1e-20))
            ctr.append(d[m].mean())
    ie = np.arange(0.0, d.max() + IMG, IMG)
    k = np.clip(np.digitize(d, ie) - 1, 0, len(ie) - 2)
    img = np.full((len(ie) - 1, t.size), np.nan)
    for i in np.unique(k):
        tr = mat[k == i].mean(axis=0)
        img[i] = tr / (np.abs(tr).max() + 1e-20)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.2, 8.0), sharex=True, gridspec_kw={"hspace": 0.06})
    clip = np.nanpercentile(np.abs(img), 98)
    te = np.r_[t[0] - 0.5 / fs, 0.5 * (t[:-1] + t[1:]), t[-1] + 0.5 / fs]
    a1.pcolormesh(te, ie, img, cmap="gray_r", vmin=-clip, vmax=clip, shading="flat", rasterized=True)
    a1.set_facecolor("0.5")
    for c, tr in zip(ctr, B):
        a2.fill_between(t, c, c + tr * BIN, where=tr > 0, color="0.35", lw=0)
        a2.plot(t, c + tr * BIN, color="k", lw=0.4)
    ylim = (0, dtop)
    for ax in (a1, a2):
        ax.set_xlim(t[0], t[-1])
        ax.set_ylim(ylim)
        ax.set_ylabel("Interstation distance (km)")
        for v, c in vref:
            for sgn in (-1, 1):
                ax.plot(sgn * np.array(ylim) / v, ylim, ls="--", c=c, lw=1.6)
    a1.text(
        0.015,
        0.035,
        "%s–%s Hz, %d pairs" % (hz(band[0]), hz(band[1]), len(d)),
        transform=a1.transAxes,
        fontsize=13,
        fontweight="bold",
        color="k",
        bbox=dict(fc="w", ec="none", alpha=0.7),
    )
    for i, (v, c) in enumerate(vref):
        a2.text(0.985, 0.035 + 0.07 * (len(vref) - 1 - i), "%.3g km/s" % v, transform=a2.transAxes, ha="right", color=c)
    a2.set_xlabel("Lag time (s)")
    return style.save(fig, out, cfg)


def _hourly(cfg, pair, bad_days):
    s1, s2 = pair.split("_")
    X, W, lab = [], [], []
    for f in sorted(glob.glob(paths.P(cfg, paths.CCF, "day_*.h5"))):
        day = os.path.basename(f)[4:11]
        if (s1, day) in bad_days or (s2, day) in bad_days:
            continue
        with h5py.File(f, "r") as h:
            pn = [q.decode() for q in h["pairs"][:]]
            if pair not in pn:
                continue
            k = pn.index(pair)
            cnt = h["count"][:, k]
            x = h["hourly"][:, k, :]
            fs, maxlag = float(h.attrs["fs"]), float(h.attrs["maxlag"])
        for hr in range(24):
            if cnt[hr] > 0:
                X.append(x[hr].astype(float))
                W.append(float(cnt[hr]))
                lab.append(day)
    if not X:
        raise RuntimeError("pair %s has no hourly correlation" % pair)
    return np.array(X), np.array(W), lab, fs, maxlag


def stack_time(cfg, stack_path, pair, out, band=None, bad_days=(), vwin=None):
    style.apply(cfg["plots"]["font_scale"])
    b, (vf, vs) = figure_band(cfg)
    band = band or b
    if vwin:
        vf, vs = vwin
    with h5py.File(stack_path, "r") as h:
        pairs = [p.decode() for p in h["pairs"][:]]
        i = pairs.index(pair)
        d = float(h["dist_km"][i])
        lon = np.r_[h["lon1"][:], h["lon2"][:]]
        lat = np.r_[h["lat1"][:], h["lat2"][:]]
        p1, p2 = (float(h["lon1"][i]), float(h["lat1"][i])), (float(h["lon2"][i]), float(h["lat2"][i]))
    X, W, lab, fs, maxlag = _hourly(cfg, pair, set(bad_days))
    nh, nt = X.shape
    t = np.arange(nt) / fs - maxlag
    sos = butter(4, list(band), btype="bandpass", fs=fs, output="sos")
    filt = lambda x: sosfiltfilt(sos, detrend(x, type="constant", axis=-1), axis=-1)
    sym = lambda y: 0.5 * (y + y[..., ::-1])
    cum = np.cumsum(X * W[:, None], axis=0) / np.cumsum(W)[:, None]
    sig = (np.abs(t) >= d / vf) & (np.abs(t) <= d / vs)
    noi = np.abs(t) > d / vs + 2.0 / band[0]
    pos = t > 0
    snr = lambda y: np.abs(y[..., sig & pos]).max(axis=-1) / (y[..., noi & pos].std(axis=-1) + 1e-20)
    XF, CF = filt(X), filt(cum)
    stack = CF[-1]
    tmax = float(min(maxlag, d / vs + 7.0 / band[0]))
    m = np.abs(t) <= tmax
    fig = plt.figure(figsize=(12.5, 6.8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1], height_ratios=[1, 4.0], hspace=0.05, wspace=0.24)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)
    ax2 = fig.add_subplot(gs[:, 1])
    ax0.plot(t[m], stack[m] / np.abs(stack[m]).max(), "k-", lw=1.0)
    for sgn in (-1, 1):
        ax0.axvspan(sgn * d / vf, sgn * d / vs, color="0.78", lw=0, zorder=0)  # signal window
        ax0.axvspan(sgn * (d / vs + 2.0 / band[0]), sgn * tmax, color="#cfe3f5", lw=0, zorder=0)  # noise window
    ax0.set_ylim(-1.15, 1.15)
    ax0.set_yticks([])
    ax0.tick_params(labelbottom=False)
    ax0.text(0.0, 1.04, "%d-hour stack" % nh, transform=ax0.transAxes, va="bottom")
    ax0.text(
        1.0,
        1.04,
        "%s   %.2f km, %s–%s Hz" % (pair.replace("_", " – "), d, hz(band[0]), hz(band[1])),
        transform=ax0.transAxes,
        ha="right",
        va="bottom",
    )
    R = XF[:, m] / (np.abs(XF[:, m]).max(axis=1, keepdims=True) + 1e-20)
    ax1.imshow(
        R,
        aspect="auto",
        cmap="RdBu_r",
        vmin=-1,
        vmax=1,
        extent=(t[m][0], t[m][-1], nh + 0.5, 0.5),
        interpolation="nearest",
        rasterized=True,
    )
    first = [k for k in range(nh) if k == 0 or lab[k] != lab[k - 1]]
    every = max(1, len(first) // 8)
    ax1.set_yticks([k + 1 for k in first][::every])
    ax1.set_yticklabels([str(int(lab[k][4:])) for k in first][::every])
    ax1.set_ylabel("Day of year")
    ax1.set_xlabel("Lag time (s)")
    ax1.set_xlim(-tmax, tmax)
    ins = ax1.inset_axes([0.012, 0.015, 0.25, 0.30])
    ins.plot(lon, lat, ".", ms=2.5, color="0.5")
    ins.plot(*zip(p1, p2), "-", color=style.RED, lw=1.4)
    ins.plot(*zip(p1, p2), "^", ms=7, mfc=style.RED, mec="k", mew=0.6)
    ins.set_aspect(1 / np.cos(np.radians(np.mean(lat))), adjustable="datalim")
    ins.set_xticks([])
    ins.set_yticks([])
    n = np.arange(1, nh + 1)
    one, acc = snr(sym(XF)), snr(sym(CF))
    ax2.loglog(n, one, ".", ms=4, color="0.5", alpha=0.35)
    ax2.loglog(n, acc, "-", lw=2.4, color=style.RED, label="Stack of the first n hours")
    ax2.plot([], [], ".", ms=8, color="0.5", label="Single hours")
    ax2.set_xlim(0.8, nh * 1.25)
    ax2.set_xlabel("Stacking length n (hours)")
    ax2.set_ylabel("SNR of the symmetric component")
    from matplotlib.ticker import FuncFormatter, NullFormatter

    for axis in (ax2.xaxis, ax2.yaxis):
        axis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
        axis.set_minor_formatter(NullFormatter())
    ax2.grid(which="major", alpha=0.3, lw=0.6)
    ax2.legend(loc="upper left", frameon=False)
    files = style.save(fig, out, cfg)
    return dict(
        files=files, pair=pair, dist_km=d, hours=nh, snr_final=float(acc[-1]), snr_single_median=float(np.median(one))
    )
