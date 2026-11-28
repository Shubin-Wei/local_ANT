"""Figures of step 2 (2_cc/figures):

ccf_<band>Hz.png    the stacked correlations against distance in sub-bands of the correlation band (at most a factor 4
                    wide each): above, an image of the traces averaged in narrow distance bins; below, wiggles averaged
                    in wider bins.  Dashed lines: the fastest and slowest velocity of data.velocity_range_kms.  Surface
                    waves appear as a moveout between them, on both sides when the noise comes from all directions.
ccf_time_<pair>.png three pairs (short, middle, long path) in the middle sub-band: the stack with its signal window
                    (between the two velocities) and noise window, the sub-stacks in time order below it (a stable
                    arrival means a stable correlation; a sub-stack that differs points at a problem of a station on
                    that day or at a strong local source), and the signal-to-noise ratio of the symmetric component of
                    the stack of the first n sub-stacks.
The lag axis ends a little after the slowest velocity reaches the longest distance shown.
"""

import math
import os

import h5py
import numpy as np
from scipy.signal import butter, sosfiltfilt

from . import style
from .style import plt

CHUNK = 2000  # pairs read and filtered at a time


def bands(fmin, fmax, ratio=4.0):
    """The correlation band cut into geometrically equal sub-bands at most `ratio` wide."""
    n = max(1, math.ceil(math.log(fmax / fmin) / math.log(ratio) - 1e-9))
    e = fmin * (fmax / fmin) ** (np.arange(n + 1) / n)
    return [(float("%.3g" % e[k]), float("%.3g" % e[k + 1])) for k in range(n)]


def lag_window(dist, v_slow, maxlag, extra=1.1):
    """Lag shown: the slowest velocity's arrival at `dist` times `extra`, rounded up to 1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6
    or 8 x 10^n."""
    t = max(extra * dist / v_slow, 1e-3)
    e = 10.0 ** math.floor(math.log10(t))
    nice = min(m * e for m in (1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10) if m * e >= t * (1 - 1e-9))
    return float(min(maxlag, nice))


def _sos(fs, band):
    lo, hi = band[0], min(band[1], 0.45 * fs)
    return butter(4, (lo, hi), btype="band", fs=fs, output="sos")


def _edges(x):
    e = np.empty(len(x) + 1)
    e[1:-1] = 0.5 * (x[:-1] + x[1:])
    e[0], e[-1] = x[0] - 0.5 * (x[1] - x[0]), x[-1] + 0.5 * (x[-1] - x[-2])
    return e


def make(cfg, step):
    style.apply()
    dpi = int(cfg["plots"]["dpi"])
    with h5py.File(step.path("stack.h5"), "r") as h:
        fs, maxlag = float(h.attrs["fs"]), float(h.attrs["maxlag"])
        fmin, fmax = float(h.attrs["fmin"]), float(h.attrs["fmax"])
        cnt, dist = h["count"][:], h["dist_km"][:]
        idx = np.where(cnt > 0)[0]
        if idx.size < 2:
            return
        names = ["%s_%s" % (a.decode(), b.decode()) for a, b in zip(h["sta1"][:], h["sta2"][:])]
        xy = np.c_[np.r_[h["lon1"][:], h["lon2"][:]], np.r_[h["lat1"][:], h["lat2"][:]]]
        v_slow, v_fast = (float(v) for v in cfg["data"]["velocity_range_kms"])
        dmax = float(dist[idx].max())
        tmax = lag_window(dmax, v_slow, maxlag)
        bs = bands(fmin, fmax)
        nlag = h["ccf"].shape[1]
        lag = (np.arange(nlag) - nlag // 2) / fs
        keep = np.abs(lag) <= tmax + 0.5 / fs
        few = idx.size < 2000
        img_e = np.linspace(0, dmax * 1.0001, 251 if not few else 101)
        wig_e = np.linspace(0, dmax * 1.0001, 51 if not few else 31)
        acc = {b: [np.zeros((len(img_e) - 1, keep.sum())), np.zeros(len(img_e) - 1), np.zeros((len(wig_e) - 1,
               keep.sum())), np.zeros(len(wig_e) - 1), np.zeros(len(wig_e) - 1)] for b in bs}  # fmt: skip
        order = idx[np.argsort(dist[idx])]
        for c0 in range(0, order.size, CHUNK):
            sel = np.sort(order[c0 : c0 + CHUNK])
            x = h["ccf"][sel].astype(np.float64)
            d = dist[sel]
            for b in bs:
                y = sosfiltfilt(_sos(fs, b), x - x.mean(axis=1, keepdims=True), axis=1)[:, keep]
                y /= np.abs(y).max(axis=1, keepdims=True) + 1e-30
                a = acc[b]
                i = np.clip(np.digitize(d, img_e) - 1, 0, len(img_e) - 2)
                np.add.at(a[0], i, y)
                np.add.at(a[1], i, 1)
                j = np.clip(np.digitize(d, wig_e) - 1, 0, len(wig_e) - 2)
                np.add.at(a[2], j, y)
                np.add.at(a[3], j, 1)
                np.add.at(a[4], j, d)
    t = lag[keep]
    minpairs = 3 if few else 10
    for b in bs:
        section(acc[b], img_e, wig_e, t, b, (v_fast, v_slow), minpairs, step.fig("ccf_%g-%gHz" % b), dpi)
    mid = bs[len(bs) // 2]
    q = idx[np.argsort(dist[idx])]
    best = cnt[q] >= 0.8 * cnt[q].max()  # among the pairs recorded about as long as the best ones
    picks = []
    for f in (0.25, 0.5, 0.9):
        target = dist[q][int(f * (q.size - 1))]
        cand = q[best] if best.any() else q
        p = int(cand[np.argmin(np.abs(dist[cand] - target))])
        if p not in picks:
            picks.append(p)
    for p in picks:
        try:
            sub_stacks(step, p, names[p], float(dist[p]), xy, fs, maxlag, mid, (v_fast, v_slow), dpi)
        except Exception as exc:  # a figure must not stop the step
            step.say("  figure of %s not drawn: %s" % (names[p], exc))


def section(a, img_e, wig_e, t, band, vel, minpairs, out, dpi):
    img = np.where(a[1][:, None] > 0, a[0] / np.maximum(a[1], 1)[:, None], np.nan)
    img /= np.nanmax(np.abs(img), axis=1, keepdims=True) + 1e-30
    ok = a[3] >= minpairs
    wig = a[2][ok] / a[3][ok][:, None]
    wig /= np.abs(wig).max(axis=1, keepdims=True) + 1e-30
    dw = a[4][ok] / a[3][ok]
    h = wig_e[1] - wig_e[0]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.6, 8.0), sharex=True, gridspec_kw={"hspace": 0.06})
    clip = np.nanpercentile(np.abs(img), 98) if np.isfinite(img).any() else 1.0
    a1.set_facecolor("0.5")
    a1.pcolormesh(_edges(t), img_e, img, cmap="gray_r", vmin=-clip, vmax=clip, shading="flat", rasterized=True)
    for c, y in zip(dw, wig):
        a2.fill_between(t, c, c + y * h, where=y > 0, color="0.35", lw=0)
        a2.plot(t, c + y * h, color="k", lw=0.4)
    ylim = (0, img_e[-1] + h)
    for ax in (a1, a2):
        ax.set_xlim(t[0], t[-1])
        ax.set_ylim(*ylim)
        ax.set_ylabel("Distance (km)")
        for v, col in zip(vel, (style.FAST, style.SLOW)):
            for sgn in (-1, 1):
                ax.plot(sgn * np.array(ylim) / v, ylim, ls="--", c=col, lw=1.6)
    a1.text(0.015, 0.03, "%g–%g Hz" % band, transform=a1.transAxes, ha="left", va="bottom", fontsize=16,
            fontweight="bold", color="k", bbox=dict(fc="white", ec="none", alpha=0.8, pad=2))  # fmt: skip
    for i, (v, col) in enumerate(zip(vel, (style.FAST, style.SLOW))):
        a2.text(0.985, 0.03 + 0.07 * (1 - i), "%g km/s" % v, transform=a2.transAxes, ha="right", va="bottom",
                fontsize=14, color=col)  # fmt: skip
    a2.set_xlabel("Lag time (s)")
    style.save(fig, out, dpi)


def _pair_rows(step, p):
    """Sub-stacks of pair p in time order -> (rows, segment counts, labels)."""
    rows, w, labels = [], [], []
    for f in sorted(f for f in os.listdir(step.path("days")) if f.endswith(".h5")):
        with h5py.File(step.path("days", f), "r") as h:
            gp = h["pair"][:]
            k = np.searchsorted(gp, p)
            if k >= gp.size or gp[k] != p:
                continue
            ns = h["count"].shape[0]
            for s in range(ns):
                n = h["count"][s, k]
                if n > 0:
                    rows.append(h["ccf"][s, k].astype(np.float64))
                    w.append(float(n))
                    labels.append(f[:-3] if ns == 1 else "%s %02dh" % (f[:-3], s * 24 // ns))
    return np.array(rows), np.array(w), labels


def sub_stacks(step, p, name, d, xy, fs, maxlag, band, vel, dpi):
    X, W, labels = _pair_rows(step, p)
    n = len(X)
    if n < 2:
        return
    v_fast, v_slow = vel
    lag = (np.arange(X.shape[1]) - X.shape[1] // 2) / fs
    sos = _sos(fs, band)
    filt = lambda x: sosfiltfilt(sos, x - x.mean(axis=-1, keepdims=True), axis=-1)  # noqa: E731
    cum = np.cumsum(X * W[:, None], axis=0) / np.cumsum(W)[:, None]
    XF, CF = filt(X), filt(cum)
    s0, s1 = d / v_fast, d / v_slow
    n0 = s1 + 2.0 / band[0]
    pos = lag >= 0
    sym = lambda y: 0.5 * (y[..., pos] + y[..., ::-1][..., pos])  # noqa: E731
    tp = lag[pos]
    sig, noi = (tp >= s0) & (tp <= s1), tp >= min(n0, 0.85 * maxlag)
    snr = lambda y: np.abs(y[..., sig]).max(axis=-1) / (y[..., noi].std(axis=-1) + 1e-30)  # noqa: E731
    one, acc = snr(sym(XF)), snr(sym(CF))
    tmax = lag_window(max(d, 1e-3), v_slow, maxlag, extra=1.6)
    tmax = max(tmax, min(maxlag, n0 + 2.0 / band[0]))
    m = np.abs(lag) <= tmax
    fig = plt.figure(figsize=(13.0, 7.0))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1], height_ratios=[1, 4.0], hspace=0.05, wspace=0.24)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)
    ax2 = fig.add_subplot(gs[:, 1])
    st = CF[-1]
    ax0.plot(lag[m], st[m] / (np.abs(st[m]).max() + 1e-30), "k-", lw=1.0)
    for sgn in (-1, 1):
        ax0.axvspan(sgn * s0, sgn * s1, color="0.8", lw=0, zorder=0)
        if n0 < tmax:
            ax0.axvspan(sgn * n0, sgn * tmax, color="#cfe3f5", lw=0, zorder=0)
    ax0.set_ylim(-1.15, 1.15)
    ax0.set_yticks([])
    ax0.tick_params(labelbottom=False)
    ax0.text(0.012, 0.92, "stack of %d" % n, transform=ax0.transAxes, ha="left", va="top", fontsize=13)
    ax0.text(0.988, 0.92, "%s, %.2f km, %g–%g Hz" % (name.replace("_", "–"), d, *band),
             transform=ax0.transAxes, ha="right", va="top", fontsize=13)  # fmt: skip
    R = XF[:, m] / (np.abs(XF[:, m]).max(axis=1, keepdims=True) + 1e-30)
    ax1.imshow(R, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1, extent=(lag[m][0], lag[m][-1], n + 0.5, 0.5),
               interpolation="nearest", rasterized=True)  # fmt: skip
    every = max(1, int(math.ceil(n / 12)))
    ax1.set_yticks(np.arange(1, n + 1)[::every])
    ax1.set_yticklabels(labels[::every], fontsize=10)
    ax1.set_xlabel("Lag time (s)")
    ax1.set_xlim(-tmax, tmax)
    ins = ax1.inset_axes([0.012, 0.015, 0.22, 0.28])
    kx = np.cos(np.radians(xy[:, 1].mean()))
    u = np.unique(np.round(xy, 6), axis=0)
    ins.plot(u[:, 0], u[:, 1], ".", ms=2, color="0.55")
    half = len(xy) // 2
    p1, p2 = xy[p], xy[half + p]
    ins.plot(*zip(p1, p2), "-", color=style.RED, lw=1.3)
    ins.plot(*zip(p1, p2), "^", ms=7, mfc=style.RED, mec="k", mew=0.5)
    ins.set_aspect(1 / kx, adjustable="datalim")
    ins.set_xticks([])
    ins.set_yticks([])
    k = np.arange(1, n + 1)
    ax2.loglog(k, one, ".", ms=5, color="0.5", alpha=0.5, label="single sub-stack")
    ax2.loglog(k, acc, "-", lw=2.4, color=style.RED, label="stack of the first n")
    ax2.set_xlim(0.8, n * 1.25)
    ax2.set_xlabel("Number of sub-stacks n")
    ax2.set_ylabel("SNR of the symmetric component")
    ax2.grid(which="major", alpha=0.3, lw=0.6)
    style.plain_log(ax2.xaxis)
    style.plain_log(ax2.yaxis)
    ax2.legend(loc="upper left")
    style.save(fig, step.fig("ccf_time_" + name), dpi)
