"""Cross-correlation figures.

  moveout     all pairs of the final product against distance in a band: top, the traces averaged in 0.05 km
              distance bins as a grey image; bottom, wiggles averaged in 0.25 km bins (>= 10 pairs); two reference
              moveouts (4.0 and 1.5 km/s) on both lag sides
  stack_time  one pair: the stack, the hourly correlations in time order (each row scaled to its peak; hours of
              listed station-days left out), and the SNR of the symmetric component for the stack of the first n
              hours against n (signal |t| = d/4 ... d/1.2 s, noise |t| > d/1.2 + 2 s)
"""
import glob
import os

import h5py
import numpy as np
from scipy.signal import butter, detrend, sosfiltfilt

from . import style
from .style import plt


def moveout(cfg, stack_path, out, band=(1.0, 5.0), maxlag=None, vref=((4.0, "#ff1a1a"), (1.5, "#00b4ff"))):
    style.apply(cfg["plots"]["font_scale"])
    with h5py.File(stack_path, "r") as h:
        fs, full = float(h.attrs["fs"]), float(h.attrs["maxlag"])
        use = np.where(h["nseg"][:] > 0)[0]
        d, S = h["dist_km"][:][use], h["stack"][:][use].astype(float)
    o = np.argsort(d, kind="stable")
    d, S = d[o], S[o]
    maxlag = maxlag or min(full, float(np.ceil(d.max() / 1.5 + 2)))
    tax = np.linspace(-full, full, S.shape[1])
    keep = np.abs(tax) <= maxlag + 0.5 / fs
    t = tax[keep]
    mat = sosfiltfilt(butter(4, list(band), btype="bandpass", fs=fs, output="sos"), detrend(S, type="constant", axis=1), axis=1)[:, keep]
    mat /= np.abs(mat).max(axis=1, keepdims=True) + 1e-20
    BIN, IMG = max(0.25, d.max() / 50), max(0.05, d.max() / 250)
    edges = np.arange(0.0, d.max() + BIN, BIN)
    ib = np.clip(np.digitize(d, edges) - 1, 0, len(edges) - 2)
    ctr, B = [], []
    for i in range(len(edges) - 1):
        m = ib == i
        if m.sum() >= 10:
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
    ylim = (0, d.max() + BIN)
    for ax in (a1, a2):
        ax.set_xlim(t[0], t[-1])
        ax.set_ylim(ylim)
        ax.set_ylabel("Interstation distance (km)")
        for v, c in vref:
            for sgn in (-1, 1):
                ax.plot(sgn * np.array(ylim) / v, ylim, ls="--", c=c, lw=1.6)
    a1.text(0.015, 0.035, "%g–%g Hz, %d pairs" % (band[0], band[1], len(d)), transform=a1.transAxes, fontsize=13,
            fontweight="bold", color="k", bbox=dict(fc="w", ec="none", alpha=0.7))
    for i, (v, c) in enumerate(vref):
        a2.text(0.985, 0.035 + 0.07 * (len(vref) - 1 - i), "%g km/s" % v, transform=a2.transAxes, ha="right", color=c)
    a2.set_xlabel("Lag time (s)")
    return style.save(fig, out, cfg)


def _hourly(cfg, pair, bad_days):
    s1, s2 = pair.split("_")
    X, W, lab = [], [], []
    for f in sorted(glob.glob(os.path.join(cfg["project"]["output_dir"], "02_ccf", "day_*.h5"))):
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


def stack_time(cfg, stack_path, pair, out, band=(1.0, 4.0), bad_days=()):
    style.apply(cfg["plots"]["font_scale"])
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
    sig = (np.abs(t) >= d / 4.0) & (np.abs(t) <= d / 1.2)
    noi = np.abs(t) > d / 1.2 + 2.0
    pos = t > 0
    snr = lambda y: np.abs(y[..., sig & pos]).max(axis=-1) / (y[..., noi & pos].std(axis=-1) + 1e-20)
    XF, CF = filt(X), filt(cum)
    stack = CF[-1]
    tmax = float(min(maxlag, d / 1.2 + 7.0))
    m = np.abs(t) <= tmax
    fig = plt.figure(figsize=(12.5, 6.8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1], height_ratios=[1, 4.0], hspace=0.05, wspace=0.24)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)
    ax2 = fig.add_subplot(gs[:, 1])
    ax0.plot(t[m], stack[m] / np.abs(stack[m]).max(), "k-", lw=1.0)
    for sgn in (-1, 1):
        ax0.axvspan(sgn * d / 4.0, sgn * d / 1.2, color="0.78", lw=0, zorder=0)
        ax0.axvspan(sgn * (d / 1.2 + 2.0), sgn * tmax, color="#cfe3f5", lw=0, zorder=0)
    ax0.set_ylim(-1.15, 1.15)
    ax0.set_yticks([])
    ax0.tick_params(labelbottom=False)
    ax0.text(0.012, 0.9, "%d-hour stack" % nh, transform=ax0.transAxes, va="top")
    ax0.text(0.988, 0.9, "%.2f km, %g–%g Hz" % (d, band[0], band[1]), transform=ax0.transAxes, ha="right", va="top")
    R = XF[:, m] / (np.abs(XF[:, m]).max(axis=1, keepdims=True) + 1e-20)
    ax1.imshow(R, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1, extent=(t[m][0], t[m][-1], nh + 0.5, 0.5),
               interpolation="nearest", rasterized=True)
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
    return dict(files=files, pair=pair, dist_km=d, hours=nh, snr_final=float(acc[-1]), snr_single_median=float(np.median(one)))
