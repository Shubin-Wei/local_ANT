"""Dispersion figures.

  pair_example  one pair: its symmetric cross-correlation with the velocity window marked, the phase-velocity image
                (real part of the analytic signal at t = d / c, phase constant removed, 0-1 shown) and the group-velocity
                image (squared envelope) with the accepted picks
  all_curves    every exported phase-velocity curve in grey, mean and standard deviation (>= 20 curves), period axis on
                top
  pair_figures  the pair example for pairs spread evenly over distance and for the conspicuous ones - the most
                rejected points and the largest deviation from the median curve - (05_dispersion/pair_figures:
                spread_*, rejected_*, deviating_*), for checking the picks by eye without looking only at good results
"""
import os

import h5py
import numpy as np
from scipy.signal import butter, sosfiltfilt

from . import style
from .style import plt
from .. import disp as D
from ..export import read_curves

XT = [m * 10.0 ** e for e in range(-3, 4) for m in (1, 2, 5)]       # 1-2-5 ticks for any band


def _freq_axis(ax, f):
    ax.set_xscale("log")
    t = [x for x in XT if f[0] * 0.99 <= x <= f[-1] * 1.01]
    if len(t) < 3:                              # a narrow band: every 1-2-3-4-5-6-8 x 10^k inside it
        fine = [m * 10.0 ** k for k in range(-3, 4) for m in (1, 2, 3, 4, 5, 6, 8)]
        t = [x for x in fine if f[0] * 0.99 <= x <= f[-1] * 1.01]
    ax.set_xticks(t)
    ax.set_xticklabels(["%g" % x for x in t])
    ax.minorticks_off()
    ax.set_xlabel("Frequency (Hz)")
    ax.set_xlim(f[0], f[-1])


def _images(cfg, disp_path, stack_path, i):
    with h5py.File(disp_path, "r") as h:
        freqs = h["freqs"][:]
        vel = h["vel_axis"][:]
        okc, vc, vg, dist = (h[k][i] for k in ("ok_c", "vc", "vg", "dist_km"))
        alpha = float(h.attrs["alpha"])
    with h5py.File(stack_path, "r") as h:
        x = h["stack"][i].astype(float)
        fs, maxlag = float(h.attrs["fs"]), float(h.attrs["maxlag"])
    D.VEL = vel
    an, gimg, _ = D.symmetric_images(x[None, :], np.array([dist]), fs, maxlag, freqs, alpha)
    an, gimg = an[0], gimg[0]
    gimg = (gimg / (gimg.max(axis=1, keepdims=True) + 1e-30)) ** 2
    # cycle branches of the measured phase (phase-matched filter, lant.disp): c_n = d / (t_p + n / f), t_p = d / vc;
    # the pick is the branch n = 0.  (The narrow-band phase read at t = d / c, the classical phase image, differs from
    # the measurement where the dispersion is strong, so it is not drawn.)
    with h5py.File(disp_path, "r") as h:
        vc_all = h["vc"][i].astype(float)
        cref = h["c_ref"][:].astype(float) if h["c_ref"].ndim == 1 else h["c_ref"][i].astype(float)
        support = h["ref_support"][:].astype(bool) if "ref_support" in h else np.ones(len(freqs), bool)
    tp = np.where(okc, dist / vc_all, np.nan)
    branches = {}
    for n in (-3, -2, -1, 1, 2, 3):
        t_n = tp + n / freqs
        branches[n] = np.where(t_n > 0, dist / np.where(t_n > 0, t_n, 1), np.nan)
    cref = np.where(support, cref, np.nan)                 # the reference only where the F-J picks support it
    w = sosfiltfilt(butter(4, [freqs[0], min(freqs[-1], 0.45 * fs)], btype="bandpass", fs=fs, output="sos"), x - x.mean())
    c = len(w) // 2
    n = min(int(round(dist / (0.8 * vel[0]) * fs)) + 1, c)
    sym = 0.5 * (w[c:c + n] + w[c::-1][:n])
    sym /= np.abs(sym).max() + 1e-30
    return dict(freqs=freqs, vel=vel, okc=okc, vc=vc, vg=vg, dist=float(dist), branches=branches, cref=cref, gimg=gimg,
                tt=np.arange(n) / fs, sym=sym)


def pair_example(cfg, disp_path, stack_path, pair, out, vwin=None):
    style.apply(cfg["plots"]["font_scale"])
    if vwin is None:
        from .ccf import figure_band
        vwin = figure_band(cfg)[1]
    with h5py.File(disp_path, "r") as h:
        i = [p.decode() for p in h["pairs"][:]].index(pair)
        cref = h["c_ref"][i].astype(float) if h["c_ref"].ndim == 2 else h["c_ref"][:].astype(float)
        fr = h["freqs"][:]
    a = _images(cfg, disp_path, stack_path, i)
    # frequencies where the path is at least one wavelength of the reference (below, the images carry no phase
    # velocity information); and not above the first accepted frequency
    okf = np.isfinite(cref) & (cref > 0) & (a["dist"] * fr / np.where(cref > 0, cref, np.inf) >= 1.0)
    f_lo = fr[okf].min() if okf.any() else a["freqs"][0]
    if a["okc"].any():
        f_lo = min(f_lo, a["freqs"][a["okc"]].min())
    keepf = a["freqs"] >= f_lo * (1 - 1e-9)
    fig = plt.figure(figsize=(11.0, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 2.7], hspace=0.42, wspace=0.26)
    ax = fig.add_subplot(gs[0, :])
    t1, t2 = a["dist"] / vwin[0], a["dist"] / vwin[1]
    inw = (a["tt"] >= t1) & (a["tt"] <= t2)
    ax.plot(a["tt"], a["sym"], color="0.55", lw=0.9)
    ax.plot(a["tt"][inw], a["sym"][inw], color="#e8201a", lw=1.5)
    for tv, v, col in ((t1, vwin[0], "#d0189a"), (t2, vwin[1], "#00a0c8")):
        ax.axvline(tv, color=col, ls=":", lw=1.8, label="%.3g km/s" % v)
    ax.set_xlim(0, a["tt"][-1])
    ax.set_ylim(-1.1, 1.1)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.set_yticks([-1, 0, 1])
    ax.legend(loc="upper right", frameon=False, ncol=2)
    ax.text(0.99, 0.06, "%s   %.2f km" % (pair.replace("_", " – "), a["dist"]), transform=ax.transAxes, ha="right", va="bottom")
    # velocity axis: the measured phase and group velocities of this pair and the reference, with a margin
    vals = np.r_[a["vc"][a["okc"]], a["vg"][a["okc"]], a["cref"][keepf][np.isfinite(a["cref"][keepf])]]
    vlo, vhi = a["vel"][0], a["vel"][-1]
    if vals.size:
        vlo = max(vlo, np.floor(10 * 0.75 * vals.min()) / 10)
        vhi = min(vhi, np.ceil(10 * 1.3 * vals.max()) / 10)
    ax = fig.add_subplot(gs[1, 0])                    # phase: cycle branches, reference, picks
    fk = a["freqs"][keepf]
    for n, cn in a["branches"].items():
        ax.plot(a["freqs"], cn, "-", color="0.7", lw=0.8, zorder=1)
    ax.plot([], [], "-", color="0.7", lw=0.8, label="Other numbers of cycles")
    ax.plot(fk, a["cref"][keepf], "--", color="#e8201a", lw=1.2, zorder=2, label="Reference (F-J)")
    ax.plot(a["freqs"][a["okc"]], a["vc"][a["okc"]], "o", ms=4, mfc="k", mec="k", zorder=3, label="Phase velocity")
    _freq_axis(ax, fk)
    ax.set_ylim(vlo, vhi)
    ax.set_ylabel("Phase velocity (km/s)")
    ax.legend(loc="best", frameon=False, fontsize="small")
    ax = fig.add_subplot(gs[1, 1])                    # group: envelope image and the measured arrival
    pm = ax.pcolormesh(fk, a["vel"], a["gimg"][keepf].T, shading="gouraud", cmap=style.ENERGY, vmin=0, vmax=1, rasterized=True)
    ax.plot(a["freqs"][a["okc"]], a["vg"][a["okc"]], "o", ms=4, mfc="k", mec="k")
    _freq_axis(ax, fk)
    ax.set_ylim(vlo, vhi)
    ax.set_ylabel("Group velocity (km/s)")
    cb = fig.colorbar(pm, cax=ax.inset_axes([1.04, 0, 0.04, 1]))
    cb.set_label("Normalised energy")
    cb.set_ticks([0, 0.5, 1])
    return style.save(fig, out, cfg)


def choose_pairs(disp_path, n=3, min_points=12):
    """n pairs with long accepted phase curves, at the 20 / 50 / 80 % quantiles of their distances (short, middle,
    long), ordered by distance; with the accepted band of each: [(pair, (f_lo, f_hi))]."""
    with h5py.File(disp_path, "r") as h:
        okc = h["ok_c"][:]
        dist = h["dist_km"][:]
        freqs = h["freqs"][:]
        pairs = [p.decode() for p in h["pairs"][:]]
    cnt = okc.sum(1)
    if cnt.max() < 1:
        raise RuntimeError("no pair has an accepted phase-velocity curve")
    cand = np.where(cnt >= min(min_points, cnt.max()))[0]
    cand = cand[np.argsort(dist[cand])]
    qs = [0.5] if n == 1 else np.linspace(0.2, 0.8, n)
    pick = sorted({int(cand[min(len(cand) - 1, int(q * len(cand)))]) for q in qs}, key=lambda j: dist[j])
    out = []
    for j in pick:
        f = freqs[okc[j]]
        lo, hi = float(f.min()), float(f.max())
        if hi / lo < 2:                               # at least an octave around the accepted band
            c = np.sqrt(lo * hi)
            lo, hi = c / np.sqrt(2), c * np.sqrt(2)
        from .ccf import nice_hz
        lo, hi = nice_hz(lo), nice_hz(hi)
        if hi <= lo:
            hi = nice_hz(2 * lo)
        out.append((pairs[j], (lo, hi)))
    return out


def choose_pair(disp_path, min_points=12):
    """A representative pair: among those with a long accepted phase curve, the one with the median distance."""
    with h5py.File(disp_path, "r") as h:
        okc = h["ok_c"][:]
        dist = h["dist_km"][:]
        pairs = [p.decode() for p in h["pairs"][:]]
    n = okc.sum(1)
    if n.max() < 1:
        raise RuntimeError("no pair has an accepted phase-velocity curve")
    cand = np.where(n >= min(min_points, n.max()))[0]
    j = cand[np.argsort(dist[cand])[len(cand) // 2]]
    return pairs[j]


def all_curves(cfg, curve_dir, out):
    style.apply(cfg["plots"]["font_scale"])
    curves = read_curves(curve_dir, "C")
    if not curves:
        raise RuntimeError("no phase-velocity curve in %s" % curve_dir)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    acc = {}
    for c in curves:
        f = 1 / c["T"]
        o = np.argsort(f)
        ax.plot(f[o], c["v"][o], "-", c="0.6", lw=0.3, alpha=0.45, rasterized=True)
        for t, v in zip(c["T"], c["v"]):
            acc.setdefault(round(float(t), 6), []).append(v)
    T = np.array(sorted(acc))
    full = np.array([len(acc[t]) >= 20 for t in T])
    mean = np.array([np.mean(acc[t]) for t in T])
    std = np.array([np.std(acc[t]) for t in T])
    f = 1 / T
    o = np.argsort(f)
    f, mean, std, full = f[o], mean[o], std[o], full[o]
    b = (np.arange(len(f)) % 4 == 2) & full
    ax.plot([], [], "-", c="0.6", lw=1.5, label="Station pairs (%d curves)" % len(curves))
    ax.errorbar(f[b], mean[b], yerr=std[b], fmt="none", ecolor="k", elinewidth=1.4, capsize=4, capthick=1.4, zorder=5)
    ax.plot(f[full], mean[full], "-", c="k", lw=2.6, zorder=6, label="Mean and standard deviation")
    _freq_axis(ax, np.array([f.min(), f.max()]))
    ax.set_ylabel("Phase velocity (km/s)")
    sec = ax.secondary_xaxis("top", functions=(lambda x: 1 / np.clip(x, 1e-3, None), lambda x: 1 / np.clip(x, 1e-3, None)))
    sec.set_xticks([t for t in XT[::-1] if 1 / f.max() * 0.99 <= t <= 1 / f.min() * 1.01])
    sec.set_xticklabels(["%g" % t for t in sec.get_xticks()])
    sec.minorticks_off()
    sec.set_xlabel("Period (s)")
    ax.legend(loc="upper right", frameon=False)
    ax.text(0.03, 0.05, "Rayleigh-wave phase velocity", transform=ax.transAxes, fontsize=14)
    return style.save(fig, out, cfg)


def pair_figures(cfg, disp_path, out_dir, n):
    """Single-pair figures: n pairs with data spread evenly over distance (the same every run), plus max(3, n // 10)
    with the most rejected candidate points and as many with the largest median deviation of their accepted phase
    velocities from the per-frequency median."""
    import glob
    from .. import paths
    os.makedirs(out_dir, exist_ok=True)
    for old in glob.glob(os.path.join(out_dir, "*")):
        os.remove(old)
    with h5py.File(disp_path, "r") as h:
        pairs = [p.decode() for p in h["pairs"][:]]
        okc = h["ok_c"][:]
        pre = h["pre_c"][:]
        vc = np.where(okc, h["vc"][:], np.nan)
        snr = h["snr"][:]
        dist = h["dist_km"][:]
    have = np.where(np.isfinite(snr).any(1))[0]
    have = have[np.argsort(dist[have], kind="stable")]
    k = min(n if n > 0 else len(have), len(have))
    pick = {}
    for j in np.unique(np.round(np.linspace(0, len(have) - 1, k)).astype(int)) if k else []:
        pick[int(have[j])] = "spread"
    m = max(3, n // 10)
    rej = (pre & ~okc).sum(1)
    for i in np.argsort(-rej)[:m]:
        if rej[i] > 0:
            pick.setdefault(int(i), "rejected")
    import warnings
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(vc, axis=0)
        dev = np.nanmedian(np.abs(vc / med[None, :] - 1), axis=1)
    dev = np.where(np.isfinite(dev), dev, -1)
    for i in np.argsort(-dev)[:m]:
        if dev[i] > 0:
            pick.setdefault(int(i), "deviating")
    stack = paths.stack_product(cfg)
    sub = dict(cfg, plots=dict(cfg["plots"], formats=["png"], dpi=120))
    done = {}
    for i, kind in pick.items():
        pair_example(sub, disp_path, stack, pairs[i], os.path.join(out_dir, "%s_%s" % (kind, pairs[i])))
        done[kind] = done.get(kind, 0) + 1
    return dict(figures=done, directory=out_dir)
