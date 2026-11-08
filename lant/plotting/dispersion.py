"""Dispersion figures.

  pair_example  one pair: its symmetric cross-correlation with the velocity window marked, the phase-velocity image
                (real part of the analytic signal at t = d / c, phase constant removed, 0-1 shown) and the group-velocity
                image (squared envelope) with the accepted picks
  all_curves    every exported phase-velocity curve in grey, mean and standard deviation (>= 20 curves), period axis on
                top
  pair_figures  the pair example for a sample of pairs (06_dispersion/pair_figures), for checking the picks by eye
"""
import os

import h5py
import numpy as np
from scipy.signal import butter, sosfiltfilt

from . import style
from .style import plt
from .. import disp as D
from .. import signal as S
from ..export import read_curves

XT = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20]


def _freq_axis(ax, f):
    ax.set_xscale("log")
    t = [x for x in XT if f[0] * 0.99 <= x <= f[-1] * 1.01]
    ax.set_xticks(t)
    ax.set_xticklabels(["%g" % x for x in t])
    ax.minorticks_off()
    ax.set_xlabel("Frequency (Hz)")
    ax.set_xlim(f[0], f[-1])


def _images(cfg, disp_path, stack_path, i):
    from ..config import phase_constant
    with h5py.File(disp_path, "r") as h:
        freqs = h["freqs"][:]
        vel = h["vel_axis"][:]
        okc, okg, vc, vg, dist = (h[k][i] for k in ("ok_c", "ok_g", "vc", "vg", "dist_km"))
        alpha = float(h.attrs["alpha"])
    with h5py.File(stack_path, "r") as h:
        x = h["stack"][i].astype(float)
        fs, maxlag = float(h.attrs["fs"]), float(h.attrs["maxlag"])
    D.VEL = vel
    bim = D.branch_images(x[None, :], np.array([dist]), fs, maxlag, freqs, alpha)
    gimg = bim[2][1][0]
    gimg = (gimg / (gimg.max(axis=1, keepdims=True) + 1e-30)) ** 2
    an = bim[2][0][0]
    t = np.broadcast_to(dist / vel, (len(freqs), len(vel))).copy()
    cimg = np.real(S.sample_at_times(an, fs, t) * np.exp(-1j * phase_constant(cfg))) / (np.abs(an).max(axis=1, keepdims=True) + 1e-30)
    w = sosfiltfilt(butter(4, [freqs[0], min(freqs[-1], 0.45 * fs)], btype="bandpass", fs=fs, output="sos"), x - x.mean())
    c = len(w) // 2
    n = min(int(round(dist / 1.0 * fs)) + 1, c)
    sym = 0.5 * (w[c:c + n] + w[c::-1][:n])
    sym /= np.abs(sym).max() + 1e-30
    return dict(freqs=freqs, vel=vel, okc=okc, okg=okg, vc=vc, vg=vg, dist=float(dist), cimg=cimg, gimg=gimg,
                tt=np.arange(n) / fs, sym=sym)


def pair_example(cfg, disp_path, stack_path, pair, out, vwin=(4.0, 1.5)):
    style.apply(cfg["plots"]["font_scale"])
    with h5py.File(disp_path, "r") as h:
        i = [p.decode() for p in h["pairs"][:]].index(pair)
    a = _images(cfg, disp_path, stack_path, i)
    fig = plt.figure(figsize=(11.0, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 2.7], hspace=0.42, wspace=0.26)
    ax = fig.add_subplot(gs[0, :])
    t1, t2 = a["dist"] / vwin[0], a["dist"] / vwin[1]
    inw = (a["tt"] >= t1) & (a["tt"] <= t2)
    ax.plot(a["tt"], a["sym"], color="0.55", lw=0.9)
    ax.plot(a["tt"][inw], a["sym"][inw], color="#e8201a", lw=1.5)
    for tv, v, col in ((t1, vwin[0], "#d0189a"), (t2, vwin[1], "#00a0c8")):
        ax.axvline(tv, color=col, ls=":", lw=1.8, label="%g km/s" % v)
    ax.set_xlim(0, a["tt"][-1])
    ax.set_ylim(-1.1, 1.1)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.set_yticks([-1, 0, 1])
    ax.legend(loc="upper right", frameon=False, ncol=2)
    ax.text(0.99, 0.06, "%s   %.2f km" % (pair.replace("_", " – "), a["dist"]), transform=ax.transAxes, ha="right", va="bottom")
    vlo, vhi = a["vel"][0], a["vel"][-1]
    for k, (img, v, ok, lab) in enumerate(((a["cimg"], a["vc"], a["okc"], "Phase velocity (km/s)"),
                                           (a["gimg"], a["vg"], a["okg"], "Group velocity (km/s)"))):
        ax = fig.add_subplot(gs[1, k])
        pm = ax.pcolormesh(a["freqs"], a["vel"], img.T, shading="gouraud", cmap=style.ENERGY, vmin=0, vmax=1, rasterized=True)
        ax.plot(a["freqs"][ok], v[ok], "o", ms=4, mfc="k", mec="k")
        _freq_axis(ax, a["freqs"])
        ax.set_ylim(vlo, vhi)
        ax.set_ylabel(lab)
    cb = fig.colorbar(pm, cax=ax.inset_axes([1.04, 0, 0.04, 1]))
    cb.set_label("Normalised energy")
    cb.set_ticks([0, 0.5, 1])
    return style.save(fig, out, cfg)


def choose_pair(disp_path, min_points=12):
    """A representative pair: among those with a long accepted phase curve, the one with the median distance."""
    with h5py.File(disp_path, "r") as h:
        okc = h["ok_c"][:]
        dist = h["dist_km"][:]
        pairs = [p.decode() for p in h["pairs"][:]]
    n = okc.sum(1)
    cand = np.where(n >= min(min_points, n.max()))[0]
    if not cand.size:
        raise RuntimeError("no pair has an accepted phase-velocity curve")
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
            acc.setdefault(round(float(t), 4), []).append(v)
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
    sec.set_xticks([t for t in (5, 2, 1, 0.5, 0.2, 0.1) if 1 / f.max() * 0.99 <= t <= 1 / f.min() * 1.01])
    sec.set_xticklabels(["%g" % t for t in sec.get_xticks()])
    sec.minorticks_off()
    sec.set_xlabel("Period (s)")
    ax.legend(loc="upper right", frameon=False)
    ax.text(0.03, 0.05, "Rayleigh-wave phase velocity", transform=ax.transAxes, fontsize=14)
    return style.save(fig, out, cfg)


def pair_figures(cfg, disp_path, out_dir, n):
    """The pair example for pairs whose images were kept by the measurement (export.curve_figures)."""
    os.makedirs(out_dir, exist_ok=True)
    with h5py.File(disp_path, "r") as h:
        pairs = [p.decode() for p in h["pairs"][:]]
        idx = h["img_index"][:]
        okc = h["ok_c"][:]
    stack = os.path.join(cfg["project"]["output_dir"], "04_notch", "stack.h5")
    done = 0
    sub = dict(cfg, plots=dict(cfg["plots"], formats=["png"], dpi=120))
    for i in idx[: (len(idx) if n < 0 else n)]:
        if not okc[i].any():
            continue
        pair_example(sub, disp_path, stack, pairs[i], os.path.join(out_dir, pairs[i]))
        done += 1
    return dict(figures=done, directory=out_dir)
