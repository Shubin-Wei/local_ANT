"""Figures of step 2 (2_disp/figures):

fj_reference.png   the array's F-J spectrum (normalised at each frequency), its picks, the reference curve used
                   and the median of the accepted phase velocities: the check that the reference follows the
                   fundamental mode
curves.png         every accepted curve, the reference, and the number of pairs per frequency
rejections.png     per frequency: accepted values and the values rejected for each reason
pair_<pair>.png    three pairs (short, middle, long path): the envelope of the filtered symmetric correlation
                   against group velocity with the picked arrival, and the measured phase velocity with its
                   jackknife error and the reference
"""

import h5py
import numpy as np

from .. import dispersion as D
from ..signal import analytic_bank, branches
from . import style
from .style import plt


def _faxis(ax, f):
    style.log_axis(ax.xaxis)
    ax.set_xlim(f[0] / 1.03, f[-1] * 1.03)
    ax.set_xlabel("Frequency (Hz)")


def make(cfg, step, ref, m, ok, freqs, dist, sta1, sta2):
    style.apply()
    dpi = int(cfg["plots"]["dpi"])
    with np.errstate(all="ignore"):
        med = np.array([np.median(m["c"][ok[:, k], k]) if ok[:, k].sum() else np.nan for k in range(len(freqs))])
    fj_reference(ref, med, cfg["data"]["velocity_range_kms"], step.fig("fj_reference"), dpi)
    curves(m["c"], ok, freqs, ref["c"], step.fig("curves"), dpi)
    rejections(m["reason"], ok, freqs, step.fig("rejections"), dpi)
    n_ok = ok.sum(axis=1)
    have = np.where(n_ok >= 3)[0]
    if have.size:
        by = have[np.argsort(dist[have])]
        for f in (0.2, 0.5, 0.85):
            part = by[int(f * (by.size - 1)) : int(f * (by.size - 1)) + max(1, by.size // 10)]
            q = int(part[np.argmax(n_ok[part])])
            pair(cfg, step, q, "%s_%s" % (sta1[q], sta2[q]), dist[q], freqs, ref, m, ok, dpi)


def fj_reference(ref, med, vrange, out, dpi):
    f, c, img = ref["freqs"], ref["fj_c"], ref["fj_image"]
    Y = img / (img.max(axis=1, keepdims=True) + 1e-30)
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    pm = ax.pcolormesh(f, c, Y.T, cmap=style.ENERGY, vmin=0, vmax=1, shading="nearest", rasterized=True)
    ax.plot(f, ref["fj_picks"], "o", ms=3.5, mfc="none", mec="k", mew=0.8, label="F-J maximum")
    ax.plot(f, ref["c"], "k-", lw=1.8, label="Reference (%s)" % ref["source"])
    ax.plot(f, med, "D", ms=3.5, color="#6a3d9a", label="Median of the measured pairs")
    for v in vrange:
        ax.axhline(v, color="0.3", lw=0.8, ls=":")
    _faxis(ax, f)
    ax.set_ylim(c[0], c[-1])
    ax.set_ylabel("Phase velocity (km/s)")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)
    fig.colorbar(pm, ax=ax, label="F-J amplitude (normalised per frequency)")
    style.save(fig, out, dpi)


def curves(c, ok, freqs, c_ref, out, dpi):
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(7.0, 6.5), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
    rows = np.where(ok.any(axis=1))[0]
    step = max(1, rows.size // 3000)
    for q in rows[::step]:
        y = np.where(ok[q], c[q], np.nan)
        ax.plot(freqs, y, color="0.35", lw=0.3, alpha=0.35)
    ax.plot(freqs, c_ref, color=style.RED, lw=1.8, label="Reference")
    ax.set_ylabel("Phase velocity (km/s)")
    ax.legend(fontsize=9)
    ax.set_title("%d pairs with a curve%s" % (rows.size, "" if step == 1 else " (every %d-th drawn)" % step),
                 fontsize=11)  # fmt: skip
    bx.plot(freqs, ok.sum(axis=0), "k.-", lw=1)
    bx.set_ylabel("Pairs")
    _faxis(bx, freqs)
    style.save(fig, out, dpi)


def rejections(reason, ok, freqs, out, dpi):
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    with_data = ~(reason == D.REASONS["no_data"]).all(axis=1)
    ax.plot(freqs, ok.sum(axis=0), "k-", lw=2.0, label="accepted")
    for name, bit in D.REASONS.items():
        n = ((reason[with_data] & bit) > 0).sum(axis=0)
        if n.any() and name != "no_data":
            ax.plot(freqs, n, lw=1.0, label=name.replace("_", " "))
    _faxis(ax, freqs)
    ax.set_ylabel("Values (pairs)")
    ax.legend(fontsize=7, ncol=2, loc="best")
    style.save(fig, out, dpi)


def pair(cfg, step, q, name, dist, freqs, ref, m, ok, dpi):
    import os

    with h5py.File(os.path.join(cfg["output_dir"], "1_cc", "stack.h5"), "r") as h:
        x = h["ccf"][q].astype(np.float64)
        fs, maxlag = float(h.attrs["fs"]), float(h.attrs["maxlag"])
    pos, neg = branches(x[None, :], fs, maxlag)
    sym = 0.5 * (pos + neg)
    env = np.abs(analytic_bank(sym, fs, freqs, float(cfg["disp"]["alpha"]))[0])
    t = np.arange(env.shape[1]) / fs
    v = np.linspace(0.5 * cfg["data"]["velocity_range_kms"][0], 1.3 * cfg["data"]["velocity_range_kms"][1], 300)
    img = np.array([np.interp(dist / v, t, e) for e in env])
    img /= img.max(axis=1, keepdims=True) + 1e-30
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.5, 4.3))
    ax.pcolormesh(freqs, v, img.T, cmap=style.ENERGY, vmin=0, vmax=1, shading="nearest", rasterized=True)
    ax.plot(freqs, ref["u"], "w--", lw=1.2, label="Reference group velocity")
    ax.plot(freqs, m["u"][q], "k.", ms=4, label="Group arrival")
    _faxis(ax, freqs)
    ax.set_ylabel("Group velocity (km/s)")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title("Envelope (normalised per frequency)", fontsize=10)
    good = ok[q]
    bx.plot(freqs, ref["c"], color=style.RED, lw=1.5, label="Reference")
    bx.plot(freqs[~good], m["c"][q][~good], "x", color="0.55", ms=4, label="Rejected")
    bx.errorbar(freqs[good], m["c"][q][good], yerr=np.nan_to_num(m["sigma"][q][good]), fmt="o", color="k", ms=3.5,
                capsize=2, label="Accepted (jackknife 1 sd)")  # fmt: skip
    _faxis(bx, freqs)
    bx.set_ylabel("Phase velocity (km/s)")
    bx.legend(fontsize=8)
    fig.suptitle("%s, %.3f km" % (name, dist), fontsize=11)
    style.save(fig, step.fig("pair_" + name), dpi)
