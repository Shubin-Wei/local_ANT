"""Figures of step 3 (3_disp/figures):

fj_reference.png      the array's F-J spectrum (normalised at each frequency), its picks, the reference curve used and
                      the median of the accepted phase velocities: the check that the reference follows the
                      fundamental mode
phase_curves.png      every accepted phase-velocity curve (grey), their median with the standard deviation, the reference
group_curves.png      the same for the group velocities
rejections.png        per frequency: the accepted phase and group velocities and the values rejected for each reason
pair_<pair>.png       three pairs (short, middle, long path): the symmetric correlation with the window between the
                      fastest and slowest velocity; the phase-velocity image (real part of the narrow-band signal at
                      t = dist / c, phase-shifted by the far-field pi/4: its crests are the 2 pi branches) and the group
                      image (squared envelope at t = dist / U), each with the accepted (black) and rejected (grey)
                      values and the reference (dashed)
"""

import os

import h5py
import numpy as np
from scipy.signal import butter, sosfiltfilt

from .. import dispersion as D
from ..signal import analytic_bank, branches, sample_at_times
from . import style
from .style import plt


def make(cfg, step, ref, m, ok, ok_u, freqs, dist, sta1, sta2):
    style.apply()
    dpi = int(cfg["plots"]["dpi"])
    vrange = [float(v) for v in cfg["data"]["velocity_range_kms"]]
    with np.errstate(all="ignore"):
        med = np.array([np.median(m["c"][ok[:, k], k]) if ok[:, k].sum() else np.nan for k in range(len(freqs))])
    fj_reference(ref, med, vrange, step.fig("fj_reference"), dpi)
    curves(m["c"], ok, freqs, ref["c"], "Phase velocity (km/s)", step.fig("phase_curves"), dpi)
    curves(m["u"], ok_u, freqs, ref["u"], "Group velocity (km/s)", step.fig("group_curves"), dpi)
    rejections(m, ok, ok_u, freqs, np.isfinite(ref["c"]), step.fig("rejections"), dpi)
    n_ok = ok.sum(axis=1)
    have = np.where(n_ok >= 3)[0]
    if have.size:
        by = have[np.argsort(dist[have])]
        for f in (0.2, 0.5, 0.85):
            part = by[int(f * (by.size - 1)) : int(f * (by.size - 1)) + max(1, by.size // 10)]
            q = int(part[np.argmax(n_ok[part])])
            try:
                pair(cfg, step, q, "%s_%s" % (sta1[q], sta2[q]), float(dist[q]), freqs, ref, m, ok, ok_u, vrange, dpi)
            except Exception as exc:  # a figure must not stop the step
                step.say("  figure of %s_%s not drawn: %s" % (sta1[q], sta2[q], exc))


def fj_reference(ref, med, vrange, out, dpi):
    f, c, img = ref["freqs"], ref["fj_c"], ref["fj_image"]
    Y = img / (img.max(axis=1, keepdims=True) + 1e-30)
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    pm = ax.pcolormesh(f, c, Y.T, cmap=style.ENERGY, vmin=0, vmax=1, shading="gouraud", rasterized=True)
    ax.plot(f, ref["fj_picks"], "o", ms=5, mfc="none", mec="k", mew=1.0, label="F-J maximum")
    ax.plot(f, ref["c"], "k-", lw=2.2, label="Reference (%s)" % ref["source"])
    ax.plot(f, med, "D", ms=4.5, mfc="#6a3d9a", mec="white", mew=0.5, label="Median of the accepted pairs")
    style.freq_axis(ax, f)
    ax.set_ylim(c[0], c[-1])
    ax.set_ylabel("Phase velocity (km/s)")
    ax.legend(loc="upper right", fontsize=11)
    cb = fig.colorbar(pm, ax=ax, pad=0.02, ticks=[0, 0.5, 1])
    cb.set_label("Normalised F-J amplitude")
    style.save(fig, out, dpi)


def curves(v, ok, freqs, v_ref, label, out, dpi):
    rows = np.where(ok.any(axis=1))[0]
    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    step = max(1, rows.size // 4000)
    for q in rows[::step]:
        y = np.where(ok[q], v[q], np.nan)
        ax.plot(freqs, y, "-", color="0.6", lw=0.3, alpha=0.45, rasterized=True)
    n = ok.sum(axis=0)
    full = n >= max(5, min(20, rows.size // 20))
    with np.errstate(all="ignore"):
        med = np.array([np.median(v[ok[:, k], k]) if full[k] else np.nan for k in range(len(freqs))])
        sd = np.array([np.std(v[ok[:, k], k]) if full[k] else np.nan for k in range(len(freqs))])
    bars = (np.arange(len(freqs)) % 4 == 2) & full
    ax.plot([], [], "-", color="0.6", lw=1.5, label="Pairs (%d%s)" % (rows.size, "" if step == 1 else
                                                                   ", every %d-th drawn" % step))  # fmt: skip
    ax.errorbar(freqs[bars], med[bars], yerr=sd[bars], fmt="none", ecolor="k", elinewidth=1.4, capsize=4, zorder=5)
    ax.plot(freqs, med, "-", color="k", lw=2.6, zorder=6, label="Median and standard deviation")
    ax.plot(freqs, v_ref, "--", color=style.RED, lw=1.6, zorder=7, label="Reference")
    style.freq_axis(ax, freqs)
    good = v[ok]
    if good.size:
        lo, hi = np.percentile(good, [0.5, 99.5])
        pad = 0.08 * (hi - lo)
        ax.set_ylim(lo - pad, hi + pad)
    ax.set_ylabel(label)
    ax.legend(loc="upper right", fontsize=11)
    style.save(fig, out, dpi)


def rejections(m, ok, ok_u, freqs, defined, out, dpi):
    """Accepted values and the values rejected for each reason, per frequency, where the reference is defined."""
    k = np.where(defined)[0]
    if not k.size:
        return
    f = freqs[k]
    with_data = ~(m["reason"] == D.REASONS["no_data"]).all(axis=1)
    fig, axs = plt.subplots(1, 2, figsize=(13.0, 4.8), sharey=True, gridspec_kw={"wspace": 0.06})
    colours = dict(zip([r for r in D.REASONS if r != "no_data"], plt.get_cmap("tab10").colors))
    for ax, key, good, title in (
        (axs[0], "reason", ok, "Phase velocity"),
        (axs[1], "reason_u", ok_u, "Group velocity"),
    ):
        ax.plot(f, good.sum(axis=0)[k], "k-", lw=2.6, label="accepted")
        for name, bit in D.REASONS.items():
            n = ((m[key][with_data] & bit) > 0).sum(axis=0)[k]
            if name != "no_data" and n.any():
                ax.plot(f, n, "-", lw=1.4, color=colours[name], label=name.replace("_", " "))
        style.freq_axis(ax, f)
        ax.text(0.02, 0.96, title, transform=ax.transAxes, ha="left", va="top", fontsize=14)
        ax.legend(loc="upper right", fontsize=10, ncol=2)
    axs[0].set_ylabel("Pairs")
    style.save(fig, out, dpi)


def pair(cfg, step, q, name, dist, freqs, ref, m, ok, ok_u, vrange, dpi):
    with h5py.File(os.path.join(cfg["output_dir"], "2_cc", "stack.h5"), "r") as h:
        x = h["ccf"][q].astype(np.float64)
        fs, maxlag = float(h.attrs["fs"]), float(h.attrs["maxlag"])
    pos, neg = branches(x[None, :], fs, maxlag)
    sym = 0.5 * (pos + neg)
    an = analytic_bank(sym, fs, freqs, float(cfg["disp"]["alpha"]))[0]  # [nf, nt]
    v = np.linspace(0.8 * vrange[0], 1.1 * vrange[1], 300)
    t = np.broadcast_to(dist / v, (len(freqs), len(v)))
    a = sample_at_times(an, fs, np.ascontiguousarray(t))
    norm = np.abs(an).max(axis=1, keepdims=True) + 1e-30
    cimg = np.real(a * np.exp(-1j * D.PHI0)) / norm
    genv = np.abs(a) / norm
    gimg = (genv / (genv.max(axis=1, keepdims=True) + 1e-30)) ** 2
    band = (freqs[0], min(freqs[-1], 0.45 * fs))
    w = sosfiltfilt(butter(4, band, btype="band", fs=fs, output="sos"), sym[0] - sym[0].mean())
    tt = np.arange(w.size) / fs
    tend = min(maxlag, 1.6 * dist / vrange[0] + 2.0 / freqs[0])
    keep = tt <= tend
    w, tt = w[keep] / (np.abs(w[keep]).max() + 1e-30), tt[keep]
    t1, t2 = dist / vrange[1], dist / vrange[0]
    inw = (tt >= t1) & (tt <= t2)

    fig = plt.figure(figsize=(12.0, 7.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 2.7], hspace=0.55, wspace=0.22)
    ax = fig.add_subplot(gs[0, :])
    ax.plot(tt, w, color="0.55", lw=0.9)
    ax.plot(np.where(inw, tt, np.nan), np.where(inw, w, np.nan), color=style.FAST, lw=1.5)
    for tv, vv, col in ((t1, vrange[1], style.FAST), (t2, vrange[0], style.SLOW)):
        ax.axvline(tv, color=col, ls=":", lw=1.8, label="%g km/s" % vv)
    ax.set_xlim(0, tt[-1])
    ax.set_ylim(-1.1, 1.1)
    ax.set_yticks([-1, 0, 1])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.legend(loc="upper right", ncol=2, fontsize=12)
    ax.text(0.01, 0.06, "%s   %.2f km   symmetric component, %.3g–%.3g Hz" % (name.replace("_", "–"), dist,
            *band), transform=ax.transAxes, ha="left", va="bottom", fontsize=12)  # fmt: skip
    for k, (img, key, sk, good, r, lab) in enumerate((
            (cimg, "c", "sigma", ok[q], ref["c"], "Phase velocity (km/s)"),
            (gimg, "u", "u_sigma", ok_u[q], ref["u"], "Group velocity (km/s)"))):  # fmt: skip
        ax = fig.add_subplot(gs[1, k])
        pm = ax.pcolormesh(freqs, v, img.T, shading="gouraud", cmap=style.ENERGY, vmin=0, vmax=1, rasterized=True)
        ax.plot(freqs, r, "--", color="0.25", lw=1.2, label="Reference")
        val = m[key][q]
        bad = ~good & np.isfinite(val)
        ax.plot(freqs[bad], val[bad], "x", color="0.45", ms=5, mew=1.0, label="Rejected")
        ax.errorbar(freqs[good], val[good], yerr=np.nan_to_num(m[sk][q][good]), fmt="o", ms=4.5, mfc="k", mec="k",
                    ecolor="k", capsize=0, label="Accepted")  # fmt: skip
        style.log_axis(ax.xaxis)
        ax.set_xlim(freqs[0], freqs[-1])
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylim(v[0], v[-1])
        ax.set_ylabel(lab)
        if k == 0:
            ax.legend(loc="upper right", fontsize=11)
    cb = fig.colorbar(pm, cax=ax.inset_axes([1.04, 0, 0.035, 1]), ticks=[0, 0.5, 1])
    cb.set_label("Normalised amplitude")
    style.save(fig, step.fig("pair_" + name), dpi)
