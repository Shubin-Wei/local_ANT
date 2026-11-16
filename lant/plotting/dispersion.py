"""Dispersion figures.

pair_example  one pair: symmetric cross-correlation with the velocity window, phase-matched coherence with picks
              and their integer-cycle alternatives, and the group-velocity squared-envelope image. The array
              F-J phase-velocity spectrum is a separate reference product, not single-pair energy.
all_curves    every exported phase-velocity curve in grey, mean and standard deviation (>= 20 curves), period axis on
              top
pair_figures  the pair example for pairs spread evenly over distance and for the conspicuous ones - the most
              rejected points and the largest deviation from the median curve - (work/05_dispersion/pair_figures:
              spread_*, rejected_*, deviating_*), for checking the picks by eye without looking only at good results
"""

import os

import h5py
import numpy as np
from scipy.signal import butter, sosfiltfilt

from .. import disp as D
from .. import signal as S
from ..export import read_curves
from . import style
from .style import plt

XT = [m * 10.0**e for e in range(-3, 4) for m in (1, 2, 5)]  # 1-2-5 ticks for any band


def _freq_axis(ax, f):
    ax.set_xscale("log")
    t = [x for x in XT if f[0] * 0.99 <= x <= f[-1] * 1.01]
    if len(t) < 3:  # a narrow band: every 1-2-3-4-5-6-8 x 10^k inside it
        fine = [m * 10.0**k for k in range(-3, 4) for m in (1, 2, 3, 4, 5, 6, 8)]
        t = [x for x in fine if f[0] * 0.99 <= x <= f[-1] * 1.01]
    ax.set_xticks(t)
    ax.set_xticklabels(["%g" % x for x in t])
    ax.minorticks_off()
    ax.set_xlabel("Frequency (Hz)")
    ax.set_xlim(f[0], f[-1])


def _phase_image(sig, fs, freqs, vel, dist, vg, u_ref, alpha, phi0):
    """Phase-only coherent power, |1 + exp(i residual)|^2 / 4, from the measured signal.

    Rebuild the measurement's phase-matched filter and evaluate its phase against d / c, without using phase
    picks or their cycle counts. This diagnostic is normalised at each frequency, not absolute signal energy
    or an acceptance score. Missing group arrivals remain blank, including pairs without any group picks.
    """
    image = np.full((len(freqs), len(vel)), np.nan)
    measured = np.isfinite(vg) & (vg > 0)
    ref_ok = np.isfinite(u_ref) & (u_ref > 0)
    if not measured.any() or not ref_ok.any():
        return image
    # Match reference_work's nearest-valid filling and process_chunk's log-frequency interpolation.
    idx = np.flatnonzero(ref_ok)
    near = idx[np.argmin(np.abs(np.arange(len(freqs))[:, None] - idx[None, :]), axis=1)]
    vg_fill = np.where(measured, vg, u_ref[near])
    po = int(round(1 / np.log2(freqs[1] / freqs[0])))
    ffine = S.log_freqs(freqs[0], freqs[-1], 4 * po)
    tg_fine = dist / np.interp(np.log(ffine), np.log(freqs), vg_fill)
    tp0 = D.matched_phase(sig, fs, ffine, tg_fine, alpha, phi0)[::4]
    residual = 2 * np.pi * freqs[:, None] * (dist / vel[None, :] - tp0[:, None])
    image[measured] = np.cos(0.5 * residual[measured]) ** 2
    return image


def _images(cfg, disp_path, stack_path, i):
    with h5py.File(disp_path, "r") as h:
        freqs = h["freqs"][:]
        vel = h["vel_axis"][:]
        okc, vc, vg, dist = (h[k][i] for k in ("ok_c", "vc", "vg", "dist_km"))
        alpha = float(h.attrs["alpha"])
        phi0 = float(h.attrs["phase_constant"])
        zone = int(h["zone"][i]) if "zone" in h else -1
        ref = h
        if zone >= 0:
            if "reference_zones" not in h or str(zone) not in h["reference_zones"]:
                raise RuntimeError(
                    "measurement has no saved reference for zone %d; rerun `lant run disp -c CONFIG` "
                    "before redrawing this pair" % zone
                )
            ref = h["reference_zones"][str(zone)]
        u_ref = ref["u_ref"][:].astype(float)
        cref = ref["c_ref"][:].astype(float) if ref["c_ref"].ndim == 1 else ref["c_ref"][i].astype(float)
        if "ref_valid" in ref:
            cref = np.where(ref["ref_valid"][:], cref, np.nan)
        support = ref["ref_support"][:].astype(bool) if "ref_support" in ref else np.ones(len(freqs), bool)
        lam_initial = h.attrs.get("lam_min")
        if lam_initial is None and cfg["project"]["output_dir"]:
            from ..manifest import load

            previous = load(cfg)["steps"].get("disp", {}).get("params", {})
            lam_initial = previous.get("dispersion.lam_min")
        lam_final = h.attrs.get("lam_c_final")
        if lam_initial is None or lam_final is None:
            raise RuntimeError(
                "measurement has no saved wavelength limits; rerun `lant run disp -c CONFIG` before redrawing"
            )
        lam_initial, lam_final = float(lam_initial), float(lam_final)
    with h5py.File(stack_path, "r") as h:
        x = h["stack"][i].astype(float)
        fs, maxlag = float(h.attrs["fs"]), float(h.attrs["maxlag"])
    D.VEL = vel
    _, gimg, _ = D.symmetric_images(x[None, :], np.array([dist]), fs, maxlag, freqs, alpha)
    gimg = gimg[0]
    gimg = (gimg / (gimg.max(axis=1, keepdims=True) + 1e-30)) ** 2
    pos, neg = S.branches(x[None, :], fs, maxlag)
    pimg = _phase_image(0.5 * (pos[0] + neg[0]), fs, freqs, vel, dist, vg, u_ref, alpha, phi0)
    # Cycle branches of the measured phase: c_n = d / (t_p + n / f), t_p = d / vc; the pick is branch n = 0.
    tp = np.where(okc, dist / vc, np.nan)
    branches = {}
    for n in (-3, -2, -1, 1, 2, 3):
        t_n = tp + n / freqs
        branches[n] = np.where(t_n > 0, dist / np.where(t_n > 0, t_n, 1), np.nan)
    reference_full = cref.copy()
    cref = np.where(support, cref, np.nan)  # the reference only where the F-J picks support it
    w = sosfiltfilt(
        butter(4, [freqs[0], min(freqs[-1], 0.45 * fs)], btype="bandpass", fs=fs, output="sos"), x - x.mean()
    )
    c = len(w) // 2
    n = min(int(round(dist / (0.8 * vel[0]) * fs)) + 1, c)
    sym = 0.5 * (w[c : c + n] + w[c::-1][:n])
    sym /= np.abs(sym).max() + 1e-30
    return dict(
        freqs=freqs,
        vel=vel,
        okc=okc,
        vc=vc,
        vg=vg,
        dist=float(dist),
        branches=branches,
        cref=cref,
        reference_full=reference_full,
        pimg=pimg,
        gimg=gimg,
        tt=np.arange(n) / fs,
        sym=sym,
        parameters=dict(
            alpha=alpha, phase_constant=phi0, lam_initial=lam_initial, lam_final=lam_final, reference_zone=zone
        ),
    )


def pair_example(cfg, disp_path, stack_path, pair, out, vwin=None):
    style.apply(cfg["plots"]["font_scale"])
    if vwin is None:
        from .ccf import figure_band

        vwin = figure_band(cfg)[1]
    with h5py.File(disp_path, "r") as h:
        i = [p.decode() for p in h["pairs"][:]].index(pair)
        fr = h["freqs"][:]
    a = _images(cfg, disp_path, stack_path, i)
    cref = a["reference_full"]
    lam_initial, lam_final = a["parameters"]["lam_initial"], a["parameters"]["lam_final"]
    # frequencies where the path is at least one wavelength of the reference (below, the images carry no phase
    # velocity information); and not above the first accepted frequency
    okf = np.isfinite(cref) & (cref > 0) & (a["dist"] * fr / np.where(cref > 0, cref, np.inf) >= 1.0)
    f_lo = fr[okf].min() if okf.any() else a["freqs"][0]
    if a["okc"].any():
        f_lo = min(f_lo, a["freqs"][a["okc"]].min())
    f_lo = fr[max(0, np.searchsorted(fr, f_lo) - 2)]  # show the approach to the wavelength gate
    keepf = a["freqs"] >= f_lo * (1 - 1e-9)
    fig = plt.figure(figsize=(12.0, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 2.7], hspace=0.42, wspace=0.48)
    ax = fig.add_subplot(gs[0, :])
    t1, t2 = a["dist"] / vwin[0], a["dist"] / vwin[1]
    inw = (a["tt"] >= t1) & (a["tt"] <= t2)
    ax.axvspan(min(t1, t2), max(t1, t2), color="0.94", zorder=0)
    ax.plot(a["tt"], a["sym"], color="0.65", lw=0.8)
    ax.plot(a["tt"][inw], a["sym"][inw], color="k", lw=1.2)
    for tv, v in ((t1, vwin[0]), (t2, vwin[1])):
        ax.axvline(tv, color="0.4", ls=":", lw=1.0, label="%.3g km/s" % v)
    ax.set_title("(a) Symmetric correlation", loc="left")
    ax.set_xlim(0, a["tt"][-1])
    ax.set_ylim(-1.1, 1.1)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.set_yticks([-1, 0, 1])
    ax.legend(loc="upper right", frameon=False, ncol=2)
    ax.text(
        0.99,
        0.06,
        "%s   %.2f km" % (pair.replace("_", " – "), a["dist"]),
        transform=ax.transAxes,
        ha="right",
        va="bottom",
    )
    # velocity axis: the measured phase and group velocities of this pair and the reference, with a margin
    vals = np.r_[a["vc"][a["okc"]], a["vg"][a["okc"]], a["cref"][keepf][np.isfinite(a["cref"][keepf])]]
    vlo, vhi = a["vel"][0], a["vel"][-1]
    if vals.size:
        vlo = max(vlo, np.floor(10 * 0.75 * vals.min()) / 10)
        vhi = min(vhi, np.ceil(10 * 1.3 * vals.max()) / 10)
    ax = fig.add_subplot(gs[1, 0])  # phase: measured-signal coherence, cycle branches, reference, picks
    fk = a["freqs"][keepf]
    pm = ax.pcolormesh(
        fk, a["vel"], a["pimg"][keepf].T, shading="gouraud", cmap=style.ENERGY,
        vmin=0, vmax=1, rasterized=True, zorder=0,
    )  # fmt: skip
    cb = fig.colorbar(pm, cax=ax.inset_axes([1.04, 0, 0.04, 1]))
    cb.set_label("Phase coherence")
    cb.set_ticks([0, 0.5, 1])
    if not np.isfinite(a["pimg"][keepf]).any():
        ax.text(0.5, 0.5, "No usable phase-matched signal", transform=ax.transAxes, ha="center", fontsize=8)
    for n, cn in a["branches"].items():
        ax.plot(a["freqs"], cn, "-", color="0.7", lw=0.8, zorder=1)
    ax.plot([], [], "-", color="0.7", lw=0.8, label="Other numbers of cycles")
    ax.plot(fk, a["cref"][keepf], "--", color=style.BLUE, lw=1.4, zorder=2, label="Reference (F-J)")
    rejected = ~a["okc"] & np.isfinite(a["vc"]) & keepf
    ax.plot(a["freqs"][rejected], a["vc"][rejected], "x", color="0.5", ms=4, label="Rejected candidate")
    ax.plot(a["freqs"][a["okc"]], a["vc"][a["okc"]], "o", ms=4, color=style.RED, mec="white", mew=0.4,
            zorder=3, label="Accepted")  # fmt: skip
    path_ok = np.isfinite(cref) & (a["dist"] * fr >= lam_initial * cref)
    ax.fill_between(
        fk,
        0,
        1,
        where=~path_ok[keepf],
        transform=ax.get_xaxis_transform(),
        color="0.9",
        alpha=0.6,
        zorder=0.5,
        label=r"Below $d/\lambda_{ref} = %g$" % lam_initial,
    )
    ax.plot(
        fk,
        a["dist"] * fk / lam_final,
        ":",
        color="k",
        lw=1,
        label=r"Final limit: $d/\lambda = %g$" % lam_final,
    )
    _freq_axis(ax, fk)
    ax.set_ylim(vlo, vhi)
    ax.set_ylabel("Phase velocity (km/s)")
    ax.set_title("(b) Phase velocity", loc="left")
    ax.legend(loc="best", frameon=True, facecolor="white", framealpha=0.85, fontsize="small")
    ax = fig.add_subplot(gs[1, 1])  # group: envelope image and the measured arrival
    pm = ax.pcolormesh(
        fk, a["vel"], a["gimg"][keepf].T, shading="gouraud", cmap=style.ENERGY, vmin=0, vmax=1, rasterized=True
    )
    ax.plot(a["freqs"][a["okc"]], a["vg"][a["okc"]], "o", ms=4, mfc="k", mec="k")
    _freq_axis(ax, fk)
    ax.set_ylim(vlo, vhi)
    ax.set_ylabel("Group velocity (km/s)")
    ax.set_title("(c) Group arrival", loc="left")
    cb = fig.colorbar(pm, cax=ax.inset_axes([1.04, 0, 0.04, 1]))
    cb.set_label("Normalised energy")
    cb.set_ticks([0, 0.5, 1])
    fig.subplots_adjust(bottom=0.14)
    fig.text(
        0.06,
        0.02,
        "Phase: normalised phase-matched coherence; group: normalised envelope energy. Colour is not a QC score.\n"
        "Grey: initial reference gate; dotted: final wavelength limit. Other QC gates apply. Array F-J: 3_fj_phase_velocity.",
        fontsize=8.5,
    )
    return dict(
        files=style.save(fig, out, cfg),
        parameters=dict(
            a["parameters"], waveform_band_hz=[float(fr[0]), float(fr[-1])], velocity_window_kms=list(vwin)
        ),
    )


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
        if hi / lo < 2:  # at least an octave around the accepted band
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
    ax.fill_between(f, mean - std, mean + std, where=full, color=style.RED, alpha=0.12, lw=0)
    ax.errorbar(f[b], mean[b], yerr=std[b], fmt="none", ecolor=style.RED, elinewidth=1.2, capsize=3, zorder=5)
    ax.plot(f[full], mean[full], "-", c=style.RED, lw=2, zorder=6, label="Mean ± std. (at least 20 pairs)")
    _freq_axis(ax, np.array([f.min(), f.max()]))
    ax.set_ylabel("Phase velocity (km/s)")
    sec = ax.secondary_xaxis(
        "top", functions=(lambda x: 1 / np.clip(x, 1e-3, None), lambda x: 1 / np.clip(x, 1e-3, None))
    )
    sec.set_xticks([t for t in XT[::-1] if 1 / f.max() * 0.99 <= t <= 1 / f.min() * 1.01])
    sec.set_xticklabels(["%g" % t for t in sec.get_xticks()])
    sec.minorticks_off()
    sec.set_xlabel("Period (s)")
    ax.legend(loc="upper right", frameon=False)
    ax.text(0.03, 0.05, "Rayleigh-wave phase velocity", transform=ax.transAxes, fontsize=14)
    return style.save(fig, out, cfg)


def _pair_figure(task):
    cfg, disp_path, stack, pair, out, kind = task
    pair_example(cfg, disp_path, stack, pair, out)
    return kind


def pair_figures(cfg, disp_path, out_dir, n):
    """Single-pair figures: n pairs with data spread evenly over distance (the same every run), plus max(3, n // 10)
    with the most rejected candidate points and as many with the largest median deviation of their accepted phase
    velocities from the per-frequency median."""
    import glob

    from .. import paths, resources
    from ..processes import pool_results
    from ..progress import track

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
    tasks = [
        (sub, disp_path, stack, pairs[i], os.path.join(out_dir, "%s_%s" % (kind, pairs[i])), kind)
        for i, kind in pick.items()
    ]
    done = {}
    if not tasks:
        return dict(figures=done, directory=out_dir)
    with h5py.File(disp_path, "r") as h:
        nf, nv = len(h["freqs"]), len(h["vel_axis"])
    estimate = D.memory_estimate(stack, 1, nf, nv, 0)
    plan = resources.plan(cfg, len(tasks), estimate, "pair diagnostic figures")
    plan["threads_each"] = 1  # rendering uses one numerical thread per process
    if plan["parallel"] == 1:
        for kind in track(map(_pair_figure, tasks), len(tasks), "pair figures"):
            done[kind] = done.get(kind, 0) + 1
    else:
        from multiprocessing import get_context

        # Matplotlib state is private; spawning is safe after threaded numerical libraries have been used.
        with get_context("spawn").Pool(plan["parallel"], initializer=resources.limit_threads, initargs=(1,)) as pool:
            for kind in track(pool_results(pool, _pair_figure, tasks), len(tasks), "pair figures"):
                done[kind] = done.get(kind, 0) + 1
            resources.note_workers(w.pid for w in pool._pool)
    return dict(figures=done, directory=out_dir)
