"""Step 05: reference phase-velocity curve.

The reference guides the group ridge, fixes the cycle count of the phase and centres the deviation screens.
  reference.source = fj (default): F-J spectrum of this array with CC-FJpy fj_noise (frequency-Bessel transform,
      Wang et al. 2019; Li et al. 2021), one frequency at a time with the pairs not longer than fj_rmax_lambda
      wavelengths of fj_c_nominal (long paths stop being coherent at high frequency and the image breaks up);
        fj_mode max_peak (areal arrays): the largest peak (assumes the fundamental mode is strongest), picks in
            notched lines dropped, cubic in log f up to fj_fit_fmax, held constant outside the fitted range;
        fj_mode track_fundamental (lines, where higher modes or off-line noise can dominate): the fundamental
            ridge is followed from an unambiguous start in track_start_band_hz (steps <= 5 % per 1/12 octave, at
            most 3 % rise with frequency), cubic in log f, made non-increasing, constant outside.
      This first-pass curve is used without deviation screens; the SECOND pass reference is the per-frequency
      median of the first-pass phase velocities (>= pass2_min_pairs pairs; gaps <= 0.25 octave interpolated
      in log f; no reference elsewhere, nothing measured there).  Written by lant.pipeline after the first pass.
  reference.source = file: a curve of another array (npz written by lant, or text "freq_Hz c_km/s"), interpolated
      in log f onto the measurement frequencies; outside its range there is no reference (no extrapolation).
      Used directly as the measurement reference (no first pass).
"""
import os
import sys

import h5py
import numpy as np

from . import paths
from .signal import log_freqs


def freq_grid(cfg):
    d = cfg["dispersion"]
    f = log_freqs(d["fmin_hz"], d["fmax_hz"], d["per_octave"])
    return f[f <= d["fmax_hz"] + 1e-9]


def import_ccfj():
    try:
        import ccfj  # noqa: F401
        return ccfj
    except ImportError:
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        sys.path.insert(0, os.path.join(here, "external", "CC-FJpy"))
        try:
            import ccfj
            return ccfj
        except ImportError as exc:
            raise RuntimeError("CC-FJpy is not installed (needed for reference.source = fj): run "
                               "`bash tools/install_external.sh ccfj`, or give a reference curve with "
                               "reference.source = file (%s)" % exc)


def fj_image(cfg, stack_path, table):
    r = cfg["reference"]
    ccfj = import_ccfj()
    with h5py.File(stack_path, "r") as h:
        x = h["stack"][:].astype(np.float64)
        dist = h["dist_km"][:]
        nseg = h["nseg"][:]
        fs, maxlag = float(h.attrs["fs"]), float(h.attrs["maxlag"])
    use = nseg > 0
    x, dist = x[use], dist[use]
    L = int(round(maxlag * fs))
    xc = np.zeros((x.shape[0], r["fj_nfft"]))
    xc[:, :L + 1] = x[:, L:]
    xc[:, -L:] = x[:, :L]
    X = np.real(np.fft.rfft(xc, axis=1))
    fall = np.fft.rfftfreq(r["fj_nfft"], 1 / fs)
    freqs = freq_grid(cfg)
    kbin = np.array([int(np.argmin(np.abs(fall - f))) for f in freqs])
    c = np.arange(r["fj_cmin"], r["fj_cmax"] + 1e-9, r["fj_dc"])
    img = np.zeros((len(freqs), len(c)))
    npair = np.zeros(len(freqs), int)
    for i, (f, kb) in enumerate(zip(freqs, kbin)):
        m = dist * f / r["fj_c_nominal"] <= r["fj_rmax_lambda"]
        npair[i] = m.sum()
        if m.sum() < 30:
            continue
        out = ccfj.fj_noise(X[m][:, kb:kb + 1], dist[m] * 1e3, c * 1e3, np.array([fall[kb]]), fstride=1, itype=1,
                            func=0, num=int(r["fj_threads"]))
        img[i] = out[:, 0]
    return freqs, c, img, npair


def _peaks(y, thr):
    return np.where((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]) & (y[1:-1] >= thr))[0] + 1


def pick_max(cfg, freqs, c, img, table):
    r = cfg["reference"]
    cfj = np.full(len(freqs), np.nan)
    for i in range(len(freqs)):
        if not img[i].any():
            continue
        j = int(np.argmax(img[i]))
        if 0 < j < len(c) - 1:
            cfj[i] = c[j]
    good = np.isfinite(cfj) & (freqs <= r["fj_fit_fmax"])
    for f0, hw in table:
        good &= np.abs(freqs - f0) > hw
    if good.sum() < 5:
        raise RuntimeError("the F-J image has only %d usable picks below %.1f Hz: check the stack (too few pairs?) or "
                           "borrow a reference (reference.source = file)" % (good.sum(), r["fj_fit_fmax"]))
    x = np.log(freqs)
    p = np.polyfit(x[good], np.log(cfj[good]), 3)
    xr = np.minimum(np.log(freqs), np.log(r["fj_fit_fmax"]))
    curve = np.exp(np.polyval(p, np.maximum(xr, x[good].min())))
    return cfj, good, curve


def pick_track(cfg, freqs, c, img, table, min_peak=0.5, min_track=0.3, step=0.05, max_miss=3, cmin=2.0, rise=0.03):
    s0, s1 = cfg["reference"]["track_start_band_hz"]
    Y = img / (img.max(axis=1, keepdims=True) + 1e-30)
    okc = (c >= cmin) & (c <= 0.97 * c.max())
    start = None
    for strict in (True, False):
        for i in np.where((freqs >= s0) & (freqs <= s1))[0]:
            k = _peaks(Y[i], min_peak)
            k = k[okc[k]]
            if k.size and np.argmax(Y[i][k]) == 0 and (k.size == 1 or not strict):
                start = i
                break
        if start is not None:
            break
    if start is None:
        raise RuntimeError("F-J image: no unambiguous fundamental-mode start in %.2f-%.2f Hz "
                           "(reference.track_start_band_hz); borrow a reference (reference.source = file)" % (s0, s1))
    pick = np.full(len(freqs), np.nan)
    k = _peaks(Y[start], min_peak)
    pick[start] = c[k[okc[k]][0]]
    for direction in (1, -1):
        last, miss, i = pick[start], 0, start + direction
        lowest = highest = last
        while 0 <= i < len(freqs) and miss < max_miss:
            k = _peaks(Y[i], min_track)
            k = k[okc[k]]
            if k.size:
                j = k[np.argmin(np.abs(c[k] - last))]
                rel = c[j] / last - 1
                mono = (c[j] <= lowest * (1 + rise)) if direction > 0 else (c[j] >= highest * (1 - rise))
                if abs(rel) <= step and mono:
                    pick[i] = c[j]
                    last = c[j]
                    miss = 0
                    lowest, highest = min(lowest, last), max(highest, last)
                else:
                    miss += 1
            else:
                miss += 1
            i += direction
    use = np.isfinite(pick)
    for f0, hw in table:
        use &= np.abs(freqs - f0) > hw
    x = np.log(freqs)
    pf = np.polyfit(x[use], np.log(pick[use]), 3)
    lo, hi = x[use].min(), x[use].max()
    curve = np.minimum.accumulate(np.exp(np.polyval(pf, np.clip(x, lo, hi))))
    return pick, use, curve


def from_file(cfg, path):
    freqs = freq_grid(cfg)
    if path.endswith(".npz"):
        z = np.load(path)
        f, c = np.asarray(z["freqs"], float), np.asarray(z["c_ref"], float)
        if "ref_valid" in z.files:
            c = np.where(z["ref_valid"], c, np.nan)
    else:
        a = np.loadtxt(path, ndmin=2)
        f, c = a[:, 0], a[:, 1]
        if np.median(f) > 0 and np.all(np.diff(f) < 0):
            f, c = f[::-1], c[::-1]
    m = np.isfinite(c) & (c > 0)
    f, c = f[m], c[m]
    if len(f) < 2 or np.any(np.diff(f) <= 0):
        raise RuntimeError("%s: frequencies must be increasing (first column frequency in Hz, second c in km/s)" % path)
    if c.min() < 0.05 or c.max() > 10:
        raise RuntimeError("%s: velocities %.3f-%.3f look wrong (km/s expected)" % (path, c.min(), c.max()))
    inside = (freqs >= f.min() * (1 - 1e-6)) & (freqs <= f.max() * (1 + 1e-6))
    cr = np.full(len(freqs), np.nan)
    cr[inside] = np.exp(np.interp(np.log(freqs[inside]), np.log(f), np.log(c)))
    return freqs, cr, inside


def pass2(cfg, ref_in, disp_path, out, report):
    """Second-pass reference: per-frequency median of the accepted first-pass phase velocities."""
    r = cfg["reference"]
    ref = dict(np.load(ref_in))
    f = ref["freqs"]
    c1 = np.asarray(ref["c_ref"], float)
    with h5py.File(disp_path, "r") as h:
        assert np.allclose(h["freqs"][:], f)
        ok = h["ok_c"][:]
        vc = np.where(ok, h["vc"][:], np.nan)
    cnt = ok.sum(axis=0)
    with np.errstate(all="ignore"):
        med = np.nanmedian(vc, axis=0)
        q16, q84 = np.nanpercentile(vc, 16, axis=0), np.nanpercentile(vc, 84, axis=0)
    have = cnt >= r["pass2_min_pairs"]
    if not have.any():
        raise RuntimeError("second-pass reference: no frequency has >= %d first-pass phase velocities (max %d); lower "
                           "reference.pass2_min_pairs or check the first pass" % (r["pass2_min_pairs"], cnt.max()))
    lf = np.log2(f)
    step = np.median(np.diff(lf))
    kind = np.where(have, 3, 0)
    c = np.where(have, med, np.nan)
    idx = np.where(have)[0]
    for i0, i1 in zip(idx[:-1], idx[1:]):
        nmiss = i1 - i0 - 1
        if nmiss > 0 and nmiss * step <= r["pass2_max_gap_octave"] + 1e-9:
            kind[i0 + 1:i1] = 4
    gap = kind == 4
    c[gap] = np.exp(np.interp(lf[gap], lf[have], np.log(med[have])))
    valid = kind > 0
    out_ref = dict(freqs=f, c_ref=np.where(valid, c, np.nan), ref_valid=valid, ref_kind=kind, median_pass1=med,
                   npairs_pass1=cnt, c_ref_pass1=c1,
                   c_ref_source=np.array("second pass: median of the first-pass phase velocities (>= %d pairs)" % r["pass2_min_pairs"]))
    np.savez(out, **out_ref)
    d = med / c1 - 1
    L = ["second-pass reference: median of the first-pass phase velocities (>= %d pairs, gaps <= %.2f octave interpolated)"
         % (r["pass2_min_pairs"], r["pass2_max_gap_octave"]),
         "  defined %.2f-%.2f Hz; median against the first-pass reference %+.2f %% to %+.2f %%"
         % (f[valid].min(), f[valid].max(), 100 * np.nanmin(d[have]), 100 * np.nanmax(d[have])), "",
         "  f(Hz)   kind          c_ref   pass1_ref   median  pairs"]
    names = {0: "none", 3: "median", 4: "interpolated"}
    for k in range(len(f)):
        L.append("  %6.3f  %-12s %6.3f   %6.3f     %6.3f  %5d" % (f[k], names[kind[k]], out_ref["c_ref"][k], c1[k], med[k], cnt[k]))
    open(report + ".txt", "w").write("\n".join(L) + "\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    ax = axs[0]
    ax.fill_between(f, q16, q84, color="0.85", label="first pass, 16-84 % of the pairs")
    ax.semilogx(f, med, ".", color="0.4", label="first-pass median")
    ax.semilogx(f, c1, "b--", lw=1, label="first-pass reference")
    ax.semilogx(f, out_ref["c_ref"], "r-", lw=1.8, label="second-pass reference")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Phase velocity (km/s)")
    ax.legend(fontsize=8, frameon=False)
    axs[1].semilogx(f, cnt, "k.-")
    axs[1].axhline(r["pass2_min_pairs"], color="r", lw=0.8)
    axs[1].set_yscale("log")
    axs[1].set_xlabel("Frequency (Hz)")
    axs[1].set_ylabel("Pairs with an accepted first-pass phase velocity")
    fig.tight_layout()
    fig.savefig(report + ".png", dpi=140)
    plt.close(fig)
    return dict(defined_hz=[float(f[valid].min()), float(f[valid].max())], frequencies_without_reference=int((~valid).sum()),
                median_vs_pass1_percent=[float(100 * np.nanmin(d[have])), float(100 * np.nanmax(d[have]))])


def run(cfg):
    from . import notch
    d = paths.mk(paths.P(cfg, paths.REFERENCE))
    r = cfg["reference"]
    if r["source"] == "file":
        freqs, c, valid = from_file(cfg, r["file"])
        np.savez(os.path.join(d, "reference.npz"), freqs=freqs, c_ref=c, ref_valid=valid,
                 c_ref_source=np.array("borrowed: %s" % r["file"]))
        return dict(source="file", file=r["file"], defined_hz=[float(freqs[valid].min()), float(freqs[valid].max())])
    table = notch.read_table(cfg).tolist()
    stack = paths.stack_product(cfg)
    freqs, c, img, npair = fj_image(cfg, stack, table)
    if r["fj_mode"] == "max_peak":
        picks, used, curve = pick_max(cfg, freqs, c, img, table)
    else:
        picks, used, curve = pick_track(cfg, freqs, c, img, table)
    np.savez(os.path.join(d, "fj.npz"), freqs=freqs, c=c, img=img, npair=npair, picks=picks, used=used)
    np.savez(os.path.join(d, "reference_pass1.npz"), freqs=freqs, c_ref=curve, ref_valid=np.ones(len(freqs), bool),
             c_ref_source=np.array("F-J (%s), rmax %.1f lambda" % (r["fj_mode"], r["fj_rmax_lambda"])))
    L = ["F-J first-pass reference (%s), pairs within %.1f wavelengths of %.2f km/s" % (r["fj_mode"], r["fj_rmax_lambda"], r["fj_c_nominal"]),
         "  f(Hz)  pairs   pick    used  curve"]
    for i in range(len(freqs)):
        L.append("  %5.2f  %5d  %s  %s  %6.3f" % (freqs[i], npair[i], "%6.3f" % picks[i] if np.isfinite(picks[i]) else "   -  ",
                                                 "y" if used[i] else "-", curve[i]))
    open(os.path.join(d, "fj.txt"), "w").write("\n".join(L) + "\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.5))
    Y = img / (img.max(axis=1, keepdims=True) + 1e-30)
    ax.pcolormesh(freqs, c, Y.T, cmap="jet", vmin=0, vmax=1, shading="nearest")
    ax.plot(freqs[used], picks[used], "o", ms=3.5, mfc="w", mec="k", mew=0.6, label="picks used")
    ax.plot(freqs, curve, "w-", lw=1.8, label="first-pass reference")
    ax.set_xscale("log")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Phase velocity (km/s)")
    ax.set_title("F-J spectrum (CC-FJpy), normalised per frequency", fontsize=10)
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.join(d, "fj.png"), dpi=150)
    plt.close(fig)
    return dict(source="fj", mode=r["fj_mode"], picks_used=int(used.sum()),
                pairs_per_frequency=[int(n) for n in npair])
