"""Quality control after the measurement, and a report of the measurement against its reference.

Wavelength limit (phase velocity, from each sample's own accepted phase velocity): path >= qc.lam_c wavelengths
(Bensen et al. 2007 use 3, Luo et al. 2015 show 1 can be enough; 2 here), then contiguous runs of >= qc.min_run
samples.  The group velocity is measured to find the arrival the phase is read at; it is not a product (no QC of
its own, not exported).
The median of the accepted phase velocities against the reference, per frequency, is reported
(median_vs_reference.txt): a systematic offset means the reference (and so the cycle count it fixes) may be off.
"""

import h5py
import numpy as np

from .signal import contiguous


def lambda_qc(cfg, path):
    """Applied in place to the measurement file; the samples it removes keep their reason bits."""
    lam_cl = cfg["dispersion"]["qc"]["lam_c"]
    min_run = cfg["dispersion"]["qc"]["min_run"]
    with h5py.File(path, "r+") as h:
        f = h["freqs"][:]
        dist = h["dist_km"][:]
        with np.errstate(invalid="ignore", divide="ignore"):
            lam_c = h["vc"][:] / f[None, :]
        ok = h["ok_c"][:]
        n0 = int(ok.sum())
        with np.errstate(invalid="ignore"):
            short = ok & ~(np.isfinite(lam_c) & (dist[:, None] >= lam_cl * lam_c))
        ok &= ~short
        for q in range(ok.shape[0]):
            if ok[q].any():
                ok[q] = contiguous(ok[q], min_run)
        dropped_run = (h["ok_c"][:] & ~short) & ~ok
        h["ok_c"][...] = ok
        if "reason" in h:
            r = h["reason"][:]
            r[short] |= REASONS["wavelength"]
            r[dropped_run] |= REASONS["short_run"]
            h["reason"][...] = r
        h.attrs["qc_lam_c"] = lam_cl
    return dict(phase_points_before=n0, phase_points=int(ok.sum()))


# reason codes of a rejected phase-velocity sample (bit mask in the "reason" dataset of the measurement files)
# (every sample that is not accepted carries at least one bit; 0 = accepted; bits 1024 and 2048 belonged to two
# criteria removed in version 1.4 and are not used)
REASONS = {
    "no_group_pick": 1,
    "snr_group": 2,
    "no_data": 4,
    "dev_pre": 8,
    "path_short_measure": 16,
    "side_lobe": 32,
    "no_anchor": 64,
    "phase_slip": 128,
    "jackknife": 512,
    "group_phase": 4096,
    "short_run": 8192,
    "wavelength": 16384,
    "no_reference": 32768,
    "off_band": 65536,
}


def reason_counts(path):
    """{reason: samples with that reason set} of a measurement file."""
    with h5py.File(path, "r") as h:
        if "reason" not in h:
            return {}
        r = h["reason"][:]
    return {k: int(((r & v) > 0).sum()) for k, v in REASONS.items()}


def median_vs_reference(path, out_txt):
    with h5py.File(path, "r") as h:
        f = h["freqs"][:]
        ok = h["ok_c"][:]
        vc = np.where(ok, h["vc"][:], np.nan)
        cref = h["c_ref"][:]
    n = ok.sum(0)
    import warnings

    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(vc, axis=0)
        dev = 100 * (med / cref - 1)
    L = ["# per frequency: accepted pairs, median phase velocity, reference, median - reference (%)"]
    for i in range(len(f)):
        L.append("%8.4f %6d %8.4f %8.4f %+7.2f" % (f[i], n[i], med[i], cref[i], dev[i]))
    open(out_txt, "w").write("\n".join(L) + "\n")
    m = (n >= 10) & np.isfinite(dev)
    return dict(
        frequencies=int(m.sum()),
        median_percent=float(np.median(dev[m])) if m.any() else None,
        max_abs_percent=float(np.max(np.abs(dev[m]))) if m.any() else None,
    )
