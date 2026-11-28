"""The reference phase-velocity curve of the measurement (disp.reference).

It does three things in lant.dispersion: it sets where the group arrival is looked for, it fixes the cycle count
of the phase, and the frequencies where it is defined are the frequencies measured.  It does not enter the measured
values.  Three sources, the common ones in the literature:
  fj     the F-J spectrum of this array (lant.fj; Wang et al. 2019): the largest value at each frequency of
         disp.fj_band_hz, a cubic in ln f through these picks; defined over the band of the picks;
  model  a 1-D Vs node model (disp.reference_file: depth_km Vs_km/s): its fundamental Rayleigh phase velocities
         (Bensen et al. 2007 and Lin et al. 2008 use a model's prediction); defined at every frequency;
  file   a curve (disp.reference_file: frequency_Hz phase_velocity_km/s), e.g. picked by hand from the F-J figure:
         interpolated in ln f between its points, not extrapolated.
The F-J spectrum is always computed and drawn, so that any reference can be checked against the array's own
spectrum.  The group velocity of the reference is U = c / (1 - d ln c / d ln f) (Bensen et al. 2007 eq. 7).
"""

import os

import numpy as np

from . import fj


def group_velocity(freqs, c):
    """U = c / (1 - d ln c / d ln f), inside each run of defined values (NaN elsewhere)."""
    from .signal import runs_of

    u = np.full(len(c), np.nan)
    for k0, k1 in runs_of(np.isfinite(c)):
        if k1 - k0 >= 2:
            u[k0:k1] = c[k0:k1] / (1.0 - np.gradient(np.log(c[k0:k1]), np.log(freqs[k0:k1])))
    return u


def from_model(path, freqs):
    from .inversion.start_model import phase_velocity, read_node_model

    z, vs = read_node_model(path)
    c = np.full(len(freqs), np.nan)
    for i0 in range(0, len(freqs), 60):  # surfdisp96 takes at most 80 periods at a time
        c[i0 : i0 + 60] = phase_velocity(z, vs, 1.0 / freqs[i0 : i0 + 60])
    c[c <= 0] = np.nan
    return c


def from_file(path, freqs):
    a = np.loadtxt(path, ndmin=2, comments="#")
    if a.shape[1] < 2 or len(a) < 2:
        raise ValueError("%s: two columns frequency_Hz phase_velocity_km/s, at least two rows" % path)
    a = a[np.argsort(a[:, 0])]
    if np.any(np.diff(a[:, 0]) <= 0) or np.any(a[:, :2] <= 0):
        raise ValueError("%s: frequencies must be distinct and positive, velocities positive" % path)
    c = np.exp(np.interp(np.log(freqs), np.log(a[:, 0]), np.log(a[:, 1])))
    c[(freqs < a[0, 0] * (1 - 1e-9)) | (freqs > a[-1, 0] * (1 + 1e-9))] = np.nan
    return c


def scan_velocities(vrange):
    """F-J velocity axis: half the slowest to 1.5 x the fastest of data.velocity_range_kms, 400 steps or finer."""
    lo, hi = 0.5 * float(vrange[0]), 1.5 * float(vrange[1])
    dc = float("%.1g" % ((hi - lo) / 400))
    return np.arange(lo, hi + dc / 2, dc)


def build(cfg, ccf, dist, fs, freqs, log=print):
    """-> dict(freqs, c, u, fj_image, fj_c, fj_picks, source)."""
    d = cfg["disp"]
    c_scan = scan_velocities(cfg["data"]["velocity_range_kms"])
    log("F-J spectrum: %d pairs, %d frequencies, velocities %.3g-%.3g km/s" % (len(dist), len(freqs), c_scan[0],
                                                                             c_scan[-1]))  # fmt: skip
    img = fj.image(ccf, dist, fs, freqs, c_scan)
    band = d["fj_band_hz"]
    picks = fj.pick(img, c_scan, freqs, band)
    vlo, vhi = (float(v) for v in cfg["data"]["velocity_range_kms"])
    picks[(picks < vlo) | (picks > vhi)] = np.nan  # outside the range expected in the area: not a pick
    if d["reference"] == "fj":
        c = fj.smooth_curve(freqs, picks)
        source = "F-J picks %s" % ("in %g-%g Hz" % tuple(band) if band else "at every frequency")
    elif d["reference"] == "model":
        c = from_model(d["reference_file"], freqs)
        source = "1-D model %s" % os.path.basename(d["reference_file"])
    else:
        c = from_file(d["reference_file"], freqs)
        source = "curve %s" % os.path.basename(d["reference_file"])
    if np.isfinite(c).sum() < 2:
        raise RuntimeError("the reference curve (%s) covers fewer than 2 measurement frequencies" % source)
    return dict(freqs=freqs, c=c, u=group_velocity(freqs, c), fj_image=img, fj_c=c_scan, fj_picks=picks,
                source=source)  # fmt: skip


def write_csv(path, ref):
    with open(path, "w") as fh:
        fh.write("# reference curve: %s\nfrequency_hz,period_s,c_ref_kms,u_ref_kms,fj_pick_kms\n" % ref["source"])
        for f, c, u, p in zip(ref["freqs"], ref["c"], ref["u"], ref["fj_picks"]):
            fh.write("%.6f,%.6f,%s,%s,%s\n" % (f, 1 / f, _v(c), _v(u), _v(p)))


def _v(x):
    return "" if not np.isfinite(x) else "%.5f" % x
