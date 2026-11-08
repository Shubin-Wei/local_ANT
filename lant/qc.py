"""Quality control after the measurement: group velocity rules, wavelength limits, redundant-path precision.

Group velocity (dispersion.group_qc; kept as a product, not inverted): SNR >= snr_g, distinct side lobe <= side_g,
jackknife <= jk_g, deviation from the reference group velocity <= dev_g, path >= lam_g PHASE wavelengths
(Bensen et al. 2007 sect. 5), and either the causal / acausal agreement (<= branch_diff, weaker branch SNR >=
snr_branch_min; areal arrays) or the jump between neighbouring frequencies (<= step_g; lines, where the noise is
too directional for the branches to agree).  Contiguous runs of >= min_run samples.
Wavelength limits (both velocities, from each sample's own phase velocity): group >= lam_g, phase >= qc.lam_c.
"""
import shutil

import h5py
import numpy as np
from scipy.spatial import cKDTree

from . import disp
from .signal import contiguous, local_xy_km


def group_qc(cfg, path_in, path_out, first_pass=False):
    g = cfg["dispersion"]["group_qc"]
    min_run = cfg["dispersion"]["qc"]["min_run"]
    shutil.copyfile(path_in, path_out)
    with h5py.File(path_out, "r+") as h:
        f = h["freqs"][:]
        dist = h["dist_km"][:]
        vg, snr, side, jk = h["vg"][:], h["snr"][:], h["side"][:], h["jk_g"][:]
        U = h["u_ref"][:][None, :]
        dev = np.abs(vg - U) / U
        with np.errstate(invalid="ignore", divide="ignore"):
            vc = h["vc"][:]
            lam = np.where(np.isfinite(vc) & (vc > 0), vc, h["c_ref"][:][None, :]) / f[None, :]
        side_off = disp.off_band_mask(f, cfg["dispersion"]["side_off_bands"])[None, :]
        with np.errstate(invalid="ignore"):
            ok = (np.isfinite(vg) & (snr >= g["snr_g"]) & ((side <= g["side_g"]) | side_off) & (jk <= g["jk_g"] * vg)
                  & np.isfinite(U) & (dist[:, None] >= g["lam_g"] * lam))
            if g["dev_g"] > 0 and not first_pass:
                ok &= dev <= g["dev_g"]
            if g["branch_diff"] > 0:
                vp, vn = h["vg_pos"][:], h["vg_neg"][:]
                bdif = 2 * np.abs(vp - vn) / (vp + vn)
                ok &= np.isfinite(bdif) & (bdif <= g["branch_diff"])
                if g["snr_branch_min"] > 0:
                    weak = np.minimum(h["snr_pos"][:], h["snr_neg"][:])
                    ok &= np.isfinite(weak) & (weak >= g["snr_branch_min"])
            if g["step_g"] > 0:
                step = np.abs(np.diff(np.log(vg), axis=1))
                bad = np.isfinite(step) & (step > g["step_g"]) & ok[:, 1:] & ok[:, :-1]
                ok[:, 1:] &= ~bad
                ok[:, :-1] &= ~bad
        for q in range(ok.shape[0]):
            if ok[q].any():
                ok[q] = contiguous(ok[q], min_run)
        n0 = int(h["ok_g"][:].sum())
        h["ok_g"][...] = ok
        h.attrs.update({"group_qc_" + k: v for k, v in g.items()})
    return dict(group_points_before=n0, group_points=int(ok.sum()))


def lambda_qc(cfg, path_in, path_out):
    lam_gl, lam_cl = cfg["dispersion"]["group_qc"]["lam_g"], cfg["dispersion"]["qc"]["lam_c"]
    min_run = cfg["dispersion"]["qc"]["min_run"]
    shutil.copyfile(path_in, path_out)
    out = {}
    with h5py.File(path_out, "r+") as h:
        f = h["freqs"][:]
        dist = h["dist_km"][:]
        with np.errstate(invalid="ignore", divide="ignore"):
            lam_c = h["vc"][:] / f[None, :]
            lam_g = np.where(np.isfinite(lam_c) & (lam_c > 0), lam_c, h["c_ref"][:][None, :] / f[None, :])
        for key, lam, limit in (("ok_g", lam_g, lam_gl), ("ok_c", lam_c, lam_cl)):
            ok = h[key][:]
            n0 = int(ok.sum())
            with np.errstate(invalid="ignore"):
                ok &= np.isfinite(lam) & (dist[:, None] >= limit * lam)
            for q in range(ok.shape[0]):
                if ok[q].any():
                    ok[q] = contiguous(ok[q], min_run)
            h[key][...] = ok
            out[key] = (n0, int(ok.sum()))
        h.attrs["qc_lam_g"] = lam_gl
        h.attrs["qc_lam_c"] = lam_cl
    return dict(group_points=out["ok_g"][1], phase_points=out["ok_c"][1])


def robust_std(x):
    return 1.4826 * np.median(np.abs(x - np.median(x)))


BANDS = [(0.5, 1.0), (1.0, 2.0), (2.0, 4.0), (4.0, 8.0), (8.0, 10.5)]


def redundant_precision(disp_path, stack_path, radius=0.5, min_dist=3.0):
    """Paths whose two ends are each within `radius` km of the other path's ends sample almost the same medium:
    the robust std of their relative velocity difference per band bounds the measurement error."""
    with h5py.File(stack_path, "r") as h:
        lon1, lat1, lon2, lat2 = (h[k][:] for k in ("lon1", "lat1", "lon2", "lat2"))
        dist = h["dist_km"][:]
    lat0, lon0 = float(np.mean(np.r_[lat1, lat2])), float(np.mean(np.r_[lon1, lon2]))
    p1 = np.c_[local_xy_km(lon1, lat1, lat0, lon0)]
    p2 = np.c_[local_xy_km(lon2, lat2, lat0, lon0)]
    with h5py.File(disp_path, "r") as h:
        freqs = h["freqs"][:]
        vg, vc, okg, okc = h["vg"][:], h["vc"][:], h["ok_g"][:], h["ok_c"][:]
    idx = np.where(dist >= min_dist)[0]
    feat = np.c_[p1[idx], p2[idx]]
    pairs = set()
    if len(idx) > 1:
        tree = cKDTree(feat)
        for f in (feat, np.c_[p2[idx], p1[idx]]):
            for a, nb in enumerate(tree.query_ball_point(f, r=radius * np.sqrt(2))):
                for b in nb:
                    if b == a:
                        continue
                    if np.hypot(*(f[a, :2] - feat[b, :2])) <= radius and np.hypot(*(f[a, 2:] - feat[b, 2:])) <= radius:
                        i, j = idx[a], idx[b]
                        pairs.add((min(i, j), max(i, j)))
    rp = sorted(pairs)
    res = {"n_redundant_pairs": len(rp), "group": {}, "phase": {}}
    lines = ["redundant path pairs (ends within %.2f km, distance >= %.1f km): %d" % (radius, min_dist, len(rp))]
    for name, V, OK in (("group", vg, okg), ("phase", vc, okc)):
        lines.append("%s velocity: relative difference between redundant paths" % name)
        for lo, hi in BANDS:
            k = (freqs >= lo) & (freqs < hi)
            diffs = []
            for i, j in rp:
                m = OK[i, k] & OK[j, k]
                if m.any():
                    vi, vj = V[i, k][m], V[j, k][m]
                    diffs.append(np.mean(2 * (vi - vj) / (vi + vj)))
            diffs = np.array(diffs)
            s = robust_std(diffs) if len(diffs) > 5 else np.nan
            res[name]["%.1f-%.1f" % (lo, hi)] = dict(n=len(diffs), robust_std=s)
            lines.append("  %4.1f-%4.1f Hz: n=%4d robust std %.3f" % (lo, hi, len(diffs), s))
    return res, "\n".join(lines) + "\n"
