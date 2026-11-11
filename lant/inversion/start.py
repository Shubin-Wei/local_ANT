"""Average curve, depth nodes and the 1-D starting model.

  average curve: per measured period the median over the pairs (robust average), kept where at least
      inversion.median_min_pairs pairs have a value; used as measured (not smoothed, not forced monotonic)
      (Jing et al. 2025: the average dispersion as the reference / starting model);
  depth nodes: inversion.depth_sequence_km, or "auto": from 0 in steps of max(a quarter of the shortest
      wavelength of the average curve, a quarter of the depth) - the shortest wavelength limits the thinnest layer
      that can be resolved and the vertical resolution of surface waves decreases with depth.  The nodes not deeper
      than half the longest wavelength are updated, the next node below is the fixed bottom node (DSurfTomo
      neither updates nor writes the last node; its velocity is a boundary assumption).  "Half the longest
      wavelength" is the empirical investigation depth of Foti et al. (2018, sect. 4.2.3).  Then a node every
      fine_step_km is added down to fine_zmax_km (linear interpolation of the 1-D model there);
  1-D model: damped and smoothed Gauss-Newton on ln Vs with a numerical Jacobian, forward = DSurfTomo's own chain
      (Brocher 2005 Vp and density, refineGrid2LayerMdl, surfdisp96) through the small program fwd1d; start
      1.087 x the phase velocity at the shortest period, smoothing of the curvature of the model itself
      ("jumping", Shaw & Orcutt 1985; Constable et al. 1987: not of the update, so the result does not drift with
      the number of iterations) on the actual node spacing (curvature_operator; weight `smooth`, independent of the
      number and spacing of the nodes), damping `damp` of the step, `iterations` steps limited to +-15 %.  DSurfTomo's own GenerateIniMOD.py sets Vs equal to the phase
      velocity (about 9 % slow) and is not used.
"""
import os
import subprocess
import tempfile

import numpy as np

from .. import export
from ..tools import external_path


def average_curve(curve_dir, tmin, tmax, min_pairs):
    """Median phase velocity per period over the pairs.  A period needs min_pairs pairs, but not more than a quarter
    of the pairs that have a curve (small arrays), and at least 5."""
    acc = {}
    curves = list(export.read_curves(curve_dir, "C", tmin, tmax))
    min_pairs = int(max(5, min(min_pairs, len(curves) // 4)))
    for c in curves:
        for t, v in zip(c["T"], c["v"]):
            acc.setdefault(round(float(t), 6), []).append(v)
    rows = [(t, float(np.median(v)), len(v)) for t, v in sorted(acc.items()) if len(v) >= min_pairs]
    if len(rows) < 3:
        raise RuntimeError("average curve: only %d periods have >= %d pairs (inversion.median_min_pairs, at most a "
                           "quarter of the %d pairs with a curve): too few measurements for a 1-D model"
                           % (len(rows), min_pairs, len(curves)))
    return np.array(rows)


def auto_sequence(avg):
    """Depth nodes from 0: steps of max(a quarter of the shortest wavelength, a quarter of the depth), rounded to two
    significant digits (and whole metres: DSurfTomo's model file keeps 3 decimals of km), one node past half the
    longest wavelength."""
    lam = avg[:, 0] * avg[:, 1]
    dz0, half = float(lam.min()) / 4, float(lam.max()) / 2
    z = [0.0]
    while z[-1] <= half + 1e-9:
        nxt = z[-1] + max(dz0, 0.25 * z[-1])
        nxt = max(round(float("%.2g" % nxt), 3), round(z[-1] + 0.001, 3))   # 2 significant digits, metres at least
        z.append(nxt)
    return z


def depth_nodes(avg, seq):
    """Updatable nodes (not deeper than half the longest wavelength) + the next node as the fixed bottom node.  When
    the sequence is not deep enough it is extended, each step 1.3 times the previous one (reported)."""
    i = int(np.argmax(avg[:, 0]))
    half = avg[i, 0] * avg[i, 1] / 2
    rule = "given"
    if seq == "auto":
        seq, rule = auto_sequence(avg), "auto: steps max(shortest wavelength / 4, depth / 4)"
    seq = [round(float(v), 3) for v in seq]
    if seq[0] != 0 or any(b <= a for a, b in zip(seq[:-1], seq[1:])):
        raise RuntimeError("inversion.depth_sequence_km must start at 0 and increase by >= 0.001 km: %s" % seq)
    extended = []
    step = (seq[-1] - seq[-2]) if len(seq) > 1 else max(half / 5, 0.001)
    while seq[-1] <= half + 1e-9:                  # extend: each step 1.3 x the previous one, at least 1 m
        step *= 1.3
        seq.append(round(seq[-1] + max(step, 0.001), 3))
        extended.append(seq[-1])
        if len(seq) > MAX_LAYERS:
            raise RuntimeError("depth_sequence_km cannot reach half the longest wavelength (%.3f km)" % half)
    upd = [v for v in seq if v <= half + 1e-9]
    below = [v for v in seq if v > half + 1e-9]
    return upd + [below[0]], dict(longest_period=float(avg[i, 0]), c=float(avg[i, 1]), half_wavelength_km=float(half),
                                  shortest_wavelength_km=float((avg[:, 0] * avg[:, 1]).min()),
                                  sequence=[float(v) for v in seq], sequence_rule=rule, sequence_extended_by=extended)


MAX_LAYERS = 200          # NL of DSurfTomo's refineGrid2LayerMdl: every node interval becomes minthk + 1 sublayers


def check_layers(nz, sublayers):
    """refineGrid2LayerMdl makes (sublayers + 1) layers per node interval and adds the half-space: <= 200 in all."""
    n = (nz - 1) * (int(sublayers) + 1) + 1
    if n > MAX_LAYERS:
        raise RuntimeError("%d depth nodes with %d sublayers make %d layers with the half-space; DSurfTomo allows %d: use "
                           "at most %d depth nodes (inversion.depth_sequence_km, fine_step_km) or fewer inversion.sublayers"
                           % (nz, sublayers, n, MAX_LAYERS, (MAX_LAYERS - 1) // (sublayers + 1) + 1))


def forward(depths, vs, periods, igr=0, minthk=3, iwave=0):
    """Phase (igr 0) or group (1) velocity of the 1-D model at the given periods, by DSurfTomo's own chain.
    surfdisp96 follows the root from one period to the next and needs increasing periods (with decreasing periods
    it can jump to another root in strongly layered models), so the periods are sorted here and the result is
    returned in the order given."""
    periods = np.asarray(periods, float)
    order = np.argsort(periods)
    out = np.empty(len(periods))
    out[order] = _forward_sorted(depths, vs, periods[order], igr, minthk, iwave)
    return out


MAX_PERIODS = 80          # surfdisp96


def _forward_sorted(depths, vs, periods, igr=0, minthk=3, iwave=0):
    check_layers(len(depths), minthk)
    if len(periods) > MAX_PERIODS:
        raise RuntimeError("%d periods: surfdisp96 (DSurfTomo, fwd1d) accepts at most %d; use fewer frequencies per "
                           "octave (dispersion.per_octave) or a narrower band" % (len(periods), MAX_PERIODS))
    exe = external_path("bin", "fwd1d")
    if not os.path.exists(exe):
        raise RuntimeError("%s not found: run `bash tools/install_external.sh dsurftomo`" % exe)
    with tempfile.TemporaryDirectory() as d:
        mf, pf = os.path.join(d, "m.txt"), os.path.join(d, "p.txt")
        with open(mf, "w") as fh:
            fh.write("%d\n" % len(depths))
            for z, v in zip(depths, vs):
                fh.write("%.6f %.6f\n" % (z, v))
        with open(pf, "w") as fh:
            fh.write("%d\n" % len(periods))
            for t in periods:
                fh.write("%.6f\n" % t)
        r = subprocess.run([exe, mf, pf, str(igr), str(minthk)] + ([str(iwave)] if iwave else []),
                           capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError("fwd1d failed: %s" % r.stderr[-500:])
        a = np.array([[float(x) for x in ln.split()] for ln in r.stdout.split("\n") if ln.strip()])
    if a.ndim != 2 or a.shape[0] != len(periods) or not np.all(np.isfinite(a[:, 1])) or np.any(a[:, 1] <= 0):
        raise RuntimeError("fwd1d returned %s values for %d periods (stderr: %s)"
                           % (a.shape, len(periods), r.stderr[-300:]))
    return a[:, 1]


def curvature_operator(z):
    """Rows R with |R m|^2 ~ integral over zeta = z / z_max of (d^2 m / d zeta^2)^2: the second derivative on the
    actual (unequal) node spacing, in units of the depth range, each row weighted by the depth interval it stands
    for - so the smoothing weight does not depend on the number or spacing of the nodes."""
    z = np.asarray(z, float)
    Z = float(z[-1] - z[0])
    R = np.zeros((len(z) - 2, len(z)))
    for i in range(1, len(z) - 1):
        h1, h2 = z[i] - z[i - 1], z[i + 1] - z[i]
        row = np.array([2 / (h1 * (h1 + h2)), -2 / (h1 * h2), 2 / (h2 * (h1 + h2))]) * Z ** 2
        R[i - 1, i - 1:i + 2] = row * np.sqrt(0.5 * (h1 + h2) / Z)
    return R


def invert1d(T, d, z, smooth=0.05, damp=0.05, iters=14, minthk=3, vs0=None, log=print):
    """Gauss-Newton on m = ln Vs.  The smoothing acts on the model itself ("jumping", Shaw & Orcutt 1985;
    Constable et al. 1987): each step solves [G; s R; damp I] m_new = [r + G m; 0; damp m], i.e. the misfit, the
    curvature of the new model (curvature_operator) and the length of the step (damping) are minimised together."""
    z = np.asarray(z, float)
    nz, nt = len(z), len(T)
    vs = np.full(nz, vs0 if vs0 else 1.087 * d[np.argmin(T)])
    L = curvature_operator(z)
    hist = []
    for it in range(iters):
        g = forward(z, vs, T, 0, minthk)
        if (~np.isfinite(g) | (g <= 0)).any():
            raise RuntimeError("1-D forward failed at periods %s" % T[~np.isfinite(g) | (g <= 0)])
        r = d - g
        G = np.zeros((nt, nz))
        for j in range(nz):
            v2 = vs.copy()
            dv = 0.01 * vs[j]
            v2[j] += dv
            G[:, j] = (forward(z, v2, T, 0, minthk) - g) / (dv / vs[j])
        A = np.vstack([G, smooth * L, damp * np.eye(nz)])
        b = np.concatenate([r, -smooth * (L @ np.log(vs)), np.zeros(nz)])          # = s R m_new = 0 for dm
        dm = np.clip(np.linalg.lstsq(A, b, rcond=None)[0], -0.15, 0.15)
        vs = vs * np.exp(dm)
        hist.append(float(np.sqrt(np.mean(r ** 2))))
        log("    1-D iteration %2d  rms %.4f km/s  Vs %.3f-%.3f" % (it + 1, hist[-1], vs.min(), vs.max()))
    g = forward(z, vs, T, 0, minthk)
    return vs, g, hist


def add_fine_nodes(z, vs, step, zmax):
    if not step:
        return np.asarray(z, float), np.asarray(vs, float)
    zz = np.unique(np.round(np.r_[z, np.arange(0, min(zmax, max(z)) + 1e-9, step)], 3))
    return zz, np.interp(zz, z, vs)


def _check_period_count(cfg, curve_dir, tmax_used):
    """Stop before the 1-D inversion if the 3-D project would carry more periods than the solver accepts."""
    from .. import export as EXPORT
    from ..tools import max_periods
    from .project import records_from_curves
    solver = cfg["inversion"]["solver"]
    curves = EXPORT.read_curves(curve_dir, "C", cfg["export"]["tmin_s"], tmax_used)
    records_from_curves(curves, max_periods(solver), solver)


def build(cfg, curve_dir, out_dir, log=print):
    """-> (z, vs, avg, info): writes average_curve.txt, vs1d.txt (and the fit) into out_dir."""
    inv = cfg["inversion"]
    e = cfg["export"]
    avg = average_curve(curve_dir, e["tmin_s"], e["tmax_s"], inv["median_min_pairs"])
    np.savetxt(os.path.join(out_dir, "average_curve.txt"), avg, fmt="%8.4f %10.4f %6d",
               header="period_s  c_km/s  pairs   (median over the pairs, >= %d pairs)" % inv["median_min_pairs"])
    _check_period_count(cfg, curve_dir, float(avg[:, 0].max()))
    z, info = depth_nodes(avg, inv["depth_sequence_km"] if inv["depth_sequence_km"] == "auto" else list(inv["depth_sequence_km"]))
    sm = inv["start_model"]
    if sm["source"] == "file":
        m = np.loadtxt(sm["file"], ndmin=2)
        if np.max(z) > m[-1, 0] + 1e-6:
            raise RuntimeError("start_model.file %s ends at %.2f km, the nodes need %.2f km" % (sm["file"], m[-1, 0], max(z)))
        vs = np.interp(z, m[:, 0], m[:, 1])
        info.update(start="file %s" % sm["file"])
        fit = None
    else:
        vs, g, hist = invert1d(avg[:, 0], avg[:, 1], z, sm["smooth"], sm["damp"], sm["iterations"], sm["minthk"], log=log)
        fit = np.c_[avg[:, 0], avg[:, 1], g]
        np.savetxt(os.path.join(out_dir, "vs1d_fit.txt"), fit, fmt="%8.4f %8.4f %8.4f", header="period_s  c_average  c_predicted")
        info.update(start="inverted", rms_km_s=hist[-1] if hist else None,
                    final_rms_km_s=float(np.sqrt(np.mean((avg[:, 1] - g) ** 2))))
    z2, vs2 = add_fine_nodes(z, vs, inv["fine_step_km"], inv["fine_zmax_km"])
    np.savetxt(os.path.join(out_dir, "vs1d.txt"), np.c_[z2, vs2], fmt="%10.4f %10.4f", header="depth_km  Vs_km/s (last node fixed)")
    info.update(depth_nodes=[float(v) for v in z2], fixed_bottom_km=float(z2[-1]), tmax_used=float(avg[:, 0].max()),
                tmin_used=float(avg[:, 0].min()))
    return z2, vs2, avg, info
