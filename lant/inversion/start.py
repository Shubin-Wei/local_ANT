"""Average curve, depth nodes and the 1-D starting model.

  average curve: per measured period the median over the pairs (robust average), kept where at least
      inversion.median_min_pairs pairs have a value; used as measured (not smoothed, not forced monotonic)
      (Jing et al. 2025; Xie et al. 2024: the average dispersion as the reference / starting model);
  depth nodes: the nodes of inversion.depth_sequence_km not deeper than half the longest wavelength of the
      average curve are updated, the next node below is the fixed bottom node (DSurfTomo neither updates nor writes
      the last node; its velocity is a boundary assumption).  "Half the longest wavelength" is the empirical
      investigation depth of Foti et al. (2018, sect. 4.2.3).  Then a node every fine_step_km is added down to
      fine_zmax_km (linear interpolation of the 1-D model there);
  1-D model: damped and smoothed Gauss-Newton on ln Vs with a numerical Jacobian, forward = DSurfTomo's own chain
      (Brocher 2005 Vp and density, refineGrid2LayerMdl, surfdisp96) through the small program fwd1d; start
      1.087 x the phase velocity at the shortest period, second-difference smoothing (weight `smooth`), damping
      `damp`, `iterations` steps limited to +-15 %.  DSurfTomo's own GenerateIniMOD.py sets Vs equal to the phase
      velocity (about 9 % slow) and is not used.
"""
import os
import subprocess
import tempfile

import numpy as np

from .. import export
from ..tools import external_path


def average_curve(curve_dir, tmin, tmax, min_pairs):
    acc = {}
    for c in export.read_curves(curve_dir, "C", tmin, tmax):
        for t, v in zip(c["T"], c["v"]):
            acc.setdefault(round(float(t), 4), []).append(v)
    rows = [(t, float(np.median(v)), len(v)) for t, v in sorted(acc.items()) if len(v) >= min_pairs]
    if len(rows) < 3:
        raise RuntimeError("average curve: only %d periods have >= %d pairs (inversion.median_min_pairs)" % (len(rows), min_pairs))
    return np.array(rows)


def depth_nodes(avg, seq):
    """Updatable nodes (not deeper than half the longest wavelength) + the next node as the fixed bottom node.  When
    the sequence is not deep enough it is extended, each step 1.3 times the previous one (reported)."""
    i = int(np.argmax(avg[:, 0]))
    half = avg[i, 0] * avg[i, 1] / 2
    seq = [float(v) for v in seq]
    extended = []
    while seq[-1] <= half + 1e-9:
        step = 1.3 * (seq[-1] - seq[-2]) if len(seq) > 1 else max(half / 5, 0.1)
        seq.append(round(seq[-1] + step, 2))
        extended.append(seq[-1])
    upd = [v for v in seq if v <= half + 1e-9]
    below = [v for v in seq if v > half + 1e-9]
    return upd + [below[0]], dict(longest_period=float(avg[i, 0]), c=float(avg[i, 1]), half_wavelength_km=float(half),
                                  sequence_extended_by=extended)


MAX_LAYERS = 200          # NL of DSurfTomo's refineGrid2LayerMdl: every node interval becomes minthk + 1 sublayers


def check_layers(nz, sublayers):
    if (nz - 1) * (int(sublayers) + 1) > MAX_LAYERS:
        raise RuntimeError("%d depth nodes with %d sublayers make %d layers; DSurfTomo allows %d: use at most %d depth nodes "
                           "(inversion.depth_sequence_km, fine_step_km) or fewer inversion.sublayers"
                           % (nz, sublayers, (nz - 1) * (sublayers + 1), MAX_LAYERS, MAX_LAYERS // (sublayers + 1) + 1))


def forward(depths, vs, periods, igr=0, minthk=3, iwave=0):
    check_layers(len(depths), minthk)
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
    return a[:, 1]


def invert1d(T, d, z, smooth=2.0, damp=0.05, iters=14, minthk=3, vs0=None, log=print):
    z = np.asarray(z, float)
    nz, nt = len(z), len(T)
    vs = np.full(nz, vs0 if vs0 else 1.087 * d[np.argmin(T)])
    L = np.zeros((nz - 2, nz))
    for i in range(nz - 2):
        L[i, i:i + 3] = (1.0, -2.0, 1.0)
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
        b = np.concatenate([r, np.zeros(nz - 2), np.zeros(nz)])
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


def build(cfg, curve_dir, out_dir, log=print):
    """-> (z, vs, avg, info): writes average_curve.txt, vs1d.txt (and the fit) into out_dir."""
    inv = cfg["inversion"]
    e = cfg["export"]
    avg = average_curve(curve_dir, e["tmin_s"], e["tmax_s"], inv["median_min_pairs"])
    np.savetxt(os.path.join(out_dir, "average_curve.txt"), avg, fmt="%8.4f %10.4f %6d",
               header="period_s  c_km/s  pairs   (median over the pairs, >= %d pairs)" % inv["median_min_pairs"])
    z, info = depth_nodes(avg, list(inv["depth_sequence_km"]))
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
