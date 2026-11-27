"""Resolution test (lant checker): the checkerboard test (Leveque et al. 1993) on the paths and periods of the
inversion.

The true model is the starting model of lant inv times 1 + a sin(pi x / s) sin(pi y / s) (a = checker.amplitude,
s = the cell size given with --size); with --flip-depth the sign of the pattern changes below that depth, to see
how anomalies of different depths are separated.  Travel times of the true model are computed for every pair and
period of the inversion (DSurfTomo's synthetic mode), Gaussian noise is added, and the data are inverted with the
same settings as the real data.  Noise: checker.noise auto takes, per period, the robust standard deviation
(1.4826 MAD) of the residuals the 3-D model leaves (so the synthetic data are about as noisy as the real data,
Rawlinson & Spakman 2016); a number is used as the standard deviation in seconds for every period.

Results in OUTPUT/4_checker_<size>km[_flip<depth>km]:
  dsurftomo/     the project (MOD.true is the true model)
  recovery.csv   per depth: nodes with data (DWS > 0), correlation and regression slope of recovered on true
                 perturbation; the slope is the fraction of the anomaly the inversion recovers at that size
"""

import os
import shutil

import numpy as np

from .. import resources
from ..step import Step, need
from . import dsurf_input
from . import project as P
from .run import DIR as INV_DIR
from .run import PROJECT, residual_table


def folder(size, flip):
    return "4_checker_%gkm%s" % (size, "" if flip is None else "_flip%gkm" % flip)


def noise_file(cfg, src, dst):
    """noise_seconds.txt: period and standard deviation (s) of the noise added to each period's travel times."""
    rows, _, _ = residual_table(src)
    v = cfg["checker"]["noise"]
    sd = rows[:, 5] if v == "auto" else np.full(len(rows), float(v))
    np.savetxt(os.path.join(dst, "noise_seconds.txt"), np.c_[rows[:, 0], np.nan_to_num(sd)], fmt="%.6f %.6g",
               header="period_s noise_sd_s")  # fmt: skip
    return sd


def run(cfg, size, flip=None):
    src = os.path.join(cfg["output_dir"], INV_DIR, PROJECT)
    need(os.path.join(src, "DSurfTomo.inMeasure.dat"), "inv")
    with Step(cfg, folder(size, flip), ["inv", "checker", "cores"]) as step:
        n_cores, how = resources.cores(cfg["cores"])
        step.say("cores: %d (%s)" % (n_cores, how))
        pdir = step.path(PROJECT)
        if os.path.isdir(pdir):
            shutil.rmtree(pdir)
        os.makedirs(pdir)
        for f in ("surfdata.dat", "MOD"):
            shutil.copy(os.path.join(src, f), pdir)
        pr = dsurf_input.read_project(src)
        c = pr.control
        g = dict(nx=c.nx, ny=c.ny, goxd=float(c.goxd), gozd=float(c.gozd), dvxd=float(c.dvxd), dvzd=float(c.dvzd))
        z = np.asarray(pr.depths, float)
        vs1d = np.asarray(pr.vs[:, 0, 0], float)
        amp = float(cfg["checker"]["amplitude"])
        pattern = P.checker_pattern(g, float(size), amp)
        sign = np.ones(len(z)) if flip is None else np.where(z < float(flip), 1.0, -1.0)
        P.write_model(os.path.join(pdir, "MOD.true"), z, vs1d, g, header=False, pattern=pattern, depth_sign=sign)
        inv = dict(cfg["inv"], vs_min=float(c.minvel), vs_max=float(c.maxvel), weight=float(c.weight),
                   damping=float(c.damp), iterations=c.maxiter, outlier_mad=float(c.threshold))  # fmt: skip
        P.write_control(os.path.join(pdir, "DSurfTomo.in"), g, c.nz, c.nsrc, inv, c.periods_rc, synthetic=1)
        sd = noise_file(cfg, src, pdir)
        step.say("checkerboard: cells %g km, +-%g %%%s; noise per period %.3g-%.3g ms (checker.noise = %s)"
                 % (size, 100 * amp, "" if flip is None else ", sign reversed below %g km" % flip,
                    1e3 * np.nanmin(sd), 1e3 * np.nanmax(sd), cfg["checker"]["noise"]))  # fmt: skip
        from . import dsurf_main

        dsurf_main.run(pdir, "DSurfTomo.in", n_cores, dws_file="dws.txt")
        rec = recovery(pdir, z, vs1d)
        np.savetxt(step.path("recovery.csv"), rec, fmt=("%.4f", "%d", "%.3f", "%.3f"), delimiter=",",
                   header="depth_km,nodes_with_data,correlation,slope", comments="")  # fmt: skip
        for row in rec:
            step.say("  depth %.3f km: %d nodes with data, correlation %.2f, slope %.2f" % tuple(row))
        step.number("median slope over the depths", round(float(np.nanmedian(rec[:, 3])), 3))
        from ..plotting import inv_figures

        inv_figures.checker(cfg, step, pdir, z, vs1d)


def recovery(pdir, z, vs1d):
    true = np.loadtxt(os.path.join(pdir, "Vs_model.real"), ndmin=2)
    syn = np.loadtxt(os.path.join(pdir, "DSurfTomo.inSyn.dat"), ndmin=2)
    dws = np.loadtxt(os.path.join(pdir, "dws.txt"), ndmin=2)[:, 3]
    rows = []
    for zk in np.unique(np.round(syn[:, 2], 4)):
        m = (np.abs(syn[:, 2] - zk) < 1e-4) & (dws > 0)
        b = np.interp(zk, z, vs1d)
        x, y = true[m, 3] / b - 1, syn[m, 3] / b - 1
        if m.sum() < 3 or np.std(x) == 0:
            rows.append((zk, int(m.sum()), np.nan, np.nan))
            continue
        rows.append((zk, int(m.sum()), float(np.corrcoef(x, y)[0, 1]), float(np.dot(x, y) / np.dot(x, x))))
    return np.array(rows)
