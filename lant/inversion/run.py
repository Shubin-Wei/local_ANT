"""Step 4 (lant inv): 1-D starting model and 3-D direct inversion of the travel times.

  0. data: inv.data = phase (the phase travel times), group (the group travel times) or both;
  1. average curves: the median velocity of each inversion period over the pairs (periods on the
     1/inv.per_octave-octave grid through 1 Hz; a period needs max(5, min(50, a quarter of the pairs)) pairs);
  2. depth nodes (inv.depths_km, or the rule of start_model.depth_nodes) and the 1-D starting model (inv.
     start_model_file, or the 1-D inversion of the average curves, start_model.invert1d); wavelengths are taken
     from the phase velocities, or from the group velocities when only these are used (velocity x period);
  3. grid: inv.grid_km, or the median wavelength of the data / 3 (lateral resolution of surface-wave tomography is
     about a wavelength; three nodes sample it);
  4. the 3-D inversion (Fang et al. 2015, DSurfTomo's method in Python, lant.inversion.dsurf_main): every iteration
     traces the rays in the current model (fast marching), computes the depth kernels, rejects travel times further
     than inv.outlier_mad x 1.4826 MAD from the median residual of their period (phase and group apart), and
     solves for the update with
     smoothing (inv.weight) and damping (inv.damping) of the model itself (lant.inversion.model_regularization).

Results in OUTPUT/4_inv:
  dsurftomo/          the DSurfTomo project (surfdata.dat, DSurfTomo.in, MOD) and its outputs (log, residuals,
                      the model of every iteration); the same files run with the original DSurfTomo
  average_curve.csv   period, velocity (phase / group), median velocity, pairs
  start_model.csv     depth, Vs of the starting model; the fit of the average curves in start_fit.csv
  vs_model.csv        lon, lat, depth, Vs, dVs (% of the starting model), DWS (data sensitivity; 0: no data)
  residuals.csv       per period and velocity: travel times, kept, rms of the residuals of the starting and the
                      final model
"""

import os
import shutil

import numpy as np

from .. import resources
from ..step import Step, need
from . import project as P
from . import start_model as S

DIR = "4_inv"
PROJECT = "dsurftomo"


def read_dispersion(path):
    return np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")


def min_pairs(n_pairs):
    return int(max(5, min(50, n_pairs // 4)))


def run(cfg):
    from .. import config as C

    C.check_inv(cfg)
    disp = need(os.path.join(cfg["output_dir"], "3_disp", "dispersion.csv"), "disp")
    with Step(cfg, DIR, ["inv", "cores"]) as step:
        inv = dict(cfg["inv"])
        n_cores, how = resources.cores(cfg["cores"])
        step.say("cores: %d (%s)" % (n_cores, how))
        tab = read_dispersion(disp)
        if tab.size == 0:
            raise RuntimeError("%s has no accepted measurement" % disp)
        recs, periods = P.records(tab, int(inv["per_octave"]), inv["data"])
        for kind, ps in periods.items():
            if len(ps) < 2:
                raise RuntimeError("fewer than 2 measured %s-velocity periods lie on the 1/%d-octave grid"
                                   % (kind, inv["per_octave"]))  # fmt: skip
            step.say("%s velocity: %d periods %.4g-%.4g s (the measured periods on the 1/%d-octave grid through 1 Hz: "
                     "inv.per_octave)" % (kind, len(ps), ps[0], ps[-1], inv["per_octave"]))  # fmt: skip
        step.say("%d pairs, %d travel times (inv.data = %s)" % (len(recs), sum(len(r["values"]) for r in recs),
                                                                inv["data"]))  # fmt: skip
        sel = P.on_grid(tab["period_s"], int(inv["per_octave"]))
        need_pairs = min_pairs(len(recs))
        avgs = {kind: S.average_curve(np.round(tab["period_s"][sel], 6), tab[P.KINDS[kind][0]][sel], need_pairs)
                for kind in periods}  # fmt: skip
        if sum(len(a) for a in avgs.values()) < 3:
            raise RuntimeError("only %d periods have %d pairs: too few for a 1-D model"
                               % (sum(len(a) for a in avgs.values()), need_pairs))  # fmt: skip
        with open(step.path("average_curve.csv"), "w") as fh:
            fh.write("period_s,velocity,v_kms,pairs\n")
            for kind, a in avgs.items():
                for t, v, n in a:
                    fh.write("%.6f,%s,%.5f,%d\n" % (t, kind, v, n))
        step.say(
            "average curves: median over the pairs, periods with >= %d pairs (max(5, min(50, pairs / 4)))" % need_pairs
        )
        lam_kind = "phase" if "phase" in avgs and len(avgs["phase"]) else "group"
        if inv["depths_km"] is not None:
            z = np.array(inv["depths_km"], float)
        else:
            z = np.array(step.derive("inv.depths_km", [round(float(v), 3) for v in S.depth_nodes(avgs[lam_kind])],
                                     "from 0 in steps of max(shortest wavelength / 4, depth / 4) to the first node past "
                                     "half the longest wavelength (the last node: fixed half-space); wavelengths "
                                     "from the %s velocities" % lam_kind))  # fmt: skip
        if inv["start_model_file"]:
            mz, mv = S.read_node_model(inv["start_model_file"])
            if z[-1] > mz[-1] + 1e-6:
                raise RuntimeError("inv.start_model_file ends at %.3f km; the nodes reach %.3f km" % (mz[-1], z[-1]))
            vs = np.interp(z, mz, mv)
            step.say("starting model: %s interpolated to the nodes" % inv["start_model_file"])
        else:
            step.say("starting model: 1-D inversion of the average curve%s" % ("s" if len(avgs) > 1 else ""))
            T = np.concatenate([a[:, 0] for a in avgs.values()])
            d = np.concatenate([a[:, 1] for a in avgs.values()])
            grp = np.concatenate([np.full(len(a), kind == "group") for kind, a in avgs.items()])
            vs, pred, _ = S.invert1d(T, d, z, log=step.say, group=grp)
            with open(step.path("start_fit.csv"), "w") as fh:
                fh.write("period_s,velocity,v_average_kms,v_start_model_kms\n")
                for t, g_, a, b in zip(T, grp, d, pred):
                    fh.write("%.6f,%s,%.5f,%.5f\n" % (t, "group" if g_ else "phase", a, b))
            step.number("1-D fit rms (km/s)", round(float(np.sqrt(np.mean((d - pred) ** 2))), 4))
        np.savetxt(step.path("start_model.csv"), np.c_[z, vs], fmt="%.4f,%.5f", header="depth_km,vs_kms", comments="")
        if inv["grid_km"] is not None:
            spacing = float(inv["grid_km"])
        else:
            lam = float(np.median([v * t for r in recs for (k, t), v in r["values"].items() if k == lam_kind]))
            spacing = step.derive("inv.grid_km", float("%.2g" % (lam / 3)), "median wavelength of the data %.4g km / 3"
                                  % lam)  # fmt: skip
        inv["vs_min"] = (
            inv["vs_min_kms"]
            if inv["vs_min_kms"] is not None
            else step.derive(
                "inv.vs_min_kms", round(0.6 * float(vs.min()), 3), "0.6 x the slowest Vs of the starting model"
            )
        )
        inv["vs_max"] = (
            inv["vs_max_kms"]
            if inv["vs_max_kms"] is not None
            else step.derive(
                "inv.vs_max_kms", round(1.4 * float(vs.max()), 3), "1.4 x the fastest Vs of the starting model"
            )
        )
        pdir = step.path(PROJECT)
        if os.path.isdir(pdir):
            shutil.rmtree(pdir)
        os.makedirs(pdir)
        info = P.write_surfdata(recs, os.path.join(pdir, "surfdata.dat"))
        g = P.grid(info["stations"], spacing)
        P.write_control(os.path.join(pdir, "DSurfTomo.in"), g, len(z), info["nsrc"], inv, periods)
        P.write_model(os.path.join(pdir, "MOD"), z, vs, g)
        step.say("grid: %d x %d nodes (latitude x longitude) every %.3g km, %d depth nodes; %d travel times"
                 % (g["nx"], g["ny"], spacing, len(z), info["travel_times"]))  # fmt: skip
        solve(pdir, n_cores, step)
        results(pdir, step.dir, z, vs, step)
        from ..plotting import inv_figures

        inv_figures.make(cfg, step, pdir)


def solve(pdir, n_cores, step):
    from . import dsurf_main

    step.say("3-D inversion (%d threads); DSurfTomo's own log: %s/DSurfTomo.in.log" % (n_cores, PROJECT))
    timings = dsurf_main.run(pdir, "DSurfTomo.in", n_cores, dws_file="dws.txt")
    step.number("iterations", len(timings))
    step.number("travel times rejected in the last iteration", timings[-1]["rejected"])


def model_table(pdir, name):
    """Inverted model (lon lat depth Vs) and the DWS of its nodes."""
    a = np.loadtxt(os.path.join(pdir, name), ndmin=2)[:, :4]
    dws = np.loadtxt(os.path.join(pdir, "dws.txt"), ndmin=2)[:, 3] if os.path.exists(os.path.join(pdir, "dws.txt")) \
        else np.full(len(a), np.nan)  # fmt: skip
    return a, dws


def residual_table(pdir):
    """Per period: travel times, kept by the rejection, rms of the kept residuals (s) of the start and final model."""
    from . import dsurf_input

    pr = dsurf_input.read_project(pdir)
    per = pr.measurements.row_periods()
    first = np.loadtxt(os.path.join(pdir, "residualFirst.dat"), ndmin=2)
    last = np.loadtxt(os.path.join(pdir, "residualLast.dat"), ndmin=2)
    rows = []
    nrc = len(pr.control.periods_rc)
    for k, t in enumerate(pr.control.periods_all, start=1):
        m = per == k
        r0 = (first[m, 2] - first[m, 1])[first[m, 5] > 0]
        r1 = (last[m, 2] - last[m, 1])[last[m, 5] > 0]
        rows.append((t, int(k > nrc), int(m.sum()), int((last[m, 5] > 0).sum()), _rms(r0), _rms(r1),
                     1.4826 * float(np.median(np.abs(r1 - np.median(r1)))) if r1.size else np.nan))  # fmt: skip
    return np.array(rows), first, last


def _rms(r):
    return float(np.sqrt(np.mean(r**2))) if r.size else np.nan


def results(pdir, out, z, vs, step):
    a, dws = model_table(pdir, "DSurfTomo.inMeasure.dat")
    base = np.interp(a[:, 2], z, vs)
    np.savetxt(os.path.join(out, "vs_model.csv"), np.c_[a, 100 * (a[:, 3] / base - 1), dws],
               fmt=("%.6f", "%.6f", "%.4f", "%.5f", "%.3f", "%.6g"), delimiter=",",
               header="lon,lat,depth_km,vs_kms,dvs_percent,dws", comments="")  # fmt: skip
    rows, first, last = residual_table(pdir)
    with open(os.path.join(out, "residuals.csv"), "w") as fh:
        fh.write("period_s,velocity,travel_times,kept_last_iteration,rms_start_s,rms_final_s,robust_sd_final_s\n")
        for r in rows:
            fh.write("%.6f,%s,%d,%d,%.6f,%.6f,%.6f\n" % (r[0], "group" if r[1] else "phase", *r[2:]))
    k0, k1 = first[:, 5] > 0, last[:, 5] > 0
    step.number("travel times", int(len(last)))
    step.number("rms residual of the starting model (s)", round(_rms(first[k0, 2] - first[k0, 1]), 5))
    step.number("rms residual of the final model (s)", round(_rms(last[k1, 2] - last[k1, 1]), 5))
    step.number("Vs range of the model (km/s)", [round(float(a[:, 3].min()), 3), round(float(a[:, 3].max()), 3)])
