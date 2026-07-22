#!/usr/bin/env python3
"""
Batch 1D dispersion inversion with evodcinv, one inversion per station pair.

All parameters are read from a YAML config file (see lant_config.yaml); nothing
is hard-coded here. For each pair the valid dispersion points (phase, group, or
both) are inverted into a 1D Vs(depth) model (CPSO), assigned to the pair
midpoint, and written out with a README-style figure. Inversions run in
parallel (spawn). Use lant_stitch.py afterwards to assemble the pseudo-2D
section.

Usage:
    cp ../configs/lant.example.yaml my_line.yaml   # then edit paths/params
    python lant_invert.py --config my_line.yaml [--max-files N]
"""

import os

# Pin native threads before importing numpy/evodcinv; parallelism is by process.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys
import glob
import argparse
import traceback
import multiprocessing as mp

import yaml
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from tqdm import tqdm
from evodcinv import EarthModel, Layer, Curve

# Per-worker config dict (populated by _init_worker in each process).
_CFG = None


def _init_worker(cfg):
    global _CFG
    _CFG = cfg
    # Route this worker's output to a per-worker log so the terminal stays clean.
    # evodcinv's progress bar writes to the underlying stderr file descriptor and
    # bypasses Python-level redirection, so redirect fds 1/2 at the OS level.
    logdir = os.path.join(cfg["outdir"], "logs")
    os.makedirs(logdir, exist_ok=True)
    logf = open(os.path.join(logdir, f"worker_{os.getpid()}.log"), "a", buffering=1)
    os.dup2(logf.fileno(), 1)
    os.dup2(logf.fileno(), 2)
    sys.stdout = logf
    sys.stderr = logf


def read_curve(path):
    """Return (period, vel, (midlon, midlat)) for valid (vel>0) points."""
    with open(path) as f:
        h1 = f.readline().split()
        h2 = f.readline().split()
    lonA, latA = float(h1[0]), float(h1[1])
    lonB, latB = float(h2[0]), float(h2[1])
    d = np.loadtxt(path, skiprows=2)
    keep = d[:, 1] > 0
    period, vel = d[keep, 0], d[keep, 1]
    return period, vel, (0.5 * (lonA + lonB), 0.5 * (latA + latB))


def vs_profile_on_grid(model_array, zmax, dz):
    """Convert evodcinv best model (thickness, vp, vs, rho) to Vs on a depth grid."""
    thickness = np.asarray(model_array)[:, 0]
    vs = np.asarray(model_array)[:, 2]
    ztop = np.concatenate(([0.0], np.cumsum(thickness)))
    zgrid = np.arange(0.0, zmax + 1e-9, dz)
    layer_idx = np.clip(np.searchsorted(ztop, zgrid, side="right") - 1, 0, len(vs) - 1)
    return zgrid, vs[layer_idx]


def _types_for(curve_type):
    """phase -> ['phase']; group -> ['group']; both -> ['phase', 'group']."""
    return ["phase", "group"] if curve_type == "both" else [curve_type]


def invert_one(task):
    index, ntot, name = task
    cfg = _CFG
    idir = cfg["input_dir"]
    try:
        # Build one Curve per requested type (CDisp=phase, GDisp=group). For
        # "both" the pair must have enough valid points in BOTH curves.
        curves, obs, midpoint, ndata = [], {}, None, 0
        for t in _types_for(cfg["curve_type"]):
            prefix = "CDisp." if t == "phase" else "GDisp."
            period, vel, midpoint = read_curve(os.path.join(idir, f"{prefix}{name}.txt"))
            if len(period) < cfg["min_points"]:
                return ("skip", name, len(period))
            curves.append(Curve(period, vel, 0, "rayleigh", t))
            obs[t] = (period, vel)
            ndata += len(period)
        midlon, midlat = midpoint

        model = EarthModel()
        for thickness, vs in cfg["layers"]:
            model.add(Layer(np.array(thickness), np.array(vs)))
        model.configure(
            optimizer="cpso", misfit="rmse",
            optimizer_args={"popsize": cfg["popsize"], "maxiter": cfg["maxiter"], "seed": 0},
        )
        res = model.invert(curves, maxrun=cfg["maxrun"])

        misfit = float(np.min(res.misfit))
        zgrid, vsgrid = vs_profile_on_grid(res.model, cfg["zmax"], cfg["dz"])

        mpath = os.path.join(cfg["outdir"], "models", f"{name}.txt")
        with open(mpath, "w") as f:
            f.write(f"# midlon {midlon:.6f} midlat {midlat:.6f} "
                    f"misfit {misfit:.5f} npts {ndata} curve {cfg['curve_type']}\n")
            f.write("# depth_km vs_kms\n")
            np.savetxt(f, np.column_stack([zgrid, vsgrid]), fmt="%.3f %.4f")

        if cfg["plot"]:
            # Figure: [Vs model | one dispersion panel per curve type | misfit].
            # Model cloud coloured by misfit, best model as a red dashed line.
            # plot_curve/plot_model only draw the *computed* model curves, so the
            # observed data is plotted manually as a black solid line (the real-data
            # equivalent of the README synthetic "True"). Real data has no true Vs
            # model, so the model panel shows only cloud + best. For "both", phase
            # and group get separate panels so the fits stay readable.
            best_args = {"color": "red", "linestyle": "--", "linewidth": 2, "label": "Best"}
            types = list(obs.keys())
            ncols = 2 + len(types)
            fig, ax = plt.subplots(1, ncols, figsize=(5 * ncols, 5), constrained_layout=True)

            # Panel 0: Vs model
            res.plot_model("vs", zmax=cfg["zmax"], show="all", ax=ax[0])
            res.plot_model("vs", zmax=cfg["zmax"], show="best", ax=ax[0], plot_args=best_args)
            ax[0].set_xlabel("Vs (km/s)")
            ax[0].set_ylabel("Depth (km)")
            ax[0].set_title(f"{name}  Vs model (misfit {misfit:.3f})")
            ax[0].legend(loc="lower left", fontsize=8)

            # One dispersion panel per curve type: observed (black) + best (red).
            for j, t in enumerate(types):
                a = ax[1 + j]
                per, v = obs[t]
                res.plot_curve(per, 0, "rayleigh", t, show="all", ax=a)
                a.plot(per, v, color="black", linewidth=2, label="Observed")
                res.plot_curve(per, 0, "rayleigh", t, show="best", ax=a, plot_args=best_args)
                # crop y-axis to observed range (model cloud can have spurious
                # group-velocity spikes that would blow up the axis)
                a.set_ylim(max(0.0, v.min() - 0.5), v.max() + 0.5)
                a.set_xlabel("Period (s)")
                a.set_ylabel(f"{t.capitalize()} velocity (km/s)")
                a.set_title(f"Rayleigh {t}")
                a.legend(loc="best", fontsize=8)

            # Last panel: misfit vs iteration
            res.plot_misfit(ax=ax[-1])
            ax[-1].set_xlabel("Iteration")
            ax[-1].set_ylabel("Misfit")
            ax[-1].set_title("Convergence")

            # shared misfit colourbar (matches the viridis_r cloud); filter out
            # the inf misfits assigned to rejected models so the range is sane.
            mf = np.asarray(res.misfits)
            mf = mf[np.isfinite(mf)]
            lo, hi = (float(mf.min()), float(mf.max())) if mf.size else (0.0, 1.0)
            if hi <= lo:
                hi = lo + 1e-6
            sm = ScalarMappable(Normalize(lo, hi), "viridis_r")
            sm.set_array([])
            fig.colorbar(sm, ax=ax.ravel().tolist(), orientation="horizontal",
                         fraction=0.04, pad=0.02, label="Misfit value")
            fig.savefig(os.path.join(cfg["outdir"], "figs", f"{name}.png"), dpi=130)
            plt.close(fig)

        return ("ok", name, midlon, midlat, misfit)
    except Exception as exc:
        print(f"ERROR {name}: {exc}", flush=True)
        traceback.print_exc(file=sys.stderr)
        if ("mpath" in locals() and os.path.isfile(mpath) and
                all(v in locals() for v in ("midlon", "midlat", "misfit"))):
            print(f"WARNING {name}: model was saved; keeping it despite "
                  f"post-inversion diagnostic failure.", flush=True)
            return ("ok", name, midlon, midlat, misfit)
        return ("err", name, str(exc))


def list_pairs(input_dir, curve_type):
    """Pair names available for the requested curve type (intersection for 'both')."""
    def names(prefix):
        return {os.path.basename(p)[len(prefix):-len(".txt")]
                for p in glob.glob(os.path.join(input_dir, prefix + "*.txt"))}
    if curve_type == "phase":
        s = names("CDisp.")
    elif curve_type == "group":
        s = names("GDisp.")
    else:
        s = names("CDisp.") & names("GDisp.")
    return sorted(s)


def build_cfg(config):
    """Flatten the YAML config into the per-worker cfg dict."""
    paths = config["paths"]
    inv = config["inversion"]
    opt = inv["optimizer"]
    return {
        "input_dir": paths["input_dir"],
        "outdir": paths["output_dir"],
        "curve_type": inv["curve_type"],
        "min_points": int(inv["min_points"]),
        "zmax": float(inv["zmax"]),
        "dz": float(inv["dz"]),
        "layers": [(layer["thickness"], layer["vs"]) for layer in inv["layers"]],
        "popsize": int(opt["popsize"]),
        "maxiter": int(opt["maxiter"]),
        "maxrun": int(opt["maxrun"]),
        "plot": bool(inv["plot"]),
        "workers": int(inv["workers"]),
    }


def main():
    p = argparse.ArgumentParser(
        description="Batch 1D evodcinv inversion. Defaults come from --config; "
                    "any flag below overrides the config value.")
    p.add_argument("--config", required=True,
                   help="path to YAML config (template: configs/lant.example.yaml).")
    p.add_argument("-n", "--max-files", type=int, default=0,
                   help="process only the first N pairs (0 = all); for testing.")
    # optional overrides of the config's 'inversion' section
    p.add_argument("--curve-type", choices=["phase", "group", "both"],
                   help="override inversion.curve_type")
    p.add_argument("--workers", type=int, help="override inversion.workers")
    p.add_argument("--no-plot", action="store_true",
                   help="disable per-pair figures (override inversion.plot)")
    p.add_argument("--min-points", type=int, help="override inversion.min_points")
    p.add_argument("--zmax", type=float, help="override inversion.zmax")
    p.add_argument("--dz", type=float, help="override inversion.dz")
    p.add_argument("--popsize", type=int, help="override optimizer.popsize")
    p.add_argument("--maxiter", type=int, help="override optimizer.maxiter")
    p.add_argument("--maxrun", type=int, help="override optimizer.maxrun")
    args = p.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    cfg = build_cfg(config)

    # apply command-line overrides (None = keep config value)
    for key, val in (("curve_type", args.curve_type), ("workers", args.workers),
                     ("min_points", args.min_points), ("zmax", args.zmax), ("dz", args.dz),
                     ("popsize", args.popsize), ("maxiter", args.maxiter),
                     ("maxrun", args.maxrun)):
        if val is not None:
            cfg[key] = val
    if args.no_plot:
        cfg["plot"] = False

    pairs = list_pairs(cfg["input_dir"], cfg["curve_type"])
    if args.max_files > 0:
        pairs = pairs[:args.max_files]
    npairs = len(pairs)
    nworkers = cfg["workers"]
    print(f"\nPairs: {npairs}   Workers: {nworkers}   Curve: {cfg['curve_type']}")

    for d in ("models", "figs", "logs"):
        os.makedirs(os.path.join(cfg["outdir"], d), exist_ok=True)
    print(f"(per-worker logs: {cfg['outdir']}/logs/worker_*.log)\n")

    tasks = [(i, npairs, pairs[i]) for i in range(npairs)]
    ok, skipped, errored = [], 0, 0
    ctx = mp.get_context("spawn")
    with ctx.Pool(nworkers, initializer=_init_worker, initargs=(cfg,)) as pool:
        pbar = tqdm(pool.imap_unordered(invert_one, tasks, chunksize=1),
                    total=npairs, desc="Inverting", unit="pair", dynamic_ncols=True)
        for r in pbar:
            if r[0] == "ok":
                ok.append(r)
            elif r[0] == "skip":
                skipped += 1
            else:
                errored += 1
            pbar.set_postfix(ok=len(ok), skip=skipped, err=errored)

    summ = os.path.join(cfg["outdir"], "summary.txt")
    with open(summ, "w") as f:
        f.write("# pair midlon midlat misfit\n")
        for _, name, midlon, midlat, misfit in sorted(ok, key=lambda r: r[1]):
            f.write(f"{name} {midlon:.6f} {midlat:.6f} {misfit:.5f}\n")
    print(f"\nDone: {len(ok)} ok, {skipped} skipped (<{cfg['min_points']} pts), {errored} errored")
    print(f"Models: {cfg['outdir']}/models/   Figures: {cfg['outdir']}/figs/")
    print(f"Summary: {summ}")


if __name__ == "__main__":
    main()
