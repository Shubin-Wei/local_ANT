"""How lant runs and evaluates the Python port of DSurfTomo (dsurf_main and its modules).

A solver project is a directory with DSurfTomo.in, MOD, surfdata.dat (and MOD.true for a checkerboard).
`python -m lant.inversion.solver --threads N` runs dsurf_main in that directory (the current one) and then writes
solver.json: the numerical identity of the implementation, the threads, the seconds of every module and iteration,
and the peak resident memory.  lant counts a project as finished only with that record (lant.inversion.run.status).
It also writes DWS.dat, the DWS of every node in the last iteration (> 0 for nonzero post-rejection Jacobian column support).

evaluate() is DSurfTomo's first iteration without the inversion, in this process: the predicted travel time of every
datum for a given model, the observed time and the weight of the rejection (inversion.outlier_mad) (the columns of
residualFirst.dat).
"""

import argparse
import hashlib
import json
import os
import sys
import time
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np

F = np.float32
DWS_FILE = "DWS.dat"
RECORD = "solver.json"
# Rough peak memory of one solver process besides its arrays (interpreter, NumPy/SciPy/Numba and compiled code), bytes.
BASE_BYTES = 400e6


def identity():
    """Packages and source bytes that decide the numbers of a solver run."""
    from . import (
        dsurf_bspline,
        dsurf_disp96,
        dsurf_fmm,
        dsurf_format,
        dsurf_forward,
        dsurf_input,
        dsurf_kernels,
        dsurf_main,
        dsurf_rays,
        dsurf_rows,
        dsurf_source,
        dsurf_system,
        least_squares,
        model_regularization,
    )

    packages = {}
    for p in ("numba", "numpy", "scipy"):
        try:
            packages[p] = version(p)
        except PackageNotFoundError:
            packages[p] = "missing"
    modules = (dsurf_bspline, dsurf_disp96, dsurf_fmm, dsurf_format, dsurf_forward, dsurf_input, dsurf_kernels,
               dsurf_main, dsurf_rays, dsurf_rows, dsurf_source, dsurf_system, least_squares,
               model_regularization)  # fmt: skip
    sources = {m.__name__: hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest() for m in modules}
    # A fixed name: run as `python -m lant.inversion.solver` this module is __main__.
    sources["lant.inversion.solver"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(implementation="lant 3-D inversion: DSurfTomo a651faf forward in Python, lant regularization (lant.inversion.dsurf_main)", packages=packages,
                sources=sources)  # fmt: skip


def memory_estimate(project):
    """Estimated peak memory (bytes) of one solver process: the receiver-time derivative fields of one block of
    sources (rows x nx x ny REAL, at most dsurf_forward.FIELD_BLOCK_BYTES plus one source), the sensitivity matrix in
    CSR and the arrays of the least-squares solve (about 72 bytes per nonzero: index pairs in int64, values in
    REAL and in double for the scaled matrix, its column order; nonzeros of a straight path: 2 x cells crossed + 4 corners
    per updated depth, plus 7 per parameter for the smoothing) and BASE_BYTES."""
    from .dsurf_forward import FIELD_BLOCK_BYTES
    from .run import read_data, read_in

    _, nx, ny, nz, periods = read_in(os.path.join(project, "DSurfTomo.in"))
    data = read_data(os.path.join(project, "surfdata.dat"), periods)
    lines = open(os.path.join(project, "DSurfTomo.in")).read().split("\n")
    dvx, dvz = (float(v) for v in lines[6].split()[:2])
    crossed = 2 * (np.abs(data[:, 0] - data[:, 2]) / dvx + np.abs(data[:, 1] - data[:, 3]) / dvz) + 4
    nonzeros = float(crossed.sum()) * (nz - 1) + 7.0 * (nx - 2) * (ny - 2) * (nz - 1)
    row_bytes = 4.0 * nx * ny
    receivers = len(np.unique(data[:, 2:4], axis=0))  # a source has at most this many rows: one block overshoots
    fields = row_bytes * min(len(data), FIELD_BLOCK_BYTES / row_bytes + receivers)  # the limit by less than this
    return BASE_BYTES + fields + 72.0 * nonzeros


def evaluate(project, depths, vs, synthetic=False, threads=1):
    """Rows of residualFirst.dat for the model (depths, vs[nz, ny, nx]) on the project's data: distance, predicted
    time, observed time, weighted predicted, weighted observed, weight of the rejection rule.  synthetic: the observed
    times come from MOD.true as in the project's own run (the same noise draw)."""
    from . import dsurf_input, dsurf_kernels, dsurf_main, dsurf_system
    from .dsurf_forward import forward

    d = Path(project)
    c = dsurf_input.read_control(d / "DSurfTomo.in")
    m = dsurf_input.read_measurements(d / c.datafile, c)
    model = dsurf_input.Project(c, m, np.asarray(depths, F), np.ascontiguousarray(vs, F))
    obst = m.observed
    if synthetic:
        true_vs = dsurf_input.read_true_model(d / "MOD.true", c.nx, c.ny, c.nz)
        noise = d / dsurf_main.NOISE_FILE
        obst = dsurf_main.synthetic_data(model, true_vs, threads, None, noise_file=noise if noise.exists() else None)
    phase = dsurf_kernels.phase_velocities(model.depths, model.vs, c.periods_rc, c.minthk, threads)
    kernels = dsurf_kernels.Kernels(phase, None, None, None)
    result = forward(model, workers=threads, kernels=kernels, derivatives=False)
    dsyn = np.zeros(len(obst), F)
    dsyn[: len(result.times)] = result.times
    weight, _ = dsurf_system.reject_mad(obst - dsyn, c.threshold, m.row_periods())
    return np.c_[m.distance, dsyn, obst, dsyn * weight, obst * weight, weight].astype(float)


def _peak_rss_gb():
    try:
        for line in open("/proc/self/status"):
            if line.startswith("VmHWM:"):
                return round(int(line.split()[1]) / 1e6, 3)
    except (OSError, ValueError):
        pass
    return None


def run(directory, threads, estimate=None):
    """dsurf_main in directory with `threads` threads, then solver.json (atomically)."""
    from . import dsurf_main

    directory = Path(directory)
    started = time.time()
    timings = dsurf_main.run(directory, "DSurfTomo.in", threads, dws_file=DWS_FILE)
    record = dict(
        identity=identity(),
        threads=threads,
        wall_seconds=round(time.time() - started, 2),
        peak_memory_gb=_peak_rss_gb(),
        memory_estimate_gb=round(estimate / 1e9, 3) if estimate else None,
        iterations=timings,
    )
    tmp = directory / (RECORD + ".tmp")
    tmp.write_text(json.dumps(record, indent=1) + "\n")
    os.replace(tmp, directory / RECORD)
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description="run the Python DSurfTomo in the current project directory")
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--estimate", type=float, default=None, help="estimated peak memory (bytes), recorded")
    a = parser.parse_args(argv)
    if a.threads < 1:
        parser.error("--threads must be at least 1")
    run(".", a.threads, a.estimate)
    sys.stdout.flush()


if __name__ == "__main__":
    main()
