"""DSurfTomo's main program in Python: the iterative inversion around the compatible forward.

Ports main.f90 of DSurfTomo a651faf (MIT; Fang et al. 2015) for Rayleigh phase
velocities, statement by statement and in the same precision (REAL where the
Fortran declares REAL, real*8 periods). Steps, with their main.f90 lines:

1. Banner on the screen (111-117); input file from the first argument, else
   asked for on stdin with DSurfTomo.in as default (120-132).
2. DSurfTomo.in (134-210): data file, grid, weight0/damp, minthk, velocity
   bounds, maxiter, spfra, the four period lists, ifsyn, noiselevel,
   threshold0 (``dsurf_input.read_control``); the log file ``<input>.log`` and
   its header (155-207). Only Rayleigh phase (kmaxRg = kmaxLc = kmaxLg = 0).
3. Measurements (217-285): grouping by period block, colatitude/longitude in
   REAL radians, REAL delsph distance and obst = dist / velocity, file order
   (``dsurf_input.read_measurements``).
4. Allocation (286-309): maxnar = spfra x dall x nx x ny x nz in REAL,
   truncated; it only sizes the arrays and stops the run when the nonzeros
   exceed it (467); it does not change any number. maxvp = (nx-2)(ny-2)(nz-1).
5. MOD (314-323): nz depths, then nz x ny records of nx Vs.
6. Synthetic test, ifsyn = 1 (326-343 and CalSurfG.f90 ``synthetic``):
   MOD.true without depth record; phase velocities of the true model
   (caldespersion) written to velmap2dRc.dat with the Fortran's index
   ``pvRc((j+1)*nx+i+1)``, which labels each value with the node one row of
   longitude before it; then the receiver times of every CalSurfG row, on a
   propagation grid diced 5 x 5 instead of 8 x 8 (gdx = gdz = 5), replace
   obst in CalSurfG order: obst = t + t * gaussian() * noiselevel. gaussian()
   draws from gfortran's random_number; that stream is NOT reproduced here
   (a Marsaglia polar draw from numpy's generator with seed 0 is used), so only
   noiselevel = 0 gives native's data.
7. Iterations 1..maxiter (348-546):
   a. CalSurfG (355-359): depth kernels of the current model, receiver times
      dsyn (CalSurfG order) and the rows of G (``dsurf_forward.forward``);
   b. cbst = obst - dsyn index by index (361-363);
   c. getpercentile q25, q75 and rejection outside threshold0 x [q25, q75]:
      weight 0, residual 0, row multiplied by the weight (365-380);
   d. DWS maximum and average to the log (382-393);
   e. residualFirst.dat (iteration 1) and residualLast.dat (last iteration):
      dist, dsyn, obst, dsyn*w, obst*w, w, list-directed (396-411);
   f. smoothing rows after the data rows (414-462), m = dall + maxvp;
   g. LSMR(m, n, ..., cbst, damp, atol = btol = 1e-6, conlim = 100,
      itnlim = 400, localSize = 10) in REAL (469-489; ``dsurf_lsmr``);
   h. mean, std and rms of the weighted residuals of this iteration's
      starting model, screen and log (491-516), min/max of dv (513-516);
   i. dv clipped to +/-0.5 km/s, added to the interior Vs above the last
      depth, Vs clipped to [Minvel, Maxvel] (518-532);
   j. ``<input>Measure.dat.iterNNN`` (534-543).
8. Final model (550-592): ``<input>Measure.dat``, or for ifsyn = 1
   Vs_model.real (true) and ``<input>Syn.dat`` (inverted); closing messages.

Added for lant, off unless asked for (``run(dws_file=...)``): the DWS of every
node in the last iteration (the column sums of step 7d), in the node order of
the model file, so that lant can tell which nodes have nonzero post-rejection data sensitivity.

Not reproduced: the boundary-ray warning that CalSurfG prints after each source
once a ray has touched the model edge (screen only), group and Love data, and
gfortran's random stream (step 6).

Depth kernels come from ``dsurf_kernels`` (DSurfTomo's own surfdisp96 and
depthkernel, ported in ``dsurf_disp96``) unless ``kernels`` supplies a function
of the model (the same Kernels arrays), e.g. for a comparison.
"""

import argparse
import math
import sys
import time
from pathlib import Path

import numba
import numpy as np

from . import dsurf_format as fmt
from . import dsurf_input, dsurf_kernels, dsurf_system
from .dsurf_forward import SYNTHETIC_DICING, forward
from .dsurf_lsmr import Matrix, lsmr

F = np.float32
# main.f90 lines 474-485: LSMR controls.
ATOL = BTOL = 1e-6
CONLIM = 100
ITNLIM = 400
LOCAL_SIZE = 10
# main.f90 lines 521-526: largest velocity change per iteration, km/s.
MAX_STEP = F(0.5)
# Seed of the stand-in noise generator (step 6); not gfortran's stream.
NOISE_SEED = 0


class Output:
    """Screen (stdout, flushed per line, as GFORTRAN_UNBUFFERED_ALL would) and the log file."""

    def __init__(self, log_path, screen=sys.stdout):
        self.log = open(log_path, "w")
        self.screen = screen

    def say(self, *lines):
        for line in lines:
            self.screen.write(line + "\n")
        self.screen.flush()

    def write(self, *lines):
        for line in lines:
            self.log.write(line + "\n")
        self.log.flush()

    def close(self):
        self.log.close()


def _periods(values):
    """'(50f7.2)' with format reversion every 50 values."""
    v = list(values)
    return [fmt.fixed_record(v[i : i + 50], 7, 2) for i in range(0, max(len(v), 1), 50)]


def _header(out, c):
    grid = [
        " model origin:latitude,longitue",
        fmt.fixed_record((c.goxd, c.gozd), 10, 5),
        " grid spacing:latitude,longitue",
        fmt.fixed_record((c.dvxd, c.dvzd), 10, 5),
        " model dimension:nx,ny,nz",
        "".join(fmt.integer(v, 5) for v in (c.nx, c.ny, c.nz)),
    ]
    out.say(*grid)
    out.write(
        "",
        "                          S U R F  T O M O",
        " Python port in local_ANT of DSurfTomo (Fang et al. 2015); report problems of this port to local_ANT",
        "",
        *grid,
    )
    title = " Rayleigh wave phase velocity used,periods:(s)"
    out.say(title, *_periods(c.periods_rc))
    out.write(title, *_periods(c.periods_rc))


def model_lines(c, depths, vs):
    """'(5f10.5)' rows lon, lat, depth, Vs of the interior nodes above the last depth (main.f90 536-542)."""
    nz, ny, nx = vs.shape
    lines = []
    for k in range(nz - 1):
        for j in range(ny - 2):
            lon = F(c.gozd) + F(j) * F(c.dvzd)
            for i in range(nx - 2):
                lat = F(c.goxd) - F(i) * F(c.dvxd)
                lines.append(fmt.fixed_record((lon, lat, depths[k], vs[k, j + 1, i + 1]), 10, 5))
    return lines


def _dws_file(path, c, depths, norm):
    """lon, lat, depth and DWS of every updated node, in the order of model_lines."""
    nvz, nvx = c.ny - 2, c.nx - 2
    k, j, i = np.meshgrid(np.arange(c.nz - 1), np.arange(nvz), np.arange(nvx), indexing="ij")
    lon = F(c.gozd) + j.ravel().astype(F) * F(c.dvzd)
    lat = F(c.goxd) - i.ravel().astype(F) * F(c.dvxd)
    rows = np.c_[lon, lat, np.asarray(depths, F)[k.ravel()], norm]
    np.savetxt(path, rows, fmt=("%10.5f", "%10.5f", "%10.5f", "%14.6e"), header="lon lat depth_km DWS")


def _write_lines(path, lines):
    Path(path).write_text("".join(line + "\n" for line in lines))


def _residual_file(path, dist, dsyn, obst, weight):
    rows = zip(dist, dsyn, obst, dsyn * weight, obst * weight, weight)
    _write_lines(path, (fmt.list_record(*r) for r in rows))


def _maxnar(c, dall):
    value = F(c.spfra) * F(dall) * F(c.nx) * F(c.ny) * F(c.nz)
    if not value < 2.0**31:
        raise SystemExit("number overflow, decrease your sparsefrac")
    return int(value)


def _gaussian(rng):
    """gaussian.f90's Marsaglia polar method in REAL (stand-in for gfortran's random_number)."""
    w = F(2.0)
    while w >= F(1.0):
        x1 = F(2.0) * F(rng.random()) - F(1.0)
        x2 = F(2.0) * F(rng.random()) - F(1.0)
        w = x1 * x1 + x2 * x2
    w = F(((F(-2.0) * F(math.log(w))) / w) ** F(0.5))
    return x1 * w


def _velmap(path, c, phase):
    nz_periods, ny, nx = phase.shape
    lines = []
    for k in range(nz_periods):
        for j in range(1, ny - 1):
            lon = F(c.gozd) + F(j - 1) * F(c.dvzd)
            for i in range(1, nx - 1):
                lat = F(c.goxd) - F(i - 1) * F(c.dvxd)
                # pvRc((j+1)*nx+i+1, k) of the Fortran: column (j+2, i+1), 1-based (longitude, latitude).
                lines.append(fmt.fixed_record((lon, lat, c.periods_rc[k], phase[k, j + 1, i]), 8, 4))
    _write_lines(path, lines)


def synthetic_data(project, true_vs, workers, directory, kernels=None):
    """CalSurfG.f90 ``synthetic``: obst of every CalSurfG row from the true model (noise as step 6).

    kernels: the function of ``run``; its phase velocities replace caldespersion's (caldespersion and
    depthkernel compute them identically).  directory None: velmap2dRc.dat is not written.
    """
    c = project.control
    if kernels is None:
        phase = dsurf_kernels.phase_velocities(project.depths, true_vs, c.periods_rc, c.minthk, workers)
    else:
        phase = kernels(project.depths, true_vs).phase
    if directory is not None:
        _velmap(Path(directory) / "velmap2dRc.dat", c, phase)
    truth = dsurf_input.Project(c, project.measurements, project.depths, true_vs)
    kernels = dsurf_kernels.Kernels(phase, None, None, None)
    times = forward(truth, workers=workers, kernels=kernels, derivatives=False, dicing=SYNTHETIC_DICING).times
    obst = project.measurements.observed.copy()
    noise = F(c.noiselevel)
    rng = np.random.default_rng(NOISE_SEED)
    for i, t in enumerate(times):
        obst[i] = t + t * _gaussian(rng) * noise
    return obst


def _set_threads(workers):
    numba.set_num_threads(max(1, min(workers, numba.config.NUMBA_NUM_THREADS)))


def run(directory=".", inputfile="DSurfTomo.in", workers=1, kernels=None, screen=sys.stdout, dws_file=None):
    """Run DSurfTomo's inversion in ``directory``; returns per-iteration seconds of each module.

    kernels: optional function (depths, vs) -> dsurf_kernels.Kernels used instead of dsurf_kernels.
    dws_file: also write the DWS of every node in the last iteration to this file of ``directory``.
    """
    directory = Path(directory)
    _set_threads(workers)
    # The native banner names its author as the contact for bugs; this port is not his program, so it says whose
    # it is instead (the rest of the log keeps DSurfTomo's layout).
    screen.write("\n                  DSurfTomo (v1.4), Python port in local_ANT\n")
    screen.write(" Problems of this port: report them to local_ANT, not to the author of DSurfTomo\n\n")
    if not (directory / inputfile).exists():
        raise SystemExit("unable to open the inputfile")
    c = dsurf_input.read_control(directory / inputfile)
    out = Output(directory / (inputfile + ".log"), screen)
    try:
        _header(out, c)
        if c.kmax_rg or c.kmax_lc or c.kmax_lg:
            raise NotImplementedError("only Rayleigh phase velocities (kmaxRg = kmaxLc = kmaxLg = 0) are supported")
        m = dsurf_input.read_measurements(directory / c.datafile, c)
        dall = len(m.observed)
        maxnar = _maxnar(c, dall)
        maxvp = (c.nx - 2) * (c.ny - 2) * (c.nz - 1)
        out.say(" Number of all measurements" + fmt.integer(dall, 7))
        depths, vs = dsurf_input.read_model(directory / "MOD", c.nx, c.ny, c.nz)
        out.say(" grid points in depth direction:(km)", *_periods(depths))
        project = dsurf_input.Project(c, m, depths, vs)
        obst = m.observed
        true_vs = None
        if c.ifsyn == 1:
            out.say(" Synthetic Test Begin")
            true_vs = dsurf_input.read_true_model(directory / "MOD.true", c.nx, c.ny, c.nz)
            obst = synthetic_data(project, true_vs, workers, directory, kernels)
        timings = [_iteration(out, directory, inputfile, project, obst, maxnar, maxvp, it, workers, kernels,
                              dws_file if it == c.maxiter else None)
                   for it in range(1, c.maxiter + 1)]  # fmt: skip
        _finish(out, directory, inputfile, project, true_vs)
    finally:
        out.close()
    return timings


def _iteration(out, directory, inputfile, project, obst, maxnar, maxvp, it, workers, kernels, dws_file=None):
    c, m = project.control, project.measurements
    dall = len(obst)
    seconds = {}
    out.say(" computing sensitivity matrix...", " Rayleigh wave phase velocity depth kernel")
    start = time.perf_counter()
    if kernels is None:
        k = dsurf_kernels.depth_kernels(project.depths, project.vs, c.periods_rc, c.minthk, workers)
    else:
        k = kernels(project.depths, project.vs)
    seconds["senk"] = time.perf_counter() - start
    result = forward(project, workers=workers, kernels=k)
    seconds["traveltime"] = result.seconds["fields_and_rays"] + result.seconds["rows"]
    start = time.perf_counter()
    dsyn = np.zeros(dall, F)  # rows CalSurfG does not visit (a repeated period block) stay 0
    dsyn[: len(result.times)] = result.times
    cbst = obst - dsyn
    q25, q75 = dsurf_system.percentiles(cbst)
    weight, cbst = dsurf_system.reject(cbst, q25, q75, c.threshold)
    g = result.jacobian
    rows = np.repeat(np.arange(g.shape[0]), np.diff(g.indptr))
    values = g.data.astype(F) * weight[rows]
    cols = g.indices.astype(np.int64)
    largest, average = dsurf_system.dws(cols, values, maxvp)
    out.write(fmt.list_record("Maximum and Average DWS values:", largest, average))
    if dws_file:
        _dws_file(directory / dws_file, c, project.depths, dsurf_system.dws_nodes(cols, values, maxvp))
    if it == 1:
        _residual_file(directory / "residualFirst.dat", m.distance, dsyn, obst, weight)
    if it == c.maxiter:
        _residual_file(directory / "residualLast.dat", m.distance, dsyn, obst, weight)
    s_rows, s_cols, s_values = dsurf_system.smoothing(c.nx, c.ny, c.nz, c.weight, dall)
    if len(values) + len(s_values) > maxnar:
        raise SystemExit("increase sparsity fraction(spfra)")
    matrix = Matrix(np.r_[rows, s_rows], np.r_[cols, s_cols], np.r_[values, s_values], dall + maxvp, maxvp)
    rhs = np.r_[cbst, np.zeros(maxvp, F)]
    seconds["assembly"] = time.perf_counter() - start
    start = time.perf_counter()
    solution = lsmr(matrix, rhs, c.damp, ATOL, BTOL, CONLIM, ITNLIM, LOCAL_SIZE)
    seconds["lsmr"] = time.perf_counter() - start
    seconds["lsmr_iterations"] = solution.itn
    start = time.perf_counter()
    mean, std, rms = dsurf_system.residual_statistics(cbst)
    numbers = fmt.fixed(mean * F(1000), 8, 1) + "ms " + fmt.fixed(F(1000) * std, 8, 2) + "ms " + fmt.fixed(rms, 8, 3)
    dv = solution.x
    variation = fmt.fixed(dv.min(), 7, 4) + fmt.fixed(dv.max(), 7, 4)
    out.say(fmt.integer(it, 2) + "th iteration...", " mean,std_devs and rms of residual after weighting: " + numbers)
    out.write(fmt.integer(it, 2) + "th iteration...", "mean,std_devs and rms of residual: " + numbers)
    out.say(" min and max velocity variation " + variation)
    out.write("min and max velocity variation " + variation)
    step = np.clip(dv, -MAX_STEP, MAX_STEP).reshape(c.nz - 1, c.ny - 2, c.nx - 2)
    inner = project.vs[:-1, 1:-1, 1:-1] + step
    project.vs[:-1, 1:-1, 1:-1] = np.minimum(np.maximum(inner, F(c.minvel)), F(c.maxvel))
    _write_lines(directory / ("%sMeasure.dat.iter%03d" % (inputfile, it)), model_lines(c, project.depths, project.vs))
    seconds["update"] = time.perf_counter() - start
    return seconds


def _finish(out, directory, inputfile, project, true_vs):
    c = project.control
    out.say(" Program finishes successfully")
    out.write(" Program finishes successfully")
    if true_vs is not None:
        syn = inputfile + "Syn.dat"
        _write_lines(directory / "Vs_model.real", model_lines(c, project.depths, true_vs))
        _write_lines(directory / syn, model_lines(c, project.depths, project.vs))
        lines = (
            " Output True velocity model to Vs_model.real",
            " Output inverted shear velocity model to " + syn.ljust(100),
        )
    else:
        name = inputfile + "Measure.dat"
        _write_lines(directory / name, model_lines(c, project.depths, project.vs))
        lines = (" Output inverted shear velocity model to " + name.ljust(100),)
    out.say(*lines)
    out.write(*lines)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="DSurfTomo's inversion (main.f90) in Python, run in the current directory"
    )
    parser.add_argument("inputfile", nargs="?", help="control file (default: asked, DSurfTomo.in)")
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="threads for the depth kernels and travel times (default: lant's rule, 80 % of the idle cores)",
    )
    a = parser.parse_args(argv)
    inputfile = a.inputfile
    if inputfile is None:
        sys.stdout.write(" input file [DSurfTomo.in(default)]:\n")
        sys.stdout.flush()
        inputfile = sys.stdin.readline().strip()
        if len(inputfile) <= 1:
            inputfile = "DSurfTomo.in"
    if a.threads is None:
        from ..resources import cores

        a.threads, how = cores("auto")
        sys.stdout.write(" threads: %d (%s)\n" % (a.threads, how))
    if a.threads < 1:
        parser.error("--threads must be at least 1")
    run(".", inputfile, a.threads)


if __name__ == "__main__":
    main()
