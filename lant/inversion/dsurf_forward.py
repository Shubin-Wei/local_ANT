"""DSurfTomo-compatible forward model: predicted times and sensitivity matrix.

Follows the forward part of CalSurfG in DSurfTomo a651faf (MIT) for fundamental
Rayleigh phase velocities: depth kernels of every grid column, then for each
period and source gridder, the refined-source fast marching, srtimes and rpaths
of every receiver, and the row assembly. Rows come in CalSurfG's order (period
block, source, receiver). The 2-D kernels are in dsurf_bspline, dsurf_fmm,
dsurf_source, dsurf_rays and dsurf_rows. Phase velocities and node kernels come
from ``dsurf_kernels`` unless supplied; a ``dsurf_kernels.LazyKernels`` is asked
for the partials of the columns each block's rays touch just before its rows
are assembled.

``real=np.float32`` reproduces the Fortran REAL arithmetic of the 2-D part;
``np.float64`` runs the same algorithm in double from the same REAL inputs.
"""

import time
from dataclasses import dataclass, field

import numba
import numpy as np
from numba import njit, prange
from scipy.sparse import csr_matrix

from . import dsurf_kernels
from .dsurf_bspline import gridder
from .dsurf_input import PI32
from .dsurf_kernels import Kernels  # noqa: F401  (re-exported for callers)
from .dsurf_rays import rpaths
from .dsurf_rows import FTOL, assemble
from .dsurf_source import source_fields, srtimes

# Fixed in CalSurfG: grid dicing, source-grid dicing and extent, Earth radius (km).
GDX = GDZ = 8
# CalSurfG.f90 ``synthetic`` computes the synthetic data on a coarser propagation grid.
SYNTHETIC_DICING = 5
SGDL = 8
SGS = 8
EARTH_KM = 6371.0
# CalSurfG switches the Vp/density relation when the deepest updated node is this deep (km).
DEEP_RELATION_KM = 35.0
# Bytes of the receiver derivative fields (fdm) held at one time.  CalSurfG keeps one field per receiver and turns it
# into its row at once; here the sources are taken in blocks whose fields fit into this, every block's rows are
# assembled before the next block (sources are independent, so the numbers are those of one block).  Engineering
# choice: memory, not numbers.
FIELD_BLOCK_BYTES = 512e6


@dataclass
class ForwardResult:
    times: np.ndarray  # (rows,) predicted time, s, CalSurfG order
    jacobian: object  # scipy CSR (rows, parameters) d time / d Vs, or None
    rows: np.ndarray  # (rows, 3): 1-based period block, 0-based source, receiver
    file_rows: np.ndarray  # (rows,) index of the surfdata.dat row behind each predicted row
    seconds: dict = field(default_factory=dict)


def depth_kernels(depths, vs, periods, minthk, workers=1, interior_only=True):
    """Phase velocities and node partials of every column (dsurf_kernels.depth_kernels)."""
    return dsurf_kernels.depth_kernels(depths, vs, periods, minthk, workers, interior_only)


def grid_constants(control, real):
    """dvx, dvz, gox, goz in radians, computed as CalSurfG (REAL pi, x = colatitude)."""
    pi = real(PI32)
    dvx = real(control.dvxd) * pi / real(180.0)
    dvz = real(control.dvzd) * pi / real(180.0)
    gox = (real(90.0) - real(control.goxd)) * pi / real(180.0)
    goz = real(control.gozd) * pi / real(180.0)
    return dvx, dvz, gox, goz


@njit(cache=True, parallel=True)
def _fields(
    F,
    sx,
    sz,
    svel,
    first,
    count,
    rx,
    rz,
    velv,
    veln,
    gox,
    goz,
    dvx,
    dvz,
    earth,
    dicing,
    derivatives,
    times,
    fdm,
    status,
):
    """All sources in parallel; dicing = (gdx, gdz, sgdl, sgs). Each source owns its fields and heap."""
    gdx, gdz, sgdl, sgs = dicing
    dnx = dvx / F(gdx)
    dnz = dvz / F(gdz)
    for s in prange(sx.size):
        k = svel[s]
        r0, r1 = first[s], first[s] + count[s]
        st, ttn, ttnr, nstsr, goxr, gozr, dnxr, dnzr = source_fields(
            F, sx[s], sz[s], velv[k], veln[k], gox, goz, dvx, dvz, earth, gdx, gdz, sgdl, sgs, rx[r0:r1], rz[r0:r1]
        )
        if st != 0:
            status[s] = st
            continue
        for r in range(first[s], first[s] + count[s]):
            st, t = srtimes(F, sx[s], sz[s], rx[r], rz[r], ttn, veln[k], gox, goz, dnx, dnz, earth)
            times[r] = t
            if st == 0 and derivatives:
                st = rpaths(
                    F,
                    sx[s],
                    sz[s],
                    rx[r],
                    rz[r],
                    fdm[r],
                    ttn,
                    veln[k],
                    gox,
                    goz,
                    dnx,
                    dnz,
                    dvx,
                    dvz,
                    gdx,
                    gdz,
                    earth,
                    ttnr,
                    nstsr,
                    goxr,
                    gozr,
                    dnxr,
                    dnzr,
                )
            if st != 0:
                status[s] = st


def source_blocks(count, row_bytes, limit=None):
    """[(s0, s1)]: consecutive sources whose rows' derivative fields together stay within `limit` bytes (default
    FIELD_BLOCK_BYTES; at least one source per block)."""
    limit = FIELD_BLOCK_BYTES if limit is None else limit
    blocks, s0, held = [], 0, 0.0
    for s, n in enumerate(count):
        if s > s0 and held + n * row_bytes > limit:
            blocks.append((s0, s))
            s0, held = s, 0.0
        held += n * row_bytes
    blocks.append((s0, len(count)))
    return blocks


def _stack_rows(pieces, ncols):
    """CSR blocks of consecutive rows -> one CSR with the same data, column indices and row order."""
    data = np.concatenate([p.data for p in pieces])
    indices = np.concatenate([p.indices for p in pieces])
    indptr = [np.zeros(1, np.int64)]
    for p in pieces:
        indptr.append(p.indptr[1:].astype(np.int64) + indptr[-1][-1])
    indptr = np.concatenate(indptr)
    return csr_matrix((data, indices, indptr), shape=(len(indptr) - 1, ncols))


def _sources(measurements):
    rows = measurements.calsurfg_rows()
    starts = [i for i in range(len(rows)) if i == 0 or tuple(rows[i, :2]) != tuple(rows[i - 1, :2])]
    keys = rows[starts, :2]
    k, s = keys[:, 0] - 1, keys[:, 1]
    first = np.asarray(starts, np.int64)
    count = np.diff(np.r_[first, len(rows)])
    kk, ss, rr = rows[:, 0] - 1, rows[:, 1], rows[:, 2]
    return rows, k, s, first, count, kk, ss, rr


def forward(project, real=np.float32, workers=1, kernels=None, derivatives=True, dicing=GDX):
    """Predicted times and (optionally) the sensitivity matrix of a dsurf_input.Project.

    dicing: propagation-grid dicing gdx = gdz (8 in CalSurfG, SYNTHETIC_DICING in its ``synthetic``).
    """
    if real not in (np.float32, np.float64):
        raise ValueError("real must be np.float32 or np.float64")
    if not isinstance(workers, int) or workers < 1:
        raise ValueError("workers must be a positive integer")
    c, m = project.control, project.measurements
    if c.kmax_rg or c.kmax_lc or c.kmax_lg:
        raise NotImplementedError("only Rayleigh phase velocities (kmaxRg = kmaxLc = kmaxLg = 0) are supported")
    seconds = {}
    start = time.perf_counter()
    if kernels is None:
        kernels = depth_kernels(project.depths, project.vs, c.periods_rc, c.minthk, workers)
    seconds["depth_kernels"] = time.perf_counter() - start
    rows, k, s, first, count, kk, ss, rr = _sources(m)
    if not len(rows):
        raise ValueError("no measurements")
    velocity_index = m.period[k, s] - 1  # CalSurfG: velf = pvRc(:, periods(srcnum, knumi))
    used = np.unique(velocity_index)
    start = time.perf_counter()
    velv = np.zeros((len(used),) + kernels.phase.shape[1:], real)
    velv[:] = kernels.phase[used].astype(np.float32)  # velv = real(pv)
    veln = np.stack([gridder(real, velv[i], dicing, dicing) for i in range(len(used))])
    dvx, dvz, gox, goz = grid_constants(c, real)
    nvz, nvx = c.ny - 2, c.nx - 2
    times = np.zeros(len(rows), real)
    status = np.zeros(len(first), np.int64)
    sx, sz, svel = m.scx[k, s].astype(real), m.scz[k, s].astype(real), np.searchsorted(used, velocity_index)
    rx, rz = m.rcx[kk, ss, rr].astype(real), m.rcz[kk, ss, rr].astype(real)
    row_bytes = (nvz + 2) * (nvx + 2) * np.dtype(real).itemsize
    blocks = source_blocks(count, row_bytes) if derivatives else [(0, len(first))]
    deep = project.depths[-2] >= np.float32(DEEP_RELATION_KM)
    pieces = []
    seconds["fields_and_rays"] = seconds["rows"] = 0.0
    previous = numba.get_num_threads()
    numba.set_num_threads(min(workers, numba.config.NUMBA_NUM_THREADS))
    try:
        for s0, s1 in blocks:
            r0, r1 = first[s0], first[s1 - 1] + count[s1 - 1]
            fdm = np.zeros((r1 - r0 if derivatives else 1, nvz + 2, nvx + 2), real)
            _fields(
                real, sx[s0:s1], sz[s0:s1], svel[s0:s1], first[s0:s1] - r0, count[s0:s1], rx[r0:r1], rz[r0:r1],
                velv, veln, gox, goz, dvx, dvz, real(EARTH_KM), (dicing, dicing, SGDL, SGS), derivatives,
                times[r0:r1], fdm, status[s0:s1],
            )  # fmt: skip
            seconds["fields_and_rays"] += time.perf_counter() - start
            bad = np.flatnonzero(status)
            if len(bad):
                i = bad[0]
                raise ValueError(
                    "period block %d source %d: %s lies outside the model (code %d), as DSurfTomo would stop"
                    % (k[i] + 1, s[i] + 1, "source" if status[i] == 1 else "receiver", status[i])
                )
            if derivatives:
                if hasattr(kernels, "ensure"):  # dsurf_kernels.LazyKernels: partials of the columns these rays touch
                    kernels.ensure((np.abs(fdm) >= real(FTOL)).any(axis=0))
                start = time.perf_counter()
                pieces.append(assemble(fdm, kk[r0:r1], project.vs, kernels.vs, kernels.vp, kernels.rho, deep, real))
                seconds["rows"] += time.perf_counter() - start
            del fdm
            start = time.perf_counter()
        jacobian = None
        if derivatives:
            ncols = nvz * nvx * (len(project.depths) - 1)
            jacobian = pieces[0] if len(pieces) == 1 else _stack_rows(pieces, ncols)
    finally:
        numba.set_num_threads(previous)
    return ForwardResult(times, jacobian, rows, m.file_row[kk, ss, rr], seconds)
