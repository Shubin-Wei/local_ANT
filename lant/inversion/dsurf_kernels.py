"""Depth kernels of every grid column for the DSurfTomo-compatible forward.

DSurfTomo a651faf (MIT), CalSurfG.f90 ``depthkernel`` and ``caldespersion``: for each (longitude, latitude) column
the phase velocity of the Brocher (2005) model and, for each node, the central differences of the phase velocity
when the node Vs, Vp or density alone changes by +/-0.5 %.  The arithmetic is DSurfTomo's own (``dsurf_disp96``, its
surfdisp96.f and refineGrid2LayerMdl statement by statement, REAL where the Fortran is REAL): the numbers equal the
native depthkernel bit for bit (tests against stored native output; docs/python_backends.md).

Less work than the Fortran loop, the same numbers where CalSurfG reads them:

- identical columns are solved once (the model is REAL, so equal columns give equal results);
- partials are computed only where CalSurfG uses them: the interior columns (the frame columns enter the
  travel-time grid through their phase velocity only) and the nodes above the last depth (the last depth is not a
  parameter).  Unused partials are left at zero;
- columns run in a thread pool; the column code is compiled with ``nogil``.
"""

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np

from . import dsurf_disp96 as D

F4 = np.float32


@dataclass
class Kernels:
    """Phase velocity (periods, ny, nx) and node partials (periods, ny, nx, nz) in double."""

    phase: np.ndarray
    vs: np.ndarray
    vp: np.ndarray
    rho: np.ndarray


def _inputs(depths, vs, periods, minthk, workers):
    vs = np.asarray(vs)
    if vs.ndim != 3:
        raise ValueError("model must be (nz, ny, nx)")
    if not isinstance(workers, int) or workers < 1:
        raise ValueError("workers must be a positive integer")
    vs = np.asarray(vs, F4)
    if not np.all(np.isfinite(vs)) or np.any(vs <= 0):
        raise ValueError("model must be finite positive Vs")
    z = np.ascontiguousarray(depths, F4)
    if z.ndim != 1 or z.size != vs.shape[0] or z[0] != 0 or np.any(np.diff(z) <= 0):
        raise ValueError("depths must start at zero and strictly increase (km), one per model depth")
    t = np.ascontiguousarray(periods, np.float64)
    return z, vs, t, F4(minthk)


def _checked(cg, periods):
    if not np.all(cg > 0) or not np.all(np.isfinite(cg)):
        bad = np.asarray(periods)[~(cg > 0)]
        raise RuntimeError(
            "surfdisp96 found no phase velocity at %s s (DSurfTomo: 'improper initial value in disper - no zero "
            "found'; it would continue with 0)" % bad
        )
    return cg


def _map(job, n, workers):
    first = job(0)  # compile / warm the Numba code before the threads enter it
    if workers == 1:
        return [first] + [job(i) for i in range(1, n)]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return [first] + list(pool.map(job, range(1, n)))


def depth_kernels(depths, vs, periods, minthk, workers=1, interior_only=True):
    """Phase velocities and node partials of every column of a (nz, ny, nx) REAL model, as CalSurfG reads them.

    interior_only=False computes partials for every column and node (depthkernel's full loop, for comparisons).
    """
    z, vs, t, mt = _inputs(depths, vs, periods, minthk, workers)
    nz, ny, nx = vs.shape
    columns = vs.reshape(nz, -1).T
    interior = np.zeros((ny, nx), bool)
    interior[1:-1, 1:-1] = True
    wanted = np.ones(ny * nx, bool) if not interior_only else interior.ravel()
    # One job per distinct (column, needs partials); the partial job also yields the phase.
    keys = np.concatenate([columns.view(np.int32), wanted[:, None].astype(np.int32)], axis=1)
    unique, first_of, inverse = np.unique(keys, axis=0, return_index=True, return_inverse=True)
    inverse = np.asarray(inverse).ravel()
    nodes = nz if not interior_only else nz - 1

    def job(i):
        col = np.ascontiguousarray(columns[first_of[i]])
        if unique[i, -1]:
            s = np.zeros((3, len(t), nz))
            cg = D.column_kernel(z, col, mt, t, nodes, s[0], s[1], s[2])
            return _checked(cg, t), s
        return _checked(D.column_phase(z, col, mt, t), t), None

    results = _map(job, len(unique), workers)
    phase = np.empty((len(unique), len(t)))
    partials = np.zeros((len(unique), 3, len(t), nz))
    for i, (c, p) in enumerate(results):
        phase[i] = c
        if p is not None:
            partials[i] = p
    phase = phase[inverse].T.reshape(len(t), ny, nx)
    partials = partials[inverse].reshape(ny, nx, 3, len(t), nz).transpose(2, 3, 0, 1, 4)
    return Kernels(phase, *(np.ascontiguousarray(partials[i]) for i in range(3)))


def phase_velocities(depths, vs, periods, minthk, workers=1):
    """Phase velocity (periods, ny, nx) only, as CalSurfG's caldespersion (synthetic data)."""
    z, vs, t, mt = _inputs(depths, vs, periods, minthk, workers)
    nz, ny, nx = vs.shape
    columns = vs.reshape(nz, -1).T
    unique, first_of, inverse = np.unique(columns.view(np.int32), axis=0, return_index=True, return_inverse=True)
    inverse = np.asarray(inverse).ravel()

    def job(i):
        return _checked(D.column_phase(z, np.ascontiguousarray(columns[first_of[i]]), mt, t), t)

    phase = np.array(_map(job, len(unique), workers))
    return phase[inverse].T.reshape(len(t), ny, nx)
