"""Depth kernels of the grid columns for the DSurfTomo-compatible forward.

The period axis is CalSurfG's period blocks: the first nrc periods are Rayleigh phase velocities (depthkernel with
igr = 0), the others Rayleigh group velocities (igr = 1); ``phase`` then holds the velocity of each block.

DSurfTomo a651faf (MIT), CalSurfG.f90 ``depthkernel`` and ``caldespersion``: for each (longitude, latitude) column
the phase velocity of the Brocher (2005) model and, for each node, the derivative of the phase velocity with respect
to its Vs.  The phase velocities are DSurfTomo's own numbers (``dsurf_disp96``, its surfdisp96.f and
refineGrid2LayerMdl statement by statement, REAL where the Fortran is REAL; bit for bit against stored native
output, docs/python_backends.md).  The derivative is lant's since 2.0.1: along the Brocher relation (Vp and density
of the node follow its Vs), two dispersion calculations per node where DSurfTomo perturbs Vs, Vp and density one at
a time (``dsurf_disp96.column_kernel``; the same sum by the chain rule).

Less work than the Fortran loop, the same numbers where CalSurfG reads them:

- identical columns are solved once (the model is REAL, so equal columns give equal results);
- partials are computed only where CalSurfG uses them: the interior columns (the frame columns enter the
  travel-time grid through their phase velocity only) and the nodes above the last depth (the last depth is not a
  parameter).  Unused partials are left at zero; ``LazyKernels`` computes them only for the columns the rays touch;
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
            "surfdisp96 found no (or no positive) velocity at %s s (DSurfTomo: 'improper initial value in disper - "
            "no zero found'; it would continue with 0)" % bad
        )
    return cg


def _velocity(z, col, mt, t, nrc):
    """Phase velocities of the first nrc periods, group velocities of the others (caldespersion, igr 0 / 1)."""
    out = np.empty(t.size)
    if nrc:
        out[:nrc] = D.column_phase(z, col, mt, t[:nrc], 0)
    if nrc < t.size:
        out[nrc:] = D.column_phase(z, col, mt, t[nrc:], 1)
    return out


def _kernel(z, col, mt, t, nrc, nodes, s):
    """Velocities and node partials (into s [3, periods, nz]) of one column: depthkernel with igr 0, then 1."""
    out = np.empty(t.size)
    if nrc:
        out[:nrc] = D.column_kernel(z, col, mt, t[:nrc], nodes, s[0][:nrc], s[1][:nrc], s[2][:nrc], 0)
    if nrc < t.size:
        out[nrc:] = D.column_kernel(z, col, mt, t[nrc:], nodes, s[0][nrc:], s[1][nrc:], s[2][nrc:], 1)
    return out


def _map(job, n, workers):
    first = job(0)  # compile / warm the Numba code before the threads enter it
    if workers == 1:
        return [first] + [job(i) for i in range(1, n)]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return [first] + list(pool.map(job, range(1, n)))


def depth_kernels(depths, vs, periods, minthk, workers=1, interior_only=True, nrc=None):
    """Velocities and node partials of every column of a (nz, ny, nx) REAL model, as CalSurfG reads them.

    nrc: the first nrc periods are phase velocities, the others group velocities (default: all phase).
    interior_only=False computes partials for every column and node (depthkernel's full loop, for comparisons).
    """
    z, vs, t, mt = _inputs(depths, vs, periods, minthk, workers)
    nrc = t.size if nrc is None else int(nrc)
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
            return _checked(_kernel(z, col, mt, t, nrc, nodes, s), t), s
        return _checked(_velocity(z, col, mt, t, nrc), t), None

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


class LazyKernels(Kernels):
    """Phase velocities of every column at once, node partials only for the columns asked for (``ensure``).

    A row of G uses the partials of a column only where the receiver derivative field of its ray is at least ftol
    (dsurf_rows); everywhere else the partials are multiplied by nothing.  The forward asks for the columns its rays
    touch, block by block, so the columns no ray reaches (on the areal waveform synthetic used to develop 2.0.1 62 % of the nodes had no
    data) cost no dispersion calculations, and G is the same as with every column's partials.  Equal columns are
    solved once, across the blocks as well.  seconds: the time spent on partials."""

    def __init__(self, depths, vs, periods, minthk, workers=1, nrc=None):
        self._z, self._vs, self._t, self._mt = _inputs(depths, vs, periods, minthk, workers)
        self._workers = workers
        self._nrc = self._t.size if nrc is None else int(nrc)
        nz, ny, nx = self._vs.shape
        phase = phase_velocities(depths, vs, periods, minthk, workers, self._nrc)
        shape = (len(self._t), ny, nx, nz)
        super().__init__(phase, np.zeros(shape), np.zeros(shape), np.zeros(shape))
        self.done = np.zeros((ny, nx), bool)
        self._cache = {}
        self.seconds = 0.0

    def ensure(self, columns):
        """Partials of the interior columns where `columns` ([ny, nx] bool) is True and not yet computed."""
        import time

        start = time.perf_counter()
        nz, ny, nx = self._vs.shape
        want = np.asarray(columns, bool).copy()
        want[[0, -1], :] = False  # the frame columns enter only through their phase velocity
        want[:, [0, -1]] = False
        want &= ~self.done
        jj, ii = np.nonzero(want)
        keys = [self._vs[:, j, i].tobytes() for j, i in zip(jj, ii)]
        todo = sorted({k for k in keys if k not in self._cache})

        def job(n):
            col = np.frombuffer(todo[n], np.float32).copy()
            s = np.zeros((3, len(self._t), nz))
            _checked(_kernel(self._z, col, self._mt, self._t, self._nrc, nz - 1, s), self._t)
            return s

        if todo:
            for key, s in zip(todo, _map(job, len(todo), self._workers)):
                self._cache[key] = s
        for j, i, key in zip(jj, ii, keys):
            s = self._cache[key]
            self.vs[:, j, i, :], self.vp[:, j, i, :], self.rho[:, j, i, :] = s[0], s[1], s[2]
        self.done |= want
        self.seconds += time.perf_counter() - start


def phase_velocities(depths, vs, periods, minthk, workers=1, nrc=None):
    """Velocity (periods, ny, nx) only, as CalSurfG's caldespersion (synthetic data); the first nrc periods phase,
    the others group velocities (default: all phase)."""
    z, vs, t, mt = _inputs(depths, vs, periods, minthk, workers)
    nrc = t.size if nrc is None else int(nrc)
    nz, ny, nx = vs.shape
    columns = vs.reshape(nz, -1).T
    unique, first_of, inverse = np.unique(columns.view(np.int32), axis=0, return_index=True, return_inverse=True)
    inverse = np.asarray(inverse).ravel()

    def job(i):
        return _checked(_velocity(z, np.ascontiguousarray(columns[first_of[i]]), mt, t, nrc), t)

    phase = np.array(_map(job, len(unique), workers))
    return phase[inverse].T.reshape(len(t), ny, nx)
