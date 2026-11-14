"""Sensitivity-matrix rows of DSurfTomo's CalSurfG (fdm -> d time / d Vs).

Numba kernels mirroring the row assembly of DSurfTomo a651faf (MIT),
CalSurfG.f90: vertices with |fdm| >= ftol (1e-4) are combined with the node
depth kernels through the empirical Brocher (2005) chain (coefficients coe_a,
coe_rho in REAL arithmetic; the alternative relation when depz(nz-1) >= 35 km),
products in double precision rounded to REAL, and entries with |row| > ftol
kept. Columns are (k-1)*nvz*nvx + (jj-1)*nvx + kk (0-based here), i.e. the
interior nodes in depth, longitude, latitude order; nodes of the last depth
and of the outer frame are fixed. Rows are returned as CSR.
"""

import numpy as np
from numba import njit, prange
from scipy.sparse import csr_matrix

FTOL = np.float32(1e-4)


@njit(cache=True)
def _chain(F, v, shallow):
    """coe_a = dVp/dVs and coe_rho = drho/dVs at Vs v, as the two Fortran branches."""
    v2 = v * v
    v3 = v2 * v
    v4 = v2 * v2
    if shallow:
        coe_a = F(2.0947) - F(0.8206) * F(2.0) * v + F(0.2683) * F(3.0) * v2 - F(0.0251) * F(4.0) * v3
        vp = F(0.9409) + F(2.0947) * v - F(0.8206) * v2 + F(0.2683) * v3 - F(0.0251) * v4
    else:
        coe_a = F(2.2110) - F(0.8984) * F(2.0) * v + F(0.2786) * F(3.0) * v2 - F(0.02412) * F(4.0) * v3
        vp = F(0.9098) + F(2.2110) * v - F(0.8984) * v2 + F(0.2786) * v3 - F(0.02412) * v4
    p2 = vp * vp
    p3 = p2 * vp
    p4 = p2 * p2
    coe_rho = coe_a * (
        F(1.6612)
        - F(0.4721) * F(2.0) * vp
        + F(0.0671) * F(3.0) * p2
        - F(0.0043) * F(4.0) * p3
        + F(0.000106) * F(5.0) * p4
    )
    return coe_a, coe_rho


@njit(cache=True)
def _row(F, fdm, vels, sen_vs, sen_vp, sen_rho, shallow, ftol, row):
    """One dense row; vels (nz, ny, nx) REAL, sen_* (ny, nx, nz) double for this period."""
    nz = vels.shape[0]
    nvz = fdm.shape[0] - 2
    nvx = fdm.shape[1] - 2
    row[:] = 0
    for jj in range(1, nvz + 1):
        for kk in range(1, nvx + 1):
            f = fdm[jj, kk]
            if abs(f) >= ftol:
                for k in range(nz - 1):
                    coe_a, coe_rho = _chain(F, vels[k, jj, kk], shallow)
                    value = (sen_vp[jj, kk, k] * np.float64(coe_a) + sen_rho[jj, kk, k] * np.float64(coe_rho)) + sen_vs[
                        jj, kk, k
                    ]
                    row[k * nvz * nvx + (jj - 1) * nvx + (kk - 1)] = F(value * np.float64(f))


@njit(cache=True, parallel=True)
def _count(F, fdm, knum, vels, sen_vs, sen_vp, sen_rho, shallow, ftol, nparam):
    counts = np.zeros(fdm.shape[0], dtype=np.int64)
    for r in prange(fdm.shape[0]):
        row = np.empty(nparam, dtype=F)
        k = knum[r]
        _row(F, fdm[r], vels, sen_vs[k], sen_vp[k], sen_rho[k], shallow, ftol, row)
        n = 0
        for i in range(nparam):
            if abs(row[i]) > ftol:
                n += 1
        counts[r] = n
    return counts


@njit(cache=True, parallel=True)
def _fill(F, fdm, knum, vels, sen_vs, sen_vp, sen_rho, shallow, ftol, nparam, indptr, indices, data):
    for r in prange(fdm.shape[0]):
        row = np.empty(nparam, dtype=F)
        k = knum[r]
        _row(F, fdm[r], vels, sen_vs[k], sen_vp[k], sen_rho[k], shallow, ftol, row)
        at = indptr[r]
        for i in range(nparam):
            if abs(row[i]) > ftol:
                indices[at] = i
                data[at] = row[i]
                at += 1


def assemble(fdm, knum, vels, sen_vs, sen_vp, sen_rho, deep_relation, real=np.float32):
    """CSR rows of all data.

    fdm (rows, nvz+2, nvx+2) REAL; knum (rows,) 0-based kernel period index;
    vels (nz, ny, nx) REAL Vs; sen_* (periods, ny, nx, nz) double node kernels;
    deep_relation: True when depz(nz-1) >= 35 km (the Fortran's second branch).
    """
    real = np.dtype(real).type
    nz = vels.shape[0]
    nparam = (fdm.shape[1] - 2) * (fdm.shape[2] - 2) * (nz - 1)
    vels = np.ascontiguousarray(vels, real)
    args = (np.ascontiguousarray(fdm, real), np.ascontiguousarray(knum, np.int64), vels)
    kernels = tuple(np.ascontiguousarray(s, np.float64) for s in (sen_vs, sen_vp, sen_rho))
    shallow = not deep_relation
    ftol = real(FTOL)
    counts = _count(real, *args, *kernels, shallow, ftol, nparam)
    indptr = np.zeros(len(counts) + 1, dtype=np.int64)
    np.cumsum(counts, out=indptr[1:])
    indices = np.empty(indptr[-1], dtype=np.int32)
    data = np.empty(indptr[-1], dtype=real)
    _fill(real, *args, *kernels, shallow, ftol, nparam, indptr, indices, data)
    return csr_matrix((data, indices, indptr), shape=(len(counts), nparam))
