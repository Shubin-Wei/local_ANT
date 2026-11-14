"""Cubic B-spline velocity grids of DSurfTomo's CalSurfG (gridder, bsplrefine).

Numba kernels mirroring DSurfTomo a651faf (MIT), CalSurfG.f90 subroutines
``gridder`` and ``bsplrefine``: the same basis polynomials, summation order and
node ranges. Every kernel takes the real type ``F`` (``np.float32`` reproduces
the Fortran REAL arithmetic; ``np.float64`` is the same algorithm in double).
Arrays use the Fortran 1-based indices: element 0 of a node axis is padding.
"""

import numpy as np
from numba import njit


@njit(cache=True)
def basis(F, u, out):
    """The four cubic B-spline weights at local coordinate u (written to out[0:4])."""
    a = F(1.0) - u
    u2 = u * u
    u3 = u2 * u
    out[0] = a * a * a / F(6.0)
    out[1] = (F(4.0) - F(6.0) * u2 + F(3.0) * u3) / F(6.0)
    out[2] = (F(1.0) + F(3.0) * u + F(3.0) * u2 - F(3.0) * u3) / F(6.0)
    out[3] = u3 / F(6.0)


@njit(cache=True)
def _basis_table(F, n):
    table = np.empty((n + 1, 4), dtype=F)
    for i in range(1, n + 2):
        u = F(n)
        u = F(i - 1) / u
        basis(F, u, table[i - 1])
    return table


@njit(cache=True)
def gridder(F, velv, gdx, gdz):
    """velv (0..nvz+1, 0..nvx+1) control values -> veln (1..nnz, 1..nnx) propagation nodes."""
    nvz = velv.shape[0] - 2
    nvx = velv.shape[1] - 2
    nnx = (nvx - 1) * gdx + 1
    nnz = (nvz - 1) * gdz + 1
    veln = np.zeros((nnz + 1, nnx + 1), dtype=F)
    ui = _basis_table(F, gdx)
    vi = _basis_table(F, gdz)
    for i in range(1, nvz):
        conz = gdz + 1 if i == nvz - 1 else gdz
        for j in range(1, nvx):
            conx = gdx + 1 if j == nvx - 1 else gdx
            for l in range(1, conz + 1):
                stz = gdz * (i - 1) + l
                for m in range(1, conx + 1):
                    stx = gdx * (j - 1) + m
                    sumi = F(0.0)
                    for i1 in range(1, 5):
                        sumj = F(0.0)
                        for j1 in range(1, 5):
                            sumj = sumj + ui[m - 1, j1 - 1] * velv[i - 2 + i1, j - 2 + j1]
                        sumi = sumi + vi[l - 1, i1 - 1] * sumj
                    veln[stz, stx] = sumi
    return veln


@njit(cache=True)
def bsplrefine(F, velv, gdx, gdz, sgdl, vnl, vnr, vnt, vnb):
    """Velocity on the refined source grid spanning coarse nodes vnl..vnr (x), vnt..vnb (z)."""
    nvz = velv.shape[0] - 2
    nvx = velv.shape[1] - 2
    nnx = (vnr - vnl) * sgdl + 1
    nnz = (vnb - vnt) * sgdl + 1
    veln = np.zeros((nnz + 1, nnx + 1), dtype=F)
    nrxr = gdx * sgdl
    nrzr = gdz * sgdl
    ui = _basis_table(F, nrxr)  # depends on the x position only
    vi = _basis_table(F, nrzr)  # depends on the z position only
    origx = (vnl - 1) * sgdl + 1
    origz = (vnt - 1) * sgdl + 1
    total = np.empty(4, dtype=F)
    for i in range(1, nvz):
        conz = nrzr + 1 if i == nvz - 1 else nrzr
        for j in range(1, nvx):
            conx = nrxr + 1 if j == nvx - 1 else nrxr
            for k in range(1, conz + 1):
                st1 = gdz * (i - 1) + (k - 1) // sgdl + 1
                if st1 < vnt or st1 > vnb:
                    continue
                st1 = nrzr * (i - 1) + k
                for l in range(1, conx + 1):
                    st2 = gdx * (j - 1) + (l - 1) // sgdl + 1
                    if st2 < vnl or st2 > vnr:
                        continue
                    st2 = nrxr * (j - 1) + l
                    for i1 in range(1, 5):
                        s = F(0.0)
                        for j1 in range(1, 5):
                            s = s + ui[l - 1, j1 - 1] * velv[i - 2 + i1, j - 2 + j1]
                        total[i1 - 1] = vi[k - 1, i1 - 1] * s
                    idm1 = st1 - origz + 1
                    idm2 = st2 - origx + 1
                    if idm1 < 1 or idm1 > nnz or idm2 < 1 or idm2 > nnx:
                        continue
                    veln[idm1, idm2] = total[0] + total[1] + total[2] + total[3]
    return veln
