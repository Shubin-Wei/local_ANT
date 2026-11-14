"""Source-grid refinement and receiver times of DSurfTomo's CalSurfG.

Numba kernels mirroring DSurfTomo a651faf (MIT), CalSurfG.f90: the asgr=1
block of the CalSurfG source loop (refined box of sgs coarse nodes around the
source, bsplrefine, refined march, resampling onto the coarse grid, narrow-band
completion, resumed coarse march) and ``srtimes`` (bilinear receiver time;
near the source 2d/(vs+vr) with d the local flat-earth distance). Arrays use
Fortran 1-based indices with element 0 as padding; ``F`` is the real type.
"""

import math

import numpy as np
from numba import njit

from .dsurf_bspline import bsplrefine
from .dsurf_fmm import bilinear, travel

RECEIVER_OUTSIDE = 2
SOURCE_CELL_OUTSIDE = 3


@njit(cache=True)
def source_fields(F, x, z, velv, veln, gox, goz, dvx, dvz, earth, gdx, gdz, sgdl, sgs, rx=None, rz=None):
    """Coarse and refined first-arrival fields of one source.

    Returns (status, ttn, ttnr, nstsr, goxr, gozr, dnxr, dnzr). The refined
    arrays are what rpaths later reads as ttnr/nstsr.

    rx, rz (lant): the receivers of this source.  The coarse march then stops once it has done every node that
    srtimes and rpaths read for them (travel's `need`, with a margin of two cell diagonals at the slowest
    velocity); the coarse field is not complete beyond, and must not be read for other receivers.
    """
    nnz = veln.shape[0] - 1
    nnx = veln.shape[1] - 1
    dnx = dvx / F(gdx)
    dnz = dvz / F(gdz)
    isx = int((x - gox) / dnx) + 1
    isz = int((z - goz) / dnz) + 1
    empty = np.zeros((1, 1), dtype=F)
    empty_i = np.zeros((1, 1), dtype=np.int64)
    if isx < 1 or isx > nnx or isz < 1 or isz > nnz:
        return 1, empty, empty, empty_i, gox, goz, dnx, dnz
    if isx == nnx:
        isx = isx - 1
    if isz == nnz:
        isz = isz - 1
    vnl = max(isx - sgs, 1)
    vnr = min(isx + sgs, nnx)
    vnt = max(isz - sgs, 1)
    vnb = min(isz + sgs, nnz)
    nrnx = (vnr - vnl) * sgdl + 1
    nrnz = (vnb - vnt) * sgdl + 1
    drnx = dvx / F(gdx * sgdl)
    drnz = dvz / F(gdz * sgdl)
    gorx = gox + dnx * F(vnl - 1)
    gorz = goz + dnz * F(vnt - 1)
    velr = bsplrefine(F, velv, gdx, gdz, sgdl, vnl, vnr, vnt, vnb)
    ttnr = np.zeros((nrnz + 1, nrnx + 1), dtype=F)
    nstsr = np.full((nrnz + 1, nrnx + 1), -1, dtype=np.int64)
    capacity = (nrnz + 1) * (nrnx + 1) + 1
    btx = np.zeros(max(capacity, (nnz + 1) * (nnx + 1) + 1), dtype=np.int64)
    btz = np.zeros_like(btx)
    status = travel(F, x, z, 1, velr, ttnr, nstsr, btx, btz, gorx, gorz, drnx, drnz, earth, vnl, vnr, vnt, vnb)
    if status != 0:
        return status, empty, empty, empty_i, gox, goz, dnx, dnz
    # Resample alive and close refined nodes onto the coarse grid.
    ttn = np.zeros((nnz + 1, nnx + 1), dtype=F)
    nsts = np.full((nnz + 1, nnx + 1), -1, dtype=np.int64)
    for k in range(1, nrnz + 1, sgdl):
        idm1 = vnt + (k - 1) // sgdl
        for l in range(1, nrnx + 1, sgdl):
            idm2 = vnl + (l - 1) // sgdl
            nsts[idm1, idm2] = nstsr[k, l]
            if nsts[idm1, idm2] >= 0:
                ttn[idm1, idm2] = ttnr[k, l]
    # Alive nodes next to far nodes become close again.
    for k in range(1, nnx + 1):
        for l in range(1, nnz + 1):
            if nsts[l, k] == 0:
                if l - 1 >= 1 and nsts[l - 1, k] == -1:
                    nsts[l, k] = 1
                if l + 1 <= nnz and nsts[l + 1, k] == -1:
                    nsts[l, k] = 1
                if k - 1 >= 1 and nsts[l, k - 1] == -1:
                    nsts[l, k] = 1
                if k + 1 <= nnx and nsts[l, k + 1] == -1:
                    nsts[l, k] = 1
    if rx is None:
        status = travel(F, x, z, 2, veln, ttn, nsts, btx, btz, gox, goz, dnx, dnz, earth, vnl, vnr, vnt, vnb)
    else:
        need = np.zeros((nnz + 1, nnx + 1), dtype=np.bool_)
        for r in range(rx.size):
            irx = int((rx[r] - gox) / dnx) + 1
            irz = int((rz[r] - goz) / dnz) + 1
            if irx < 1 or irx > nnx or irz < 1 or irz > nnz:
                continue  # srtimes reports it
            if irx == nnx:
                irx = irx - 1
            if irz == nnz:
                irz = irz - 1
            need[irz, irx] = need[irz + 1, irx] = need[irz, irx + 1] = need[irz + 1, irx + 1] = True
        vmin = np.inf
        for j in range(1, nnz + 1):
            for i in range(1, nnx + 1):
                if veln[j, i] < vmin:
                    vmin = veln[j, i]
        margin = 2.0 * float(earth) * math.sqrt(float(dnx) ** 2 + float(dnz) ** 2) / vmin
        status = travel(F, x, z, 2, veln, ttn, nsts, btx, btz, gox, goz, dnx, dnz, earth, vnl, vnr, vnt, vnb, need,
                        margin)  # fmt: skip
    return status, ttn, ttnr, nstsr, gorx, gorz, drnx, drnz


@njit(cache=True)
def srtimes(F, scx, scz, rcx1, rcz1, ttn, veln, gox, goz, dnx, dnz, earth):
    """Receiver time on the coarse grid. Returns (status, time)."""
    nnz = veln.shape[0] - 1
    nnx = veln.shape[1] - 1
    irx = int((rcx1 - gox) / dnx) + 1
    irz = int((rcz1 - goz) / dnz) + 1
    if irx < 1 or irx > nnx or irz < 1 or irz > nnz:
        return RECEIVER_OUTSIDE, F(0.0)
    sw = 0
    if irx == nnx:
        irx = irx - 1
    if irz == nnz:
        irz = irz - 1
    isx = int((scx - gox) / dnx) + 1
    isz = int((scz - goz) / dnz) + 1
    dpl = dnx * earth
    rd1 = dnz * earth * math.sin(gox)
    if rd1 < dpl:
        dpl = rd1
    rd1 = dnz * earth * math.sin(gox + F(nnx - 1) * dnx)
    if rd1 < dpl:
        dpl = rd1
    t = (scx - rcx1) * earth
    sred = t * t
    t = (scz - rcz1) * earth * math.sin(rcx1)
    sred = sred + t * t
    sred = math.sqrt(sred)
    if sred < dpl:
        sw = 1
    if isx == irx and isz == irz:
        sw = 1
    vss = np.empty((3, 3), dtype=F)
    if sw == 1:
        if isx < 1 or isx + 1 > nnx or isz < 1 or isz + 1 > nnz:
            # The Fortran reads past the grid here (a bounds-check abort).
            return SOURCE_CELL_OUTSIDE, F(0.0)
        for k in range(1, 3):
            for l in range(1, 3):
                vss[k, l] = veln[isz - 1 + l, isx - 1 + k]
        drx = (scx - gox) - F(isx - 1) * dnx
        drz = (scz - goz) - F(isz - 1) * dnz
        vels = bilinear(F, vss, drx, drz, dnx, dnz)
        for k in range(1, 3):
            for l in range(1, 3):
                vss[k, l] = veln[irz - 1 + l, irx - 1 + k]
        drx = (rcx1 - gox) - F(irx - 1) * dnx
        drz = (rcz1 - goz) - F(irz - 1) * dnz
        velr = bilinear(F, vss, drx, drz, dnx, dnz)
        trr = F(2.0) * sred / (vels + velr)
    else:
        drx = (rcx1 - gox) - F(irx - 1) * dnx
        drz = (rcz1 - goz) - F(irz - 1) * dnz
        trr = F(0.0)
        for k in range(1, 3):
            for l in range(1, 3):
                produ = (F(1.0) - abs((F(l - 1) * dnz - drz) / dnz)) * (F(1.0) - abs((F(k - 1) * dnx - drx) / dnx))
                trr = trr + ttn[irz - 1 + l, irx - 1 + k] * produ
    return 0, trr
