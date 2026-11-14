"""Fast marching on DSurfTomo's spherical propagation grid (module traveltime).

Numba kernels mirroring DSurfTomo a651faf (MIT), CalSurfG.f90: ``travel``,
``fouds2`` (mixed-order upwind stencil; CalSurfG fixes fom=1, so the unused
first-order ``fouds1`` is not ported), the binary-tree heap ``addtree``,
``downtree``, ``updtree``, and ``bilinear``. Tie order, strict comparisons,
the refined-grid exit test and the source-cell initialisation follow the
Fortran line by line. ``F`` selects the real type (np.float32 = Fortran REAL).

Arrays use Fortran 1-based indices with element 0 as padding: ttn, nsts and
veln are (nnz + 1, nnx + 1); the heap arrays btx/btz hold grid indices from
position 1. nsts: -1 far, 0 alive, > 0 position in the heap.
"""

import math

import numpy as np
from numba import njit

# Return codes of travel (0 = success).
SOURCE_OUTSIDE = 1


@njit(cache=True)
def bilinear(F, nv, dsx, dsz, dnx, dnz):
    """nv[i, j] (i = x, j = z, 1-based 2 x 2) at offset (dsx, dsz) inside the cell."""
    biv = F(0.0)
    for i in range(1, 3):
        for j in range(1, 3):
            produ = (F(1.0) - abs((F(i - 1) * dnx - dsx) / dnx)) * (F(1.0) - abs((F(j - 1) * dnz - dsz) / dnz))
            biv = biv + nv[i, j] * produ
    return biv


@njit(cache=True)
def addtree(iz, ix, ntr, ttn, nsts, btx, btz):
    ntr = ntr + 1
    nsts[iz, ix] = ntr
    btx[ntr] = ix
    btz[ntr] = iz
    tpc = ntr
    tpp = tpc // 2
    while tpp > 0:
        if ttn[iz, ix] < ttn[btz[tpp], btx[tpp]]:
            nsts[iz, ix] = tpp
            nsts[btz[tpp], btx[tpp]] = tpc
            _swap(btx, btz, tpc, tpp)
            tpc = tpp
            tpp = tpc // 2
        else:
            tpp = 0
    return ntr


@njit(cache=True)
def _swap(btx, btz, a, b):
    x, z = btx[a], btz[a]
    btx[a], btz[a] = btx[b], btz[b]
    btx[b], btz[b] = x, z


@njit(cache=True)
def downtree(ntr, ttn, nsts, btx, btz):
    """Remove the root (already marked alive by the caller); return the new size."""
    if ntr == 1:
        return 0
    nsts[btz[ntr], btx[ntr]] = 1
    btx[1] = btx[ntr]
    btz[1] = btz[ntr]
    ntr = ntr - 1
    tpp = 1
    tpc = 2 * tpp
    while tpc < ntr:
        rd1 = ttn[btz[tpc], btx[tpc]]
        rd2 = ttn[btz[tpc + 1], btx[tpc + 1]]
        if rd1 > rd2:
            tpc = tpc + 1
        rd1 = ttn[btz[tpc], btx[tpc]]
        rd2 = ttn[btz[tpp], btx[tpp]]
        if rd1 < rd2:
            nsts[btz[tpp], btx[tpp]] = tpc
            nsts[btz[tpc], btx[tpc]] = tpp
            _swap(btx, btz, tpc, tpp)
            tpp = tpc
            tpc = 2 * tpp
        else:
            tpc = ntr + 1
    if tpc == ntr:
        rd1 = ttn[btz[tpc], btx[tpc]]
        rd2 = ttn[btz[tpp], btx[tpp]]
        if rd1 < rd2:
            nsts[btz[tpp], btx[tpp]] = tpc
            nsts[btz[tpc], btx[tpc]] = tpp
            _swap(btx, btz, tpc, tpp)
    return ntr


@njit(cache=True)
def updtree(iz, ix, ttn, nsts, btx, btz):
    tpc = nsts[iz, ix]
    tpp = tpc // 2
    while tpp > 0:
        if ttn[iz, ix] < ttn[btz[tpp], btx[tpp]]:
            nsts[iz, ix] = tpp
            nsts[btz[tpp], btx[tpp]] = tpc
            _swap(btx, btz, tpc, tpp)
            tpc = tpp
            tpp = tpc // 2
        else:
            tpp = 0


@njit(cache=True)
def fouds2(F, iz, ix, ttn, nsts, veln, nnx, nnz, gox, dnx, dnz, earth):
    """Mixed-order trial time at (iz, ix) from its alive neighbours (Popovici & Sethian 2002)."""
    tsw1 = 0
    travm = F(0.0)
    slown = F(1.0) / veln[iz, ix]
    ri = earth
    risti = ri * math.sin(gox + F(ix - 1) * dnx)
    for j in (ix - 1, ix + 1):
        if j < 1 or j > nnx:
            continue
        swj = -1
        if j == ix - 1:
            j2 = j - 1
            if j2 >= 1:
                if nsts[iz, j2] == 0:
                    swj = 0
        else:
            j2 = j + 1
            if j2 <= nnx:
                if nsts[iz, j2] == 0:
                    swj = 0
        if nsts[iz, j] == 0 and swj == 0:
            swj = -1
            if ttn[iz, j] > ttn[iz, j2]:
                swj = 0
        else:
            swj = -1
        for k in (iz - 1, iz + 1):
            if k < 1 or k > nnz:
                continue
            swk = -1
            if k == iz - 1:
                k2 = k - 1
                if k2 >= 1:
                    if nsts[k2, ix] == 0:
                        swk = 0
            else:
                k2 = k + 1
                if k2 <= nnz:
                    if nsts[k2, ix] == 0:
                        swk = 0
            if nsts[k, ix] == 0 and swk == 0:
                swk = -1
                if ttn[k, ix] > ttn[k2, ix]:
                    swk = 0
            else:
                swk = -1
            swsol = 0
            a = F(1.0)
            b = F(0.0)
            c = F(0.0)
            tref = F(0.0)
            tdiv = F(1.0)
            if swj == 0:
                swsol = 1
                if swk == 0:
                    u = F(2.0) * ri * dnx
                    v = F(2.0) * risti * dnz
                    em = F(4.0) * ttn[iz, j] - ttn[iz, j2] - F(4.0) * ttn[k, ix]
                    em = em + ttn[k2, ix]
                    a = v * v + u * u
                    b = F(2.0) * em * (u * u)
                    c = (u * u) * (em * em - (slown * slown) * (v * v))
                    tref = F(4.0) * ttn[iz, j] - ttn[iz, j2]
                    tdiv = F(3.0)
                elif nsts[k, ix] == 0:
                    u = risti * dnz
                    v = F(2.0) * ri * dnx
                    em = F(3.0) * ttn[k, ix] - F(4.0) * ttn[iz, j] + ttn[iz, j2]
                    a = v * v + F(9.0) * (u * u)
                    b = F(6.0) * em * (u * u)
                    c = (u * u) * (em * em - (slown * slown) * (v * v))
                    tref = ttn[k, ix]
                    tdiv = F(1.0)
                else:
                    u = F(2.0) * ri * dnx
                    a = F(1.0)
                    b = F(0.0)
                    c = -((u * u) * (slown * slown))
                    tref = F(4.0) * ttn[iz, j] - ttn[iz, j2]
                    tdiv = F(3.0)
            elif nsts[iz, j] == 0:
                swsol = 1
                if swk == 0:
                    u = ri * dnx
                    v = F(2.0) * risti * dnz
                    em = F(3.0) * ttn[iz, j] - F(4.0) * ttn[k, ix] + ttn[k2, ix]
                    a = v * v + F(9.0) * (u * u)
                    b = F(6.0) * em * (u * u)
                    c = (u * u) * (em * em - (v * v) * (slown * slown))
                    tref = ttn[iz, j]
                    tdiv = F(1.0)
                elif nsts[k, ix] == 0:
                    u = ri * dnx
                    v = risti * dnz
                    em = ttn[k, ix] - ttn[iz, j]
                    a = u * u + v * v
                    b = -(F(2.0) * (u * u) * em)
                    c = (u * u) * (em * em - (v * v) * (slown * slown))
                    tref = ttn[iz, j]
                    tdiv = F(1.0)
                else:
                    a = F(1.0)
                    b = F(0.0)
                    c = -((slown * slown) * (ri * ri) * (dnx * dnx))
                    tref = ttn[iz, j]
                    tdiv = F(1.0)
            else:
                if swk == 0:
                    swsol = 1
                    u = F(2.0) * risti * dnz
                    a = F(1.0)
                    b = F(0.0)
                    c = -((u * u) * (slown * slown))
                    tref = F(4.0) * ttn[k, ix] - ttn[k2, ix]
                    tdiv = F(3.0)
                elif nsts[k, ix] == 0:
                    swsol = 1
                    a = F(1.0)
                    b = F(0.0)
                    c = -((slown * slown) * (risti * risti) * (dnz * dnz))
                    tref = ttn[k, ix]
                    tdiv = F(1.0)
            if swsol == 1:
                rd1 = b * b - F(4.0) * a * c
                if rd1 < F(0.0):
                    rd1 = F(0.0)
                tdsh = (-b + math.sqrt(rd1)) / (F(2.0) * a)
                trav = (tref + tdsh) / tdiv
                if tsw1 == 1:
                    travm = min(trav, travm)
                else:
                    travm = trav
                    tsw1 = 1
    ttn[iz, ix] = travm


@njit(cache=True)
def _visit(F, iz, ix, ntr, ttn, nsts, btx, btz, veln, nnx, nnz, gox, dnx, dnz, earth):
    if nsts[iz, ix] == -1:
        fouds2(F, iz, ix, ttn, nsts, veln, nnx, nnz, gox, dnx, dnz, earth)
        ntr = addtree(iz, ix, ntr, ttn, nsts, btx, btz)
    elif nsts[iz, ix] > 0:
        fouds2(F, iz, ix, ttn, nsts, veln, nnx, nnz, gox, dnx, dnz, earth)
        updtree(iz, ix, ttn, nsts, btx, btz)
    return ntr


@njit(cache=True)
def _latest(need, ttn):
    t = -np.inf
    for j in range(need.shape[0]):
        for i in range(need.shape[1]):
            if need[j, i] and ttn[j, i] > t:
                t = ttn[j, i]
    return t


@njit(cache=True)
def travel(
    F, scx, scz, urg, veln, ttn, nsts, btx, btz, gox, goz, dnx, dnz, earth, vnl, vnr, vnt, vnb, need=None, margin=0.0
):
    """First-arrival field from source (scx, scz); urg 0 plain, 1 refined grid, 2 resume.

    With urg=1 the march stops when the root reaches an edge of the refined grid
    that is not an edge of the model. As in the Fortran, the test compares vnr
    and vnb (coarse node indices) with the node counts of the grid in use.

    need (lant, not in the Fortran): nodes whose final times are read afterwards (the corners of the receivers'
    cells).  Once all of them are alive, the march goes on only while the root is not later than their latest time
    plus `margin` (s), and stops there.  The nodes it leaves undone are later than every node a ray from these
    receivers reads (rays descend in time from the receiver cell; margin covers the cells they cross), and the
    nodes it does are marched exactly as in the full march, so every value read is the same.
    """
    nnz = veln.shape[0] - 1
    nnx = veln.shape[1] - 1
    isx = int((scx - gox) / dnx) + 1
    isz = int((scz - goz) / dnz) + 1
    if isx < 1 or isx > nnx or isz < 1 or isz > nnz:
        return SOURCE_OUTSIDE
    if isx == nnx:
        isx = isx - 1
    if isz == nnz:
        isz = isz - 1
    ntr = 0
    if urg == 2:
        for i in range(1, nnx + 1):
            for j in range(1, nnz + 1):
                if nsts[j, i] > 0:
                    ntr = addtree(j, i, ntr, ttn, nsts, btx, btz)
    else:
        nsts[:, :] = -1
        vss = np.empty((3, 3), dtype=F)
        for i in range(1, 3):
            for j in range(1, 3):
                vss[i, j] = veln[isz - 1 + j, isx - 1 + i]
        dsx = (scx - gox) - F(isx - 1) * dnx
        dsz = (scz - goz) - F(isz - 1) * dnz
        vsrc = bilinear(F, vss, dsx, dsz, dnx, dnz)
        for i in range(1, 3):
            for j in range(1, 3):
                ex = dsx - F(i - 1) * dnx
                ez = dsz - F(j - 1) * dnz
                ds = math.sqrt(ex * ex + ez * ez)
                ttn[isz - 1 + j, isx - 1 + i] = F(2.0) * ds / (vss[i, j] + vsrc)
                ntr = addtree(isz - 1 + j, isx - 1 + i, ntr, ttn, nsts, btx, btz)
    remaining = 0
    cut = np.inf
    if need is not None:
        for j in range(1, nnz + 1):
            for i in range(1, nnx + 1):
                if need[j, i] and nsts[j, i] != 0:
                    remaining += 1
        if remaining == 0:
            cut = _latest(need, ttn) + margin
    while ntr > 0:
        if need is not None and ttn[btz[1], btx[1]] > cut:
            break
        if urg == 1:
            ix = btx[1]
            iz = btz[1]
            swrg = 0
            if ix == 1 and vnl != 1:
                swrg = 1
            if ix == nnx and vnr != nnx:
                swrg = 1
            if iz == 1 and vnt != 1:
                swrg = 1
            if iz == nnz and vnb != nnz:
                swrg = 1
            if swrg == 1:
                nsts[iz, ix] = 0
                break
        ix = btx[1]
        iz = btz[1]
        nsts[iz, ix] = 0
        ntr = downtree(ntr, ttn, nsts, btx, btz)
        if need is not None and need[iz, ix]:
            remaining -= 1
            if remaining == 0:
                cut = _latest(need, ttn) + margin
        for i in (ix - 1, ix + 1):
            if 1 <= i <= nnx:
                ntr = _visit(F, iz, i, ntr, ttn, nsts, btx, btz, veln, nnx, nnz, gox, dnx, dnz, earth)
        for i in (iz - 1, iz + 1):
            if 1 <= i <= nnz:
                ntr = _visit(F, i, ix, ntr, ttn, nsts, btx, btz, veln, nnx, nnz, gox, dnx, dnz, earth)
    return 0
