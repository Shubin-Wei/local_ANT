"""Ray paths and B-spline Frechet weights of DSurfTomo's CalSurfG (rpaths).

Numba kernel mirroring DSurfTomo a651faf (MIT), CalSurfG.f90 subroutine
``rpaths`` with asgr=1. The ray steps from the receiver down the first-order
traveltime gradient (refined-grid gradient while all four refined nodes are
alive) with step dpl = half the smallest coarse node spacing. Behaviour kept
from the Fortran, including what may look like defects:

* if the receiver is within 2*dpl of the source, or the ray starts in the
  refined source cell, the loop exits before any derivative is added (zero row);
* the closing segment from the last ray point to the source is not integrated;
* a segment is split at B-spline cell faces (vrat/chp); on a later sub-segment
  the previous point's basis values, computed in the old cell, are added to the
  new cell's vertices;
* end points leaving the model are clamped to its edge after the source test.

``fdm`` is (nvz + 2, nvx + 2), indices 0..nvz+1 as in the Fortran.
"""

import math

import numpy as np
from numba import njit

from .dsurf_bspline import basis

RECEIVER_OUTSIDE = 2


@njit(cache=True)
def _velocity(F, rx, rz, ipx, ipz, veln, nnx, nnz, gox, goz, dnx, dnz):
    drx = (rx - gox) - F(ipx - 1) * dnx
    drz = (rz - goz) - F(ipz - 1) * dnz
    vel = F(0.0)
    for l in range(1, 3):
        for m in range(1, 3):
            produ = F(1.0) - abs((F(m - 1) * dnz - drz) / dnz)
            produ = produ * (F(1.0) - abs((F(l - 1) * dnx - drx) / dnx))
            if ipz - 1 + m <= nnz and ipx - 1 + l <= nnx:
                vel = vel + veln[ipz - 1 + m, ipx - 1 + l] * produ
    return vel


@njit(cache=True)
def _refined_cell(rx, rz, goxr, gozr, dnxr, dnzr, nstsr):
    nnzr = nstsr.shape[0] - 1
    nnxr = nstsr.shape[1] - 1
    ipxr = int((rx - goxr) / dnxr) + 1
    ipzr = int((rz - gozr) / dnzr) + 1
    igref = 1
    if ipxr < 1 or ipxr >= nnxr:
        igref = 0
    if ipzr < 1 or ipzr >= nnzr:
        igref = 0
    if igref == 1:
        if nstsr[ipzr, ipxr] != 0 or nstsr[ipzr + 1, ipxr] != 0:
            igref = 0
        if nstsr[ipzr, ipxr + 1] != 0 or nstsr[ipzr + 1, ipxr + 1] != 0:
            igref = 0
    return ipxr, ipzr, igref


@njit(cache=True)
def _distance(F, scx, scz, rx, rz, earth):
    t = (scx - rx) * earth
    sred = t * t
    t = (scz - rz) * earth * math.sin(rx)
    sred = sred + t * t
    return math.sqrt(sred)


@njit(cache=True)
def rpaths(
    F,
    scx,
    scz,
    rcx,
    rcz,
    fdm,
    ttn,
    veln,
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
):
    """Accumulate d(time)/d(velocity vertex) of one receiver into fdm (zeroed here). Returns status."""
    nnz = veln.shape[0] - 1
    nnx = veln.shape[1] - 1
    isx = int((scx - goxr) / dnxr) + 1
    isz = int((scz - gozr) / dnzr) + 1
    dpl = dnx * earth
    rd1 = dnz * earth * math.sin(gox)
    if rd1 < dpl:
        dpl = rd1
    rd1 = dnz * earth * math.sin(gox + F(nnx - 1) * dnx)
    if rd1 < dpl:
        dpl = rd1
    dpl = F(0.5) * dpl
    fdm[:, :] = F(0.0)
    ipx = int((rcx - gox) / dnx) + 1
    ipz = int((rcz - goz) / dnz) + 1
    if ipx < 1 or ipx >= nnx or ipz < 1 or ipz >= nnz:
        return RECEIVER_OUTSIDE
    rgx = rcx
    rgz = rcz
    sw = 0
    if _distance(F, scx, scz, rgx, rgz, earth) < F(2.0) * dpl:
        sw = 1
    ipxr, ipzr, igref = _refined_cell(rcx, rcz, goxr, gozr, dnxr, dnzr, nstsr)
    if sw == 0 and igref == 1 and ipxr == isx and ipzr == isz:
        sw = 1
    vrat = np.zeros(5, dtype=F)
    chp = np.zeros(5, dtype=np.int64)
    vi = np.zeros(4, dtype=F)
    wi = np.zeros(4, dtype=F)
    vio = np.zeros(4, dtype=F)
    wio = np.zeros(4, dtype=F)
    two_earth = F(2.0) * earth
    for _ in range(nnx * nnz):
        if sw == 1:
            break
        if igref == 1:
            dtx = ttnr[ipzr, ipxr + 1] - ttnr[ipzr, ipxr]
            dtx = dtx + ttnr[ipzr + 1, ipxr + 1] - ttnr[ipzr + 1, ipxr]
            dtx = dtx / (two_earth * dnxr)
            dtz = ttnr[ipzr + 1, ipxr] - ttnr[ipzr, ipxr]
            dtz = dtz + ttnr[ipzr + 1, ipxr + 1] - ttnr[ipzr, ipxr + 1]
            dtz = dtz / (two_earth * math.sin(rgx) * dnzr)
        else:
            dtx = ttn[ipz, ipx + 1] - ttn[ipz, ipx]
            dtx = dtx + ttn[ipz + 1, ipx + 1] - ttn[ipz + 1, ipx]
            dtx = dtx / (two_earth * dnx)
            dtz = ttn[ipz + 1, ipx] - ttn[ipz, ipx]
            dtz = dtz + ttn[ipz + 1, ipx + 1] - ttn[ipz, ipx + 1]
            dtz = dtz / (two_earth * math.sin(rgx) * dnz)
        rd1 = math.sqrt(dtx * dtx + dtz * dtz)
        nx_ = rgx - dpl * dtx / (earth * rd1)
        nz_ = rgz - dpl * dtz / (earth * math.sin(rgx) * rd1)
        ipxo = ipx
        ipzo = ipz
        ipxr, ipzr, igref = _refined_cell(nx_, nz_, goxr, gozr, dnxr, dnzr, nstsr)
        ipx = int((nx_ - gox) / dnx) + 1
        ipz = int((nz_ - goz) / dnz) + 1
        sw = 0
        if _distance(F, scx, scz, nx_, nz_, earth) < F(2.0) * dpl:
            sw = 1
        if sw == 0 and igref == 1 and ipxr == isx and ipzr == isz:
            sw = 1
        if ipx < 1:
            nx_ = gox
            ipx = 1
        if ipx >= nnx:
            nx_ = gox + F(nnx - 1) * dnx
            ipx = nnx - 1
        if ipz < 1:
            nz_ = goz
            ipz = 1
        if ipz >= nnz:
            nz_ = goz + F(nnz - 1) * dnz
            ipz = nnz - 1
        _segment(
            F,
            rgx,
            rgz,
            nx_,
            nz_,
            ipx,
            ipz,
            ipxo,
            ipzo,
            fdm,
            veln,
            nnx,
            nnz,
            gox,
            goz,
            dnx,
            dnz,
            dvx,
            dvz,
            gdx,
            gdz,
            dpl,
            vrat,
            chp,
            vi,
            wi,
            vio,
            wio,
        )
        rgx = nx_
        rgz = nz_
    return 0


@njit(cache=True)
def _segment(
    F,
    rx0,
    rz0,
    rx1,
    rz1,
    ipx,
    ipz,
    ipxo,
    ipzo,
    fdm,
    veln,
    nnx,
    nnz,
    gox,
    goz,
    dnx,
    dnz,
    dvx,
    dvz,
    gdx,
    gdz,
    dpl,
    vrat,
    chp,
    vi,
    wi,
    vio,
    wio,
):
    """Trapezoidal 1/v^2 B-spline weights of one ray step, split at B-spline cell faces."""
    ivx = (ipx - 1) // gdx + 1
    ivz = (ipz - 1) // gdz + 1
    ivxo = (ipxo - 1) // gdx + 1
    ivzo = (ipzo - 1) // gdz + 1
    nhp = 0
    if ivx != ivxo:
        nhp = nhp + 1
        if ivx > ivxo:
            xi = gox + F(ivx - 1) * dvx
        else:
            xi = gox + F(ivx) * dvx
        vrat[nhp] = (xi - rx0) / (rx1 - rx0)
        chp[nhp] = 1
    if ivz != ivzo:
        nhp = nhp + 1
        if ivz > ivzo:
            zi = goz + F(ivz - 1) * dvz
        else:
            zi = goz + F(ivz) * dvz
        rd1 = (zi - rz0) / (rz1 - rz0)
        if nhp == 1:
            vrat[nhp] = rd1
            chp[nhp] = 2
        elif rd1 >= vrat[nhp - 1]:
            vrat[nhp] = rd1
            chp[nhp] = 2
        else:
            vrat[nhp] = vrat[nhp - 1]
            chp[nhp] = chp[nhp - 1]
            vrat[nhp - 1] = rd1
            chp[nhp - 1] = 2
    nhp = nhp + 1
    vrat[nhp] = F(1.0)
    chp[nhp] = 0
    vel = _velocity(F, rx0, rz0, ipxo, ipzo, veln, nnx, nnz, gox, goz, dnx, dnz)
    v = ((rx0 - gox) - F(ivxo - 1) * dvx) / dvx
    w = ((rz0 - goz) - F(ivzo - 1) * dvz) / dvz
    basis(F, v, vi)
    basis(F, w, wi)
    ivxt = ivxo
    ivzt = ivzo
    for k in range(1, nhp + 1):
        velo = vel
        vio[:] = vi
        wio[:] = wi
        if k > 1:
            if chp[k - 1] == 1:
                ivxt = ivx
            elif chp[k - 1] == 2:
                ivzt = ivz
        rigz = rz0 + vrat[k] * (rz1 - rz0)
        rigx = rx0 + vrat[k] * (rx1 - rx0)
        ipxt = int((rigx - gox) / dnx) + 1
        ipzt = int((rigz - goz) / dnz) + 1
        vel = _velocity(F, rigx, rigz, ipxt, ipzt, veln, nnx, nnz, gox, goz, dnx, dnz)
        v = ((rigx - gox) - F(ivxt - 1) * dvx) / dvx
        w = ((rigz - goz) - F(ivzt - 1) * dvz) / dvz
        basis(F, v, vi)
        basis(F, w, wi)
        if k == 1:
            dinc = vrat[k] * dpl
        else:
            dinc = (vrat[k] - vrat[k - 1]) * dpl
        for l in range(1, 5):
            for m in range(1, 5):
                rd1 = vi[m - 1] * wi[l - 1] / (vel * vel)
                rd2 = vio[m - 1] * wio[l - 1] / (velo * velo)
                rd1 = -(rd1 + rd2) * dinc / F(2.0)
                fdm[ivzt - 2 + l, ivxt - 2 + m] = rd1 + fdm[ivzt - 2 + l, ivxt - 2 + m]
