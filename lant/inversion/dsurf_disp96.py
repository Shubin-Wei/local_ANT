"""DSurfTomo's own dispersion code, statement by statement: surfdisp96.f (Rayleigh phase velocity, spherical earth,
fundamental mode) and, from CalSurfG.f90, the Brocher (2005) materials and refineGrid2LayerMdl of depthkernel /
caldespersion.

DSurfTomo a651faf (MIT); surfdisp96.f is Herrmann's surf96 (Computer Programs in Seismology) as modified there.  The
precision of every quantity follows the Fortran: the layer model, the earth flattening, the starting velocity of the
root search and the stored phase velocity are REAL (float32, with Fortran's implicit typing: a-h, o-z REAL), the
dispersion function and the root search are double precision.  REAL constants used in double expressions keep their
single-precision value (e.g. dc = dble(0.005) = 0.004999999888...).  Integer powers are multiplied as gfortran expands
them (x**3 = x**2 * x, x**4 = x**2 * x**2, x**5 = x**3 * x**2).  Only what DSurfTomo calls is ported: iwave = 2
(Rayleigh), igr = 0 (phase), mode = 1, iflsph = 1; dltar1 (Love) and the group-velocity branch are not.

Differences, both outside normal use: where the root search fails ("improper initial value in disper"), native
DSurfTomo writes a message and phase velocity 0 for the remaining periods; here the period is returned as 0 as well
and the caller stops.  getsol keeps del1st in a SAVE variable that the OpenMP threads of depthkernel share; here it is
local to the call (the value of a serial run).
"""

import math

import numpy as np
from numba import njit

F4 = np.float32
F8 = np.float64
NL = 200  # layers of the refined model (refineGrid2LayerMdl, surfdisp96)

# REAL literals, as gfortran stores them
_0P95, _0P90, _0P01 = F4(0.95), F4(0.90), F4(0.01)
_ONE, _TWO, _FOUR, _HALF = F4(1.0), F4(2.0), F4(4.0), F4(0.5)
_RHOEXP = F4(-2.275)
_DC = F8(F4(0.005))  # dc = dble(ddc0), ddc0 = 0.005 REAL
_ONEA = F8(F4(1.5))  # onea = dble(sone0)
_ONE_D = 1.0e-2  # one = 1.0d-2
_S_D = F8(F4(0.01))  # nevill: s1 = 0.01*ss1, a REAL constant in a double expression
_TWOPI = 2.0 * 3.141592653589793
_AR = 6370.0
# CalSurfG: Brocher (2005) Vp(Vs) and density(Vp), REAL
_B = tuple(F4(x) for x in (0.9409, 2.0947, 0.8206, 0.2683, 0.0251))
_R = tuple(F4(x) for x in (1.6612, 0.4721, 0.0671, 0.0043, 0.000106))


@njit(cache=True, nogil=True)
def brocher(vs):
    """vpz, rhoz of depthkernel (REAL) from vsz."""
    n = vs.size
    vp = np.empty(n, F4)
    rho = np.empty(n, F4)
    b0, b1, b2, b3, b4 = _B
    r0, r1, r2, r3, r4 = _R
    for k in range(n):
        v = vs[k]
        v2 = v * v
        v3 = v2 * v
        v4 = v2 * v2
        p = b0 + b1 * v - b2 * v2 + b3 * v3 - b4 * v4
        vp[k] = p
        p2 = p * p
        p3 = p2 * p
        p4 = p2 * p2
        p5 = p3 * p2
        rho[k] = r0 * p - r1 * p2 + r2 * p3 - r3 * p4 + r4 * p5
    return vp, rho


@njit(cache=True, nogil=True)
def refine(minthk0, dep, vp, vs, rho):
    """refineGrid2LayerMdl -> (rmax, rthk, rvp, rvs, rrho), REAL arrays of length NL."""
    rthk = np.zeros(NL, F4)
    rvp = np.zeros(NL, F4)
    rvs = np.zeros(NL, F4)
    rrho = np.zeros(NL, F4)
    mmax = dep.size
    k = 0
    for i in range(mmax - 1):
        thk = dep[i + 1] - dep[i]
        minthk = thk / minthk0
        ns = int((thk + F4(1.0e-4)) / minthk) + 1
        newthk = thk / F4(ns)
        two_ns = F4(2 * ns)
        for j in range(1, ns + 1):
            m = F4(2 * j - 1)
            rthk[k] = newthk
            rvp[k] = vp[i] + m * (vp[i + 1] - vp[i]) / two_ns
            rvs[k] = vs[i] + m * (vs[i + 1] - vs[i]) / two_ns
            rrho[k] = rho[i] + m * (rho[i + 1] - rho[i]) / two_ns
            k += 1
    rthk[k] = F4(0.0)
    rvp[k] = vp[mmax - 1]
    rvs[k] = vs[mmax - 1]
    rrho[k] = rho[mmax - 1]
    return k + 1, rthk, rvp, rvs, rrho


@njit(cache=True, nogil=True)
def _sign(x):
    return math.copysign(1.0, x)  # dsign(1.0d0, x)


@njit(cache=True, nogil=True)
def _var(p, q, ra, rb, wvno, xka, xkb, dpth):
    exa = 0.0
    a0 = 1.0
    pex = 0.0
    sex = 0.0
    w = x = cosp = 0.0
    if wvno < xka:
        sinp = math.sin(p)
        w = sinp / ra
        x = -ra * sinp
        cosp = math.cos(p)
    elif wvno == xka:
        cosp = 1.0
        w = dpth
        x = 0.0
    elif wvno > xka:
        pex = p
        fac = 0.0
        if p < 16:
            fac = math.exp(-2.0 * p)
        cosp = (1.0 + fac) * 0.5
        sinp = (1.0 - fac) * 0.5
        w = sinp / ra
        x = ra * sinp
    y = z = cosq = 0.0
    if wvno < xkb:
        sinq = math.sin(q)
        y = sinq / rb
        z = -rb * sinq
        cosq = math.cos(q)
    elif wvno == xkb:
        cosq = 1.0
        y = dpth
        z = 0.0
    elif wvno > xkb:
        sex = q
        fac = 0.0
        if q < 16:
            fac = math.exp(-2.0 * q)
        cosq = (1.0 + fac) * 0.5
        sinq = (1.0 - fac) * 0.5
        y = sinq / rb
        z = rb * sinq
    exa = pex + sex
    a0 = 0.0
    if exa < 60.0:
        a0 = math.exp(-exa)
    cpcq = cosp * cosq
    cpy = cosp * y
    cpz = cosp * z
    cqw = cosq * w
    cqx = cosq * x
    xy = x * y
    xz = x * z
    wy = w * y
    wz = w * z
    # (var then scales cosq, y and z by exp(sex - pex); dltar4 does not use them afterwards)
    return w, cosp, exa, a0, cpcq, cpy, cpz, cqw, cqx, xy, xz, wy, wz


@njit(cache=True, nogil=True)
def _dnka(ca, wvno2, gam, gammk, rho, a0, cpcq, cpy, cpz, cqw, cqx, xy, xz, wy, wz):
    one, two = 1.0, 2.0
    gamm1 = gam - one
    twgm1 = gam + gamm1
    gmgmk = gam * gammk
    gmgm1 = gam * gamm1
    gm1sq = gamm1 * gamm1
    rho2 = rho * rho
    a0pq = a0 - cpcq
    ca[0, 0] = cpcq - two * gmgm1 * a0pq - gmgmk * xz - wvno2 * gm1sq * wy
    ca[0, 1] = (wvno2 * cpy - cqx) / rho
    ca[0, 2] = -(twgm1 * a0pq + gammk * xz + wvno2 * gamm1 * wy) / rho
    ca[0, 3] = (cpz - wvno2 * cqw) / rho
    ca[0, 4] = -(two * wvno2 * a0pq + xz + wvno2 * wvno2 * wy) / rho2
    ca[1, 0] = (gmgmk * cpz - gm1sq * cqw) * rho
    ca[1, 1] = cpcq
    ca[1, 2] = gammk * cpz - gamm1 * cqw
    ca[1, 3] = -wz
    ca[1, 4] = ca[0, 3]
    ca[3, 0] = (gm1sq * cpy - gmgmk * cqx) * rho
    ca[3, 1] = -xy
    ca[3, 2] = gamm1 * cpy - gammk * cqx
    ca[3, 3] = ca[1, 1]
    ca[3, 4] = ca[0, 1]
    ca[4, 0] = -(two * gmgmk * gm1sq * a0pq + gmgmk * gmgmk * xz + gm1sq * gm1sq * wy) * rho2
    ca[4, 1] = ca[3, 0]
    ca[4, 2] = -(gammk * gamm1 * twgm1 * a0pq + gam * gammk * gammk * xz + gamm1 * gm1sq * wy) * rho
    ca[4, 3] = ca[1, 0]
    ca[4, 4] = ca[0, 0]
    t = -two * wvno2
    ca[2, 0] = t * ca[4, 2]
    ca[2, 1] = t * ca[3, 2]
    ca[2, 2] = a0 + two * (cpcq - ca[0, 0])
    ca[2, 3] = t * ca[1, 2]
    ca[2, 4] = t * ca[0, 2]


@njit(cache=True, nogil=True)
def _dltar4(wvno, omga, d, a, b, rho, mmax, llw, e, ee, ca):
    """Rayleigh-wave dispersion function (Dunkin's compound matrices), double precision."""
    omega = omga
    if omega < 1.0e-4:
        omega = 1.0e-4
    wvno2 = wvno * wvno
    xka = omega / F8(a[mmax - 1])
    xkb = omega / F8(b[mmax - 1])
    wvnop = wvno + xka
    wvnom = abs(wvno - xka)
    ra = math.sqrt(wvnop * wvnom)
    wvnop = wvno + xkb
    wvnom = abs(wvno - xkb)
    rb = math.sqrt(wvnop * wvnom)
    t = F8(b[mmax - 1]) / omega
    gammk = 2.0 * t * t
    gam = gammk * wvno2
    gamm1 = gam - 1.0
    rho1 = F8(rho[mmax - 1])
    e[0] = rho1 * rho1 * (gamm1 * gamm1 - gam * gammk * ra * rb)
    e[1] = -rho1 * ra
    e[2] = rho1 * (gamm1 - gammk * ra * rb)
    e[3] = rho1 * rb
    e[4] = wvno2 - ra * rb
    for m in range(mmax - 2, llw - 2, -1):
        xka = omega / F8(a[m])
        xkb = omega / F8(b[m])
        t = F8(b[m]) / omega
        gammk = 2.0 * t * t
        gam = gammk * wvno2
        wvnop = wvno + xka
        wvnom = abs(wvno - xka)
        ra = math.sqrt(wvnop * wvnom)
        wvnop = wvno + xkb
        wvnom = abs(wvno - xkb)
        rb = math.sqrt(wvnop * wvnom)
        dpth = F8(d[m])
        rho1 = F8(rho[m])
        p = ra * dpth
        q = rb * dpth
        w, cosp, exa, a0, cpcq, cpy, cpz, cqw, cqx, xy, xz, wy, wz = _var(p, q, ra, rb, wvno, xka, xkb, dpth)
        _dnka(ca, wvno2, gam, gammk, rho1, a0, cpcq, cpy, cpz, cqw, cqx, xy, xz, wy, wz)
        for i in range(5):
            cr = 0.0
            for j in range(5):
                cr = cr + e[j] * ca[j, i]
            ee[i] = cr
        # normc(ee, exa)
        t1 = 0.0
        for i in range(5):
            if abs(ee[i]) > t1:
                t1 = abs(ee[i])
        if t1 < 1.0e-40:
            t1 = 1.0
        for i in range(5):
            e[i] = ee[i] / t1
    if llw != 1:
        xka = omega / F8(a[0])
        wvnop = wvno + xka
        wvnom = abs(wvno - xka)
        ra = math.sqrt(wvnop * wvnom)
        dpth = F8(d[0])
        rho1 = F8(rho[0])
        p = ra * dpth
        znul = 1.0e-05
        w, cosp = _var(p, znul, ra, znul, wvno, xka, znul, dpth)[:2]
        w0 = -rho1 * w
        return cosp * e[0] + w0 * e[1]
    return e[0]


@njit(cache=True, nogil=True)
def _half(c1, c2, omega, d, a, b, rho, mmax, llw, e, ee, ca):
    c3 = 0.5 * (c1 + c2)
    wvno = omega / c3
    return c3, _dltar4(wvno, omega, d, a, b, rho, mmax, llw, e, ee, ca)


@njit(cache=True, nogil=True)
def _nevill(t, c1, c2, del1, del2, d, a, b, rho, mmax, llw, e, ee, ca, x, y):
    omega = _TWOPI / t
    c3, del3 = _half(c1, c2, omega, d, a, b, rho, mmax, llw, e, ee, ca)
    nev = 1
    nctrl = 1
    m = 1
    while True:
        nctrl = nctrl + 1
        if nctrl >= 100:
            break
        if c3 < min(c1, c2) or c3 > max(c1, c2):
            nev = 0
            c3, del3 = _half(c1, c2, omega, d, a, b, rho, mmax, llw, e, ee, ca)
        s13 = del1 - del3
        s32 = del3 - del2
        if _sign(del3) * _sign(del1) < 0.0:
            c2 = c3
            del2 = del3
        else:
            c1 = c3
            del1 = del3
        if abs(c1 - c2) <= 1.0e-6 * c1:
            break
        if _sign(s13) != _sign(s32):
            nev = 0
        ss1 = abs(del1)
        s1 = _S_D * ss1
        ss2 = abs(del2)
        s2 = _S_D * ss2
        if s1 > ss2 or s2 > ss1 or nev == 0:
            c3, del3 = _half(c1, c2, omega, d, a, b, rho, mmax, llw, e, ee, ca)
            nev = 1
            m = 1
        else:
            if nev == 2:
                x[m] = c3
                y[m] = del3
            else:
                x[0] = c1
                y[0] = del1
                x[1] = c2
                y[1] = del2
                m = 1
            fail = False
            for kk in range(1, m + 1):
                j = m - kk  # 0-based j of the Fortran j = m-kk+1
                denom = y[m] - y[j]
                if abs(denom) < 1.0e-10 * abs(y[m]):
                    fail = True
                    break
                x[j] = (-y[j] * x[j + 1] + y[m] * x[j]) / denom
            if fail:
                c3, del3 = _half(c1, c2, omega, d, a, b, rho, mmax, llw, e, ee, ca)
                nev = 1
                m = 1
            else:
                c3 = x[0]
                wvno = omega / c3
                del3 = _dltar4(wvno, omega, d, a, b, rho, mmax, llw, e, ee, ca)
                nev = 2
                m = m + 1
                if m > 10:
                    m = 10
    return c3


@njit(cache=True, nogil=True)
def _getsol(t1, c1, clow, dc, cm, betmx, ifirst, del1st, d, a, b, rho, mmax, llw, e, ee, ca, x, y):
    """-> (iret, c1, del1st)."""
    omega = _TWOPI / t1
    wvno = omega / c1
    del1 = _dltar4(wvno, omega, d, a, b, rho, mmax, llw, e, ee, ca)
    if ifirst == 1:
        del1st = del1
    plmn = _sign(del1st) * _sign(del1)
    if ifirst == 1:
        idir = 1
    elif plmn >= 0.0:
        idir = 1
    else:
        idir = -1
    while True:
        if idir > 0:
            c2 = c1 + dc
        else:
            c2 = c1 - dc
        if c2 <= clow:
            idir = 1
            c1 = clow
            continue
        omega = _TWOPI / t1
        wvno = omega / c2
        del2 = _dltar4(wvno, omega, d, a, b, rho, mmax, llw, e, ee, ca)
        if _sign(del1) != _sign(del2):
            cn = _nevill(t1, c1, c2, del1, del2, d, a, b, rho, mmax, llw, e, ee, ca, x, y)
            c1 = cn
            if c1 > F8(betmx):
                return -1, c1, del1st
            return 1, c1, del1st
        c1 = c2
        del1 = del2
        if c1 < cm:
            return -1, c1, del1st
        if c1 >= (F8(betmx) + dc):
            return -1, c1, del1st


@njit(cache=True, nogil=True)
def _gtsolh(a, b):
    c = _0P95 * b
    for _ in range(5):
        gamma = b / a
        kappa = c / b
        k2 = kappa * kappa
        gk = gamma * kappa
        gk2 = gk * gk
        fac1 = math.sqrt(_ONE - gk2)
        fac2 = math.sqrt(_ONE - k2)
        tk = _TWO - k2
        fr = tk * tk - _FOUR * fac1 * fac2
        frp = -(_FOUR * tk * kappa) + _FOUR * fac2 * gamma * gamma * kappa / fac1 + _FOUR * fac1 * kappa / fac2
        frp = frp / b
        c = c - fr / frp
    return c


@njit(cache=True, nogil=True)
def surfdisp96(thkm, vpm, vsm, rhom, nlayer, t, cg):
    """Rayleigh phase velocities cg (double, REAL-rounded) of the layer model at the periods t; 0 after a failure.
    surfdisp96(thkm, vpm, vsm, rhom, nlayer, iflsph=1, iwave=2, mode=1, igr=0, kmax, t, cg)."""
    mmax = nlayer
    d = np.empty(NL, F4)
    a = np.empty(NL, F4)
    b = np.empty(NL, F4)
    rho = np.empty(NL, F4)
    rtp = np.empty(NL, F4)
    btp = np.empty(NL, F4)
    for i in range(mmax):
        b[i] = vsm[i]
        a[i] = vpm[i]
        d[i] = thkm[i]
        rho[i] = rhom[i]
    llw = 1
    if b[0] <= F4(0.0):
        llw = 2
    # sphere(0, 0, ...)
    d[mmax - 1] = _ONE
    dr = 0.0
    r0 = _AR
    for i in range(mmax):
        rtp[i] = rho[i]
    for i in range(mmax):
        dr = dr + F8(d[i])
        r1 = _AR - dr
        z0 = _AR * math.log(_AR / r0)
        z1 = _AR * math.log(_AR / r1)
        d[i] = F4(z1 - z0)
        tmp = (_AR + _AR) / (r0 + r1)
        a[i] = F4(F8(a[i]) * tmp)
        b[i] = F4(F8(b[i]) * tmp)
        btp[i] = F4(tmp)
        r0 = r1
    dhalf = d[mmax - 1]
    d[mmax - 1] = F4(0.0)
    betmx = F4(-1.0e20)
    betmn = F4(1.0e20)
    jmn = 0
    jsol = 1
    for i in range(mmax):
        if b[i] > _0P01 and b[i] < betmn:
            betmn = b[i]
            jmn = i
            jsol = 1
        elif b[i] <= _0P01 and a[i] < betmn:
            betmn = a[i]
            jmn = i
            jsol = 0
        if b[i] > betmx:
            betmx = b[i]
    # ifunc = 2: sphere(2, 1, ...)
    d[mmax - 1] = dhalf
    for i in range(mmax):
        rho[i] = rtp[i] * btp[i] ** _RHOEXP
    d[mmax - 1] = F4(0.0)
    if jsol == 0:
        cc1 = betmn
    else:
        cc1 = _gtsolh(a[jmn], b[jmn])
    cc1 = _0P95 * cc1
    cc1 = _0P90 * cc1
    cc = F8(cc1)
    dc = abs(_DC)
    cm = cc
    kmax = t.size
    c = np.zeros(kmax, F8)
    e = np.empty(5, F8)
    ee = np.empty(5, F8)
    ca = np.empty((5, 5), F8)
    x = np.empty(20, F8)
    y = np.empty(20, F8)
    del1st = 0.0
    c1 = cc
    for k in range(kmax):
        t1 = t[k]
        if k == 0:
            c1 = cc
            clow = cc
            ifirst = 1
        else:
            ifirst = 0
            c1 = c[k - 1] - _ONEA * dc
            clow = cm
        iret, c1, del1st = _getsol(t1, c1, clow, dc, cm, betmx, ifirst, del1st, d, a, b, rho, mmax, llw, e, ee, ca,
                                   x, y)  # fmt: skip
        if iret == -1:
            for i in range(k, kmax):
                cg[i] = 0.0
            return
        c[k] = c1
        cg[k] = F8(F4(c[k]))


@njit(cache=True, nogil=True)
def column_phase(dep, vs, minthk, t):
    """caldespersion for one column: phase velocities at t of the REAL node model (dep, vs)."""
    vp, rho = brocher(vs)
    rmax, rthk, rvp, rvs, rrho = refine(minthk, dep, vp, vs, rho)
    cg = np.zeros(t.size, F8)
    surfdisp96(rthk, rvp, rvs, rrho, rmax, t, cg)
    return cg


@njit(cache=True, nogil=True)
def column_kernel(dep, vsz, minthk, t, nodes, sen_vs, sen_vp, sen_rho):
    """depthkernel for one column: phase velocities (returned) and the partials of the first `nodes` nodes, written
    into sen_*[period, node] (double).  nodes = nz gives DSurfTomo's full loop; CalSurfG reads nz - 1."""
    vpz, rhoz = brocher(vsz)
    cg = np.zeros(t.size, F8)
    rmax, rthk, rvp, rvs, rrho = refine(minthk, dep, vpz, vsz, rhoz)
    surfdisp96(rthk, rvp, rvs, rrho, rmax, t, cg)
    vsm = vsz.copy()
    vpm = vpz.copy()
    rhom = rhoz.copy()
    cg1 = np.zeros(t.size, F8)
    cg2 = np.zeros(t.size, F8)
    hv = _HALF * _0P01  # 0.5*dlnVs, REAL
    for i in range(nodes):
        vsm[i] = vsz[i] - hv * vsz[i]
        rmax, rthk, rvp, rvs, rrho = refine(minthk, dep, vpm, vsm, rhom)
        surfdisp96(rthk, rvp, rvs, rrho, rmax, t, cg1)
        vsm[i] = vsz[i] + hv * vsz[i]
        rmax, rthk, rvp, rvs, rrho = refine(minthk, dep, vpm, vsm, rhom)
        surfdisp96(rthk, rvp, rvs, rrho, rmax, t, cg2)
        vsm[i] = vsz[i]
        den = F8(_0P01 * vsz[i])
        for nn in range(t.size):
            sen_vs[nn, i] = (cg2[nn] - cg1[nn]) / den
        vpm[i] = vpz[i] - hv * vpz[i]
        rmax, rthk, rvp, rvs, rrho = refine(minthk, dep, vpm, vsm, rhom)
        surfdisp96(rthk, rvp, rvs, rrho, rmax, t, cg1)
        vpm[i] = vpz[i] + hv * vpz[i]
        rmax, rthk, rvp, rvs, rrho = refine(minthk, dep, vpm, vsm, rhom)
        surfdisp96(rthk, rvp, rvs, rrho, rmax, t, cg2)
        vpm[i] = vpz[i]
        den = F8(_0P01 * vpz[i])
        for nn in range(t.size):
            sen_vp[nn, i] = (cg2[nn] - cg1[nn]) / den
        rhom[i] = rhoz[i] - hv * rhoz[i]
        rmax, rthk, rvp, rvs, rrho = refine(minthk, dep, vpm, vsm, rhom)
        surfdisp96(rthk, rvp, rvs, rrho, rmax, t, cg1)
        rhom[i] = rhoz[i] + hv * rhoz[i]
        rmax, rthk, rvp, rvs, rrho = refine(minthk, dep, vpm, vsm, rhom)
        surfdisp96(rthk, rvp, rvs, rrho, rmax, t, cg2)
        rhom[i] = rhoz[i]
        den = F8(_0P01 * rhoz[i])
        for nn in range(t.size):
            sen_rho[nn, i] = (cg2[nn] - cg1[nn]) / den
    return cg
