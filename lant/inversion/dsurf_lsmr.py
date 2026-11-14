"""LSMR of DSurfTomo in single precision, with its sparse product and BLAS.

DSurfTomo a651faf (MIT) solves with lsmrModule.f90 (Fong & Saunders 2011,
LSMR, version of 2010 with local reorthogonalization) where lsmrDataModule
sets ``dp = selected_real_kind(4)``: every vector and scalar is REAL, and so
are its dnrm2/dscal (lsmrblas.f90, the reference BLAS loops) and the product
aprod.f90 (y += A x and x += A' y, one nonzero after another in storage order).
scipy.sparse.linalg.lsmr is double precision and has no local
reorthogonalization, so it does not reproduce these iterates; this module
follows the Fortran statement by statement in float32:

- aprod: the result of each row (mode 1) or column (mode 2) is accumulated in
  the storage order of its nonzeros, as the Fortran loop over all nonzeros does;
  rows and columns are independent, so they are computed in parallel without
  changing a single sum;
- dnrm2: the scaled sum of squares of the reference BLAS, sequential;
- dot_product of the local reorthogonalization: sequential REAL sum;
- d2norm, the plane rotations, the norm estimates and the stopping tests in the
  order and association of lsmrModule.f90.

``nout`` (the print unit) is not set in main.f90; nothing is printed here.
"""

import math
from dataclasses import dataclass

import numpy as np
from numba import njit, prange

F = np.float32


@dataclass
class Matrix:
    """Nonzeros in DSurfTomo storage order: 0-based rows, columns, REAL values; m x n."""

    rows: np.ndarray
    cols: np.ndarray
    values: np.ndarray
    m: int
    n: int

    def __post_init__(self):
        self.rows = np.ascontiguousarray(self.rows, np.int64)
        self.cols = np.ascontiguousarray(self.cols, np.int64)
        self.values = np.ascontiguousarray(self.values, F)
        if len(self.rows) and (np.any(np.diff(self.rows) < 0)):
            raise ValueError("nonzeros must be grouped by row in increasing row order, as DSurfTomo stores them")
        if len(self.rows) and (self.rows.max() >= self.m or self.cols.max() >= self.n or self.cols.min() < 0):
            raise ValueError("nonzero index outside the matrix")
        self.row_start = np.searchsorted(self.rows, np.arange(self.m + 1)).astype(np.int64)
        # A stable sort keeps, within each column, the storage order of its nonzeros.
        self.by_col = np.argsort(self.cols, kind="stable").astype(np.int64)
        self.col_start = np.searchsorted(self.cols[self.by_col], np.arange(self.n + 1)).astype(np.int64)


@njit(cache=True, parallel=True)
def _times(row_start, cols, values, x, y):
    """y = y + A x (aprod mode 1)."""
    for r in prange(row_start.size - 1):
        acc = y[r]
        for k in range(row_start[r], row_start[r + 1]):
            acc = acc + values[k] * x[cols[k]]
        y[r] = acc


@njit(cache=True, parallel=True)
def _times_transpose(col_start, by_col, rows, values, y, x):
    """x = x + A' y (aprod mode 2)."""
    for c in prange(col_start.size - 1):
        acc = x[c]
        for i in range(col_start[c], col_start[c + 1]):
            k = by_col[i]
            acc = acc + values[k] * y[rows[k]]
        x[c] = acc


@njit(cache=True)
def aprod_serial(mode, rows, cols, values, x, y):
    """aprod.f90 literally: one loop over the nonzeros (reference for the parallel form)."""
    for k in range(values.size):
        if mode == 1:
            y[rows[k]] = y[rows[k]] + values[k] * x[cols[k]]
        else:
            x[cols[k]] = x[cols[k]] + values[k] * y[rows[k]]


@njit(cache=True)
def dnrm2(x):
    """Reference-BLAS dnrm2 in REAL (lsmrblas.f90)."""
    n = x.size
    if n < 1:
        return F(0.0)
    if n == 1:
        return abs(x[0])
    scale = F(0.0)
    ssq = F(1.0)
    for i in range(n):
        if x[i] != F(0.0):
            absxi = abs(x[i])
            if scale < absxi:
                t = scale / absxi
                ssq = F(1.0) + ssq * (t * t)
                scale = absxi
            else:
                t = absxi / scale
                ssq = ssq + t * t
    return scale * F(math.sqrt(ssq))


@njit(cache=True)
def _scale(a, x):
    for i in range(x.size):
        x[i] = a * x[i]


@njit(cache=True)
def _d2norm(a, b):
    scale = abs(a) + abs(b)
    if scale == F(0.0):
        return F(0.0)
    p = a / scale
    q = b / scale
    return scale * F(math.sqrt(p * p + q * q))


@njit(cache=True)
def _dot(a, b):
    s = F(0.0)
    for i in range(a.size):
        s = s + a[i] * b[i]
    return s


@njit(cache=True)
def _lsmr(row_start, cols, col_start, by_col, rows, values, m, n, b, damp, atol, btol, conlim, itnlim, local_size):
    one = F(1.0)
    zero = F(0.0)
    x = np.zeros(n, F)
    u = b.copy()
    v = np.zeros(n, F)
    local_vecs = min(local_size, m, n)
    local_v = np.zeros((max(local_vecs, 1), n), F)
    alpha = zero
    beta = dnrm2(u)
    if beta > zero:
        _scale(one / beta, u)
        _times_transpose(col_start, by_col, rows, values, u, v)
        alpha = dnrm2(v)
    if alpha > zero:
        _scale(one / alpha, v)
    normar = alpha * beta
    stats = np.zeros(6, F)  # normA, condA, normr, normAr, normx, unused
    if normar == zero:
        return x, 0, 0, stats
    local_ortho = False
    local_pointer = 0
    queue_full = False
    if local_vecs > 0:
        local_pointer = 1
        local_ortho = True
        local_v[0, :] = v
    itn = 0
    zetabar = alpha * beta
    alphabar = alpha
    rho = one
    rhobar = one
    cbar = one
    sbar = zero
    h = v.copy()
    hbar = np.zeros(n, F)
    betadd = beta
    betad = zero
    rhodold = one
    tautildeold = zero
    thetatilde = zero
    zeta = zero
    d = zero
    norma2 = alpha * alpha
    maxrbar = zero
    minrbar = F(1e30)
    normb = beta
    istop = 0
    ctol = zero
    if conlim > zero:
        ctol = one / conlim
    normr = beta
    normar = alpha * beta
    norma = zero
    conda = zero
    normx = zero
    while True:
        itn += 1
        _scale(-alpha, u)
        _times(row_start, cols, values, v, u)
        beta = dnrm2(u)
        if beta > zero:
            _scale(one / beta, u)
            if local_ortho:  # localVEnqueue: store the old v
                if local_pointer < local_vecs:
                    local_pointer += 1
                else:
                    local_pointer = 1
                    queue_full = True
                local_v[local_pointer - 1, :] = v
            _scale(-beta, v)
            _times_transpose(col_start, by_col, rows, values, u, v)
            if local_ortho:  # localVOrtho
                limit = local_vecs if queue_full else local_pointer
                for c in range(limit):
                    dd = _dot(v, local_v[c])
                    for i in range(n):
                        v[i] = v[i] - dd * local_v[c, i]
            alpha = dnrm2(v)
            if alpha > zero:
                _scale(one / alpha, v)
        alphahat = _d2norm(alphabar, damp)
        chat = alphabar / alphahat
        shat = damp / alphahat
        rhoold = rho
        rho = _d2norm(alphahat, beta)
        c = alphahat / rho
        s = beta / rho
        thetanew = s * alpha
        alphabar = c * alpha
        rhobarold = rhobar
        zetaold = zeta
        thetabar = sbar * rho
        rhotemp = cbar * rho
        rhobar = _d2norm(cbar * rho, thetanew)
        cbar = cbar * rho / rhobar
        sbar = thetanew / rhobar
        zeta = cbar * zetabar
        zetabar = -sbar * zetabar
        f_hbar = thetabar * rho / (rhoold * rhobarold)
        for i in range(n):
            hbar[i] = h[i] - f_hbar * hbar[i]
        f_x = zeta / (rho * rhobar)
        for i in range(n):
            x[i] = x[i] + f_x * hbar[i]
        f_h = thetanew / rho
        for i in range(n):
            h[i] = v[i] - f_h * h[i]
        betaacute = chat * betadd
        betacheck = -shat * betadd
        betahat = c * betaacute
        betadd = -s * betaacute
        thetatildeold = thetatilde
        rhotildeold = _d2norm(rhodold, thetabar)
        ctildeold = rhodold / rhotildeold
        stildeold = thetabar / rhotildeold
        thetatilde = stildeold * rhobar
        rhodold = ctildeold * rhobar
        betad = -stildeold * betad + ctildeold * betahat
        tautildeold = (zetaold - thetatildeold * tautildeold) / rhotildeold
        taud = (zeta - thetatilde * tautildeold) / rhodold
        d = d + betacheck * betacheck
        t = betad - taud
        normr = F(math.sqrt(d + t * t + betadd * betadd))
        norma2 = norma2 + beta * beta
        norma = F(math.sqrt(norma2))
        norma2 = norma2 + alpha * alpha
        maxrbar = max(maxrbar, rhobarold)
        if itn > 1:
            minrbar = min(minrbar, rhobarold)
        conda = max(maxrbar, rhotemp) / min(minrbar, rhotemp)
        normar = abs(zetabar)
        normx = dnrm2(x)
        test1 = normr / normb
        test2 = normar / (norma * normr)
        test3 = one / conda
        t1 = test1 / (one + norma * normx / normb)
        rtol = btol + atol * norma * normx / normb
        if itn >= itnlim:
            istop = 7
        if one + test3 <= one:
            istop = 6
        if one + test2 <= one:
            istop = 5
        if one + t1 <= one:
            istop = 4
        if test3 <= ctol:
            istop = 3
        if test2 <= atol:
            istop = 2
        if test1 <= rtol:
            istop = 1
        if istop != 0:
            break
    if damp > zero and istop == 2:
        istop = 3
    stats[0], stats[1], stats[2], stats[3], stats[4] = norma, conda, normr, normar, normx
    return x, istop, itn, stats


@dataclass
class Result:
    x: np.ndarray
    istop: int
    itn: int
    norm_a: float
    cond_a: float
    norm_r: float
    norm_ar: float
    norm_x: float


def lsmr(matrix, b, damp, atol, btol, conlim, itnlim, local_size):
    """DSurfTomo's LSMR call (all REAL); b has matrix.m entries."""
    b = np.ascontiguousarray(b, F)
    if b.shape != (matrix.m,):
        raise ValueError("right-hand side must have %d entries" % matrix.m)
    x, istop, itn, s = _lsmr(
        matrix.row_start, matrix.cols, matrix.col_start, matrix.by_col, matrix.rows, matrix.values,
        matrix.m, matrix.n, b, F(damp), F(atol), F(btol), F(conlim), int(itnlim), int(local_size),
    )  # fmt: skip
    return Result(x, int(istop), int(itn), *(float(v) for v in s[:5]))
