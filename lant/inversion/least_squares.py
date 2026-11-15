"""The least-squares solve of each 3-D iteration (lant 2.0.1, in place of DSurfTomo's single-precision LSMR call).

LSMR (Fong & Saunders 2011, SIAM J. Sci. Comput. 33, 2950-2971; scipy.sparse.linalg.lsmr) in double precision on
the column-scaled matrix A D, D = diag(1 / ||a_j||), x = D y:

- double precision: DSurfTomo runs LSMR in REAL with atol = btol = 1e-6.  Its stopping tests compare quantities
  down to about 1e-6 relative while REAL resolves about 1e-7, so they are seldom met: on lant's areal synthetic
  example four of six iterations ended at the iteration limit (400).  A Krylov solve stopped by its iteration
  count is itself a regularization (Hansen 1998, Rank-Deficient and Discrete Ill-Posed Problems, ch. 6), whose
  strength depends on the size of the problem, on top of the explicit one (model_regularization), so that the
  smoothing weight no longer meant the same on every data set.
- column scaling: the travel-time columns (s per km/s) and the smoothing rows differ by orders of magnitude between
  shallow and deep nodes; scaling every column to unit length is the simple preconditioning that LSQR's authors
  advise (Paige & Saunders 1982, ACM TOMS 8, 43-71) and does not change the least-squares solution, only how fast
  the iterations reach it.
- tolerances of the authors (atol = btol = 1e-6, conlim = 1e8) and their iteration limit, min(m, n); no damping.
"""

from dataclasses import dataclass

import numpy as np
from numba import njit, prange
from scipy.sparse.linalg import LinearOperator, lsmr

ATOL = BTOL = 1e-6
CONLIM = 1e8


@njit(cache=True, parallel=True)
def _times(row_start, cols, values, x, y):
    """y = A x."""
    for r in prange(row_start.size - 1):
        acc = 0.0
        for k in range(row_start[r], row_start[r + 1]):
            acc += values[k] * x[cols[k]]
        y[r] = acc


@njit(cache=True, parallel=True)
def _times_transpose(col_start, by_col, rows, values, y, x):
    """x = A' y."""
    for c in prange(col_start.size - 1):
        acc = 0.0
        for i in range(col_start[c], col_start[c + 1]):
            k = by_col[i]
            acc += values[k] * y[rows[k]]
        x[c] = acc


@dataclass
class Solution:
    x: np.ndarray  # float64, n
    istop: int  # scipy/Fong & Saunders codes: 1, 2 converged, 3 condition limit, 7 iteration limit
    itn: int


def solve(rows, cols, values, m, n, b, atol=ATOL, btol=BTOL, conlim=CONLIM, maxiter=None):
    """min ||A x - b|| for the m x n matrix of nonzeros (rows, cols, values) (any order)."""
    rows = np.asarray(rows, np.int64)
    cols = np.asarray(cols, np.int64)
    norm = np.sqrt(np.bincount(cols, weights=np.asarray(values, np.float64) ** 2, minlength=n))
    d = np.where(norm > 0, 1.0 / np.where(norm > 0, norm, 1.0), 1.0)
    vals = np.asarray(values, np.float64) * d[cols]
    if rows.size and np.any(rows[1:] < rows[:-1]):
        order = np.argsort(rows, kind="stable")
        r_rows, r_cols, r_vals = rows[order], cols[order], vals[order]
    else:  # already grouped by row (dsurf_main: data rows, then the smoothing rows)
        r_rows, r_cols, r_vals = rows, cols, vals
    row_start = np.searchsorted(r_rows, np.arange(m + 1)).astype(np.int64)
    by_col = np.argsort(r_cols, kind="stable").astype(np.int64)
    col_start = np.searchsorted(r_cols[by_col], np.arange(n + 1)).astype(np.int64)

    def matvec(x):
        y = np.empty(m)
        _times(row_start, r_cols, r_vals, np.ascontiguousarray(x, np.float64).ravel(), y)
        return y

    def rmatvec(y):
        x = np.empty(n)
        _times_transpose(col_start, by_col, r_rows, r_vals, np.ascontiguousarray(y, np.float64).ravel(), x)
        return x

    op = LinearOperator((m, n), matvec=matvec, rmatvec=rmatvec, dtype=np.float64)
    y, istop, itn = lsmr(op, np.asarray(b, np.float64), atol=atol, btol=btol, conlim=conlim, maxiter=maxiter)[:3]
    return Solution(y * d, int(istop), int(itn))
