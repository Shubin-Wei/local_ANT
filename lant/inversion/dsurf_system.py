"""The linear system of one DSurfTomo iteration, as main.f90 builds it (REAL arithmetic).

DSurfTomo a651faf (MIT), main.f90 lines 361-467 and getpercentile.f90:

- residual ``cbst(i) = obst(i) - dsyn(i)`` index by index: obst in file order,
  dsyn in CalSurfG order (they agree when surfdata.dat is grouped by period in
  increasing order, as lant writes it);
- getpercentile: heap sort of the residuals, q25 = sorted(int(0.25 N)),
  q75 = sorted(int(0.75 N)) (1-based, REAL product truncated);
- rejection: a residual below threshold x q25 or above threshold x q75 gets
  weight 0 and residual 0; the nonzeros of its row are multiplied by the weight
  (they stay in the matrix as zeros);
- DWS: the sum of |G| over each parameter column, its maximum and mean (log only);
- smoothing rows, one per parameter in (depth, longitude, latitude) order: on
  the faces of the parameter box only the diagonal 2 x weight, inside the
  7-point Laplacian 6 x weight and -1 x weight for the latitude -/+, longitude
  -/+ and depth -/+ neighbours, in that storage order; right-hand side 0;
- residual statistics after LSMR: REAL sums of the weighted residuals and their
  squares, and the dnrm2-based rms of the log line.
"""

import math

import numpy as np
from numba import njit

from .dsurf_lsmr import dnrm2

F = np.float32


def percentiles(residuals):
    """getpercentile.f90: q25 and q75 of the REAL residuals."""
    n = len(residuals)
    sorted_ = np.sort(np.asarray(residuals, F), kind="stable")
    i25, i75 = int(F(0.25) * F(n)), int(F(0.75) * F(n))
    if i25 < 1:
        # RA(0) is outside the array; the native build (-fbounds-check) stops there.
        raise ValueError("getpercentile needs at least 4 travel times, got %d" % n)
    return sorted_[i25 - 1], sorted_[i75 - 1]


def reject(residuals, q25, q75, threshold):
    """main.f90 rejection: (weights, residuals with rejected entries set to 0)."""
    r = np.asarray(residuals, F).copy()
    out = (r < F(q25) * F(threshold)) | (r > F(q75) * F(threshold))
    weights = np.ones(len(r), F)
    weights[out] = F(0.0)
    r[out] = F(0.0)
    return weights, r


MAD_SCALE = 1.4826  # MAD x 1.4826 estimates the standard deviation of normal data (Rousseeuw & Croux 1993)
# Rows a period needs for its own median and MAD; fewer: the nearest periods are pooled with it until there are
# this many.  The MAD scale of n normal values has a standard error of about 1.16 sigma / sqrt(n) (its asymptotic
# efficiency is 37 %, Rousseeuw & Croux 1993), about 25 % at 20 rows.
MAD_MIN_ROWS = 20


def robust_centre_scale(residuals, groups=None, min_rows=MAD_MIN_ROWS):
    """Median and 1.4826 MAD of the residuals, per group (period) when groups is given -> two arrays, one value per
    residual.  A group with fewer than min_rows residuals takes the nearest groups (in group order) as well."""
    r = np.asarray(residuals, np.float64)
    if groups is None:
        med = np.median(r) if r.size else 0.0
        return np.full(r.size, med), np.full(r.size, MAD_SCALE * np.median(np.abs(r - med)) if r.size else 0.0)
    g = np.asarray(groups)
    keys, inverse, counts = np.unique(g, return_inverse=True, return_counts=True)
    centre, scale = np.zeros(r.size), np.zeros(r.size)
    for i in range(len(keys)):
        lo = hi = i
        while counts[lo : hi + 1].sum() < min_rows and (lo > 0 or hi < len(keys) - 1):
            # widen towards the side whose next group is nearer in group order; on a tie both
            below = keys[i] - keys[lo - 1] if lo > 0 else np.inf
            above = keys[hi + 1] - keys[i] if hi < len(keys) - 1 else np.inf
            if below <= above:
                lo -= 1
            if above <= below:
                hi += 1
        pool = r[(inverse >= lo) & (inverse <= hi)]
        med = np.median(pool)
        centre[inverse == i] = med
        scale[inverse == i] = MAD_SCALE * np.median(np.abs(pool - med))
    return centre, scale


def reject_mad(residuals, k, groups=None, min_rows=MAD_MIN_ROWS):
    """lant's rejection (replaces main.f90's): a travel-time residual further than k x 1.4826 x MAD from the median
    gets weight 0 and residual 0 (median absolute deviation; Rousseeuw & Croux 1993, JASA 88, 1273-1283; k = 3,
    'very conservative', Leys et al. 2013, J. Exp. Soc. Psychol. 49, 764-766).  Centred on the median and
    symmetric: DSurfTomo's rule (outside threshold x [q25, q75] of the signed residuals) is neither, and on a
    synthetic checkerboard with 1 % noise it dropped 12 % of the data, mostly the pattern's own signal.

    groups: the period of each residual.  The rule assumes the residuals share one distribution, but a phase error
    becomes a travel-time error proportional to the period (t = -(phi - phi0) / omega): on lant's areal synthetic
    example the robust spread rose from 7 ms at 0.21 s to 19 ms at 0.84 s, and one MAD of all periods took 77 of
    its 96 rejections from the longest third of the periods.  The median and the MAD are therefore taken per period
    (robust_centre_scale).  -> (weights, residuals with rejected entries set to 0)."""
    r = np.asarray(residuals, F).copy()
    centre, scale = robust_centre_scale(r, groups, min_rows)
    out = (np.abs(r.astype(np.float64) - centre) > k * scale) & (scale > 0)
    weights = np.ones(len(r), F)
    weights[out] = F(0.0)
    r[out] = F(0.0)
    return weights, r


@njit(cache=True)
def _dws(cols, values, n):
    norm = np.zeros(n, np.float32)
    for i in range(values.size):
        norm[cols[i]] = norm[cols[i]] + abs(values[i])
    total = np.float32(0.0)
    largest = np.float32(0.0)
    for i in range(n):
        total = total + norm[i]
        if norm[i] > largest:
            largest = norm[i]
    return largest, total / np.float32(n)


def dws(cols, values, n):
    """Maximum and average DWS (sum of |G| per column), REAL as main.f90."""
    largest, average = _dws(np.ascontiguousarray(cols, np.int64), np.ascontiguousarray(values, F), n)
    return F(largest), F(average)


@njit(cache=True)
def _column_sums(cols, values, n):
    norm = np.zeros(n, np.float32)
    for i in range(values.size):
        norm[cols[i]] = norm[cols[i]] + abs(values[i])
    return norm


def dws_nodes(cols, values, n):
    """DWS of every column (sum of |G|, REAL, as _dws), for lant: > 0 for nonzero post-rejection Jacobian column support."""
    return _column_sums(np.ascontiguousarray(cols, np.int64), np.ascontiguousarray(values, F), n)


def smoothing(nx, ny, nz, weight, first_row):
    """Smoothing rows (rows, cols, values) in main.f90 storage order; rows start at first_row (0-based)."""
    nvx, nvz = nx - 2, ny - 2
    w = F(weight)
    diag_face, diag_inside, off = F(2.0) * w, F(6.0) * w, F(-1.0) * w
    k, j, i = np.meshgrid(np.arange(nz - 1), np.arange(nvz), np.arange(nvx), indexing="ij")
    k, j, i = k.ravel(), j.ravel(), i.ravel()
    node = k * nvz * nvx + j * nvx + i
    face = (i == 0) | (i == nvx - 1) | (j == 0) | (j == nvz - 1) | (k == 0) | (k == nz - 2)
    counts = np.where(face, 1, 7)
    rows = np.repeat(first_row + np.arange(node.size), counts)
    cols = np.empty(counts.sum(), np.int64)
    values = np.empty(counts.sum(), F)
    start = np.r_[0, np.cumsum(counts)[:-1]]
    cols[start] = node
    values[start] = np.where(face, diag_face, diag_inside)
    inside = start[~face]
    n_in = node[~face]
    for offset, step in enumerate((-1, 1, -nvx, nvx, -nvz * nvx, nvz * nvx), 1):
        cols[inside + offset] = n_in + step
        values[inside + offset] = off
    return rows, cols, values


@njit(cache=True)
def _sums(r):
    s = np.float32(0.0)
    s2 = np.float32(0.0)
    for i in range(r.size):
        s = s + r[i]
        s2 = s2 + r[i] * r[i]
    return s, s2


def residual_statistics(residuals):
    """mean, std_devs and rms of main.f90's log line (REAL; std is NaN if the variance rounds below 0)."""
    r = np.ascontiguousarray(residuals, F)
    n = F(len(r))
    s, s2 = _sums(r)
    mean = F(s) / n
    variance = F(s2) / n - mean * mean
    std = F(math.sqrt(variance)) if variance >= 0 else F(np.nan)
    rms = F(dnrm2(r)) / F(math.sqrt(n))
    return mean, std, rms
