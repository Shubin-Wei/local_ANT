"""Compiled inner loops of the cc step: the linear detrend of a segment (detrend_linear) and the pair products of
the segment spectra of one time bin.

For every pair p = (iu[p], ju[p]) and frequency f:  acc[p, f] = sum over the segments k, in their order, of
conj(F[k, i, f]) * F[k, j, f], each product formed in single precision (as the spectra are stored) and added in
double precision.  This is the sum that the cc step formed segment by segment with NumPy before; done here per time
bin, a block of frequencies at a time, so that the spectra of the block stay in the cache, without the gathered copies
of the spectra for every pair, and spread over threads (the kernel releases the interpreter lock).
"""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
from numba import njit

FREQ_BLOCK = 16  # frequencies per block: the spectra of a block, segments x stations x 16 x 8 bytes, stay in cache


@njit(cache=True, nogil=True)
def _pair_sums(spec, good, nk, iu, ju, p0, p1, acc):
    nb = spec.shape[2]
    for f0 in range(0, nb, FREQ_BLOCK):
        f1 = min(nb, f0 + FREQ_BLOCK)
        for p in range(p0, p1):
            i = iu[p]
            j = ju[p]
            for f in range(f0, f1):
                acc[p, f] = 0.0
            for k in range(nk):
                if good[k, i] and good[k, j]:  # a product with a zeroed spectrum adds nothing
                    for f in range(f0, f1):
                        acc[p, f] += np.complex128(np.conj(spec[k, i, f]) * spec[k, j, f])


def pair_sums(spec, good, nk, iu, ju, acc, threads=1):
    """acc[p] = sum_k<nk conj(spec[k, iu[p]]) * spec[k, ju[p]] where good[k] holds both stations, else 0.
    spec complex64 [K, nsta, nfreq], good bool [K, nsta], acc complex128 [npair, nfreq] (written)."""
    npair = iu.size
    iu = np.ascontiguousarray(iu, np.int64)
    ju = np.ascontiguousarray(ju, np.int64)
    if threads <= 1 or npair < 2:
        _pair_sums(spec, good, nk, iu, ju, 0, npair, acc)
        return acc
    edges = np.linspace(0, npair, 4 * threads + 1).astype(np.int64)  # more pieces than threads: even load
    with ThreadPoolExecutor(threads) as ex:
        for f in [ex.submit(_pair_sums, spec, good, nk, iu, ju, a, b, acc) for a, b in zip(edges[:-1], edges[1:])]:
            f.result()
    return acc


@njit(cache=True, nogil=True)
def detrend_linear(X):
    """scipy.signal.detrend(X, axis=1, type="linear") of a 2-D float64 array, bit for bit: the same model (columns
    (1..N)/N and 1), the same LAPACK least squares (xGELSD through NumPy's lstsq in Numba) and the same subtraction,
    without the Python overhead and without holding the GIL (the segments of an hour run in threads)."""
    ns, N = X.shape
    A = np.ones((N, 2))
    for i in range(N):
        A[i, 0] = (i + 1) / N
    B = np.ascontiguousarray(X.T)
    coef = np.linalg.lstsq(A, B)[0]
    return (B - A @ coef).T
