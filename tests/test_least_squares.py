"""The double-precision, column-scaled LSMR solve of each 3-D iteration."""

import numpy as np
import pytest

pytest.importorskip("numba")
from lant.inversion import least_squares as LS  # noqa: E402


def _problem(seed, m=300, n=40, spread=1e4):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(m, n)) * (rng.random((m, n)) < 0.2)
    a[np.arange(n), np.arange(n)] += 1.0  # full column rank
    a *= np.logspace(0, np.log10(spread), n)[None, :]  # columns of very different size, as shallow and deep nodes
    b = rng.normal(size=m)
    r, c = np.nonzero(a)
    return a, r, c, a[r, c], b


def test_solution_is_the_least_squares_solution_despite_badly_scaled_columns():
    a, r, c, v, b = _problem(0)
    s = LS.solve(r, c, v, a.shape[0], a.shape[1], b)
    x = np.linalg.lstsq(a, b, rcond=None)[0]
    assert s.istop in (1, 2)
    assert np.abs(s.x - x).max() <= 1e-4 * np.abs(x).max()  # atol = 1e-6 times the condition of the scaled matrix


def test_nonzeros_in_any_order_give_the_same_solution():
    a, r, c, v, b = _problem(1)
    p = np.random.default_rng(2).permutation(len(v))
    s1 = LS.solve(r, c, v, a.shape[0], a.shape[1], b)
    s2 = LS.solve(r[p], c[p], v[p], a.shape[0], a.shape[1], b)
    assert np.allclose(s1.x, s2.x, rtol=1e-6)


def test_empty_columns_stay_zero():
    a, r, c, v, b = _problem(3)
    keep = c != 5
    s = LS.solve(r[keep], c[keep], v[keep], a.shape[0], a.shape[1], b)
    assert s.x[5] == 0
