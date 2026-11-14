"""Single-precision LSMR of DSurfTomo: storage-order products and convergence to the least-squares solution."""

import numpy as np
import pytest

pytest.importorskip("numba")
from lant.inversion import dsurf_lsmr as L  # noqa: E402


def random_matrix(rng, m=60, n=25, density=0.3):
    mask = rng.random((m, n)) < density
    rows, cols = np.nonzero(mask)
    values = rng.standard_normal(len(rows)).astype(np.float32)
    return L.Matrix(rows, cols, values, m, n)


def test_parallel_products_equal_the_literal_aprod_loop():
    rng = np.random.default_rng(1)
    a = random_matrix(rng)
    x = rng.standard_normal(a.n).astype(np.float32)
    y = rng.standard_normal(a.m).astype(np.float32)
    y1, y2 = y.copy(), y.copy()
    L._times(a.row_start, a.cols, a.values, x, y1)
    L.aprod_serial(1, a.rows, a.cols, a.values, x, y2)
    np.testing.assert_array_equal(y1, y2)
    x1, x2 = x.copy(), x.copy()
    L._times_transpose(a.col_start, a.by_col, a.rows, a.values, y, x1)
    L.aprod_serial(2, a.rows, a.cols, a.values, x2, y)
    np.testing.assert_array_equal(x1, x2)


def test_dnrm2_is_the_scaled_sum_of_squares():
    v = np.array([3.0, 0.0, -4.0], np.float32)
    assert L.dnrm2(v) == np.float32(5.0)
    assert L.dnrm2(np.array([-2.5], np.float32)) == np.float32(2.5)
    assert L.dnrm2(np.zeros(0, np.float32)) == 0


@pytest.mark.parametrize("damp", [0.0, 0.5])
def test_lsmr_solves_damped_least_squares(damp):
    rng = np.random.default_rng(2)
    a = random_matrix(rng)
    b = rng.standard_normal(a.m).astype(np.float32)
    r = L.lsmr(a, b, damp, 1e-6, 1e-6, 100, 400, 10)
    dense = np.zeros((a.m, a.n))
    dense[a.rows, a.cols] = a.values
    full = np.vstack([dense, damp * np.eye(a.n)])
    expected = np.linalg.lstsq(full, np.r_[b, np.zeros(a.n)], rcond=None)[0]
    assert r.x.dtype == np.float32 and 1 <= r.itn <= 400 and r.istop in (1, 2, 3)
    np.testing.assert_allclose(r.x, expected, atol=1e-3 * np.abs(expected).max())


def test_zero_right_hand_side_returns_zero():
    a = random_matrix(np.random.default_rng(3))
    r = L.lsmr(a, np.zeros(a.m, np.float32), 1.0, 1e-6, 1e-6, 100, 400, 10)
    assert r.itn == 0 and not r.x.any()


def test_matrix_rejects_rows_out_of_storage_order():
    with pytest.raises(ValueError, match="grouped by row"):
        L.Matrix([1, 0], [0, 0], [1.0, 1.0], 2, 1)
