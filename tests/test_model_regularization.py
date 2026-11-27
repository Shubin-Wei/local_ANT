"""lant's regularization of the 3-D inversion: the graph Laplacian of the parameter box, smoothing and damping of
the model itself (not of the update)."""

import numpy as np
import pytest

pytest.importorskip("numba")
from lant.inversion import model_regularization  # noqa: E402


def test_laplacian_is_the_graph_laplacian_of_the_box():
    rows, cols, vals = model_regularization.laplacian(3, 4, 2)
    n = 3 * 4 * 2
    L = np.zeros((n, n))
    np.add.at(L, (rows, cols), vals)
    assert np.allclose(L, L.T) and np.allclose(L.sum(1), 0)  # symmetric, constants in the null space
    assert L[0, 0] == 3  # a corner node has three neighbours
    inner = 1 * 3 * 4 + 1 * 3 + 1
    assert L[inner, inner] == 5 and np.count_nonzero(L[inner]) == 6
    assert np.linalg.matrix_rank(L) == n - 1


def test_rows_hold_the_model_not_the_update():
    m0 = np.ones((2, 2, 2))
    reg = model_regularization.Regularization(2, 2, 2, m0, weight=2.0)
    m = m0.copy()
    m[0, 0, 0] += 0.1
    rows, cols, vals, rhs = reg.rows(m, first_row=10)
    assert rows.min() == 10 and rows.max() == 10 + 8 - 1 and reg.row_count() == 8
    s = 1 / np.sqrt(8)  # the mean over the 8 parameters
    assert rhs[0] == pytest.approx(-2.0 * s * 3 * 0.1)
    assert np.allclose(reg.rows(m0, 0)[3], 0)  # the starting model itself costs nothing


def test_damping_rows_pull_towards_the_starting_model():
    m0 = np.ones((2, 2, 2))
    reg = model_regularization.Regularization(2, 2, 2, m0, weight=1.0, damping=0.5)
    m = m0 + 0.2
    rows, cols, vals, rhs = reg.rows(m, first_row=0)
    assert reg.row_count() == 16 and rows.max() == 15
    damp = rows >= 8
    s = 1 / np.sqrt(8)
    assert np.allclose(vals[damp], 0.5 * s) and np.allclose(rhs[8:], -0.5 * s * 0.2, atol=1e-7)
    # the least-squares update of these rows alone brings the model back to m0
    A = np.zeros((16, 8))
    np.add.at(A, (rows, cols), vals)
    dv = np.linalg.lstsq(A[8:], rhs[8:], rcond=None)[0]
    assert np.allclose(m.ravel() + dv, 1.0, atol=1e-6)
