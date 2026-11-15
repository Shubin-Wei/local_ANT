"""lant's regularization (inversion.regularization: lant): the graph Laplacian of the parameter box and the rows that
smooth the model itself (jumping)."""

import numpy as np
import pytest

pytest.importorskip("numba")
from lant.inversion import model_regularization  # noqa: E402
from tests.util import solver_project  # noqa: E402


def test_laplacian_is_the_graph_laplacian_of_the_box():
    rows, cols, vals = model_regularization.laplacian(3, 4, 2)
    n = 3 * 4 * 2
    L = np.zeros((n, n))
    np.add.at(L, (rows, cols), vals)
    assert np.allclose(L, L.T) and np.allclose(L.sum(1), 0)  # symmetric, constants in the null space
    assert L[0, 0] == 3  # a corner node has three neighbours (DSurfTomo: only 2 w on the diagonal)
    inner = 1 * 3 * 4 + 1 * 3 + 1  # k=1 is a face (two depths); interior in i, j
    assert L[inner, inner] == 5 and np.count_nonzero(L[inner]) == 6
    assert np.linalg.matrix_rank(L) == n - 1


def test_rows_hold_the_model_not_the_update():
    m0 = np.ones((2, 2, 2))
    reg = model_regularization.Regularization(2, 2, 2, m0, weight=2.0)
    m = m0.copy()
    m[0, 0, 0] += 0.1  # a perturbed model: the right-hand sides pull it back
    rows, cols, vals, rhs = reg.rows(m, first_row=10)
    assert rows.min() == 10 and rows.max() == 10 + 8 - 1 and reg.row_count() == 8
    s = 1 / np.sqrt(8)  # the mean over the 8 parameters
    assert rhs[0] == pytest.approx(-2.0 * s * 3 * 0.1)
    assert np.allclose(reg.rows(m0, 0)[3], 0)  # the starting model itself costs nothing


def test_inversion_runs_and_records_the_solver_stop(tmp_path):
    from lant.inversion import solver

    d = solver_project(tmp_path / "p", velocity=[2.2, 2.7, 3.1])
    record = solver.run(d, 1)
    assert "Program finishes successfully" in (d / "DSurfTomo.in.log").read_text()
    assert all("lsmr_istop" in it and "rejected" in it for it in record["iterations"])
