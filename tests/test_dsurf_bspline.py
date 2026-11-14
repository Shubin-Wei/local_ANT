"""B-spline propagation and refined grids of the DSurfTomo-compatible forward."""

import numpy as np
import pytest

pytest.importorskip("numba")
from lant.inversion import dsurf_bspline as bs  # noqa: E402


@pytest.mark.parametrize("real", [np.float32, np.float64])
def test_basis_partition_of_unity(real):
    out = np.empty(4, real)
    for u in np.linspace(0, 1, 9, dtype=real):
        bs.basis(real, u, out)
        assert abs(out.sum() - 1) < 10 * np.finfo(real).eps and out.dtype == real
    bs.basis(real, real(0), out)
    assert np.allclose(out, [1 / 6, 4 / 6, 1 / 6, 0])


def test_gridder_reproduces_linear_field_and_shape():
    nvz, nvx = 4, 5
    i, j = np.indices((nvz + 2, nvx + 2))
    velv = (2 + 0.1 * i + 0.05 * j).astype(np.float64)
    veln = bs.gridder(np.float64, velv, 8, 8)
    assert veln.shape == ((nvz - 1) * 8 + 2, (nvx - 1) * 8 + 2)
    # Node (stz, stx) lies at control coordinate (1 + (stz-1)/8, 1 + (stx-1)/8): cubic B-splines reproduce lines.
    stz, stx = np.indices(veln.shape)
    expected = 2 + 0.1 * (1 + (stz - 1) / 8) + 0.05 * (1 + (stx - 1) / 8)
    assert np.allclose(veln[1:, 1:], expected[1:, 1:], atol=1e-12)


def test_refined_grid_matches_coarse_nodes_exactly():
    rng = np.random.default_rng(1)
    velv = (2 + rng.random((6, 7))).astype(np.float32)
    coarse = bs.gridder(np.float32, velv, 8, 8)
    vnl, vnr, vnt, vnb = 3, 19, 1, 17
    fine = bs.bsplrefine(np.float32, velv, 8, 8, 8, vnl, vnr, vnt, vnb)
    assert fine.shape == ((vnb - vnt) * 8 + 2, (vnr - vnl) * 8 + 2)
    # Every 8th refined node coincides with a coarse node; same basis values and summation order.
    np.testing.assert_array_equal(fine[1::8, 1::8], coarse[vnt : vnb + 1, vnl : vnr + 1])
