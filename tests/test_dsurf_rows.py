"""Row assembly of the DSurfTomo-compatible sensitivity matrix."""

import numpy as np
import pytest

pytest.importorskip("numba")
from lant.inversion import dsurf_rows as rows  # noqa: E402


def materials(vs):
    """Brocher (2005) Vp and density and their derivatives with respect to Vs (double), the reference."""
    vs = np.asarray(vs, float)
    vp = 0.9409 + 2.0947 * vs - 0.8206 * vs**2 + 0.2683 * vs**3 - 0.0251 * vs**4
    rho = 1.6612 * vp - 0.4721 * vp**2 + 0.0671 * vp**3 - 0.0043 * vp**4 + 0.000106 * vp**5
    dvp = 2.0947 - 2 * 0.8206 * vs + 3 * 0.2683 * vs**2 - 4 * 0.0251 * vs**3
    drho = dvp * (1.6612 - 2 * 0.4721 * vp + 3 * 0.0671 * vp**2 - 4 * 0.0043 * vp**3 + 5 * 0.000106 * vp**4)
    return np.column_stack((vp, vs, rho)), dvp, drho


def test_chain_matches_brocher_derivatives():
    vs = np.array([0.8, 1.5, 2.5, 3.5])
    _, dvp, drho = materials(vs)
    for v, a, r in zip(vs, dvp, drho):
        coe_a, coe_rho = rows._chain(np.float64, v, True)
        assert coe_a == pytest.approx(a, rel=1e-12) and coe_rho == pytest.approx(r, rel=1e-12)
    deep_a, _ = rows._chain(np.float64, 3.5, False)
    assert deep_a != pytest.approx(dvp[-1])


def test_cutoffs_order_and_double_products():
    nz, ny, nx = 3, 4, 5
    nvz, nvx = ny - 2, nx - 2
    fdm = np.zeros((2, nvz + 2, nvx + 2), np.float32)
    fdm[0, 1, 2] = -0.5
    fdm[0, 2, 1] = 5e-5  # below ftol: the vertex is skipped
    fdm[0, 0, 1] = -1.0  # outer frame vertex: never a parameter
    fdm[1, 2, 3] = 2e-4
    vels = np.full((nz, ny, nx), 2.0, np.float32)
    shape = (1, ny, nx, nz)
    sen_vs, sen_vp, sen_rho = np.full(shape, 0.3), np.full(shape, 0.1), np.full(shape, 0.05)
    for s in (sen_vs, sen_vp, sen_rho):
        s[..., 1] = 1e-5  # second depth: |row| below ftol after assembly
    g = rows.assemble(fdm, np.array([0, 0]), vels, sen_vs, sen_vp, sen_rho, False)
    coe_a, coe_rho = rows._chain(np.float32, np.float32(2.0), True)
    value = np.float32((0.1 * np.float64(coe_a) + 0.05 * np.float64(coe_rho) + 0.3) * -0.5)
    column = (1 - 1) * nvx + (2 - 1)  # jj = 1, kk = 2
    assert g.dtype == np.float32 and g.shape == (2, (nz - 1) * nvz * nvx)
    assert g[0].indices.tolist() == [column] and g[0, column] == value
    # Row 1: 2e-4 * ~0.38 is below ftol for every depth.
    assert g[1].nnz == 0
