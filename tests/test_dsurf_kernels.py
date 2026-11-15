"""DSurfTomo depth kernels against the native depthkernel (tests/data/native_depthkernel.npz, made
with the unmodified Fortran objects), and only the work CalSurfG reads."""

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("numba")
from lant.inversion import dsurf_kernels as K  # noqa: E402

REF = np.load(Path(__file__).parent / "data" / "native_depthkernel.npz")


def _native_total():
    """dc/dVs + dc/dVp dVp/dVs + dc/drho drho/dVs of the stored native partials, with CalSurfG's coefficients (its
    shallow branch: the Brocher relation the model is built with)."""
    v = np.transpose(REF["vs"].astype(np.float64), (1, 2, 0))[None]  # (1, ny, nx, nz)
    coe_a = 2.0947 - 0.8206 * 2 * v + 0.2683 * 3 * v**2 - 0.0251 * 4 * v**3
    vp = 0.9409 + 2.0947 * v - 0.8206 * v**2 + 0.2683 * v**3 - 0.0251 * v**4
    coe_rho = coe_a * (1.6612 - 0.4721 * 2 * vp + 0.0671 * 3 * vp**2 - 0.0043 * 4 * vp**3 + 0.000106 * 5 * vp**4)
    return REF["sen_vs"] + REF["sen_vp"] * coe_a + REF["sen_rho"] * coe_rho


def test_native_phase_bit_for_bit_and_its_kernel_combination():
    """Phase velocities as native depthkernel bit for bit; the derivative along the Brocher relation equals the
    native combination of the three partials to the finite-difference error (at most 0.1 % of the largest value of
    each period; measured 0.008-0.03 %)."""
    k = K.depth_kernels(REF["depths"], REF["vs"], REF["periods"], REF["minthk"], workers=2, interior_only=False)
    np.testing.assert_array_equal(k.phase, REF["phase"])
    native = _native_total()
    for p in range(native.shape[0]):
        assert np.abs(k.vs[p] - native[p]).max() <= 1e-3 * np.abs(native[p]).max()
    assert not k.vp.any() and not k.rho.any()  # the whole derivative is in k.vs


def test_interior_only_skips_unused_partials():
    args = (REF["depths"], REF["vs"], REF["periods"], REF["minthk"])
    full = K.depth_kernels(*args, interior_only=False)
    used = K.depth_kernels(*args)
    np.testing.assert_array_equal(used.phase, full.phase)
    np.testing.assert_array_equal(used.vs[:, 1:-1, 1:-1, :-1], full.vs[:, 1:-1, 1:-1, :-1])
    assert not used.vs[:, 0].any() and not used.vs[..., -1].any()
    np.testing.assert_array_equal(K.phase_velocities(*args, workers=2), full.phase)


def test_identical_columns_are_solved_once_with_the_same_numbers():
    vs = np.repeat(REF["vs"][:, :1, :1], 3, axis=1).repeat(4, axis=2)  # every column the same
    k = K.depth_kernels(REF["depths"], vs, REF["periods"], REF["minthk"], interior_only=False)
    np.testing.assert_array_equal(k.phase, np.broadcast_to(REF["phase"][:, :1, :1], k.phase.shape))
    one = K.depth_kernels(REF["depths"], REF["vs"][:, :1, :1], REF["periods"], REF["minthk"], interior_only=False)
    np.testing.assert_array_equal(k.vs, np.broadcast_to(one.vs, k.vs.shape))
