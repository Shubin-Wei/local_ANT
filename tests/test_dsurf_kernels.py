"""DSurfTomo depth kernels: the native depthkernel's numbers bit for bit (tests/data/native_depthkernel.npz, made
with the unmodified Fortran objects), and only the work CalSurfG reads."""

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("numba")
from lant.inversion import dsurf_kernels as K  # noqa: E402

REF = np.load(Path(__file__).parent / "data" / "native_depthkernel.npz")


def test_identical_to_native_depthkernel_bit_for_bit():
    k = K.depth_kernels(REF["depths"], REF["vs"], REF["periods"], REF["minthk"], workers=2, interior_only=False)
    np.testing.assert_array_equal(k.phase, REF["phase"])
    for name in ("vs", "vp", "rho"):
        np.testing.assert_array_equal(getattr(k, name), REF["sen_" + name])


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
    np.testing.assert_array_equal(k.vs, np.broadcast_to(REF["sen_vs"][:, :1, :1], k.vs.shape))
