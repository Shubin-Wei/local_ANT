"""Ray tracing and Frechet weights of the DSurfTomo-compatible forward (rpaths)."""

import numpy as np
import pytest

pytest.importorskip("numba")
from lant.inversion import dsurf_bspline as bs  # noqa: E402
from lant.inversion import dsurf_rays as rays  # noqa: E402
from lant.inversion import dsurf_source as src  # noqa: E402

F = np.float32
EARTH = F(6371)
VELOCITY = 2.0


def ray_weights(fx, fz, nvz=8, nvx=7):
    velv = np.full((nvz + 2, nvx + 2), VELOCITY, F)
    veln = bs.gridder(F, velv, 8, 8)
    gox, goz = F(np.radians(60.0)), F(np.radians(100.0))
    dvx = dvz = F(np.radians(0.01))
    dnx, dnz = dvx / F(8), dvz / F(8)
    sx, sz = gox + F(13.37) * dnx, goz + F(20.81) * dnz
    status, ttn, ttnr, nstsr, goxr, gozr, dnxr, dnzr = src.source_fields(
        F, sx, sz, velv, veln, gox, goz, dvx, dvz, EARTH, 8, 8, 8, 8
    )
    rx, rz = gox + F(fx) * dnx, goz + F(fz) * dnz
    fdm = np.ones((nvz + 2, nvx + 2), F)
    status = rays.rpaths(
        F,
        sx,
        sz,
        rx,
        rz,
        fdm,
        ttn,
        veln,
        gox,
        goz,
        dnx,
        dnz,
        dvx,
        dvz,
        8,
        8,
        EARTH,
        ttnr,
        nstsr,
        goxr,
        gozr,
        dnxr,
        dnzr,
    )
    d = 6371 * np.hypot(float(sx) - float(rx), (float(sz) - float(rz)) * np.sin(float(rx)))
    dpl = 0.5 * min(float(dnx) * 6371, float(dnz) * 6371 * np.sin(float(gox + 48 * dnx)))
    return status, fdm, d, dpl


def test_far_ray_integrates_path_except_closing_segment():
    status, fdm, d, dpl = ray_weights(41.3, 52.7)
    assert status == 0
    # Uniform 1/v^2 and a partition of unity: -sum(fdm) * v^2 is the integrated length, which
    # stops when the ray comes within 2*dpl of the source (the closing segment is not added).
    length = -float(fdm.sum()) * VELOCITY**2
    assert d - 3 * dpl < length < d


def test_receiver_near_source_gives_zero_weights():
    status, fdm, d, dpl = ray_weights(13.9, 21.1)
    assert status == 0 and d < 2 * dpl and not np.any(fdm)


def test_receiver_on_last_node_is_outside():
    status, fdm, d, dpl = ray_weights(48.0, 10.0)
    assert status == rays.RECEIVER_OUTSIDE
