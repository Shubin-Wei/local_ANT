"""Refined source fields and receiver times of the DSurfTomo-compatible forward."""

import numpy as np
import pytest

pytest.importorskip("numba")
from lant.inversion import dsurf_bspline as bs  # noqa: E402
from lant.inversion import dsurf_source as src  # noqa: E402

F = np.float32
EARTH = F(6371)


def uniform(nvz=8, nvx=7, velocity=2.0):
    velv = np.full((nvz + 2, nvx + 2), velocity, F)
    veln = bs.gridder(F, velv, 8, 8)
    gox, goz = F(np.radians(60.0)), F(np.radians(100.0))
    dvx = dvz = F(np.radians(0.01))
    return velv, veln, gox, goz, dvx, dvz


def distance(scx, scz, rcx, rcz):
    return 6371 * np.hypot(float(scx) - float(rcx), (float(scz) - float(rcz)) * np.sin(float(rcx)))


def test_refined_source_field_is_complete_and_accurate():
    velv, veln, gox, goz, dvx, dvz = uniform()
    dnx, dnz = dvx / F(8), dvz / F(8)
    sx, sz = gox + F(13.37) * dnx, goz + F(20.81) * dnz
    status, ttn, ttnr, nstsr, goxr, gozr, dnxr, dnzr = src.source_fields(
        F, sx, sz, velv, veln, gox, goz, dvx, dvz, EARTH, 8, 8, 8, 8
    )
    assert status == 0 and ttnr.shape == (8 * 16 + 2, 8 * 16 + 2)
    assert dnxr == dvx / F(64) and np.any(nstsr[1:, 1:] == 0)
    errors = []
    for fx, fz in [(1.2, 3.3), (40.7, 50.1), (5.5, 2.25), (13.6, 20.9)]:
        rx, rz = gox + F(fx) * dnx, goz + F(fz) * dnz
        status, t = src.srtimes(F, sx, sz, rx, rz, ttn, veln, gox, goz, dnx, dnz, EARTH)
        assert status == 0
        errors.append(float(t) * 2.0 / distance(sx, sz, rx, rz) - 1)
    # With the 1/64 source refinement the times are within a fraction of a per cent of d/v.
    assert np.max(np.abs(errors)) < 5e-3


def test_near_source_rule_uses_mean_velocity():
    velv, veln, gox, goz, dvx, dvz = uniform(velocity=2.5)
    dnx, dnz = dvx / F(8), dvz / F(8)
    sx, sz = gox + F(13.2) * dnx, goz + F(20.2) * dnz
    ttn = np.zeros_like(veln)  # never read on this branch
    rx, rz = gox + F(13.7) * dnx, goz + F(20.6) * dnz  # same coarse cell as the source
    status, t = src.srtimes(F, sx, sz, rx, rz, ttn, veln, gox, goz, dnx, dnz, EARTH)
    assert status == 0 and abs(float(t) - distance(sx, sz, rx, rz) / 2.5) < 1e-6


def test_receiver_outside_is_reported():
    velv, veln, gox, goz, dvx, dvz = uniform()
    dnx, dnz = dvx / F(8), dvz / F(8)
    status, _ = src.srtimes(F, gox, goz, gox - F(3) * dnx, goz, veln, veln, gox, goz, dnx, dnz, EARTH)
    assert status == src.RECEIVER_OUTSIDE


def test_march_stopped_after_the_receivers_reads_the_same_values():
    """With the receivers given, the coarse march stops early; receiver times and ray weights are unchanged."""
    from lant.inversion import dsurf_bspline as bs
    from lant.inversion import dsurf_rays as rays

    F = np.float32
    earth = F(6371)
    rng = np.random.default_rng(3)
    nvz, nvx = 12, 11
    velv = (2.0 + 0.3 * rng.random((nvz + 2, nvx + 2))).astype(F)
    veln = bs.gridder(F, velv, 8, 8)
    gox, goz = F(np.radians(60.0)), F(np.radians(100.0))
    dvx = dvz = F(np.radians(0.01))
    dnx, dnz = dvx / F(8), dvz / F(8)
    sx, sz = gox + F(30.3) * dnx, goz + F(35.7) * dnz
    rx = gox + np.array([20.4, 45.1, 33.3], F) * dnx
    rz = goz + np.array([41.2, 28.6, 52.9], F) * dnz
    full = src.source_fields(F, sx, sz, velv, veln, gox, goz, dvx, dvz, earth, 8, 8, 8, 8)
    cut = src.source_fields(F, sx, sz, velv, veln, gox, goz, dvx, dvz, earth, 8, 8, 8, 8, rx, rz)
    assert np.count_nonzero(cut[1]) < np.count_nonzero(full[1])  # the march did stop early
    for r in range(rx.size):
        t = [src.srtimes(F, sx, sz, rx[r], rz[r], f[1], veln, gox, goz, dnx, dnz, earth) for f in (full, cut)]
        assert t[0] == t[1]
        w = []
        for f in (full, cut):
            fdm = np.zeros((nvz + 2, nvx + 2), F)
            st = rays.rpaths(F, sx, sz, rx[r], rz[r], fdm, f[1], veln, gox, goz, dnx, dnz, dvx, dvz, 8, 8, earth,
                             f[2], f[3], f[4], f[5], f[6], f[7])  # fmt: skip
            w.append((st, fdm))
        assert w[0][0] == w[1][0] and np.array_equal(w[0][1], w[1][1])
