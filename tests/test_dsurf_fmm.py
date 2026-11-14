"""Heap and fast marching of the DSurfTomo-compatible forward."""

import numpy as np
import pytest

pytest.importorskip("numba")
from lant.inversion import dsurf_fmm as fmm  # noqa: E402


def test_heap_pops_in_order_and_tracks_positions():
    rng = np.random.default_rng(0)
    nz, nx = 6, 7
    ttn = rng.random((nz + 1, nx + 1)).astype(np.float32)
    nsts = np.full((nz + 1, nx + 1), -1, np.int64)
    btx = np.zeros(nz * nx + 1, np.int64)
    btz = np.zeros_like(btx)
    ntr = 0
    for iz in range(1, nz + 1):
        for ix in range(1, nx + 1):
            ntr = fmm.addtree(iz, ix, ntr, ttn, nsts, btx, btz)
    # Lower one value and update its position, as fouds2 + updtree do.
    ttn[3, 4] = -1
    fmm.updtree(3, 4, ttn, nsts, btx, btz)
    popped = []
    while ntr > 0:
        for position in range(1, ntr + 1):
            assert nsts[btz[position], btx[position]] == position
        iz, ix = btz[1], btx[1]
        popped.append(ttn[iz, ix])
        nsts[iz, ix] = 0
        ntr = fmm.downtree(ntr, ttn, nsts, btx, btz)
    assert len(popped) == nz * nx and popped[0] == -1
    assert np.all(np.diff(popped) >= 0)


def _plain_field(real, velocity, nz=41, nx=33):
    veln = np.full((nz + 1, nx + 1), velocity, real)
    ttn = np.zeros_like(veln)
    nsts = np.full(veln.shape, -1, np.int64)
    btx = np.zeros(veln.size + 1, np.int64)
    btz = np.zeros_like(btx)
    gox, goz = real(np.radians(59.9)), real(np.radians(100.0))
    dnx = dnz = real(np.radians(0.001))
    return veln, ttn, nsts, btx, btz, gox, goz, dnx, dnz


@pytest.mark.parametrize("real", [np.float32, np.float64])
def test_uniform_medium_times_follow_distance(real):
    veln, ttn, nsts, btx, btz, gox, goz, dnx, dnz = _plain_field(real, 2.0)
    scx, scz = gox + real(16.3) * dnx, goz + real(20.6) * dnz
    status = fmm.travel(real, scx, scz, 0, veln, ttn, nsts, btx, btz, gox, goz, dnx, dnz, real(6371), 1, 1, 1, 1)
    assert status == 0 and np.all(nsts[1:, 1:] == 0)
    # Far from the source the front is nearly plane: the gradient along x approaches 1/v. (Without
    # source refinement, urg=0, absolute times near the source carry errors of several per cent.)
    row = ttn[21, 1:]
    slope = (row[0] - row[8]) / (8 * float(dnx) * 6371)
    assert np.all(np.diff(row[:17]) < 0) and np.all(np.diff(row[18:]) > 0)
    assert abs(slope * 2.0 - 1) < 0.01


def test_source_outside_is_reported():
    veln, ttn, nsts, btx, btz, gox, goz, dnx, dnz = _plain_field(np.float32, 2.0)
    status = fmm.travel(
        np.float32, gox - dnx * 2, goz, 0, veln, ttn, nsts, btx, btz, gox, goz, dnx, dnz, np.float32(6371), 1, 1, 1, 1
    )
    assert status == fmm.SOURCE_OUTSIDE


def test_refined_march_stops_at_inner_edge():
    veln, ttn, nsts, btx, btz, gox, goz, dnx, dnz = _plain_field(np.float32, 2.0, 17, 17)
    scx, scz = gox + np.float32(8.5) * dnx, goz + np.float32(8.5) * dnz
    # The refined box is not at the model edge on any side: the march ends at its first edge node.
    status = fmm.travel(
        np.float32, scx, scz, 1, veln, ttn, nsts, btx, btz, gox, goz, dnx, dnz, np.float32(6371), 5, 30, 5, 30
    )
    assert status == 0
    alive = nsts[1:, 1:] == 0
    assert 0 < alive.sum() < alive.size and np.any(nsts[1:, 1:] > 0)
