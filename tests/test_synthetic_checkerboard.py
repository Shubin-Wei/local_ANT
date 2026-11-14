"""Independent straight-path checks and the saved synthetic truth."""

import numpy as np

from lant import synthetic_checkerboard as checker

_trapezoid = getattr(np, "trapezoid", None) or np.trapz  # NumPy 2 renamed trapz (removed in 2.4)


def test_axis_path_has_the_integral_of_a_half_sine():
    positive, negative = checker.path_weights(1.5, 1.5, [np.pi / 2], 9.0, 3.0)
    assert np.allclose(positive, 3 / np.pi, atol=1e-10)
    assert np.allclose(negative, 0, atol=1e-10)


def test_incoming_paths_match_dense_quadrature():
    angles = np.r_[np.arange(0, 2 * np.pi, np.pi / 8), 1e-14, np.pi / 2 + 1e-14]
    got = checker.path_weights(4.2, 6.3, angles, 9.0, 3.0)
    s = np.linspace(0, 20, 100001)
    for i, az in enumerate(angles):
        x, y = 4.2 - s * np.sin(az), 6.3 - s * np.cos(az)
        inside = (x >= 0) & (x <= 9) & (y >= 0) & (y <= 9)
        values = np.sin(np.pi * x / 3) * np.sin(np.pi * y / 3) * inside
        assert abs(got[0][i] - _trapezoid(np.maximum(values, 0), s)) < 0.001
        assert abs(got[1][i] - _trapezoid(np.maximum(-values, 0), s)) < 0.001


def test_ray_outside_square_has_no_checkerboard_delay():
    p, n = checker.path_weights(-1.0, 2.0, [0.0, -np.pi / 2], 9.0, 3.0)
    assert p[0] == n[0] == 0
    assert p[1] + n[1] > 0  # the other incoming ray crosses the square before reaching this station


def test_truth_includes_both_signs_and_background_outside_square():
    desc = dict(km_per_degree=[1, 1], origin_lon=0, origin_lat=0, size_km=9, amplitude=0.06, cell_km=3)
    v = checker.model(desc, [-1, 1.5, 4.5], [1.5], [0, 1], [0, 1], [2, 3])
    assert np.allclose(v[:, 0], [[2, 3], [2.12, 3.18], [1.88, 2.82]])
