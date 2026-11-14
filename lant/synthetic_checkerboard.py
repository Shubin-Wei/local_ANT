"""Checkerboard truth and straight-path integrals for synthetic noise (not the inversion forward model).

The local slowness is interpolated between the background and the +/- amplitude 1-D profiles, as in the
synthetic correlation example. Plane waves acquire this delay while crossing the square array domain. This
kinematic approximation omits refraction, scattering and amplitude focusing. The truth and its limits are written
beside the generated records; no processing or quality-control rule is changed to recover it.
"""

import numpy as np

AMPLITUDE = 0.06
CELLS = 3  # cell width = one third of the synthetic square; an illustrative engineering choice
NQUAD = 128
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(NQUAD)


def pattern(x, y, cell):
    return np.sin(np.pi * np.asarray(x) / cell) * np.sin(np.pi * np.asarray(y) / cell)


def path_weights(x, y, az, size, cell):
    """Integrals (km) of the positive and negative checker pattern on each incoming wave's path to (x,y).
    Integrate only the intersection with [0,size]^2, behind the station. Azimuths are radians clockwise from north.
    Splitting by sign permits the same three-profile slowness interpolation as examples/synthetic/make_stack.py.
    """
    az = np.atleast_1d(az)
    lower, upper = np.zeros_like(az), np.full_like(az, np.inf)
    for position, direction in ((x, np.sin(az)), (y, np.cos(az))):
        parallel = np.abs(direction) < 1e-12
        safe = np.where(parallel, 1.0, direction)
        a, b = (position - size) / safe, position / safe
        outside = position < 0 or position > size
        lower = np.maximum(lower, np.where(parallel, np.inf if outside else -np.inf, np.minimum(a, b)))
        upper = np.minimum(upper, np.where(parallel, -np.inf if outside else np.inf, np.maximum(a, b)))
    valid = np.isfinite(lower) & np.isfinite(upper) & (upper > lower)
    lower = np.where(valid, lower, 0.0)
    length = np.where(valid, upper - lower, 0.0)
    s = lower[:, None] + 0.5 * length[:, None] * (_NODES + 1)
    values = pattern(x - s * np.sin(az[:, None]), y - s * np.cos(az[:, None]), cell)
    positive = 0.5 * length * (np.maximum(values, 0) @ _WEIGHTS)
    negative = 0.5 * length * (np.maximum(-values, 0) @ _WEIGHTS)
    return positive, negative


def model(description, lon, lat, depths, background_depths, background_vs):
    """3-D Vs truth at requested nodes, using the coordinates and amplitude saved with the raw example."""
    kx, ky = description["km_per_degree"]
    x = (np.asarray(lon) - description["origin_lon"]) * kx
    y = (np.asarray(lat) - description["origin_lat"]) * ky
    inside = (
        (x[:, None] >= 0)
        & (x[:, None] <= description["size_km"])
        & (y[None, :] >= 0)
        & (y[None, :] <= description["size_km"])
    )
    anomaly = description["amplitude"] * pattern(x[:, None], y[None, :], description["cell_km"]) * inside
    vs = np.interp(depths, background_depths, background_vs)
    return (1 + anomaly[:, :, None]) * vs
