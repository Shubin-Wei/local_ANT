"""Geometry of the stations: distances on the WGS84 ellipsoid, local coordinates in km."""

import numpy as np

EARTH_RADIUS_KM = 6371.0  # the sphere DSurfTomo computes its distances on (sphere_km)
WGS84_A_KM, WGS84_F = 6378.137, 1.0 / 298.257223563


def geodesic_km(lon1, lat1, lon2, lat2):
    """Distance (km) and azimuth (deg) from point 1 to point 2 on the WGS84 ellipsoid (Vincenty 1975, inverse
    problem).  Station coordinates are geodetic, and on a sphere the distance would be off by up to 0.6 % with a
    sign that depends on the azimuth of the path."""
    lon1, lat1, lon2, lat2 = np.broadcast_arrays(*[np.asarray(v, float) for v in (lon1, lat1, lon2, lat2)])
    a, f = WGS84_A_KM, WGS84_F
    b = a * (1 - f)
    U1, U2 = np.arctan((1 - f) * np.tan(np.radians(lat1))), np.arctan((1 - f) * np.tan(np.radians(lat2)))
    L = np.radians((lon2 - lon1 + 180.0) % 360.0 - 180.0)
    s1, c1, s2, c2 = np.sin(U1), np.cos(U1), np.sin(U2), np.cos(U2)
    lam = L.copy()
    with np.errstate(divide="ignore", invalid="ignore"):
        for _ in range(30):
            sl, cl = np.sin(lam), np.cos(lam)
            ss = np.hypot(c2 * sl, c1 * s2 - s1 * c2 * cl)
            cs = s1 * s2 + c1 * c2 * cl
            sig = np.arctan2(ss, cs)
            sa = np.where(ss > 0, c1 * c2 * sl / ss, 0.0)
            ca2 = 1 - sa**2
            c2m = np.where(ca2 > 0, cs - 2 * s1 * s2 / ca2, 0.0)
            C = f / 16 * ca2 * (4 + f * (4 - 3 * ca2))
            new = L + (1 - C) * f * sa * (sig + C * ss * (c2m + C * cs * (-1 + 2 * c2m**2)))
            done = np.max(np.abs(new - lam), initial=0.0) < 1e-13
            lam = new
            if done:
                break
        u2 = ca2 * (a * a - b * b) / (b * b)
        A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
        B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
        dsig = B * ss * (c2m + B / 4 * (cs * (-1 + 2 * c2m**2) - B / 6 * c2m * (-3 + 4 * ss**2) * (-3 + 4 * c2m**2)))
        dist = np.where(ss > 0, b * A * (sig - dsig), 0.0)
        az = np.degrees(np.arctan2(c2 * np.sin(lam), c1 * s2 - s1 * c2 * np.cos(lam)))
    if dist.ndim == 0:
        return float(dist), float(az % 360.0)
    return dist, az % 360.0


def sphere_km(lon1, lat1, lon2, lat2):
    """Great-circle distance (km) on a sphere of 6371 km: the distance DSurfTomo computes from the coordinates it
    is given.  Only for reproducing DSurfTomo's own travel times (travel time = this distance / velocity)."""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(np.asarray(lon2, float) - np.asarray(lon1, float))
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2))


def km_per_degree(lat0):
    """(km per degree of longitude, km per degree of latitude) on the WGS84 ellipsoid at latitude lat0."""
    e2 = WGS84_F * (2 - WGS84_F)
    w = np.sqrt(1 - e2 * np.sin(np.radians(lat0)) ** 2)
    return (np.radians(1.0) * WGS84_A_KM / w * np.cos(np.radians(lat0)), np.radians(1.0) * WGS84_A_KM * (1 - e2) / w**3)


def local_xy_km(lon, lat, lat0=None, lon0=None):
    """Local x (east) / y (north) in km around (lon0, lat0) (defaults: the mean position), with the WGS84 radii of
    curvature at lat0; for geometry and figures, not for the distances of the measurement (geodesic_km)."""
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    lat0 = float(np.mean(lat)) if lat0 is None else lat0
    lon0 = float(np.mean(lon)) if lon0 is None else lon0
    kx, ky = km_per_degree(lat0)
    return (lon - lon0) * kx, (lat - lat0) * ky


def nearest_neighbour_km(lon, lat):
    """Distance of every station to its nearest neighbour (km), with a k-d tree on local x/y."""
    from scipy.spatial import cKDTree

    x, y = local_xy_km(lon, lat)
    if len(x) < 2:
        return np.array([np.nan] * len(x))
    d, _ = cKDTree(np.c_[x, y]).query(np.c_[x, y], k=2)
    return d[:, 1]


def nice_number(x, mode="nearest", steps=(1.0, 2.0, 2.5, 5.0)):
    """Round x > 0 to the 1-2-2.5-5 x 10^k sequence (mode nearest | up | down)."""
    if not x > 0:
        return float(x)
    e = np.floor(np.log10(x))
    cand = np.array(sorted({s * 10.0 ** (e + k) for k in (-1, 0, 1) for s in steps}))
    if mode == "up":
        return float(cand[cand >= x * (1 - 1e-9)][0])
    if mode == "down":
        return float(cand[cand <= x * (1 + 1e-9)][-1])
    return float(cand[np.argmin(np.abs(np.log(cand / x)))])
