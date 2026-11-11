"""Shared numerics: frequency axes, branches, Gaussian filter bank, ridge tracking, peaks.

These functions are used unchanged by the reference, the dispersion and the figure code, so that every place
applies exactly the same operation.
"""

import numpy as np

EARTH_RADIUS_KM = 6371.0  # the sphere DSurfTomo computes its distances on (sphere_km)
WGS84_A_KM, WGS84_F = 6378.137, 1.0 / 298.257223563


def log_freqs(fmin=0.5, fmax=10.0, per_octave=12):
    """Log-spaced frequencies fmin * 2**(k / per_octave) up to fmax (inclusive within rounding)."""
    n = int(round(np.log2(fmax / fmin) * per_octave)) + 1
    return fmin * 2 ** (np.arange(n) / per_octave)


def branches(ccf, fs, maxlag):
    """Split two-sided CCFs [npair, 2L+1] into causal / acausal (time reversed) branches."""
    L = int(round(maxlag * fs))
    pos = ccf[:, L:]
    neg = ccf[:, : L + 1][:, ::-1]
    return pos, neg


def analytic_bank(x, fs, freqs, alpha, taper_s=1.0):
    """Narrow-band analytic signals with a Gaussian filter of constant relative width.

    x      [npair, nt] one-sided branch (t = 0 .. maxlag)
    alpha  scalar or [npair]:  H = exp(-alpha ((f - fc) / fc)^2)
    return complex [npair, nfreq, nt]
    """
    npair, nt = x.shape
    w = np.ones(nt)
    ntap = int(taper_s * fs)
    w[-ntap:] = 0.5 + 0.5 * np.cos(np.pi * np.arange(ntap) / ntap)
    nfft = int(2 ** np.ceil(np.log2(max(4 * nt, 4096))))
    X = np.fft.fft(x * w, n=nfft, axis=1)
    f = np.fft.fftfreq(nfft, 1 / fs)
    pos = f > 0
    alpha = np.broadcast_to(np.asarray(alpha, float), (npair,))[:, None]
    out = np.empty((npair, len(freqs), nt), np.complex64)
    for k, fc in enumerate(freqs):
        H = np.zeros((npair, nfft))
        H[:, pos] = 2.0 * np.exp(-alpha * ((f[pos] - fc) / fc) ** 2)
        out[:, k] = np.fft.ifft(X * H, axis=1)[:, :nt]
    return out


def sample_at_times(a, fs, times):
    """Linear interpolation of a [..., nt] at fractional sample times [..., m]."""
    nt = a.shape[-1]
    pos = np.clip(times * fs, 0, nt - 1.000001)
    i0 = np.floor(pos).astype(int)
    fr = pos - i0
    lo = np.take_along_axis(a, i0, axis=-1)
    hi = np.take_along_axis(a, i0 + 1, axis=-1)
    return lo * (1 - fr) + hi * fr


def ridge_dp(img, vel, prior=None, prior_tol=0.3, smooth=0.02, valid=None):
    """Viterbi ridge through img [nf, nv] (larger = better) on a velocity grid.

    Penalises velocity jumps between neighbouring frequencies (relative units) and departures from a prior
    curve beyond prior_tol (relative)."""
    nf, nv = img.shape
    score = np.log(np.clip(img, 1e-3, None))
    if prior is not None:  # NaN prior: no pull at that frequency
        with np.errstate(invalid="ignore"):
            rel = np.abs(vel[None, :] - prior[:, None]) / prior[:, None]
        score = score - 20.0 * np.nan_to_num(np.clip(rel - prior_tol, 0, None))
    if valid is not None:
        score = np.where(valid[:, None], score, 0.0)
    dv = np.abs(vel[:, None] - vel[None, :]) / vel[None, :]
    trans = -((dv / smooth) ** 2)
    acc = score[0].copy()
    back = np.zeros((nf, nv), int)
    for k in range(1, nf):
        tot = acc[None, :] + trans
        back[k] = np.argmax(tot, axis=1)
        acc = tot[np.arange(nv), back[k]] + score[k]
    path = np.zeros(nf, int)
    path[-1] = int(np.argmax(acc))
    for k in range(nf - 1, 0, -1):
        path[k - 1] = back[k, path[k]]
    return path


def refine_peak(profile, idx):
    """Parabolic sub-grid refinement of a peak; returns a fractional index."""
    if idx <= 0 or idx >= len(profile) - 1:
        return float(idx)
    y0, y1, y2 = profile[idx - 1], profile[idx], profile[idx + 1]
    den = y0 - 2 * y1 + y2
    if den >= 0:
        return float(idx)
    return idx + 0.5 * (y0 - y2) / den


def local_peak(profile, idx, halfwin):
    lo = max(0, idx - halfwin)
    hi = min(len(profile), idx + halfwin + 1)
    return lo + int(np.argmax(profile[lo:hi]))


def runs_of(mask):
    """Contiguous runs [(start, stop), ...] of True values."""
    out = []
    k = 0
    while k < len(mask):
        if mask[k]:
            j = k
            while j < len(mask) and mask[j]:
                j += 1
            out.append((k, j))
            k = j
        else:
            k += 1
    return out


def gradient_in_runs(y, x):
    """dy/dx inside each run of finite y (central differences, one-sided at the ends of a run); nan elsewhere and
    on runs of a single sample."""
    out = np.full(len(y), np.nan)
    for k0, k1 in runs_of(np.isfinite(y)):
        if k1 - k0 >= 2:
            out[k0:k1] = np.gradient(y[k0:k1], x[k0:k1])
    return out


def contiguous(mask, min_run):
    """Keep only runs of at least min_run True values."""
    out = np.zeros_like(mask)
    for k, j in runs_of(mask):
        if j - k >= min_run:
            out[k:j] = True
    return out


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


def common_rates(rates, fmin_rate=1.0):
    """Whole-number processing rates (Hz) reachable from every sampling rate by an integer decimation factor,
    largest first (whole numbers keep an integer number of samples per second and per hour)."""
    rates = [float(r) for r in rates]
    lo = min(rates)
    out = []
    for q in range(1, int(lo / fmin_rate) + 1):
        f = lo / q
        if abs(f - round(f)) < 1e-6 and all(abs(r / f - round(r / f)) < 1e-6 for r in rates):
            out.append(float(round(f)))
    return out


def octave_bands(freqs):
    """Contiguous one-octave bands covering freqs (the last one stretched to the highest frequency)."""
    f0, f1 = float(np.min(freqs)), float(np.max(freqs))
    edges = [f0]
    while edges[-1] * 2 < f1 * (1 - 1e-9):
        edges.append(edges[-1] * 2)
    if len(edges) > 1 and f1 / edges[-1] < 2**0.5:  # a short last band joins the previous one
        edges.pop()
    edges.append(f1 * (1 + 1e-9))
    return [(round(a, 6), round(b, 6)) for a, b in zip(edges[:-1], edges[1:])]


def band_key(lo, hi):
    return "%.4g-%.4g" % (lo, hi)


def nearest_neighbour_km(lon, lat):
    """Distance of every station to its nearest neighbour (km), with a k-d tree on local x/y."""
    from scipy.spatial import cKDTree

    x, y = local_xy_km(lon, lat)
    if len(x) < 2:
        return np.array([np.nan] * len(x))
    d, _ = cKDTree(np.c_[x, y]).query(np.c_[x, y], k=2)
    return d[:, 1]


def local_xy_km(lon, lat, lat0=None, lon0=None):
    """Local x (east) / y (north) in km around (lon0, lat0) (defaults: the mean position), with the WGS84 radii of
    curvature at lat0; for geometry and figures, not for the distances of the measurement (geodesic_km)."""
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    lat0 = float(np.mean(lat)) if lat0 is None else lat0
    lon0 = float(np.mean(lon)) if lon0 is None else lon0
    kx, ky = km_per_degree(lat0)
    return (lon - lon0) * kx, (lat - lat0) * ky
