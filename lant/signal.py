"""Shared numerics: frequency axes, branches, notch, Gaussian filter bank, ridge tracking, peaks.

These functions are used unchanged by the cross-correlation, the notch and the dispersion steps, so that every
step applies exactly the same operation (for example the same notch gain in the stack and in the measurement).
"""
import numpy as np

EARTH_RADIUS_KM = 6371.0      # spherical distances, as DSurfTomo computes them internally


def log_freqs(fmin=0.5, fmax=10.0, per_octave=12):
    """Log-spaced frequencies fmin * 2**(k / per_octave) up to fmax (inclusive within rounding)."""
    n = int(round(np.log2(fmax / fmin) * per_octave)) + 1
    return fmin * 2 ** (np.arange(n) / per_octave)


def branches(ccf, fs, maxlag):
    """Split two-sided CCFs [npair, 2L+1] into causal / acausal (time reversed) branches."""
    L = int(round(maxlag * fs))
    pos = ccf[:, L:]
    neg = ccf[:, :L + 1][:, ::-1]
    return pos, neg


def notch_lines(x, fs, lines, half_width=0.06):
    """Remove narrow spectral lines from time series [..., nt] with a raised-cosine gain around each line."""
    if not len(lines):
        return x
    nt = x.shape[-1]
    nfft = int(2 ** np.ceil(np.log2(nt * 4)))
    X = np.fft.rfft(x, n=nfft, axis=-1)
    f = np.fft.rfftfreq(nfft, 1 / fs)
    g = np.ones_like(f)
    for fl in lines:
        d = np.abs(f - fl) / half_width
        g *= np.where(d < 1, 0.5 - 0.5 * np.cos(np.pi * d), 1.0)
    return np.fft.irfft(X * g, n=nfft, axis=-1)[..., :nt]


def notch_table(x, fs, table):
    """Apply a notch table [[centre_Hz, half_width_Hz], ...] line by line (each with its own half-width)."""
    x = np.asarray(x, np.float64)
    for f0, hw in np.atleast_2d(table):
        x = notch_lines(x, fs, [f0], half_width=hw)
    return x


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
    if prior is not None:
        rel = np.abs(vel[None, :] - prior[:, None]) / prior[:, None]
        score = score - 20.0 * np.clip(rel - prior_tol, 0, None)
    if valid is not None:
        score = np.where(valid[:, None], score, 0.0)
    dv = np.abs(vel[:, None] - vel[None, :]) / vel[None, :]
    trans = -(dv / smooth) ** 2
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


def best_lag(a, b, maxs):
    """Lag (samples) maximising |corr(a(t), b(t+lag))|, and that correlation."""
    best, bl = 0.0, 0
    na = np.sqrt((a * a).sum()) + 1e-30
    for s in range(-maxs, maxs + 1):
        bb = np.roll(b, -s)
        c = (a * bb).sum() / (na * (np.sqrt((bb * bb).sum()) + 1e-30))
        if abs(c) > abs(best):
            best, bl = c, s
    return bl, best


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


def contiguous(mask, min_run):
    """Keep only runs of at least min_run True values."""
    out = np.zeros_like(mask)
    for k, j in runs_of(mask):
        if j - k >= min_run:
            out[k:j] = True
    return out


def geodesic_km(lon1, lat1, lon2, lat2):
    """Great-circle distance (km, sphere of 6371 km) and azimuth (deg) from point 1 to point 2."""
    r = EARTH_RADIUS_KM
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(np.asarray(lon2) - np.asarray(lon1))
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    dist = 2 * r * np.arcsin(np.sqrt(a))
    az = np.degrees(np.arctan2(np.sin(dl) * np.cos(p2),
                               np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl)))
    return dist, az % 360.0


def local_xy_km(lon, lat, lat0=None, lon0=None):
    """Equirectangular x (east) / y (north) in km around (lon0, lat0) (defaults: the mean position)."""
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    lat0 = float(np.mean(lat)) if lat0 is None else lat0
    lon0 = float(np.mean(lon)) if lon0 is None else lon0
    k = np.radians(1.0) * EARTH_RADIUS_KM
    return (lon - lon0) * k * np.cos(np.radians(lat0)), (lat - lat0) * k
