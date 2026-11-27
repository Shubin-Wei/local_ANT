"""Narrow-band filtering and helpers shared by the reference curve, the dispersion measurement and the figures."""

import numpy as np


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


def refine_peak(profile, idx):
    """Parabolic sub-grid refinement of a peak; returns a fractional index."""
    if idx <= 0 or idx >= len(profile) - 1:
        return float(idx)
    y0, y1, y2 = profile[idx - 1], profile[idx], profile[idx + 1]
    den = y0 - 2 * y1 + y2
    if den >= 0:
        return float(idx)
    return idx + 0.5 * (y0 - y2) / den


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
