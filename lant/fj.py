"""Array phase-velocity spectrum by the frequency-Bessel (F-J) transform (Wang, Wu & Chen 2019).

For a laterally uniform medium and isotropic noise the cross-correlation spectrum of two stations r apart is
C(r, w) ~ J0(w r / c(w)) (Aki 1957).  The F-J spectrum

    I(w, c) = integral C(r, w) J0(w r / c) r dr

is largest where c is the phase velocity of a mode at w.  As in the open-source CC-FJpy (Li et al. 2021) for noise
data, C is the real part of the spectrum of the stacked correlation, and the integral over the pairs is the
trapezoidal rule over the distances in increasing order; double precision.  Every pair with data is used.
"""

import numpy as np
from scipy.special import j0


def spectrum_at(x, lag, f):
    """Real part of the spectrum of the two-sided correlations x [npair, nlag] (lag in s, centred on 0) at the
    frequencies f: the cosine transform (the sine part is the antisymmetric part, not used)."""
    return x.astype(np.float64) @ np.cos(2 * np.pi * np.outer(lag, f))


def fj_integral(X, r, k):
    """I(k) = sum over the pairs of X J0(k r) r dr with the trapezoidal weights of the sorted distances r."""
    o = np.argsort(r)
    r, X = np.asarray(r, np.float64)[o], np.asarray(X, np.float64)[o]
    w = np.empty(len(r))
    w[1:-1] = 0.5 * (r[2:] - r[:-2])
    w[0], w[-1] = 0.5 * (r[1] - r[0]), 0.5 * (r[-1] - r[-2])
    a = X * r * w
    out = np.zeros(len(k))
    for i0 in range(0, len(r), 2000):  # blocks of pairs: bounded memory
        out += a[i0 : i0 + 2000] @ j0(np.outer(r[i0 : i0 + 2000], k))
    return out


def image(ccf, dist, fs, freqs, c):
    """|I(f, c)| [nfreq, nc] of the correlations ccf [npair, 2L+1] at distances dist (km)."""
    L = (ccf.shape[1] - 1) // 2
    lag = np.arange(-L, L + 1) / fs
    X = spectrum_at(ccf, lag, freqs)
    out = np.zeros((len(freqs), len(c)))
    for i, f in enumerate(freqs):
        out[i] = np.abs(fj_integral(X[:, i], dist, 2 * np.pi * f / c))
    return out


def pick(img, c, freqs, band=None):
    """Largest value of each frequency inside `band` (Hz; None: every frequency); a peak on the edge of the
    velocity scan is not a ridge and is not picked.  -> picks (NaN where none)."""
    out = np.full(len(freqs), np.nan)
    for i, f in enumerate(freqs):
        if band is not None and not (band[0] * (1 - 1e-9) <= f <= band[1] * (1 + 1e-9)):
            continue
        j = int(np.argmax(img[i]))
        if 0 < j < len(c) - 1 and img[i, j] > 0:
            out[i] = c[j]
    return out


def smooth_curve(freqs, picks, degree=3):
    """Cubic polynomial of ln c in ln f through the picks (fewer points: lower degree); defined over the range of
    the picks only (NaN outside: no extrapolation)."""
    ok = np.isfinite(picks)
    if ok.sum() < 3:
        raise RuntimeError("the F-J spectrum gave %d picks; at least 3 are needed for a reference curve" % ok.sum())
    x = np.log(freqs)
    p = np.polyfit(x[ok], np.log(picks[ok]), min(degree, ok.sum() - 2))
    out = np.exp(np.polyval(p, x))
    out[(freqs < freqs[ok].min() * (1 - 1e-9)) | (freqs > freqs[ok].max() * (1 + 1e-9))] = np.nan
    return out
