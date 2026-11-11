"""F-J spectrum of an array (frequency-Bessel transform, Wang et al. 2019) and the picking of its fundamental-mode
ridge, for the reference step (lant/reference.py).

The spectrum (fj_spectrum, double precision) is computed on a scan grid from the lowest usable frequency (the
longest pair = 2 wavelengths at the slowest velocity of data.velocity_range_kms) to dispersion.fmax_hz, velocities
fj_cmin - fj_cmax; one frequency at a time with the pairs not longer than fj_rmax_lambda wavelengths of fj_c_nominal
(auto: of the curve picked in a first pass that takes every pair).
A pick is "clear" when the second-highest peak is below 0.7 of the largest and >= 30 pairs were used;
  fj_mode max_peak (areal arrays): the largest peak (assumes the fundamental mode is strongest); a cubic in log f
      fitted to the longest run of clear picks that follow each other smoothly (fj_fit_fmax, when given, caps the
      run), held at its end values outside;
  fj_mode track_fundamental (lines, where higher modes or off-line noise can dominate): the fundamental ridge is
      followed from a start in track_start_band_hz (auto: the tracking is tried from every octave and the start
      that follows the ridge over the most frequencies is taken, the slower ridge between equally long ones; at
      each frequency the strongest peak whose change from the last pick is <= 5 % per 1/12-octave step since that
      pick (a limit on d ln c / d ln f), at most 3 % rise with frequency; picks only inside
      data.velocity_range_kms and, on a line, only where the wavelength is at least twice the mean station
      spacing: the spatial Nyquist limit of a sampled line, Landau 1967, beyond which the spectrum has aliased
      ridges), cubic in log f, made non-increasing (reference.normal_dispersion), constant outside; fewer than 5
      frequencies: an error that says what to do.
"""

import math

import h5py
import numpy as np

from .signal import log_freqs


def scan_grid(cfg):
    """The F-J scan frequencies: from the lowest usable frequency (correlation.fmin_hz undone by the filter margin)
    to dispersion.fmax_hz, on the same grid."""
    import math

    from .auto import FILTER_EDGE, _grid_up

    d = cfg["dispersion"]
    dlt = math.sqrt(FILTER_EDGE / float(d["alpha"]))
    f0 = _grid_up(cfg["correlation"]["fmin_hz"] / (1 - dlt) * (1 - 1e-6), int(d["per_octave"]))
    f = log_freqs(f0, d["fmax_hz"], d["per_octave"])
    return f[f <= d["fmax_hz"] + 1e-9]


def fj_nfft(value, L):
    """reference.fj_nfft; auto: a power of 2 >= 4 x the two-sided lag window, at least 8192."""
    if value == "auto":
        return 1 << int(np.ceil(np.log2(max(8192, 4 * (2 * L + 1)))))
    if value < 2 * L + 1:
        raise RuntimeError("reference.fj_nfft %d is shorter than the lag window (%d samples)" % (value, 2 * L + 1))
    return int(value)


def min_fj_pairs(n_pairs):
    """Pairs a frequency needs for an F-J spectrum: 30, for small arrays a quarter of the pairs, at least 10."""
    return int(min(30, max(10, n_pairs // 4)))


def fj_spectrum(X, r, k):
    """Frequency-Bessel spectrum I(k) = integral of C(r) J0(k r) r dr (Wang et al. 2019, eq. 4) for the real
    spectral values X [npair] at the distances r [npair] (any order) and the wavenumbers k [nk]: the trapezoidal
    rule over the sorted distances, in double precision (the sum stays exact for any number of pairs)."""
    from scipy.special import j0

    o = np.argsort(r)
    r, X = np.asarray(r, np.float64)[o], np.asarray(X, np.float64)[o]
    w = np.empty(len(r))
    w[1:-1] = 0.5 * (r[2:] - r[:-2])
    w[0], w[-1] = 0.5 * (r[1] - r[0]), 0.5 * (r[-1] - r[-2])
    a = X * r * w
    out = np.zeros(len(k))
    for i0 in range(0, len(r), 2000):  # pairs in blocks: memory stays small for any array
        out += a[i0 : i0 + 2000] @ j0(np.outer(r[i0 : i0 + 2000], k))
    return out


def fj_image(cfg, stack_path, freqs=None, select=None, c_limit=None):
    """F-J image [nfreq, nc] of the pairs of the stack (select: boolean mask over all pairs of the file)."""
    global MIN_FJ_PAIRS
    r = cfg["reference"]
    with h5py.File(stack_path, "r") as h:
        dist = h["dist_km"][:]
        use = h["nseg"][:] > 0
        if select is not None:
            use &= select
        idx = np.where(use)[0]
        x = h["stack"][idx].astype(np.float64) if len(idx) < len(use) else h["stack"][:].astype(np.float64)
        fs, maxlag = float(h.attrs["fs"]), float(h.attrs["maxlag"])
    dist = dist[idx]
    L = int(round(maxlag * fs))
    nfft = fj_nfft(r["fj_nfft"], L)
    freqs = scan_grid(cfg) if freqs is None else freqs
    fall = np.fft.rfftfreq(nfft, 1 / fs)
    kbin = np.array([int(np.argmin(np.abs(fall - f))) for f in freqs])
    lag = np.arange(-L, L + 1) / fs
    # the real part of the spectrum at the scan frequencies only (the CCFs are real: a cosine transform)
    X = x @ np.cos(2 * np.pi * np.outer(lag, fall[kbin]))
    c = np.arange(r["fj_cmin"], r["fj_cmax"] + 1e-9, r["fj_dc"])
    img = np.zeros((len(freqs), len(c)))
    npair = np.zeros(len(freqs), int)
    MIN_FJ_PAIRS = min_fj_pairs(len(idx))
    # velocity that turns fj_rmax_lambda into a distance: a number, or (frequencies, curve) of an earlier pass
    cl = (
        np.full(len(freqs), float(r["fj_c_nominal"]))
        if c_limit is None
        else np.interp(np.log(freqs), np.log(c_limit[0]), c_limit[1])
    )
    for i, (f, kb) in enumerate(zip(freqs, kbin)):
        m = dist * f / cl[i] <= r["fj_rmax_lambda"]
        npair[i] = m.sum()
        if m.sum() < MIN_FJ_PAIRS:
            continue
        img[i] = np.abs(fj_spectrum(X[m, i], dist[m], 2 * np.pi * fall[kb] / c))
    return freqs, c, img, npair


def _peaks(y, thr):
    return np.where((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]) & (y[1:-1] >= thr))[0] + 1


CLEAR_RATIO = 0.7  # second-highest F-J peak below this fraction of the largest: a clear pick


MIN_FJ_PAIRS = 30  # set by fj_image for the array at hand (min_fj_pairs)


def fj_picks(img, cgrid, npair):
    """Largest peak per frequency, the contrast 1 - second / largest peak, and whether the pick is clear."""
    nf = img.shape[0]
    pick = np.full(nf, np.nan)
    contrast = np.zeros(nf)
    for i in range(nf):
        y = img[i]
        if not y.any() or npair[i] < MIN_FJ_PAIRS:
            continue
        pk = np.where((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]))[0] + 1
        if not pk.size:
            continue
        hs = np.sort(y[pk])[::-1]
        j = pk[np.argmax(y[pk])]
        pick[i] = cgrid[j]
        contrast[i] = 1.0 - (hs[1] / hs[0] if hs.size > 1 and hs[0] > 0 else 0.0)
    clear = np.isfinite(pick) & (contrast >= 1 - CLEAR_RATIO)
    return pick, contrast, clear


def trusted_run(pick, clear, max_jump=0.10, max_gap=2):
    """Indices of the longest run of clear picks that follow each other smoothly (gaps of at most max_gap
    frequencies, |d ln c| <= max_jump per step): scattered picks outside the ridge are left out."""
    k = np.where(clear)[0]
    runs, cur = [], [k[0]] if k.size else []
    for a, b in zip(k[:-1], k[1:]):
        if b - a <= max_gap + 1 and abs(np.log(pick[b] / pick[a])) <= max_jump * (b - a):
            cur.append(b)
        else:
            runs.append(cur)
            cur = [b]
    if cur:
        runs.append(cur)
    return np.array(max(runs, key=len), int) if runs else np.zeros(0, int)


def pick_max(cfg, freqs, c, img, npair):
    r = cfg["reference"]
    cfj, contrast, clear = fj_picks(img, c, npair)
    clear &= (cfj > c[0]) & (cfj < c[-1])  # a peak on the edge of the scan is not a ridge
    if r["fj_fit_fmax"] != "auto":
        clear &= freqs <= r["fj_fit_fmax"] * (1 + 1e-9)
    k = trusted_run(cfj, clear)
    if k.size < 5:
        raise RuntimeError(
            "the F-J image has only %d clear picks that follow each other (work/04_reference/fj.png): check the "
            "stack (too few pairs or days?) or borrow a reference (reference.source = file)" % k.size
        )
    good = np.zeros(len(freqs), bool)
    good[k] = True
    x = np.log(freqs)
    p = np.polyfit(x[good], np.log(cfj[good]), 3 if k.size >= 8 else 2)
    curve = np.exp(np.polyval(p, np.clip(x, x[good].min(), x[good].max())))
    return cfj, good, curve, None


MIN_TRACK = 5  # a reference curve needs the ridge at this many frequencies at least (as pick_max)


def alias_wavelength_km(cfg, geom):
    """Shortest wavelength the F-J spectrum of a line resolves: twice the mean station spacing along the line
    (aperture / (stations - 1)).  Spatial Nyquist limit of a sampled line; for irregular spacing the average
    sampling density is what counts (Landau 1967).  None for areal arrays."""
    if cfg["project"].get("array_type") != "line":
        return None
    g = geom or {}
    if not g.get("aperture_km") or int(g.get("stations", 0)) < 2:
        return None
    return 2.0 * float(g["aperture_km"]) / (int(g["stations"]) - 1)


def track_start(cfg, freqs, c, img, npair, lam_min_km=None):
    """reference.track_start_band_hz; auto: the tracking is tried from every octave (in steps of a third of an
    octave) and the start that lets the ridge be followed over the most frequencies is taken; between equally long
    tracks the slower one (at a given frequency the fundamental mode is the slowest mode)."""
    r = cfg["reference"]
    if r["track_start_band_hz"] != "auto":
        return list(r["track_start_band_hz"])
    best = None
    lo = float(freqs[npair >= MIN_FJ_PAIRS].min()) if (npair >= MIN_FJ_PAIRS).any() else float(freqs[0])
    while lo < freqs[-1]:
        try:
            pick, use, _, _ = pick_track(cfg, freqs, c, img, start_band=[lo, 2 * lo], need=1, lam_min_km=lam_min_km)
            key = (int(use.sum()), -float(np.mean(pick[use])))
            if best is None or key > best[0]:
                best = (key, [math.floor(lo * 1e4) / 1e4, math.ceil(2 * lo * 1e4) / 1e4])  # rounded outwards
        except RuntimeError:
            pass
        lo *= 2 ** (1 / 3)
    if best is None:
        raise RuntimeError(
            "F-J image: no ridge to start from at any frequency (work/04_reference/fj.png); borrow a "
            "reference (reference.source = file)"
        )
    return best[1]


def pick_track(
    cfg,
    freqs,
    c,
    img,
    min_peak=0.5,
    min_track=0.3,
    step=0.05,
    max_miss=3,
    rise=0.03,
    start_band=None,
    need=None,
    lam_min_km=None,
):
    s0, s1 = start_band or cfg["reference"]["track_start_band_hz"]
    Y = img / (img.max(axis=1, keepdims=True) + 1e-30)
    okc = (c >= 1.03 * c.min()) & (c <= 0.97 * c.max())  # peaks on the edges of the scan are not ridges
    v_lo, v_hi = [float(v) for v in cfg["data"]["velocity_range_kms"]]
    okc &= (c >= v_lo) & (c <= v_hi)  # the scan is wider than the range of the area only
    # to show what lies beside it
    normal = bool(cfg["reference"].get("normal_dispersion", True))
    ok_at = (lambda i: okc) if not lam_min_km else (lambda i: okc & (c >= lam_min_km * freqs[i]))  # wavelength limit
    start = None
    for strict in (True, False):
        for i in np.where((freqs >= s0) & (freqs <= s1))[0]:
            k = _peaks(Y[i], min_peak)
            k = k[ok_at(i)[k]]
            if k.size and np.argmax(Y[i][k]) == 0 and (k.size == 1 or not strict):
                start = i
                break
        if start is not None:
            break
    if start is None:
        raise RuntimeError(
            "F-J image: no unambiguous fundamental-mode start in %.2f-%.2f Hz "
            "(reference.track_start_band_hz); borrow a reference (reference.source = file)" % (s0, s1)
        )
    pick = np.full(len(freqs), np.nan)
    k = _peaks(Y[start], min_peak)
    pick[start] = c[k[ok_at(start)[k]][0]]
    for direction in (1, -1):
        last, miss, i = pick[start], 0, start + direction
        lowest = highest = last
        while 0 <= i < len(freqs) and miss < max_miss:
            k = _peaks(Y[i], min_track)
            k = k[ok_at(i)[k]]
            # the limit is on the slope d ln c / d ln f: `step` per frequency step, so it grows with the number of
            # steps since the last pick; within it the strongest peak (the ridge is the energy maximum)
            if not normal:
                mono = np.ones(k.size, bool)
            elif direction > 0:
                mono = c[k] <= lowest * (1 + rise)
            else:
                mono = c[k] >= highest * (1 - rise)
            k = k[mono & (np.abs(c[k] / last - 1) <= step * (miss + 1))]
            if k.size:
                j = k[np.argmax(Y[i][k])]
                pick[i] = last = c[j]
                miss = 0
                lowest, highest = min(lowest, last), max(highest, last)
            else:
                miss += 1
            i += direction
    use = np.isfinite(pick)
    if use.sum() < (MIN_TRACK if need is None else need):
        raise RuntimeError(
            "the F-J ridge could be followed over only %d frequencies from %.3g-%.3g Hz (at least %d "
            "are needed for a curve): look at work/04_reference/fj.png; say where the fundamental mode is "
            "clear (reference.track_start_band_hz), narrow data.velocity_range_kms to the velocities "
            "of the area, or borrow a curve (reference.source = file)" % (use.sum(), s0, s1, MIN_TRACK)
        )
    x = np.log(freqs)
    pf = np.polyfit(x[use], np.log(pick[use]), 3 if use.sum() >= 8 else min(2, use.sum() - 1))
    lo, hi = x[use].min(), x[use].max()
    raw = np.exp(np.polyval(pf, np.clip(x, lo, hi)))
    if not normal:
        return pick, use, raw, None
    curve = np.minimum.accumulate(raw)
    flat = np.abs(curve / raw - 1) > 0.01  # where forcing normal dispersion changed the fit by > 1 %
    note = None
    if flat.any():
        note = (
            "the fitted curve rises with frequency at %.3g-%.3g Hz and was held flat there (reference."
            "normal_dispersion); if the F-J ridge really rises there (a velocity reversal), set it to false"
            % (freqs[flat].min(), freqs[flat].max())
        )
    return pick, use, curve, note
