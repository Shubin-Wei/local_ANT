"""The UTC-locked 1-s periodic signal (the timing pulse of some loggers, Takagi et al. 2015): its detection in the
clean step (fold test of every station-day, summary per station in 01_clean/pulse_1s.csv) and its removal in the cc
step (correlation.remove_1s_comb).
"""
import csv
from collections import defaultdict

import numpy as np

PULSE_P = 1e-3            # significance level of the 1-s fold test of one station-day


def fold_1s(pieces, fs):
    """Test for a signal that repeats every UTC second (the timing pulse of some loggers, Takagi et al. 2015).
    pieces: [(index of the first sample counted from a whole UTC second, samples)] of gap-free record at a
    whole-number rate fs.  Every whole second has its mean removed and the seconds are averaged sample by sample
    (folding modulo 1 s).  For records without such a signal the variance of the folded waveform is the variance
    of the record / number of seconds; their ratio F follows an F distribution with (nu, nu (M - 1)) degrees of
    freedom (M seconds).  For white noise nu = p - 1 (p samples per second); neighbouring samples of real noise are
    correlated, which leaves fewer independent values: nu = (p - 1) / sum of the squared autocorrelation of the
    record inside a second (the equivalent degrees of freedom of a sum of squares of correlated values), otherwise
    coloured noise passes the test far more often than PULSE_P.  -> dict(ratio=F, power=share of the record's
    variance in the folded waveform, seconds=M, limit=F at the level PULSE_P, dof=nu), or None when it cannot be
    tested."""
    from scipy.stats import f as fdist
    p = int(round(fs))
    if p < 4 or abs(fs - p) > 1e-6:
        return None
    acc, sq, m = np.zeros(p), 0.0, 0
    spec = np.zeros(p // 2 + 1)
    for k0, d in pieces:
        off = (-int(k0)) % p
        n = (len(d) - off) // p
        if n < 1:
            continue
        b = np.asarray(d[off:off + n * p], np.float64).reshape(n, p)
        b = b - b.mean(axis=1, keepdims=True)
        acc += b.sum(axis=0)
        sq += float((b * b).sum())
        for i in range(0, n, 2000):                  # mean power spectrum of the seconds -> their autocorrelation
            spec += (np.abs(np.fft.rfft(b[i:i + 2000], axis=1)) ** 2).sum(axis=0)
        m += n
    if m < 60 or sq <= 0:
        return None
    tmpl = acc / m
    between = float((tmpl * tmpl).sum()) * m / (p - 1)             # variance between the phases x M
    within = (sq - float((tmpl * tmpl).sum()) * m) / ((m - 1) * (p - 1))
    if within <= 0:
        return None
    rho = np.fft.irfft(spec, p)
    nu = max(1.0, (p - 1) / float(((rho / rho[0]) ** 2).sum())) if rho[0] > 0 else float(p - 1)
    return dict(ratio=between / within, power=float((tmpl * tmpl).sum()) * m / sq, seconds=m, dof=round(nu, 1),
                limit=float(fdist.isf(PULSE_P, nu, nu * (m - 1))))


def pulse_summary(stats, path):
    """Stations with a 1-s periodic signal: the fold test (fold_1s) is passed in more than half of their days.
    Writes 01_clean/pulse_1s.csv (station, days tested, days with the signal, median F, median limit, median share
    of the variance)."""
    by = defaultdict(list)
    for s in stats:
        if "pulse_ratio" in s:
            by[s["station"]].append((s["pulse_ratio"], s["pulse_limit"], s["pulse_power"]))
    found, mx = [], 0.0
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["station", "days_tested", "days_detected", "median_F", "median_limit", "median_variance_share"])
        for sta in sorted(by):
            a = np.array(by[sta])
            nd = int((a[:, 0] > a[:, 1]).sum())
            w.writerow([sta, len(a), nd, "%.3f" % np.median(a[:, 0]), "%.3f" % np.median(a[:, 1]), "%.3g" % np.median(a[:, 2])])
            if nd > 0.5 * len(a):
                found.append(sta)
                mx = max(mx, float(np.median(a[:, 2])))
    return dict(tested=len(by), stations=found, max_power=mx)


def remove_comb(y, ok, fs):
    """Subtract the UTC-locked 1-s periodic component hour by hour (hours with < 50 % valid samples untouched)."""
    p = int(round(fs))
    ph_all = np.arange(y.size) % p
    for h in range(24):
        a, b = int(h * 3600 * fs), int((h + 1) * 3600 * fs)
        m = ok[a:b]
        if m.mean() < 0.5:
            continue
        seg = y[a:b].astype(np.float64)
        ph = ph_all[a:b]
        fold = np.bincount(ph[m], weights=seg[m], minlength=p)
        cnt = np.bincount(ph[m], minlength=p)
        tmpl = fold / np.maximum(cnt, 1)
        tmpl -= tmpl.mean()
        y[a:b] = (seg - tmpl[ph]).astype(np.float32)
    return y
