"""The measurement against known answers.

In a laterally uniform medium with phase velocity c(f) and isotropic noise, the cross-correlation spectrum of two
stations r apart is J0(2 pi f r / c(f)) (Aki 1957): the two-station measurement and the F-J spectrum must give c(f)
back.  These are exact correlations, not data: they check that the formulas are implemented right.
"""

import numpy as np
import pytest
from scipy.special import j0

from lant import dispersion as D
from lant import fj, reference
from lant.signal import branches, log_freqs

C_TRUE = lambda f: 2.2 + 0.9 / (1 + (np.asarray(f) / 2.0) ** 1.5)  # km/s, normal dispersion
FS, MAXLAG = 50.0, 14.0


def j0_correlations(dist, band=(0.4, 16.0)):
    L = int(MAXLAG * FS)
    nfft = 2 * L
    f = np.fft.rfftfreq(nfft, 1 / FS)
    m = (f >= band[0]) & (f <= band[1])
    tap = np.ones(m.sum())
    tap[:20] = 0.5 - 0.5 * np.cos(np.pi * np.arange(20) / 20)
    tap[-20:] = tap[:20][::-1]
    out = np.zeros((len(dist), 2 * L + 1))
    for q, d in enumerate(dist):
        S = np.zeros(len(f))
        S[m] = j0(2 * np.pi * f[m] * d / C_TRUE(f[m])) * tap
        x = np.fft.irfft(S, nfft)
        out[q] = np.r_[x[-L:], x[: L + 1]]
    return out / np.abs(out).max(axis=1, keepdims=True)


def settings(**kw):
    p = dict(alpha=5.0, group_window=0.3, snr_min=5.0, min_wavelengths=2.0, reference_error=0.1,
             max_jackknife_rel=None, slip_cycles=0.3, min_points=6)  # fmt: skip
    p.update(kw)
    return p


def measure(dist, ccf, freqs, c_ref, p, subsets=None):
    out = []
    for q in range(len(dist)):
        pos, neg = branches(ccf[q : q + 1], FS, MAXLAG)
        sym = 0.5 * (pos + neg)[0]
        out.append(D.measure_pair(sym, np.zeros((0, sym.size)) if subsets is None else subsets[q], dist[q], FS, freqs,
                                  c_ref, reference.group_velocity(freqs, c_ref), p))  # fmt: skip
    return out


def true_group_velocity(f):
    """U = c / (1 - (f / c) dc/df) of C_TRUE."""
    ff = np.exp(np.linspace(np.log(0.3), np.log(30.0), 20000))
    c = C_TRUE(ff)
    return np.interp(f, ff, c / (1 - ff / c * np.gradient(c, ff)))


def test_group_velocity_of_exact_correlations():
    """Envelope maximum at the instantaneous frequency (FTAN): within 1 % (95 %) of the exact group velocity; the
    median bias of the wide Gaussian filter (alpha = 5) is below 0.7 %."""
    dist = np.linspace(1.5, 11.0, 40)
    freqs = log_freqs(1.0, 10.0, 12)
    res = measure(dist, j0_correlations(dist), freqs, C_TRUE(freqs), settings())
    u = np.array([r["u"] for r in res])
    ok = np.array([D.accepted(r["reason_u"], r["u"]) for r in res])
    err = u / true_group_velocity(freqs)[None, :] - 1
    assert ok.sum() > 0.5 * ok.size
    assert abs(np.median(err[ok])) < 0.007
    assert np.percentile(np.abs(err[ok]), 95) < 0.01


def test_group_times_follow_the_instantaneous_frequency():
    """Values at the instantaneous frequencies are moved back onto the grid, never across a gap."""
    freqs = np.array([1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
    f_inst = np.array([1.02, 1.12, 1.22, np.nan, 1.42, 1.52])
    tg = np.array([2.0, 2.1, 2.2, np.nan, 2.4, 2.5])
    out = D.group_times(freqs, f_inst, tg)
    assert np.isnan(out[0]) and np.isnan(out[3]) and np.isnan(out[4])  # below the first run, in the gap
    assert np.allclose(out[1:3], np.interp(np.log(freqs[1:3]), np.log(f_inst[:3]), tg[:3]))
    assert np.isclose(out[5], np.interp(np.log(1.5), np.log(f_inst[4:]), tg[4:]))


def test_phase_velocity_of_exact_correlations():
    """Far field (>= 5 wavelengths): a few tenths of a percent; 3-5 wavelengths: within 1.5 % (the far-field phase
    pi/4 is an approximation of J0 there)."""
    dist = np.linspace(1.5, 11.0, 40)
    freqs = log_freqs(1.0, 10.0, 12)
    res = measure(dist, j0_correlations(dist), freqs, C_TRUE(freqs), settings())
    c = np.array([r["c"] for r in res])
    ok = np.array([D.accepted(r["reason"], r["c"]) for r in res])
    lam = dist[:, None] * freqs[None, :] / C_TRUE(freqs)[None, :]
    err = c / C_TRUE(freqs)[None, :] - 1
    assert ok.sum() > 0.3 * ok.size
    far, mid = ok & (lam >= 5), ok & (lam >= 3) & (lam < 5)
    assert abs(np.median(err[ok])) < 0.003
    assert np.percentile(np.abs(err[far]), 95) < 0.004
    assert np.percentile(np.abs(err[mid]), 95) < 0.015


def test_a_wrong_reference_within_the_assumed_error_gives_the_same_cycle():
    """The cycle count is set where the path is at most 1 / (2 e) wavelengths: a reference 8 % off (e = 0.1) gives
    the same phase velocities."""
    dist = np.array([2.0, 4.0, 6.0])
    freqs = log_freqs(1.0, 8.0, 12)
    ccf = j0_correlations(dist)
    a = measure(dist, ccf, freqs, C_TRUE(freqs), settings())
    b = measure(dist, ccf, freqs, 1.08 * C_TRUE(freqs), settings(group_window=0.4))
    for ra, rb in zip(a, b):
        both = D.accepted(ra["reason"], ra["c"]) & D.accepted(rb["reason"], rb["c"])
        assert both.sum() >= 6 and np.allclose(ra["c"][both], rb["c"][both], rtol=2e-3)


def test_path_too_long_for_the_cycle_count_is_not_measured():
    dist = np.array([12.0])
    freqs = log_freqs(4.0, 10.0, 12)  # 12 km at >= 4 Hz is > 5 wavelengths everywhere
    r = measure(dist, j0_correlations(dist), freqs, C_TRUE(freqs), settings())[0]
    assert not D.accepted(r["reason"], r["c"]).any()
    assert ((r["reason"] & D.REASONS["cycle_ambiguous"]) > 0).any()


def test_jackknife_of_identical_substacks_is_zero_and_rejection_works():
    dist = np.array([4.0])
    freqs = log_freqs(1.0, 8.0, 12)
    ccf = j0_correlations(dist)
    pos, neg = branches(ccf, FS, MAXLAG)
    sym = 0.5 * (pos + neg)
    r = measure(dist, ccf, freqs, C_TRUE(freqs), settings(max_jackknife_rel=0.02), [np.repeat(sym, 4, axis=0)])[0]
    ok = D.accepted(r["reason"], r["c"])
    assert ok.sum() >= 6 and np.nanmax(r["sigma"][ok]) < 1e-9
    noisy = sym + 0.3 * np.random.default_rng(1).standard_normal((4, sym.shape[1]))
    r2 = measure(dist, ccf, freqs, C_TRUE(freqs), settings(max_jackknife_rel=1e-6), [noisy])[0]
    assert ((r2["reason"] & D.REASONS["jackknife"]) > 0).any()


def test_one_curve_per_pair():
    dist = np.array([5.0])
    freqs = log_freqs(1.0, 10.0, 12)
    r = measure(dist, j0_correlations(dist), freqs, C_TRUE(freqs), settings())[0]
    ok = D.accepted(r["reason"], r["c"])
    k = np.where(ok)[0]
    assert k.size and np.all(np.diff(k) == 1)  # one consecutive run


def test_fj_spectrum_peaks_on_the_curve():
    rng = np.random.default_rng(31)
    x, y = rng.uniform(-4, 4, 20), rng.uniform(-4, 4, 20)
    iu, ju = np.triu_indices(20, 1)
    dist = np.hypot(x[iu] - x[ju], y[iu] - y[ju])
    freqs = log_freqs(1.0, 10.0, 12)
    c = np.arange(1.0, 4.5, 0.005)
    img = fj.image(j0_correlations(dist), dist, FS, freqs, c)
    picks = fj.pick(img, c, freqs)
    ok = np.isfinite(picks)
    assert ok.sum() >= 20 and np.abs(picks[ok] / C_TRUE(freqs[ok]) - 1).max() < 0.015
    r = np.linspace(0.05, 10, 400)
    k = np.linspace(0.5, 6, 1101)
    assert abs(k[np.argmax(fj.fj_integral(j0(3.0 * r), r, k))] - 3.0) < 0.02  # a single J0 peaks at its k


def test_group_velocity_of_the_reference():
    f = log_freqs(0.5, 8.0, 12)
    assert np.allclose(reference.group_velocity(f, np.full(len(f), 2.0)), 2.0)
    c = C_TRUE(f)
    u = reference.group_velocity(f, c)
    assert np.all(u[1:-1] < c[1:-1])  # normal dispersion: group slower than phase


def test_matched_phase_reads_a_delayed_pulse():
    """A pulse delayed by t0 has phase travel time t0 (+ the pi/4 convention) at every frequency."""
    fs, t0 = 50.0, 3.0
    t = np.arange(int(14 * fs)) / fs
    f = log_freqs(1.0, 8.0, 12)
    sig = np.zeros(t.size)
    sig[int(t0 * fs)] = 1.0
    tp = D.matched_phase(sig[None, :], fs, f, np.full(len(f), t0), 5.0, phi0=0.0)[0]
    T = 1 / f
    assert np.allclose((tp - t0 + T / 2) % T - T / 2, 0, atol=1e-6)


@pytest.mark.parametrize("bad", [0.0, -1.0])
def test_smooth_curve_needs_picks(bad):
    with pytest.raises(RuntimeError):
        fj.smooth_curve(np.array([1.0, 2.0, 3.0]), np.array([np.nan, np.nan, 2.0 + bad]))


def test_delete_one_stacks_of_a_pair_with_one_sub_stack():
    """A pair whose data all come from one sub-stack has no delete-one stack (and so no jackknife error)."""
    nlag = 2 * int(MAXLAG * FS) + 1
    x = np.random.default_rng(0).standard_normal(nlag)
    subs = np.array([x, np.zeros(nlag), np.zeros(nlag)])
    sets = D.delete_one_stacks(10 * x, 10.0, subs, np.array([10.0, 0.0, 0.0]), FS, MAXLAG, int(MAXLAG * FS) + 1)
    assert sets.shape == (0, int(MAXLAG * FS) + 1)
    sets = D.delete_one_stacks(10 * x + 5 * x, 15.0, np.array([x, x, np.zeros(nlag)]), np.array([10.0, 5.0, 0.0]), FS,
                               MAXLAG, int(MAXLAG * FS) + 1)  # fmt: skip
    assert sets.shape[0] == 2 and np.allclose(sets[0], sets[1])
