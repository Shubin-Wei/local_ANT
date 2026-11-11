"""Numerical benchmarks of the measurement core against known answers (no external program needed).

The correlations of a laterally uniform medium with phase velocity c(f) are J0(2 pi f d / c(f)); the two-station
measurement and the F-J spectrum must give c(f) back.  The cross-coherence of a record with its delayed copy must
peak at the positive delay.
"""

import h5py
import numpy as np
from scipy.special import j0

from lant import config, signal

C_TRUE = lambda f: 2.2 + 0.9 / (1 + (np.asarray(f) / 2.0) ** 1.5)  # km/s, normal dispersion


def _stack(path, nsta=14, fs=50.0, maxlag=14.0, band=(0.4, 16.0), nblk=4):
    rng = np.random.default_rng(31)
    x, y = rng.uniform(-4, 4, nsta), rng.uniform(-4, 4, nsta)
    iu, ju = np.triu_indices(nsta, 1)
    dist = np.hypot(x[iu] - x[ju], y[iu] - y[ju])
    L = int(maxlag * fs)
    nfft = 2 * L
    f = np.fft.rfftfreq(nfft, 1 / fs)
    m = (f >= band[0]) & (f <= band[1])
    tap = np.ones(m.sum())
    tap[:20] = 0.5 - 0.5 * np.cos(np.pi * np.arange(20) / 20)
    tap[-20:] = tap[:20][::-1]
    st = np.zeros((len(dist), 2 * L + 1), np.float32)
    for q, d in enumerate(dist):
        S = np.zeros(len(f))
        S[m] = j0(2 * np.pi * f[m] * d / C_TRUE(f[m])) * tap
        xx = np.fft.irfft(S, nfft)
        st[q] = np.r_[xx[-L:], xx[: L + 1]]
    st /= np.abs(st).max(axis=1, keepdims=True)
    daily = np.stack([st + 0.002 * rng.standard_normal(st.shape).astype(np.float32) for _ in range(nblk)])
    with h5py.File(path, "w") as h:
        h["stack"] = daily.mean(0)
        h["daily"] = daily
        h["daily_count"] = np.full((nblk, len(dist)), 100, np.int32)
        h["nseg"] = np.full(len(dist), 100 * nblk)
        h["pairs"] = np.array([("S%02d_S%02d" % (i, j)).encode() for i, j in zip(iu, ju)])
        h["days"] = np.array([("202400%d" % (b + 1)).encode() for b in range(nblk)])
        h["dist_km"] = dist
        h.attrs.update(fs=fs, maxlag=maxlag, fmin=band[0], fmax=band[1])
    return dist


def _cfg():
    cfg = config.defaults()
    cfg["data"]["velocity_range_kms"] = [1.5, 3.5]
    cfg["correlation"].update(fmin_hz=0.4, fmax_hz=16.0)
    cfg["dispersion"].update(fmin_hz=1.0, fmax_hz=10.0, vmin=1.2, vmax=4.2, dv=0.01, njack=4, nproc=1)
    cfg["reference"].update(fj_cmin=1.0, fj_cmax=4.5, fj_dc=0.005, fj_c_nominal=2.6)
    return cfg


def test_two_station_phase_velocity(tmp_path):
    """Phase velocities of J0 correlations against c(f).  (Between 2 and 3 wavelengths the far-field phase pi/4
    is itself in error by up to about 3 %: the reason for the wavelength limit qc.lam_c.)"""
    from lant import disp, reference, qc

    stack = str(tmp_path / "stack.h5")
    dist = _stack(stack)
    cfg = _cfg()
    f = reference.freq_grid(cfg)
    ref = str(tmp_path / "reference.npz")
    np.savez(ref, freqs=f, c_ref=C_TRUE(f), ref_valid=np.ones(len(f), bool), ref_support=np.ones(len(f), bool))
    out = str(tmp_path / "disp.h5")
    s = disp.measure(cfg, stack, ref, out)
    qc.lambda_qc(cfg, out)
    with h5py.File(out, "r") as h:
        ok, vc = h["ok_c"][:], h["vc"][:]
    err = vc / C_TRUE(f)[None, :] - 1
    lam = dist[:, None] * f[None, :] / C_TRUE(f)[None, :]  # path length in wavelengths
    assert ok.sum() > 0.25 * s["pairs"] * len(f)
    far, mid = ok & (lam >= 5), ok & (lam >= 3) & (lam < 5)
    assert abs(np.median(err[ok])) < 0.003
    assert np.percentile(np.abs(err[far]), 95) < 0.004  # far field: a few tenths of a percent
    assert np.percentile(np.abs(err[mid]), 95) < 0.015  # 3-5 wavelengths: the far-field phase pi/4
    per_f = np.array([np.median(err[far[:, k], k]) for k in range(len(f)) if far[:, k].sum() >= 10])
    assert np.abs(per_f).max() < 0.004  # no frequency-dependent bias


def test_fj_spectrum_picks_the_curve(tmp_path):
    """The F-J image of J0 correlations peaks on c(f): picks within 1.5 %."""
    from lant import fj

    stack = str(tmp_path / "stack.h5")
    _stack(stack, nsta=20)
    cfg = _cfg()
    f = signal.log_freqs(1.0, 10.0, 12)
    freqs, c, img, npair = fj.fj_image(cfg, stack, f)
    picks, used, curve, _ = fj.pick_max(cfg, freqs, c, img, npair)
    assert used.sum() >= 20
    assert np.abs(picks[used] / C_TRUE(freqs[used]) - 1).max() < 0.015
    # the integral itself: a single J0 of wavenumber k0 peaks at k0
    r = np.linspace(0.05, 10, 400)
    k = np.linspace(0.5, 6, 1101)
    assert abs(k[np.argmax(fj.fj_spectrum(j0(3.0 * r), r, k))] - 3.0) < 0.02


def test_correlation_sign(monkeypatch):
    """Station B records what station A recorded 0.4 s earlier: the cross-coherence A-B peaks at +0.4 s."""
    from lant import cc, data as D

    fs, n = 50.0, int(86400 * 50)
    rng = np.random.default_rng(32)
    a = rng.standard_normal(n).astype(np.float32)
    recs = {"XX.A": a, "XX.B": np.roll(a, 20)}
    monkeypatch.setattr(D, "read_day", lambda rows, t0, fs_: (recs[rows].copy(), np.ones(n, bool)))
    opt = dict(
        fs=fs,
        seg_len=40.0,
        seg_step=20.0,
        maxlag=5.0,
        fmin=0.5,
        fmax=20.0,
        smooth_bins=1,
        max_over_std=10.0,
        decomb=False,
        pair_chunk=100,
    )
    names, hourly, count, _, _ = cc.day_ccf({"XX.A": "XX.A", "XX.B": "XX.B"}, 0.0, opt, nbin=24)
    x = hourly[:2].mean(0)[0]
    lag = (np.arange(len(x)) - len(x) // 2) / fs
    assert names == ["XX.A", "XX.B"] and abs(lag[np.argmax(x)] - 0.4) < 1e-9 and count[0, 0] > 100


def test_distance_on_the_ellipsoid():
    """geodesic_km against ObsPy (WGS84) at several latitudes and azimuths; the sphere differs by up to 0.6 %."""
    from obspy.geodetics import gps2dist_azimuth
    from lant.signal import geodesic_km, sphere_km, local_xy_km

    rng = np.random.default_rng(1)
    la1, lo1 = rng.uniform(-70, 70, 200), rng.uniform(-180, 180, 200)
    la2, lo2 = la1 + rng.uniform(-0.2, 0.2, 200), lo1 + rng.uniform(-0.2, 0.2, 200)
    d, az = geodesic_km(lo1, la1, lo2, la2)
    ref = np.array([gps2dist_azimuth(a, b, c, e)[:2] for a, b, c, e in zip(la1, lo1, la2, lo2)])
    assert np.abs(d / (ref[:, 0] / 1e3) - 1).max() < 1e-8
    assert np.abs((az - ref[:, 1] + 180) % 360 - 180).max() < 1e-6
    assert 0.001 < np.abs(sphere_km(lo1, la1, lo2, la2) / d - 1).max() < 0.007
    x, y = local_xy_km(np.r_[100.0, 100.03], np.r_[30.0, 30.02], 30.01, 100.015)
    assert abs(np.hypot(x[1] - x[0], y[1] - y[0]) / geodesic_km(100.0, 30.0, 100.03, 30.02)[0] - 1) < 1e-6
