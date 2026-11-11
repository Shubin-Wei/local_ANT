"""The reference step and the F-J ridge (lant/reference.py, lant/fj.py): ridge tracking, the normal-dispersion switch,
borrowed curves."""

import os

import numpy as np

from lant import config
from util import tmp_dir


def test_track_steep_ridge_and_alias_limit():
    """track_fundamental (1) takes the strongest peak inside the window, not a weak one closer to the last pick;
    (2) bridges a single drop of 8 % in one step (the limit is on the slope, so it widens with the steps since the
    last pick); (3) on a line keeps only wavelengths of at least twice the mean station spacing."""
    from lant import fj as FJ

    cfg = config.defaults("line")
    cfg["data"]["velocity_range_kms"] = [0.5, 4.0]
    cfg["reference"]["track_start_band_hz"] = [1.0, 1.1]
    freqs = 2.0 ** (np.arange(0, 30) / 12)  # 1-5.3 Hz, 1/12 octave
    c = np.linspace(0.3, 5.0, 2000)
    gauss = lambda v, a=1.0: a * np.exp(-(((c[None, :] - v[:, None]) / 0.01) ** 2))
    # (1) ridge -4 % per step; a weak peak where the ridge was one step before (0 % change: nearer)
    ridge = 2.5 * 0.96 ** np.arange(30)
    img = gauss(ridge) + gauss(np.r_[ridge[0], ridge[:-1]], 0.4)
    pick, use, _, _ = FJ.pick_track(cfg, freqs, c, img)
    assert use.sum() == 30 and np.max(np.abs(pick / ridge - 1)) < 0.005
    # (2) ridge -2 % per step with one step of -8 %
    ridge = 2.5 * 0.98 ** np.arange(30) * np.where(np.arange(30) >= 10, 0.94, 1.0)
    pick, use, _, _ = FJ.pick_track(cfg, freqs, c, gauss(ridge))
    assert use.sum() >= 29 and freqs[use].max() > 5 and np.nanmax(np.abs(pick / ridge - 1)) < 0.005
    # (3) the spatial Nyquist limit of a line
    geom = dict(aperture_km=10.0, stations=21)  # mean spacing 0.5 km -> wavelength >= 1 km
    lam = FJ.alias_wavelength_km(cfg, geom)
    assert abs(lam - 1.0) < 1e-12
    pick, use, _, _ = FJ.pick_track(cfg, freqs, c, gauss(ridge), lam_min_km=lam)
    assert (
        np.all(pick[use] >= lam * freqs[use]) and use.sum() >= 5 and freqs[use].max() < 1.9
    )  # at 1.89 Hz the ridge (1.88 km/s) is below lam x f
    assert FJ.alias_wavelength_km(dict(cfg, project=dict(cfg["project"], array_type="surface")), geom) is None


def test_normal_dispersion_switch():
    """track_fundamental: a rising ridge is held flat (with a note) unless reference.normal_dispersion is false."""
    from lant import fj as FJ

    cfg = config.defaults("surface")
    freqs = 2.0 ** (np.arange(0, 25) / 12)
    c = np.linspace(0.5, 3.0, 300)
    ridge = 1.0 + 0.2 * np.log2(freqs)  # rises with frequency
    img = np.exp(-(((c[None, :] - ridge[:, None]) / 0.03) ** 2))
    cfg["reference"]["track_start_band_hz"] = [1.0, 1.5]
    _, use, curve, note = FJ.pick_track(cfg, freqs, c, img, step=0.1, rise=0.5)
    assert note and np.all(np.diff(curve) <= 1e-12)
    cfg["reference"]["normal_dispersion"] = False
    _, use, curve, note = FJ.pick_track(cfg, freqs, c, img)
    assert note is None and use.sum() > 20 and np.max(np.abs(curve / ridge - 1)) < 0.02


def test_borrowed_curve():
    """A reference curve borrowed from another project: valid only where the source is, and a held (unsupported)
    stretch of the source stays unsupported."""
    from lant.reference import borrowed_curve

    d = tmp_dir()
    np.savez(
        os.path.join(d, "r.npz"),
        freqs=np.array([1.0, 2, 3, 4, 5]),
        c_ref=np.array([3.0, 2.9, 2.8, 2.7, 2.6]),
        ref_valid=np.array([True, True, False, True, True]),
    )
    cb, valid, sup, origin, _ = borrowed_curve(os.path.join(d, "r.npz"), np.array([1.5, 2.5, 3.0, 4.5, 6.0]))
    assert list(valid) == [True, False, False, True, False]
    np.savez(
        os.path.join(d, "s.npz"),
        freqs=np.array([1.0, 2, 3, 4]),
        c_ref=np.array([3.0, 2.9, 2.8, 2.8]),
        ref_valid=np.ones(4, bool),
        ref_support=np.array([True, True, True, False]),
    )
    cb, valid, sup, origin, _ = borrowed_curve(os.path.join(d, "s.npz"), np.array([1.5, 3.0, 3.5]))
    assert list(valid) == [True, True, True] and list(sup) == [True, True, False]
