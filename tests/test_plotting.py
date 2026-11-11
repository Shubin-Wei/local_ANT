"""Choices the figures make from the data: band edges and bands of the moveout figures, their distance and velocity
window, the example pairs, the relief tiles; and the checkerboard noise from the fit of the main model."""
import os

import h5py
import numpy as np

from lant import config
from lant.plotting import ccf, dispersion, relief


def test_band_edges_and_bands():
    """Band edges at whole Hz from 2 Hz, 0.5 Hz from 0.5 Hz, 0.1 Hz below; the bands join, are at least 2/3 octave
    wide, at most four, and cover the rounded measurement band."""
    assert [ccf.nice_hz(f) for f in (13.6, 2.4, 1.26, 0.63, 0.31, 0.04)] == [14.0, 2.0, 1.5, 0.5, 0.3, 0.1]
    for lo, hi in ((0.63, 14.3), (1.0, 1.3), (2.2, 47.0), (0.12, 3.1)):
        b = ccf.figure_bands(lo, hi)
        assert b[0][0] == ccf.nice_hz(lo) and b[-1][1] == ccf.nice_hz(hi) and 1 <= len(b) <= 4
        assert all(b[k][1] == b[k + 1][0] for k in range(len(b) - 1))
        assert len(b) == 1 or all(np.log2(h / l) >= 2 / 3 - 1e-9 for l, h in b)
    assert ccf.figure_bands(0.63, 14.3) == [(0.5, 1.0), (1.0, 2.0), (2.0, 5.0), (5.0, 14.0)]


def _disp(path, n=40):
    rng = np.random.default_rng(1)
    freqs = np.geomspace(1.0, 10.0, 30)
    dist = np.sort(rng.uniform(0.5, 10.0, n))
    ok = np.zeros((n, 30), bool)
    for j, d in enumerate(dist):
        ok[j, (freqs * d / 2.0 >= 1.0) & (np.arange(30) < 10 + j % 15)] = True       # at least one wavelength
    vc = np.tile(2.0 - 0.05 * np.log(freqs), (n, 1)).astype(np.float32)
    with h5py.File(path, "w") as h:
        h["freqs"], h["dist_km"], h["ok_c"], h["vc"], h["vg"] = freqs, dist, ok, vc, (0.8 * vc).astype(np.float32)
        h["pairs"] = np.array([("P%02d_Q%02d" % (j, j)).encode() for j in range(n)])
    return dist, ok, freqs


def test_moveout_window_and_example_pairs(tmp_path):
    """The moveout window of a band from the measurements in it (distance, slow and fast velocity on 0.1 km/s); the
    example pairs at the 20 / 50 / 80 % distance quantiles, ordered by distance, with rounded bands."""
    p = str(tmp_path / "disp.h5")
    dist, ok, freqs = _disp(p)
    w = ccf.band_window(p, (2.0, 5.0))
    m = (freqs >= 2.0) & (freqs <= 5.0)
    assert abs(w["dmax"] - np.percentile(dist[ok[:, m].any(1)], 99)) < 1e-9
    assert w["vslow"] <= 0.8 * 1.95 and w["vfast"] >= 2.0 - 0.05 * np.log(2.0) - 1e-6
    assert ccf.band_window(p, (20.0, 30.0)) is None                  # nothing measured in the band
    pairs = dispersion.choose_pairs(p)
    d = {("P%02d_Q%02d" % (j, j)): dist[j] for j in range(len(dist))}
    assert len(pairs) == 3 and [d[n] for n, _ in pairs] == sorted(d[n] for n, _ in pairs)
    assert all(ccf.nice_hz(lo) == lo and ccf.nice_hz(hi) == hi and hi > lo for _, (lo, hi) in pairs)


def test_relief_tiles_and_limits(tmp_path):
    """SRTM tile names on both sides of the equator and of Greenwich; no automatic relief for wide or polar maps,
    none when switched off, a missing grid file named."""
    assert relief.tile_names(99.5, 100.5, 29.5, 30.2) == ["N29E099", "N29E100", "N30E099", "N30E100"]
    assert relief.tile_names(-0.5, 0.5, -0.5, 0.5) == ["S01W001", "S01E000", "N00W001", "N00E000"]
    assert relief.tile_names(100.0, 101.0, 30.0, 31.0) == ["N30E100"]                # edges on whole degrees
    cfg = config.defaults()
    cfg["plots"]["topography"] = "auto"
    assert relief.grid(cfg, 100, 104, 30, 31, str(tmp_path))[0] is None
    assert "60 N" in relief.grid(cfg, 10, 11, 65, 66, str(tmp_path))[1]
    cfg["plots"]["topography"] = None
    assert relief.grid(cfg, 100, 101, 30, 31, str(tmp_path)) == (None, "plots.topography is none")
    cfg["plots"]["topography"] = str(tmp_path / "none.nc")
    assert "not found" in relief.grid(cfg, 100, 101, 30, 31, str(tmp_path))[1]


def test_checkerboard_noise_from_the_main_fit(tmp_path, monkeypatch):
    """inversion.noise auto: the relative rms of the residuals that DSurfTomo's rejection rule keeps; computed again
    only when the main model changes."""
    from lant import pipeline
    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path)
    main = tmp_path / "06_inversion" / "main"
    main.mkdir(parents=True)
    (main / "DSurfTomo.inMeasure.dat").write_text("model")
    obs = np.array([1.0, 2.0, 1.0, 4.0])
    pred = np.array([0.98, 2.04, 1.5, 4.0])
    R = np.c_[np.zeros(4), pred, obs, np.zeros(4), np.zeros(4), [1, 1, 0, 1]]      # the third datum rejected
    calls = []
    monkeypatch.setattr(pipeline.RUN, "evaluate", lambda cfg, d: calls.append(d) or (None, R, None))
    nz = pipeline.checkerboard_noise(cfg, str(main))
    assert abs(nz["value"] - np.sqrt((0.02 ** 2 + 0.02 ** 2 + 0) / 3)) < 1e-12 and nz["n_kept"] == 3 and nz["n"] == 4
    assert pipeline.checkerboard_noise(cfg, str(main))["value"] == nz["value"] and len(calls) == 1     # kept
    os.utime(main / "DSurfTomo.inMeasure.dat", (1e9, 1e9))
    pipeline.checkerboard_noise(cfg, str(main))
    assert len(calls) == 2                                            # a new main model: computed again
