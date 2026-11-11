"""Signal and geometry helpers (lant.signal)."""

import numpy as np

from lant import signal


def test_signal_helpers():
    """Common sampling rates, octave bands, logarithmic frequencies, runs, distances on the ellipsoid and sphere."""
    assert signal.common_rates([500.0, 200.0])[:3] == [100.0, 50.0, 25.0]
    assert signal.octave_bands(signal.log_freqs(0.5, 10.08))[0] == (0.5, 1.0)
    f = signal.log_freqs(0.5, 10.5, 12)
    assert len(f[f <= 10.5]) == 53 and abs(f[12] - 1.0) < 1e-12
    m = np.array([1, 1, 0, 1, 1, 1, 1, 1, 1, 0], bool)
    assert signal.contiguous(m, 6).tolist() == [0, 0, 0, 1, 1, 1, 1, 1, 1, 0]
    d, az = signal.geodesic_km(100.0, 30.0, 100.0, 31.0)
    assert abs(d - 110.861) < 0.001 and abs(az) < 1e-6  # WGS84 (111.195 on a sphere of 6371 km)
    assert abs(signal.sphere_km(100.0, 30.0, 100.0, 31.0) - 111.195) < 0.001


def test_gradient_in_runs():
    """The slope used by the smoothness and group-phase checks exists at the ends of a run (one-sided), not nan."""
    from lant.signal import gradient_in_runs

    y = np.array([np.nan, 1.0, 2.0, 4.0, np.nan, 5.0, np.nan])
    g = gradient_in_runs(y, np.arange(7.0))
    assert np.allclose(g[1:4], [1.0, 1.5, 2.0]) and np.isnan(g[[0, 4, 5, 6]]).all()
