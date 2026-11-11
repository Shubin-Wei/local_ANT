"""Rules that derive parameters from the data (lant.auto)."""
import numpy as np

from lant import config


def test_parameters_from_velocity_and_reference():
    """lant.auto: the correlation parameters from the velocity range, the band and velocity axis from the reference."""
    from lant import auto
    cfg = config.defaults("surface")
    cfg["_user"] = {}
    cfg["data"]["velocity_range_kms"] = [1.5, 3.5]
    geom = dict(max_distance_km=14.0, p90_distance_km=9.5, median_distance_km=6.0)
    d = {k: v["value"] for k, v in auto.from_velocity(cfg, geom, [200.0]).items()}
    a = np.sqrt(np.log(2) / 5.0)
    f_lo = 2 ** (np.ceil(12 * np.log2(2 * 1.5 / 14.0) - 1e-9) / 12)          # longest pair = 2 wavelengths at 1.5 km/s
    assert d["correlation.fs"] == 50.0 and abs(d["dispersion.fmax_hz"] - 2 ** (46 / 12)) < 1e-6
    assert abs(d["correlation.fmin_hz"] - round(f_lo * (1 - a), 6)) < 1e-6
    assert d["correlation.maxlag_s"] >= 14.0 / (0.6 * 1.5) + 4 / f_lo
    assert d["correlation.seg_len_s"] >= 2 * d["correlation.maxlag_s"] and d["correlation.seg_step_s"] == d["correlation.seg_len_s"] / 2
    assert d["reference.fj_cmin"] <= 0.75 and d["reference.fj_cmax"] >= 5.25
    assert abs(d["export.tmin_s"] - 1 / d["dispersion.fmax_hz"]) < 1e-9
    cfg["dispersion"]["fmax_hz"] = d["dispersion.fmax_hz"]
    f = 2.0 ** (np.arange(-40, 46) / 12)
    c = np.exp(np.interp(np.log(f), np.log([0.1, 1, 10, 40]), np.log([3.4, 3.0, 2.0, 1.8])))
    u = c / (1 - np.gradient(np.log(c), np.log(f)))
    dec, warn = auto.from_reference(cfg, f, c, u, (0.5, 10.0), geom)
    r = {k: v["value"] for k, v in dec.items()}
    want = f[np.nonzero(f * 14.0 < 2.0 * c)[0].max() + 1]
    assert abs(r["dispersion.fmin_hz"] - want) < 1e-6
    band = (f >= want) & (f <= r["dispersion.fmin_hz"] * 1e9)
    assert r["dispersion.vmin"] <= 0.6 * u[band & (f <= d["dispersion.fmax_hz"])].min() and r["dispersion.vmax"] >= 1.3 * 3.0
    assert not warn
    cfg["data"]["velocity_range_kms"] = [0.2, 0.8]                         # measured velocities far outside it
    assert auto.from_reference(cfg, f, c, u, (0.5, 10.0), geom)[1]
