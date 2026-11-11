"""The stack step: excluded station-days and the station diagnosis."""
import os

import numpy as np

from lant import config


def test_block_exclusion(tmp_path):
    """An entry for a sub-day block (YYYYDDD_hh) removes that block; a whole-day entry removes all its blocks."""
    import h5py
    from lant.stack import exclude_station_days
    p = str(tmp_path / "s.h5")
    with h5py.File(p, "w") as h:
        h["stack"] = np.ones((1, 5), np.float32)
        h["daily"] = np.ones((2, 1, 5), np.float32)
        h["daily_count"] = np.full((2, 1), 10, np.int32)
        h["nseg"] = np.array([20])
        h["pairs"] = np.array([b"A_B"])
        h["days"] = np.array([b"2024001_00", b"2024001_12"])
        for k in ("az_deg", "dist_km", "lat1", "lat2", "lon1", "lon2"):
            h[k] = np.zeros(1)
    assert exclude_station_days(p, [("A", "2024001_00")], str(tmp_path / "a.h5"))["dropped_segments"] == 10
    assert exclude_station_days(p, [("A", "2024001")], str(tmp_path / "b.h5"))["dropped_segments"] == 20


def test_station_diagnosis(tmp_path):
    """A station whose records are 0.03 s late in one sub-stack is the only candidate; the thresholds come from the
    distribution; with too few station sub-stacks nothing is judged."""
    import h5py
    from lant import stack as SK
    rng = np.random.default_rng(21)
    fs, nl, nblk, names = 50.0, 401, 12, ["S%02d" % i for i in range(8)]
    t = (np.arange(nl) - nl // 2) / fs
    pairs, dist = [], []
    for i in range(8):
        for j in range(i + 1, 8):
            pairs.append("%s_%s" % (names[i], names[j]))
            dist.append(0.5 * (j - i))
    wav = lambda t0, d: np.exp(-((np.abs(t - t0) - d / 2.0) / 0.08) ** 2) * np.cos(2 * np.pi * 5 * (np.abs(t - t0) - d / 2.0))
    daily = np.zeros((nblk, len(pairs), nl), np.float32)
    for k, (p, d) in enumerate(zip(pairs, dist)):
        for b in range(nblk):
            shift = 0.0
            if b == 4 and "S03" in p.split("_"):
                shift = 0.03 if p.startswith("S03") else -0.03
            daily[b, k] = wav(shift, d) + 0.05 * rng.standard_normal(nl)
    f = str(tmp_path / "s.h5")
    with h5py.File(f, "w") as h:
        h["daily"] = daily
        h["stack"] = daily.mean(0)
        h["daily_count"] = np.full((nblk, len(pairs)), 10, np.int32)
        h["nseg"] = np.full(len(pairs), 10 * nblk)
        h["pairs"] = np.array([p.encode() for p in pairs])
        h["days"] = np.array([("2024001_%02d" % (2 * b)).encode() for b in range(nblk)])
        h["dist_km"] = np.array(dist)
        h.attrs.update(fs=fs, fmin=2.0, fmax=10.0)
    cfg = config.defaults()
    cfg["data"]["velocity_range_kms"] = [1.0, 3.0]
    cfg["correlation"].update(fmin_hz=2.0, fmax_hz=10.0)
    flags, s = SK.station_qc(cfg, f, str(tmp_path / "q"))
    assert list(flags) == [("S03", "2024001_08")] and "shift" in flags[("S03", "2024001_08")]
    z = np.load(tmp_path / "q" / "station_qc.npz")
    assert abs(abs(z["shift"][3, 4]) - 0.03) < 0.004 and s["limits"]["n"] == 96 and s["checked"]
    assert os.path.exists(tmp_path / "q" / "station_qc.png")
    c, lag = SK._shift_corr(wav(0, 1.0), wav(0.013, 1.0), 100, 301, 10)
    assert c > 0.95 and abs(lag / fs - 0.013) < 0.002                  # between samples (0.02 s)
    assert SK.outlier_limits(np.zeros(5), np.zeros(5), fs) is None
