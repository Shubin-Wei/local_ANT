"""The stack step: excluded station-days (the station diagnosis: test_station_qc.py)."""

import numpy as np


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


def test_hour_sum_is_the_numpy_form_bit_for_bit():
    """The compiled hour sum gives the numbers of x = float64(h) * n; total[g] += x; daily[g] += float32(x)."""
    from lant.stack import _add_hour

    rng = np.random.default_rng(1)
    hourly = rng.standard_normal((6, 50)).astype(np.float32)
    cnt = rng.integers(0, 40, 6).astype(np.float64)
    gi = np.array([4, 0, 2, 7, 5, 1])
    total = rng.standard_normal((8, 50))
    daily = rng.standard_normal((8, 50)).astype(np.float32)
    t_ref, d_ref = total.copy(), daily.copy()
    x = hourly.astype(np.float64) * cnt[:, None]
    t_ref[gi] += x
    d_ref[gi] += x.astype(np.float32)
    _add_hour(hourly, cnt, gi, total, daily)
    assert np.array_equal(total, t_ref) and np.array_equal(daily, d_ref)


def test_rerun_after_confirming_exclusions_reuses_the_stack_of_every_hour(tmp_path, monkeypatch):
    """The run after `lant confirm stack --exclude-candidates` does not stack the hours again when the product is
    the stack of every hour with the same inputs; it removes the candidates from it."""
    import h5py

    from lant import config, paths
    from lant import stack as S

    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path)
    built = []

    def fake_build(cfg, out, exclude_hours=None, block_hours=None):
        built.append(out)
        with h5py.File(out, "w") as h:
            h["stack"] = np.ones((1, 5), np.float32)
            h["daily"] = np.ones((2, 1, 5), np.float32)
            h["daily_count"] = np.full((2, 1), 10, np.int32)
            h["nseg"] = np.array([20])
            h["pairs"] = np.array([b"A_B"])
            h["days"] = np.array([b"2024001_00", b"2024001_12"])
            for k in ("az_deg", "dist_km", "lat1", "lat2", "lon1", "lon2"):
                h[k] = np.zeros(1)
        return dict(pairs=1)

    monkeypatch.setattr(S, "build_stack", fake_build)
    monkeypatch.setattr(S, "all_hours_key", lambda cfg: "same inputs")
    monkeypatch.setattr(S, "station_qc", lambda cfg, raw, d: ({("A", "2024001_00"): "shift"}, {}))
    first = S.run(cfg)
    assert len(built) == 1 and first["product"].startswith("the stack of every hour")
    cfg["stack"]["exclude_candidates"] = True
    second = S.run(cfg)
    assert len(built) == 1 and second["all_hours"]["reused_from_last_run"]
    assert second["product"]["dropped_segments"] == 10
    with h5py.File(paths.stack_product(cfg)) as h:
        assert h["nseg"][0] == 10
    cfg["stack"]["exclude_candidates"] = False  # the product is no longer the stack of every hour: built again
    S.run(cfg)
    assert len(built) == 2


def test_sub_stack_chunks_stay_small_and_memory_is_planned():
    """1000 stations at 50 Hz and 30 s lag: one chunk of every pair would be 6 GB (older HDF5 refuses 4 GiB); the
    chunks hold whole correlations of as many pairs as fit into 64 MiB.  The estimate grows with the sub-stacks."""
    from lant import stack

    npair, nt = 499500, 3001
    c = stack.daily_chunks(npair, nt)
    assert c[0] == 1 and c[2] == nt and 4 * c[1] * nt <= stack.CHUNK_BYTES and c[1] > 1000
    assert stack.daily_chunks(10, 101) == (1, 10, 101)
    one, many = stack.memory_estimate(npair, nt, 1, 30), stack.memory_estimate(npair, nt, 24, 720)
    assert 2.4e10 < one < many
