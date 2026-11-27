"""Station table, record index, station matching, reading a station-day, configuration and lant init."""

import numpy as np
import pytest

from lant import config, init, records, stations


def test_station_table_forms(tmp_path):
    p = tmp_path / "a.csv"
    p.write_text("network,station,longitude,latitude,elevation\nXX,A1,100.0,30.0,12\nXX,A2,100.01,30.01,\n")
    t = stations.read_table(str(p))
    assert set(t) == {"XX.A1", "XX.A2"} and t["XX.A1"]["elev"] == 12 and np.isnan(t["XX.A2"]["elev"])
    p.write_text("# comment\nsta lon lat\nB1 359.0 -10\nB2 0.5 -10.1\n")
    t = stations.read_table(str(p))
    assert set(t) == {"B1", "B2"} and t["B1"]["lon"] == -1.0 and t["B1"]["network"] == ""
    p.write_text("name x y\nA 1 2\nB 1 2\n")
    with pytest.raises(stations.StationError, match="must name the columns"):
        stations.read_table(str(p))
    p.write_text("sta lon lat\nA 1 2\nA 1 3\n")
    with pytest.raises(stations.StationError, match="twice"):
        stations.read_table(str(p))


def _row(net, sta, cha="DPZ", fs=100.0, start=0.0, end=3600.0, loc=""):
    return dict(file="f", size=1, mtime=1, network=net, station=sta, location=loc, channel=cha, fs=fs, start=start,
                end=end, cmpinc="")  # fmt: skip


def test_matching_is_exact_and_reports_what_does_not_match():
    rows = [_row("XX", "A"), _row("XX", "B"), _row("YY", "B"), _row("XX", "C", cha="DPE"), _row("XX", "D")]
    rows += [_row("XX", "A", cha="DPZ", fs=500.0, end=600.0)]
    table = {"XX.A": dict(network="XX", station="A"), "B": dict(network="", station="B"),
             "XX.C": dict(network="XX", station="C"), "E": dict(network="", station="E")}  # fmt: skip
    found, notes = records.match(rows, table, "*Z")
    assert set(found) == {"XX.A"}  # B: two networks share the code; C: no vertical channel; E: no records
    assert all(r["fs"] == 100.0 for r in found["XX.A"])  # the sensor with most data
    text = " ".join(notes)
    assert "networks XX, YY" in text and "not in the table" in text and "XX.D" in text and "without *Z" in text


def _mseed(path, sta, data, fs, t0, net="XX", cha="DPZ"):
    import obspy

    tr = obspy.Trace(np.asarray(data, np.float32))
    tr.stats.update(dict(network=net, station=sta, channel=cha, sampling_rate=fs, starttime=obspy.UTCDateTime(t0)))
    tr.write(str(path), format="MSEED")


def test_index_and_station_day(tmp_path):
    raw = tmp_path / "raw"
    (raw / "sub").mkdir(parents=True)
    fs, t0 = 20.0, 1.7e9 - (1.7e9 % 86400)
    rng = np.random.default_rng(0)
    x = rng.standard_normal(int(3600 * fs))
    _mseed(raw / "a.mseed", "A", x[: len(x) // 2], fs, t0)
    _mseed(raw / "sub" / "a2.mseed", "A", x[len(x) // 2 + 2000 :], fs, t0 + (len(x) // 2 + 2000) / fs)  # a gap
    (raw / "notes.txt").write_text("not seismic data")
    idx = str(tmp_path / "out" / "record_index.csv")
    rows, bad = records.build_index(str(raw), idx, 1, log=lambda s: None)
    assert len(rows) == 2 and len(bad) == 1 and bad[0][0].endswith("notes.txt")
    rows2, _ = records.build_index(str(raw), idx, 1, log=lambda s: None)  # read again from the index
    assert [r["start"] for r in rows2] == [r["start"] for r in rows]
    y, ok = records.read_day(rows, t0, 5.0)  # valid only beyond the guard of max(1 s, 50 samples) = 10 s
    assert (
        y.size == 86400 * 5 and ok[60:1000].all() and not ok[:40].any() and not ok[int(1900 * 5)] and ok[int(2400 * 5)]
    )
    assert not ok[int(3700 * 5) :].any()
    assert records.days({"A": rows}) == [t0]
    with pytest.raises(records.DataError, match="whole multiple"):
        records.read_day(rows, t0, 3.0)


def test_configuration_rules(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text("output_dir: out\ncc:\n  fs_hz: 50\n  fmax_hz: 30\n")
    cfg = config.load(str(p))
    assert cfg["output_dir"] == str(tmp_path / "out") and cfg["cc"]["overlap"] == 0.5
    with pytest.raises(config.ConfigError, match="fill in"):
        config.check_cc(cfg)
    cfg["data"].update(raw_dir="r", station_file="s")
    cfg["cc"].update(fmin_hz=0.5, maxlag_s=20, segment_s=60)
    with pytest.raises(config.ConfigError, match="0.4 x cc.fs_hz"):
        config.check_cc(cfg)
    cfg["cc"]["fmax_hz"] = 15
    cfg["cc"]["segment_s"] = 30
    with pytest.raises(config.ConfigError, match="circular"):
        config.check_cc(cfg)
    p.write_text("cc:\n  fs_hzz: 50\n")
    with pytest.raises(config.ConfigError, match="did you mean 'fs_hz'"):
        config.load(str(p))


def test_init_suggestions_and_written_file(tmp_path):
    rng = np.random.default_rng(3)
    table = {"S%d" % k: dict(lon=100 + 0.05 * a, lat=30 + 0.05 * b) for k, (a, b) in enumerate(rng.random((30, 2)))}
    sug, geom = init.suggest(table, [100.0, 250.0], (0.3, 3.0))
    v = {k: x for k, (x, _) in sug.items()}
    assert v["cc.fs_hz"] in (10.0, 25.0, 50.0) and v["cc.fmax_hz"] <= 0.4 * v["cc.fs_hz"] + 1e-9
    assert v["cc.segment_s"] >= 2 * v["cc.maxlag_s"] and 0 < v["cc.fmin_hz"] < v["cc.fmax_hz"]
    path = tmp_path / "c.yaml"
    values = {"output_dir": str(tmp_path), "data.velocity_range_kms": [0.3, 3.0]}
    values.update(v)
    init.write(str(path), values, {k: r for k, (_, r) in sug.items()})
    cfg = config.load(str(path))
    assert cfg["cc"]["fs_hz"] == v["cc.fs_hz"] and cfg["data"]["velocity_range_kms"] == [0.3, 3.0]
    assert "lant init:" in path.read_text() and set(cfg) - {"_file"} == set(config.defaults())
