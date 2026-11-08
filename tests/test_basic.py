"""Fast tests (about half a minute): configuration, cleaning of mixed formats, numerics, export and project files.
The full chain is tested with the synthetic example (examples/synthetic, docs/tutorial.md section 1.4).
    python -m pytest tests
"""
import os
import textwrap

import numpy as np
import pytest

from lant import config, signal, export
from lant.inversion import project


def write_cfg(tmp_path, text):
    p = tmp_path / "cfg.yaml"
    p.write_text(textwrap.dedent(text))
    return str(p)


def test_config_defaults_and_errors(tmp_path):
    p = write_cfg(tmp_path, """
        project: {output_dir: out}
        data: {raw_dir: raw, station_file: s.csv}
        """)
    cfg = config.load(p)
    assert cfg["correlation"]["fs"] == "auto" and cfg["correlation"]["seg_step_s"] == 30.0
    assert cfg["project"]["output_dir"] == str(tmp_path / "out")
    p = write_cfg(tmp_path, """
        project: {output_dir: out, array_type: line}
        data: {raw_dir: raw}
        """)
    cfg = config.load(p)
    assert cfg["stack"]["station_hour_qc"]["enabled"] and cfg["inversion"]["weight"] == 1.5
    p = write_cfg(tmp_path, """
        project: {output_dir: out}
        data: {raw_dir: raw, chanel: "*Z"}
        """)
    with pytest.raises(config.ConfigError, match="unknown configuration key 'data.chanel'"):
        config.load(p)
    p = write_cfg(tmp_path, """
        project: {output_dir: out}
        data: {raw_dir: raw}
        correlation: {seg_len_s: 30, maxlag_s: 20}
        """)
    with pytest.raises(config.ConfigError, match="seg_len_s"):
        config.load(p)


def _make_raw(tmp_path):
    import obspy
    rng = np.random.default_rng(0)
    t0 = obspy.UTCDateTime(2024, 1, 1)
    raw = tmp_path / "raw"
    for i, name in enumerate(["A01", "A02", "A03"]):
        tr = obspy.Trace(rng.standard_normal(int(3 * 3600 * 100)).astype(np.float32),
                         header=dict(network="", station=name, channel="HHZ", sampling_rate=100.0,
                                     starttime=t0 + (0.004 if i == 2 else 0.0003)))
        d = raw / ("deep" if i else "flat") / "x"
        d.mkdir(parents=True, exist_ok=True)
        if i == 0:
            tr.stats.network = "ZZ"
            tr.write(str(d / "rec1.mseed"), format="MSEED", encoding="FLOAT32")
            tr.copy().write(str(d / "dup.mseed"), format="MSEED", encoding="FLOAT32")        # exact duplicate
        else:
            tr.write(str(d / ("%s.sac" % name)), format="SAC")
    (raw / "notes.txt").write_text("not data")
    (tmp_path / "sta.csv").write_text("network,station,latitude,longitude,elevation\n"
                                      "ZZ,A01,30.0,100.0,10\nZZ,A02,30.01,100.0,20\nZZ,A03,30.0,100.01,\n")
    return raw


def test_clean_mixed_formats(tmp_path):
    from lant import clean, auto
    _make_raw(tmp_path)
    p = write_cfg(tmp_path, """
        project: {output_dir: out}
        data: {raw_dir: raw, station_file: sta.csv}
        clean: {nproc: 2}
        """)
    cfg = config.load(p)
    s = clean.run(cfg)
    assert s["stations"] == 3 and s["formats"] == {"MSEED": 2, "SAC": 2}
    assert s["names"]["matched"] == 3 and s["names"]["renamed"] == 2      # SAC without network matched by station code
    assert s["totals"]["duplicates_or_identical_overlaps"] == 1
    assert s["totals"]["shifted"] == 3          # 0.03 sample off the grid: start moved (A01 twice, A02)
    assert s["totals"]["interpolated"] == 1     # 0.4 sample off the grid: Lanczos interpolation (A03)
    assert os.path.exists(tmp_path / "out" / "01_clean" / "archive" / "ZZ.A02" / "ZZ.A02..HHZ.2024001.mseed")
    cfg = auto.apply(cfg)
    assert cfg["correlation"]["fs"] == 50.0 and cfg["data"]["day_start_utc_hour"] == 0


def test_station_tables_and_names(tmp_path):
    from lant import stations as ST
    tables = {
        "a.csv": "network,station,channel,latitude,longitude,elevation\nYY,1001,HHZ,35.1234,105.4321,306\nYY,1002,HHZ,35.2,105.5,483\n",
        "b.txt": "YY1001 105.4321 35.1234 306\nYY1002 105.5 35.2 483\n",
        "c.txt": "台站名\t经度\t纬度\nYY.1001\t105°25'55.56\"E\t35°7'24.24\"N\nYY.1002\t105.5\t35.2\n",
        "d.txt": "S01 100.04 30.07 335\nS02 100.01 30.08 517\n",
    }
    for name, text in tables.items():
        f = tmp_path / name
        f.write_text(text)
        rows, how = ST.read_table(str(f))
        assert len(rows) == 2 and abs(rows[0]["lat"] - (30.07 if name == "d.txt" else 35.1234)) < 1e-3, (name, rows, how)
    rows, _ = ST.read_table(str(tmp_path / "b.txt"))
    mapping, _, rep = ST.match_names(rows, {("YY", "1001"), ("", "1002")})
    assert mapping == {("YY", "1001"): "YY.1001", ("", "1002"): "YY.1002"}
    mapping, _, rep = ST.match_names([dict(network=None, station="A01", lat=1, lon=2, elev=0)], {("X", "A01"), ("Y", "A01")})
    assert not mapping and rep["ambiguous"]
    assert not ST._is_num("S01") and ST._is_num("105.43E") and abs(ST.parse_coord("40:30:00") - 40.5) < 1e-9


def test_numerics():
    f = signal.log_freqs(0.5, 10.5, 12)
    assert len(f[f <= 10.5]) == 53 and abs(f[12] - 1.0) < 1e-12
    m = np.array([1, 1, 0, 1, 1, 1, 1, 1, 1, 0], bool)
    assert signal.contiguous(m, 6).tolist() == [0, 0, 0, 1, 1, 1, 1, 1, 1, 0]
    d, az = signal.geodesic_km(100.0, 30.0, 100.0, 31.0)
    assert abs(d - 111.195) < 0.01 and abs(az) < 1e-6
    x = np.sin(2 * np.pi * 5.0 * np.arange(4000) / 50.0)[None, :]
    y = signal.notch_lines(x, 50.0, [5.0], half_width=0.1)
    assert np.abs(y[:, 1000:3000]).max() < 0.05


def test_one_segment():
    ok = np.array([0, 1, 1, 0, 1, 1, 1, 0], bool)
    snr = np.array([0, 50, 50, 0, 5, 5, 5, 0], float)
    assert export.longest_run(ok, snr).tolist() == [0, 0, 0, 0, 1, 1, 1, 0]
    ok = np.array([1, 1, 0, 1, 1], bool)
    snr = np.array([5, 5, 0, 9, 9], float)
    assert export.longest_run(ok, snr).tolist() == [0, 0, 0, 1, 1]      # tie: higher mean SNR


def test_project_files(tmp_path):
    curves = [dict(name="A_B", lat1=30.0, lon1=100.0, lat2=30.05, lon2=100.02, T=np.array([0.5, 1.0]), v=np.array([2.5, 2.7])),
              dict(name="A_C", lat1=30.0, lon1=100.0, lat2=30.02, lon2=100.06, T=np.array([1.0]), v=np.array([2.8]))]
    recs, periods = project.records_from_curves(curves)
    assert periods == [0.5, 1.0]
    info = project.write_surfdata(recs, str(tmp_path / "surfdata.dat"))
    lines = (tmp_path / "surfdata.dat").read_text().splitlines()
    assert lines[0] == "# 30.000000 100.000000 1 2 0" and info["travel_times"] == 3
    g = project.grid(info["stations"], (0.005, 0.005), refine=2)
    assert g["dvxd"] == 0.0025 and (g["nx"] - 1) % 2 == 0
    lat, lon = project.node_coords(g)
    assert lat[1] > 30.05 and lon[1] < 100.0              # stations inside the boundary ring
    z = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    s = project.depth_signs(z, [0, 0.2, 0.4, 0.6], [0.3])
    assert np.allclose(s, [1, 1, 1, 0, -1, -1, -1])
