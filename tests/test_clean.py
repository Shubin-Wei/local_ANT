"""The clean step: input formats, names, sample grid, gaps, polarity, responses, processing days."""
import os

import numpy as np
import pytest

from lant import config
from util import REQ, write_cfg, make_raw


def test_clean_mixed_formats(tmp_path):
    """miniSEED and SAC in nested folders: duplicates dropped, records off the sample grid shifted or interpolated,
    names matched, the clean-step rules applied."""
    from lant import clean, auto
    make_raw(tmp_path)
    p = write_cfg(tmp_path, """
        project: {output_dir: out, array_type: surface}
        data: {raw_dir: raw, station_file: sta.csv, %s}
        clean: {nproc: 2}
        """ % REQ)
    cfg = config.load(p)
    s = clean.run(cfg)
    assert s["stations"] == 3 and s["formats"] == {"MSEED": 2, "SAC": 2}
    assert s["names"]["matched"] == 3 and s["names"]["renamed"] == 2      # SAC without network matched by station code
    assert s["totals"]["duplicates_or_identical_overlaps"] == 1
    assert s["totals"]["shifted"] == 3          # 0.03 sample off the grid: start moved (A01 twice, A02)
    assert s["totals"]["interpolated"] == 1     # 0.4 sample off the grid: Lanczos interpolation (A03)
    assert os.path.exists(tmp_path / "out" / "work" / "01_clean" / "archive" / "ZZ.A02" / "ZZ.A02..HHZ.2024001.mseed")
    cfg = auto.apply(cfg)
    assert cfg["data"]["day_start_utc_hour"] == 0 and cfg["stack"]["block_hours"] == 1   # 3 h of record
    assert cfg["correlation"]["fs"] == 50.0     # from the velocity range and the rates (fs_max 50)
    assert auto.unresolved(cfg, auto.CLEAN_KEYS) == [] and auto.unresolved(cfg, auto.REFERENCE_KEYS)
    assert s["response"]["removed"] is False and s["array_type"]["type"] == "surface"


def test_clean_long_names_nan_polarity(tmp_path):
    """Station names longer than miniSEED allows are archived under codes and read back; NaN samples are gaps; a
    downward vertical sensor (SAC cmpinc 180) is turned over."""
    import obspy
    from lant import clean, data as D
    t0 = obspy.UTCDateTime(2024, 1, 1)
    raw = tmp_path / "raw"
    raw.mkdir()
    rng = np.random.default_rng(0)
    base = rng.standard_normal(3600 * 50).astype(np.float32)
    for i, name in enumerate(["NODE0001", "NODE0002", "NODE0003"]):
        x = base.copy() if i < 2 else -base.copy()
        if i == 1:
            x[1000:1100] = np.nan
        tr = obspy.Trace(x, header=dict(network="", station=name, channel="DPZ", sampling_rate=50.0, starttime=t0))
        tr.stats.sac = dict(cmpinc=180.0 if i == 2 else 0.0, cmpaz=0.0)
        tr.write(str(raw / ("%s.sac" % name)), format="SAC")
    (tmp_path / "sta.csv").write_text("station,latitude,longitude\nNODE0001,30.0,100.0\nNODE0002,30.01,100.0\n"
                                      "NODE0003,30.0,100.01\n")
    p = write_cfg(tmp_path, "project: {output_dir: out, array_type: surface}\ndata: {raw_dir: raw, station_file: sta.csv, "
                            "%s}\nclean: {nproc: 1}\n" % REQ)
    cfg = config.load(p)
    s = clean.run(cfg)
    assert s["totals"]["nonfinite_samples"] == 100 and s["totals"]["polarity_reversed"] == 1
    rows = D.processing_rows(cfg)
    assert {r["station"] for r in rows} == {"XX.NODE0001", "XX.NODE0002", "XX.NODE0003"}
    by = {}
    for r in rows:
        by.setdefault(r["station"], []).append(r)
    y1, ok1 = D.read_day(by["XX.NODE0001"], t0.timestamp, 25.0)
    y3, ok3 = D.read_day(by["XX.NODE0003"], t0.timestamp, 25.0)
    y2, ok2 = D.read_day(by["XX.NODE0002"], t0.timestamp, 25.0)
    m = ok1 & ok3
    assert m.sum() > 1000 and np.allclose(y1[m], y3[m], atol=1e-3 * np.abs(y1[m]).max())   # turned over: equal
    assert np.isfinite(y2).all() and not ok2[500:560].any() and ok2[2000:2100].all()


def test_other_input_formats(tmp_path):
    """GSE2 records and a StationXML station table go through the clean step like miniSEED and a CSV table."""
    import obspy
    from obspy.core.inventory import Inventory, Network, Station
    from lant import pipeline
    raw = tmp_path / "raw"
    raw.mkdir()
    rng = np.random.default_rng(11)
    stas = []
    for k, code in enumerate(("G01", "G02", "G03")):
        tr = obspy.Trace((1000 * rng.standard_normal(100 * 7200)).astype(np.int32),
                         dict(network="ZZ", station=code, channel="HHZ", sampling_rate=100.0,
                              starttime=obspy.UTCDateTime(2024, 1, 1)))
        tr.write(str(raw / ("%s.gse" % code)), format="GSE2")
        stas.append(Station(code, 30.0 + 0.01 * (k % 2), 100.0 + 0.01 * (k // 2), 0.0))
    Inventory(networks=[Network("ZZ", stations=stas)]).write(str(tmp_path / "sta.xml"), format="STATIONXML")
    p = write_cfg(tmp_path, "project: {output_dir: out}\ndata: {raw_dir: raw, station_file: sta.xml, %s}\nclean: {nproc: 1}" % REQ)
    s = pipeline.run_step(config.load(p), "clean")
    assert s["formats"] == {"GSE2": 3} and s["stations"] == 3 and s["station_hours"] == 6.0


def test_segy_input(tmp_path):
    """SEG-Y without station codes: one file per station, the record cut into traces of 300 s (SEG-Y holds at most
    32767 samples per trace), named from the file name and joined by the clean step; a file with one trace per
    channel at the same times is not named and is reported."""
    import obspy
    from lant import pipeline
    raw = tmp_path / "raw"
    raw.mkdir()
    rng = np.random.default_rng(12)
    t0 = obspy.UTCDateTime(2024, 1, 1)
    rows = []
    for k, code in enumerate(("S01", "S02", "S03")):
        st = obspy.Stream([obspy.Trace(rng.standard_normal(30000).astype(np.float32),
                                       dict(sampling_rate=100.0, starttime=t0 + 300 * j)) for j in range(24)])
        st.write(str(raw / ("%s.sgy" % code)), format="SEGY", data_encoding=5)
        rows.append("%s,%.3f,%.3f,0" % (code, 30.0 + 0.01 * (k % 2), 100.0 + 0.01 * (k // 2)))
    st = obspy.Stream([obspy.Trace(rng.standard_normal(30000).astype(np.float32), dict(sampling_rate=100.0, starttime=t0))
                       for _ in range(3)])
    st.write(str(raw / "S01_shot.sgy"), format="SEGY", data_encoding=5)
    (tmp_path / "sta.csv").write_text("station,lat,lon,elev\n" + "\n".join(rows) + "\n")
    p = write_cfg(tmp_path, "project: {output_dir: out}\ndata: {raw_dir: raw, station_file: sta.csv, %s}\nclean: {nproc: 1}" % REQ)
    s = pipeline.run_step(config.load(p), "clean")
    assert s["formats"] == {"SEGY": 75} and s["stations"] == 3 and s["station_hours"] == 6.0
    rep = open(tmp_path / "out" / "work" / "01_clean" / "report.txt").read()
    assert "72 records without a station code in the header: named" in rep and "1 files hold several records" in rep
    # without the archive the later steps read the raw files with the code of the file (empty), not the table's
    import json
    from lant import data as D
    p2 = write_cfg(tmp_path, "project: {output_dir: out2}\ndata: {raw_dir: raw, station_file: sta.csv, %s}\n"
                             "clean: {nproc: 1, write_archive: false}" % REQ)
    pipeline.run_step(config.load(p2), "clean")
    idx = json.load(open(tmp_path / "out2" / "work" / "01_clean" / "archive_index.json"))["traces"]
    rows = [t for t in idx if t["station"] == "XX.S02"]
    x, ok = D.read_day(rows, t0.timestamp, 50.0)
    assert len(rows) == 24 and x is not None and ok.sum() > 0.95 * 2 * 3600 * 50 - 200


def test_snap_to_grid():
    """Records off the sample grid: the same placement for the archive and for direct reading (ceil when resampled)."""
    import obspy
    from lant.data import snap_to_grid
    t0 = 1.6e9
    for frac, kind in ((0.0, "on_grid"), (0.03, "shifted"), (0.2, "interpolated"), (-0.4, "interpolated")):
        tr = obspy.Trace(np.sin(np.arange(2000) / 7.0), header=dict(sampling_rate=100.0))
        tr.stats.starttime = obspy.UTCDateTime(t0 + (10 + frac) / 100.0)
        assert snap_to_grid(tr, t0, 100.0) == kind
        off = (tr.stats.starttime.timestamp - t0) * 100.0
        assert abs(off - round(off)) < 1e-3                  # float timestamps: ~1e-5 sample
        if kind == "interpolated":                        # resampled from the first grid point inside the record
            assert off >= 10 + frac - 1e-3


def test_fold_1s():
    """A waveform repeating every second is detected by the fold test; noise alone is not."""
    from lant.pulse import fold_1s
    rng = np.random.default_rng(13)
    fs, n = 100.0, 100 * 3600
    noise = rng.standard_normal(n)
    r0 = fold_1s([(37, noise)], fs)
    pulse = np.zeros(100)
    pulse[10:13] = 0.2
    r1 = fold_1s([(0, noise + np.tile(pulse, 3600))], fs)
    assert r0["ratio"] < r0["limit"] < r1["ratio"] and r1["seconds"] == 3600 and r0["seconds"] == 3599
    assert fold_1s([(0, noise)], 62.5) is None


def test_response_removal():
    """Two sensors with different responses record the same ground motion: after lant.response.remove the two
    records agree in phase inside the pre-filter band."""
    import obspy
    from obspy.core.inventory import Inventory, Network, Station, Channel
    from obspy.core.inventory.response import Response
    from lant import response as R
    fs, n = 100.0, 60000
    t0 = obspy.UTCDateTime(2024, 1, 1)
    rng = np.random.default_rng(14)
    g = rng.standard_normal(n)
    sos = __import__("scipy.signal", fromlist=["x"]).butter(4, [1, 10], "band", fs=fs, output="sos")
    g = __import__("scipy.signal", fromlist=["x"]).sosfiltfilt(sos, g)
    chans, recs = [], []
    for code, f0 in (("A", 5.0), ("B", 0.1)):                 # 5-Hz geophone, broadband-like sensor
        w0 = 2 * np.pi * f0
        resp = Response.from_paz(zeros=[0j, 0j], poles=[-0.7 * w0 + 0.714j * w0, -0.7 * w0 - 0.714j * w0], stage_gain=100.0)
        tr = obspy.Trace(g.copy(), dict(network="XX", station=code, channel="DPZ", sampling_rate=fs, starttime=t0))
        f = np.fft.rfftfreq(n, 1 / fs)
        h = resp.get_evalresp_response_for_frequencies(f[1:], output="VEL")
        spec = np.fft.rfft(g)
        spec[1:] *= h
        tr.data = np.fft.irfft(spec, n)
        recs.append(tr)
        chans.append(Station(code, 30.0, 100.0, 0.0, channels=[Channel("DPZ", "", 30.0, 100.0, 0.0, 0.0, response=resp,
                                                                       sample_rate=fs)]))
    inv = Inventory(networks=[Network("XX", stations=chans)])
    raw = np.angle(np.fft.rfft(recs[0].data)[3000:3600] * np.conj(np.fft.rfft(recs[1].data)[3000:3600]))
    assert np.abs(raw).mean() > 0.3                              # the raw records differ in phase around 5 Hz
    assert R.lookup(inv, "", "A", "", "DPZ", t0.timestamp + 1) == "XX.A..DPZ"       # header without a network code
    for tr in recs:
        R.remove(tr, inv, "XX.%s..DPZ" % tr.stats.station, [0.5, 1.0, 10.0, 12.5])
    a, b = (np.fft.rfft(tr.data[5000:55000]) for tr in recs)
    f = np.fft.rfftfreq(50000, 1 / fs)
    m = (f > 1.5) & (f < 9)
    assert np.abs(np.angle(a[m] * np.conj(b[m]))).max() < 0.02


def test_processing_day_labels():
    """Processing-day labels (data.day_label) and the processing-day windows."""
    import obspy
    from lant import data as D
    cfg = config.defaults("surface")
    t0 = obspy.UTCDateTime("2024-01-01T16:00:00").timestamp          # Beijing midnight of 2 January
    assert D.day_label(t0, cfg) == "2024001"                         # utc_start: UTC date of the start
    cfg["data"]["day_label"] = "UTC+8"
    assert D.day_label(t0, cfg) == "2024002" and D.label_date(t0, cfg) == "2024-01-02"
    assert D.label_offset_s("UTC-05:30") == -5.5 * 3600 and D.label_offset_s(8) == 8 * 3600
    with pytest.raises(ValueError):
        D.label_offset_s("+8")
    tr = dict(start=obspy.UTCDateTime("2023-12-31T16:00:00").timestamp, end=obspy.UTCDateTime("2024-01-03T16:00:00").timestamp)
    cfg["data"]["day_start_utc_hour"] = 16
    cfg["data"]["start_date"] = cfg["data"]["end_date"] = "2024-01-02"
    days = D.day_windows([tr], cfg)
    assert list(days) == ["2024002"] and days["2024002"][0] == t0       # one date, one processing day
    cfg["data"]["day_label"] = "utc_start"
    days = D.day_windows([tr], cfg)
    assert list(days) == ["2024002"] and days["2024002"][0] == t0 + 86400
