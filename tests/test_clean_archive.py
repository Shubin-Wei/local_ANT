"""Reusing a cleaned station-day must not hide same-second changes in the samples."""

import os

import numpy as np
import obspy

from lant import clean_archive, record_index


def test_same_second_raw_edit_rebuilds_archive(tmp_path):
    path = tmp_path / "record.sac"
    t0 = obspy.UTCDateTime(2024, 1, 1)
    trace = obspy.Trace(
        np.arange(100, dtype=np.float32),
        header=dict(network="XX", station="A01", channel="HHZ", sampling_rate=10.0, starttime=t0),
    )
    trace.write(str(path), format="SAC")
    stamp = 1_700_000_000_100_000_000
    os.utime(path, ns=(stamp, stamp))
    rows = record_index._scan_one((str(path), "SAC"))[1]
    options = dict(
        archive=str(tmp_path / "archive"),
        format="SAC",
        min_len=0.0,
        encoding="auto",
        codes={"XX.A01": ("XX", "A01", "", "HHZ")},
    )
    task = ("XX.A01", rows, t0.timestamp, options)
    _, original, _ = clean_archive._clean_station_day(task)
    assert clean_archive._clean_station_day(task)[2]["reused"] == 1
    trace.data += 5
    size = path.stat().st_size
    trace.write(str(path), format="SAC")
    os.utime(path, ns=(stamp + 100_000_000, stamp + 100_000_000))
    assert path.stat().st_size == size
    _, rebuilt, stats = clean_archive._clean_station_day(task)
    assert not stats.get("reused")
    assert rebuilt[0]["file"] == original[0]["file"]
    np.testing.assert_array_equal(obspy.read(rebuilt[0]["file"])[0].data, trace.data)
