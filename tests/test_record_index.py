"""Header-cache reuse must skip work without retaining same-second stale metadata."""

import json
import os

import numpy as np
import obspy
import pytest

from lant import config, record_index


def fixture(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    path = raw / "record.sac"
    trace = obspy.Trace(
        np.arange(100, dtype=np.float32),
        header=dict(network="XX", station="A01", channel="HHZ", sampling_rate=10.0),
    )
    trace.write(str(path), format="SAC")
    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path / "out")
    cfg["data"]["raw_dir"] = str(raw)
    cfg["resources"]["cores"] = 1
    return cfg, path, trace


def test_unchanged_cache_opens_no_pool_and_does_not_rewrite(tmp_path, monkeypatch):
    import multiprocessing

    cfg, _, _ = fixture(tmp_path)
    first = record_index.build_index(cfg)
    cache = tmp_path / "out/work/01_clean/header_cache.json"
    stamp = cache.stat().st_mtime_ns

    def unexpected(*args, **kwargs):
        pytest.fail("unchanged headers must not start a pool or read raw records")

    monkeypatch.setattr(multiprocessing, "Pool", unexpected)
    monkeypatch.setattr(record_index, "_scan_one", unexpected)
    again = record_index.build_index(cfg)
    assert again["traces"][0]["sta"] == first["traces"][0]["sta"] == "A01"
    assert again["n_files"] == first["n_files"] == 1
    assert cache.stat().st_mtime_ns == stamp


@pytest.mark.parametrize("legacy", [False, True])
def test_same_second_same_size_header_edit_is_rescanned(tmp_path, legacy):
    cfg, path, trace = fixture(tmp_path)
    first_ns = 1_700_000_000_100_000_000
    os.utime(path, ns=(first_ns, first_ns))
    assert record_index.build_index(cfg)["traces"][0]["sta"] == "A01"
    size = path.stat().st_size
    if legacy:
        cache = tmp_path / "out/work/01_clean/header_cache.json"
        payload = json.loads(cache.read_text())
        payload["files"][str(path)][0] = [size, first_ns // 1_000_000_000]
        cache.write_text(json.dumps(payload))
    trace.stats.station = "B02"
    trace.write(str(path), format="SAC")
    second_ns = first_ns + 100_000_000
    os.utime(path, ns=(second_ns, second_ns))
    assert path.stat().st_size == size
    assert record_index.build_index(cfg)["traces"][0]["sta"] == "B02"


def test_files_of_stations_outside_the_table_are_read_only_at_their_ends(tmp_path):
    """A miniSEED file whose first and last records are of a station no rule can match to the table keeps only those
    two records in the index (header_only); when the table names the station, the file is read in full."""
    cfg, _, _ = fixture(tmp_path)
    raw = tmp_path / "raw"
    for sta in ("A01", "Z99"):
        tr = obspy.Trace(np.arange(3000, dtype=np.int32), header=dict(network="XX", station=sta, channel="HHZ",
                                                                       sampling_rate=10.0))  # fmt: skip
        obspy.Stream([tr]).write(str(raw / ("%s.mseed" % sta)), format="MSEED", reclen=512, encoding="STEIM2")
    table = tmp_path / "stations.csv"
    table.write_text("network,station,latitude,longitude,elevation\nXX,A1,30.0,100.0,0\n")  # A1 = A01 (rule 5)
    cfg["data"]["station_file"] = str(table)
    index = record_index.build_index(cfg)
    by = {}
    for t in index["traces"]:
        by.setdefault(t["sta"], []).append(t)
    assert index["header_only"] == 1
    assert all(t.get("header_only") for t in by["Z99"]) and not any(t.get("header_only") for t in by["A01"])
    assert by["Z99"][0]["end"] == pytest.approx(by["A01"][0]["end"])  # the extent from the last record
    table.write_text("network,station,latitude,longitude,elevation\nXX,A1,30.0,100.0,0\nXX,Z99,30.1,100.1,0\n")
    again = record_index.build_index(cfg)
    assert again["header_only"] == 0 and not any(t.get("header_only") for t in again["traces"])


def test_read_error_is_not_kept_in_the_header_cache(tmp_path, monkeypatch):
    """A file that could not be read (OSError: disk, network share) is scanned again by the next run; a file that is
    not seismic data stays cached as such."""
    cfg, path, _ = fixture(tmp_path)
    original = record_index.read_quiet

    def unreadable(*args, **kwargs):
        raise OSError(5, "Input/output error")

    monkeypatch.setattr(record_index, "read_quiet", unreadable)
    try:
        record_index.build_index(cfg)
    except Exception:
        pass  # no readable record at all may stop the index; the cache is what matters here
    cache = json.loads((tmp_path / "out/work/01_clean/header_cache.json").read_text())
    assert str(path) not in cache["files"]
    monkeypatch.setattr(record_index, "read_quiet", original)
    assert record_index.build_index(cfg)["traces"][0]["sta"] == "A01"
