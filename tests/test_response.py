"""Response-file edits invalidate both an in-memory inventory and clean-input provenance."""

import os

import obspy
from obspy.core.inventory import Inventory, Network

from lant import config, manifest, response


def test_response_edit_reloads_inventory_and_invalidates_clean_input(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "record.txt").write_text("raw record placeholder")
    table = tmp_path / "stations.csv"
    table.write_text("station,latitude,longitude\nA,30,100\n")
    path = tmp_path / "response.xml"
    Inventory(networks=[Network("XX")], source="source-a").write(str(path), format="STATIONXML")
    cfg = config.defaults()
    cfg["data"].update(raw_dir=str(raw), station_file=str(table), response=str(path))
    before = manifest.raw_inputs(cfg)
    first = response.load(str(path))
    assert response.load(str(path)) is first
    stamp, size = path.stat().st_mtime_ns, path.stat().st_size
    # Preserve mtime and size: ctime/file identity must still invalidate these caches.
    path.write_text(path.read_text().replace("source-a", "source-b"))
    os.utime(path, ns=(stamp, stamp))
    assert path.stat().st_size == size
    assert response.load(str(path)).source == "source-b"
    assert manifest.raw_inputs(cfg) != before
    assert isinstance(response.load(str(path)), obspy.Inventory)
