"""Native GMT API operation validation and failed/cross-filesystem delivery."""

import pytest

from lant.plotting import gmt


def test_api_argv_preserves_label_spaces_without_shell_quotes():
    pytest.importorskip("pygmt")
    from lant.plotting.pygmt_backend import _call

    calls = []

    class Session:
        def call_module(self, module, args):
            calls.append((module, args))

    _call(Session(), "gmt basemap -R0/1/0/1 -JX10c -Bxa+l'Distance (km)'", {})
    assert calls == [("basemap", ["-R0/1/0/1", "-JX10c", "-Bxa+lDistance (km)"])]


def test_pygmt_failure_preserves_old_image(tmp_path, monkeypatch):
    from lant.plotting import pygmt_backend

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    out = tmp_path / "old"
    (tmp_path / "old.png").write_bytes(b"old complete")

    def failure(*args, **kwargs):
        raise RuntimeError("render failed")

    monkeypatch.setattr(pygmt_backend, "render", failure)
    with pytest.raises(gmt.GMTError, match="render failed"):
        gmt.run(str(scratch), [], [], str(out), formats=["png"], backend="pygmt")
    assert (tmp_path / "old.png").read_bytes() == b"old complete"


def test_pygmt_copies_before_promoting_any_format(tmp_path, monkeypatch):
    from lant.plotting import pygmt_backend

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    out = tmp_path / "old"
    for suffix in ("png", "pdf"):
        (tmp_path / ("old." + suffix)).write_bytes(b"old complete")
        (scratch / ("lant_figure." + suffix)).write_bytes(b"new complete")
    monkeypatch.setattr(pygmt_backend, "render", lambda *a, **k: None)
    import shutil

    copy = shutil.copy2

    def fail_second(source, target):
        if str(source).endswith(".pdf"):
            raise OSError("copy failed")
        return copy(source, target)

    monkeypatch.setattr(gmt.shutil, "copy2", fail_second)
    with pytest.raises(OSError, match="copy failed"):
        gmt.run(str(scratch), [], [], str(out), formats=["png", "pdf"], backend="pygmt")
    assert all((tmp_path / ("old." + s)).read_bytes() == b"old complete" for s in ("png", "pdf"))


def test_relief_mosaic_preserves_coordinates_units_and_nodata(tmp_path):
    import numpy as np

    rasterio = pytest.importorskip("rasterio")
    xr = pytest.importorskip("xarray")
    from rasterio.transform import from_origin

    from lant.plotting.relief import _mosaic_python

    inputs = []
    for i, values in enumerate(([[1, 2], [3, -9999]], [[5, 6], [7, 8]])):
        path = tmp_path / ("tile%d.tif" % i)
        with rasterio.open(
            path,
            "w",
            driver="GTiff",
            width=2,
            height=2,
            count=1,
            dtype="int16",
            crs="EPSG:4326",
            transform=from_origin(100 + i, 31, 0.5, 0.5),
            nodata=-9999,
        ) as dst:
            dst.write(np.array(values, dtype="int16"), 1)
        inputs.append(path)
    target = str(tmp_path / "mosaic.nc")
    _mosaic_python(inputs, 100, 102, 30, 31, target)
    with xr.open_dataset(target) as grid:
        assert grid.node_offset == 1 and grid.z.units == "m"
        assert np.array_equal(grid.lon, [100.25, 100.75, 101.25, 101.75])
        assert np.array_equal(grid.lat, [30.25, 30.75])
        assert np.isnan(grid.z.values[0, 1])
        assert np.array_equal(grid.z.values[1], [1, 2, 5, 6])


def test_pdf_map_coordinates_and_api_labels(tmp_path):
    import shutil
    import subprocess

    pytest.importorskip("pygmt")
    if not shutil.which("pdftotext"):
        pytest.skip("PDF text extraction is optional")
    from lant.plotting.maps import _geo_format

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    target = tmp_path / "map"
    statement = "gmt basemap -R100/100.08/30/30.08 -JM12c -BWSne -Bxa0.05 -Bya0.025 " + _geo_format(0.05, 0.025)
    gmt.run(
        str(scratch),
        [],
        [statement, "gmt colorbar -Cviridis -Dx0c/-2c+w10c/0.3c+h -Bxaf+l'Distance (km)'"],
        str(target),
        formats=["pdf"],
        backend="pygmt",
    )
    text = subprocess.check_output(["pdftotext", "-layout", str(target) + ".pdf", "-"], text=True)
    assert "30.075" in text and "30.75" not in text
    assert "Distance (km)" in text and "’Distance" not in text
