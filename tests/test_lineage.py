"""Parameter provenance survives config edits, failures and partial redraws."""

import json
from pathlib import Path

import pytest
import yaml

from lant import config, lineage, manifest, pipeline


def project(tmp_path):
    file = tmp_path / "config.yaml"
    file.write_text(
        "project: {output_dir: out, array_type: surface}\n"
        "data: {raw_dir: raw, station_file: stations.csv, min_trace_length_s: 5}\n"
    )
    return file, config.load(str(file))


def test_actual_parameters_are_frozen_before_running_and_retained_after_failure(tmp_path, monkeypatch):
    file, cfg = project(tmp_path)
    monkeypatch.setattr(manifest, "gate", lambda *args: None)
    monkeypatch.setattr(manifest, "raw_inputs", lambda *args: "fixture")

    def clean(frozen):
        file.write_text(file.read_text().replace("min_trace_length_s: 5", "min_trace_length_s: 9"))
        coverage = Path(frozen["project"]["output_dir"]) / "work/01_clean/coverage.png"
        coverage.parent.mkdir(parents=True, exist_ok=True)
        coverage.write_bytes(b"test image")
        return {}

    monkeypatch.setitem(pipeline.FUNCS, "clean", clean)
    pipeline.run_step(cfg, "clean")
    actual = lineage.parameters(config.load(str(file)), "results/figures/0_data_coverage.png")
    assert actual["execution"]["params"]["data.min_trace_length_s"] == 5
    assert actual["execution"]["effective"]["data"]["min_trace_length_s"] == 5
    difference = lineage.parameters(config.load(str(file)), "clean")["differences"]
    assert difference["data.min_trace_length_s"] == {"used": 5, "next": 9}

    def failed(cfg):
        raise RuntimeError("read failed")

    monkeypatch.setitem(pipeline.FUNCS, "clean", failed)
    with pytest.raises(RuntimeError, match="read failed"):
        pipeline.run_step(config.load(str(file)), "clean")
    retained = lineage.parameters(config.load(str(file)), "results/figures/0_data_coverage.png")
    assert retained["origin"] == actual["origin"]
    assert retained["execution"]["state"] == "done"
    root = Path(cfg["project"]["output_dir"])
    assert len(list((root / ".lant/runs").glob("*.json"))) == 2
    assert manifest.load(cfg)["steps"]["clean"]["state"] == "failed"
    index = yaml.safe_load((root / "results/parameters.yaml").read_text())
    assert index["products"]["results/figures/0_data_coverage.png"]["record"] == actual["origin"]["record"]


def test_crlf_configuration_is_not_reported_as_edited_during_execution(tmp_path):
    file, _ = project(tmp_path)
    file.write_bytes(file.read_bytes().replace(b"\n", b"\r\n"))
    cfg = config.load(str(file))
    lineage.begin(cfg, "figures", {})
    lineage.finish(cfg, "figures", "done")
    record = json.loads((Path(cfg["project"]["output_dir"]) / cfg["_run_record"]).read_text())
    assert record["configuration_edited_during_run"] is False


def test_each_redrawn_image_keeps_its_own_execution_and_upstream(tmp_path):
    _, cfg = project(tmp_path)
    root = Path(cfg["project"]["output_dir"])
    images = root / "results/figures"
    images.mkdir(parents=True)
    record = lineage.begin(cfg, "figures", {"which": ["stations", "ccf_moveout"]})
    for name in ("stations.png", "moveout.png"):
        (images / name).write_bytes(b"first")
        lineage.mark(cfg, "results/figures/" + name, "figures", {"band_hz": [1, 3]})
    lineage.finish(cfg, "figures", "done")
    first = cfg["_run_record"]
    # A later partial redraw replaces only the successful image; a failed image retains its old record.
    lineage.begin(cfg, "figures", {"which": ["ccf_moveout"]})
    (images / "moveout.png").write_bytes(b"second")
    lineage.mark(cfg, "results/figures/moveout.png", "figures", {"band_hz": [2, 4]})
    lineage.finish(cfg, "figures", "partial", {"failed": ["other renderer failed"]})
    catalog = lineage.index(cfg)
    assert catalog["products"]["results/figures/stations.png"]["record"] == first
    assert catalog["products"]["results/figures/moveout.png"]["record"] == cfg["_run_record"]
    assert catalog["products"]["results/figures/moveout.png"]["parameters"]["band_hz"] == [2, 4]
    assert json.loads((root / first).read_text())["id"] == record["id"]
    (images / "moveout.png").write_bytes(b"manually edited after publication")
    queried = lineage.parameters(cfg, "results/figures/moveout.png")
    assert "changed after publication" in queried["note"]


def test_result_queries_normalise_paths_and_reject_paths_outside_the_project(tmp_path):
    _, cfg = project(tmp_path)
    root = Path(cfg["project"]["output_dir"])
    image = root / "results/figures/example.png"
    image.parent.mkdir(parents=True)
    image.write_bytes(b"figure")
    lineage.begin(cfg, "figures", {})
    lineage.mark(cfg, "results/figures/example.png", "figures")
    lineage.finish(cfg, "figures", "done")
    for path in (str(image), "./results/figures/example.png", "results/figures/../figures/example.png"):
        assert lineage.parameters(cfg, path)["product"] == "results/figures/example.png"
    with pytest.raises(ValueError, match="inside the project"):
        lineage.parameters(cfg, "../outside.png")


def test_coverage_distinguishes_input_paths_from_retained_paths(tmp_path, monkeypatch):
    from lant.plotting import coverage, style

    source = tmp_path / "residuals.csv"
    source.write_text(
        "src_lat,src_lon,rcv_lat,rcv_lon,period_s,kept_by_rule\n30,100,30.1,100.1,1,1\n30,100,30.2,100.2,1,0\n"
    )
    monkeypatch.setattr(style, "save", lambda fig, *args: fig)
    result = coverage.path_coverage(config.defaults(), str(source), "unused")
    fig = result["files"]
    try:
        assert "retained 1 / 2" in fig.axes[0].get_title()
        assert len(fig.axes[0].collections[0].get_segments()) == 2
        assert len(fig.axes[0].collections[1].get_segments()) == 1
        assert "not bent rays" in result["parameters"]["path_representation"]
    finally:
        style.plt.close(fig)


def test_explicit_band_is_not_rejected_by_an_unused_automatic_band_rule():
    from lant import auto

    cfg = config.defaults()
    user = {
        "correlation": {"fs": 10.0, "fmin_hz": 1.0, "fmax_hz": 4.0},
        "dispersion": {"fmin_hz": 2.0, "fmax_hz": 3.0, "qc": {"lam_c": 1.0}},
    }
    config._merge(cfg, user)
    cfg["_user"] = user
    params = auto.from_velocity(cfg, {"max_distance_km": 0.05}, [100.0])
    assert params["correlation.fmin_hz"]["value"] == 1.0
    assert params["dispersion.fmax_hz"]["value"] == 3.0
    assert "given in the configuration" in params["correlation.fmin_hz"]["reason"]


def test_display_band_is_validated_against_saved_computation(tmp_path):
    import h5py

    from lant.plotting.ccf import validate_display

    path = tmp_path / "stack.h5"
    with h5py.File(path, "w") as h:
        h.attrs.update(fs=25.0, maxlag=30.0, fmin=1.0, fmax=10.0)
    validate_display(str(path), [2, 5], 20)
    with pytest.raises(ValueError, match="saved correlation band"):
        validate_display(str(path), [0.5, 5])
    with pytest.raises(ValueError, match="saved lag window"):
        validate_display(str(path), [2, 5], 40)
