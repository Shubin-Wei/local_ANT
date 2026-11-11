"""Delivery keeps complete products independently of working files."""
import json
from pathlib import Path

import pytest

from lant import config, delivery, inspect, paths


def project(tmp_path):
    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path / "project")
    return cfg


def test_status_does_not_initialize_output(tmp_path):
    cfg = project(tmp_path)
    inspect.status(cfg)
    assert not Path(cfg["project"]["output_dir"]).exists()


def test_publish_model_and_validation_survive_work_removal(tmp_path):
    import shutil
    cfg = project(tmp_path)
    delivery.initialize(cfg)
    main = Path(paths.P(cfg, paths.INV, "main"))
    main.mkdir(parents=True)
    (main / "DSurfTomo.inMeasure.dat").write_text("model\n")
    (main / "surfdata.dat").write_text("input\n")
    check = main.parent / "checkerboards" / "check1"
    check.mkdir(parents=True)
    (check / "recovered.dat").write_text("recovery\n")
    delivery.publish(cfg, "run")
    root = Path(cfg["project"]["output_dir"])
    shutil.rmtree(str(root / "work"))
    assert (root / "results/model/DSurfTomo.inMeasure.dat").read_text() == "model\n"
    assert (root / "results/model/surfdata.dat").read_text() == "input\n"
    assert (root / "results/validation/checkerboards/check1/recovered.dat").exists()
    assert json.loads((root / ".lant/layout.json").read_text())["version"] == 1


def test_failed_publication_preserves_previous_tree(tmp_path, monkeypatch):
    source, target = tmp_path / "source", tmp_path / "target"
    source.mkdir()
    target.mkdir()
    (target / "model").write_text("previous")
    monkeypatch.setattr(delivery.shutil, "copytree", lambda *a: (_ for _ in ()).throw(OSError("full disk")))
    with pytest.raises(OSError, match="full disk"):
        delivery._replace_tree(source, target)
    assert (target / "model").read_text() == "previous"


def test_moving_project_keeps_clean_state_and_archive_readable(tmp_path):
    from lant import data, manifest, pipeline, report
    from util import make_raw
    make_raw(tmp_path)
    root = tmp_path / "project"
    root.mkdir()
    file = root / "config.yaml"
    file.write_text("project: {output_dir: ., array_type: surface}\n"
                    "data: {raw_dir: %s, station_file: %s, velocity_range_kms: [0.5, 3.0]}\n"
                    "clean: {nproc: 1}\n" % (tmp_path / "raw", tmp_path / "sta.csv"))
    pipeline.run_step(config.load(str(file)), "clean")
    moved = tmp_path / "moved"
    root.rename(moved)
    cfg = config.load(str(moved / "config.yaml"))
    assert manifest.status(report.resolved(cfg))["clean"]["state"] == "done"
    assert all(Path(r["file"]).exists() for r in data.processing_rows(cfg))
