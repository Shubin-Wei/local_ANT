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


def test_publish_model_and_validation_survive_work_removal(tmp_path, monkeypatch):
    import shutil

    from lant.inversion import run as RUN

    monkeypatch.setattr(RUN, "status", lambda d: "done")
    cfg = project(tmp_path)
    delivery.initialize(cfg)
    main = Path(paths.P(cfg, paths.INV, "main"))
    main.mkdir(parents=True)
    (main / "DSurfTomo.inMeasure.dat").write_text("model\n")
    (main / "surfdata.dat").write_text("input\n")
    check = main.parent / "checkerboards" / "check1"
    check.mkdir(parents=True)
    (check / "recovered.dat").write_text("recovery\n")
    (main.parent / "evaluation").mkdir()
    (main.parent / "evaluation" / "evaluation.json").write_text("{}\n")
    delivery.publish(cfg, "evaluate")
    root = Path(cfg["project"]["output_dir"])
    shutil.rmtree(str(root / "work"))
    assert (root / "results/model/DSurfTomo.inMeasure.dat").read_text() == "model\n"
    assert (root / "results/model/surfdata.dat").read_text() == "input\n"
    assert (root / "results/validation/checkerboards/check1/recovered.dat").exists()
    assert json.loads((root / ".lant/layout.json").read_text())["version"] == 1


def test_failed_publication_preserves_previous_tree(tmp_path, monkeypatch):
    cfg = project(tmp_path)
    source, target = tmp_path / "source", Path(cfg["project"]["output_dir"]) / "results/model"
    source.mkdir()
    target.mkdir(parents=True)
    (target / "model").write_text("previous")
    import lant.artifacts

    monkeypatch.setattr(lant.artifacts, "load", lambda cfg: {"results/model/model": [0, 0]})
    monkeypatch.setattr(delivery.shutil, "copytree", lambda *a, **k: (_ for _ in ()).throw(OSError("full disk")))
    with pytest.raises(OSError, match="full disk"):
        delivery._swap(cfg, target, lambda new: delivery._copy_into(source, new))
    assert (target / "model").read_text() == "previous"
    assert sorted(p.name for p in target.parent.iterdir()) == ["model"]


def _evaluated_project(tmp_path, monkeypatch):
    """A project whose main model and one of two checkerboards have run (solver status faked)."""
    from lant import artifacts
    from lant.inversion import run as RUN

    cfg = project(tmp_path)
    delivery.initialize(cfg)
    inv = Path(paths.P(cfg, paths.INV))
    for d in ("main", "checkerboards/cb_done", "checkerboards/cb_unrun", "evaluation"):
        (inv / d).mkdir(parents=True)
    (inv / "main/DSurfTomo.in").write_text("input 1\n")
    (inv / "main/DSurfTomo.inMeasure.dat").write_text("model 1\n")
    (inv / "checkerboards/cb_done/DSurfTomo.inSyn.dat").write_text("recovered\n")
    (inv / "checkerboards/cb_unrun/DSurfTomo.in").write_text("not run\n")
    (inv / "evaluation/evaluation.json").write_text("{}\n")
    (inv / "vs1d.txt").write_text("start 1\n")
    finished = {str(inv / "main"), str(inv / "checkerboards/cb_done")}
    monkeypatch.setattr(RUN, "status", lambda d: "done" if str(d) in finished else "new")
    before = artifacts.snapshot(cfg)
    delivery.publish(cfg, "evaluate")
    artifacts.register(cfg, before)
    return cfg, inv, finished


def test_rebuilding_the_inversion_keeps_the_evaluated_model(tmp_path, monkeypatch):
    import shutil
    from lant import artifacts

    cfg, inv, finished = _evaluated_project(tmp_path, monkeypatch)
    results = Path(cfg["project"]["output_dir"]) / "results"
    assert sorted(p.name for p in (results / "validation/checkerboards").iterdir()) == ["cb_done"]
    # the configuration snapshot is written in staging but must resolve from its final place
    assert "output_dir: ../.." in (results / "reproduce/effective.yaml").read_text()
    # --rebuild or a removed main directory: the new project is built, has not run, and its run is interrupted
    before = artifacts.snapshot(cfg)
    shutil.rmtree(str(inv / "main"))
    (inv / "main").mkdir()
    (inv / "main/DSurfTomo.in").write_text("input 2\n")
    finished.discard(str(inv / "main"))
    for step in ("inversion", "run"):
        delivery.publish(cfg, step)
    with pytest.raises(RuntimeError, match="retained model is kept"):
        delivery.publish(cfg, "evaluate")
    artifacts.register(cfg, before)
    assert (results / "model/DSurfTomo.inMeasure.dat").read_text() == "model 1\n"
    assert (results / "model/DSurfTomo.in").read_text() == "input 1\n"
    assert (results / "reproduce/vs1d.txt").read_text() == "start 1\n"


def test_report_names_the_run_of_each_retained_tree(tmp_path, monkeypatch):
    from lant import manifest

    cfg, inv, _ = _evaluated_project(tmp_path, monkeypatch)
    root = Path(cfg["project"]["output_dir"])
    marks = json.loads((root / ".lant/published.json").read_text())
    assert marks["results/model"]["step"] == "evaluate" and marks["results/model"]["run"] == 1
    # evaluate has run again since (run 2) without being published
    monkeypatch.setattr(manifest, "status", lambda *a, **k: {"evaluate": dict(state="done", run=2, time="t")})
    delivery.write_report(cfg)
    text = (root / "REPORT.md").read_text()
    assert "results/model/" in text and "OLDER than the latest run" in text


def test_atomic_text_uses_ordinary_permissions(tmp_path):
    import os

    old = os.umask(0o022)
    try:
        delivery.atomic_text(tmp_path / "REPORT.md", "report\n")
        assert (tmp_path / "REPORT.md").stat().st_mode & 0o777 == 0o644
        os.chmod(str(tmp_path / "REPORT.md"), 0o664)
        delivery.atomic_text(tmp_path / "REPORT.md", "again\n")
        assert (tmp_path / "REPORT.md").stat().st_mode & 0o777 == 0o664
    finally:
        os.umask(old)


def test_publication_refuses_to_remove_unregistered_retained_files(tmp_path, monkeypatch):
    from lant.inversion import run as RUN

    monkeypatch.setattr(RUN, "status", lambda d: "done")
    cfg = project(tmp_path)
    root = Path(cfg["project"]["output_dir"])
    source = root / "work/06_inversion/main"
    source.mkdir(parents=True)
    (source / "model").write_text("new")
    target = root / "results/model"
    target.mkdir(parents=True)
    (target / "user-note.txt").write_text("keep")
    with pytest.raises(RuntimeError, match="unregistered retained files"):
        delivery.publish(cfg, "evaluate")
    assert (target / "user-note.txt").read_text() == "keep"


def test_moving_project_keeps_clean_state_and_archive_readable(tmp_path):
    from lant import data, manifest, pipeline, report
    from util import make_raw

    make_raw(tmp_path)
    root = tmp_path / "project"
    root.mkdir()
    file = root / "config.yaml"
    file.write_text(
        "project: {output_dir: ., array_type: surface}\n"
        "data: {raw_dir: %s, station_file: %s, velocity_range_kms: [0.5, 3.0]}\n"
        "clean: {nproc: 1}\n" % (tmp_path / "raw", tmp_path / "sta.csv")
    )
    pipeline.run_step(config.load(str(file)), "clean")
    moved = tmp_path / "moved"
    root.rename(moved)
    cfg = config.load(str(moved / "config.yaml"))
    assert manifest.status(report.resolved(cfg))["clean"]["state"] == "done"
    assert all(Path(r["file"]).exists() for r in data.processing_rows(cfg))


def test_report_says_when_a_retained_result_is_no_longer_current(tmp_path, monkeypatch):
    from lant import manifest

    cfg, inv, _ = _evaluated_project(tmp_path, monkeypatch)
    monkeypatch.setattr(manifest, "status", lambda *a, **k: {"evaluate": dict(state="out of date", run=1, time="t")})
    delivery.write_report(cfg)
    text = (Path(cfg["project"]["output_dir"]) / "REPORT.md").read_text()
    assert "that step is now out of date" in text and "OLDER" not in text
