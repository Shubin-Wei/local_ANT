"""Delivery keeps complete products independently of working files."""

import json
from pathlib import Path

import pytest

from lant import assistant_rules, config, delivery, inspect, paths


def project(tmp_path):
    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path / "project")
    return cfg


def test_status_does_not_initialize_output(tmp_path):
    cfg = project(tmp_path)
    inspect.status(cfg)
    assert not Path(cfg["project"]["output_dir"]).exists()


MODEL_1 = "  99.97642  30.10929   0.50000   2.48578\n  99.97642  30.10929   0.00000   2.28578\n"


def _main(inv, model, tag):
    """A main project as DSurfTomo leaves it (inputs and the final model)."""
    main = inv / "main"
    main.mkdir(parents=True, exist_ok=True)
    for name in ("DSurfTomo.in", "surfdata.dat", "MOD"):
        (main / name).write_text("%s %s\n" % (name, tag))
    if model:
        (main / "DSurfTomo.inMeasure.dat").write_text(model)
    return main


def test_results_hold_readable_products_without_solver_projects(tmp_path, monkeypatch):
    import shutil

    from lant.inversion import run as RUN

    monkeypatch.setattr(RUN, "status", lambda d: "done")
    cfg = project(tmp_path)
    delivery.initialize(cfg)
    inv = Path(paths.P(cfg, paths.INV))
    _main(inv, MODEL_1, "1")
    check = inv / "checkerboards" / "check1"
    check.mkdir(parents=True)
    (check / "DSurfTomo.inSyn.dat").write_text("recovery\n")
    (inv / "evaluation").mkdir()
    (inv / "evaluation" / "evaluation.json").write_text("{}\n")
    delivery.publish(cfg, "evaluate")
    root = Path(cfg["project"]["output_dir"])
    shutil.rmtree(str(root / "work"))
    # depth first, values exactly as DSurfTomo wrote them
    assert (root / "results/model/vs_model.csv").read_text() == (
        "longitude_deg,latitude_deg,depth_km,vs_km_s\n99.97642,30.10929,0.00000,2.28578\n"
        "99.97642,30.10929,0.50000,2.48578\n"
    )
    assert sorted(p.name for p in (root / "results/model").iterdir()) == ["vs_model.csv"]
    assert sorted(p.name for p in (root / "results/validation").iterdir()) == ["evaluation.json"]
    assert not (root / "results/reproduce").exists()  # the solver files stay in work/06_inversion/main
    assert "vs_model.csv" in (root / "README.md").read_text()
    assert json.loads((root / ".lant/layout.json").read_text())["version"] == delivery.LAYOUT_VERSION


def test_dispersion_results_leave_measurements_and_pair_figures_in_work(tmp_path):
    cfg = project(tmp_path)
    disp = Path(paths.P(cfg, paths.DISP))
    (disp / "pair_figures").mkdir(parents=True)
    (disp / "pair_figures/spread_A_B.png").write_bytes(b"png")
    for name in ("curves.csv", "disp.h5", "median_vs_reference.txt"):
        (disp / name).write_text(name)
    delivery.publish(cfg, "disp")
    kept = Path(cfg["project"]["output_dir"]) / "results/dispersion"
    assert sorted(p.name for p in kept.iterdir()) == ["curves.csv"]
    assert (disp / "disp.h5").is_file() and (disp / "median_vs_reference.txt").is_file()


def test_step_products_are_available_before_inversion_and_keep_their_origins(tmp_path, monkeypatch):
    import yaml

    from lant import artifacts, manifest

    cfg = project(tmp_path)
    root = Path(cfg["project"]["output_dir"])
    delivery.initialize(cfg)
    for step, figures in delivery.STEP_FIGURES.items():
        before = artifacts.snapshot(cfg)
        for folder, name, target in figures:
            source = Path(paths.P(cfg, folder, name))
            source.parent.mkdir(parents=True, exist_ok=True)
            source.write_bytes(target.encode())
        delivery.publish(cfg, step)
        artifacts.register(cfg, before)
    # A stale zones figure must not be exported when zoned references are disabled.
    assert not (root / "results/figures/3_reference_zones.png").exists()
    assert (root / "results/figures/3_fj_phase_velocity.png").read_bytes() == b"3_fj_phase_velocity.png"
    assert not Path(paths.P(cfg, paths.REFERENCE, "fj.png")).exists()  # moved, not copied
    assert not (root / "results/model").exists()
    effective = yaml.safe_load((root / "results/parameters.yaml").read_text())
    assert effective["project"]["output_dir"] == ".."
    assert not (root / "results/parameters").exists()
    assert effective["reference"] == cfg["reference"]
    state = {s: dict(state="done", run=1, time="t") for s in delivery.STEP_FIGURES}
    state["reference"]["run"] = 2  # the retained spectrum is from the previous reference run
    monkeypatch.setattr(manifest, "status", lambda *a, **k: state)
    delivery.write_report(cfg)
    report = (root / "REPORT.md").read_text()
    assert "results/parameters.yaml" in report
    assert "![3_fj_phase_velocity](results/figures/3_fj_phase_velocity.png)" in report
    assert report.count("OLDER than the latest run 2") == 2  # spectrum and the parameters (last written by reference)


def test_borrowed_reference_does_not_publish_a_stale_fj_spectrum(tmp_path):
    cfg = project(tmp_path)
    cfg["reference"]["source"] = "file"
    source = Path(paths.P(cfg, paths.REFERENCE, "fj.png"))
    source.parent.mkdir(parents=True)
    source.write_bytes(b"old spectrum")
    delivery.publish(cfg, "reference")
    root = Path(cfg["project"]["output_dir"])
    assert not (root / "results/figures/3_fj_phase_velocity.png").exists()
    assert (root / "results/parameters.yaml").exists()


def test_step_publication_protects_unknown_files_and_survives_copy_failure(tmp_path, monkeypatch):
    from lant import artifacts

    cfg = project(tmp_path)
    root = Path(cfg["project"]["output_dir"])
    source = Path(paths.P(cfg, "01_clean", "coverage.png"))
    source.parent.mkdir(parents=True)
    source.write_bytes(b"first")
    before = artifacts.snapshot(cfg)
    delivery.publish(cfg, "clean")
    artifacts.register(cfg, before)
    source.write_bytes(b"new")
    marks = (root / ".lant/published.json").read_bytes()
    replace = delivery.os.replace

    def failed_move(a, b):
        if str(a) == str(source):
            Path(b).write_bytes(b"partial")
            raise OSError("copy failed")
        return replace(a, b)

    monkeypatch.setattr(delivery.os, "replace", failed_move)
    with pytest.raises(OSError, match="copy failed"):
        delivery.publish(cfg, "clean")
    target = root / "results/figures/0_data_coverage.png"
    assert target.read_bytes() == b"first"
    assert (root / ".lant/published.json").read_bytes() == marks
    assert not list(target.parent.glob(".deliver-*"))
    monkeypatch.setattr(delivery.os, "replace", replace)
    foreign = root / "results/parameters.yaml"
    before = artifacts.snapshot(cfg)
    foreign.unlink()
    artifacts.register(cfg, before)  # no longer lant's: what is put there now is the user's
    foreign.write_text("user file\n")
    with pytest.raises(RuntimeError, match="unregistered retained files"):
        delivery.publish(cfg, "cc")
    assert foreign.read_text() == "user file\n"


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
        delivery._swap(cfg, target, lambda new: delivery.shutil.copytree(str(source), str(new / "x")))
    assert (target / "model").read_text() == "previous"
    assert sorted(p.name for p in target.parent.iterdir()) == ["model"]


def _evaluated_project(tmp_path, monkeypatch):
    """A project whose main model and one of two checkerboards have run (solver status faked)."""
    from lant import artifacts
    from lant.inversion import run as RUN

    cfg = project(tmp_path)
    delivery.initialize(cfg)
    inv = Path(paths.P(cfg, paths.INV))
    for d in ("checkerboards/cb_done", "checkerboards/cb_unrun", "evaluation"):
        (inv / d).mkdir(parents=True)
    _main(inv, MODEL_1, "1")
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
    # --rebuild or a removed main directory: the new project is built, has not run, and its run is interrupted
    before = artifacts.snapshot(cfg)
    shutil.rmtree(str(inv / "main"))
    _main(inv, None, "2")
    finished.discard(str(inv / "main"))
    for step in ("inversion", "run"):
        delivery.publish(cfg, step)
    with pytest.raises(RuntimeError, match="retained model is kept"):
        delivery.publish(cfg, "evaluate")
    artifacts.register(cfg, before)
    assert "2.28578" in (results / "model/vs_model.csv").read_text()


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
    from util import make_raw

    from lant import data, manifest, pipeline, report

    make_raw(tmp_path)
    root = tmp_path / "project"
    root.mkdir()
    file = root / "config.yaml"
    file.write_text(
        "project: {output_dir: ., array_type: surface}\n"
        "data: {raw_dir: %s, station_file: %s, velocity_range_kms: [0.5, 3.0]}\n"
        "resources: {cores: 1}\n" % (tmp_path / "raw", tmp_path / "sta.csv")
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


def test_report_opens_with_state_numbers_judgements_and_figures(tmp_path, monkeypatch):
    from lant import manifest

    cfg = project(tmp_path)
    root = Path(cfg["project"]["output_dir"])
    delivery.initialize(cfg)
    (root / "results/figures").mkdir(parents=True)
    for name in ("1_stations", "6_checkerboard_cb_2", "6_checkerboard_cb_1", "5_depth_slices"):
        (root / "results/figures" / (name + ".png")).write_bytes(b"png")
    (root / ".lant/summary").mkdir(parents=True)
    (root / ".lant/summary/clean.json").write_text('{"stations": 12, "days": 3}')
    state = {
        "clean": dict(state="done", run=1, time="t"),
        "stack": dict(state="to confirm", run=1, time="t", reason="1 item"),
    }
    monkeypatch.setattr(manifest, "status", lambda *a, **k: state)
    manifest._save(cfg, {"steps": {"stack": {"to_confirm": ["remove the 3 candidate sub-stacks?"]}}})
    delivery.write_report(cfg)
    text = (root / "REPORT.md").read_text()
    order = [text.index(h) for h in ("## Key numbers", "## For you to judge", "## Main figures", "## Steps")]
    assert order == sorted(order)
    assert "12 stations, 3 processing days" in text and "remove the 3 candidate sub-stacks?" in text
    assert "your decision" in text
    figures = [line for line in text.splitlines() if line.startswith("![")]
    assert figures == [
        "![1_stations](results/figures/1_stations.png)",
        "![5_depth_slices](results/figures/5_depth_slices.png)",
        "![6_checkerboard_cb_1](results/figures/6_checkerboard_cb_1.png)",
    ]


def test_a_new_project_carries_the_rules_for_assistants(tmp_path):
    cfg = project(tmp_path)
    root = Path(cfg["project"]["output_dir"])
    cfg["_config_file"] = str(root / "config.yaml")
    delivery.initialize(cfg)
    rules = (root / "AGENTS.md").read_text()
    assert "scratch/" in rules and "lant note" in rules and "lant status -c config.yaml --json" in rules
    assert not (root / "CLAUDE.md").exists()
    assert not (root / "GEMINI.md").exists()
    assert "scratch/" in (root / "README.md").read_text()


@pytest.mark.parametrize(
    "host,filename", [("gemini", "GEMINI.md"), ("claude", "CLAUDE.md"), ("codex", "AGENTS.md"), (None, "AGENTS.md")]
)
def test_only_the_rules_file_of_the_chosen_assistant_is_written(tmp_path, host, filename):
    """One file with the full rules text, from the one template; an existing file is never overwritten."""
    cfg = project(tmp_path)
    cfg["project"]["ai_assistant"] = host
    root = Path(cfg["project"]["output_dir"])
    cfg["_config_file"] = str(root / "config.yaml")
    delivery.initialize(cfg)
    assert (root / filename).read_text() == assistant_rules.PROJECT_AGENTS.replace("CONFIG", "config.yaml")
    others = {"AGENTS.md", "CLAUDE.md", "GEMINI.md"} - {filename}
    assert not any((root / o).exists() for o in others)
    (root / filename).write_text("User instructions\n")
    delivery.initialize(cfg)
    assert (root / filename).read_text() == "User instructions\n"


def test_main_preview_keeps_model_even_if_plotting_fails(tmp_path, monkeypatch):
    from lant import artifacts, data
    from lant.inversion import preview, run

    cfg = project(tmp_path)
    delivery.initialize(cfg)
    main = _main(Path(paths.P(cfg, paths.INV)), MODEL_1, "preview")
    monkeypatch.setattr(run, "status", lambda _: "done")
    monkeypatch.setattr(data, "used_stations", lambda _: (_ for _ in ()).throw(RuntimeError("missing fixture")))
    before = artifacts.snapshot(cfg)
    preview.publish(cfg, str(main))
    artifacts.register(cfg, before)
    target = Path(cfg["project"]["output_dir"]) / "results/main_preview"
    assert (target / "DSurfTomo.inMeasure.dat").read_text() == MODEL_1
    note = json.loads((target / "status.json").read_text())
    assert not note["scientific_acceptance"] and note["figure_status"] == "failed"
    assert "missing fixture" in note["figure_error"]
    assert (target / "effective.yaml").is_file()
    (target / "unknown.txt").write_text("keep")
    with pytest.raises(RuntimeError, match="unregistered"):
        preview.publish(cfg, str(main))
    assert (target / "unknown.txt").read_text() == "keep"


def test_earlier_copies_in_results_are_retired(tmp_path, monkeypatch):
    """results/reproduce and the per-step results/parameters of earlier versions go when they are unchanged."""
    from lant import artifacts
    from lant.inversion import run as RUN

    monkeypatch.setattr(RUN, "status", lambda d: "done")
    cfg = project(tmp_path)
    delivery.initialize(cfg)
    root = Path(cfg["project"]["output_dir"])
    inv = Path(paths.P(cfg, paths.INV))
    _main(inv, MODEL_1, "1")
    (inv / "evaluation").mkdir()
    (inv / "evaluation" / "evaluation.json").write_text("{}\n")
    before = artifacts.snapshot(cfg)
    for rel in ("results/reproduce/MOD", "results/parameters/clean.yaml", "results/parameters/cc.yaml"):
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        (root / rel).write_text("old copy\n")
    artifacts.register(cfg, before)
    (root / "results/parameters/cc.yaml").write_text("edited\n")  # changed by hand: the tree stays
    delivery.publish(cfg, "evaluate")
    assert not (root / "results/reproduce").exists()
    assert (root / "results/parameters/cc.yaml").read_text() == "edited\n"
    assert (root / "results/parameters.yaml").is_file()


def test_main_preview_is_retired_by_the_evaluated_model_unless_changed(tmp_path, monkeypatch):
    from lant import artifacts, data
    from lant.inversion import preview, run

    cfg = project(tmp_path)
    delivery.initialize(cfg)
    inv = Path(paths.P(cfg, paths.INV))
    main = _main(inv, MODEL_1, "1")
    (inv / "evaluation").mkdir()
    (inv / "evaluation" / "evaluation.json").write_text("{}\n")
    monkeypatch.setattr(run, "status", lambda _: "done")
    monkeypatch.setattr(data, "used_stations", lambda _: (_ for _ in ()).throw(RuntimeError("no figure")))
    root = Path(cfg["project"]["output_dir"])

    def published(step_fn):
        before = artifacts.snapshot(cfg)
        step_fn()
        artifacts.register(cfg, before)

    published(lambda: preview.publish(cfg, str(main)))
    (root / "results/main_preview/status.json").write_text("edited by hand\n")
    published(lambda: delivery.publish(cfg, "evaluate"))
    assert (root / "results/main_preview/status.json").read_text() == "edited by hand\n"  # changed: kept
    published(lambda: preview.publish(cfg, str(main)))
    published(lambda: delivery.publish(cfg, "evaluate"))
    assert not (root / "results/main_preview").exists() and (root / "results/model/vs_model.csv").is_file()
    assert "results/main_preview" not in json.loads((root / ".lant/published.json").read_text())


def test_step_times_come_from_the_latest_completed_run(tmp_path):
    """The steps table of REPORT.md shows the wall time of the latest completed run of each step."""
    import json

    from lant import delivery, records

    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path)
    log = Path(records.rec_dir(cfg, create=True)) / "run_log.jsonl"
    with open(log, "w") as fh:
        for e in (
            dict(step="cc", event="done", seconds=50.0),
            dict(step="cc", event="failed", seconds=3.0),
            dict(step="cc", event="done", seconds=136.0),
            dict(step="run", event="done", seconds=1256.2),
        ):
            fh.write(json.dumps(e) + "\n")
        fh.write("not json\n")
    assert delivery.step_seconds(cfg) == {"cc": 136.0, "run": 1256.2}
    assert [delivery._duration(s) for s in (None, 52.0, 1256.2, 9000.0)] == ["", "52 s", "21 min", "2.5 h"]
