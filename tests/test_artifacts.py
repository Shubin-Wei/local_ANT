"""Cleanup never adopts unknown files or follows links outside the project."""

import json
from pathlib import Path

import pytest

from lant import artifacts, config, records


def setup_project(tmp_path):
    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path)
    records.rec_dir(cfg, create=True)
    return cfg


def test_audit_and_preview_preserve_unknown_and_modified_files(tmp_path):
    cfg = setup_project(tmp_path)
    archive = tmp_path / "work/01_clean/archive"
    archive.mkdir(parents=True)
    unknown = archive / "user.txt"
    unknown.write_text("keep")
    before = artifacts.snapshot(cfg)
    cache = archive / "station.mseed"
    cache.write_text("registered")
    artifacts.register(cfg, before)
    plan = artifacts.clean(cfg, ["archive"])
    assert [f["path"] for f in plan["files"]] == ["work/01_clean/archive/station.mseed"]
    assert cache.exists() and unknown.exists()
    assert "work/01_clean/archive/user.txt" in artifacts.audit(cfg)["unregistered"]
    cache.write_text("modified")
    assert artifacts.clean(cfg, ["archive"])["protected"] == ["work/01_clean/archive/station.mseed"]


def test_cleanup_rejects_parent_symlink(tmp_path):
    cfg = setup_project(tmp_path / "project")
    root = Path(cfg["project"]["output_dir"])
    before = artifacts.snapshot(cfg)
    cache = root / "work/02_ccf/day.h5"
    cache.parent.mkdir(parents=True)
    cache.write_text("original")
    artifacts.register(cfg, before)
    outside = tmp_path / "outside"
    cache.parent.rename(outside)
    cache.parent.symlink_to(outside, target_is_directory=True)
    with pytest.raises(RuntimeError, match="symlinks"):
        artifacts.clean(cfg, ["cc"])
    assert (outside / "day.h5").read_text() == "original"


def test_cleanup_honors_lock_and_retains_deliverables(tmp_path, monkeypatch):
    from lant import manifest

    cfg = setup_project(tmp_path)
    cfg["data"]["raw_dir"] = str(tmp_path)
    cache = tmp_path / "work/02_ccf/day_2024001.h5"
    cache.parent.mkdir(parents=True)
    before = artifacts.snapshot(cfg)
    cache.write_text("cache")
    model = tmp_path / "results/model/model.dat"
    model.parent.mkdir(parents=True)
    model.write_text("model")
    artifacts.register(cfg, before)
    manifest._save(cfg, {"steps": {"cc": {}}})
    monkeypatch.setattr(manifest, "status", lambda *a, **k: {"evaluate": {"state": "done"}})
    import lant.run_report

    monkeypatch.setattr(lant.run_report, "write_report", lambda cfg: None)
    artifacts.clean(cfg, ["cc"], execute=True)
    assert not cache.exists() and model.read_text() == "model"
    assert manifest.load(cfg)["steps"]["cc"]["cache_cleaned"]


def test_unsafe_inventory_paths_cannot_become_cache():
    assert not artifacts.is_cache("../work/02_ccf/file", ["cc"])
    assert not artifacts.is_cache("/work/02_ccf/file", ["cc"])
    assert not artifacts.is_cache("results/model/file", ["cc", "archive"])


def test_step_cannot_overwrite_an_unregistered_figure(tmp_path):
    cfg = setup_project(tmp_path)
    figure = tmp_path / "results/figures/2_ccf_moveout_custom.png"
    figure.parent.mkdir(parents=True)
    figure.write_bytes(b"user figure")
    with pytest.raises(RuntimeError, match="protected"):
        with records.Step(cfg, "figures"):
            pytest.fail("the exporter must not start")
    assert figure.read_bytes() == b"user figure"


def test_interrupted_step_is_recovered_instead_of_blocking_the_rerun(tmp_path):
    import signal
    import subprocess
    import sys

    from lant import manifest

    child = (
        "import os, sys, time\n"
        "from lant import config, records\n"
        "cfg = config.defaults(); cfg['project']['output_dir'] = sys.argv[1]\n"
        "with records.Step(cfg, 'cc'):\n"
        "    d = os.path.join(sys.argv[1], 'work', '02_ccf'); os.makedirs(d, exist_ok=True)\n"
        "    open(os.path.join(d, 'day_2024001.h5'), 'w').write('partial')\n"
        "    sys.__stdout__.write('ready\\n'); sys.__stdout__.flush()\n"
        "    time.sleep(60)\n"
    )
    process = subprocess.Popen(
        [sys.executable, "-c", child, str(tmp_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        cwd=str(Path(__file__).parents[1]),
    )
    assert process.stdout.readline().strip() == "ready"
    process.send_signal(signal.SIGKILL)  # no cleanup at all (out of memory, kill -9)
    process.wait()
    cfg = setup_project(tmp_path)
    assert (tmp_path / ".lant/pending.json").exists()
    from lant import inspect

    shown = inspect.status(cfg)
    assert shown["steps"]["cc"]["state"] == "failed" and "interrupted" in shown["steps"]["cc"]["reason"]
    with records.Step(cfg, "cc"):  # was refused: unregistered files of the step
        pass
    assert "work/02_ccf/day_2024001.h5" in artifacts.load(cfg)
    assert not (tmp_path / ".lant/pending.json").exists()
    events = [json.loads(line)["event"] for line in (tmp_path / ".lant/run_log.jsonl").read_text().splitlines()]
    assert "interrupted" in events
    assert manifest.load(cfg)["steps"]["cc"]["state"] == "failed"  # recorded at recovery, before this run
    assert artifacts.audit(cfg)["unregistered"] == []


def test_status_during_a_write_shows_the_running_step(tmp_path):
    from lant import inspect
    from lant.provenance import Lock

    cfg = setup_project(tmp_path)
    (tmp_path / "work/02_ccf").mkdir(parents=True)
    (tmp_path / "work/02_ccf/day_2024001.h5").write_text("being written")
    with Lock(records.rec_dir(cfg), "cc"):
        st = inspect.status(cfg)
    assert st["steps"]["cc"]["state"] == "running" and st["next"]["command"] is None
    assert st["artifacts"]["unregistered"] == [] and st["artifacts"]["deferred"]
    assert "work/02_ccf/day_2024001.h5" in inspect.status(cfg)["artifacts"]["unregistered"]


def test_snapshot_tolerates_files_renamed_during_the_scan(tmp_path, monkeypatch):
    cfg = setup_project(tmp_path)
    (tmp_path / "work").mkdir()
    (tmp_path / "work/kept.txt").write_text("kept")
    (tmp_path / "work/x.tmp").write_text("renamed by a writer")
    real = Path.lstat
    monkeypatch.setattr(
        Path,
        "lstat",
        lambda self: (_ for _ in ()).throw(FileNotFoundError(self)) if self.name == "x.tmp" else real(self),
    )
    assert list(artifacts.snapshot(cfg)) == ["work/kept.txt"]


def test_cleanup_only_takes_the_cache_files_the_steps_write():
    assert artifacts.is_cache("work/02_ccf/day_2024001.h5", ["cc"])
    assert artifacts.is_cache("work/01_clean/archive/ZZ.A01/ZZ.A01..HHZ.2024001.mseed.done.json", ["archive"])
    assert not artifacts.is_cache("work/02_ccf/check.py", ["cc"])  # e.g. written there during a run
    assert not artifacts.is_cache("work/01_clean/archive/notes.txt", ["archive"])


def test_archived_copy_keeps_the_audit_clean(tmp_path):
    import shutil

    cfg = setup_project(tmp_path / "a")
    before = artifacts.snapshot(cfg)
    (tmp_path / "a/results").mkdir()
    (tmp_path / "a/results/model.dat").write_text("model")
    artifacts.register(cfg, before)
    shutil.copytree(str(tmp_path / "a"), str(tmp_path / "b"))  # copy2: keeps modification times
    copy = setup_project(tmp_path / "b")
    assert artifacts.audit(copy)["modified"] == []


def test_stray_files_at_the_top_of_the_project_are_reported(tmp_path):
    cfg = setup_project(tmp_path)
    cfg["_config_file"] = str(tmp_path / "config.yaml")
    (tmp_path / "config.yaml").write_text("project: {output_dir: .}\n")  # the user's own files
    (tmp_path / "notes.md").write_text("# Notes\n")
    (tmp_path / "check_dispersion.py").write_text("print(1)\n")  # what an assistant may leave behind
    (tmp_path / "tmp_plots").mkdir()
    (tmp_path / "tmp_plots/a.png").write_bytes(b"png")
    (tmp_path / "scratch").mkdir()
    (tmp_path / "scratch/extract.py").write_text("print(2)\n")  # the place for temporary files
    found = artifacts.audit(cfg)
    assert found["unregistered"] == ["check_dispersion.py", "tmp_plots/a.png"]
    assert found["scratch"] == dict(files=1, bytes=9)


@pytest.mark.parametrize("step", ["run", "evaluate"])
def test_step_protects_unknown_forward_cache_files(tmp_path, step):
    cfg = setup_project(tmp_path)
    entry = tmp_path / ".lant/forward_cache/user.npz"
    entry.parent.mkdir()
    entry.write_bytes(b"user data")
    with pytest.raises(RuntimeError, match="unregistered files"):
        with records.Step(cfg, step):
            pytest.fail("must not overwrite an unknown cache entry")
    assert entry.read_bytes() == b"user data"


def test_forward_cleanup_is_optional_and_preserves_other_state(tmp_path, monkeypatch):
    from lant import manifest, run_report

    cfg = setup_project(tmp_path)
    cfg["data"]["raw_dir"] = str(tmp_path)
    folder = tmp_path / ".lant/forward_cache"
    folder.mkdir()
    before = artifacts.snapshot(cfg)
    (folder / "registered.npz").write_bytes(b"cached")
    artifacts.register(cfg, before)
    (folder / "unknown.npz").write_bytes(b"user")
    manifest._save(cfg, {"steps": {"cc": {"state": "done"}}})
    monkeypatch.setattr(manifest, "status", lambda *a, **kw: {"evaluate": {"state": "done"}})
    monkeypatch.setattr(run_report, "write_report", lambda _: None)
    plan = artifacts.clean(cfg, ["forward"])
    assert [v["path"] for v in plan["files"]] == [".lant/forward_cache/registered.npz"]
    artifacts.clean(cfg, ["forward"], execute=True)
    assert not (folder / "registered.npz").exists()
    assert (folder / "unknown.npz").read_bytes() == b"user"
    assert "cache_cleaned" not in manifest.load(cfg)["steps"]["cc"]


def test_pruning_removes_only_unread_registered_solver_files_of_check_projects(tmp_path):
    from lant import config, paths

    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path / "project")
    root = Path(cfg["project"]["output_dir"])
    inv = Path(paths.P(cfg, paths.INV))
    cb, main = inv / "checkerboards" / "cb_1km", inv / "main"
    names = ["DSurfTomo.inMeasure.dat.iter001", "residualFirst.dat", "residualLast.dat", "velmap2dRc.dat"]
    keep = ["DSurfTomo.in", "DSurfTomo.in.log", "DSurfTomo.inSyn.dat", "Vs_model.real", "MOD", "run.log"]
    with artifacts.writing(cfg, "run"):
        for d in (cb, main):
            d.mkdir(parents=True, exist_ok=True)
            for n in names + keep:
                (d / n).write_text(n)
    (cb / "velmap2dRc.dat").write_text("changed by hand")
    out = artifacts.prune_check_projects(cfg, [str(cb), str(main)])
    assert out["files"] == 3 and out["kept"] == [(cb / "velmap2dRc.dat").relative_to(root).as_posix()]
    assert sorted(p.name for p in cb.iterdir()) == sorted(keep + ["velmap2dRc.dat"])
    assert sorted(p.name for p in main.iterdir()) == sorted(names + keep)  # the main project is never pruned
