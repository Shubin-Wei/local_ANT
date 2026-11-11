"""Cleanup never adopts unknown files or follows links outside the project."""
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
    plan = artifacts.clean(cfg, ["cc"])
    assert plan["protected"] == ["work/02_ccf/day.h5"]
    assert (outside / "day.h5").read_text() == "original"


def test_cleanup_honors_lock_and_retains_deliverables(tmp_path, monkeypatch):
    from lant import manifest
    cfg = setup_project(tmp_path)
    cfg["data"]["raw_dir"] = str(tmp_path)
    cache = tmp_path / "work/02_ccf/day.h5"
    cache.parent.mkdir(parents=True)
    before = artifacts.snapshot(cfg)
    cache.write_text("cache")
    model = tmp_path / "results/model/model.dat"
    model.parent.mkdir(parents=True)
    model.write_text("model")
    artifacts.register(cfg, before)
    manifest._save(cfg, {"steps": {"cc": {}}})
    monkeypatch.setattr(manifest, "status", lambda *a, **k: {"evaluate": {"state": "done"}})
    monkeypatch.setattr(artifacts, "write_report", lambda cfg: None, raising=False)
    import lant.delivery
    monkeypatch.setattr(lant.delivery, "write_report", lambda cfg: None)
    artifacts.clean(cfg, ["cc"], execute=True)
    assert not cache.exists() and model.read_text() == "model"
    assert manifest.load(cfg)["steps"]["cc"]["cache_cleaned"]


def test_unsafe_inventory_paths_cannot_become_cache():
    assert not artifacts.is_cache("../work/02_ccf/file", ["cc"])
    assert not artifacts.is_cache("/work/02_ccf/file", ["cc"])
    assert not artifacts.is_cache("results/model/file", ["cc", "archive"])
