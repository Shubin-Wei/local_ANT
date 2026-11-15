"""Long tasks cannot silently mix imported and subsequently edited source files."""

import json
import os
import shutil
import subprocess
import sys

import pytest

from lant import config, records, runtime_identity


def test_source_drift_detects_edits_additions_deletions_and_ignores_mtime(tmp_path, monkeypatch):
    first = tmp_path / "solver.py"
    first.write_text("old")
    second = tmp_path / "lazy.py"
    second.write_text("lazy")
    monkeypatch.setattr(runtime_identity, "PACKAGE", tmp_path)
    monkeypatch.setattr(runtime_identity, "SOURCES", runtime_identity.source_files())
    os.utime(first, (1, 1))
    runtime_identity.assert_unchanged()
    first.write_text("new")
    os.utime(first, (1, 1))  # identical size/mtime does not hide changed implementation
    second.unlink()
    (tmp_path / "added.py").write_text("new module")
    assert runtime_identity.changed_sources() == ["added.py", "lazy.py", "solver.py"]
    with pytest.raises(RuntimeError, match="source changed after process startup"):
        runtime_identity.assert_unchanged()


def test_worker_rejects_a_different_source_snapshot(tmp_path):
    package = tmp_path / "lant"
    package.mkdir()
    for name in ("__init__.py", "runtime_identity.py"):
        shutil.copyfile(runtime_identity.PACKAGE / name, package / name)
    result = subprocess.run(
        [sys.executable, "-c", "import lant"],
        cwd=tmp_path,
        env=dict(os.environ, PYTHONPATH=str(tmp_path), LANT_EXPECTED_SOURCE_SHA256="different"),
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0 and "worker source differs" in result.stderr


def test_step_records_frozen_runtime_and_fails_if_source_changes(tmp_path, monkeypatch):
    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path)
    original = runtime_identity.source_files()
    changed = dict(original, **{"lazy.py": "edited"})
    with pytest.raises(RuntimeError, match="source changed"):
        with records.Step(cfg, "run"):
            monkeypatch.setattr(runtime_identity, "source_files", lambda: changed)
    log = [json.loads(line) for line in (tmp_path / ".lant/run_log.jsonl").read_text().splitlines()]
    assert log[0]["runtime"]["source_sha256"] == runtime_identity.SOURCE_SHA256
    assert log[-1]["event"] == "failed" and "source changed" in log[-1]["error"]


def test_source_drift_does_not_hide_original_failure(tmp_path, monkeypatch):
    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path)
    with pytest.raises(ValueError, match="actual failure"):
        with records.Step(cfg, "run"):
            monkeypatch.setattr(runtime_identity, "source_files", lambda: {})
            raise ValueError("actual failure")


def test_code_identity_ignores_docstrings_and_comments_only():
    from lant.runtime_identity import _code_only

    a = 'def f(x):\n    """Twice x."""\n    return 2 * x  # double\n'
    b = 'def f(x):\n    """Two times x, reworded."""\n\n    # a comment\n    return 2 * x\n'
    c = 'def f(x):\n    """Twice x."""\n    return 3 * x\n'
    assert _code_only(a) == _code_only(b) != _code_only(c)


def test_step_record_carries_the_code_of_its_step(monkeypatch):
    """A processing step recorded by other code (another lant version) is out of date: its numerics differ."""
    from lant import config, manifest, runtime_identity

    cfg = config.defaults()
    for step in runtime_identity.STEP_SOURCES:
        assert manifest.numerical_dependencies(cfg, step)["code"] == runtime_identity.code_digest(step)
    before = manifest.numerical_dependencies(cfg, "disp")
    monkeypatch.setitem(runtime_identity._CODE, "disp", "0" * 16)
    assert manifest.numerical_dependencies(cfg, "disp") != before
