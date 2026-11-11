"""Resource decisions use effective limits and leave scientific configuration untouched."""

from collections import namedtuple

import pytest

from lant import config, resources


def test_effective_memory_includes_nested_cgroup_limit(tmp_path):
    proc, group = tmp_path / "proc", tmp_path / "cgroup"
    (proc / "self").mkdir(parents=True)
    (proc / "meminfo").write_text("MemAvailable: 8000000 kB\n")
    (proc / "self/cgroup").write_text("0::/job/task\n")
    (group / "job/task").mkdir(parents=True)
    (group / "job/memory.max").write_text("2000000000")
    (group / "job/memory.current").write_text("1500000000")
    (group / "memory.max").write_text("max")
    assert resources.available_memory(proc, group) == 500000000


def test_preflight_rejects_shortage_without_changing_priors(tmp_path, monkeypatch):
    cfg = config.defaults()
    prior = list(cfg["data"]["velocity_range_kms"])
    estimate = dict(memory_per_process_gb=2, disk_gb=1)
    monkeypatch.setattr(resources, "available_memory", lambda: 1000000000)
    with pytest.raises(RuntimeError, match="independent evidence"):
        resources.preflight(cfg, estimate, tmp_path)
    assert cfg["data"]["velocity_range_kms"] == prior
    monkeypatch.setattr(resources, "available_memory", lambda: None)
    usage = namedtuple("usage", "total used free")
    monkeypatch.setattr(resources.shutil, "disk_usage", lambda d: usage(100, 90, 10))
    with pytest.raises(RuntimeError, match="additional storage"):
        resources.preflight(cfg, estimate, tmp_path)
