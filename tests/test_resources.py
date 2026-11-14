"""Resource decisions use effective limits and leave scientific configuration untouched."""

import math
import os
from collections import namedtuple

import pytest

from lant import config, resources


def test_worker_thread_limit_does_not_change_numba_pool_ceiling():
    """A fresh compilation after a loaded runtime changes its active mask must not reload a different ceiling."""
    import os
    import subprocess
    import sys

    script = """
import os
import numba
from lant import resources
numba.set_num_threads(2)
ceiling = numba.config.NUMBA_NUM_THREADS
resources.limit_threads(1)
assert numba.get_num_threads() == 1
assert int(os.environ['NUMBA_NUM_THREADS']) == ceiling
@numba.njit
def cold_kernel(x):
    return x + 1
assert cold_kernel(2) == 3
"""
    subprocess.run(
        [sys.executable, "-c", script],
        env=dict(os.environ, NUMBA_NUM_THREADS="2"),
        check=True,
        capture_output=True,
        text=True,
        timeout=30,
    )


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


def test_inactive_page_cache_counts_as_available(tmp_path):
    proc, group = tmp_path / "proc", tmp_path / "cgroup"
    (proc / "self").mkdir(parents=True)
    (proc / "meminfo").write_text("MemAvailable: 80000000 kB\n")
    (proc / "self/cgroup").write_text("0::/job\n")
    (group / "job").mkdir(parents=True)
    (group / "job/memory.max").write_text("8000000000")
    (group / "job/memory.current").write_text("7900000000")  # mostly cache from reading the raw files
    (group / "job/memory.stat").write_text("anon 900000000\ninactive_file 6000000000\n")
    assert resources.available_memory(proc, group) == 6100000000


def test_clean_stage_shortages_are_reported_not_raised(tmp_path, monkeypatch):
    monkeypatch.setattr(resources, "available_memory", lambda: 1000000000)
    found = resources.shortages(dict(memory_per_process_gb=2, disk_gb=0), tmp_path)
    assert len(found) == 1 and "one cc process" in found[0]


def test_cores_80_percent_of_the_idle_cpus_or_as_given(monkeypatch):
    monkeypatch.setattr(resources, "usable_cpus", lambda: 40)
    monkeypatch.setattr(resources.os, "getloadavg", lambda: (10.6, 0, 0))
    assert resources.cores("auto")[0] == 23  # floor(0.8 x (40 - 10.6))
    assert resources.cores("auto", own_load=8.0)[0] == 29  # this run's own load does not count against it
    monkeypatch.setattr(resources.os, "getloadavg", lambda: (45.0, 0, 0))
    assert resources.cores("auto")[0] == 1
    assert resources.cores(3) == (3, "resources.cores = 3 (given)")


def test_usable_cpus_follow_a_cgroup_quota(tmp_path, monkeypatch):
    proc, group = tmp_path / "proc", tmp_path / "cgroup"
    (proc / "self").mkdir(parents=True)
    (proc / "self/cgroup").write_text("0::/job\n")
    (group / "job").mkdir(parents=True)
    (group / "job/cpu.max").write_text("250000 100000\n")
    monkeypatch.setattr(resources.os, "cpu_count", lambda: 40)
    monkeypatch.setattr(resources.os, "sched_getaffinity", lambda pid: set(range(16)))
    assert resources.usable_cpus(proc, group) == 2
    (group / "job/cpu.max").write_text("max 100000\n")
    assert resources.usable_cpus(proc, group) == 16


def test_own_load_decays_like_the_load_average(monkeypatch):
    clock = [1000.0]
    monkeypatch.setattr(resources.time, "monotonic", lambda: clock[0])
    own = resources.OwnLoad()
    job = own.start(8)
    clock[0] += 600
    assert own.estimate() == pytest.approx(8, rel=1e-3)
    own.stop(job)
    clock[0] += 60
    assert own.estimate() == pytest.approx(8 * math.exp(-1), rel=1e-3)


def test_plan_limits_parallel_work_to_half_the_memory(monkeypatch):
    cfg = config.defaults()
    monkeypatch.setattr(resources, "usable_cpus", lambda: 20)
    monkeypatch.setattr(resources.os, "getloadavg", lambda: (0.0, 0, 0))
    monkeypatch.setattr(resources, "available_memory", lambda: 10e9)
    p = resources.plan(cfg, 30, 1.2e9, "days")
    # 0.8 x 20 cores; 4 x 1.2 GB within half of 10 GB
    assert (p["cores"], p["parallel"], p["threads_each"], p["limited_by_memory"]) == (16, 4, 4, True)
    p = resources.plan(cfg, 3, None, "projects")
    assert (p["parallel"], p["threads_each"], p["limited_by_memory"]) == (3, 5, False)  # 16 // 3
    with pytest.raises(RuntimeError, match="needs about 6.00 GB"):
        resources.plan(cfg, 1, 6e9, "projects")


def test_step_usage_records_plans_and_worker_peaks():
    cfg = config.defaults()
    cfg["resources"]["cores"] = 2
    resources.begin_step(cfg)
    resources.plan(cfg, 5, None, "pair batches")
    resources.note_workers([os.getpid()])
    u = resources.usage()
    assert u["at_start"]["cores"] == 2 and u["plans"][0]["parallel"] == 2
    assert u["worker_peak_gb"] > 0 and u["main_peak_gb"] > 0
