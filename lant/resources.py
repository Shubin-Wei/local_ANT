"""Computing resources of a step: cores, memory and disk (Linux); never adjust scientific priors.

Nobody picks the number of cores: resources.cores = auto (the default) takes CORE_FRACTION (80 %) of what the machine
leaves idle, max(1, floor(0.8 x (usable CPUs - 1-minute load average))), where the usable CPUs are the CPU count limited by the
affinity mask and a cgroup CPU quota.  An integer in the configuration is used as given.  The number is taken again
at the start of every step and of every solver project; work already running keeps its threads.  Work that runs in
parallel is limited further so that its estimated peak memory stays within MEMORY_FRACTION of the available memory;
the rest waits in the queue.
"""

import math
import os
import shutil
import time
from pathlib import Path
from typing import Optional

# Share of the effective available memory that the estimated peak of parallel work may take.
MEMORY_FRACTION = 0.5
# Share of the idle CPUs that auto takes (user's choice, 2026-11-13; was 0.5): the rest is left for others on a
# shared machine and for the 1-minute load average lagging behind work that has just started.
CORE_FRACTION = 0.8
# Time constant of the Linux 1-minute load average (s): an exponential average updated every 5 s.
LOAD_TAU_S = 60.0


def available_memory(proc: Path = Path("/proc"), cgroup: Path = Path("/sys/fs/cgroup")) -> Optional[int]:
    limits = []
    try:
        fields = dict(line.split(":", 1) for line in (proc / "meminfo").read_text().splitlines())
        limits.append(int(fields["MemAvailable"].split()[0]) * 1024)
    except (OSError, KeyError, ValueError):
        try:
            limits.append(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_AVPHYS_PAGES"))
        except (ValueError, OSError):
            pass
    # Include limits on every ancestor: a job may live below the cgroup mount root.
    roots = [cgroup]
    try:
        for line in (proc / "self/cgroup").read_text().splitlines():
            hierarchy, controllers, relative = line.split(":", 2)
            if hierarchy == "0" or "memory" in controllers.split(","):
                base = cgroup if hierarchy == "0" else cgroup / "memory"
                path = base / relative.lstrip("/")
                if ".." not in path.parts:
                    roots += [path] + [p for p in path.parents if p == base or base in p.parents]
    except (OSError, ValueError):
        pass
    for root in roots:
        for maximum, usage, cache in (
            ("memory.max", "memory.current", "inactive_file"),
            ("memory.limit_in_bytes", "memory.usage_in_bytes", "total_inactive_file"),
        ):
            try:
                cap = int((root / maximum).read_text())
                used = int((root / usage).read_text())
            except (OSError, ValueError):
                continue
            # The usage counts the page cache; its inactive part is reclaimed on demand (as MemAvailable counts it).
            limits.append(max(0, cap - max(0, used - _stat(root / "memory.stat", cache))))
    return min(limits) if limits else None


def _stat(path, key):
    try:
        for line in path.read_text().splitlines():
            name, _, value = line.partition(" ")
            if name == key:
                return int(value)
    except (OSError, ValueError):
        pass
    return 0


def check_disk(directory, required_bytes: float) -> int:
    free = shutil.disk_usage(directory).free
    if required_bytes > free:
        raise RuntimeError(
            "estimated additional storage %.2f GB exceeds %.2f GB free on %s; "
            "free space or choose another output filesystem" % (required_bytes / 1e9, free / 1e9, directory)
        )
    return free


def shortages(estimate, directory, archive_bytes=0):
    """What does not fit now: one cc process into MEMORY_FRACTION of the effective memory, the full new output on
    disk."""
    found = []
    per_process = estimate["memory_per_process_gb"] * 1e9
    available = available_memory()
    if available is not None and per_process > MEMORY_FRACTION * available:
        found.append(
            "one cc process needs about %.2f GB; effective available memory is %.2f GB "
            "(including readable cgroup limits). Free memory or use a larger worker. "
            "Scientific velocity priors must only change with independent evidence."
            % (per_process / 1e9, available / 1e9)
        )
    # Full new output is a conservative bound and includes the transient replacement of existing files.
    try:
        check_disk(directory, estimate["disk_gb"] * 1e9 + archive_bytes)
    except RuntimeError as exc:
        found.append(str(exc))
    return found


def preflight(cfg, estimate, directory, archive_bytes=0):
    """Before cc allocates its workers: refuse what cannot fit."""
    found = shortages(estimate, directory, archive_bytes)
    if found:
        raise RuntimeError("; ".join(found))


# ------------------------------------------------------------------------------------------------ cores
def usable_cpus(proc: Path = Path("/proc"), cgroup: Path = Path("/sys/fs/cgroup")) -> int:
    """os.cpu_count() limited by the affinity mask and by a cgroup CPU quota (v2 cpu.max, v1 cfs quota)."""
    n = os.cpu_count() or 1
    try:
        n = min(n, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        pass
    roots = [cgroup]
    try:
        for line in (proc / "self/cgroup").read_text().splitlines():
            hierarchy, controllers, relative = line.split(":", 2)
            if hierarchy == "0" or "cpu" in controllers.split(","):
                base = cgroup if hierarchy == "0" else cgroup / "cpu"
                path = base / relative.lstrip("/")
                if ".." not in path.parts:
                    roots += [path] + [p for p in path.parents if p == base or base in p.parents]
    except (OSError, ValueError):
        pass
    for root in roots:
        quota = None
        try:
            fields = (root / "cpu.max").read_text().split()
            if fields[0] != "max":
                quota = int(fields[0]) / int(fields[1])
        except (OSError, ValueError, IndexError, ZeroDivisionError):
            try:
                q, period = int((root / "cpu.cfs_quota_us").read_text()), int((root / "cpu.cfs_period_us").read_text())
                if q > 0 and period > 0:
                    quota = q / period
            except (OSError, ValueError):
                pass
        if quota:
            n = min(n, max(1, math.floor(quota)))
    return max(1, n)


def cores(setting="auto", own_load=0.0):
    """-> (cores, how).  setting: resources.cores (auto or an integer).  own_load: the part of the load average
    that is this run's own work (OwnLoad.estimate), which does not count against it."""
    if setting != "auto":
        return int(setting), "resources.cores = %d (given)" % int(setting)
    total = usable_cpus()
    load = os.getloadavg()[0]
    others = max(0.0, load - own_load)
    n = max(1, math.floor(CORE_FRACTION * (total - others)))
    own = ", without %.1f of its own" % own_load if own_load > 0.05 else ""
    return n, "auto: %g x (%d usable CPUs - 1-minute load %.1f%s)" % (CORE_FRACTION, total, others, own)


class OwnLoad:
    """This run's own contribution to the 1-minute load average, so that work started after other work of the same
    step is not counted as somebody else's load.  A job of n busy threads from t0 to t1 adds
    n (exp(-(now - t1) / tau) - exp(-(now - t0) / tau)) (the exponential average of the kernel)."""

    def __init__(self):
        self.jobs = []

    def start(self, threads):
        job = [threads, time.monotonic(), None]
        self.jobs.append(job)
        return job

    @staticmethod
    def stop(job):
        job[2] = time.monotonic()

    def estimate(self):
        now = time.monotonic()
        return sum(
            n * (math.exp(-(now - (end if end is not None else now)) / LOAD_TAU_S) - math.exp(-(now - t0) / LOAD_TAU_S))
            for n, t0, end in self.jobs
        )


def plan(cfg, tasks, per_task_bytes=None, what="processes"):
    """Cores of this step (taken now) and how many of `tasks` run at the same time: processes = min(tasks, cores,
    as many as fit into MEMORY_FRACTION of the available memory), threads = cores // processes.  per_task_bytes:
    estimated peak memory of one task (None: not estimated).  Raises when one task alone does not fit.
    -> dict, also kept for the step report (usage())."""
    n, how = cores(cfg["resources"]["cores"])
    available = available_memory()
    fit = None
    if per_task_bytes and available:
        fit = math.floor(MEMORY_FRACTION * available / per_task_bytes)
        if fit < 1:
            raise RuntimeError(
                "one of the %s needs about %.2f GB (estimate); %.0f %% of the %.2f GB of effectively available "
                "memory is %.2f GB: free memory or use a machine with more"
                % (
                    what,
                    per_task_bytes / 1e9,
                    100 * MEMORY_FRACTION,
                    available / 1e9,
                    MEMORY_FRACTION * available / 1e9,
                )
            )
    parallel = max(1, min(int(tasks), n, fit or n))
    out = dict(
        what=what,
        cores=n,
        cores_rule=how,
        tasks=int(tasks),
        parallel=parallel,
        threads_each=max(1, n // parallel),
        memory_available_gb=round(available / 1e9, 3) if available else None,
        memory_fraction=MEMORY_FRACTION,
        memory_each_estimate_gb=round(per_task_bytes / 1e9, 3) if per_task_bytes else None,
        limited_by_memory=bool(fit is not None and fit < min(int(tasks), n)),
    )
    _usage["plans"].append(out)
    return out


def thread_environment(threads):
    """Environment of a child process that may use `threads` threads (OpenMP, BLAS, Numba)."""
    t = str(int(threads))
    return dict(OMP_NUM_THREADS=t, OPENBLAS_NUM_THREADS=t, MKL_NUM_THREADS=t, NUMBA_NUM_THREADS=t)


def limit_threads(threads):
    """Initializer of pool workers: their OpenMP/BLAS/Numba threads.  The variables reach libraries that start in the
    worker; loaded Numba keeps its import-time pool ceiling and only its active thread mask changes."""
    import sys

    environment = thread_environment(threads)
    numba = sys.modules.get("numba")
    if numba is not None:
        # The compiler reloads environment settings even for a serial kernel. Changing this ceiling after
        # set_num_threads launches the pool fails on a cold compile; the active mask is the runtime control.
        environment["NUMBA_NUM_THREADS"] = str(numba.config.NUMBA_NUM_THREADS)
    os.environ.update(environment)
    try:
        if numba is not None:
            numba.set_num_threads(min(int(threads), numba.config.NUMBA_NUM_THREADS))
    except (ValueError, AttributeError):
        pass


# ------------------------------------------------------------------------------------------------ usage of a step
_usage = dict(plans=[], workers_peak_kb=0, solver=[], cores=None)


def _status_kb(pid, key):
    try:
        for line in open("/proc/%s/status" % pid):
            if line.startswith(key + ":"):
                return int(line.split()[1])
    except (OSError, ValueError, IndexError):
        pass
    return None


def begin_step(cfg):
    """At the start of a step: its cores, a new record of what it uses, and the peak resident memory of this process
    reset (Linux clear_refs 5)."""
    n, how = cores(cfg["resources"]["cores"])
    _usage.update(plans=[], workers_peak_kb=0, solver=[], cores=dict(cores=n, rule=how))
    try:
        with open("/proc/self/clear_refs", "w") as fh:
            fh.write("5")
    except OSError:
        pass


def note_workers(pids):
    """Peak resident memory of worker processes that are still alive (read before a pool closes)."""
    for pid in pids:
        kb = _status_kb(pid, "VmHWM")
        if kb:
            _usage["workers_peak_kb"] = max(_usage["workers_peak_kb"], kb)


def note_solver(record):
    _usage["solver"].append(record)


def usage():
    """What the step used: the plans (cores, parallel tasks, threads, memory estimate) and the measured peaks."""
    main = _status_kb("self", "VmHWM")
    out = dict(
        at_start=_usage["cores"], plans=list(_usage["plans"]), main_peak_gb=round(main / 1e6, 3) if main else None
    )
    if _usage["workers_peak_kb"]:
        out["worker_peak_gb"] = round(_usage["workers_peak_kb"] / 1e6, 3)
    if _usage["solver"]:
        out["solver_projects"] = list(_usage["solver"])
    return out
