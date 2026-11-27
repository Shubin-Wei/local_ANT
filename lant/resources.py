"""Cores and memory of a step (Linux).

cores: auto takes 80 % of the CPUs the machine leaves idle, max(1, floor(0.8 x (CPUs - 1-minute load average))), so
that a shared machine keeps room for others; a number in the configuration is used as given.  Work that runs in
parallel is also limited so that its estimated memory stays within half of the available memory.
"""

import math
import os

CORE_FRACTION = 0.8
MEMORY_FRACTION = 0.5


def usable_cpus():
    n = os.cpu_count() or 1
    try:
        n = min(n, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        pass
    return n


def cores(value="auto"):
    """-> (number of cores, how it was chosen)."""
    if value != "auto":
        n = int(value)
        if n < 1:
            raise ValueError("cores must be auto or a positive number")
        return n, "as configured"
    cpus = usable_cpus()
    try:
        load = os.getloadavg()[0]
    except OSError:
        load = 0.0
    n = max(1, int(math.floor(CORE_FRACTION * (cpus - load))))
    return n, "auto: 80 %% of the idle CPUs (%d CPUs, load %.1f)" % (cpus, load)


def available_memory():
    """MemAvailable of /proc/meminfo in bytes (None when unknown)."""
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    return None


def parallel(n_tasks, n_cores, bytes_per_task=None):
    """How many tasks run at once: at most the cores and the tasks, and within MEMORY_FRACTION of the memory."""
    n = max(1, min(n_tasks, n_cores))
    mem = available_memory()
    if bytes_per_task and mem:
        n = max(1, min(n, int(MEMORY_FRACTION * mem // bytes_per_task)))
    return n


def limit_threads(n=1):
    """Worker processes: BLAS libraries use n threads each (the processes share the cores).  Numba's thread count
    is set with numba.set_num_threads, not through the environment (it may not change once its threads run)."""
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = str(n)
