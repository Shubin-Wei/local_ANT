"""Effective Linux memory limits and storage estimates; never adjust scientific priors."""

import os
import shutil
from pathlib import Path


def available_memory(proc=Path("/proc"), cgroup=Path("/sys/fs/cgroup")):
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
        for maximum, usage in (("memory.max", "memory.current"), ("memory.limit_in_bytes", "memory.usage_in_bytes")):
            try:
                cap = int((root / maximum).read_text())
                used = int((root / usage).read_text())
                limits.append(max(0, cap - used))
            except (OSError, ValueError):
                pass
    return min(limits) if limits else None


def check_disk(directory, required_bytes):
    free = shutil.disk_usage(directory).free
    if required_bytes > free:
        raise RuntimeError(
            "estimated additional storage %.2f GB exceeds %.2f GB free on %s; "
            "free space or choose another output filesystem" % (required_bytes / 1e9, free / 1e9, directory)
        )
    return free


def preflight(cfg, estimate, directory, archive_bytes=0):
    per_process = estimate["memory_per_process_gb"] * 1e9
    available = available_memory()
    if available is not None and per_process > 0.7 * available:
        raise RuntimeError(
            "one cc process needs about %.2f GB; effective available memory is %.2f GB "
            "(including readable cgroup limits). Free memory or use a larger worker. "
            "Scientific velocity priors must only change with independent evidence."
            % (per_process / 1e9, available / 1e9)
        )
    # Full new output is a conservative bound and includes the transient replacement of existing files.
    check_disk(directory, estimate["disk_gb"] * 1e9 + archive_bytes)
    return available
