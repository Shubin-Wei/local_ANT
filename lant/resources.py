"""Effective Linux memory limits and storage estimates; never adjust scientific priors."""

import os
import shutil
from pathlib import Path
from typing import Optional


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
    """What does not fit now: one cc process into 70 % of the effective memory, the full new output on disk."""
    found = []
    per_process = estimate["memory_per_process_gb"] * 1e9
    available = available_memory()
    if available is not None and per_process > 0.7 * available:
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
