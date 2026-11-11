"""Registered output files and conservative, preview-first cache cleanup."""

import json
import os
from pathlib import Path, PurePosixPath

from . import paths, records
from .delivery import atomic_text

MANAGED = ("work", "results", ".lant")


def check_roots(cfg):
    root = Path(cfg["project"]["output_dir"])
    for name in MANAGED:
        if (root / name).is_symlink():
            raise RuntimeError("managed directory must not be a symlink: %s" % (root / name))


def snapshot(cfg):
    root = Path(cfg["project"]["output_dir"])
    found = {}
    for name in MANAGED:
        for directory, dirs, files in os.walk(str(root / name), followlinks=False):
            # Include links as unknown files, but never traverse them.
            files += [d for d in dirs if (Path(directory) / d).is_symlink()]
            dirs[:] = [d for d in dirs if not (Path(directory) / d).is_symlink()]
            for file in files:
                p = Path(directory) / file
                st = p.lstat()
                found[p.relative_to(root).as_posix()] = [st.st_size, st.st_mtime_ns, st.st_ctime_ns]
    return found


def load(cfg):
    p = Path(records.rec_dir(cfg)) / "artifacts.json"
    return json.loads(p.read_text()) if p.exists() else {}


def register(cfg, before):
    """Record files created or changed by a serialized command, not pre-existing unknown files."""
    previous = load(cfg)
    after = snapshot(cfg)
    changed = {p: sig for p, sig in after.items() if before.get(p) != sig}
    previous = {p: sig for p, sig in previous.items() if p in after}
    previous.update(changed)
    previous.pop(".lant/artifacts.json", None)
    atomic_text(Path(records.rec_dir(cfg)) / "artifacts.json", json.dumps(previous, indent=1) + "\n")


def audit(cfg):
    current, owned = snapshot(cfg), load(cfg)
    unknown = sorted(p for p in current if p not in owned and p not in (".lant/artifacts.json", ".lant/lock"))
    missing = sorted(p for p in owned if p not in current)
    changed = sorted(p for p in owned if p in current and owned[p] != current[p] and p != ".lant/lock")
    return dict(unregistered=unknown, missing=missing, modified=changed)


def is_cache(relative, selection):
    p = PurePosixPath(relative)
    if p.is_absolute() or ".." in p.parts:
        return False
    prefixes = {"archive": ("work", "01_clean", "archive"), "cc": ("work", "02_ccf")}
    return any(p.parts[: len(prefixes[s])] == prefixes[s] for s in selection)


def cleanup_plan(cfg, selection):
    root = Path(cfg["project"]["output_dir"])
    owned = load(cfg)
    files, changed = [], []
    for rel, stamp in sorted(owned.items()):
        if not is_cache(rel, selection):
            continue
        p = root / rel
        if not p.exists() and not p.is_symlink():
            continue
        safe = not any(parent.is_symlink() for parent in [p] + list(p.parents) if parent != root)
        st = p.lstat()
        if not safe or not p.is_file() or [st.st_size, st.st_mtime_ns, st.st_ctime_ns] != stamp:
            changed.append(rel)
        else:
            files.append(dict(path=rel, bytes=st.st_size))
    return dict(
        files=files,
        bytes=sum(p["bytes"] for p in files),
        protected=changed,
        impact="Retained results stay readable. Rebuilding cc requires clean archives; rebuilding stack requires cc files.",
    )


def clean(cfg, selection, execute=False):
    """Delete registered, unchanged caches only after scientific products are complete."""
    from . import manifest, report
    from .delivery import write_report
    from .provenance import Lock

    check_roots(cfg)
    if not execute:
        return cleanup_plan(cfg, selection)
    if not Path(records.rec_dir(cfg)).is_dir():
        raise RuntimeError("no initialized project to clean")
    with Lock(records.rec_dir(cfg)):
        before = snapshot(cfg)
        state = manifest.status(report.resolved(cfg))
        if state["evaluate"]["state"] != "done":
            raise RuntimeError("cleanup requires a completed, current evaluation; finish the run first")
        if not Path(cfg["data"]["raw_dir"]).is_dir():
            raise RuntimeError("raw data are unavailable; caches cannot be declared rebuildable")
        plan = cleanup_plan(cfg, selection)
        if plan["protected"]:
            raise RuntimeError("modified files or links are protected: %s" % ", ".join(plan["protected"]))
        root = Path(cfg["project"]["output_dir"])
        # Record intent before deletion, so interruption cannot leave status claiming all caches exist.
        m = manifest.load(cfg)
        for step in ("clean" if s == "archive" else "cc" for s in selection):
            m["steps"][step]["cache_cleaned"] = True
        manifest._save(cfg, m)
        try:
            for item in plan["files"]:
                (root / item["path"]).unlink()
            records.log(cfg, event="cache_cleaned", selection=selection, files=len(plan["files"]), bytes=plan["bytes"])
            write_report(cfg)
        finally:
            register(cfg, before)
        return plan


def required(cfg, step):
    """Essential working products; path storage is always relative to the project."""
    table = {
        "clean": [paths.P(cfg, "01_clean", "archive_index.json"), paths.P(cfg, "01_clean", "stations.csv")],
        "cc": [paths.P(cfg, paths.CCF)],
        "stack": [paths.stack_product(cfg)],
        "reference": [paths.P(cfg, paths.REFERENCE, "reference.npz")],
        "disp": [paths.P(cfg, paths.DISP, "disp.h5")],
        "inversion": [paths.P(cfg, paths.INV, "main", "DSurfTomo.in")],
        "run": [paths.P(cfg, paths.INV, "main", "DSurfTomo.inMeasure.dat")],
        "evaluate": [paths.P(cfg, paths.INV, "evaluation", "evaluation.json")],
    }
    return [os.path.relpath(p, cfg["project"]["output_dir"]) for p in table.get(step, [])]
