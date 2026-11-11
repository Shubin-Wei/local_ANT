"""Registered output files and conservative, preview-first cache cleanup."""

import contextlib
import fnmatch
import json
import os
from pathlib import Path, PurePosixPath

from . import paths, records
from .delivery import atomic_text

MANAGED = ("work", "results", ".lant")
INTERNAL = (".lant/artifacts.json", ".lant/lock", ".lant/pending.json")  # bookkeeping of the writer itself
# Cache files by what the steps write; anything else in these directories is never cleaned.
CACHES = {
    "archive": (("work", "01_clean", "archive"), ("*.mseed", "*.done.json")),
    "cc": (("work", "02_ccf"), ("day_*.h5",)),
}


def check_roots(cfg):
    root = Path(cfg["project"]["output_dir"])
    for name in MANAGED:
        if (root / name).is_symlink():
            raise RuntimeError("managed directory must not be a symlink: %s" % (root / name))
        for directory, dirs, files in os.walk(str(root / name), followlinks=False):
            for entry in dirs + files:
                candidate = Path(directory) / entry
                if candidate.is_symlink():
                    raise RuntimeError("managed paths must not be symlinks: %s" % candidate)


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
                try:
                    st = p.lstat()
                except FileNotFoundError:  # renamed or removed by a writer after the listing
                    continue
                found[p.relative_to(root).as_posix()] = signature(st)
    return found


def signature(st):
    """Size and modification time: kept by archiving copies (cp -a, rsync -a), changed by edits.  Not a content
    hash, and not proof of who wrote a file."""
    return [st.st_size, st.st_mtime_ns]


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
    for p in INTERNAL:
        previous.pop(p, None)
    atomic_text(Path(records.rec_dir(cfg)) / "artifacts.json", json.dumps(previous, indent=1) + "\n")


def audit(cfg: dict) -> dict:
    current, owned = snapshot(cfg), load(cfg)
    unknown = sorted(p for p in current if p not in owned and p not in INTERNAL)
    missing = sorted(p for p in owned if p not in current)
    changed = sorted(p for p in owned if p in current and owned[p] != current[p] and p not in INTERNAL)
    return dict(unregistered=unknown, missing=missing, modified=changed)


def _pending(cfg):
    return Path(records.rec_dir(cfg)) / "pending.json"


def interrupted(cfg):
    """Read-only: the step of a write command that ended without its cleanup, before the next one recovers it."""
    try:
        return json.loads(_pending(cfg).read_text()).get("step") or "unknown"
    except (OSError, ValueError):
        return None


def recover(cfg):
    """A write command that ended without its cleanup (killed, out of memory, lost terminal) left pending.json.
    Holding the lock proves it is gone: register what it wrote, as its normal end would have, and record its step
    as failed, so that status does not keep showing the step's earlier run.  -> the step, or None."""
    p = _pending(cfg)
    if not p.exists():
        return None
    try:
        pending = json.loads(p.read_text())
    except ValueError:  # written atomically, so this is not expected; adopt nothing rather than guess
        pending = {}
    p.unlink()
    register(cfg, pending.get("before", snapshot(cfg)))
    step = pending.get("step")
    from . import manifest

    if step in manifest._steps():
        manifest.record(
            cfg,
            step,
            "failed",
            error="interrupted: the command (pid %s) ended without finishing; "
            "its partial outputs were kept and are replaced by a new run" % pending.get("pid"),
        )
    records.log(cfg, event="interrupted", step=step, pid=pending.get("pid"), started=pending.get("started"))
    return step


@contextlib.contextmanager
def writing(cfg, what):
    """The one way a command writes to a project: the project lock, recovery of an interrupted earlier command,
    and the registration of every file this command creates or changes (also when it fails)."""
    from .provenance import Lock

    check_roots(cfg)
    rec = records.rec_dir(cfg, create=True)
    with Lock(rec, what) as lock:
        if lock.nested:
            yield
            return
        before = snapshot(cfg)  # taken first, so that the records written by the recovery are registered too
        recovered = recover(cfg)
        if recovered:
            import sys

            print(
                "[lant] the earlier %s command was interrupted; its outputs are registered and the step is marked "
                "failed" % recovered,
                file=sys.stderr,
                flush=True,
            )
        atomic_text(_pending(cfg), json.dumps(dict(step=what, pid=os.getpid(), started=records._now(), before=before)))
        try:
            yield
        finally:
            _pending(cfg).unlink()
            register(cfg, before)


def protect_step_outputs(cfg, step):
    """Do not let name-based exporters or figure cleanup overwrite an unregistered file."""
    directories = {
        "clean": "01_clean",
        "cc": paths.CCF,
        "stack": paths.STACK,
        "reference": paths.REFERENCE,
        "disp": paths.DISP,
        "inversion": paths.INV,
        "run": paths.INV,
        "evaluate": paths.INV,
        "figures": paths.FIG,
    }
    prefix = os.path.relpath(paths.P(cfg, directories[step]), cfg["project"]["output_dir"]) + "/"
    conflicts = [p for p in audit(cfg)["unregistered"] if p.startswith(prefix)]
    if conflicts:
        raise RuntimeError(
            "unregistered files in this step's output directory are protected; relocate them before "
            "rerunning: " + ", ".join(conflicts[:5])
        )


def is_cache(relative, selection):
    p = PurePosixPath(relative)
    if p.is_absolute() or ".." in p.parts:
        return False
    return any(
        p.parts[: len(CACHES[s][0])] == CACHES[s][0]
        and any(fnmatch.fnmatch(p.name, pattern) for pattern in CACHES[s][1])
        for s in selection
    )


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
        if not safe or not p.is_file() or signature(st) != stamp:
            changed.append(rel)
        else:
            files.append(dict(path=rel, bytes=st.st_size))
    return dict(
        files=files,
        bytes=sum(p["bytes"] for p in files),
        protected=changed,
        impact="Retained results stay readable. Rebuilding cc requires clean archives; rebuilding stack, and redrawing "
        "the stack-time figures (2_ccf_stack_time_*), requires cc files.",
    )


def clean(cfg: dict, selection, execute: bool = False) -> dict:
    """Delete registered, unchanged caches only after scientific products are complete."""
    from . import manifest, report
    from .delivery import write_report

    check_roots(cfg)
    if not execute:
        return cleanup_plan(cfg, selection)
    if not Path(records.rec_dir(cfg)).is_dir():
        raise RuntimeError("no initialized project to clean")
    with writing(cfg, "clean-work"):
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
        for item in plan["files"]:
            (root / item["path"]).unlink()
        records.log(cfg, event="cache_cleaned", selection=selection, files=len(plan["files"]), bytes=plan["bytes"])
        write_report(cfg)
        return plan


def required(cfg, step):
    """Essential working products; path storage is always relative to the project."""
    table = {
        "clean": [paths.P(cfg, "01_clean", "archive_index.json"), paths.P(cfg, "01_clean", "stations.csv")],
        "cc": [paths.P(cfg, paths.CCF)],
        "stack": [paths.stack_product(cfg)],
        "reference": [paths.P(cfg, paths.REFERENCE, "reference.npz")],
        "disp": [paths.P(cfg, paths.DISP, "disp.h5"), paths.P(cfg, paths.DISP, "curves.csv")],
        "inversion": [paths.P(cfg, paths.INV, "main", "DSurfTomo.in")],
        "run": [paths.P(cfg, paths.INV, "main", "DSurfTomo.inMeasure.dat")],
        "evaluate": [paths.P(cfg, paths.INV, "evaluation", "evaluation.json")],
    }
    selected = table.get(step, [])
    if step == "cc":
        selected += [str(p) for p in Path(paths.P(cfg, paths.CCF)).glob("day_*.h5")]
    if step == "clean":
        selected += [str(p) for p in Path(paths.P(cfg, "01_clean", "archive")).rglob("*.mseed")]
    return [os.path.relpath(p, cfg["project"]["output_dir"]) for p in selected]
