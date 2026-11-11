"""Project entry points and retained scientific deliverables; no computation lives here."""

import json
import os
import shutil
import tempfile
from pathlib import Path

from . import paths, records

LAYOUT_VERSION = 1
RULES = [
    "Raw data are read-only. Use lant commands to produce managed outputs.",
    "Read REPORT.md, the step report and its evidence before continuing.",
    "The user decides scientific judgements; assistants must not confirm on their own.",
    "Do not loosen QC or alter scientific priors for more points or prettier figures.",
    "Append notes with lant note; do not create scripts or extra logs inside managed directories.",
    "Configuration changes require a recorded reason; --yes is not scientific authorization.",
]


def _umask():
    mask = os.umask(0)
    os.umask(mask)
    return mask


def atomic_text(path, text):
    """Replace a small text file in one step, with the permissions an ordinary write would give (mkstemp: 0600)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".write-", dir=str(path.parent))
    try:
        mode = path.stat().st_mode & 0o777 if path.exists() else 0o666 & ~_umask()
        os.fchmod(fd, mode)
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(text)
        os.replace(tmp, str(path))
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def initialize(cfg: dict) -> None:
    root = Path(cfg["project"]["output_dir"])
    marker = root / ".lant" / "layout.json"
    if marker.exists():
        if json.loads(marker.read_text())["version"] != LAYOUT_VERSION:
            raise RuntimeError("unsupported project layout version")
    else:
        atomic_text(marker, json.dumps({"version": LAYOUT_VERSION}) + "\n")
    readme = root / "README.md"
    if not readme.exists():
        atomic_text(
            readme,
            "# lant project\n\nStart with [REPORT.md](REPORT.md).\n\n"
            "- `results/`: retained models, dispersion, figures, validation and reproduction inputs.\n"
            "- `work/`: processing files. Only explicitly selected registered caches may be cleaned.\n"
            "- `.lant/`: state, step reports, configuration history, provenance and rotating logs.\n"
            "- `notes.md`: optional append-only notes and scientific decisions.\n\n"
            "## Operating rules\n\n" + "\n".join("- " + r for r in RULES) + "\n",
        )


def _swap(cfg, target, fill):
    """Fill a new tree beside `target`, then put it in place; the previous complete tree stays if anything fails.
    Unregistered files in the old tree are never replaced."""
    from .artifacts import load

    root = Path(cfg["project"]["output_dir"])
    target = Path(target)
    if target.is_symlink():
        raise RuntimeError("refusing to replace a symlink: %s" % target)
    owned = load(cfg)
    unknown = [
        p.relative_to(root).as_posix()
        for p in target.rglob("*")
        if p.is_file() and p.relative_to(root).as_posix() not in owned
    ]
    if unknown:
        raise RuntimeError("unregistered retained files are protected: " + ", ".join(unknown[:5]))
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".deliver-", dir=str(target.parent)))
    new, old = staging / "new", staging / "previous"
    try:
        new.mkdir()  # ordinary permissions; the mkdtemp directory itself is private
        fill(new)
        if target.exists():
            os.rename(str(target), str(old))
        try:
            os.rename(str(new), str(target))
        except BaseException:
            if old.exists():
                os.rename(str(old), str(target))
            raise
    finally:
        shutil.rmtree(str(staging))


def _copy_into(source, new):
    shutil.copytree(str(source), str(new), dirs_exist_ok=True)


def _mark(cfg, target, step):
    """Which run of which step a retained tree belongs to (manifest.record gives the run this number next)."""
    from . import manifest

    p = Path(records.rec_dir(cfg)) / "published.json"
    marks = json.loads(p.read_text()) if p.exists() else {}
    run = int(manifest.load(cfg)["steps"].get(step, {}).get("run", 0)) + 1
    marks[target] = dict(step=step, run=run, time=records._now())
    atomic_text(p, json.dumps(marks, indent=1) + "\n")


def publish(cfg, step):
    """Replace retained products only by complete ones: the dispersion table after disp; the model, its checks
    and its reproduction inputs together after evaluate (the main model has run and been evaluated).  Building
    or rerunning the inversion leaves the previous complete model in results/ until then."""
    from .inversion import run as RUN

    dest = Path(cfg["project"]["output_dir"]) / "results"
    inv = Path(paths.P(cfg, paths.INV))
    if step == "disp":
        _swap(cfg, dest / "dispersion", lambda new: _copy_into(paths.P(cfg, paths.DISP), new))
        _mark(cfg, "results/dispersion", step)
        return
    if step != "evaluate":
        return
    if RUN.status(str(inv / "main")) != "done":
        raise RuntimeError("the main model has not finished; the retained model is kept")

    def validation(new):
        for name in ("checkerboards", "controls"):
            for project in sorted((inv / name).iterdir()) if (inv / name).is_dir() else []:
                if project.is_dir() and RUN.status(str(project)) == "done":  # unrun checks are not results
                    _copy_into(project, new / name / project.name)
        _copy_into(inv / "evaluation", new / "evaluation")

    def reproduce(new):
        for p in inv.iterdir():
            if p.is_file():
                shutil.copy2(str(p), str(new / p.name))
        stations = Path(paths.P(cfg, "01_clean", "stations.csv"))
        if stations.exists():
            shutil.copy2(str(stations), str(new / "stations.csv"))
        from . import config

        # Small and useful without the raw files; relative paths as seen from its final place.
        config.dump(cfg, str(new / "effective.yaml"), base=str(dest / "reproduce"))

    for name, fill in (
        ("model", lambda new: _copy_into(inv / "main", new)),
        ("validation", validation),
        ("reproduce", reproduce),
    ):
        _swap(cfg, dest / name, fill)
        _mark(cfg, "results/" + name, step)  # each tree is labelled as soon as it is in place


def write_report(cfg: dict) -> None:
    from . import manifest, report

    root = Path(cfg["project"]["output_dir"])
    state = manifest.status(report.resolved(cfg), check_inputs=False)
    lines = [
        "# Run report",
        "",
        "State belongs to the latest recorded run; retained files may be from an earlier complete run.",
        "Use `lant status -c CONFIG --json` to recheck current inputs and configuration.",
        "",
        "| Step | State | Detail |",
        "|---|---|---|",
    ]
    for step, info in state.items():
        detail = ("Caches cleaned; rebuilding needs upstream work. " if info.get("cache_cleaned") else "") + info.get(
            "reason", ""
        )
        detail = detail.replace("|", "\\|").replace("\n", " ")
        link = ".lant/steps/%s.md" % step
        label = "[%s](%s)" % (step, link) if (root / link).exists() else step
        lines.append("| %s | %s | %s |" % (label, info["state"], detail))
    lines += ["", "## Retained deliverables", ""]
    marks_path = Path(records.rec_dir(cfg)) / "published.json"
    marks = json.loads(marks_path.read_text()) if marks_path.exists() else {}
    trees = [p for p in (root / "results").iterdir() if p.is_dir()] if (root / "results").is_dir() else []
    for tree in sorted(p for p in trees if not p.name.startswith(".")):  # not staging left by an interruption
        rel = tree.relative_to(root).as_posix()
        mark = marks.get(rel) or (dict(step="figures") if tree.name == "figures" else {})
        made = state.get(mark.get("step"), {})
        if not mark:
            origin = "origin not recorded"
        elif "run" not in mark:
            origin = "written by step %s" % mark["step"]
        else:
            origin = "from run %d of step %s (%s)" % (mark["run"], mark["step"], mark["time"])
            if made.get("run") != mark["run"]:
                origin += "; OLDER than the latest run %s of that step (%s)" % (made.get("run"), made.get("state"))
            elif made.get("state") != "done":
                origin += "; that step is now %s: earlier steps changed and a new result is not complete yet" % (
                    made.get("state", "not run")
                )
        files = sorted(p for p in tree.rglob("*") if p.is_file())
        lines.append("- [%s/](%s): %d files, %s" % (rel, rel, len(files), origin))
        shown = files if len(files) <= 25 else [p for p in files if p.parent == tree][:25]
        for p in shown:
            link = p.relative_to(root).as_posix()
            lines.append("  - [%s](%s)" % (p.relative_to(tree).as_posix(), link))
        if len(shown) < len(files):
            lines.append("  - ... %d more files in the directory" % (len(files) - len(shown)))
    lines += [
        "",
        "## Judgements and limits",
        "",
        "Read the linked step reports for key numbers, parameters, diagnostic figures and scientific judgements.",
        "Jackknife uncertainty is conditional, not a calibrated total error. Review checkerboards with the model.",
        "Directional-noise and higher-mode limitations remain those of the documented method.",
    ]
    conf = cfg.get("_confirmed") or {}
    if conf:
        lines += ["", "## Confirmation record", "", "Invocation metadata is not identity verification.", ""]
        for step, entry in conf.items():
            lines.append(
                "- %s: %s; reason: %s; invocation: %s"
                % (
                    step,
                    entry.get("_at", "unknown"),
                    entry.get("_reason", "not recorded"),
                    entry.get("_invocation", "not recorded"),
                )
            )
    atomic_text(root / "REPORT.md", "\n".join(lines) + "\n")
