"""Project entry points and retained scientific deliverables; no computation lives here."""

import json
import os
import shutil
import tempfile
from pathlib import Path

from . import paths, records

# 2: lant 2.0 (Python inversion, automatic resources); earlier projects are not continued (old_project_message).
LAYOUT_VERSION = 2
RULES = [
    "Raw data are read-only. Use lant commands to produce managed outputs.",
    "Read REPORT.md, the step report and its evidence before continuing.",
    "The user decides scientific judgements; assistants must not confirm on their own.",
    "Do not loosen QC or alter scientific priors for more points or prettier figures.",
    "Write no files of your own into this project except notes (lant note) and temporary files under scratch/; "
    "lant status lists every other file it did not write.",
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


def old_project_message(root, version):
    """What to do with an output directory written by an earlier lant (no compatibility)."""
    return (
        "%s was written by an earlier version of lant (project layout %s; this version writes layout %d) and is not "
        "continued: run the chain again with this version into a new output directory (set project.output_dir in the "
        "configuration to a new directory, then lant run all -c CONFIG); the old directory can stay for comparison"
        % (root, version if version is not None else "1.4 or older", LAYOUT_VERSION)
    )


def initialize(cfg: dict) -> None:
    root = Path(cfg["project"]["output_dir"])
    marker = root / ".lant" / "layout.json"
    if marker.exists():
        version = json.loads(marker.read_text())["version"]
        if version != LAYOUT_VERSION:
            raise RuntimeError(old_project_message(root, version))
    elif (root / ".lant" / "manifest.json").exists():  # a project of lant 1.4 or older (no layout record)
        raise RuntimeError(old_project_message(root, None))
    else:
        atomic_text(marker, json.dumps({"version": LAYOUT_VERSION}) + "\n")
    readme = root / "README.md"
    if not readme.exists():
        atomic_text(readme, PROJECT_README + "\n".join("- " + r for r in RULES) + "\n")
    from .assistant_rules import write

    write(cfg, root)


PROJECT_README = """# lant project

Start with [REPORT.md](REPORT.md): state of every step, key numbers, what is left for you to judge, the main
figures and links to everything below.

| Directory | Content |
|---|---|
| `results/model/` | `vs_model.csv`: the 3-D shear-velocity model |
| `results/dispersion/` | `curves.csv` (phase velocities with QC and inversion flags); all measurements with their diagnostics are in `work/05_dispersion/` |
| `results/figures/` | maps, sections, dispersion, checks |
| `results/parameters.yaml` | the configuration with every value resolved, and the record of your confirmations (latest completed step) |
| `results/validation/` | `evaluation.json` (residuals, checkerboard recovery, control inversions), travel-time residual tables |
| `results/main_preview/` | main model, section and parameters available during the run, before checks finish; removed when the evaluated model is published |
| `work/` | processing files of every step, kept so that nothing has to be computed again; the solver files to rerun the inversion are in `work/06_inversion/main/` (`DSurfTomo.in`, `surfdata.dat`, `MOD`); caches can be removed with `lant clean-work` |
| `.lant/` | state, step reports, configuration history, provenance and logs (written by lant only) |
| `notes.md` | your notes and decisions (`lant note "text" -c CONFIG`) |
| `scratch/` | temporary files of people and assistants (optional; not checked, safe to delete) |
| `AGENTS.md`, `CLAUDE.md` or `GEMINI.md` | the rules for AI assistants working in this directory: one file, for the assistant chosen with project.ai_assistant (AGENTS.md when none is chosen) |

`vs_model.csv` columns: `longitude_deg`, `latitude_deg` (WGS84), `depth_km` (below the surface of the model),
`vs_km_s`.  The values are those DSurfTomo wrote (`work/06_inversion/main/DSurfTomo.inMeasure.dat`), sorted by depth, latitude
and longitude.  Nodes far from the ray paths keep the starting model: read the model with the checkerboard figures.

## Operating rules

"""


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


def _mark(cfg, target, step):
    """Which run of which step a retained tree belongs to (manifest.record gives the run this number next)."""
    from . import manifest

    p = Path(records.rec_dir(cfg)) / "published.json"
    marks = json.loads(p.read_text()) if p.exists() else {}
    run = int(manifest.load(cfg)["steps"].get(step, {}).get("run", 0)) + 1
    marks[target] = dict(step=step, run=run, time=records._now())
    atomic_text(p, json.dumps(marks, indent=1) + "\n")


STEP_FIGURES = {
    "clean": (("01_clean", "coverage.png", "0_data_coverage.png"),),
    "stack": ((paths.STACK, "station_qc/station_qc.png", "2_station_qc.png"),),
    "reference": (
        (paths.REFERENCE, "fj.png", "3_fj_phase_velocity.png"),
        (paths.REFERENCE, "zones.png", "3_reference_zones.png"),
    ),
}


def _publish_file(cfg, relative, step, fill):
    """Atomically replace one owned file; a failed copy leaves the previous result intact."""
    from .artifacts import load

    target = Path(cfg["project"]["output_dir"]) / relative
    if target.is_symlink() or (target.exists() and relative not in load(cfg)):
        raise RuntimeError("unregistered retained files are protected: " + relative)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".deliver-", dir=str(target.parent))
    try:
        with os.fdopen(fd, "wb") as fh:
            os.fchmod(fh.fileno(), target.stat().st_mode & 0o777 if target.exists() else 0o666 & ~_umask())
        fill(temporary)
        os.replace(temporary, target)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    _mark(cfg, relative, step)


def publish(cfg, step):
    """Publish completed step figures and parameters, dispersion after disp, and evaluated model products.
    Rebuilding an inversion keeps the previous evaluated model until a new evaluation is complete."""
    from . import config

    dest = Path(cfg["project"]["output_dir"]) / "results"
    inv = Path(paths.P(cfg, paths.INV))
    if step == "disp":
        # the measurements with their diagnostics (disp.h5, median_vs_reference.txt) and the per-pair figures stay in
        # work/05_dispersion
        def dispersion(new):
            shutil.copy2(paths.P(cfg, paths.DISP, "curves.csv"), str(new / "curves.csv"))

        _swap(cfg, dest / "dispersion", dispersion)
        _mark(cfg, "results/dispersion", step)
    elif step == "evaluate":
        _publish_evaluation(cfg, dest, inv)
    for folder, name, target in STEP_FIGURES.get(step, ()):
        if step == "stack" and not cfg["stack"]["station_qc"]["enabled"]:
            continue
        if step == "reference" and (
            cfg["reference"]["source"] != "fj" or (name == "zones.png" and cfg["reference"]["zone_km"] is None)
        ):
            continue
        source = Path(paths.P(cfg, folder, name))
        if source.is_file():  # moved: the figure is kept once, in results/figures
            _publish_file(cfg, "results/figures/" + target, step, lambda new: os.replace(str(source), new))
    _publish_file(cfg, "results/parameters.yaml", step, lambda new: config.dump(cfg, new, base=str(dest)))
    _retire(cfg, dest / "parameters")  # per-step copies of earlier versions; each step's values are in the manifest


def _publish_evaluation(cfg, dest, inv):
    from .inversion import run as RUN

    step = "evaluate"
    if RUN.status(str(inv / "main")) != "done":
        raise RuntimeError("the main model has not finished; the retained model is kept")

    def model(new):
        write_model_table(inv / "main" / "DSurfTomo.inMeasure.dat", new / "vs_model.csv")

    def validation(new):  # the numbers of the checks; their recovered models are in the figures and in work/
        for p in sorted((inv / "evaluation").iterdir()):
            shutil.copy2(str(p), str(new / p.name))

    for name, fill in (("model", model), ("validation", validation)):
        _swap(cfg, dest / name, fill)
        _mark(cfg, "results/" + name, step)  # each tree is labelled as soon as it is in place
    _retire(cfg, dest / "main_preview")  # superseded by the evaluated model
    _retire(cfg, dest / "reproduce")  # copies of work/ written by earlier versions


def _retire(cfg, target):
    """Remove a retained tree that a newer one supersedes, only when every file in it is one lant wrote and has not
    been changed since; otherwise it stays as it is."""
    from .artifacts import load, signature

    root = Path(cfg["project"]["output_dir"])
    target = Path(target)
    if not target.is_dir() or target.is_symlink():
        return False
    owned = load(cfg)
    files = [p for p in target.rglob("*") if not p.is_dir() or p.is_symlink()]
    if any(p.is_symlink() or owned.get(p.relative_to(root).as_posix()) != signature(p.lstat()) for p in files):
        return False
    shutil.rmtree(str(target))
    marks = Path(records.rec_dir(cfg)) / "published.json"
    if marks.exists():
        m = json.loads(marks.read_text())
        rel = target.relative_to(root).as_posix()
        kept = {k: v for k, v in m.items() if k != rel and not k.startswith(rel + "/")}
        if kept != m:
            atomic_text(marks, json.dumps(kept, indent=1) + "\n")
    return True


MODEL_COLUMNS = ("longitude_deg", "latitude_deg", "depth_km", "vs_km_s")


def write_model_table(model, target):
    """DSurfTomo's "lon lat depth vs" rows with a header, sorted by depth, latitude, longitude.  The numbers are
    copied as DSurfTomo wrote them (no conversion, no rounding)."""
    rows = [line.split() for line in Path(model).read_text().splitlines() if line.strip()]
    if any(len(r) != 4 for r in rows):
        raise RuntimeError("%s: expected four columns (lon lat depth vs)" % model)
    rows.sort(key=lambda r: (float(r[2]), float(r[1]), float(r[0])))
    Path(target).write_text(",".join(MODEL_COLUMNS) + "\n" + "".join(",".join(r) + "\n" for r in rows))


def _origin(mark, state):
    if not mark:
        return "origin not recorded"
    if "run" not in mark:
        return "written by step %s" % mark["step"]
    made = state.get(mark["step"], {})
    origin = "from run %d of step %s (%s)" % (mark["run"], mark["step"], mark["time"])
    if made.get("run") != mark["run"]:
        origin += "; OLDER than the latest run %s of that step (%s)" % (made.get("run"), made.get("state"))
    elif made.get("state") != "done":
        origin += "; that step is now %s; read its report before use" % made.get("state", "not run")
    return origin


def write_report(cfg: dict) -> None:
    from . import manifest, report

    root = Path(cfg["project"]["output_dir"])
    marks_path = Path(records.rec_dir(cfg)) / "published.json"
    marks = json.loads(marks_path.read_text()) if marks_path.exists() else {}
    resolved = report.resolved(cfg)
    state = manifest.status(resolved, check_inputs=False)
    command, why = manifest.next_action(resolved, state)
    lines = ["# Run report", "", "**%s.**" % (why[0].upper() + why[1:].rstrip(".")) if why else ""]
    if command:
        who = " (your decision: read the items below first)" if command.startswith("lant confirm") else ""
        lines += ["", "Next%s: `%s`" % (who, command)]
    lines += ["", "*The state of the latest recorded run; `lant status -c CONFIG` also rechecks the raw inputs.*"]
    numbers = key_numbers(cfg)
    if numbers:
        lines += ["", "## Key numbers", ""] + ["- " + n for n in numbers]
    recorded = manifest.load(cfg)["steps"]
    waiting = [
        (s, recorded.get(s, {}).get("to_confirm") or []) for s, info in state.items() if info["state"] == "to confirm"
    ]
    if waiting:
        lines += ["", "## For you to judge", ""]
        for step, items in waiting:
            lines += [
                "Step **%s** (details and figures: [.lant/steps/%s.md](.lant/steps/%s.md)):" % (step, step, step),
                "",
            ]
            lines += ["- " + item for item in items] + [""]
        lines.append('Decide, then record it: `lant confirm STEP -c CONFIG --reason "your decision and its evidence"`.')
    if (root / "results/parameters.yaml").is_file():
        rel = "results/parameters.yaml"
        lines += [
            "",
            "## Parameters",
            "",
            "- [%s](%s): every value resolved, with your confirmations; %s. The values each step used are in its "
            "report (.lant/steps/)." % (rel, rel, _origin(marks.get(rel), state)),
        ]
    figures = main_figures(root)
    if figures:
        lines += ["", "## Main figures", ""]
        for f, caption in figures:
            rel = f.relative_to(root).as_posix()
            lines += ["**%s** (`%s`)" % (caption, f.name), "", "![%s](%s)" % (f.stem, rel), ""]
            if rel in marks:
                lines += [_origin(marks[rel], state), ""]
    lines += ["", "## Steps", "", "| Step | State | Time | Detail |", "|---|---|---|---|"]
    seconds = step_seconds(cfg)
    for step, info in state.items():
        detail = ("Caches cleaned; rebuilding needs upstream work. " if info.get("cache_cleaned") else "") + info.get(
            "reason", ""
        )
        detail = detail.replace("|", "\\|").replace("\n", " ")
        link = ".lant/steps/%s.md" % step
        label = "[%s](%s)" % (step, link) if (root / link).exists() else step
        when = _duration(seconds.get(step))  # of its latest completed run, also when it is out of date now
        lines.append("| %s | %s | %s | %s |" % (label, info["state"], when, detail))
    lines += ["", "## Retained results", ""]
    trees = [p for p in (root / "results").iterdir() if p.is_dir()] if (root / "results").is_dir() else []
    for tree in sorted(p for p in trees if not p.name.startswith(".")):  # not staging left by an interruption
        rel = tree.relative_to(root).as_posix()
        origin = "the main ones are shown above" if tree.name == "figures" else _origin(marks.get(rel), state)
        files = sorted(p for p in tree.rglob("*") if p.is_file())
        lines.append("- [%s/](%s): %d files, %s" % (rel, rel, len(files), origin))
        if tree.name == "figures":
            continue  # the main ones are shown above; the folder holds them all
        shown = files if len(files) <= 25 else [p for p in files if p.parent == tree][:25]
        for p in shown:
            link = p.relative_to(root).as_posix()
            lines.append("  - [%s](%s)" % (p.relative_to(tree).as_posix(), link))
        if len(shown) < len(files):
            lines.append("  - ... %d more files in the directory" % (len(files) - len(shown)))
    lines += [
        "",
        "## Limits",
        "",
        "Every parameter, its value and what it rests on are in the linked step reports.",
        "Jackknife uncertainty is conditional, not a calibrated total error. Review checkerboards with the model.",
        "Directional-noise and higher-mode limitations remain those of the documented method.",
    ]
    conf = cfg.get("_confirmed") or {}
    if conf:
        lines += ["", "## Confirmation record", ""]
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


FIGURES = (
    ("0_data_coverage", "Data coverage after cleaning"),
    ("1_stations", "Stations"),
    ("2_station_qc", "Station diagnostics (candidates require your decision)"),
    ("3_fj_phase_velocity", "Array F-J phase-velocity spectrum and reference picks"),
    ("3_reference_zones", "Local reference curves"),
    ("3_dispersion_curves", "Phase-velocity curves used by the inversion"),
    ("4_start_model", "1-D starting model and its fit"),
    ("5_depth_slices", "Vs at depth"),
    ("5_sections", "Vs sections"),
    ("6_checkerboard_*", "Checkerboard test (the first; all are in results/figures)"),
    ("6_recovery", "Checkerboard recovery with depth"),
    ("6_residuals", "Travel-time residuals of the starting and final model"),
)


def main_figures(root):
    """(figure, caption) a reader looks at first, in reading order (the first checkerboard stands for all)."""
    folder = Path(root) / "results" / "figures"
    found = []
    for name, caption in FIGURES:
        found += [(f, caption) for f in sorted(folder.glob(name + ".png"))[:1]]
    return found


def key_numbers(cfg):
    """A few numbers that tell what the run produced, from the step summaries (absent steps are left out)."""

    def summary(step):
        return records.read_summary(cfg, step) or {}

    def count(n, word):
        return "%d %s%s" % (n, word, "" if n == 1 else "s")

    out = []
    clean = summary("clean")
    if clean.get("stations"):
        out.append("%s, %s" % (count(clean["stations"], "station"), count(clean.get("days", 0), "processing day")))
    export = summary("disp").get("export") or {}
    if "phase_curves" in export:
        out.append(
            "%s (%s)" % (count(export["phase_curves"], "phase-velocity curve"), count(export["phase_points"], "point"))
        )
    residuals = summary("evaluate").get("residuals") or {}
    start = next((r for r in residuals.get("iterations") or [] if r.get("model") == "start"), None)
    final = residuals.get("final") or {}
    if final.get("rms_ms") is not None:
        out.append(
            "residual rms of all %d travel times %s%.1f ms"
            % (final.get("n", 0), "%.1f -> " % start["rms_ms"] if start else "", final["rms_ms"])
        )
    project = summary("inversion").get("project") or {}
    if project.get("weight") is not None:
        out.append("smoothing weight %g, damping %g" % (project["weight"], project["damp"]))
    checks = summary("evaluate").get("checkerboards") or {}
    if checks:
        out.append(count(len(checks), "checkerboard test"))
    folder = Path(cfg["project"]["output_dir"]) / "results" / "figures"
    shown = sorted(folder.glob("*.png")) if folder.is_dir() else []
    if shown:
        out.append("%s in results/figures" % count(len(shown), "figure"))
    return out


def step_seconds(cfg):
    """Wall time (s) of the latest completed run of every step, from .lant/run_log.jsonl."""
    path = Path(records.rec_dir(cfg)) / "run_log.jsonl"
    out = {}
    if path.exists():
        for line in path.read_text().splitlines():
            try:
                e = json.loads(line)
            except ValueError:
                continue
            if e.get("event") == "done" and e.get("step") and e.get("seconds") is not None:
                out[e["step"]] = float(e["seconds"])
    return out


def _duration(seconds):
    if seconds is None:
        return ""
    if seconds < 90:
        return "%.0f s" % seconds
    if seconds < 5400:
        return "%.0f min" % (seconds / 60)
    return "%.1f h" % (seconds / 3600)


def finish_lines(cfg):
    """What a finished chain tells the user: where to look, and what it produced."""
    root = Path(cfg["project"]["output_dir"])
    lines = ["  report:  %s" % (root / "REPORT.md"), "  results: %s/" % (root / "results")]
    numbers = key_numbers(cfg)
    if numbers:
        lines.append("  " + "; ".join(numbers))
    return lines


def append_note(cfg, line):
    """Append one line to notes.md under its own small lock, not the project lock: notes can be written while a
    step runs.  Appending keeps every earlier line; notes.md is the user's file, not a managed output."""
    import fcntl

    root = Path(cfg["project"]["output_dir"])
    rec = Path(records.rec_dir(cfg, create=True))
    file = root / "notes.md"
    with open(str(rec / "notes.lock"), "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        new = not file.exists()
        with open(str(file), "a", encoding="utf-8") as fh:
            fh.write(("# Notes and decisions\n" if new else "") + "\n- %s\n" % line)
    return file
