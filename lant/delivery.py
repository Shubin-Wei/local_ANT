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
        atomic_text(readme, PROJECT_README + "\n".join("- " + r for r in RULES) + "\n")


PROJECT_README = """# lant project

Start with [REPORT.md](REPORT.md): state of every step, key numbers, what is left for you to judge, the main
figures and links to everything below.

| Directory | Content |
|---|---|
| `results/model/` | `vs_model.csv`: the 3-D shear-velocity model |
| `results/dispersion/` | `curves.csv` (phase velocities with QC and inversion flags), `disp.h5` (all measurements and diagnostics), `median_vs_reference.txt` |
| `results/figures/` | maps, sections, dispersion, checks |
| `results/validation/` | `evaluation.json` (residuals, checkerboard recovery, control inversions), travel-time residual tables |
| `results/reproduce/` | effective configuration, station table, starting model and the solver input and output files |
| `work/` | processing files of every step; caches can be removed with `lant clean-work` |
| `.lant/` | state, step reports, configuration history, provenance and logs (written by lant only) |
| `notes.md` | your notes and decisions (`lant note "text" -c CONFIG`) |

`vs_model.csv` columns: `longitude_deg`, `latitude_deg` (WGS84), `depth_km` (below the surface of the model),
`vs_km_s`.  The values are those DSurfTomo wrote (`reproduce/DSurfTomo.inMeasure.dat`), sorted by depth, latitude
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


def publish(cfg, step):
    """Replace retained products only by complete ones: the dispersion table after disp; the model, its checks
    and its reproduction inputs together after evaluate (the main model has run and been evaluated).  Building
    or rerunning the inversion leaves the previous complete model in results/ until then."""
    from .inversion import run as RUN

    dest = Path(cfg["project"]["output_dir"]) / "results"
    inv = Path(paths.P(cfg, paths.INV))
    if step == "disp":

        def dispersion(new):  # the per-pair diagnostic figures stay in work/05_dispersion/pair_figures
            for name in ("curves.csv", "disp.h5", "median_vs_reference.txt"):
                if os.path.exists(paths.P(cfg, paths.DISP, name)):
                    shutil.copy2(paths.P(cfg, paths.DISP, name), str(new / name))

        _swap(cfg, dest / "dispersion", dispersion)
        _mark(cfg, "results/dispersion", step)
        return
    if step != "evaluate":
        return
    if RUN.status(str(inv / "main")) != "done":
        raise RuntimeError("the main model has not finished; the retained model is kept")

    def model(new):
        write_model_table(inv / "main" / "DSurfTomo.inMeasure.dat", new / "vs_model.csv")

    def validation(new):  # the numbers of the checks; their recovered models are in the figures and in work/
        for p in sorted((inv / "evaluation").iterdir()):
            shutil.copy2(str(p), str(new / p.name))

    def reproduce(new):
        for p in inv.iterdir():
            if p.is_file():
                shutil.copy2(str(p), str(new / p.name))
        for name in ("DSurfTomo.in", "surfdata.dat", "MOD", "DSurfTomo.inMeasure.dat"):
            shutil.copy2(str(inv / "main" / name), str(new / name))
        stations = Path(paths.P(cfg, "01_clean", "stations.csv"))
        if stations.exists():
            shutil.copy2(str(stations), str(new / "stations.csv"))
        from . import config

        # Small and useful without the raw files; relative paths as seen from its final place.
        config.dump(cfg, str(new / "effective.yaml"), base=str(dest / "reproduce"))

    for name, fill in (("model", model), ("validation", validation), ("reproduce", reproduce)):
        _swap(cfg, dest / name, fill)
        _mark(cfg, "results/" + name, step)  # each tree is labelled as soon as it is in place


MODEL_COLUMNS = ("longitude_deg", "latitude_deg", "depth_km", "vs_km_s")


def write_model_table(model, target):
    """DSurfTomo's "lon lat depth vs" rows with a header, sorted by depth, latitude, longitude.  The numbers are
    copied as DSurfTomo wrote them (no conversion, no rounding)."""
    rows = [line.split() for line in Path(model).read_text().splitlines() if line.strip()]
    if any(len(r) != 4 for r in rows):
        raise RuntimeError("%s: expected four columns (lon lat depth vs)" % model)
    rows.sort(key=lambda r: (float(r[2]), float(r[1]), float(r[0])))
    Path(target).write_text(",".join(MODEL_COLUMNS) + "\n" + "".join(",".join(r) + "\n" for r in rows))


def write_report(cfg: dict) -> None:
    from . import manifest, report

    root = Path(cfg["project"]["output_dir"])
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
    figures = main_figures(root)
    if figures:
        lines += ["", "## Main figures", ""]
        for f, caption in figures:
            rel = f.relative_to(root).as_posix()
            lines += ["**%s** (`%s`)" % (caption, f.name), "", "![%s](%s)" % (f.stem, rel), ""]
    lines += ["", "## Steps", "", "| Step | State | Detail |", "|---|---|---|"]
    for step, info in state.items():
        detail = ("Caches cleaned; rebuilding needs upstream work. " if info.get("cache_cleaned") else "") + info.get(
            "reason", ""
        )
        detail = detail.replace("|", "\\|").replace("\n", " ")
        link = ".lant/steps/%s.md" % step
        label = "[%s](%s)" % (step, link) if (root / link).exists() else step
        lines.append("| %s | %s | %s |" % (label, info["state"], detail))
    lines += ["", "## Retained results", ""]
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
        "## Limits",
        "",
        "Every parameter, its value and what it rests on are in the linked step reports.",
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


FIGURES = (
    ("1_stations", "Stations"),
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
    checks = summary("evaluate").get("checkerboards") or {}
    if checks:
        out.append(count(len(checks), "checkerboard test"))
    drawn = summary("figures").get("drawn") or {}
    if drawn:
        out.append(count(len(drawn), "figure"))
    return out


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
