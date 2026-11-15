"""REPORT.md, the run report at the top of a project: state and next command, key numbers, items for the user to
judge, the main figures and the retained results with the run each comes from (lant.delivery publishes them)."""

import json
from pathlib import Path

import numpy as np

from . import records
from .delivery import atomic_text


def best_checkerboard(evaluation):
    """(name, median slope over depth) of the horizontal checkerboard recovered best, or None."""
    import json

    path = Path(evaluation)
    if not path.exists():
        return None
    checks = json.loads(path.read_text()).get("checkerboards") or {}
    scores = {n: float(np.median([r["slope"] for r in rows])) for n, rows in checks.items()
              if isinstance(rows, list) and rows and n != "cb_vertical"}  # fmt: skip
    return max(scores.items(), key=lambda x: x[1]) if scores else None


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
    (
        "6_checkerboard_*",
        "Checkerboard test recovered best (all tests are in results/figures; numbers in "
        "results/validation/checkerboard_recovery.csv)",
    ),
    ("6_recovery", "Checkerboard recovery with depth"),
    ("6_residuals", "Travel-time residuals of the starting and final model"),
)


def main_figures(root):
    """(figure, caption) a reader looks at first, in reading order; of the checkerboards the one recovered best
    (highest median slope over depth) stands for all."""
    folder = Path(root) / "results" / "figures"
    best = best_checkerboard(Path(root) / "results" / "validation" / "evaluation.json")
    found = []
    for name, caption in FIGURES:
        files = sorted(folder.glob(name + ".png"))
        if name == "6_checkerboard_*" and best and (folder / ("6_checkerboard_%s.png" % best[0])).exists():
            files = [folder / ("6_checkerboard_%s.png" % best[0])]
        found += [(f, caption) for f in files[:1]]
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
    inversion = summary("inversion")
    project = inversion.get("project") or {}
    periods = project.get("periods") or []
    if periods:
        out.append(
            "%s used by the inversion, %.3g-%.3g s" % (count(len(periods), "period"), min(periods), max(periods))
        )
    grid = inversion.get("grid") or {}
    if grid.get("spacing_reason"):
        out.append("grid: %s" % grid["spacing_reason"])
    if project.get("weight") is not None:
        out.append("smoothing weight %g" % project["weight"])
    checks = summary("evaluate").get("checkerboards") or {}
    if checks:
        out.append(count(len(checks), "checkerboard test"))
    for name, rows in sorted(checks.items()):
        rows = [r for r in rows if isinstance(r, dict)] if isinstance(rows, list) else []
        if rows:
            top = max(rows, key=lambda r: r["slope"])
            out.append("checkerboard %s: largest recovery %.0f %% of the anomaly at %g km (r %.2f)"
                       % (name, 100 * top["slope"], top["depth_km"], top["r"]))  # fmt: skip
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
