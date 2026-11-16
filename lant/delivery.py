"""Project entry points and retained scientific deliverables; no computation lives here (REPORT.md: run_report)."""

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
| `results/parameters.yaml` | actual parameters and source execution records of retained data and individual figures; edit config.yaml for the next command |
| `results/validation/` | `evaluation.json` (residuals, checkerboard recovery, control inversions), travel-time residual tables |
| `results/main_preview/` | main model, section and parameters available during the run, before checks finish; removed when the evaluated model is published |
| `work/` | processing files of every step, kept so that nothing has to be computed again; the solver files to rerun the inversion are in `work/06_inversion/main/` (`DSurfTomo.in`, `surfdata.dat`, `MOD`); caches can be removed with `lant clean-work` |
| `.lant/` | state, step reports, configuration history, provenance and logs (written by lant only) |
| `notes.md` | your notes and decisions (`lant note "text" -c CONFIG`) |
| `scratch/` | temporary files of people and assistants (optional; not checked, safe to delete) |
| `AGENTS.md`, `CLAUDE.md` or `GEMINI.md` | the rules for AI assistants working in this directory: one file, for the assistant chosen with project.ai_assistant (AGENTS.md when none is chosen) |

`vs_model.csv` columns: `longitude_deg`, `latitude_deg` (WGS84), `depth_km` (below the surface of the model),
`vs_km_s`, `dws` (summed data sensitivity at the node in the last iteration: 0 where no ray passes, empty where the
inversion does not update the node).  `results/validation/checkerboard_recovery.csv`: per checkerboard and depth the
recovered fraction of the anomaly (slope) and the pattern correlation r.  The values are those DSurfTomo wrote (`work/06_inversion/main/DSurfTomo.inMeasure.dat`), sorted by depth, latitude
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
    from .lineage import mark

    mark(cfg, target, step)


STEP_FIGURES = {
    "clean": (("01_clean", "coverage.png", "0_data_coverage.png"),),
    "stack": ((paths.STACK, "station_qc/station_qc.png", "2_station_qc.png"),),
    "reference": (
        (paths.REFERENCE, "fj.png", "3_fj_phase_velocity.png"),
        (paths.REFERENCE, "zones.png", "3_reference_zones.png"),
    ),
    "inversion": ((paths.INV, "start_model.png", "4_start_model.png"),),
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
    if relative != "results/parameters.yaml":
        _mark(cfg, relative, step)


def publish(cfg, step):
    """Publish completed step figures and parameters, dispersion after disp, and evaluated model products.
    Rebuilding an inversion keeps the previous evaluated model until a new evaluation is complete."""
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
    if not cfg.get("_run_record"):  # direct callers without a pipeline record
        from .lineage import publish_index

        publish_index(cfg)
    _retire(cfg, dest / "parameters")  # per-step copies of earlier versions; each step's values are in the manifest


def _publish_evaluation(cfg, dest, inv):
    from .inversion import run as RUN

    step = "evaluate"
    if RUN.status(str(inv / "main")) != "done":
        raise RuntimeError("the main model has not finished; the retained model is kept")

    def model(new):
        write_model_table(inv / "main" / "DSurfTomo.inMeasure.dat", new / "vs_model.csv", inv / "main" / "DWS.dat")

    def validation(new):  # the numbers of the checks; their recovered models are in the figures and in work/
        for p in sorted((inv / "evaluation").iterdir()):
            shutil.copy2(str(p), str(new / p.name))
        write_recovery_table(inv / "evaluation" / "evaluation.json", new / "checkerboard_recovery.csv")

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


def write_model_table(model, target, dws=None):
    """The solver's "lon lat depth vs" rows with a header, sorted by depth, latitude, longitude; the numbers are copied
    as written (no conversion, no rounding).  dws: the solver's DWS.dat (lon lat depth DWS of the updated nodes, last
    iteration): a column `dws`, the summed sensitivity of the kept data at the node (0: no ray; empty: a node the
    inversion does not update, the outer ring and the deepest depth).  Read the model where the data reach, together
    with the checkerboards (Rawlinson & Spakman 2016)."""
    rows = [line.split() for line in Path(model).read_text().splitlines() if line.strip()]
    if any(len(r) != 4 for r in rows):
        raise RuntimeError("%s: expected four columns (lon lat depth vs)" % model)
    rows.sort(key=lambda r: (float(r[2]), float(r[1]), float(r[0])))
    columns = list(MODEL_COLUMNS)
    if dws is not None and Path(dws).exists():
        key = lambda a, b, c: (round(float(a), 4), round(float(b), 4), round(float(c), 4))  # noqa: E731
        sens = {}
        for line in Path(dws).read_text().splitlines():
            if line.strip() and not line.startswith("#"):
                a, b, c, v = line.split()
                sens[key(a, b, c)] = "%.4g" % float(v)
        rows = [r + [sens.get(key(*r[:3]), "")] for r in rows]
        columns.append("dws")
    Path(target).write_text(",".join(columns) + "\n" + "".join(",".join(r) + "\n" for r in rows))


def write_recovery_table(evaluation, target):
    """checkerboard_recovery.csv: test, depth_km, slope (fraction of the anomaly recovered), r (pattern
    correlation), nodes, from the evaluation (nodes with data sensitivity)."""
    import json

    path = Path(evaluation)
    checks = json.loads(path.read_text()).get("checkerboards") or {} if path.exists() else {}
    lines = ["test,depth_km,slope,r,nodes"]
    for name, rows in sorted(checks.items()):
        for r in rows if isinstance(rows, list) else []:
            lines.append("%s,%g,%.3f,%.3f,%d" % (name, r["depth_km"], r["slope"], r["r"], r["n"]))
    Path(target).write_text("\n".join(lines) + "\n")


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
