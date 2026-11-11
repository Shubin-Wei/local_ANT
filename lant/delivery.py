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


def atomic_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".write-", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(text)
        os.replace(tmp, str(path))
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def initialize(cfg):
    root = Path(cfg["project"]["output_dir"])
    marker = root / ".lant" / "layout.json"
    if marker.exists():
        if json.loads(marker.read_text())["version"] != LAYOUT_VERSION:
            raise RuntimeError("unsupported project layout version")
    else:
        atomic_text(marker, json.dumps({"version": LAYOUT_VERSION}) + "\n")
    readme = root / "README.md"
    if not readme.exists():
        atomic_text(readme, "# lant project\n\nStart with [REPORT.md](REPORT.md).\n\n"
                    "- `results/`: retained models, dispersion, figures, validation and reproduction inputs.\n"
                    "- `work/`: processing files. Only explicitly selected registered caches may be cleaned.\n"
                    "- `.lant/`: state, step reports, configuration history, provenance and rotating logs.\n"
                    "- `notes.md`: optional append-only notes and scientific decisions.\n\n"
                    "## Operating rules\n\n" + "\n".join("- " + r for r in RULES) + "\n")


def _replace_tree(source, target):
    """Build before replacing; retain the old complete tree if copying or replacement fails."""
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.is_symlink():
        raise RuntimeError("refusing to replace a symlink: %s" % target)
    staging = Path(tempfile.mkdtemp(prefix=".deliver-", dir=str(target.parent)))
    new, old = staging / "new", staging / "previous"
    try:
        shutil.copytree(str(source), str(new))
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


def publish(cfg, step):
    """Retain complete scientific products before cache cleanup becomes available."""
    root = Path(cfg["project"]["output_dir"])
    dest = root / "results"
    inv = Path(paths.P(cfg, paths.INV))
    if step == "disp":
        source = Path(paths.P(cfg, paths.DISP))
        target = dest / "dispersion"
        _replace_tree(source, target)
    if step in ("inversion", "run", "evaluate"):
        if (inv / "main").exists():
            _replace_tree(inv / "main", dest / "model")
        for name in ("checkerboards", "controls", "evaluation"):
            if (inv / name).exists():
                _replace_tree(inv / name, dest / "validation" / name)
        reproduce = dest / "reproduce"
        reproduce.mkdir(parents=True, exist_ok=True)
        for p in inv.iterdir():
            if p.is_file():
                atomic_text(reproduce / p.name, p.read_text())
    # Effective configuration and station geometry are small and useful without the raw files.
    effective = Path(records.rec_dir(cfg)) / "config.yaml"
    if effective.exists():
        from . import config
        target = dest / "reproduce" / "effective.yaml"
        target.parent.mkdir(parents=True, exist_ok=True)
        config.dump(cfg, str(target) + ".tmp")
        os.replace(str(target) + ".tmp", str(target))
    stations = Path(paths.P(cfg, "01_clean", "stations.csv"))
    if stations.exists():
        atomic_text(dest / "reproduce" / "stations.csv", stations.read_text())


def write_report(cfg):
    from . import manifest, report
    root = Path(cfg["project"]["output_dir"])
    state = manifest.status(report.resolved(cfg), check_inputs=False)
    lines = ["# Run report", "", "State belongs to the latest recorded run; retained files may be from an earlier complete run.",
             "Use `lant status -c CONFIG --json` to recheck current inputs and configuration.", "",
             "| Step | State | Detail |", "|---|---|---|"]
    for step, info in state.items():
        detail = ("Caches cleaned; rebuilding needs upstream work. " if info.get("cache_cleaned") else "") + info.get("reason", "")
        detail = detail.replace("|", "\\|").replace("\n", " ")
        link = ".lant/steps/%s.md" % step
        label = "[%s](%s)" % (step, link) if (root / link).exists() else step
        lines.append("| %s | %s | %s |" % (label, info["state"], detail))
    lines += ["", "## Retained deliverables", ""]
    for p in sorted((root / "results").rglob("*")):
        if p.is_file():
            rel = p.relative_to(root).as_posix()
            lines.append("- [%s](%s)" % (rel, rel))
    lines += ["", "## Judgements and limits", "",
              "Read the linked step reports for key numbers, parameters, diagnostic figures and scientific judgements.",
              "Jackknife uncertainty is conditional, not a calibrated total error. Review checkerboards with the model.",
              "Directional-noise and higher-mode limitations remain those of the documented method."]
    conf = cfg.get("_confirmed") or {}
    if conf:
        lines += ["", "## Confirmation record", "", "Invocation metadata is not identity verification.", ""]
        for step, entry in conf.items():
            lines.append("- %s: %s; reason: %s; invocation: %s" %
                         (step, entry.get("_at", "unknown"), entry.get("_reason", "not recorded"),
                          entry.get("_invocation", "not recorded")))
    atomic_text(root / "REPORT.md", "\n".join(lines) + "\n")
