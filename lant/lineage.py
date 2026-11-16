"""Frozen execution records and a single parameter index for the retained results.

An editable configuration describes the next command. A record describes one
execution, including defaults and resolved values; publication points to that
record even after a failed rerun or a partial redraw. History contains metadata,
not copies of large scientific products.
"""

import copy
import datetime
import hashlib
import json
import os
import uuid
from pathlib import Path

import yaml

from . import records
from .delivery import atomic_text


def _public(cfg):
    return records.strict({k: v for k, v in cfg.items() if not k.startswith("_")})


def begin(cfg, step, options):
    from . import manifest, provenance, report

    stamp = datetime.datetime.now().strftime("%Y%m%dT%H%M%S_%f")
    name = "%s_%s_%s" % (step, stamp, uuid.uuid4().hex[:8])
    cfg["_run_record"] = ".lant/runs/%s.json" % name
    previous = manifest.load(cfg)["steps"]
    upstream = {
        s: dict(run=r["run"], record=r.get("record"), params=r.get("params", {}))
        for s, r in previous.items()
        if s in manifest._steps()[: manifest._steps().index(step)]
    }
    entry = dict(
        id=name,
        step=step,
        state="running",
        started=records._now(),
        configuration=cfg.get("_config_file"),
        requested=_public(cfg),
        configuration_sha256=cfg.get("_config_sha256"),
        params=report.used_params(cfg, step),
        options=options,
        upstream=upstream,
        code=provenance.record(records.rec_dir(cfg)),
        confirmations=copy.deepcopy(cfg.get("_confirmed", {})),
    )
    _write(cfg, entry)
    return entry


def _write(cfg, entry):
    target = Path(cfg["project"]["output_dir"]) / cfg["_run_record"]
    atomic_text(target, json.dumps(records.strict(entry), indent=1, ensure_ascii=False) + "\n")


def finish(cfg, step, state, summary=None, error=None):
    from . import report

    root = Path(cfg["project"]["output_dir"])
    entry = json.loads((root / cfg["_run_record"]).read_text())
    if entry["state"] != "running":
        return  # a completed record is never rewritten by error handling
    entry.update(
        state=state,
        finished=records._now(),
        effective=_public(cfg),
        params=report.used_params(cfg, step),
        error=error,
    )
    if summary:
        entry["figure_parameters"] = summary.get("figure_parameters", {})
        entry["failed"] = summary.get("failed", [])
    published = root / ".lant/published.json"
    marks = json.loads(published.read_text()) if published.exists() else {}
    entry["outputs"] = {p: m.get("signature") for p, m in marks.items() if m.get("record") == cfg["_run_record"]}
    source = cfg.get("_config_file")
    if source and cfg.get("_config_sha256"):
        entry["configuration_edited_during_run"] = (
            not Path(source).is_file() or hashlib.sha256(Path(source).read_bytes()).hexdigest() != cfg["_config_sha256"]
        )
    _write(cfg, entry)


def mark(cfg, relative, step, parameters=None):
    """Associate a successfully installed output with its execution, independently of other outputs."""
    mark_many(cfg, [relative], step, parameters)


def mark_many(cfg, relatives, step, parameters=None):
    """One index update even for thousands of cleaned station-day files."""
    from . import artifacts, manifest

    root = Path(cfg["project"]["output_dir"])
    path = root / ".lant/published.json"
    marks = json.loads(path.read_text()) if path.exists() else {}
    origin = dict(
        step=step,
        run=int(manifest.load(cfg)["steps"].get(step, {}).get("run", 0)) + 1,
        time=records._now(),
        record=cfg.get("_run_record"),
        parameters=parameters or {},
    )
    for relative in relatives:
        marks[relative] = dict(origin)
        target = root / relative
        if target.is_file():
            marks[relative]["signature"] = artifacts.signature(target.stat())
    atomic_text(path, json.dumps(marks, indent=1) + "\n")


def index(cfg):
    """Portable, human-readable provenance of existing published products, not the editable configuration."""
    from . import manifest

    root = Path(cfg["project"]["output_dir"])
    path = root / ".lant/published.json"
    marks = json.loads(path.read_text()) if path.exists() else {}
    products = {p: m for p, m in marks.items() if p != "results/parameters.yaml" and (root / p).exists()}
    runs = {}

    def include(relative):
        if not relative or relative in runs:
            return
        file = root / relative
        if not file.is_file():
            return
        entry = json.loads(file.read_text())
        runs[relative] = entry
        for upstream in entry.get("upstream", {}).values():
            include(upstream.get("record"))

    steps = manifest.load(cfg)["steps"]
    for m in products.values():
        include(m.get("record"))
    for r in steps.values():
        include(r.get("record"))
    return dict(
        description="Actual parameters of retained products. Edit the project config.yaml for the next command. "
        "Records preserve parameters, not historical data files. Paths are relative to the project root.",
        products=products,
        steps=steps,
        runs=runs,
    )


def publish_index(cfg):
    from .delivery import _publish_file

    contents = yaml.safe_dump(records.strict(index(cfg)), sort_keys=False, allow_unicode=True)
    _publish_file(cfg, "results/parameters.yaml", "figures", lambda path: Path(path).write_text(contents))


def parameters(cfg, target):
    """Read-only query: next settings versus the last execution, or one retained product's actual settings."""
    from . import auto, manifest, report

    if target in manifest._steps():
        current = auto.apply(copy.deepcopy(cfg))
        previous = manifest.load(cfg)["steps"].get(target, {})
        planned = report.used_params(current, target)
        actual = previous.get("params", {})
        changed = {
            k: dict(used=actual.get(k), next=planned.get(k))
            for k in sorted(set(actual) | set(planned))
            if actual.get(k) != planned.get(k)
        }
        return dict(step=target, next=planned, last_execution=previous, differences=changed)
    root = Path(cfg["project"]["output_dir"]).resolve()
    candidate = Path(target) if os.path.isabs(target) else root / target
    try:
        relative = candidate.resolve().relative_to(root).as_posix()
    except ValueError as exc:
        raise ValueError("result path must be inside the project output directory: %s" % target) from exc
    catalog = index(cfg)
    matches = [
        (p, m) for p, m in catalog["products"].items() if relative == p or relative.startswith(p.rstrip("/") + "/")
    ]
    if not matches:
        raise ValueError(
            "no recorded result for %s; use a step name or a path listed in results/parameters.yaml" % target
        )
    path, origin = max(matches, key=lambda item: len(item[0]))
    entry = catalog["runs"].get(origin.get("record"), {})
    note = None if entry else "This result predates execution records; its exact historical settings are unavailable."
    if origin.get("signature"):
        from .artifacts import signature

        if not (root / relative).is_file() or signature((root / relative).stat()) != origin["signature"]:
            note = "The file changed after publication; this record describes the earlier file, not verified current contents."
    return dict(
        product=relative,
        origin=origin,
        execution=entry,
        note=note,
    )
