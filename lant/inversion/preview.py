"""Publish a completed main model before slower recovery experiments finish."""

import copy
import json
import shutil
from pathlib import Path

from .. import config, data, delivery, records
from . import run


def publish(cfg, project):
    if run.status(project) != "done":
        raise RuntimeError("cannot publish an unfinished main model")
    target = Path(cfg["project"]["output_dir"]) / "results/main_preview"

    def fill(new):
        from ..plotting import figures, sections

        delivery.write_model_table(Path(project) / "DSurfTomo.inMeasure.dat", new / "vs_model.csv")
        for name in ("DSurfTomo.in", "surfdata.dat", "MOD", "DSurfTomo.inMeasure.dat"):
            shutil.copy2(Path(project) / name, new / name)
        config.dump(cfg, str(new / "effective.yaml"), base=str(target))
        note = dict(
            status="main_model_complete_checks_pending",
            scientific_acceptance=False,
            meaning="Main inversion complete; recovery/control evaluation may still be pending.",
        )
        try:
            drawing = copy.deepcopy(cfg)
            stations = data.used_stations(cfg)
            figures.resolve(drawing, str(new / "DSurfTomo.inMeasure.dat"), stations)
            draw = sections.line_section if cfg["project"]["array_type"] == "line" else sections.surface_sections
            draw(drawing, project, stations, str(new / "sections"))
            note["figure_status"] = "complete"
        except Exception as exc:
            # A missing plotting dependency must not discard a valid solver result or stop other projects.
            note["figure_status"] = "failed"
            note["figure_error"] = "%s: %s" % (type(exc).__name__, exc)
            records.log(cfg, event="main_preview_figure_failed", error=note["figure_error"])
        (new / "status.json").write_text(json.dumps(note, indent=2) + "\n")

    with records.phase(cfg, "main_preview"):
        delivery._swap(cfg, target, fill)
        delivery._mark(cfg, "results/main_preview", "run")
    records.log(cfg, event="main_preview_published", path=str(target))
    print("  Main model preview (checks may still be pending): %s" % target, flush=True)
