"""The standard figure set -> <output>/results/figures (each can be drawn alone: `lant plot <name>`).

  1_stations             station map over coloured relief, excluded stations as red crosses (GMT)
  2_ccf_moveout_<band>   cross-correlations against distance, one figure per band (the measurement band cut at
                         round frequencies); distances shown up to the last distance bin with enough pairs
  2_ccf_stack_time_<k>   three pairs (short, middle, long): hourly correlations and SNR against stacking length
  3_dispersion_pair_<k>  the same pairs: waveform, phase and group images with the accepted picks
  3_dispersion_curves    all phase-velocity curves, mean and standard deviation
  4_start_model          average curve, 1-D model fit and the starting model
  4_kernels              illustrative 1-D sensitivity (surfdisp96 finite differences; not the inversion Jacobian)
  5_depth_slices         areal arrays: Vs depth slices over grey relief, with the section lines (GMT)
  5_slices_3d            areal arrays: five slices in a perspective view, the relief above them (GMT)
  5_sections             Vs sections: plots.sections (areal) or along the line (line) (GMT)
  6_checkerboard_*       checkerboard input and recovery: maps (areal), sections along the line (line) (GMT)
  6_recovery             recovery slope and correlation with depth
  6_residuals            travel-time residuals of the starting and final model
  6_path_coverage        retained and input source-receiver connections at representative periods
Figures of the processing steps are written by the steps themselves (work/01_clean/coverage.png,
work/03_stack/station_qc/station_qc.png, work/04_reference/fj.png) and moved by delivery to results/figures
as 0_data_coverage, 2_station_qc and 3_fj_phase_velocity when the corresponding step completes.
"""

import json
import os
import time

from .. import paths
from . import ccf, coverage, dispersion, inversion, maps, sections
from . import gmt as G

NAMES = [
    "stations",
    "ccf_moveout",
    "ccf_stack_time",
    "dispersion_pair",
    "dispersion_curves",
    "start_model",
    "kernels",
    "depth_slices",
    "slices_3d",
    "sections",
    "checkerboards",
    "recovery",
    "residuals",
    "path_coverage",
]
# Groups with numbered or banded copies: discard obsolete copies only after successful rendering.
MULTI = {
    "ccf_moveout": "2_ccf_moveout*",
    "ccf_stack_time": "2_ccf_stack_time*",
    "dispersion_pair": "3_dispersion_pair*",
    "checkerboards": "6_checkerboard_*",
    "depth_slices": "5_depth_slices*",
    "slices_3d": "5_slices_3d*",
}


def resolve(cfg, model, stations):
    """plots.vs_range / slice_depths_km / contour_step given as "auto": from the final model (nodes within twice the
    median nearest-station distance of a station, updated depths)."""
    import numpy as np

    from ..inversion.run import read_model
    from ..pipeline import recovery_radius
    from ..signal import local_xy_km, nice_number

    p = cfg["plots"]
    notes = {}
    if not os.path.exists(model):
        return notes
    lon, lat, dep, v = read_model(model)
    if p["slice_depths_km"] == "auto":  # the model file holds the updated depths only
        upd = [float(z) for z in dep]
        k = np.unique(np.round(np.linspace(0, len(upd) - 1, min(8, len(upd)))).astype(int))
        p["slice_depths_km"] = [upd[i] for i in k]
        notes["slice_depths_km"] = "up to 8 of the updated depth nodes, evenly chosen"
    if p["vs_range"] == "auto":
        LO, LA = np.meshgrid(lon, lat, indexing="ij")
        sl = np.array([stations[n]["lon"] for n in stations])
        sb = np.array([stations[n]["lat"] for n in stations])
        lat0, lon0 = float(np.mean(sb)), float(np.mean(sl))
        gx, gy = local_xy_km(LO.ravel(), LA.ravel(), lat0, lon0)
        sx, sy = local_xy_km(sl, sb, lat0, lon0)
        from scipy.spatial import cKDTree

        near = cKDTree(np.c_[sx, sy]).query(np.c_[gx, gy])[0] <= recovery_radius(cfg)
        vals = v.reshape(-1, len(dep))[near].ravel()
        vals = vals[np.isfinite(vals)]
        lo, hi = np.percentile(vals, [2, 98])
        step = nice_number((hi - lo) / 10, "down") if hi > lo else 0.1
        p["vs_range"] = [round(float(np.floor(lo / step) * step), 6), round(float(np.ceil(hi / step) * step), 6)]
        notes["vs_range"] = "2nd-98th percentile %.3f-%.3f km/s of the model near the stations" % (lo, hi)
    if p["contour_step"] == "auto":
        lo, hi = p["vs_range"]
        p["contour_step"] = nice_number((hi - lo) / 10, "nearest")
        notes["contour_step"] = "about a tenth of the colour range"
    return notes


class _Figures:
    """What every group of figures needs: the configuration, the input files, the record of what was drawn."""

    def __init__(self, cfg, which):
        from .. import data as D

        self.cfg, self.todo = cfg, which or NAMES
        self.out = paths.mk(paths.P(cfg, paths.FIG))
        self.stations = D.used_stations(cfg)
        self.stack = paths.stack_product(cfg)
        self.disp = paths.P(cfg, paths.DISP, "disp.h5")
        self.main = paths.P(cfg, paths.INV, "main")
        self.model = os.path.join(self.main, "DSurfTomo.inMeasure.dat")
        ev_path = paths.P(cfg, paths.INV, "evaluation", "evaluation.json")
        self.ev = json.load(open(ev_path)) if os.path.exists(ev_path) else None
        self.line = cfg["project"]["array_type"] == "line"
        self.done, self.skipped, self.failed = {}, {}, set()
        self.seconds = {}
        self.figure_parameters = {}
        self.excluded = {}
        try:
            allst = D.read_canonical(D.canonical_table(cfg["project"]["output_dir"]))
            self.excluded = {k: v for k, v in allst.items() if k not in self.stations}
        except Exception:
            pass

    def path(self, name):
        return os.path.join(self.out, name)

    def want(self, n, need=()):
        if n not in self.todo:
            return False
        miss = [p for p in need if not os.path.exists(p)]
        if miss:
            self.skipped[n] = "missing %s" % ", ".join(miss)
            self.failed.add(n)  # an unavailable input must not retire the previous complete figures
            return False
        return True

    def attempt(self, n, fn, parameters=None):
        import tempfile

        destination = self.out
        started = time.perf_counter()
        try:
            with tempfile.TemporaryDirectory(prefix=".draw-", dir=destination) as staging:
                self.out = staging
                result = fn()  # file list, dict(files=..., numbers...), or None: nothing to draw for these data
                produced = result.get("files", []) if isinstance(result, dict) else (result or [])
                if not produced:
                    self.skipped[n] = "nothing to draw"
                for path in produced:
                    if not os.path.isfile(path) or not os.path.getsize(path):
                        raise RuntimeError("figure output is missing or empty: %s" % path)
                    if os.path.commonpath([staging, os.path.abspath(path)]) != staging:
                        raise RuntimeError("figure output escaped its staging directory")
                targets = []
                for path in produced:
                    target = os.path.join(destination, os.path.basename(path))
                    os.replace(path, target)
                    targets.append(target)
                if targets:
                    self.done[n] = dict(result, files=targets) if isinstance(result, dict) else targets
                    if hasattr(self, "figure_parameters"):
                        actual = dict(parameters or {})
                        if isinstance(result, dict):
                            actual.update(result.get("parameters", {}))
                        self.figure_parameters[n] = dict(plots=self.cfg["plots"].copy(), **actual)
            from ..progress import detail

            detail("%d figure tasks completed" % len(self.done))
        except Exception as exc:
            self.skipped[n] = "%s: %s" % (type(exc).__name__, exc)
            self.failed.add(n)
            print("  figure %s not drawn: %s" % (n, exc), flush=True)
        finally:
            self.out = destination
            # Some low-level callers construct _Figures without __init__.
            if hasattr(self, "seconds"):
                self.seconds[n] = time.perf_counter() - started

    def remove_old(self):
        """Remove obsolete variants after the group was drawn again.  A group with a failed figure keeps its previous
        ones; a variant left out for the present data (a band without measurements) is obsolete and removed."""
        import glob

        drawn = {path for r in self.done.values() for path in (r["files"] if isinstance(r, dict) else r)}
        for name, pattern in MULTI.items():
            failed = any(
                key == name
                or key.startswith(name)
                or key.startswith(pattern.rstrip("*"))
                or (name == "checkerboards" and key.startswith("checkerboard_"))
                for key in self.failed
            )
            if name in self.todo and not failed:
                for file in glob.glob(os.path.join(self.out, pattern)):
                    if file not in drawn:
                        os.remove(file)


def _station_map(F):
    if F.want("stations"):
        metadata = os.path.join(F.main, "project.json")
        count = json.load(open(metadata)).get("stations_with_data") if os.path.isfile(metadata) else None
        F.attempt(
            "stations",
            lambda: maps.station_map(
                F.cfg, F.stations, F.path("1_stations"), excluded=F.excluded, inversion_count=count
            ),
        )


def _moveout(F):
    """One moveout figure per band; distance and velocity window from the measurements in the band."""
    if not F.want("ccf_moveout", [F.stack]):
        return
    band, _ = ccf.figure_band(F.cfg)
    explicit = F.cfg["plots"]["ccf_bands_hz"]
    lag = F.cfg["plots"]["ccf_maxlag_s"]
    for lo, hi in explicit or ccf.figure_bands(*band):
        if not explicit:
            lo, hi = ccf.clamp_default_band(F.stack, [lo, hi])
        nm = "2_ccf_moveout_%.12g-%.12gHz" % (lo, hi)
        win = ccf.band_window(F.disp, (lo, hi)) if os.path.exists(F.disp) else None
        if os.path.exists(F.disp) and not win and not explicit:
            F.skipped[nm] = "no phase velocity was measured in this band"
            continue
        kw = dict(dmax=win["dmax"], vref=((win["vfast"], "#ff1a1a"), (win["vslow"], "#00b4ff"))) if win else {}
        F.attempt(
            nm,
            lambda lo=lo, hi=hi, nm=nm, kw=kw: ccf.moveout(F.cfg, F.stack, F.path(nm), band=(lo, hi), maxlag=lag, **kw),
            parameters=dict(band_hz=[lo, hi], maxlag_s=lag),
        )


def _example_pairs(F):
    """The same example pairs (plots.pair_example, or short / middle / long) for the stacking and dispersion figures."""
    requested = {"ccf_stack_time", "dispersion_pair"} & set(F.todo)
    available = [F.want(name, [F.stack, F.disp]) for name in sorted(requested)]
    if not any(available):
        return
    chosen = dispersion.choose_pairs(F.disp)
    given = F.cfg["plots"]["pair_example"]
    if given:
        given = [given] if isinstance(given, str) else list(given)
        bands = dict(dispersion.choose_pairs(F.disp, n=10**6, min_points=1))
        chosen = [(p, bands.get(p)) for p in given]

    def vwin(band):
        win = ccf.band_window(F.disp, band, min_samples=1) if band else None
        return (win["vfast"], win["vslow"]) if win else None

    if F.want("ccf_stack_time", [F.stack, F.disp]):
        bad = []
        p = paths.P(F.cfg, paths.STACK, "removed.txt")
        if os.path.exists(p):
            bad = [tuple(l.split()) for l in open(p) if len(l.split()) == 2 and not l.startswith("#")]
        for k, (pair, band) in enumerate(chosen, 1):
            band = F.cfg["plots"]["ccf_stack_band_hz"] or ccf.clamp_default_band(
                F.stack, band or ccf.figure_band(F.cfg)[0]
            )
            F.attempt(
                "ccf_stack_time_%d" % k,
                lambda pair=pair, band=band, k=k, vw=vwin(band): ccf.stack_time(
                    F.cfg,
                    F.stack,
                    pair,
                    F.path("2_ccf_stack_time_%d" % k),
                    band=band,
                    bad_days=bad,
                    vwin=vw,
                    maxlag=F.cfg["plots"]["ccf_maxlag_s"],
                ),
                parameters=dict(pair=pair, band_hz=band),
            )
    if F.want("dispersion_pair", [F.stack, F.disp]):
        for k, (pair, band) in enumerate(chosen, 1):
            F.attempt(
                "dispersion_pair_%d" % k,
                lambda pair=pair, k=k, vw=vwin(band): dispersion.pair_example(
                    F.cfg, F.disp, F.stack, pair, F.path("3_dispersion_pair_%d" % k), vwin=vw
                ),
                parameters=dict(pair=pair),
            )


def _curves_and_start(F):
    curves = paths.P(F.cfg, paths.DISP, "curves.csv")
    if F.want("dispersion_curves", [curves]):
        F.attempt("dispersion_curves", lambda: dispersion.all_curves(F.cfg, curves, F.path("3_dispersion_curves")))
    vs1d = paths.P(F.cfg, paths.INV, "vs1d.txt")
    if F.want("start_model", [vs1d]):
        F.attempt(
            "start_model", lambda: inversion.start_model(F.cfg, paths.P(F.cfg, paths.INV), F.path("4_start_model"))
        )
    if F.want("kernels", [vs1d]):
        F.attempt("kernels", lambda: inversion.kernels(F.cfg, paths.P(F.cfg, paths.INV), F.path("4_kernels")))


def _model(F):
    """Depth slices and the 3-D view (areal arrays), sections (both)."""
    cfg = F.cfg
    if F.line:
        for n in ("depth_slices", "slices_3d"):
            if n in F.todo:
                F.skipped[n] = "line array: the model is a narrow strip along the line, drawn as the section 5_sections"
    else:
        traces = (
            [(s["name"], s["start"], s["end"]) for s in cfg["plots"]["sections"]] if cfg["plots"]["sections"] else None
        )
        if F.want("depth_slices", [F.model]):
            F.attempt(
                "depth_slices",
                lambda: maps.depth_slices(cfg, F.model, F.stations, F.path("5_depth_slices"), traces=traces),
            )
        if F.want("slices_3d", [F.model]):
            F.attempt("slices_3d", lambda: maps.slices_3d(cfg, F.model, F.stations, F.path("5_slices_3d")))
    if F.want("sections", [F.model]):
        if F.line:
            F.attempt("sections", lambda: sections.line_section(cfg, F.main, F.stations, F.path("5_sections")))
        elif cfg["plots"]["sections"]:
            F.attempt("sections", lambda: sections.surface_sections(cfg, F.main, F.stations, F.path("5_sections")))
        else:
            F.skipped["sections"] = "plots.sections is empty (give start and end points of the sections)"


def _checkerboards(F):
    """Every checkerboard that has run: maps (areal; vertical tests in their own layout), sections (line)."""
    if "checkerboards" not in F.todo:
        return
    cfg = F.cfg
    from ..pipeline import check_dirs

    root = paths.P(cfg, paths.INV, "checkerboards")
    # the checkerboards of the present configuration (as evaluated), not every directory that happens to be there
    names = sorted(os.path.basename(d) for d in check_dirs(cfg) if os.path.dirname(d) == root)
    vert = {c["name"] for c in cfg["checks"]["vertical"]} if isinstance(cfg["checks"]["vertical"], list) else set()
    vert |= {n for n in names if n.startswith("cb_vertical")}
    for name in names:
        d = os.path.join(root, name)
        if not os.path.exists(os.path.join(d, "DSurfTomo.inSyn.dat")):
            continue
        if F.line:
            F.attempt(
                "checkerboard_" + name,
                lambda d=d, name=name: sections.line_checkerboard(
                    cfg, d, F.stations, F.path("6_checkerboard_" + name), label=name
                ),
            )
        else:
            F.attempt(
                "checkerboard_" + name,
                lambda d=d, name=name: maps.checkerboard(
                    cfg, d, F.stations, F.path("6_checkerboard_" + name), vertical=name in vert
                ),
            )


def _evaluation(F):
    if F.ev and "recovery" in F.todo:
        F.attempt("recovery", lambda: inversion.recovery(F.cfg, F.ev, F.path("6_recovery")))
    res_csv = paths.P(F.cfg, paths.INV, "evaluation", "residuals_final.csv")
    if F.ev and F.want("residuals", [res_csv]):
        F.attempt("residuals", lambda: inversion.residuals(F.cfg, res_csv, F.ev, F.path("6_residuals")))
    if F.want("path_coverage", [res_csv]):
        F.attempt("path_coverage", lambda: coverage.path_coverage(F.cfg, res_csv, F.path("6_path_coverage")))


GROUPS = [_station_map, _moveout, _example_pairs, _curves_and_start, _model, _checkerboards, _evaluation]


def make_all(cfg, which=None):
    """Draw the figures `which` (names of NAMES; all by default) -> dict(drawn, not_drawn, failed, plot_parameters).
    failed: figures whose drawing raised ("name: reason"); the step is then incomplete (lant.pipeline.run_step),
    their previous images are kept.  A figure with nothing to draw is not_drawn only."""
    F = _Figures(cfg, which)
    if not hasattr(F, "figure_parameters"):
        F.figure_parameters = {}
    unknown = set(F.todo) - set(NAMES)
    if unknown:
        raise ValueError("unknown figure name(s): %s; available: %s" % (", ".join(sorted(unknown)), ", ".join(NAMES)))
    notes = resolve(cfg, F.model, F.stations)
    G.RELIEF_NOTES.clear()
    from .figure_batch import render

    results = render(cfg, F.todo) if len(F.todo) > 1 else None
    if results is None:
        for group in GROUPS:
            group(F)
    else:
        for result in results:
            F.done.update(result["drawn"])
            F.skipped.update(result["not_drawn"])
            F.failed.update(result["failed"])
            F.seconds.update(result["figure_seconds"])
            F.figure_parameters.update(result.get("figure_parameters", {}))
            G.RELIEF_NOTES.extend(result["relief_notes"])
    F.remove_old()
    if cfg.get("_run_record"):
        from ..lineage import mark

        for name, result in F.done.items():
            for file in result.get("files", []) if isinstance(result, dict) else result:
                mark(cfg, os.path.relpath(file, cfg["project"]["output_dir"]), "figures", F.figure_parameters.get(name))
    if G.RELIEF_NOTES:
        F.skipped["relief"] = "; ".join(dict.fromkeys(G.RELIEF_NOTES))
    return dict(
        drawn=F.done,
        not_drawn=F.skipped,
        failed=["%s: %s" % (n, F.skipped.get(n, "failed")) for n in sorted(F.failed)],
        figure_seconds=F.seconds,
        plot_parameters={k: dict(value=cfg["plots"][k], rule=v) for k, v in notes.items()},
        figure_parameters=F.figure_parameters,
    )
