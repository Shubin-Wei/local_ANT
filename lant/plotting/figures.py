"""The standard figure set -> <output>/figures (each figure can also be drawn alone: `lant plot <name>`).

  1_stations           station map (GMT; relief if plots.topography is set)
  2_ccf_moveout        cross-correlations against distance (1-5 Hz)
  2_ccf_stack_time     one pair: hourly correlations and SNR against stacking length
  3_dispersion_pair    one pair: waveform, phase and group images with the picks (plots.pair_example or chosen)
  3_dispersion_curves  all phase-velocity curves, mean and standard deviation
  4_start_model        average curve, 1-D model fit and the starting model
  4_kernels            depth sensitivity of the phase velocity for the starting model (fwd1d finite differences)
  5_depth_slices       Vs depth slices (GMT)
  5_sections           Vs sections: plots.sections (areal) or along the line (line)
  6_checkerboard_*     checkerboard input and recovery (GMT)
  6_recovery           recovery slope and correlation with depth
  6_residuals          travel-time residuals of the final model, rms against period and iteration
Figures of the processing steps are written by the steps themselves (01_clean/coverage.png,
03_stack/station_day_qc/station_day_qc.png, 04_notch/notch_report.png, 05_reference/fj.png and reference_pass2.png).
"""
import json
import os

from .. import paths
from . import ccf, dispersion, inversion, maps, sections

NAMES = ["stations", "ccf_moveout", "ccf_stack_time", "dispersion_pair", "dispersion_curves", "start_model", "kernels",
         "depth_slices", "sections", "checkerboards", "recovery", "residuals"]


def make_all(cfg, which=None):
    from .. import data as D
    out = paths.mk(paths.P(cfg, paths.FIG))
    stations = D.used_stations(cfg)
    stack = paths.stack_product(cfg)
    disp = paths.P(cfg, paths.DISP, "disp_final.h5")
    main = paths.P(cfg, paths.INV, "main")
    ev_path = paths.P(cfg, paths.INV, "evaluation", "evaluation.json")
    ev = json.load(open(ev_path)) if os.path.exists(ev_path) else None
    done, skipped = {}, {}
    todo = which or NAMES

    def want(n, need=()):
        if n not in todo:
            return False
        miss = [p for p in need if not os.path.exists(p)]
        if miss:
            skipped[n] = "missing %s" % ", ".join(miss)
            return False
        return True

    def attempt(n, fn):
        try:
            done[n] = fn()
        except Exception as exc:
            skipped[n] = "%s: %s" % (type(exc).__name__, exc)
            print("  figure %s not drawn: %s" % (n, exc), flush=True)

    if want("stations"):
        attempt("stations", lambda: maps.station_map(cfg, stations, os.path.join(out, "1_stations")))
    if want("ccf_moveout", [stack]):
        attempt("ccf_moveout", lambda: ccf.moveout(cfg, stack, os.path.join(out, "2_ccf_moveout")))
    pair = cfg["plots"]["pair_example"]
    if os.path.exists(disp) and not pair and ({"ccf_stack_time", "dispersion_pair"} & set(todo)):
        pair = dispersion.choose_pair(disp)
    if want("ccf_stack_time", [stack, disp]):
        bad = []
        p = paths.P(cfg, paths.STACK, "exclude_station_days.txt")
        if os.path.exists(p):
            bad = [tuple(l.split()[:2]) for l in open(p) if l.strip() and not l.startswith("#")]
        attempt("ccf_stack_time", lambda: ccf.stack_time(cfg, stack, pair, os.path.join(out, "2_ccf_stack_time"), bad_days=bad))
    if want("dispersion_pair", [stack, disp]):
        attempt("dispersion_pair", lambda: dispersion.pair_example(cfg, disp, stack, pair, os.path.join(out, "3_dispersion_pair")))
    curves = paths.P(cfg, paths.DISP, "curves")
    if want("dispersion_curves", [curves]):
        attempt("dispersion_curves", lambda: dispersion.all_curves(cfg, curves, os.path.join(out, "3_dispersion_curves")))
    if want("start_model", [paths.P(cfg, paths.INV, "vs1d.txt")]):
        attempt("start_model", lambda: inversion.start_model(cfg, paths.P(cfg, paths.INV), os.path.join(out, "4_start_model")))
    if want("kernels", [paths.P(cfg, paths.INV, "vs1d.txt")]):
        attempt("kernels", lambda: inversion.kernels(cfg, paths.P(cfg, paths.INV), os.path.join(out, "4_kernels")))
    model = os.path.join(main, "DSurfTomo.inMeasure.dat")
    if want("depth_slices", [model]):
        attempt("depth_slices", lambda: maps.depth_slices(cfg, model, stations, os.path.join(out, "5_depth_slices")))
    if want("sections", [model]):
        if cfg["project"]["array_type"] == "line":
            attempt("sections", lambda: sections.line_section(cfg, main, stations, os.path.join(out, "5_sections")))
        elif cfg["plots"]["sections"]:
            attempt("sections", lambda: sections.surface_sections(cfg, main, stations, os.path.join(out, "5_sections")))
        else:
            skipped["sections"] = "plots.sections is empty (give start and end points of the sections)"
    if "checkerboards" in todo:
        root = paths.P(cfg, paths.INV, "checkerboards")
        vert = {c["name"] for c in cfg["checks"]["vertical"]} if isinstance(cfg["checks"]["vertical"], list) else set()
        vert |= {n for n in (os.listdir(root) if os.path.isdir(root) else []) if n.startswith("cb_vertical")}
        for name in sorted(os.listdir(root)) if os.path.isdir(root) else []:
            d = os.path.join(root, name)
            if os.path.exists(os.path.join(d, "DSurfTomo.inSyn.dat")):
                attempt("checkerboard_" + name, lambda d=d, name=name: maps.checkerboard(
                    cfg, d, stations, os.path.join(out, "6_checkerboard_" + name), vertical=name in vert))
    if ev and "recovery" in todo:
        attempt("recovery", lambda: inversion.recovery(cfg, ev, os.path.join(out, "6_recovery")))
    res_csv = paths.P(cfg, paths.INV, "evaluation", "residuals_final.csv")
    if ev and want("residuals", [res_csv]):
        attempt("residuals", lambda: inversion.residuals(cfg, res_csv, ev, os.path.join(out, "6_residuals")))
    return dict(drawn=done, not_drawn=skipped)
