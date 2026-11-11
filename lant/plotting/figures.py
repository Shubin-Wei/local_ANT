"""The standard figure set -> <output>/figures (each figure can also be drawn alone: `lant plot <name>`).

  1_stations             station map over coloured relief, excluded stations as red crosses (GMT)
  2_ccf_moveout_<band>   cross-correlations against distance, one figure per band (the measurement band cut at
                         round frequencies); distances shown up to the last distance bin with enough pairs
  2_ccf_stack_time_<k>   three pairs (short, middle, long): hourly correlations and SNR against stacking length
  3_dispersion_pair_<k>  the same pairs: waveform, phase and group images with the accepted picks
  3_dispersion_curves    all phase-velocity curves, mean and standard deviation
  4_start_model          average curve, 1-D model fit and the starting model
  4_kernels              depth sensitivity of the phase velocity for the starting model (fwd1d finite differences)
  5_depth_slices         areal arrays: Vs depth slices over grey relief, with the section lines (GMT)
  5_slices_3d            areal arrays: five slices in a perspective view, the relief above them (GMT)
  5_sections             Vs sections: plots.sections (areal) or along the line (line) (GMT)
  6_checkerboard_*       checkerboard input and recovery: maps (areal), sections along the line (line) (GMT)
  6_recovery             recovery slope and correlation with depth
  6_residuals            travel-time residuals of the starting and final model
Figures of the processing steps are written by the steps themselves (work/01_clean/coverage.png,
work/03_stack/station_qc/station_qc.png, work/04_reference/fj.png).
"""
import json
import os

from .. import paths
from . import ccf, dispersion, inversion, maps, sections
from . import gmt as G

NAMES = ["stations", "ccf_moveout", "ccf_stack_time", "dispersion_pair", "dispersion_curves", "start_model", "kernels",
         "depth_slices", "slices_3d", "sections", "checkerboards", "recovery", "residuals"]
# figures that exist in several numbered or banded copies: removed before they are drawn again
MULTI = {"ccf_moveout": "2_ccf_moveout*", "ccf_stack_time": "2_ccf_stack_time*", "dispersion_pair": "3_dispersion_pair*",
         "checkerboards": "6_checkerboard_*", "depth_slices": "5_depth_slices*", "slices_3d": "5_slices_3d*"}


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
    if p["slice_depths_km"] == "auto":                # the model file holds the updated depths only
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
        self.done, self.skipped = {}, {}
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
            return False
        return True

    def attempt(self, n, fn):
        try:
            self.done[n] = fn()
        except Exception as exc:
            self.skipped[n] = "%s: %s" % (type(exc).__name__, exc)
            print("  figure %s not drawn: %s" % (n, exc), flush=True)

    def remove_old(self):
        import glob
        for n, pat in MULTI.items():
            if n in self.todo:
                for f in glob.glob(os.path.join(self.out, pat)):
                    os.remove(f)


def _station_map(F):
    if F.want("stations"):
        F.attempt("stations", lambda: maps.station_map(F.cfg, F.stations, F.path("1_stations"), excluded=F.excluded))


def _moveout(F):
    """One moveout figure per band; distance and velocity window from the measurements in the band."""
    if not F.want("ccf_moveout", [F.stack]):
        return
    band, _ = ccf.figure_band(F.cfg)
    for lo, hi in ccf.figure_bands(*band):
        nm = "2_ccf_moveout_%s-%sHz" % (ccf.hz(lo), ccf.hz(hi))
        win = ccf.band_window(F.disp, (lo, hi)) if os.path.exists(F.disp) else None
        if os.path.exists(F.disp) and not win:
            F.skipped[nm] = "no phase velocity was measured in this band"
            continue
        kw = dict(dmax=win["dmax"], vref=((win["vfast"], "#ff1a1a"), (win["vslow"], "#00b4ff"))) if win else {}
        F.attempt(nm, lambda lo=lo, hi=hi, nm=nm, kw=kw: ccf.moveout(F.cfg, F.stack, F.path(nm), band=(lo, hi), **kw))


def _example_pairs(F):
    """The same example pairs (plots.pair_example, or short / middle / long) for the stacking and dispersion figures."""
    if not (os.path.exists(F.disp) and {"ccf_stack_time", "dispersion_pair"} & set(F.todo)):
        return
    chosen = dispersion.choose_pairs(F.disp)
    given = F.cfg["plots"]["pair_example"]
    if given:
        given = [given] if isinstance(given, str) else list(given)
        bands = dict(dispersion.choose_pairs(F.disp, n=10 ** 6, min_points=1))
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
            F.attempt("ccf_stack_time_%d" % k, lambda pair=pair, band=band, k=k, vw=vwin(band): ccf.stack_time(
                F.cfg, F.stack, pair, F.path("2_ccf_stack_time_%d" % k), band=band, bad_days=bad, vwin=vw))
    if F.want("dispersion_pair", [F.stack, F.disp]):
        for k, (pair, band) in enumerate(chosen, 1):
            F.attempt("dispersion_pair_%d" % k, lambda pair=pair, k=k, vw=vwin(band): dispersion.pair_example(
                F.cfg, F.disp, F.stack, pair, F.path("3_dispersion_pair_%d" % k), vwin=vw))


def _curves_and_start(F):
    curves = paths.P(F.cfg, paths.DISP, "curves.csv")
    if F.want("dispersion_curves", [curves]):
        F.attempt("dispersion_curves", lambda: dispersion.all_curves(F.cfg, curves, F.path("3_dispersion_curves")))
    vs1d = paths.P(F.cfg, paths.INV, "vs1d.txt")
    if F.want("start_model", [vs1d]):
        F.attempt("start_model", lambda: inversion.start_model(F.cfg, paths.P(F.cfg, paths.INV), F.path("4_start_model")))
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
        traces = [(s["name"], s["start"], s["end"]) for s in cfg["plots"]["sections"]] if cfg["plots"]["sections"] else None
        if F.want("depth_slices", [F.model]):
            F.attempt("depth_slices", lambda: maps.depth_slices(cfg, F.model, F.stations, F.path("5_depth_slices"),
                                                                traces=traces))
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
    root = paths.P(cfg, paths.INV, "checkerboards")
    names = sorted(os.listdir(root)) if os.path.isdir(root) else []
    vert = {c["name"] for c in cfg["checks"]["vertical"]} if isinstance(cfg["checks"]["vertical"], list) else set()
    vert |= {n for n in names if n.startswith("cb_vertical")}
    for name in names:
        d = os.path.join(root, name)
        if not os.path.exists(os.path.join(d, "DSurfTomo.inSyn.dat")):
            continue
        if F.line:
            F.attempt("checkerboard_" + name, lambda d=d, name=name: sections.line_checkerboard(
                cfg, d, F.stations, F.path("6_checkerboard_" + name), label=name))
        else:
            F.attempt("checkerboard_" + name, lambda d=d, name=name: maps.checkerboard(
                cfg, d, F.stations, F.path("6_checkerboard_" + name), vertical=name in vert))


def _evaluation(F):
    if F.ev and "recovery" in F.todo:
        F.attempt("recovery", lambda: inversion.recovery(F.cfg, F.ev, F.path("6_recovery")))
    res_csv = paths.P(F.cfg, paths.INV, "evaluation", "residuals_final.csv")
    if F.ev and F.want("residuals", [res_csv]):
        F.attempt("residuals", lambda: inversion.residuals(F.cfg, res_csv, F.ev, F.path("6_residuals")))


GROUPS = [_station_map, _moveout, _example_pairs, _curves_and_start, _model, _checkerboards, _evaluation]


def make_all(cfg, which=None):
    """Draw the figures `which` (names of NAMES; all by default) -> dict(drawn, not_drawn, plot_parameters)."""
    F = _Figures(cfg, which)
    notes = resolve(cfg, F.model, F.stations)
    F.remove_old()
    G.RELIEF_NOTES.clear()
    for group in GROUPS:
        group(F)
    if G.RELIEF_NOTES:
        F.skipped["relief"] = "; ".join(G.RELIEF_NOTES)
    return dict(drawn=F.done, not_drawn=F.skipped,
                plot_parameters={k: dict(value=cfg["plots"][k], rule=v) for k, v in notes.items()})
