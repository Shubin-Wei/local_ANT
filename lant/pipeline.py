"""The steps, in order.  Each step reads the outputs of the earlier ones, so it can be rerun alone once those exist.

  clean      01_clean      raw data of any layout / format -> index, checks, cleaned day archive
  cc         02_ccf        hourly cross-coherence, one file per processing day
  stack      03_stack      stack, station-day (areal) or station-hour (line) screening, cleaned stack
  notch      04_notch      notch table by rule, final cross-correlation product stack.h5
  reference  05_reference  first-pass reference (F-J of this array) or a borrowed curve
  disp1      06_dispersion first pass (no deviation screens) -> second-pass reference (median)   [F-J source only]
  disp2      06_dispersion formal measurement, QC, precision, export of the curves, single-pair figures
  inversion  07_inversion  average curve, depth nodes, 1-D starting model, main project, checkerboards, controls
  run        07_inversion  run DSurfTomo / pDSurfTomo on the main project, then the check projects
  evaluate   07_inversion  residuals of the final model and of each iteration, checkerboard recovery, control differences
  figures    figures       the standard figure set (lant.plotting)
"""
import json
import os

import numpy as np

from . import paths, records
from . import clean as CLEAN
from . import cc as CC
from . import stack as STACK
from . import notch as NOTCH
from . import reference as REF
from . import disp as DISP
from . import qc as QC
from . import export as EXPORT
from .inversion import start as START, project as PROJ, run as RUN

STEPS = ["clean", "cc", "stack", "notch", "reference", "disp1", "disp2", "inversion", "run", "evaluate", "figures"]


def _need(path, step):
    if not os.path.exists(path):
        raise RuntimeError("%s not found: run `lant run %s` first" % (path, step))


def step_clean(cfg):
    return CLEAN.run(cfg)


def step_cc(cfg):
    _need(paths.P(cfg, "01_clean", "archive_index.json"), "clean")
    return CC.run(cfg)


def step_stack(cfg):
    return STACK.run(cfg)


def step_notch(cfg):
    return NOTCH.run(cfg)


def step_reference(cfg):
    _need(paths.stack_product(cfg), "notch")
    return REF.run(cfg)


def ref_final(cfg):
    return paths.P(cfg, paths.REFERENCE, "reference.npz")


def _measure_chain(cfg, ref, tag, first_pass):
    d = paths.mk(paths.P(cfg, paths.DISP))
    stack = paths.stack_product(cfg)
    raw = os.path.join(d, "disp_%s_raw.h5" % tag)
    s = {"measure": DISP.measure(cfg, stack, ref, raw, first_pass=first_pass,
                                 save_images=0 if first_pass else max(0, int(cfg["export"]["curve_figures"])))}
    q = os.path.join(d, "disp_%s_qc.h5" % tag)
    s["group_qc"] = QC.group_qc(cfg, raw, q, first_pass=first_pass)
    fin = os.path.join(d, "disp_%s.h5" % tag)
    s["lambda_qc"] = QC.lambda_qc(cfg, q, fin)
    os.remove(q)          # the raw measurement is kept: `lant run` of the QC alone does not need to re-measure
    return fin, s


def step_disp1(cfg):
    r = cfg["reference"]
    if r["source"] == "file":
        return dict(skipped="reference.source is file: the borrowed curve is the measurement reference")
    ref1 = paths.P(cfg, paths.REFERENCE, "reference_pass1.npz")
    _need(ref1, "reference")
    if not r["pass2"]:
        import shutil
        shutil.copyfile(ref1, ref_final(cfg))
        return dict(skipped="reference.pass2 is false: the F-J curve is the measurement reference")
    fin, s = _measure_chain(cfg, ref1, "pass1", first_pass=True)
    s["pass2"] = REF.pass2(cfg, ref1, fin, ref_final(cfg), paths.P(cfg, paths.REFERENCE, "reference_pass2"))
    return s


def step_disp2(cfg):
    ref = ref_final(cfg)
    _need(ref, "disp1" if cfg["reference"]["source"] == "fj" else "reference")
    fin, s = _measure_chain(cfg, ref, "final", first_pass=False)
    d = paths.P(cfg, paths.DISP)
    prec, txt = QC.redundant_precision(fin, paths.stack_product(cfg))
    open(os.path.join(d, "precision.txt"), "w").write(txt)
    s["precision"] = prec
    s["export"] = EXPORT.export(cfg, fin, paths.stack_product(cfg), d, prec)
    n = int(cfg["export"]["curve_figures"])
    if n:
        from .plotting import dispersion as PD
        s["curve_figures"] = PD.pair_figures(cfg, fin, os.path.join(d, "pair_figures"), n)
    return s


def inv_dir(cfg, *p):
    return paths.P(cfg, paths.INV, *p)


def step_inversion(cfg):
    curves_dir = paths.P(cfg, paths.DISP, "curves")
    _need(curves_dir, "disp2")
    d = paths.mk(inv_dir(cfg))
    main = inv_dir(cfg, "main")
    if RUN.status(main) in ("done", "unfinished", "stale"):
        raise RuntimeError("%s holds DSurfTomo output: move it away (or delete it) before building a new project" % main)
    z, vs, avg, info = START.build(cfg, curves_dir, d)
    e = cfg["export"]
    curves = EXPORT.read_curves(curves_dir, "C", e["tmin_s"], info["tmax_used"])
    meta = PROJ.build_main(cfg, curves, z, vs, main)
    checks = PROJ.build_checks(cfg, main, d)
    return dict(start=info, project={k: v for k, v in meta.items() if k != "grid"}, grid=meta["grid"], checks=checks)


def check_dirs(cfg):
    out = []
    for sub in ("checkerboards", "controls"):
        p = inv_dir(cfg, sub)
        if os.path.isdir(p):
            out += [os.path.join(p, n) for n in sorted(os.listdir(p))]
    return out


def step_run(cfg, what="all"):
    main = inv_dir(cfg, "main")
    _need(os.path.join(main, "DSurfTomo.in"), "inversion")
    dirs = ([main] if what in ("main", "all") else []) + (check_dirs(cfg) if what in ("checks", "all") else [])
    res = RUN.run_projects(cfg, dirs)
    if res["failed"]:
        raise RuntimeError("%d project(s) failed:\n  " % len(res["failed"])
                           + "\n  ".join("%s: %s" % (d, res["results"][d]) for d in res["failed"]))
    return res


def step_evaluate(cfg):
    from . import data as D
    main = inv_dir(cfg, "main")
    if RUN.status(main) != "done":
        raise RuntimeError("%s has not finished: run `lant run run` first" % main)
    out = paths.mk(inv_dir(cfg, "evaluation"))
    s = {"residuals": RUN.residuals(cfg, main, os.path.join(out, "residuals_final"), cfg["checks"]["iteration_residuals"])}
    stations = D.used_stations(cfg)
    s["checkerboards"] = {}
    for d in check_dirs(cfg):
        name = os.path.basename(d)
        if os.path.dirname(d).endswith("checkerboards") and RUN.status(d) == "done":
            s["checkerboards"][name] = RUN.checker_recovery(d, stations)
    s["controls"] = {}
    base = RUN.read_model(os.path.join(main, "DSurfTomo.inMeasure.dat"))
    for d in check_dirs(cfg):
        if os.path.dirname(d).endswith("controls") and RUN.status(d) == "done":
            m = RUN.read_model(os.path.join(d, "DSurfTomo.inMeasure.dat"))
            diff = m[3] - base[3]
            s["controls"][os.path.basename(d)] = dict(
                rms_km_s=float(np.sqrt(np.nanmean(diff ** 2))), max_abs_km_s=float(np.nanmax(np.abs(diff))),
                rms_by_depth={"%.2f" % z: float(np.sqrt(np.nanmean(diff[:, :, k] ** 2))) for k, z in enumerate(base[2])})
    json.dump(s, open(os.path.join(out, "evaluation.json"), "w"), indent=1)
    return s


def step_figures(cfg, which=None):
    from .plotting import figures as F
    return F.make_all(cfg, which)


FUNCS = {n: globals()["step_" + n] for n in STEPS}


def run_step(cfg, name, **kw):
    from . import auto
    if name not in FUNCS:
        raise ValueError("unknown step %s (steps: %s)" % (name, ", ".join(STEPS)))
    if name != "clean":
        cfg = auto.apply(cfg)
        left = auto.unresolved(cfg)
        if left:
            raise RuntimeError("parameters %s are 'auto' and not resolved yet: run `lant run clean` first" % ", ".join(left))
    with records.Step(cfg, name):
        s = FUNCS[name](cfg, **kw)
        records.write_summary(cfg, name, s)
    return s
