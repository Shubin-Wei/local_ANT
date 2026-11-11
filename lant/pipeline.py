"""The steps, in order.  Each step reads the outputs of the earlier ones, so it can be rerun alone once those exist.

  clean      01_clean      raw data of any layout / format -> index, checks, cleaned day archive; parameters from the
                           geometry and the velocity range                                         [confirm clean]
  cc         02_ccf        hourly cross-coherence, one file per processing day
  stack      03_stack      stack and sub-stacks; diagnosis: candidate station sub-stacks;
                           only the listed (or confirmed) ones are removed                         [confirm stack]
  reference  04_reference  F-J reference of this array (or a borrowed curve), measurement band and velocity axis
                                                                                                   [confirm reference]
  disp       05_dispersion two-station phase velocities, QC, export of the curves, single-pair figures
  inversion  06_inversion  average curve, depth nodes, 1-D starting model, main project, checkerboards, controls
                                                                                                   [confirm inversion]
  run        06_inversion  precheck, then DSurfTomo / pDSurfTomo on the main project and the check projects
  evaluate   06_inversion  residuals of the final model, checkerboard recovery, control differences
  figures    figures       the standard figure set (lant.plotting)
A step after a confirmation point does not run before `lant confirm <step>` (lant.confirm).
"""
import json
import shutil
import os

import numpy as np

from . import paths, records
from . import clean as CLEAN
from . import cc as CC
from . import stack as STACK
from . import reference as REF
from . import disp as DISP
from . import qc as QC
from . import export as EXPORT
from .inversion import start as START, project as PROJ, run as RUN

STEPS = ["clean", "cc", "stack", "reference", "disp", "inversion", "run", "evaluate", "figures"]


def _need(path, step):
    if not os.path.exists(path):
        raise RuntimeError("%s not found: run `lant run %s` first" % (path, step))


def step_clean(cfg):
    return CLEAN.run(cfg)


def step_cc(cfg):
    _need(paths.P(cfg, "01_clean", "archive_index.json"), "clean")
    return CC.run(cfg)


def step_stack(cfg):
    _need(paths.P(cfg, paths.CCF), "cc")
    return STACK.run(cfg)


def step_reference(cfg):
    _need(paths.stack_product(cfg), "stack")
    return REF.run(cfg)


def ref_final(cfg):
    return paths.P(cfg, paths.REFERENCE, "reference.npz")


def curve_figure_count(cfg):
    """export.curve_figures: 0 none, n > 0 that many pairs spread over distance (plus the conspicuous ones), -1 every pair."""
    n = int(cfg["export"]["curve_figures"])
    return 10 ** 9 if n < 0 else n


def step_disp(cfg):
    ref = ref_final(cfg)
    _need(ref, "reference")
    d = paths.mk(paths.P(cfg, paths.DISP))
    fin = os.path.join(d, "disp.h5")
    s = {"measure": DISP.measure(cfg, paths.stack_product(cfg), ref, fin)}
    s["lambda_qc"] = QC.lambda_qc(cfg, fin)
    s["rejections"] = QC.reason_counts(fin)
    s["median_vs_reference"] = QC.median_vs_reference(fin, os.path.join(d, "median_vs_reference.txt"))
    s["export"] = EXPORT.export(cfg, fin, paths.stack_product(cfg), d)
    n = int(cfg["export"]["curve_figures"])
    if n:
        from .plotting import dispersion as PD
        s["curve_figures"] = PD.pair_figures(cfg, fin, os.path.join(d, "pair_figures"), n)
    return s


def inv_dir(cfg, *p):
    return paths.P(cfg, paths.INV, *p)


def step_inversion(cfg):
    curves_dir = paths.P(cfg, paths.DISP, "curves")
    _need(curves_dir, "disp")
    d = paths.mk(inv_dir(cfg))
    main = inv_dir(cfg, "main")
    z, vs, avg, info = START.build(cfg, curves_dir, d)
    inv = cfg["inversion"]
    from . import auto
    from .signal import nice_number
    steps = (1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8)
    rule = {"inversion.depth_sequence_km": dict(value=[float(x) for x in START.auto_sequence(avg)],
                                                reason="steps max(shortest wavelength / 4, depth / 4) to half the longest "
                                                       "wavelength (average curve)"),
            "inversion.vmin": dict(value=round(nice_number(0.6 * float(np.min(vs)), "down", steps), 4),
                                   reason="0.6 x the slowest Vs of the starting model %.3f km/s" % float(np.min(vs))),
            "inversion.vmax": dict(value=round(nice_number(1.4 * float(np.max(vs)), "up", steps), 4),
                                   reason="1.4 x the fastest Vs of the starting model %.3f km/s" % float(np.max(vs)))}
    auto.write(cfg, rule, None)                     # the rule's values, proposed at the confirmation point
    if inv["depth_sequence_km"] == "auto":
        inv["depth_sequence_km"] = info["sequence"]
    for k in ("vmin", "vmax"):
        if inv[k] == "auto":
            inv[k] = rule["inversion." + k]["value"]
            info[k + "_rule"] = rule["inversion." + k]["reason"]
    info.update(vmin=inv["vmin"], vmax=inv["vmax"])
    e = cfg["export"]
    curves = EXPORT.read_curves(curves_dir, "C", e["tmin_s"], info["tmax_used"])
    shutil.rmtree(main + ".new", ignore_errors=True)
    meta = PROJ.build_main(cfg, curves, z, vs, main + ".new")
    main_state = PROJ.settle(main + ".new", main)  # a main project that has run is kept if its inputs are the same
    if auto.noise_is_auto(cfg):                     # set in the run step from the fit of the (new) main model
        auto.reset(cfg, [("inversion", "noise")])
        noise = None
    else:
        noise = float(inv["noise"])
    checks = PROJ.build_checks(cfg, main, d, noise=noise)
    projects = [("checkerboards", n) for n, _ in checks["built"] if not n.startswith(("reg_", "start_"))] + \
               [("checkerboards", n) for n in checks["deferred"]] + \
               [("controls", n) for n, _ in checks["built"] if n.startswith(("reg_", "start_"))]
    json.dump(dict(projects=projects), open(inv_dir(cfg, "checks.json"), "w"), indent=1)
    return dict(start=info, project={k: v for k, v in meta.items() if k != "grid"}, main=main_state, grid=meta["grid"],
                checks=checks, short_paths=short_path_note(curves))


SHORT_PATH_KM = 0.2     # empirical: DSurfTomo keeps coordinates and distances in single precision; on two test
                        # projects its distance differed from the double-precision one by 0.3-1.1 % near 80 m and by
                        # 4.4 % at 15 m (depends on position and orientation)


def short_path_note(curves):
    """Shortest station pair going into DSurfTomo, and how many pairs are shorter than SHORT_PATH_KM."""
    from .signal import geodesic_km
    if not curves:
        return {}
    d = np.array([geodesic_km(c["lon1"], c["lat1"], c["lon2"], c["lat2"])[0] for c in curves], float)
    n = int((d < SHORT_PATH_KM).sum())
    out = dict(shortest_pair_km=round(float(d.min()), 4), pairs_shorter_than_km={str(SHORT_PATH_KM): n})
    if n:
        out["note"] = ("%d of %d pairs are shorter than %.1f km: DSurfTomo computes their distances in single "
                       "precision, an error that reached 0.3-1.1 %% near 0.08 km in tests" % (n, len(d), SHORT_PATH_KM))
        print("  note: " + out["note"], flush=True)
    return out


def check_dirs(cfg):
    """The check projects built for the present configuration (06_inversion/checks.json), not every directory that
    happens to be there."""
    p = inv_dir(cfg, "checks.json")
    if not os.path.exists(p):
        return []
    return [inv_dir(cfg, sub, n) for sub, n in json.load(open(p))["projects"] if os.path.isdir(inv_dir(cfg, sub, n))]


def checkerboard_noise(cfg, main):
    """inversion.noise auto: the relative rms of the residuals of the final main model that DSurfTomo's rejection
    rule keeps (one forward run; kept in 06_inversion/checkerboard_noise.json while the main model is unchanged)."""
    p = inv_dir(cfg, "checkerboard_noise.json")
    model = os.path.join(main, "DSurfTomo.inMeasure.dat")
    stamp = os.path.getmtime(model)
    if os.path.exists(p):
        old = json.load(open(p))
        if old.get("model_mtime") == stamp:
            return old
    _, R, _ = RUN.evaluate(cfg, main)
    kept = R[:, 5] != 0
    rel = (R[kept, 2] - R[kept, 1]) / R[kept, 2]
    out = dict(value=float(np.sqrt(np.mean(rel ** 2))), n_kept=int(kept.sum()), n=int(len(R)), model_mtime=stamp,
               reason="relative rms of the %d residuals of the final main model kept by DSurfTomo's rejection rule "
                      "(of %d)" % (int(kept.sum()), len(R)))
    json.dump(out, open(p, "w"), indent=1)
    return out


def _run_batch(cfg, dirs, main):
    invalid = precheck_projects(cfg, [d for d in dirs if RUN.status(d) == "new"])
    if main in invalid:
        raise RuntimeError("the main project cannot be inverted: %s (06_inversion/main/precheck.json); check the starting "
                           "model and the data" % "; ".join(invalid[main]))
    res = RUN.run_projects(cfg, [d for d in dirs if d not in invalid])
    res["invalid"] = invalid
    return res


def step_run(cfg, what="all"):
    """Main project and controls side by side; the checkerboards with them when inversion.noise is a number, after
    the main project when it is auto (their noise follows from the fit of the main model)."""
    from . import auto
    main = inv_dir(cfg, "main")
    _need(os.path.join(main, "DSurfTomo.in"), "inversion")
    want_checks = what in ("checks", "all")
    pending = []
    if want_checks and os.path.exists(inv_dir(cfg, "checks.json")):
        pending = [n for sub, n in json.load(open(inv_dir(cfg, "checks.json")))["projects"]
                   if sub == "checkerboards" and auto.noise_is_auto(cfg)]
    first = [d for d in check_dirs(cfg) if os.path.basename(d) not in pending] if want_checks else []
    res = _run_batch(cfg, ([main] if what in ("main", "all") else []) + first, main)
    if pending:
        if RUN.status(main) != "done":
            raise RuntimeError("the checkerboards need the finished main project (inversion.noise: auto): run "
                               "`lant run run` first")
        nz = checkerboard_noise(cfg, main)
        print("  checkerboard noise %.4f: %s" % (nz["value"], nz["reason"]), flush=True)
        auto.write(cfg, {"inversion.noise": dict(value=round(nz["value"], 4), reason=nz["reason"])}, None)
        built = PROJ.build_checks(cfg, main, inv_dir(cfg), noise=nz["value"], part="checkerboards")
        res2 = _run_batch(cfg, [inv_dir(cfg, "checkerboards", n) for n, _ in built["built"]], main)
        res = dict(results={**res["results"], **res2["results"]}, failed=res["failed"] + res2["failed"],
                   invalid={**res["invalid"], **res2["invalid"]})
        res["checkerboard_noise"] = {k: nz[k] for k in ("value", "reason")}
    elif want_checks and check_dirs(cfg):
        res["checkerboard_noise"] = dict(value=float(cfg["inversion"]["noise"]), reason="inversion.noise in the configuration")
    for d, why in res["invalid"].items():
        print("  %s INVALID, not run: %s" % (d, "; ".join(why)), flush=True)
    if res["failed"]:
        raise RuntimeError("%d project(s) failed:\n  " % len(res["failed"])
                           + "\n  ".join("%s: %s" % (d, res["results"][d]) for d in res["failed"]))
    return res


def precheck_projects(cfg, dirs):
    """lant.inversion.run.precheck for every project about to run -> {dir: reasons} of the invalid ones."""
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=int(cfg["inversion"]["njobs"])) as ex:
        res = list(ex.map(lambda d: (d, RUN.precheck(cfg, d)), dirs))
    for d, r in res:
        st = r["stats"]
        print("  precheck %s: median residual %.2f ms, kept by the rejection rule %.1f %% -> %s"
              % (d, st["median_ms"], 100 * st["kept_fraction"], "ok" if r["valid"] else "INVALID"), flush=True)
    return {d: r["reasons"] for d, r in res if not r["valid"]}


def precheck_state(d):
    """(valid, reasons) from <project>/precheck.json, or (None, []) when the project was not checked."""
    p = os.path.join(d, "precheck.json")
    if not os.path.exists(p):
        return None, []
    r = json.load(open(p))
    return r["valid"], r["reasons"]


def recovery_radius(cfg):
    """Checkerboard statistics use the nodes within twice the median nearest-station distance of a station."""
    from . import auto
    nn = auto.load(cfg).get("geometry", {}).get("median_nearest_km")
    if not nn:
        raise RuntimeError(".lant/auto_parameters.json has no geometry: run `lant run clean` first")
    return round(2 * float(nn), 3)


def step_evaluate(cfg):
    from . import data as D
    main = inv_dir(cfg, "main")
    if RUN.status(main) != "done":
        raise RuntimeError("%s has not finished: run `lant run run` first" % main)
    out = paths.mk(inv_dir(cfg, "evaluation"))
    s = {"residuals": RUN.residuals(cfg, main, os.path.join(out, "residuals_final"), cfg["checks"]["iteration_residuals"])}
    stations = D.used_stations(cfg)
    s["invalid"] = {os.path.basename(d): precheck_state(d)[1] for d in check_dirs(cfg) if precheck_state(d)[0] is False}
    s["checkerboards"] = {}
    for d in check_dirs(cfg):
        name = os.path.basename(d)
        if os.path.dirname(d).endswith("checkerboards") and RUN.status(d) == "done":
            s["checkerboards"][name] = RUN.checker_recovery(d, stations, edge_km=recovery_radius(cfg))
    s["controls"] = {}
    base = RUN.read_model(os.path.join(main, "DSurfTomo.inMeasure.dat"))
    for d in check_dirs(cfg):
        if os.path.dirname(d).endswith("controls") and RUN.status(d) == "done":
            m = RUN.read_model(os.path.join(d, "DSurfTomo.inMeasure.dat"))
            if any(a.shape != b.shape or not np.allclose(a, b) for a, b in zip(m[:3], base[:3])) or \
                    np.isnan(m[3]).any() or np.isnan(base[3]).any():
                raise RuntimeError("%s: the model grid differs from the main project (or has missing nodes)" % d)
            diff = m[3] - base[3]
            s["controls"][os.path.basename(d)] = dict(
                rms_km_s=float(np.sqrt(np.nanmean(diff ** 2))), max_abs_km_s=float(np.nanmax(np.abs(diff))),
                rms_by_depth={"%.2f" % z: float(np.sqrt(np.nanmean(diff[:, :, k] ** 2))) for k, z in enumerate(base[2])},
                residual_std_ms=RUN.logged_residuals(d), main_residual_std_ms=RUN.logged_residuals(main))
    json.dump(s, open(os.path.join(out, "evaluation.json"), "w"), indent=1)
    return s


def step_figures(cfg, which=None):
    from .plotting import figures as F
    return F.make_all(cfg, which)


FUNCS = {n: globals()["step_" + n] for n in STEPS}


def run_step(cfg, name, **kw):
    """Run one step: only when every earlier step stands (lant.manifest.gate); the run is written to the manifest."""
    import time
    from . import auto, confirm, manifest, provenance, report
    if name not in FUNCS:
        raise ValueError("unknown step %s (steps: %s)" % (name, ", ".join(STEPS)))
    if name in confirm.VALUE_KEYS and manifest.upstream_changed(cfg, name):
        cfg.get("_confirmed", {}).pop(name, None)      # confirmed for the results of other inputs: the rule decides anew
    manifest.gate(report.resolved(cfg), name)
    cfg = dict(cfg)
    if name != "clean":
        cfg = auto.apply(cfg)
        left = auto.unresolved(cfg, auto.CLEAN_KEYS)
        if left:
            raise RuntimeError("parameters %s are 'auto' and not resolved yet: run `lant run clean -c %s` first"
                               % (", ".join(left), cfg["_config_file"]))
        if name in ("disp", "inversion", "run", "evaluate"):
            left = auto.unresolved(cfg, auto.REFERENCE_KEYS)
            if left:
                raise RuntimeError("parameters %s are 'auto' and not resolved yet: run `lant run reference -c %s` first"
                                   % (", ".join(left), cfg["_config_file"]))
    with records.Step(cfg, name) as st:
        code = provenance.record(records.rec_dir(cfg))
        inputs = manifest.raw_inputs(cfg) if name == "clean" else None
        t0 = time.time()
        try:
            s = FUNCS[name](cfg, **kw)
            records.write_summary(cfg, name, s)
            st.failed = list(s.get("failed") or []) if isinstance(s, dict) else []
            dropped = name in confirm.VALUE_KEYS and name not in cfg.get("_confirmed", {})
            rcfg = report.resolved(cfg, without=name if dropped else None)
            items = confirm.judgements(rcfg, name, s) if name in confirm.CONFIRM_STEPS else []
            from .delivery import publish
            if not st.failed:
                publish(rcfg, name)
            rec = manifest.record(
                cfg, name, "partial" if st.failed else "done", seconds=round(time.time() - t0, 1),
                params=report.used_params(rcfg, name), to_confirm=items, basis=confirm.basis(rcfg, name, items),
                error="%d items failed, e.g. %s" % (len(st.failed), "; ".join(st.failed[:3])) if st.failed else None,
                code=code, inputs=inputs)
            try:
                from .delivery import write_report
                report.write_step_report(rcfg, name, s, rec)
                write_report(rcfg)
            except Exception as exc:
                records.log(cfg, step=name, event="report_failed", error=str(exc))
                import sys
                print("[lant] report failed: %s; see run_log.jsonl" % exc, file=sys.stderr)
        except BaseException as exc:
            manifest.record(cfg, name, "failed", seconds=round(time.time() - t0, 1),
                            error="%s: %s" % (type(exc).__name__, exc), code=code, inputs=inputs)
            raise
        if st.failed:
            raise records.StepIncomplete(
                "step %s: %d items failed (listed in the step summary); repair or exclude failed inputs, then "
                "run lant run all -c %s" % (name, len(st.failed), cfg["_config_file"]))
    return s
