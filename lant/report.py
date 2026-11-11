"""Step reports and the plan: what a step did (in order, with counts), the key numbers, every parameter it used
with its value, what the value rests on (lant.params) and where it comes from (written by the user, confirmed,
derived by a rule with the numbers of the present data, or the default), what to look at, what is left for the user
to judge, and what to change - kept in <output>/.lant/steps/<step>.md (printed with -v).
`lant plan` shows the same parameter table for every step before anything runs.

Nothing here changes a result.
"""
import json
import os

from . import auto
from . import config as C
from . import records, paths

# The configuration each step reads (sections or single keys).  A change in any of them makes the step and every
# later step out of date.
STEP_KEYS = {
    "clean": ["project.array_type", "data", "clean", "correlation.fs_max", "correlation.fs", "correlation.fmin_hz", "correlation.fmax_hz",
              "correlation.maxlag_s", "correlation.seg_len_s", "dispersion.alpha", "dispersion.per_octave",
              "dispersion.fmax_hz", "stack.block_hours", "inversion.grid_spacing_deg", "inversion.grid_refine"],
    "cc": ["correlation"],
    "stack": ["stack"],
    "reference": ["reference", "dispersion.fmin_hz", "dispersion.fmax_hz", "dispersion.per_octave", "dispersion.vmin",
                  "dispersion.vmax", "dispersion.dv", "dispersion.qc.lam_c"],
    "disp": ["dispersion", "export"],
    "inversion": ["inversion", "checks", "export.tmin_s", "export.tmax_s"],
    "run": ["inversion.solver", "inversion.threads", "inversion.njobs", "inversion.noise", "checks.precheck"],
    "evaluate": ["checks"],
    "figures": ["plots"],
}

# The parameters one usually changes after looking at a step, and what they do.
TUNE = {
    "clean": [("data.velocity_range_kms", "the approximate velocities: lag window, band, F-J scan"),
              ("data.exclude_stations", "leave stations out (instrument or coupling problems)"),
              ("data.rename", "match header names to the station table by hand"),
              ("data.channel", "choose the channel when the names do not tell the vertical one"),
              ("data.start_date / end_date", "process only part of the record"),
              ("correlation.fs_max", "the highest processing rate, which sets the top of the band")],
    "cc": [("correlation.max_over_std", "segment rejection (transients)"),
           ("correlation.remove_1s_comb", "the 1-s timing-pulse removal (Takagi et al. 2015)")],
    "stack": [("lant confirm stack --exclude-candidates", "remove the candidate station sub-stacks"),
              ("stack.exclude_station_days", "station-days of your own, also ranges [station, \"YYYYDDD-\"]")],
    "reference": [("reference.fj_mode", "max_peak / track_fundamental when the picks jump to another ridge"),
                  ("reference.track_start_band_hz", "where the fundamental ridge is unambiguous"),
                  ("reference.fj_fit_fmax", "cap the fitted picks"),
                  ("reference.source: file", "a borrowed reference curve"),
                  ("dispersion.fmin_hz / fmax_hz", "the measured band")],
    "disp": [("export.tmin_s / tmax_s", "periods passed on")],
    "inversion": [("inversion.depth_sequence_km", "depth nodes"), ("inversion.weight / damp", "regularisation"),
                  ("inversion.start_model", "a starting model of your own"), ("checks.*", "checkerboards and controls")],
    "run": [("inversion.solver", "pdsurftomo: parallel, same model"), ("inversion.threads / njobs", "speed and memory")],
    "evaluate": [("checks.iteration_residuals", "residuals of every iteration (one forward run each)")],
    "figures": [("plots.*", "colour range, slice depths, sections, topography, font size")],
}

LOOK_AT = {
    "clean": ["work/01_clean/report.txt", "work/01_clean/coverage.png", "work/01_clean/stations.csv", "work/01_clean/processing_days.txt"],
    "cc": ["work/02_ccf/"],
    "stack": ["work/03_stack/candidates.txt", "work/03_stack/station_qc/station_qc.png", "work/03_stack/station_qc/station_qc.txt"],
    "reference": ["work/04_reference/fj.png", "work/04_reference/fj.txt"],
    "disp": ["work/05_dispersion/pair_figures/", "work/05_dispersion/median_vs_reference.txt"],
    "inversion": ["work/06_inversion/average_curve.txt", "work/06_inversion/vs1d_fit.txt", "work/06_inversion/vs1d.txt"],
    "run": ["work/06_inversion/main/run.log", "work/06_inversion/main/precheck.json"],
    "evaluate": ["work/06_inversion/evaluation/evaluation.json"],
    "figures": ["results/figures/"],
}

RERUN = "after a change, run this step again and every step after it (`lant run all --from %s -c ...`)"


def _get(cfg, dotted):
    cur = cfg
    for k in dotted.split("."):
        if not isinstance(cur, dict) or k not in cur:
            return None
        cur = cur[k]
    return cur


def _flatten(prefix, v, out):
    if isinstance(v, dict):
        for k, x in v.items():
            _flatten(prefix + "." + k if prefix else k, x, out)
    else:
        out[prefix] = v
    return out


def resolved(cfg, without=None):
    """The configuration as it now stands in its file, with every "auto" value resolved (confirmed values, then the
    rules' values of .lant/auto_parameters.json): the same view for the recorded and the present parameters.
    without: a step whose confirmation is left out (lant.pipeline.run_step)."""
    if not cfg.get("_config_file"):
        return cfg
    new = C.load(cfg["_config_file"])
    if without:
        new["_confirmed"].pop(without, None)
    new["_verbose"] = cfg.get("_verbose", False)
    return auto.apply(new)


def used_params(cfg, step):
    """{dotted key: value} of the configuration the step reads (after the "auto" values are resolved), in the form
    the run log keeps (JSON: dates as text)."""
    out = {}
    for k in STEP_KEYS.get(step, []):
        _flatten(k, _get(cfg, k), out)
    for k, setting in auto.DERIVED_LATER.get(step, {}).items():
        if k in out:
            out[k] = setting(cfg)
    root = cfg["project"]["output_dir"]
    def portable(value):
        if isinstance(value, str) and os.path.isabs(value):
            try:
                if os.path.commonpath([root, value]) == root:
                    return "@project/" + os.path.relpath(value, root)
            except ValueError:
                pass
        if isinstance(value, list):
            return [portable(v) for v in value]
        if isinstance(value, dict):
            return {k: portable(v) for k, v in value.items()}
        return value
    return portable(json.loads(json.dumps(records.strict(out), default=str)))


def source(cfg, key):
    """'config' (written by the user), 'confirmed (<step>)', 'auto: <rule>' (proposed, not yet confirmed) or
    'default'."""
    u = cfg.get("_user") or {}
    cur = u
    for k in key.split("."):
        if not isinstance(cur, dict) or k not in cur:
            cur = None
            break
        cur = cur[k]
    if cur is not None and not C.is_auto(cur):
        return "config"
    from . import confirm as CF
    if key in CF.valid_values(cfg):
        return "confirmed (%s)" % key.split(".")[0]
    dec = auto.load(cfg).get("parameters", {})
    if key in dec:
        return "auto: " + dec[key]["reason"]
    return "default"


def _round(v):
    if isinstance(v, float):
        return float("%.6g" % v)
    if isinstance(v, (list, tuple)):
        return [_round(x) for x in v]
    if isinstance(v, dict):
        return {k: _round(x) for k, x in v.items()}
    return v


def _short(v, n=70):
    s = json.dumps(_round(v), ensure_ascii=False) if not isinstance(v, str) else v
    return s if len(s) <= n else s[:n - 3] + "..."


def _numbers(summary, n=14):
    """A few top-level numbers of the step summary."""
    flat = _flatten("", summary or {}, {})
    out = []
    skip = ("day_median", "auto_parameters", "parameters.", "messages", "score_percentiles",
            "rms_by_depth", "drawn.", "not_drawn.", "plot_parameters")
    for k, v in flat.items():
        if any(x in k for x in skip) or k.count(".") > 2:
            continue
        if isinstance(v, bool) or isinstance(v, (int, float)) or (isinstance(v, str) and len(v) < 60):
            out.append((k, v))
        if len(out) >= n:
            break
    return out


def did(cfg, step, s):
    """What the step did, in the order it did it, with the counts (from the step summary)."""
    s = s or {}
    f = _DID.get(step)
    return f(cfg, lambda *k: _dig(s, *k)) if f else []


def _did_clean(cfg, g):
    L = []
    t = g("totals") or {}
    L.append("read the headers of %s candidate files (%s not seismic data or unreadable); formats %s"
             % (g("files"), g("unreadable"), _short(g("formats"))))
    n = g("names") or {}
    L.append("matched the header names to the station table: %s matched (%s under another spelling), %s ambiguous, "
             "%s table stations without data, %s records without a station"
             % (n.get("matched"), n.get("renamed"), len(n.get("ambiguous") or []),
                len(n.get("table_without_data") or []), len(n.get("data_without_table") or [])))
    L.append("chose the vertical channel and one sensor / rate per station: %s stations, raw rates %s Hz"
             % (g("stations"), _short(g("sampling_rates"))))
    L.append("cut the records into processing days starting at %02d h UTC (labels: %s): %s days, %s station-days"
             % (int(g("day_start_utc_hour") or 0), g("day_label"), g("days"), g("station_days")))
    if t:
        L.append("put %s records on the UTC sample grid: %s moved by < 0.05 sample, %s interpolated (Lanczos); "
                 "removed %s duplicates / identical overlaps; left out %.1f s of overlaps with differing samples "
                 "and %s pieces too short; %s non-finite samples became gaps; %s records flipped (sensor pointing "
                 "down); gaps are not filled (%s gaps)"
                 % (t.get("traces_in"), t.get("shifted"), t.get("interpolated"),
                    t.get("duplicates_or_identical_overlaps"), float(t.get("conflicting_overlap_s") or 0),
                    t.get("dropped_short"), t.get("nonfinite_samples"), t.get("polarity_reversed"), t.get("gaps")))
        L.append("wrote the day archive: %s station-days, %.0f station-hours; no filtering"
                 % (g("archive_station_days"), float(g("station_hours") or 0)))
    L.append("derived the parameters marked 'rule' below from the geometry, the records and "
             "data.velocity_range_kms")
    return L


def _did_cc(cfg, g):
    L = []
    c = cfg["correlation"]
    L.append("per station-day: zero-phase low-pass and decimation to %s Hz%s"
             % (c["fs"], "; 1-s periodic waveform subtracted hour by hour (%s)" % (
                 g("pulse_removed") if isinstance(g("pulse_removed"), str) else "%d stations where the clean step "
                 "detected it" % len(g("pulse_removed"))) if g("pulse_removed") else
                "; no 1-s periodic waveform subtracted (none detected)" if c["remove_1s_comb"] == "auto" else ""))
    L.append("segments of %s s every %s s; a segment is dropped when it holds a gap or its peak exceeds %s x the "
             "standard deviation of the station-day; linear detrend, 5 %% cosine taper"
             % (c["seg_len_s"], c["seg_step_s"], c["max_over_std"]))
    L.append("cross-coherence conj(Fi) Fj / (|Fi| |Fj|) in %.4g-%.4g Hz (quarter-octave cosine edges), averaged per "
             "hour, lags +-%s s" % (c["fmin_hz"], c["fmax_hz"], c["maxlag_s"]))
    L.append("%s processing days: %s computed, %s computed again (settings or input files changed), %s reused "
             "unchanged; %s days skipped (fewer than 2 stations); %s station-days failed to read"
             % (g("days"), g("computed"), g("recomputed"), g("reused"), len(g("skipped_days") or []),
                len(g("failed") or [])))
    return L


def _did_stack(cfg, g):
    L = []
    a = g("all_hours") or {}
    L.append("stacked every hour, weighted by its segments: %s pairs (%s with data), %s segments; sub-stacks of %s h "
             "for the jackknife: %s" % (a.get("pairs"), a.get("pairs_with_data"), a.get("segments"),
                                        cfg["stack"]["block_hours"], a.get("blocks")))
    q = g("station_qc")
    if q:
        lim = q.get("limits")
        L.append("diagnosis: every sub-stack of a pair against the stack of its other sub-stacks (correlation and "
                 "time shift), per station the median over its pairs with the nearest stations: %s station "
                 "sub-stacks judged, %s not judged (too few pairs)" % (q.get("evaluated"), q.get("not_judged")))
        if lim:
            L.append("candidates = outliers of the two distributions (median and 1.4826 MAD, Chauvenet's "
                     "criterion for N = %s): atanh(correlation) more than %.2f sigma (%.3f) below the array median of the "
                     "sub-stack, or a time shift beyond %.2f sigma (%.4f s): %s candidates"
                     % (lim["n"], lim["k_corr"], lim["k_corr"] * lim["corr_sigma"], lim["k_shift"],
                        lim["shift_max"], q.get("candidates", 0)))
        else:
            L.append("NOT CHECKED: too few station sub-stacks to set thresholds from their distribution")
    r = g("removed") or {}
    L.append("removed %s station sub-stacks and %s station-hours (your lists%s); product: %s"
             % (r.get("station_blocks", 0), r.get("station_hours", 0),
                " and the confirmed candidates" if cfg["stack"]["exclude_candidates"] else "; candidates are kept",
                _short(g("product"), 80) if isinstance(g("product"), str) else "the stack without them"))
    return L


def _did_reference(cfg, g):
    L = []
    if g("source") == "fj":
        npf = g("pairs_per_frequency") or [0]
        L.append("frequency-Bessel spectrum of the stack, one frequency at a time with the pairs within %s "
                 "wavelengths (%s-%s pairs per frequency), velocities %s-%s km/s in steps of %s"
                 % (cfg["reference"]["fj_rmax_lambda"], min(npf), max(npf), cfg["reference"]["fj_cmin"],
                    cfg["reference"]["fj_cmax"], cfg["reference"]["fj_dc"]))
        L.append("picked the ridge (%s): %s picks used, %.3g-%.3g Hz; fitted a cubic in log f, held constant outside"
                 % (g("mode"), g("picks_used"), *(g("picked_hz") or [float("nan")] * 2)))
    else:
        L.append("read the borrowed curve %s: defined at %s Hz; %s measurement frequencies have no reference"
                 % (g("file"), _short(g("defined_hz")), g("frequencies_without_reference")))
    L.append("set the measurement band %s Hz and the velocity axis %s km/s from the curve"
             % (_short(g("band_hz")), _short(g("velocity_axis_kms"))))
    return L


def _did_disp(cfg, g):
    L = []
    d = cfg["dispersion"]
    m = g("measure") or {}
    L.append("symmetric component of every pair, Gaussian filters (alpha %s) at %s frequencies per octave, "
             "%.4g-%.4g Hz: %s pairs" % (d["alpha"], d["per_octave"], d["fmin_hz"], d["fmax_hz"], m.get("pairs")))
    L.append("group arrival by ridge tracking guided by the reference; phase read at the group arrival (phase "
             "constant %s); cycle count fixed where the path is %s-%s reference wavelengths, then followed in "
             "frequency" % (d["phase_constant"], d["lam_min"], d["anchor_max_lambda"]))
    L.append("delete-a-group jackknife over %s groups of sub-stacks" % m.get("njack"))
    rc = g("rejections") or {}
    L.append("quality control: %s samples accepted by the measurement; rejected samples by reason (a sample can "
             "have several): %s" % (m.get("phase_points"), ", ".join("%s %d" % kv for kv in rc.items() if kv[1]) or "-"))
    q = g("lambda_qc") or {}
    L.append("wavelength rule (path >= %s of the sample's own wavelengths) and runs of >= %s samples: %s -> %s"
             % (d["qc"]["lam_c"], d["qc"]["min_run"], q.get("phase_points_before"), q.get("phase_points")))
    e = g("export") or {}
    L.append("exported %s curves with %s points (periods %.4g-%.4g s; one segment per pair: %s points of other "
             "segments left out)" % (e.get("phase_curves"), e.get("phase_points"), cfg["export"]["tmin_s"],
                                     cfg["export"]["tmax_s"], e.get("points_dropped_by_one_segment")))
    return L


def _did_inversion(cfg, g):
    L = []
    p, st, gr = g("project") or {}, g("start") or {}, g("grid") or {}
    L.append("average curve of the array (periods with >= %s pairs), 1-D inversion for the starting model"
             % cfg["inversion"]["median_min_pairs"])
    L.append("main project: %s curves, %s travel times, %s periods (%.4g-%.4g s), %s depth nodes down to %s km"
             % (p.get("curves"), p.get("travel_times"), len(p.get("periods") or []),
                min(p.get("periods") or [0]), max(p.get("periods") or [0]), len(p.get("depths") or []),
                (p.get("depths") or [None])[-1]))
    L.append("grid %s x %s nodes, spacing %s x %s deg; weight %s, damp %s; Vs bounds %s-%s km/s"
             % (gr.get("nx"), gr.get("ny"), gr.get("dvxd"), gr.get("dvzd"), p.get("weight"), p.get("damp"),
                st.get("vmin"), st.get("vmax")))
    L.append("check projects built: %s" % (", ".join(n for n, _ in (g("checks", "built") or [])) or "none"))
    sp = g("short_paths") or {}
    if sp:
        L.append("shortest pair %s km%s" % (sp.get("shortest_pair_km"), "; " + sp["note"] if sp.get("note") else ""))
    return L


def _did_run(cfg, g):
    L = []
    res = g("results") or {}
    L.append("precheck of every project (DSurfTomo's first-iteration rejection rule on the starting model), then "
             "%s: %s projects finished, %s invalid (not run), %s failed"
             % (cfg["inversion"]["solver"], sum(1 for v in res.values() if v == "ok"), len(g("invalid") or {}),
                len(g("failed") or [])))
    nz = g("checkerboard_noise") or {}
    if nz:
        L.append("Gaussian noise of the checkerboard travel times: %.2f %% of each time (%s)%s"
                 % (100 * nz["value"], nz["reason"], "; the checkerboards ran after the main project"
                    if auto.noise_is_auto(cfg) else ""))
    return L


def _did_evaluate(cfg, g):
    import numpy as np
    L = []
    f = g("residuals", "final") or {}
    L.append("travel-time residuals of the final model: %s times, rms %.2f ms (%.2f %%), %.0f %% kept by the "
             "rejection rule" % (f.get("n"), f.get("rms_ms", float("nan")), f.get("rel_rms_percent", float("nan")),
                                 100 * f.get("kept_fraction", float("nan"))))
    s0 = next((r for r in (g("residuals", "iterations") or []) if r.get("model") == "start"), {})
    if s0.get("rms_ms") is not None:
        L.append("the same for the starting model: rms %.2f ms (%.2f %%): the 3-D inversion lowered the residual "
                 "by %.1f %%" % (s0["rms_ms"], s0.get("rel_rms_percent", float("nan")),
                                 100 * (1 - f.get("rms_ms", float("nan")) / s0["rms_ms"])))
    for name, rows in sorted((g("checkerboards") or {}).items()):
        ok = [r for r in rows if r.get("slope") is not None and np.isfinite(r["slope"])]
        if ok:
            best = max(ok, key=lambda r: r["slope"])
            L.append("checkerboard %s: recovered amplitude / input amplitude (regression slope) %.2f at best (%.2f "
                     "km depth), %.2f-%.2f over the depths; pattern correlation %.2f-%.2f.  An anomaly of this "
                     "size in the model is damped by about that factor"
                     % (name, best["slope"], best["depth_km"], min(r["slope"] for r in ok),
                        max(r["slope"] for r in ok), min(r["r"] for r in ok), max(r["r"] for r in ok)))
    for name, c in sorted((g("controls") or {}).items()):
        if name.startswith("start_"):
            rc, rm = c.get("residual_std_ms") or [], c.get("main_residual_std_ms") or []
            L.append("control %s (another starting model): final models differ by %.3f km/s rms, %.3f km/s at "
                     "most%s" % (name, c.get("rms_km_s", float("nan")), c.get("max_abs_km_s", float("nan")),
                                 ".  Weighted residual (standard deviation) of the control from the first "
                                 "to the last iteration: %.1f to %.1f ms; of the main inversion: %.1f to %.1f ms "
                                 "(a control that still fits worse has not finished converging: part of the "
                                 "difference is then the starting model)" % (rc[0], rc[-1], rm[0], rm[-1])
                                 if rc and rm else ""))
        else:
            L.append("control %s: final models differ by %.3f km/s rms, %.3f km/s at most"
                     % (name, c.get("rms_km_s", float("nan")), c.get("max_abs_km_s", float("nan"))))
    return L


def _did_figures(cfg, g):
    L = []
    L.append("drew %s figures; not drawn: %s" % (len(g("drawn") or {}), _short(g("not_drawn") or "none", 100)))
    return L


_DID = {"clean": _did_clean, "cc": _did_cc, "stack": _did_stack, "reference": _did_reference, "disp": _did_disp,
        "inversion": _did_inversion, "run": _did_run, "evaluate": _did_evaluate, "figures": _did_figures}


def _dig(d, *keys):
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    return d


def param_rows(cfg, keys):
    """[(key, value, kind, where the value comes from, what it is)] for the report and `lant plan`."""
    from . import auto, params
    dec = auto.load(cfg).get("parameters", {})
    later = {".".join(x) for x in auto.KEYS} | {"inversion.depth_sequence_km", "inversion.vmin", "inversion.vmax"}
    rows = []
    for k in keys:
        v = _get(cfg, k)
        kind, what = params.describe(k)
        src = source(cfg, k)
        if src == "config":
            how = "written in the configuration"
        elif src.startswith("confirmed"):
            how = "confirmed by you (%s)" % src[11:-1]
        elif src.startswith("auto: "):
            how = "derived: " + dec[k]["reason"]
        elif C.is_auto(v) and k in later:
            how = "not derived yet (decided from the data when the chain reaches it)"
        else:
            how = "default"
        rows.append((k, v, kind, how, what))
    return rows


def _param_lines(rows):
    L = []
    for k, v, kind, how, what in rows:
        L.append("  %-40s = %s" % (k, _short(v, 60)))
        L.append("      %-11s %s" % ("[%s]" % kind, what))
        if how not in ("default",):
            L.append("      %-11s %s" % ("", how))
    return L


def step_report(cfg, step, summary, rec=None):
    """Text of the report of one step; printed and written to .lant/steps/<step>.md."""
    out_dir = cfg["project"]["output_dir"]
    cf = cfg.get("_config_file", "CONFIG")
    L = ["", "=" * 100, "step %s: report  (.lant/steps/%s.md)" % (step, step), "-" * 100]
    done = did(cfg, step, summary)
    if done:
        L.append("WHAT THIS STEP DID")
        L += ["  %d. %s" % (i + 1, t) for i, t in enumerate(done)]
    nums = _numbers(summary)
    if nums:
        L.append("KEY NUMBERS (all of them: .lant/summary/%s.json)" % step)
        L += ["  %-46s %s" % (k, _short(v)) for k, v in nums]
    keys = list(used_params(cfg, step))
    L.append("PARAMETERS THIS STEP USED  (value; [what it rests on] meaning; where the value comes from)")
    L += _param_lines(param_rows(cfg, keys))
    look = [f for f in LOOK_AT.get(step, []) if os.path.exists(paths.P(cfg, *f.rstrip("/").split("/")))]
    if look:
        L.append("LOOK AT")
        L += ["  " + paths.P(cfg, *f.rstrip("/").split("/")) for f in look]
    items = (rec or {}).get("to_confirm") or []
    if items:
        from . import confirm as CF
        L.append("FOR YOU TO JUDGE: the chain stops here until you confirm")
        L += ["  - " + c for c in items]
        prop = CF.proposals(cfg, step)
        if prop:
            L.append("  `lant confirm %s -c %s` accepts these values (or write your own in the configuration and run "
                     "the step again):" % (step, cf))
            L += ["    %-42s %s" % (k, _short(v)) for k, v in prop.items()]
        else:
            L.append("  when it is right: lant confirm %s -c %s" % (step, cf))
    elif step in C.CONFIRM_STEPS:
        L.append("FOR YOU TO JUDGE: nothing; the chain goes on")
    tune = TUNE.get(step) or []
    if tune:
        L.append("TO CHANGE THE RESULT")
        L += ["  %-46s %s" % k for k in tune]
        L.append("  then: lant run all -c %s   (it runs this step and every later one again)" % cf)
    L.append("=" * 100)
    return "\n".join(L)


def write_step_report(cfg, step, summary, rec=None):
    txt = step_report(cfg, step, summary, rec)
    d = os.path.join(records.rec_dir(cfg), "steps")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "%s.md" % step), "w") as fh:
        fh.write("```\n" + txt.strip("\n") + "\n```\n")
    if cfg.get("_verbose"):
        import sys
        print(txt, file=sys.stderr, flush=True)
    return txt


def plan(cfg):
    """What the chain will do with this configuration, step by step, before anything runs."""
    from . import confirm as CF
    from .pipeline import STEPS
    L = ["plan for %s  (array type %s, output %s)" % (cfg.get("_config_file", "CONFIG"), cfg["project"]["array_type"],
                                                    cfg["project"]["output_dir"]),
         "values marked 'not derived yet' are decided from the data by a step of the chain; every "
         "step writes its report to .lant/steps/<step>.md", ""]
    for i, step in enumerate(STEPS):
        L.append("%d. %s - %s" % (i + 1, step, STEP_DOC[step]))
        if step in CF.CONFIRM_STEPS:
            L.append("   stops for your judgement: %s" % JUDGED[step])
        L += ["   " + x[2:] for x in _param_lines(param_rows(cfg, list(used_params(cfg, step))))
              if True]
        L.append("")
    return "\n".join(L)


STEP_DOC = {
    "clean": "read every raw file, match the stations, put the records on the UTC sample grid, write day files; "
             "derive the parameters that follow from the geometry and the velocity range",
    "cc": "hourly cross-coherence of every station pair",
    "stack": "stack the hours (sub-stacks for the jackknife); diagnose station sub-stacks that differ from the rest",
    "reference": "reference phase-velocity curve of the array (F-J spectrum) or a borrowed curve; measurement band",
    "disp": "two-station phase velocities with quality control; export of the curves",
    "inversion": "average curve, 1-D starting model, DSurfTomo project, checkerboard and control projects",
    "run": "DSurfTomo / pDSurfTomo on the main project and the check projects",
    "evaluate": "residuals, checkerboard recovery, control differences",
    "figures": "the standard figure set",
}
JUDGED = {
    "clean": "only when station names, coordinates or channel orientations could not be settled from the data",
    "stack": "only when the diagnosis lists candidate station sub-stacks",
    "reference": "do the F-J picks follow the fundamental mode? (not asked for a borrowed curve)",
    "inversion": "starting model, deepest depth node, velocity bounds",
}
