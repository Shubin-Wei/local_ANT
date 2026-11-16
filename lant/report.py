"""Step reports and the plan: what a step did (in order, with counts), the key numbers, every parameter it used
with its value, what the value rests on (lant.params) and where it comes from (written by the user, confirmed,
derived by a rule with the numbers of the present data, or the default), what to look at, what is left for the user
to judge, and what to change - kept in <output>/.lant/steps/<step>.md (printed with -v).
`lant plan` shows the same parameter table for every step before anything runs.

Nothing here changes a result.
"""

import copy
import json
import os

from . import auto, paths, records
from . import config as C
from .step_descriptions import _short, did

# The configuration each step reads (sections or single keys).  A change in any of them makes the step and every
# later step out of date.
STEP_KEYS = {
    "clean": [
        "project.array_type",
        "data",
        "clean",
        "correlation.fs_max",
        "correlation.fs",
        "correlation.fmin_hz",
        "correlation.fmax_hz",
        "correlation.maxlag_s",
        "correlation.seg_len_s",
        "dispersion.alpha",
        "dispersion.per_octave",
        "dispersion.fmax_hz",
        "stack.block_hours",
        "inversion.grid_spacing_deg",
    ],
    "cc": ["correlation"],
    "stack": ["stack"],
    "reference": [
        "reference",
        "dispersion.fmin_hz",
        "dispersion.fmax_hz",
        "dispersion.per_octave",
        "dispersion.vmin",
        "dispersion.vmax",
        "dispersion.dv",
        "dispersion.qc.lam_c",
    ],
    "disp": ["dispersion", "export"],
    "inversion": ["inversion", "checks", "export.tmin_s", "export.tmax_s"],
    "run": ["inversion.noise", "checks.precheck"],
    "evaluate": ["checks"],
    "figures": ["plots"],
}

# Small input files named in the configuration: their content, not only their path, belongs to the step's
# parameters, so that a file rewritten in place makes the step out of date (the raw records, station table and
# responses are checked by lant.manifest.raw_inputs).
CONTENT_KEYS = ("reference.file", "inversion.start_model.file")


def content_digest(path):
    import hashlib

    try:
        with open(path, "rb") as fh:
            return hashlib.sha256(fh.read()).hexdigest()
    except OSError as exc:
        return "unreadable: %s" % exc.strerror


# The parameters one usually changes after looking at a step, and what they do.
TUNE = {
    "clean": [
        ("data.velocity_range_kms", "the approximate velocities: lag window, band, F-J scan"),
        ("data.exclude_stations", "leave stations out (instrument or coupling problems)"),
        ("data.rename", "match header names to the station table by hand"),
        ("data.channel", "choose the channel when the names do not tell the vertical one"),
        ("data.start_date / end_date", "process only part of the record"),
        ("correlation.fs_max", "the highest processing rate, which sets the top of the band"),
    ],
    "cc": [
        ("correlation.max_over_std", "segment rejection (transients)"),
        ("correlation.remove_1s_comb", "the 1-s timing-pulse removal (Takagi et al. 2015)"),
    ],
    "stack": [
        ("lant confirm stack --exclude-candidates", "remove the candidate station sub-stacks"),
        ("stack.exclude_station_days", 'station-days of your own, also ranges [station, "YYYYDDD-"]'),
    ],
    "reference": [
        ("reference.fj_mode", "max_peak / track_fundamental when the picks jump to another ridge"),
        ("reference.track_start_band_hz", "where the fundamental ridge is unambiguous"),
        ("reference.fj_fit_fmax", "cap the fitted picks"),
        ("reference.source: file", "a borrowed reference curve"),
        ("dispersion.fmin_hz / fmax_hz", "the measured band"),
    ],
    "disp": [("export.tmin_s / tmax_s", "periods passed on")],
    "inversion": [
        ("inversion.depth_sequence_km", "depth nodes"),
        ("inversion.weight", "regularisation"),
        ("inversion.start_model", "a starting model of your own"),
        ("checks.*", "checkerboards and controls"),
    ],
    "run": [("inversion.weight", "regularisation (then from the inversion step)")],
    "evaluate": [("checks.iteration_residuals", "residuals of every iteration (one forward run each)")],
    "figures": [("plots.*", "colour range, slice depths, sections, topography, font size")],
}

LOOK_AT = {
    "clean": [
        "work/01_clean/report.txt",
        "results/figures/0_data_coverage.png",
        "work/01_clean/stations.csv",
        "work/01_clean/processing_days.txt",
    ],
    "cc": ["work/02_ccf/"],
    "stack": [
        "work/03_stack/candidates.txt",
        "results/figures/2_station_qc.png",
        "work/03_stack/station_qc/station_qc.txt",
    ],
    "reference": ["results/figures/3_fj_phase_velocity.png", "work/04_reference/fj.txt"],
    "disp": ["work/05_dispersion/pair_figures/", "work/05_dispersion/median_vs_reference.txt"],
    "inversion": [
        "work/06_inversion/average_curve.txt",
        "work/06_inversion/vs1d_fit.txt",
        "work/06_inversion/vs1d.txt",
    ],
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
    if cfg.get("_frozen"):
        new = copy.deepcopy(cfg)
        if without:
            new.get("_confirmed", {}).pop(without, None)
        return auto.apply(new)
    if not cfg.get("_config_file"):
        return cfg
    new = C.load(cfg["_config_file"])
    if without:
        new["_confirmed"].pop(without, None)
    new["_verbose"] = cfg.get("_verbose", False)
    new["_position"] = cfg.get("_position", "")
    return auto.apply(new)


def used_params(cfg, step):
    """{dotted key: value} of the configuration the step reads (after the "auto" values are resolved), in the form
    the run log keeps (JSON: dates as text)."""
    out = {}
    for k in STEP_KEYS.get(step, []):
        value = cfg.get("_plot_settings", cfg["plots"]) if step == "figures" and k == "plots" else _get(cfg, k)
        _flatten(k, value, out)
    for k, setting in auto.DERIVED_LATER.get(step, {}).items():
        if k in out:
            out[k] = setting(cfg)
    for k in CONTENT_KEYS:
        if isinstance(out.get(k), str):
            out[k + ".sha256"] = content_digest(out[k])
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


def _numbers(summary, n=14, root=None):
    """A few top-level numbers of the step summary; keys that hold a path of the project are shown relative to its
    root."""
    flat = _flatten("", summary or {}, {})
    if root:
        prefix = os.path.abspath(root).rstrip(os.sep) + os.sep
        flat = {k.replace(prefix, ""): v for k, v in flat.items()}
    out = []
    skip = (
        "day_median",
        "auto_parameters",
        "parameters.",
        "messages",
        "score_percentiles",
        "rms_by_depth",
        "drawn.",
        "not_drawn.",
        "plot_parameters",
        "resources.",
    )
    for k, v in flat.items():
        if any(x in k for x in skip) or k.count(".") > 2:
            continue
        if isinstance(v, bool) or isinstance(v, (int, float)) or (isinstance(v, str) and len(v) < 60):
            out.append((k, v))
        if len(out) >= n:
            break
    return out


def resource_lines(r):
    """COMPUTING RESOURCES of a step: the cores and the memory it used (summary["resources"], lant.resources)."""
    if not r:
        return []
    L = []
    start = r.get("at_start") or {}
    if start:
        L.append("cores at the start: %d (%s)" % (start["cores"], start["rule"]))
    for p in r.get("plans") or []:
        mem = "memory not estimated"
        if p.get("memory_each_estimate_gb") is not None:
            mem = "estimated %.2f GB each, %.0f %% of %.1f GB available%s" % (
                p["memory_each_estimate_gb"], 100 * p["memory_fraction"], p["memory_available_gb"] or 0,
                "; fewer at a time because of the memory" if p["limited_by_memory"] else "")  # fmt: skip
        L.append("%s: %d at a time of %d, %d thread(s) each, from %d cores; %s"
                 % (p["what"], p["parallel"], p["tasks"], p["threads_each"], p["cores"], mem))  # fmt: skip
    for q in r.get("solver_projects") or []:
        L.append("solver %s: %d thread(s), %.1f s, peak memory %s GB (estimate %s GB)"
                 % (q["project"], q["threads"], q["wall_seconds"], q["peak_memory_gb"], q["memory_estimate_gb"]))  # fmt: skip
    if r.get("worker_peak_gb") is not None:
        L.append("peak memory of one worker process: %.2f GB" % r["worker_peak_gb"])
    if r.get("main_peak_gb") is not None:
        L.append("peak memory of the lant process during the step: %.2f GB" % r["main_peak_gb"])
    return L


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
    cf = cfg.get("_config_file", "CONFIG")
    L = ["", "=" * 100, "step %s: report  (.lant/steps/%s.md)" % (step, step), "-" * 100]
    done = did(cfg, step, summary)
    if done:
        L.append("WHAT THIS STEP DID")
        L += ["  %d. %s" % (i + 1, t) for i, t in enumerate(done)]
    nums = _numbers(summary, root=cfg["project"]["output_dir"])
    if nums:
        L.append("KEY NUMBERS (all of them: .lant/summary/%s.json)" % step)
        L += ["  %-46s %s" % (k, _short(v)) for k, v in nums]
    res = resource_lines((summary or {}).get("resources"))
    if res:
        L.append("COMPUTING RESOURCES")
        L += ["  " + t for t in res]
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
            L.append(
                "  `lant confirm %s -c %s` accepts these values (or write your own in the configuration and run "
                "the step again):" % (step, cf)
            )
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


def plan(cfg, steps=None):
    """What the chain will do with this configuration, step by step, before anything runs."""
    from . import confirm as CF
    from .pipeline import STEPS

    L = [
        "plan for %s  (array type %s, output %s)"
        % (cfg.get("_config_file", "CONFIG"), cfg["project"]["array_type"], cfg["project"]["output_dir"]),
        "values marked 'not derived yet' are decided from the data by a step of the chain; every "
        "step writes its report to .lant/steps/<step>.md",
        "",
    ]
    from . import manifest

    previous = manifest.load(cfg)["steps"]
    for i, step in enumerate(steps or STEPS):
        L.append("%d. %s - %s" % (i + 1, step, STEP_DOC[step]))
        if step in CF.CONFIRM_STEPS:
            L.append("   stops for your judgement: %s" % JUDGED[step])
        L += ["   " + x[2:] for x in _param_lines(param_rows(cfg, list(used_params(cfg, step)))) if True]
        old = previous.get(step)
        if old:
            changed = [
                (k, old["params"].get(k), v)
                for k, v in used_params(cfg, step).items()
                if old["params"].get(k) != v and k not in manifest.RESOURCE_KEYS
            ]
            L.append(
                "   compared with recorded run %s: %s"
                % (old["run"], "changed (existing results keep their old settings)" if changed else "unchanged")
            )
            L += ["     %s: %s -> %s" % (k, _short(a), _short(b)) for k, a, b in changed]
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
    "run": "the Python DSurfTomo on the main project and the check projects",
    "evaluate": "residuals, checkerboard recovery, control differences",
    "figures": "the standard figure set",
}
JUDGED = {
    "clean": "only when station names, coordinates or channel orientations could not be settled from the data",
    "stack": "only when the diagnosis lists candidate station sub-stacks",
    "reference": "do the F-J picks follow the fundamental mode? (not asked for a borrowed curve)",
    "inversion": "starting model, deepest depth node, velocity bounds",
}
