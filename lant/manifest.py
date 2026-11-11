""".lant/manifest.json: the one record of what each step produced, from what, and whether it still stands.

Per step (its latest run):
    state       done | partial (some inputs failed) | failed
    run         counter, raised every time the step runs
    time, seconds
    params      the configuration values the step read (after the "auto" values were resolved)
    upstream    {earlier step: its run counter when this step ran}
    inputs      clean only: fingerprint of the raw directory listing and of the station table
    code        id of .lant/provenance/<id>.json (code version, packages, third-party programs)
    to_confirm  the scientific judgements the step leaves to the user (empty: nothing to confirm)
    basis       fingerprint of the results those judgements are about; `lant confirm` writes it into the
                configuration, so a confirmation holds exactly as long as the results stay the same
    error       failed / partial: what went wrong

`status` derives the state of every step from this file alone, and `gate` lets a step run only when every earlier
step stands:
    not run       no record
    failed        the last run raised an error
    partial       the last run finished but some inputs failed
    out of date   its parameters changed, its inputs changed (clean), an earlier step was run again since, or an
                  earlier step does not stand
    to confirm    done, with judgements not yet confirmed (or confirmed for other results)
    done          stands
"""
import datetime as _dt
import hashlib
import json
import os

from . import records


# keys that set only the computing resources (processes, threads, task size, disk margin): changing them does not
# change any result, so a step does not become out of date (and its confirmation stands)
RESOURCE_KEYS = {"clean.nproc", "clean.min_free_disk_factor", "correlation.nproc", "correlation.threads_per_process",
                 "dispersion.nproc", "dispersion.chunk", "inversion.threads", "inversion.njobs"}


def _steps():
    from .pipeline import STEPS
    return STEPS


def path(cfg):
    return os.path.join(records.rec_dir(cfg), "manifest.json")


def load(cfg):
    p = path(cfg)
    if os.path.exists(p):
        return json.load(open(p))
    return {"steps": {}}


def _save(cfg, m):
    p = path(cfg)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p + ".tmp", "w") as fh:
        json.dump(records.strict(m), fh, indent=1, default=str, ensure_ascii=False, allow_nan=False)
    os.replace(p + ".tmp", p)


def raw_inputs(cfg):
    """Fingerprint of what the clean step reads: every candidate file under data.raw_dir (relative path, size,
    modification time) and the content of the station table."""
    from . import data as D
    h = hashlib.sha1()
    d = cfg["data"]
    try:
        files = D.find_files(d["raw_dir"], d["file_patterns"], d["ignore_patterns"])
    except Exception as exc:
        return "unreadable: %s" % exc
    for f in sorted(files):
        try:
            st = os.stat(f)
            h.update(repr((os.path.relpath(f, d["raw_dir"]), st.st_size, st.st_mtime_ns)).encode())
        except OSError:
            h.update(repr((f, "missing")).encode())
    try:
        with open(d["station_file"], "rb") as fh:
            h.update(fh.read())
    except OSError:
        h.update(b"no station file")
    return h.hexdigest()[:16]


def record(cfg, step, state, seconds=None, params=None, to_confirm=None, basis=None, error=None, code=None,
           inputs=None):
    """Write the record of a run of `step` (replacing that of its previous run)."""
    m = load(cfg)
    steps = m["steps"]
    order = _steps()
    prev = steps.get(step, {})
    steps[step] = dict(
        state=state, run=int(prev.get("run", 0)) + 1, time=_dt.datetime.now().isoformat(timespec="seconds"),
        seconds=seconds, params=params or {},
        upstream={s: steps[s]["run"] for s in order[:order.index(step)] if s in steps},
        inputs=inputs, code=code, to_confirm=list(to_confirm or []), basis=basis, error=error)
    from .artifacts import required
    steps[step]["outputs"] = required(cfg, step)
    _save(cfg, m)
    return steps[step]


def upstream_changed(cfg, step):
    """True when an earlier step was run again since `step` last ran (or `step` has no record)."""
    steps = load(cfg)["steps"]
    r = steps.get(step)
    if r is None:
        return True
    order = _steps()
    return any(s in steps and r["upstream"].get(s) != steps[s]["run"] for s in order[:order.index(step)])


def confirmed(cfg, step, rec):
    """Is the record's basis the one written into the configuration by `lant confirm`?"""
    if not rec.get("to_confirm"):
        return True
    c = (cfg.get("_confirmed") or {}).get(step) or {}
    return c.get("_basis") == rec.get("basis")


def status(cfg, check_inputs=True):
    """{step: {state, time, reason}} in the order of the steps, from the manifest and the present configuration
    (cfg: resolved, lant.report.resolved)."""
    from . import report
    steps = load(cfg)["steps"]
    order = _steps()
    out = {}
    broken = None                                   # the first earlier step that does not stand
    for s in order:
        r = steps.get(s)
        if r is None:
            out[s] = dict(state="not run")
            if s != "figures":
                broken = broken or (s, "has not run")
            continue
        e = dict(state=r["state"], time=r["time"], run=r["run"])
        if r.get("cache_cleaned"):
            e["cache_cleaned"] = True
        if r["state"] in ("failed", "partial"):
            e["reason"] = r.get("error") or ""
            broken = broken or (s, "failed" if r["state"] == "failed" else "is incomplete")
            out[s] = e
            continue
        now = report.used_params(cfg, s)
        changed = sorted(k for k in set(now) | set(r["params"]) if now.get(k) != r["params"].get(k)
                         and k not in RESOURCE_KEYS)
        again = [u for u in order[:order.index(s)] if u in steps and r["upstream"].get(u) != steps[u]["run"]]
        missing = [p for p in r.get("outputs", []) if not os.path.exists(os.path.join(cfg["project"]["output_dir"], p))]
        why = None
        if missing:
            why = "required products missing: " + ", ".join(missing)
        elif changed:
            why = "parameters changed since it ran: %s%s" % (", ".join(changed[:8]), " ..." if len(changed) > 8 else "")
        elif s == "clean" and check_inputs and r.get("inputs") and raw_inputs(cfg) != r["inputs"]:
            why = "the raw files or the station table changed since it ran"
        elif again:
            why = "step %s was run again since" % again[-1]
        elif broken:
            why = "an earlier step (%s) %s" % broken
        if why:
            e.update(state="out of date", reason=why)
            broken = broken or (s, "is out of date")
        elif not confirmed(cfg, s, r):
            e.update(state="to confirm", reason="; ".join(r["to_confirm"]))
            broken = broken or (s, "waits for confirmation")
        out[s] = e
    return out


def next_action(cfg, st):
    """(command, why) of the next thing to do, or (None, "everything is done")."""
    c = cfg.get("_config_file", "CONFIG")
    for s, e in st.items():
        if e["state"] == "done":
            continue
        if e["state"] == "to confirm":
            return ("lant confirm %s -c %s" % (s, c),
                    "step %s leaves judgements to you: read .lant/steps/%s.md first" % (s, s))
        if e["state"] == "partial":
            return ("lant run all -c %s" % c, "step %s is incomplete (%s): repair or exclude what failed, then run "
                                               "again" % (s, e.get("reason", "")))
        if e["state"] == "failed":
            return ("lant run all -c %s" % c, "step %s failed (%s): fix the cause, then run again" % (s, e.get("reason", "")))
        return ("lant run all -c %s" % c, "step %s %s%s" % (s, "has not run yet" if e["state"] == "not run" else
                                                               "is " + e["state"], ": " + e["reason"] if e.get("reason") else ""))
    return (None, "everything is done")


class NotCurrent(RuntimeError):
    """An earlier step does not stand: the step asked for must not read its outputs."""


def gate(cfg, step, st=None):
    """Raise unless every earlier step stands (figures: earlier steps that never ran are left out)."""
    from .confirm import NeedsConfirmation
    st = st or status(cfg)
    order = _steps()
    c = cfg.get("_config_file", "CONFIG")
    needed = {"cc": "clean", "stack": "cc"}.get(step)
    if needed and st[needed].get("cache_cleaned"):
        raise NotCurrent("%s caches were cleaned; rebuild with lant run all --from clean -c %s" % (needed, c))
    for s in order[:order.index(step)]:
        e = st[s]
        if e["state"] == "done" or (step == "figures" and e["state"] == "not run"):
            continue
        if e["state"] == "to confirm":
            raise NeedsConfirmation(
                "step %s waits for your confirmation of %s:\n  - %s\nread .lant/steps/%s.md, then\n"
                "  lant confirm %s -c %s\n  lant run all -c %s" % (step, s, e["reason"].replace("; ", "\n  - "), s, s, c, c))
        raise NotCurrent("step %s cannot run: step %s is %s%s.  Run\n  lant run all -c %s\n(it continues from the "
                         "first step that does not stand)" % (step, s, e["state"],
                                                              " (%s)" % e["reason"] if e.get("reason") else "", c))
