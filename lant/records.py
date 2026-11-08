"""Records kept in the user's output directory (never in the source tree):

    <output>/records/config.yaml            the full configuration of the latest run (defaults filled in)
    <output>/records/config_history/        earlier versions, kept whenever the configuration changed
    <output>/records/run_log.jsonl          one line per step: start, end, status, duration, outputs, messages
    <output>/records/summary/<step>.json    the key numbers of each step (what `lant inspect` reports)
    <output>/records/decisions.md           parameter choices and their reasons (written by the user or the AI)
"""
import datetime as _dt
import hashlib
import json
import os
import shutil
import time
import traceback

from . import config as C

DECISIONS_TEMPLATE = """# Parameter decisions

One entry per decision: what was chosen, from which evidence (file or figure in this output directory), what
the alternatives were, and who decided (user / AI).  Keep earlier entries; add a new one when a value changes.

| date | step | parameter | value | evidence | decided by |
|---|---|---|---|---|---|
"""


def rec_dir(cfg):
    d = os.path.join(cfg["project"]["output_dir"], "records")
    os.makedirs(os.path.join(d, "summary"), exist_ok=True)
    return d


def _now():
    return _dt.datetime.now().isoformat(timespec="seconds")


def save_config(cfg):
    d = rec_dir(cfg)
    cur = os.path.join(d, "config.yaml")
    tmp = cur + ".tmp"
    C.dump(cfg, tmp)
    if os.path.exists(cur):
        if open(cur, "rb").read() == open(tmp, "rb").read():
            os.remove(tmp)
            return cur
        hist = os.path.join(d, "config_history")
        os.makedirs(hist, exist_ok=True)
        stamp = _dt.datetime.fromtimestamp(os.path.getmtime(cur)).strftime("%Y%m%dT%H%M%S")
        shutil.move(cur, os.path.join(hist, "config_%s.yaml" % stamp))
    os.replace(tmp, cur)
    dec = os.path.join(d, "decisions.md")
    if not os.path.exists(dec):
        open(dec, "w").write(DECISIONS_TEMPLATE)
    return cur


def config_hash(cfg, sections):
    sub = {s: cfg[s] for s in sections}
    return hashlib.sha1(json.dumps(sub, sort_keys=True, default=str).encode()).hexdigest()[:12]


def log(cfg, **entry):
    entry.setdefault("time", _now())
    with open(os.path.join(rec_dir(cfg), "run_log.jsonl"), "a") as fh:
        fh.write(json.dumps(entry, default=str, ensure_ascii=False) + "\n")


def write_summary(cfg, step, data):
    p = os.path.join(rec_dir(cfg), "summary", "%s.json" % step)
    with open(p, "w") as fh:
        json.dump(data, fh, indent=1, default=_json_default, ensure_ascii=False)
    return p


def read_summary(cfg, step):
    p = os.path.join(rec_dir(cfg), "summary", "%s.json" % step)
    return json.load(open(p)) if os.path.exists(p) else None


def _json_default(o):
    import numpy as np
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


class Step:
    """Context manager around one pipeline step: logs start / end / failure with the duration."""

    def __init__(self, cfg, name):
        self.cfg, self.name = cfg, name

    def __enter__(self):
        save_config(self.cfg)
        self.t = time.time()
        log(self.cfg, step=self.name, event="start", config=self.cfg.get("_config_file"))
        print("[lant] step %s ..." % self.name, flush=True)
        return self

    def __exit__(self, et, ev, tb):
        dt = round(time.time() - self.t, 1)
        if et is None:
            log(self.cfg, step=self.name, event="done", seconds=dt)
            print("[lant] step %s done (%.0f s)" % (self.name, dt), flush=True)
        else:
            log(self.cfg, step=self.name, event="failed", seconds=dt, error="%s: %s" % (et.__name__, ev),
                traceback=traceback.format_exception(et, ev, tb)[-3:])
            print("[lant] step %s FAILED after %.0f s: %s" % (self.name, dt, ev), flush=True)
        return False
