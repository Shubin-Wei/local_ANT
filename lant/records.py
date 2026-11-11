"""Records kept in the user's output directory (never in the source tree):

<output>/.lant/config.yaml            the full configuration of the latest run (defaults filled in)
<output>/.lant/config_history/        earlier versions, kept whenever the configuration changed
<output>/.lant/run_log.jsonl          one line per step: start, end, status, duration, outputs, messages
<output>/.lant/summary/<step>.json    the key numbers of each step (what `lant inspect` reports)
<output>/notes.md                     parameter choices and their reasons (appended through lant)
<output>/.lant/provenance/<id>.json   code version and environment of the runs (lant/provenance.py)
<output>/.lant/lock                   the process writing to the output directory
"""

import datetime as _dt
import json
import os
import shutil
import time
import traceback

from . import config as C


def rec_dir(cfg, create=False):
    d = os.path.join(cfg["project"]["output_dir"], ".lant")
    if create:
        os.makedirs(os.path.join(d, "summary"), exist_ok=True)
    return d


def _now():
    return _dt.datetime.now().isoformat(timespec="seconds")


def save_config(cfg):
    d = rec_dir(cfg, create=True)
    cur = os.path.join(d, "config.yaml")
    tmp = cur + ".tmp"
    C.dump(cfg, tmp)
    if os.path.exists(cur):
        if open(cur, "rb").read() == open(tmp, "rb").read():
            os.remove(tmp)
            return cur
        hist = os.path.join(d, "config_history")
        os.makedirs(hist, exist_ok=True)
        stamp = _dt.datetime.fromtimestamp(os.path.getmtime(cur)).strftime("%Y%m%dT%H%M%S_%f")
        dst = os.path.join(hist, "config_%s.yaml" % stamp)
        k = 1
        while os.path.exists(dst):  # never overwrite an earlier version
            dst = os.path.join(hist, "config_%s_%d.yaml" % (stamp, k))
            k += 1
        shutil.move(cur, dst)
    os.replace(tmp, cur)
    return cur


def log(cfg, **entry):
    entry.setdefault("time", _now())
    with open(os.path.join(rec_dir(cfg, create=True), "run_log.jsonl"), "a") as fh:
        fh.write(json.dumps(strict(entry), default=str, ensure_ascii=False, allow_nan=False) + "\n")


def strict(o):
    """JSON-safe copy: numpy types to Python, NaN / inf to None (standard JSON has no NaN)."""
    import math

    import numpy as np

    if isinstance(o, dict):
        return {str(k): strict(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [strict(v) for v in o]
    if isinstance(o, np.ndarray):
        return strict(o.tolist())
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (float, np.floating)):
        return float(o) if math.isfinite(o) else None
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def write_summary(cfg, step, data):
    p = os.path.join(rec_dir(cfg, create=True), "summary", "%s.json" % step)
    with open(p + ".tmp", "w") as fh:
        json.dump(strict(data), fh, indent=1, default=str, ensure_ascii=False, allow_nan=False)
    os.replace(p + ".tmp", p)
    return p


def read_summary(cfg, step):
    p = os.path.join(rec_dir(cfg), "summary", "%s.json" % step)
    return json.load(open(p)) if os.path.exists(p) else None


class StepIncomplete(RuntimeError):
    """A step finished, but some of its inputs failed (unreadable files, ...): the chain stops here."""


class Step:
    """Own the write lock, detailed stdout log and progress for an entire step."""

    def __init__(self, cfg, name):
        self.cfg, self.name = cfg, name

    def __enter__(self):
        import contextlib
        import logging
        import sys
        from logging.handlers import RotatingFileHandler

        from .progress import Progress
        from .provenance import record

        self.context = contextlib.ExitStack()
        try:
            from .artifacts import protect_step_outputs, writing

            self.context.enter_context(writing(self.cfg, self.name))
            protect_step_outputs(self.cfg, self.name)
            from .delivery import initialize

            initialize(self.cfg)
            save_config(self.cfg)
            self.t = time.time()
            log(
                self.cfg,
                step=self.name,
                event="start",
                config=self.cfg.get("_config_file"),
                provenance=record(rec_dir(self.cfg)),
            )
            handler = RotatingFileHandler(
                os.path.join(rec_dir(self.cfg), "lant.log"), maxBytes=10 * 1024 * 1024, backupCount=3, encoding="utf-8"
            )
            self.context.callback(handler.close)
            verbose = self.cfg.get("_verbose", False)

            class Output:
                def write(self, text):
                    if text.strip():
                        handler.emit(logging.LogRecord("lant", logging.INFO, "", 0, text.rstrip(), (), None))
                    if verbose:
                        sys.stderr.write(text)
                    return len(text)

                def flush(self):
                    handler.flush()

            self.context.enter_context(contextlib.redirect_stdout(Output()))
            print("[lant] step %s started" % self.name, flush=True)
            position = self.cfg.get("_position")  # "6/9" within lant run all
            self.context.enter_context(Progress("[%s] %s" % (position, self.name) if position else self.name))
            return self
        except BaseException:
            self.context.close()
            raise

    def __exit__(self, et, ev, tb):
        import sys

        dt = round(time.time() - self.t, 1)
        failed = getattr(self, "failed", None)
        event = "failed" if et else ("partial" if failed else "done")
        extra = dict(error="%s: %s" % (et.__name__, ev), traceback=traceback.format_exception(et, ev, tb)) if et else {}
        if failed:
            extra.update(failed=failed, n_failed=len(failed))
        try:
            log(self.cfg, step=self.name, event=event, seconds=dt, **extra)
            print("[lant] step %s %s (%.0fs)" % (self.name, event, dt), flush=True)
        finally:
            self.context.__exit__(et, ev, tb)  # registers the outputs and releases the lock
        print(
            "[lant] %s %s (%.0fs); details: %s" % (self.name, event, dt, os.path.join(rec_dir(self.cfg), "lant.log")),
            file=sys.stderr,
            flush=True,
        )
        return False
