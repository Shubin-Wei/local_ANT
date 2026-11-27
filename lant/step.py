"""What every step leaves in its folder besides its results, so that nothing it did is hidden:

params.yaml   the settings it used (the sections of the configuration) and the values it derived, each with the
              rule that gave it
log.txt       everything it printed, with the time
summary.txt   the key numbers at the end
figures/      its figures
"""

import datetime
import os
import sys
import time

import yaml

from . import __version__


class _Tee:
    def __init__(self, stream, fh):
        self.stream, self.fh = stream, fh

    def write(self, text):
        self.stream.write(text)
        self.fh.write(text)
        return len(text)

    def flush(self):
        self.stream.flush()
        self.fh.flush()


class Step:
    """with Step(cfg, "1_cc", ["data", "cc"]) as s: ...  s.dir, s.fig(name), s.derive(key, value, rule), s.say()"""

    def __init__(self, cfg, folder, sections):
        self.cfg = cfg
        self.dir = os.path.join(cfg["output_dir"], folder)
        self.sections = sections
        self.derived = {}
        self.numbers = {}

    def __enter__(self):
        os.makedirs(os.path.join(self.dir, "figures"), exist_ok=True)
        self._log = open(os.path.join(self.dir, "log.txt"), "a")
        self._out = sys.stdout
        sys.stdout = _Tee(sys.stdout, self._log)
        self.t0 = time.time()
        self.say("lant %s %s  (configuration %s)" % (__version__, os.path.basename(self.dir), self.cfg["_file"]))
        self._write_params()
        return self

    def __exit__(self, kind, value, tb):
        dt = time.time() - self.t0
        if kind is None:
            self.numbers["run time (s)"] = round(dt, 1)
            with open(os.path.join(self.dir, "summary.txt"), "w") as fh:
                for k, v in self.numbers.items():
                    fh.write("%s: %s\n" % (k, v))
            self.say("done in %.0f s; results in %s" % (dt, self.dir))
        else:
            self.say("STOPPED after %.0f s: %s" % (dt, value))
        sys.stdout = self._out
        self._log.close()
        return False

    def say(self, text):
        stamp = datetime.datetime.now().strftime("%H:%M:%S")
        print("[%s] %s" % (stamp, text), flush=True)

    def fig(self, name):
        return os.path.join(self.dir, "figures", name)

    def path(self, *parts):
        return os.path.join(self.dir, *parts)

    def derive(self, key, value, rule):
        """A value the step computed from the data: kept in params.yaml with its rule, and printed."""
        self.derived[key] = dict(value=value, rule=rule)
        self.say("%s = %s   (%s)" % (key, value, rule))
        self._write_params()
        return value

    def number(self, key, value):
        self.numbers[key] = value

    def _write_params(self):
        out = {s: self.cfg[s] for s in self.sections}
        out["derived"] = self.derived
        out["lant_version"] = __version__
        with open(os.path.join(self.dir, "params.yaml"), "w") as fh:
            yaml.safe_dump(_plain(out), fh, sort_keys=False, allow_unicode=True)


def _plain(x):
    """numpy scalars and tuples -> plain YAML types."""
    import numpy as np

    if isinstance(x, dict):
        return {k: _plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_plain(v) for v in x]
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, np.ndarray):
        return x.tolist()
    return x


def need(path, step):
    """A result of an earlier step, or a clear message which step to run first."""
    if not os.path.exists(path):
        raise FileNotFoundError("%s not found: run `lant %s` first" % (path, step))
    return path
