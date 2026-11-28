"""The configuration file: defaults, loading and checks.

The defaults are the values in config_template.yaml, which is also the documented list of every setting: a key that
is not in the template is an error (a typing mistake cannot silently fall back to a default).  Relative paths are
taken relative to the configuration file.  Each step checks only the keys it uses (`require`), so a configuration
can be filled in step by step.
"""

import copy
import difflib
import os

import yaml

TEMPLATE = os.path.join(os.path.dirname(__file__), "config_template.yaml")
PATH_KEYS = [
    ("output_dir",),
    ("data", "raw_dir"),
    ("data", "station_file"),
    ("data", "response"),
    ("disp", "reference_file"),
    ("inv", "start_model_file"),
]


class ConfigError(ValueError):
    pass


def template_text():
    with open(TEMPLATE) as fh:
        return fh.read()


def defaults():
    return yaml.safe_load(template_text())


def _merge(base, over, where=""):
    if not isinstance(over, dict):
        raise ConfigError("'%s' must be a mapping (key: value)" % (where.rstrip(".") or "the file"))
    for k, v in over.items():
        if k not in base:
            close = difflib.get_close_matches(str(k), [str(b) for b in base], n=1, cutoff=0.7)
            raise ConfigError(
                "unknown key '%s%s'%s (the template lant/config_template.yaml lists every key)"
                % (where, k, "; did you mean '%s'?" % close[0] if close else "")
            )
        if isinstance(base[k], dict):
            _merge(base[k], v, where + k + ".")
        else:
            base[k] = v
    return base


def get(cfg, key):
    """cfg value of a dotted key such as 'clean.fs_hz'."""
    v = cfg
    for part in key.split("."):
        v = v[part]
    return v


def load(path):
    """Configuration file -> dict (defaults filled in, paths absolute)."""
    if not os.path.exists(path):
        raise ConfigError("configuration file %s not found (write one with `lant init`)" % path)
    with open(path) as fh:
        try:
            user = yaml.safe_load(fh) or {}
        except yaml.YAMLError as exc:
            raise ConfigError("%s is not valid YAML: %s" % (path, exc)) from exc
    cfg = _merge(defaults(), user)
    base = os.path.dirname(os.path.abspath(path))
    for key in PATH_KEYS:
        d = cfg
        for part in key[:-1]:
            d = d[part]
        if d[key[-1]] is not None:
            d[key[-1]] = os.path.normpath(os.path.join(base, os.path.expanduser(str(d[key[-1]]))))
    cfg["_file"] = os.path.abspath(path)
    return cfg


def require(cfg, *keys):
    """Stop with one message listing every required key that is still empty."""
    missing = [k for k in keys if get(cfg, k) is None]
    if missing:
        raise ConfigError(
            "%s: fill in %s (the comments in the file say what each one is; `lant init` suggests values)"
            % (cfg.get("_file", "configuration"), ", ".join(missing))
        )


def _check_band(cfg):
    f0, f1 = float(cfg["cc"]["fmin_hz"]), float(cfg["cc"]["fmax_hz"])
    if not 0 < f0 < f1:
        raise ConfigError("cc.fmin_hz must be positive and below cc.fmax_hz")
    return f0, f1


def check_fs(fs, fmax, where="clean.fs_hz"):
    if fmax > 0.4 * fs:
        raise ConfigError(
            "cc.fmax_hz %g Hz is above 0.4 x %s = %g Hz: the anti-alias low-pass (0.8 x Nyquist) cuts it"
            % (fmax, where, 0.4 * fs)
        )


def check_clean(cfg):
    """Settings of lant clean; the correlation band is needed only to remove responses (the pre-filter)."""
    require(cfg, "output_dir", "data.raw_dir", "data.station_file", "clean.fs_hz")
    if float(cfg["clean"]["fs_hz"]) <= 0:
        raise ConfigError("clean.fs_hz must be positive")
    if cfg["data"]["response"]:
        require(cfg, "cc.fmin_hz", "cc.fmax_hz")
        check_fs(float(cfg["clean"]["fs_hz"]), _check_band(cfg)[1])


def check_cc(cfg):
    """Consistency of the correlation settings, with the reason of each rule (the sampling rate is that of the
    cleaned days, checked by lant cc)."""
    require(cfg, "output_dir", "cc.fmin_hz", "cc.fmax_hz", "cc.maxlag_s", "cc.segment_s")
    c = cfg["cc"]
    f0, _ = _check_band(cfg)
    lag, seg = float(c["maxlag_s"]), float(c["segment_s"])
    if seg < 2 * lag:
        raise ConfigError(
            "cc.segment_s %g s is shorter than 2 x cc.maxlag_s: the correlation of a segment computed with the FFT "
            "is circular and its lags would wrap around" % seg
        )
    if seg * f0 < 2:
        raise ConfigError("cc.segment_s %g s holds fewer than 2 periods of cc.fmin_hz" % seg)
    if not 0 <= float(c["overlap"]) < 1:
        raise ConfigError("cc.overlap must be in [0, 1)")
    if 24 % int(c["substack_hours"]) or int(c["substack_hours"]) < 1:
        raise ConfigError("cc.substack_hours must divide 24 (1, 2, 3, 4, 6, 8, 12 or 24)")


def check_disp(cfg):
    require(cfg, "output_dir", "data.velocity_range_kms")
    v = cfg["data"]["velocity_range_kms"]
    if not (isinstance(v, (list, tuple)) and len(v) == 2 and 0 < float(v[0]) < float(v[1])):
        raise ConfigError("data.velocity_range_kms must be [slowest, fastest] in km/s")
    d = cfg["disp"]
    if d["reference"] not in ("fj", "model", "file"):
        raise ConfigError("disp.reference must be fj, model or file")
    if d["reference"] != "fj" and not d["reference_file"]:
        raise ConfigError("disp.reference = %s needs disp.reference_file" % d["reference"])
    if not 0 < float(d["group_window"]) < 1:
        raise ConfigError("disp.group_window must be between 0 and 1 (a fraction of the reference group velocity)")
    if not 0 < float(d["reference_error"]) <= 0.5:
        raise ConfigError("disp.reference_error must be in (0, 0.5]")
    if int(d["per_octave"]) < 1 or int(d["min_points"]) < 2:
        raise ConfigError("disp.per_octave must be >= 1 and disp.min_points >= 2")


def snapshot(cfg, sections):
    """The parts of the configuration a step used (for its params.yaml)."""
    out = {}
    for s in sections:
        out[s] = copy.deepcopy(cfg[s])
    return out
