"""Project configuration: defaults, loading, validation.

A project is described by one YAML file (see config/surface.example.yaml and config/line.example.yaml).  Keys
that are not given take the defaults below; the defaults are the values of the method as published for the
dense nodal array studies (docs/parameters.md gives each value, its basis and when to change it).  A key that
does not exist in the defaults is an error, so that a typing mistake cannot silently fall back to a default.

`array_type: line` applies LINE_OVERRIDES on top of the defaults before the user's file is read.
"""
import copy
import os

import yaml

DEFAULTS = {
    "project": {
        "name": "my_array",
        "array_type": "auto",             # auto (from the station layout) | surface | line
        "output_dir": None,               # required
    },
    "data": {
        "raw_dir": None,                  # required: directory searched recursively; any layout, any file names
        "format": "auto",                 # auto (ObsPy detects MSEED, SAC, GSE2, SEG2, SEISAN, ...) or one ObsPy format name
        "file_patterns": ["*"],           # file-name patterns to try (default: every file)
        "ignore_patterns": ["*.txt", "*.csv", "*.xml", "*.json", "*.yaml", "*.yml", "*.py", "*.sh", "*.png",
                            "*.jpg", "*.pdf", "*.log", "*.md", "*.zip", "*.gz", "*.tar", ".*"],
        "station_file": None,             # required: station table in any common form (lant/stations.py), or StationXML
        "channel": "auto",                # auto: the vertical component of each station; or a pattern such as "*Z"
        "network_default": None,          # optional override: network code for headers without one
        "rename": {},                     # optional override: {"NET.STA of the headers": "NET.STA of the table"}
        "day_start_utc_hour": "auto",     # auto (the hour most records start) | 0-23: processing days start at this UTC hour
        "min_trace_length_s": 0.0,        # pieces of record shorter than this (s) are not used
        "exclude_stations": [],           # station names (NET.STA) left out entirely; or a file name
        "field_note_stations": [],        # stations with a field note (moved / dug out): the persistence rule of stack.station_day_qc
        "start_date": None,               # optional YYYY-MM-DD limits of the processing days
        "end_date": None,
    },
    "clean": {
        "write_archive": True,            # write the cleaned day files (01_clean/archive); False = later steps read the raw files
        "nproc": 4,
        "encoding": "FLOAT32",            # miniSEED encoding of the archive (FLOAT32 keeps integer and float data exactly
                                          # up to 2**24 counts; use INT32 or STEIM2 for integer data to save space)
        "min_free_disk_factor": 1.2,      # refuse to write the archive if free disk < factor x the estimated size
    },
    "correlation": {
        "fs": "auto",                     # Hz after decimation; auto: largest of 50, 40, 25, 20, 10, 5 dividing every raw rate
        "seg_len_s": 60.0,
        "seg_step_s": 30.0,
        "maxlag_s": "auto",               # auto: max distance / 1.2 km/s + 2 longest periods, rounded up to 10 s
        "fmin_hz": 0.2,
        "fmax_hz": "auto",                # auto: min(12.5, 0.4 fs)
        "smooth_bins": 1,                 # spectral smoothing of the coherence denominator (1 = none)
        "max_over_std": 10.0,             # segment rejection: peak / station-day std
        "remove_1s_comb": True,           # hourly removal of a UTC-locked 1-s periodic signal (nodal GPS pulse)
        "nproc": 4,
        "threads_per_process": 2,
    },
    "stack": {
        "block_hours": 24,                # length of the sub-stacks used for the jackknife (24 = days)
        "station_day_qc": {
            "enabled": True,
            "band_hz": [1.0, 4.0],
            "dist_km": [2.0, 10.0],
            "corr_drop": 0.25,
            "shift_s": 0.04,
            "min_partners": 10,
            "min_segments": 200,
            "field_drop": 0.05,
        },
        "station_hour_qc": {
            "enabled": False,
            "band_hz": [1.0, 4.0],
            "maxlag_s": 10.0,
            "max_dist_km": 6.0,
            "min_partners": 3,
            "threshold": 0.4,
        },
        "exclude_station_days": [],       # extra [station, YYYYDDD] entries, or a file "station day" per line
    },
    "notch": {
        "enabled": True,
        "band_hz": [0.5, 12.0],
        "raw_nsta": 20,                   # stations sampled for the raw-waveform line search
        "raw_ndays": 4,                   # days sampled (evenly spaced over the record)
        "raw_days": [],                   # explicit YYYYDDD list instead of raw_ndays
        "raw_min_frac": 0.15,
        "db": 6.0,
        "seed": 0,
        "table_file": None,               # use this notch table instead of the rule (freq_Hz half_width_Hz)
    },
    "reference": {
        "source": "fj",                   # fj: F-J of this array (CC-FJpy) | file: borrowed curve (reference.file)
        "file": None,                     # npz written by lant, or a text file "freq_Hz c_km/s"
        "fj_mode": "max_peak",            # max_peak (areal arrays) | track_fundamental (lines)
        "fj_rmax_lambda": 5.0,
        "fj_c_nominal": 2.5,
        "fj_cmin": 1.2,
        "fj_cmax": 4.0,
        "fj_dc": 0.005,
        "fj_nfft": 8192,
        "fj_fit_fmax": 6.0,
        "fj_threads": 8,
        "track_start_band_hz": [0.6, 1.2],
        "pass2": True,                    # second pass: median of the first-pass phase velocities
        "pass2_min_pairs": 50,
        "pass2_max_gap_octave": 0.25,
    },
    "dispersion": {
        "fmin_hz": "auto",                # auto: p90 distance = 2 wavelengths at 2.5 km/s, down to the 1/12-octave grid
        "fmax_hz": "auto",                # auto: 0.84 x correlation.fmax_hz
        "per_octave": 12,
        "vmin": 1.2,
        "vmax": 4.0,
        "dv": 0.01,
        "alpha": 5.0,
        "phase_constant": "pi/4",
        "min_dist_km": 0.3,
        "lam_min": 1.0,                   # paths shorter than this many reference wavelengths are not measured
        "anchor_max_lambda": 3.5,         # cycle count anchored only where the path is 1.0-3.5 wavelengths
        "njack": 5,                       # jackknife groups; "auto" = blocks with >= 25 % of the fullest block
        "min_jackknife_groups": 3,
        "nproc": 8,
        "chunk": 20,
        "side_off_bands": "",             # bands 'lo-hi,...' where the side-lobe rule is not applied
        "qc": {
            "snr_c": 6.0, "jk_c_rel": 0.02, "dev_c": 0.20, "smooth_c": 0.03, "ucons_c": 0.25,
            "side_c": 1.0, "dev_pre": 0.30, "lam_c": 2.0, "min_run": 6,
        },
        "group_qc": {
            "snr_g": 8.0, "side_g": 0.9, "jk_g": 0.03, "dev_g": 0.25, "lam_g": 3.0,
            "branch_diff": 0.20, "snr_branch_min": 4.0, "step_g": 0.0,
        },
    },
    "export": {
        "tmin_s": 0.09,
        "tmax_s": "auto",                 # auto: 1 / dispersion.fmin_hz
        "one_segment": True,
        "curve_figures": 200,             # number of pairs drawn as single-pair figures (0 = none, -1 = all)
    },
    "inversion": {
        "wave": "rayleigh_phase",
        "median_min_pairs": 50,
        "depth_sequence_km": [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.9, 2.4, 3.0],
        "fine_step_km": 0.1,              # extra depth nodes every fine_step_km down to fine_zmax_km (0 = none)
        "fine_zmax_km": 1.5,
        "grid_spacing_deg": "auto",       # auto: refined spacing ~ half the median nearest-station distance
        "grid_refine": 2,                 # every grid interval divided by this after the extent is fixed
        "line_strip_cells": None,         # lines: cells kept on each side of the stations across the line
        "start_model": {
            "source": "invert",           # invert: from this array's median curve | file: two columns depth Vs
            "file": None,
            "smooth": 2.0, "damp": 0.05, "iterations": 14, "minthk": 3,
        },
        "weight": 1.0,
        "damp": 1.0,
        "sublayers": 3,
        "vmin": 1.5,
        "vmax": 4.5,
        "iterations": 10,
        "sparsity": 0.05,
        "threshold": 3.0,
        "noise": 0.02,
        "solver": "dsurftomo",            # dsurftomo | pdsurftomo
        "threads": 4,
        "njobs": 1,                       # projects run at the same time
    },
    "checks": {
        # horizontal checkerboards: cell size (half wavelength) north-south and east-west in km
        # auto: cells of 0.12, 0.18 and 0.25 x the array aperture
        "checkerboards": "auto",
        # depth-varying checkerboards: the lateral pattern changes sign at each of these depths
        # auto: cell 0.18 x aperture, sign reversals at 0.21, 0.46, 0.71 x the deepest updated node
        "vertical": "auto",
        "amplitude": 0.4,
        "controls": ["reg_w0.5_d0.5", "reg_w1.5_d1.5", "start_x0.95", "start_x1.05"],
        "min_rays_per_node": 5,
        "iteration_residuals": True,      # evaluate also the model of every iteration (one forward run each)
    },
    "plots": {
        "formats": ["png", "pdf"],
        "dpi": 300,
        "topography": None,               # GMT grid file or a GMT remote grid such as @earth_relief_01s; None = no relief
        "vs_range": [1.9, 3.6],
        "vs_cpt": "seis",
        "slice_depths_km": [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.5],
        "sections": "auto",               # auto: two sections through the array centre along its principal axes; or
                                          # [{name: A, start: [lon, lat], end: [lon, lat]}]
        "section_zmax_km": None,
        "contour_step": 0.2,
        "font_scale": 1.0,
        "pair_example": None,             # pair name drawn in the dispersion example figure (default: chosen automatically)
    },
}

LINE_OVERRIDES = {
    "data": {"min_trace_length_s": 600.0},
    "stack": {"block_hours": 12,
              "station_day_qc": {"enabled": False},
              "station_hour_qc": {"enabled": True}},
    "reference": {"fj_mode": "track_fundamental", "pass2_min_pairs": 20},
    "dispersion": {"njack": "auto", "qc": {"dev_c": 0.0},
                   "group_qc": {"branch_diff": 0.0, "snr_branch_min": 0.0, "step_g": 0.03}},
    "inversion": {"median_min_pairs": 20,
                  "depth_sequence_km": [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.9, 2.4, 3.0],
                  "fine_step_km": 0.0, "grid_refine": 1, "line_strip_cells": 1.0, "grid_spacing_deg": [0.005, 0.005],
                  "weight": 1.5, "damp": 1.5, "sparsity": 0.30},
    "checks": {"controls": ["reg_w0.75_d1.5", "reg_w3_d1.5", "reg_w1.5_d0.75", "reg_w1.5_d3", "start_x0.95", "start_x1.05"]},
}

REQUIRED = [("project", "output_dir"), ("data", "raw_dir")]
PATH_KEYS = [("project", "output_dir"), ("data", "raw_dir"), ("data", "station_file"), ("reference", "file"),
             ("notch", "table_file"), ("inversion", "start_model", "file"), ("plots", "topography")]


class ConfigError(ValueError):
    pass


def _merge(base, over, where=""):
    for k, v in over.items():
        if k not in base:
            raise ConfigError("unknown configuration key '%s%s'" % (where, k))
        if isinstance(base[k], dict) and isinstance(v, dict):
            _merge(base[k], v, where + k + ".")
        else:
            base[k] = v
    return base


def _get(cfg, keys):
    for k in keys:
        cfg = cfg[k]
    return cfg


def _set(cfg, keys, value):
    for k in keys[:-1]:
        cfg = cfg[k]
    cfg[keys[-1]] = value


def defaults(array_type="surface"):
    cfg = copy.deepcopy(DEFAULTS)
    if array_type == "line":
        _merge(cfg, copy.deepcopy(LINE_OVERRIDES))
    cfg["project"]["array_type"] = array_type
    return cfg


def load(path):
    """Read a project YAML file -> full configuration dict (defaults filled in, paths absolute).  Values "auto" are
    resolved from the data by lant.auto after the clean step (array_type included: until then the areal defaults
    are used)."""
    path = os.path.abspath(path)
    with open(path) as fh:
        user = yaml.safe_load(fh) or {}
    if not isinstance(user, dict):
        raise ConfigError("%s: the file must be a YAML mapping" % path)
    array_type = (user.get("project") or {}).get("array_type", "auto")
    if array_type not in ("surface", "line", "auto"):
        raise ConfigError("project.array_type must be 'auto', 'surface' or 'line', not %r" % array_type)
    cfg = build(user, array_type, os.path.dirname(path))
    cfg["_config_file"] = path
    cfg["_user"] = copy.deepcopy(user)
    validate(cfg)
    return cfg


def build(user, array_type, base):
    cfg = _merge(defaults("line" if array_type == "line" else "surface"), copy.deepcopy(user))
    cfg["project"]["array_type"] = array_type
    for keys in REQUIRED:
        if not _get(cfg, keys):
            raise ConfigError("missing required key '%s'" % ".".join(keys))
    for keys in PATH_KEYS:
        v = _get(cfg, keys)
        if isinstance(v, str) and v and not v.startswith("@") and v != "auto":
            _set(cfg, keys, os.path.normpath(os.path.join(base, os.path.expanduser(v))))
    for keys in (("data", "exclude_stations"), ("data", "field_note_stations"), ("stack", "exclude_station_days")):
        v = _get(cfg, keys)
        if isinstance(v, str):
            p = os.path.normpath(os.path.join(base, os.path.expanduser(v)))
            if not os.path.exists(p):
                raise ConfigError("%s: file %s not found" % (".".join(keys), p))
            rows = [l.split("#")[0].split() for l in open(p)]
            rows = [r for r in rows if r]
            _set(cfg, keys, [r if keys[1] == "exclude_station_days" else r[0] for r in rows])
    return cfg


def rebuild(cfg, array_type):
    """The configuration for a given array type: its defaults, then the user's own values (used once the type is known)."""
    new = build(cfg["_user"], array_type, os.path.dirname(cfg["_config_file"]))
    new["_config_file"], new["_user"] = cfg["_config_file"], cfg["_user"]
    return new


def is_auto(v):
    return isinstance(v, str) and v == "auto"


def validate(cfg):
    """Checks that do not need the data (values still "auto" are checked after they are resolved)."""
    c, d = cfg["correlation"], cfg["dispersion"]
    num = lambda *v: not any(is_auto(x) for x in v)
    if num(c["seg_len_s"], c["maxlag_s"]) and c["seg_len_s"] < 2 * c["maxlag_s"]:
        raise ConfigError("correlation.seg_len_s must be at least 2 * maxlag_s (circular correlation)")
    h = cfg["data"]["day_start_utc_hour"]
    if num(h) and not (isinstance(h, int) and 0 <= h < 24):
        raise ConfigError("data.day_start_utc_hour must be 'auto' or an integer 0-23")
    if 24 % cfg["stack"]["block_hours"]:
        raise ConfigError("stack.block_hours must divide 24")
    if num(c["fmax_hz"], c["fs"]) and c["fmax_hz"] >= c["fs"] / 2:
        raise ConfigError("correlation.fmax_hz must be below the Nyquist frequency fs/2")
    if num(d["fmax_hz"], c["fmax_hz"]) and d["fmax_hz"] > c["fmax_hz"]:
        raise ConfigError("dispersion.fmax_hz (%.2f) is above the correlation band (%.2f)" % (d["fmax_hz"], c["fmax_hz"]))
    if d["njack"] != "auto" and (not isinstance(d["njack"], int) or d["njack"] < 2):
        raise ConfigError("dispersion.njack must be an integer >= 2 or 'auto'")
    if cfg["reference"]["source"] not in ("fj", "file"):
        raise ConfigError("reference.source must be 'fj' or 'file'")
    if cfg["reference"]["source"] == "file" and not cfg["reference"]["file"]:
        raise ConfigError("reference.source is 'file' but reference.file is empty")
    if cfg["reference"]["fj_mode"] not in ("max_peak", "track_fundamental"):
        raise ConfigError("reference.fj_mode must be 'max_peak' or 'track_fundamental'")
    sm = cfg["inversion"]["start_model"]
    if sm["source"] not in ("invert", "file") or (sm["source"] == "file" and not sm["file"]):
        raise ConfigError("inversion.start_model.source must be 'invert', or 'file' with start_model.file")
    if cfg["inversion"]["solver"] not in ("dsurftomo", "pdsurftomo"):
        raise ConfigError("inversion.solver must be 'dsurftomo' or 'pdsurftomo'")
    seq = cfg["inversion"]["depth_sequence_km"]
    if sorted(seq) != list(seq) or seq[0] != 0:
        raise ConfigError("inversion.depth_sequence_km must start at 0 and increase")
    return cfg


def dump(cfg, path):
    out = {k: v for k, v in cfg.items() if not k.startswith("_")}
    with open(path, "w") as fh:
        yaml.safe_dump(out, fh, sort_keys=False, allow_unicode=True)


def phase_constant(cfg):
    import numpy as np
    v = cfg["dispersion"]["phase_constant"]
    if isinstance(v, str):
        v = v.replace(" ", "")
        table = {"pi/4": np.pi / 4, "0": 0.0, "-pi/4": -np.pi / 4}
        if v not in table:
            raise ConfigError("dispersion.phase_constant: use a number (rad) or pi/4")
        return table[v]
    return float(v)
