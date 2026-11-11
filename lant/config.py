"""Project configuration: defaults, loading, validation.

A project is described by one YAML file (config/example.yaml; `lant init` writes one).  Keys that are not given
take the defaults below.  Only three paths are required: the raw data, the station table and the output directory.
The array type is read from the station geometry (project.array_type = auto), the day labels are UTC dates unless
data.day_label says otherwise, the velocity range has a wide default, and instrument responses are removed when
response files are given (data.response).  Quantities that follow from the data (bands, lag and window lengths,
velocity axis, depth nodes, grid, ...) are "auto": derived by written rules (lant/auto.py) and shown with their
derivation in the step reports.  The scientific judgements a step leaves to the user are confirmed with
`lant confirm <step>` (lant/confirm.py), which records them in the `confirmed:` block of the configuration file.
A value written in the configuration always wins over a rule.  What every default rests on is in lant/params.py.
A key that does not exist in the defaults is an error, so that a typing mistake cannot silently fall back to a
default.

A line array (project.array_type line, or auto on an elongated array) applies LINE_OVERRIDES on top of the
defaults before the user's file is read.
"""
import copy
import os
import math
import re

import yaml

DEFAULTS = {
    "project": {
        "name": "my_array",
        "array_type": "auto",             # auto: from the station geometry (array_type_of) | surface (areal array) | line
        "output_dir": None,               # required
    },
    "data": {
        "raw_dir": None,                  # required: directory searched recursively; any layout, any file names
        "format": "auto",                 # auto (ObsPy detects MSEED, SAC, GSE2, SEG2, SEISAN, ...) or one ObsPy format name
        "file_patterns": ["*"],           # file-name patterns to try (default: every file)
        "ignore_patterns": ["*.txt", "*.csv", "*.xml", "*.json", "*.yaml", "*.yml", "*.py", "*.sh", "*.png",
                            "*.jpg", "*.pdf", "*.log", "*.md", "*.zip", "*.gz", "*.tar", ".*"],
        "station_file": None,             # required: station table in any common form (lant/stations.py), or StationXML
        "component": "Z",                 # correlations of this component pair: Z = vertical (ZZ, Rayleigh waves); the only
                                          # one supported so far (RR / TT need three components and orientations)
        "channel": "auto",                # auto: the vertical channel of each station (see lant.data.choose_channel); or a
                                          # channel pattern such as "*3" when the names do not tell
        "network_default": None,          # optional override: network code for headers without one
        "rename": {},                     # optional override: {"NET.STA of the headers": "NET.STA of the table"}
        "day_start_utc_hour": "auto",     # auto (the hour most records start) | 0-23: processing days start at this UTC hour
        "day_label": "utc_start",         # label YYYYDDD of a day: utc_start (UTC date of its start) | a fixed
                                          # offset "UTC+8" / "UTC-05:30": the date of its start there (no daylight saving)
        "velocity_range_kms": [0.1, 5.0], # approximate Rayleigh-wave phase velocities of the area [slow, fast] (km/s).
                                          # The default covers soft soil to crystalline rock; a narrower range, e.g.
                                          # [0.2, 1.5] on soft sediments or [1.5, 3.5] on hard rock, shortens the lag
                                          # window and the computing time
        "response": None,                 # instrument responses: a StationXML / dataless SEED / RESP file, or a directory
                                          # of such files; given: the responses are removed in the clean step (stations
                                          # of different models can then be used together); not given: the records are
                                          # used as they are, which is exact only when every station has the same model
        "line": None,                     # station tables with a line column: the line of this project
        "min_trace_length_s": 0.0,        # pieces of record shorter than this (s) are not used
        "exclude_stations": [],           # station names (NET.STA) left out entirely; or a file name
        "start_date": None,               # optional YYYY-MM-DD limits of the processing days
        "end_date": None,
    },
    "clean": {
        "write_archive": True,            # write the cleaned day files (work/01_clean/archive); False = later steps read the raw files
        "nproc": 4,
        "encoding": "FLOAT32",            # miniSEED encoding of the archive (FLOAT32 keeps integer and float data exactly
                                          # up to 2**24 counts; use INT32 or STEIM2 for integer data to save space)
        "min_free_disk_factor": 1.2,      # refuse to write the archive if free disk < factor x the estimated size
    },
    "correlation": {
        "fs_max": 50.0,                   # Hz: highest processing rate; it sets the top of the band, 0.4 fs_max / 1.37
                                          # (14.6 Hz for 50 Hz): raise it for small, high-frequency arrays
        "fs": "auto",                     # Hz after decimation; auto: the highest whole-number rate <= fs_max reachable
                                          # from every raw rate
        "seg_len_s": "auto",              # auto: max(2 maxlag, 10 longest periods)
        "seg_step_s": "auto",             # auto: seg_len / 2
        "maxlag_s": "auto",               # auto: max distance / (0.6 x slowest velocity) + 4 longest periods
        "fmin_hz": "auto",                # auto: longest pair = 2 wavelengths at the slowest velocity, widened to where
        "fmax_hz": "auto",                #       the Gaussian filter of the edge frequency has fallen to half; top 0.4 fs
        "max_distance_km": None,          # only pairs up to this distance are correlated (None: all); for arrays of
                                          # many hundred stations, where all pairs do not fit into memory
        "smooth_bins": 1,                 # spectral smoothing of the coherence denominator (1 = none)
        "max_over_std": 10.0,             # segment rejection: peak / station-day std
        "remove_1s_comb": "auto",         # removal of a UTC-locked 1-s periodic signal (logger timing pulses,
                                          # Takagi et al. 2015): subtract the waveform folded modulo 1 s, hour by hour.
                                          # auto: for the stations where the clean step detects it | true: all | false
        "nproc": 4,
        "threads_per_process": 2,
    },
    "stack": {
        "block_hours": "auto",            # sub-stack length for the jackknife; auto: the longest of 24, 12, 8, 6, 4, 3, 2, 1 h
                                          # giving the median station at least 8 sub-stacks
        "station_qc": {                   # diagnosis only: candidates are listed, nothing is removed automatically
            "enabled": True,
            "partners": 20,               # each station is compared through its pairs with this many nearest stations
            "min_partners": 5,            # a station sub-stack is judged from at least this many pairs
        },
        "exclude_station_days": [],       # [station, YYYYDDD] (one day or sub-stack YYYYDDD_hh), [station, "YYYYDDD-"]
                                          # (from that day on) or [station, "YYYYDDD-YYYYDDD"]; or a file "station day" per line
        "exclude_station_hours": [],      # [station, YYYYDDD, hour] entries, or a file "station day hour" per line
        "exclude_candidates": False,      # true: also remove the station sub-stacks the diagnosis lists
                                          # (set by `lant confirm stack --exclude-candidates`)
    },
    "reference": {
        "source": "fj",                   # fj: F-J of this array | file: borrowed curve (reference.file)
        "file": None,                     # npz written by lant, or a text file "freq_Hz c_km/s"
        "zone_km": None,                  # None: one reference curve for the whole array | size (km) of square zones:
                                          # one F-J curve per zone from the pairs whose midpoint lies in it
        "fj_mode": "max_peak",            # max_peak (areal arrays) | track_fundamental (lines)
        "normal_dispersion": True,        # track_fundamental: c may not rise with frequency (fundamental Rayleigh
                                          # mode without velocity reversal); false where the ridge really rises
        "fj_rmax_lambda": 5.0,
        "fj_c_nominal": "auto",           # auto: the curve of a first pass with every pair (selects the pairs
                                          # within fj_rmax_lambda wavelengths in the second pass)
        "fj_cmin": "auto",                # auto: 0.5 x slowest / 1.5 x fastest of data.velocity_range_kms
        "fj_cmax": "auto",
        "fj_dc": "auto",                  # auto: the 1-2-2.5-5 step not above fj_cmin / 200
        "fj_nfft": "auto",                # auto: power of 2 >= 4 x the lag window, at least 8192
        "fj_fit_fmax": "auto",            # auto: the top of the longest run of clear F-J picks
        "track_start_band_hz": "auto",    # auto: the start (every 1/3 octave) followed over the most frequencies
    },
    "dispersion": {
        "fmin_hz": "auto",                # auto (reference step): the longest pair = qc.lam_c wavelengths of the F-J curve
        "fmax_hz": "auto",                # auto: the highest the processing rate allows (correlation.fs_max); the QC
                                          #   decides how far the data really reach
        "per_octave": 12,
        "vmin": "auto",                   # auto (reference step): 0.6 x the slowest F-J group / phase velocity in the band
        "vmax": "auto",                   # auto: 1.3 x the fastest
        "dv": "auto",                     # auto: the 1-2-2.5-5 step not above vmin / 100
        "alpha": 5.0,
        "phase_constant": "pi/4",
        "min_dist_km": 0.0,               # pairs shorter than this are not measured (lam_min rules the short paths)
        "lam_min": 1.0,                   # paths shorter than this many reference wavelengths are not measured
        "anchor_max_lambda": 3.5,         # cycle count anchored only where the path is 1.0-3.5 wavelengths
        "njack": 5,                       # jackknife groups; "auto" = blocks with >= 25 % of the fullest block
        "min_jackknife_groups": 3,
        "nproc": 8,
        "chunk": 20,
        "side_off_bands": "",             # bands 'lo-hi,...' where the side-lobe rule is not applied
        "qc": {                           # phase velocity; every threshold of the measurement is here (none hidden)
            "snr_c": 6.0, "jk_c_rel": 0.02, "ucons_c": 0.25,
            "side_c": 1.0, "dev_pre": 0.30, "lam_c": 2.0, "min_run": 6,
            "band_edge": 0.5,             # amplitude spectrum at the frequency at least this fraction of its mean over
                                          # the filter band (0.5 = -6 dB; 0 = not applied)
            "snr_g": 4.0,                 # first gate: the largest SNR of the frequency must reach this
        },
    },
    "export": {
        "tmin_s": "auto",                 # auto: 1 / dispersion.fmax_hz
        "tmax_s": "auto",                 # auto: 1 / dispersion.fmin_hz
        "one_segment": True,
        "curve_figures": 30,              # single-pair figures: this many pairs spread over distance, plus the most
                                          # rejected points and the largest deviation from the median (0 = none)
    },
    "inversion": {
        "wave": "rayleigh_phase",
        "median_min_pairs": 50,
        "depth_sequence_km": "auto",      # auto: first step a quarter of the shortest wavelength, then steps of a quarter
                                          # of the depth, down past half the longest wavelength
        "fine_step_km": 0.0,              # extra depth nodes every fine_step_km down to fine_zmax_km (0 = none)
        "fine_zmax_km": 0.0,
        "grid_spacing_deg": "auto",       # auto: refined spacing ~ half the median nearest-station distance
        "grid_refine": 2,                 # every grid interval divided by this after the extent is fixed
        "line_strip_cells": None,         # lines: cells kept on each side of the stations across the line
        "start_model": {
            "source": "invert",           # invert: from this array's median curve | file: two columns depth Vs
            "file": None,
            "smooth": 0.01, "damp": 0.05, "iterations": 14, "minthk": 3,   # smooth: curvature of the model (jumping)
        },
        "weight": 1.0,
        "damp": 1.0,
        "sublayers": 3,
        "vmin": "auto",                   # auto: 0.6 x the slowest Vs of the starting model
        "vmax": "auto",                   # auto: 1.4 x the fastest
        "iterations": 10,
        "sparsity": "auto",               # auto: twice the estimated filling of DSurfTomo's matrix, at least 0.05
        "threshold": 3.0,
        "noise": "auto",                  # checkerboard noise; auto: relative rms of the kept residuals of the main model
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
        "amplitude": 0.10,                # checkerboard perturbation, a fraction of the 1-D Vs of each depth (+-10 %)
        # reg_w<w>_d<d>: other regularisation; start_grad<f>: starting model x(1 + f) at the surface to x(1 - f) at the
        # deepest updatable node, centred on the data (start_x<f>, a uniform factor, fails with DSurfTomo's rejection
        # rule: the precheck marks it invalid)
        "controls": ["start_grad0.05"],          # add e.g. "reg_w0.5_d0.5", "reg_w2_d2" for other regularisation
        "precheck": {"min_kept": 0.5, "max_median_iqr": 1.0},   # a project is invalid below / beyond these
        "iteration_residuals": False,     # evaluate also the model of every iteration (one forward run each)
    },
    "plots": {
        "formats": ["png"],
        "dpi": 300,
        "topography": "auto",             # auto: SRTM 1" tiles from the GMT data server; a grid file; null = no relief
        "vs_range": "auto",               # auto: 2nd-98th percentile of the model under the stations, to 0.1 km/s
        "vs_cpt": "seis",
        "slice_depths_km": "auto",        # auto: up to 8 updated depth nodes, evenly chosen
        "sections": "auto",               # auto: two sections through the array centre along its principal axes; or
                                          # [{name: A, start: [lon, lat], end: [lon, lat]}]
        "section_zmax_km": None,
        "contour_step": "auto",           # auto: about a tenth of the colour range, 1-2-2.5-5 step
        "font_scale": 1.0,
        "pair_example": None,             # pair(s) of the single-pair figures (default: three, short / middle / long)
    },
}

LINE_OVERRIDES = {
    "data": {"min_trace_length_s": 600.0},
    "reference": {"fj_mode": "track_fundamental"},
    "dispersion": {"njack": "auto"},
    "inversion": {"median_min_pairs": 20, "grid_refine": 1, "line_strip_cells": 1.0,
                  "weight": 1.5, "damp": 1.5},
}

REQUIRED = [("project", "output_dir"), ("data", "raw_dir"), ("data", "station_file")]
REQUIRED_HELP = {
    "project.output_dir": "the directory the results are written to",
    "data.raw_dir": "the directory with the raw records (searched recursively)",
    "data.station_file": "the station table (names and coordinates)",
}
# project.array_type = auto: elongation of the station cloud = sqrt(ratio of the eigenvalues of the covariance of
# the station positions).  Empirical limits: a square or round areal array is near 1, a straight line with a few
# percent of scatter is above 10; in between the choice is shown for confirmation at the clean step.
ELONGATION_SURFACE, ELONGATION_LINE, ELONGATION_SPLIT = 3.0, 10.0, 5.5
CONFIRM_STEPS = ("clean", "stack", "reference", "inversion")
PATH_KEYS = [("project", "output_dir"), ("data", "raw_dir"), ("data", "station_file"), ("data", "response"),
             ("reference", "file"),
             ("inversion", "start_model", "file"), ("plots", "topography")]


class ConfigError(ValueError):
    pass


FREE_MAPPINGS = {"data.rename"}          # keys whose value is a mapping with arbitrary keys (not a section)
REMOVED = {                              # keys of earlier versions: what to do instead
    "notch": "the notch step was removed in 1.2 (no literature basis); delete the 'notch' section - frequencies "
             "disturbed by spectral lines are rejected by the dispersion QC",
    "reference.pass2": "the second pass was removed in 1.2 (it is self-consistent, not an independent check); "
                       "delete the key",
    "reference.pass2_min_pairs": "removed in 1.2 with reference.pass2; delete the key",
    "reference.pass2_max_gap_octave": "removed in 1.2 with reference.pass2; delete the key",
    "stack.station_day_qc": "replaced in 1.3 by stack.station_qc (one diagnosis for areal arrays and lines, "
                            "thresholds from the distribution of the data); delete the section",
    "stack.station_hour_qc": "replaced in 1.3 by stack.station_qc (one diagnosis for areal arrays and lines, "
                             "thresholds from the distribution of the data); delete the section",
    "reference.fj_threads": "removed in 1.3 (the F-J spectrum is computed inside lant; CC-FJpy is no longer used); "
                            "delete the key",
    "data.instrument": "removed in 1.3: give the response files in data.response (they are then removed), or "
                       "nothing when every station has the same instrument model; delete the key",
    "dispersion.qc.width_g": "removed in 1.2 (it was never applied); delete the key",
    "dispersion.qc.dev_c": "removed in 1.4 (it changed nothing in the tests with noise and a second arrival, "
                           "docs/validation.md); delete the key",
    "dispersion.qc.smooth_c": "removed in 1.4 (it changed nothing in the tests with noise and a second arrival, "
                              "docs/validation.md); delete the key",
}


def _merge(base, over, where=""):
    if not isinstance(over, dict):
        raise ConfigError("'%s' must be a mapping (key: value), not %r" % (where.rstrip(".") or "the file", over))
    for k, v in over.items():
        if k not in base:
            if (where + k) in REMOVED:
                raise ConfigError("configuration key '%s%s': %s" % (where, k, REMOVED[where + k]))
            import difflib
            close = difflib.get_close_matches(str(k), [str(b) for b in base], n=1, cutoff=0.7)
            raise ConfigError("unknown configuration key '%s%s'%s" % (where, k, " (did you mean '%s%s'?)" % (
                where, close[0]) if close else "; all keys: lant init --full"))
        if (where + k) in FREE_MAPPINGS:
            if not isinstance(v, dict) or not all(isinstance(a, str) and isinstance(b, str) for a, b in v.items()):
                raise ConfigError("%s%s must map names to names, e.g. {\"XX.001\": \"YY.A01\"}" % (where, k))
            if len(set(v.values())) != len(v):
                raise ConfigError("%s%s: two names are renamed to the same name" % (where, k))
            base[k] = dict(v)
        elif isinstance(base[k], dict):
            _merge(base[k], v, where + k + ".")
        else:
            base[k] = v
    return base


class _UniqueKeyLoader(yaml.SafeLoader):
    """YAML loader that refuses a key written twice in the same mapping (the second would silently win)."""


def _no_duplicates(loader, node, deep=False):
    seen = set()
    for kn, _ in node.value:
        k = loader.construct_object(kn, deep=deep)
        if k in seen:
            raise ConfigError("key %r is written twice (line %d)" % (k, kn.start_mark.line + 1))
        seen.add(k)
    return loader.construct_mapping(node, deep)


_UniqueKeyLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _no_duplicates)


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
    """Read a project YAML file -> full configuration dict (defaults filled in, paths absolute).  The `confirmed:`
    block (written by `lant confirm`) is kept apart in cfg["_confirmed"]: {step: {dotted key: value, "_at": time,
    "_basis": fingerprint of what was looked at}}; lant.auto.apply uses its values where the configuration says
    "auto" (a value written in the configuration wins)."""
    path = os.path.abspath(path)
    with open(path) as fh:
        user = yaml.load(fh, Loader=_UniqueKeyLoader) or {}
    if not isinstance(user, dict):
        raise ConfigError("%s: the file must be a YAML mapping" % path)
    confirmed = user.pop("confirmed", None) or {}
    if not isinstance(confirmed, dict) or any(k not in CONFIRM_STEPS or not isinstance(v, dict)
                                              for k, v in confirmed.items()):
        raise ConfigError("the confirmed: block must map the steps %s to their accepted values (written by `lant "
                          "confirm`)" % ", ".join(CONFIRM_STEPS))
    line = (user.get("data") or {}).get("line")
    array_type = (user.get("project") or {}).get("array_type", "auto")
    if array_type not in ("surface", "line", "auto"):
        raise ConfigError("project.array_type must be auto, surface (areal array) or line, not %r" % (array_type,))
    how = dict(given=True)
    if array_type == "auto":
        sf = (user.get("data") or {}).get("station_file")
        array_type, how = array_type_of(os.path.normpath(os.path.join(os.path.dirname(path), os.path.expanduser(sf)))
                                        if isinstance(sf, str) and sf else None, line)
    cfg = build(user, array_type, os.path.dirname(path))
    cfg["_array_type"] = how
    cfg["_config_file"] = path
    cfg["_user"] = copy.deepcopy(user)
    cfg["_confirmed"] = copy.deepcopy(confirmed)
    for step, vals in confirmed.items():
        for key in vals:
            if not key.startswith("_"):
                try:
                    _get(cfg, key.split("."))
                except (KeyError, TypeError):
                    raise ConfigError("confirmed.%s: unknown key %s" % (step, key))
    validate(cfg)
    marker = os.path.join(cfg["project"]["output_dir"], ".lant", "layout.json")
    if os.path.exists(marker):
        import json
        from .delivery import LAYOUT_VERSION
        try:
            with open(marker) as fh:
                version = json.load(fh).get("version")
        except (ValueError, AttributeError) as exc:
            raise ConfigError("invalid layout record %s: %s" % (marker, exc))
        if version != LAYOUT_VERSION:
            raise ConfigError("unsupported project layout version: %r" % version)
    return cfg


def array_type_of(station_file, line=None):
    """(surface | line, how) from the station table: the elongation of the station cloud against the limits
    ELONGATION_*.  how: dict(given=False, elongation, clear, text).  An unreadable table gives surface (the clean
    step reports the table's problem)."""
    try:
        from . import stations as ST
        from .signal import local_xy_km
        import numpy as np
        rows, _ = ST.read_table(station_file)
        if line:
            rows = [r for r in rows if str(r.get("line")) == str(line)]
        x, y = local_xy_km([r["lon"] for r in rows], [r["lat"] for r in rows])
        if len(rows) < 3:
            return "line", dict(given=False, elongation=None, clear=True, text="fewer than 3 stations")
        w = np.sort(np.linalg.eigvalsh(np.cov(np.vstack([x, y]))))
        e = float(np.sqrt(w[1] / max(w[0], 1e-30)))
    except Exception as exc:
        return "surface", dict(given=False, elongation=None, clear=True, text="station table not read (%s)" % exc)
    kind = "line" if e >= ELONGATION_SPLIT else "surface"
    clear = e <= ELONGATION_SURFACE or e >= ELONGATION_LINE
    return kind, dict(given=False, elongation=round(min(e, 1e6), 2), clear=bool(clear),
                      text="the station cloud is %.1f times longer than wide (%d stations): %s array%s"
                           % (min(e, 1e6), len(rows), "a line" if kind == "line" else "an areal",
                              "" if clear else " - between the clear cases (below %g areal, above %g line)"
                              % (ELONGATION_SURFACE, ELONGATION_LINE)))


def build(user, array_type, base):
    cfg = _merge(defaults("line" if array_type == "line" else "surface"), copy.deepcopy(user))
    cfg["project"]["array_type"] = array_type
    for keys in REQUIRED:
        if _get(cfg, keys) in (None, "", []):
            k = ".".join(keys)
            raise ConfigError("missing required key '%s'%s" % (k, ": " + REQUIRED_HELP[k] if k in REQUIRED_HELP else ""))
    for keys in PATH_KEYS:
        v = _get(cfg, keys)
        if isinstance(v, str) and v and not v.startswith("@") and v != "auto":
            _set(cfg, keys, os.path.normpath(os.path.join(base, os.path.expanduser(v))))
    for keys in (("data", "exclude_stations"), ("stack", "exclude_station_days"), ("stack", "exclude_station_hours")):
        v = _get(cfg, keys)
        if isinstance(v, str):
            p = os.path.normpath(os.path.join(base, os.path.expanduser(v)))
            if not os.path.exists(p):
                raise ConfigError("%s: file %s not found" % (".".join(keys), p))
            rows = [l.split("#")[0].split() for l in open(p)]
            rows = [r for r in rows if r]
            _set(cfg, keys, rows if keys[0] == "stack" else [r[0] for r in rows])
    return cfg


def np_isfinite(v):
    import math
    return math.isfinite(v)


def is_auto(v):
    return isinstance(v, str) and v == "auto"


POSITIVE = [("correlation", "fs"), ("correlation", "seg_len_s"), ("correlation", "seg_step_s"), ("correlation", "maxlag_s"),
            ("correlation", "fmin_hz"), ("correlation", "fmax_hz"), ("correlation", "max_over_std"),
            ("dispersion", "fmin_hz"), ("dispersion", "fmax_hz"), ("dispersion", "vmin"), ("dispersion", "vmax"),
            ("dispersion", "dv"), ("dispersion", "alpha"), ("dispersion", "per_octave"), ("export", "tmin_s"),
            ("export", "tmax_s"), ("inversion", "vmin"), ("inversion", "vmax"), ("inversion", "iterations"),
            ("inversion", "sublayers"), ("inversion", "grid_refine"), ("checks", "amplitude"), ("correlation", "fs_max"),
            ("stack", "block_hours"), ("inversion", "median_min_pairs")]
INTEGERS = [("dispersion", "per_octave"), ("dispersion", "min_jackknife_groups"), ("dispersion", "qc", "min_run"),
            ("dispersion", "nproc"), ("dispersion", "chunk"), ("correlation", "nproc"),
            ("correlation", "threads_per_process"), ("correlation", "smooth_bins"),
            ("inversion", "iterations"), ("inversion", "sublayers"), ("inversion", "grid_refine"),
            ("inversion", "threads"), ("inversion", "njobs"), ("inversion", "median_min_pairs"),
            ("inversion", "start_model", "iterations"), ("inversion", "start_model", "minthk")]
ORDERED = [(("correlation", "fmin_hz"), ("correlation", "fmax_hz")), (("dispersion", "fmin_hz"), ("dispersion", "fmax_hz")),
           (("dispersion", "vmin"), ("dispersion", "vmax")), (("inversion", "vmin"), ("inversion", "vmax")),
           (("export", "tmin_s"), ("export", "tmax_s"))]


def validate(cfg):
    """Checks that do not need the data (values still "auto" are checked after they are resolved)."""
    c, d = cfg["correlation"], cfg["dispersion"]
    num = lambda *v: not any(is_auto(x) for x in v)
    for keys in POSITIVE:
        v = _get(cfg, keys)
        if is_auto(v):
            continue
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not np_isfinite(v) or v <= 0:
            raise ConfigError("%s must be a positive number, not %r" % (".".join(keys), v))
    for a, b in ORDERED:
        va, vb = _get(cfg, a), _get(cfg, b)
        if num(va, vb) and not va < vb:
            raise ConfigError("%s (%r) must be below %s (%r)" % (".".join(a), va, ".".join(b), vb))
    for keys in (("reference", "track_start_band_hz"),):
        v = _get(cfg, keys)
        if not is_auto(v) and not (isinstance(v, (list, tuple)) and len(v) == 2 and 0 <= v[0] < v[1]):
            raise ConfigError("%s must be [low, high] with 0 <= low < high, not %r" % (".".join(keys), v))
    vr = cfg["data"]["velocity_range_kms"]
    if not (isinstance(vr, (list, tuple)) and len(vr) == 2 and all(isinstance(x, (int, float)) and not isinstance(x, bool)
                                                                   for x in vr) and 0 < vr[0] < vr[1] < 20):
        raise ConfigError("data.velocity_range_kms must be [slow, fast] in km/s with 0 < slow < fast, not %r" % (vr,))
    nz = cfg["inversion"]["noise"]
    if not is_auto(nz) and (isinstance(nz, bool) or not isinstance(nz, (int, float)) or not np_isfinite(nz) or nz < 0):
        raise ConfigError("inversion.noise must be auto or a number >= 0 (a fraction of each travel time), not %r" % (nz,))
    if cfg["correlation"]["remove_1s_comb"] not in ("auto", True, False):
        raise ConfigError("correlation.remove_1s_comb must be auto, true or false")
    if cfg["data"]["response"] and not cfg["clean"]["write_archive"]:
        raise ConfigError("data.response needs clean.write_archive: true (the responses are removed when the cleaned "
                          "archive is written)")
    if not isinstance(d["per_octave"], int):
        raise ConfigError("dispersion.per_octave must be a whole number")
    if cfg["data"]["component"] != "Z":
        raise ConfigError("data.component %r: only Z (vertical-component correlations, Rayleigh waves) is supported "
                          "so far" % cfg["data"]["component"])
    if num(c["seg_len_s"], c["maxlag_s"]) and c["seg_len_s"] < 2 * c["maxlag_s"]:
        raise ConfigError("correlation.seg_len_s must be at least 2 * maxlag_s (circular correlation)")
    try:
        from .data import label_offset_s
        label_offset_s(cfg["data"].get("day_label"))
    except ValueError as exc:
        raise ConfigError(str(exc))
    for k in ("start_date", "end_date"):
        v = cfg["data"][k]
        if v is not None:
            v = cfg["data"][k] = str(v)        # YAML reads an unquoted 2024-05-01 as a date object
        if v is not None and not re.fullmatch(r"\d{4}-\d{2}-\d{2}", v):
            raise ConfigError("data.%s must be YYYY-MM-DD (the date of a processing-day label)" % k)
    h = cfg["data"]["day_start_utc_hour"]
    if num(h) and not (isinstance(h, int) and 0 <= h < 24):
        raise ConfigError("data.day_start_utc_hour must be 'auto' or an integer 0-23")
    bh = cfg["stack"]["block_hours"]
    if num(bh) and (isinstance(bh, bool) or not isinstance(bh, int) or 24 % bh):
        raise ConfigError("stack.block_hours must be a whole number of hours that divides 24 (1, 2, 3, 4, 6, 8, 12, 24), "
                          "not %r" % (bh,))
    for keys in INTEGERS:
        v = _get(cfg, keys)
        if num(v) and (isinstance(v, bool) or not isinstance(v, int) or v < 1):
            raise ConfigError("%s must be a whole number >= 1, not %r" % (".".join(keys), v))
    for k, v in d["qc"].items():
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not np_isfinite(v) or v < 0:
            raise ConfigError("dispersion.qc.%s must be a number >= 0, not %r" % (k, v))
    if num(d["alpha"]) and d["alpha"] <= math.log(2):
        raise ConfigError("dispersion.alpha must be above ln 2 = 0.693 (the Gaussian filter of the band edge must "
                          "reach half gain at a positive frequency), not %r" % d["alpha"])
    if not (isinstance(d["lam_min"], (int, float)) and d["lam_min"] > 0):
        raise ConfigError("dispersion.lam_min must be a positive number of wavelengths")
    if d["anchor_max_lambda"] < d["lam_min"] + 1.5:
        raise ConfigError("dispersion.anchor_max_lambda (%r) must be at least dispersion.lam_min + 1.5 (%r): the cycle "
                          "count is anchored where the path is lam_min to anchor_max_lambda wavelengths, a window of at "
                          "least 1.5 wavelengths" % (d["anchor_max_lambda"], d["lam_min"] + 1.5))
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
    if cfg["inversion"]["wave"] != "rayleigh_phase":
        raise ConfigError("inversion.wave %r: only rayleigh_phase (Rayleigh-wave phase velocities) is supported"
                          % cfg["inversion"]["wave"])
    if cfg["inversion"]["solver"] not in ("dsurftomo", "pdsurftomo"):
        raise ConfigError("inversion.solver must be 'dsurftomo' or 'pdsurftomo'")
    seq = cfg["inversion"]["depth_sequence_km"]
    if num(seq) and (seq[0] != 0 or any(round(b, 3) <= round(a, 3) for a, b in zip(seq[:-1], seq[1:]))):
        raise ConfigError("inversion.depth_sequence_km must start at 0 and increase by at least 0.001 km (DSurfTomo's "
                          "model file keeps 3 decimals of km)")
    return cfg


def dump(cfg, path, base=None):
    """base: the directory the file will be read from, when it is written elsewhere first (default: its own)."""
    out = copy.deepcopy({k: v for k, v in cfg.items() if not k.startswith("_")})
    root = cfg["project"]["output_dir"]
    base = base or os.path.dirname(os.path.abspath(path))
    for keys in PATH_KEYS:
        value = _get(out, keys)
        if isinstance(value, str) and os.path.isabs(value) and os.path.commonpath([root, value]) == root:
            _set(out, keys, os.path.relpath(value, base))
    with open(path, "w") as fh:
        yaml.safe_dump(out, fh, sort_keys=False, allow_unicode=True)
    if cfg.get("_confirmed"):
        from .confirm import _write_block
        _write_block(path, cfg["_confirmed"])


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
