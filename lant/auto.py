"""Parameters given as "auto": proposed by written rules, confirmed by the user.

The values and the numbers they came from are written to <output>/.lant/auto_parameters.json.  At the four
confirmation points (after clean, stack, reference and inversion) `lant confirm <step>` writes the accepted values
into the `confirmed:` block of the configuration file; a later step does not run before that.  Order of precedence:
a value written in the configuration, then a confirmed value, then the rule.  docs/parameters.md gives the reasons.

Clean step, from the station geometry and the records (infer):
  data.day_start_utc_hour        the UTC hour at which most records longer than 12 h start (0 when they are spread)
  stack.block_hours              the longest of 24, 12, 8, 6, 4, 3, 2, 1 h giving the median station >= 8 sub-stacks
  inversion.grid_spacing_deg     half the median nearest-station distance (at least the longest extent / 150), on
                                 the 1-2-2.5-5 sequence of degrees, times grid_refine for the coarse grid
  checks.checkerboards           cells of 0.12, 0.18, 0.25 x the aperture (lines: 0.165, 0.24, 0.34 x the line length,
                                 along the line), to 0.1 km
  checks.vertical                set by the inversion step from the depth nodes
  plots.sections                 areal: two sections through the centre of the stations along the principal axes
Clean step, from the velocity range the user gives (data.velocity_range_kms) and the sampling rates (from_velocity):
  correlation.fs / fmin_hz / fmax_hz / maxlag_s / seg_len_s / seg_step_s, dispersion.fmax_hz, export.tmin_s,
  reference.fj_cmin / fj_cmax / fj_dc / fj_c_nominal
Reference step, from the F-J curve of this array or the borrowed curve (from_reference):
  dispersion.fmin_hz / vmin / vmax / dv, export.tmax_s; reference.fj_fit_fmax and track_start_band_hz from the
  F-J image itself (lant.fj)
"""

import json
import math
import os
from collections import Counter

import numpy as np

from . import config as C
from .signal import common_rates, geodesic_km, km_per_degree, local_xy_km, nearest_neighbour_km, nice_number

GEOMETRY_KEYS = [
    ("data", "day_start_utc_hour"),
    ("stack", "block_hours"),
    ("inversion", "grid_spacing_deg"),
    ("checks", "checkerboards"),
    ("checks", "vertical"),
    ("plots", "sections"),
]
VELOCITY_KEYS = [
    ("correlation", "fs"),
    ("correlation", "fmin_hz"),
    ("correlation", "fmax_hz"),
    ("correlation", "maxlag_s"),
    ("correlation", "seg_len_s"),
    ("correlation", "seg_step_s"),
    ("dispersion", "fmax_hz"),
    ("export", "tmin_s"),
    ("reference", "fj_c_nominal"),
    ("reference", "fj_cmin"),
    ("reference", "fj_cmax"),
    ("reference", "fj_dc"),
]
CLEAN_KEYS = GEOMETRY_KEYS + VELOCITY_KEYS
REFERENCE_KEYS = [
    ("dispersion", "fmin_hz"),
    ("dispersion", "vmin"),
    ("dispersion", "vmax"),
    ("dispersion", "dv"),
    ("export", "tmax_s"),
]
KEYS = CLEAN_KEYS + REFERENCE_KEYS
LATER = [("checks", "vertical")]  # resolved by the inversion step
FILTER_EDGE = math.log(2.0)  # the Gaussian filter of the band-edge frequency has fallen to half here


def path(cfg):
    return os.path.join(cfg["project"]["output_dir"], ".lant", "auto_parameters.json")


def geometry(stations):
    names = sorted(stations)
    lon = np.array([stations[n]["lon"] for n in names])
    lat = np.array([stations[n]["lat"] for n in names])
    iu, ju = np.triu_indices(len(names), 1)
    d, _ = geodesic_km(lon[iu], lat[iu], lon[ju], lat[ju])
    nn = nearest_neighbour_km(lon, lat)
    x, y = local_xy_km(lon, lat)
    w, v = np.linalg.eigh(np.cov(np.c_[x, y].T))
    u_major = v[:, np.argmax(w)]
    s = x * u_major[0] + y * u_major[1]
    return dict(
        n=len(names),
        d=d,
        nn=nn,
        elongation=float(np.sqrt(w.max() / max(w.min(), 1e-12))),
        aperture=float(s.max() - s.min()),
        axes=v,
        lon=lon,
        lat=lat,
        x=x,
        y=y,
        extent=float(max(x.max() - x.min(), y.max() - y.min())),
    )


def geometry_summary(g):
    return dict(
        stations=g["n"],
        elongation=g["elongation"],
        aperture_km=g["aperture"],
        p90_distance_km=float(np.percentile(g["d"], 90)),
        max_distance_km=float(g["d"].max()),
        median_distance_km=float(np.median(g["d"])),
        median_nearest_km=float(np.nanmedian(g["nn"])),
    )


def user_value(cfg, keys):
    """The value the user wrote in the configuration, or None ("auto" is the same as leaving the key out)."""
    u = cfg.get("_user") or {}
    for k in keys:
        if not isinstance(u, dict) or k not in u:
            return None
        u = u[k]
    return None if C.is_auto(u) else u


def noise_is_auto(cfg):
    """inversion.noise left at auto: derived by the run step from the fit of the main model (lant.pipeline)."""
    return user_value(cfg, ["inversion", "noise"]) is None


# derived by a later step than the one that reads them: the earlier step records the setting ("auto" or the user's
# value), not the derived value, so that deriving it does not make the earlier step out of date
DERIVED_LATER = {
    "inversion": {
        "inversion.noise": lambda cfg: "auto" if noise_is_auto(cfg) else user_value(cfg, ["inversion", "noise"])
    }
}


class _Put:
    """Collects {dotted key: {value, reason}}; a value given in the configuration wins."""

    def __init__(self, cfg, out=None):
        self.cfg, self.out = cfg, out if out is not None else {}

    def __call__(self, keys, value, reason):
        uv = user_value(self.cfg, keys)
        if uv is not None and not C.is_auto(uv):
            self.out[".".join(keys)] = dict(value=uv, reason="given in the configuration")
        else:
            self.out[".".join(keys)] = dict(value=value, reason=reason)
        return self.out[".".join(keys)]["value"]


def record_hours(traces):
    """Median over the stations of the hours of record (union of the record intervals: overlapping or repeated
    files are counted once)."""
    iv = {}
    for t in traces:
        iv.setdefault(t.get("station") or (t["net"], t["sta"]), []).append((t["start"], t["end"]))
    tot = []
    for lst in iv.values():
        lst.sort()
        h, cur0, cur1 = 0.0, None, None
        for a, b in lst:
            if cur1 is None or a > cur1:
                if cur1 is not None:
                    h += cur1 - cur0
                cur0, cur1 = a, b
            else:
                cur1 = max(cur1, b)
        if cur1 is not None:
            h += cur1 - cur0
        tot.append(h / 3600)
    return float(np.median(tot)) if tot else 0.0


def infer(cfg, stations, traces):
    """Clean step: the parameters that follow from the station geometry and the records."""
    g = geometry(stations)
    put = _Put(cfg)
    at = cfg["project"]["array_type"]
    long_rec = [t for t in traces if t["end"] - t["start"] >= 12 * 3600]
    import obspy

    hours = Counter(obspy.UTCDateTime(t["start"] + 30).hour for t in long_rec)
    if hours and hours.most_common(1)[0][1] >= 0.6 * len(long_rec):
        h, n = hours.most_common(1)[0]
        put(
            ("data", "day_start_utc_hour"),
            int(h),
            "%d of %d records longer than 12 h start at %02d h UTC" % (n, len(long_rec), h),
        )
    else:
        put(("data", "day_start_utc_hour"), 0, "no common start hour of the long records (%s): UTC days" % dict(hours))
    rh = record_hours(traces)
    bh = next((b for b in (24, 12, 8, 6, 4, 3, 2) if rh / b >= 8), 1)
    put(
        ("stack", "block_hours"),
        bh,
        "median record %.1f h per station: %d h gives %.0f sub-stacks (>= 8 wanted)" % (rh, bh, rh / bh),
    )
    nn_med = float(np.nanmedian(g["nn"]))
    want_km = max(nn_med / 2, g["extent"] / 150)
    fine = nice_number(want_km / 111.195)
    r = cfg["inversion"]["grid_refine"]
    put(
        ("inversion", "grid_spacing_deg"),
        [round(fine * r, 6)] * 2,
        "max(half the median nearest-station distance %.3f km, extent %.1f km / 150) = %.3f km -> %g deg (x grid_refine %d)"
        % (nn_med, g["extent"], want_km, fine, r),
    )
    ap = g["aperture"]
    frac = (0.165, 0.24, 0.34) if at == "line" else (0.12, 0.18, 0.25)
    cells = [float("%.2g" % (f * ap)) for f in frac]  # two significant digits: any aperture
    put(
        ("checks", "checkerboards"),
        [{"name": "cb_%d_%gkm" % (i + 1, c), "cell_km": [c, 0.0 if at == "line" else c]} for i, c in enumerate(cells)],
        "%s x the aperture %.1f km" % (" / ".join("%g" % f for f in frac), ap),
    )
    put(("checks", "vertical"), "auto", "set by the inversion step from the depth nodes")
    if at == "surface":
        lat0, lon0 = float(np.mean(g["lat"])), float(np.mean(g["lon"]))
        kx, ky = km_per_degree(lat0)
        secs = []
        for name, u in zip("AB", (g["axes"][:, 1], g["axes"][:, 0])):
            s = g["x"] * u[0] + g["y"] * u[1]
            ends = []
            for sv in (s.min(), s.max()):
                ends.append([round(lon0 + sv * u[0] / kx, 5), round(lat0 + sv * u[1] / ky, 5)])
            secs.append({"name": name, "start": ends[0], "end": ends[1]})
        put(("plots", "sections"), secs, "through the centre of the stations along the principal axes (A: long axis)")
    else:
        put(("plots", "sections"), [], "line arrays: section along the line")
    gs = geometry_summary(g)
    gs["record_hours_median"] = rh
    return put.out, gs


def _grid_down(f, per_octave):
    """The 1/per_octave-octave grid through 1 Hz, rounded down / up."""
    return 2 ** (math.floor(per_octave * math.log2(f) + 1e-9) / per_octave)


def _grid_up(f, per_octave):
    return 2 ** (math.ceil(per_octave * math.log2(f) - 1e-9) / per_octave)


def from_velocity(cfg, geom, raw_rates):
    """Clean step, part 2: the parameters needed before the correlations, from the approximate velocity range the
    user gives (data.velocity_range_kms), the geometry and the sampling rates.  Each is a theory relation with a
    margin (the margins are project values, docs/parameters.md):
      lowest frequency  f_lo = 2 v_slow / max distance (the longest pair is 2 wavelengths at the slowest velocity);
      correlation.fs    the highest whole-number rate <= correlation.fs_max reachable from every raw rate;
      band top          dispersion.fmax_hz = 0.4 fs / (1 + sqrt(ln 2 / alpha)) on the 1/per_octave grid;
      correlation band  widened to where the Gaussian filter of the edge frequency falls to half;
      maxlag            max distance / (0.6 v_slow) (group velocities are slower than phase velocities) + 4 periods
                        of f_lo; segments max(2 maxlag, 10 periods of f_lo), half overlap;
      F-J scan          0.5 v_slow to 1.5 v_fast, step <= 1/200 of the slowest (at most 3000 steps), nominal
                        velocity the geometric mean."""
    put = _Put(cfg)
    v_lo, v_hi = (float(x) for x in cfg["data"]["velocity_range_kms"])
    v_mid = math.sqrt(v_lo * v_hi)
    po = int(cfg["dispersion"]["per_octave"])
    dlt = math.sqrt(FILTER_EDGE / float(cfg["dispersion"]["alpha"]))
    dmax = float(geom["max_distance_km"])
    cap = float(cfg["correlation"]["fs_max"])
    rates = [r for r in common_rates(raw_rates) if 1 <= r <= cap]
    if not rates:
        raise RuntimeError(
            "no whole-number processing rate <= correlation.fs_max %g Hz is reachable from %s Hz" % (cap, raw_rates)
        )
    fs = put(
        ("correlation", "fs"),
        float(max(rates)),
        "the highest whole-number rate <= correlation.fs_max %g Hz reachable from %s Hz" % (cap, raw_rates),
    )
    fmax = put(
        ("dispersion", "fmax_hz"),
        float(_grid_down(0.4 * fs / (1 + dlt) * (1 + 1e-9), po)),
        "0.4 x %g Hz / (1 + sqrt(ln 2 / %g)) = %.4g Hz, taken down to the 1/%d-octave grid: the top the "
        "processing rate allows; the QC decides how far the data really reach"
        % (fs, cfg["dispersion"]["alpha"], 0.4 * fs / (1 + dlt), po),
    )
    f_lo = _grid_up(2 * v_lo / dmax * (1 - 1e-9), po)
    if f_lo >= fmax:
        raise RuntimeError(
            "the longest pair (%.3f km) is shorter than 2 wavelengths at %.3g km/s up to %.3g Hz: raise "
            "correlation.fs_max or check data.velocity_range_kms" % (dmax, v_lo, fmax)
        )
    put(
        ("correlation", "fmin_hz"),
        round(f_lo * (1 - dlt), 6),
        "%.4g Hz (the longest pair %.3f km = 2 wavelengths at the slowest velocity %.3g km/s) x (1 - sqrt(ln 2 "
        "/ alpha)): the Gaussian filter of that frequency at half gain" % (f_lo, dmax, v_lo),
    )
    put(
        ("correlation", "fmax_hz"),
        round(min(fmax * (1 + dlt), 0.4 * fs), 6),
        "dispersion.fmax_hz x (1 + sqrt(ln 2 / alpha)), not above 0.4 fs",
    )
    put(("export", "tmin_s"), 1 / fmax, "1 / dispersion.fmax_hz = 1 / %.4g Hz" % fmax)
    fine = (1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8)
    ml = put(
        ("correlation", "maxlag_s"),
        nice_number(dmax / (0.6 * v_lo) + 4 / f_lo, "up", fine),
        "max distance %.3f km / (0.6 x %.3g km/s) + 4 periods of %.4g Hz (2 to clear the arrival, 2 of noise), "
        "rounded up" % (dmax, v_lo, f_lo),
    )
    sl = put(
        ("correlation", "seg_len_s"),
        nice_number(max(2 * ml, 10 / f_lo), "up", fine),
        "max(2 x maxlag (circular correlation), 10 periods of %.4g Hz), rounded up" % f_lo,
    )
    put(("correlation", "seg_step_s"), sl / 2, "half the segment of %g s (50 %% overlap)" % sl)
    steps = (1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8)
    fj0 = put(
        ("reference", "fj_cmin"),
        round(nice_number(0.5 * v_lo, "down", steps), 4),
        "0.5 x the slowest velocity %.3g km/s, rounded down" % v_lo,
    )
    put(
        ("reference", "fj_cmax"),
        round(nice_number(1.5 * v_hi, "up", steps), 4),
        "1.5 x the fastest velocity %.3g km/s, rounded up" % v_hi,
    )
    fj1 = put.out["reference.fj_cmax"]["value"]
    put(
        ("reference", "fj_dc"),
        nice_number(max(fj0 / 200, (fj1 - fj0) / 3000), "down"),
        "the 1-2-2.5-5 step not above max(fj_cmin / 200 = %.3g, scan width / 3000 = %.3g) km/s"
        % (fj0 / 200, (fj1 - fj0) / 3000),
    )
    put(
        ("reference", "fj_c_nominal"),
        round(v_mid, 4),
        "not used unless given: the reference step limits the pairs with the curve of a first pass that takes every pair "
        "(shown here: sqrt(%.3g x %.3g) km/s)" % (v_lo, v_hi),
    )
    return put.out


def from_reference(cfg, freqs, c, u, picked_hz, geom):
    """Reference step: the measurement band and velocity axis from the F-J curve c(f), U(f) of this array (or the
    borrowed curve): dispersion.fmin_hz where the longest pair is qc.lam_c wavelengths; vmin / vmax 0.6 x the
    slowest and 1.3 x the fastest of c and U in the band; dv the 1-2-2.5-5 step not above vmin / 100.  The fitted
    range picked_hz (the F-J picks) is reported against data.velocity_range_kms.  -> (decisions, warnings)"""
    put = _Put(cfg)
    po = int(cfg["dispersion"]["per_octave"])
    f, c, u = np.asarray(freqs), np.asarray(c), np.asarray(u)
    lam = float(cfg["dispersion"]["qc"]["lam_c"])
    dmax = float(geom["max_distance_km"])
    fmax = float(cfg["dispersion"]["fmax_hz"])
    short = f * dmax < lam * c * (1 - 1e-9)
    i0 = int(np.nonzero(short)[0].max()) + 1 if short.any() else 0
    if i0 >= len(f):
        raise RuntimeError(
            "the longest pair (%.3f km) is shorter than %g wavelengths at every frequency up to %.3g Hz"
            % (dmax, lam, f[-1])
        )
    fmin = put(
        ("dispersion", "fmin_hz"),
        float(_grid_up(f[i0] * (1 - 1e-9), po)),
        "the longest pair %.3f km = %g wavelengths of the reference curve (%.3f km/s); below, no pair is long "
        "enough" % (dmax, lam, float(c[i0])),
    )
    put(("export", "tmax_s"), 1 / fmin, "1 / dispersion.fmin_hz = 1 / %.4g Hz" % fmin)
    band = (f >= fmin * (1 - 1e-6)) & (f <= fmax * (1 + 1e-6))
    cb, ub = c[band], u[band]
    vlo, vhi = float(min(cb.min(), ub.min())), float(max(cb.max(), ub.max()))
    steps = (1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8)
    vmin = put(
        ("dispersion", "vmin"),
        round(nice_number(0.6 * vlo, "down", steps), 4),
        "0.6 x the slowest reference velocity %.3f km/s in %.3g-%.3g Hz (c %.3f-%.3f, U %.3f-%.3f)"
        % (vlo, fmin, fmax, cb.min(), cb.max(), ub.min(), ub.max()),
    )
    put(
        ("dispersion", "vmax"),
        round(nice_number(1.3 * vhi, "up", steps), 4),
        "1.3 x the fastest reference velocity %.3f km/s" % vhi,
    )
    put(
        ("dispersion", "dv"),
        nice_number(vmin / 100, "down"),
        "the 1-2-2.5-5 step not above vmin / 100 = %.3g km/s" % (vmin / 100),
    )
    warn = []
    v_lo, v_hi = (float(x) for x in cfg["data"]["velocity_range_kms"])
    if picked_hz is not None:
        m = (f >= picked_hz[0] * (1 - 1e-9)) & (f <= picked_hz[1] * (1 + 1e-9))
        cm = c[m] if m.any() else c
        if cm.min() < 0.8 * v_lo or cm.max() > 1.25 * v_hi:
            warn.append(
                "the F-J phase velocities %.3f-%.3f km/s lie outside data.velocity_range_kms %s by more than 20 %%: "
                "the lag window and the band were set for that range; correct the range and rerun from clean"
                % (cm.min(), cm.max(), cfg["data"]["velocity_range_kms"])
            )
    return put.out, warn


def write(cfg, decisions, geom, extra=None):
    """Write (or update) .lant/auto_parameters.json: parameters merged, geometry and extra blocks replaced."""
    p = path(cfg)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    old = json.load(open(p)) if os.path.exists(p) else {}
    par = old.get("parameters", {})
    par.update(decisions)
    new = dict(old, parameters=par)
    if geom is not None:
        new["geometry"] = geom
    for k, v in (extra or {}).items():
        new[k] = v
    json.dump(new, open(p, "w"), indent=1, ensure_ascii=False, default=str)
    return p


def reset(cfg, keys):
    """Remove derived values (when the step that derives them runs again)."""
    p = path(cfg)
    if not os.path.exists(p):
        return
    d = json.load(open(p))
    for k in keys:
        d.get("parameters", {}).pop(".".join(k), None)
    json.dump(d, open(p, "w"), indent=1, ensure_ascii=False, default=str)


def load(cfg):
    p = path(cfg)
    return json.load(open(p)) if os.path.exists(p) else {}


def apply(cfg):
    """The configuration with the "auto" values resolved: a value written by the user stays; else the confirmed
    value while its confirmation holds (lant.confirm.valid_values); else the rule's value from
    .lant/auto_parameters.json."""
    from .confirm import valid_values

    dec = load(cfg).get("parameters", {})
    conf = valid_values(cfg)
    for key in set(dec) | set(conf):
        keys = key.split(".")
        cur = cfg
        try:
            for k in keys[:-1]:
                cur = cur[k]
        except (KeyError, TypeError):
            continue
        if key in conf and user_value(cfg, keys) is None:  # confirmed, and not written by the user
            cur[keys[-1]] = conf[key]
        elif C.is_auto(cur.get(keys[-1])) and key in dec and not C.is_auto(dec[key]["value"]):
            cur[keys[-1]] = dec[key]["value"]
    C.validate(cfg)
    return cfg


def restore(cfg, keys):
    """Put the given keys back to the user's value, or to the default ("auto"), before their rule is applied again
    (a step that re-derives them must not keep the values of its previous run)."""
    base = C.defaults(cfg["project"]["array_type"])
    for ks in keys:
        uv = user_value(cfg, ks)
        cur, b = cfg, base
        for k in ks[:-1]:
            cur, b = cur[k], b[k]
        cur[ks[-1]] = uv if uv is not None else b[ks[-1]]
    return cfg


def unresolved(cfg, keys=KEYS):
    """Keys still "auto" (other than those resolved by later steps)."""
    left = []
    for ks in keys:
        cur = cfg
        for k in ks:
            cur = cur[k]
        if C.is_auto(cur) and ks not in LATER:
            left.append(".".join(ks))
    return left


def vertical_tests(cfg, depths):
    """checks.vertical "auto": from the depth nodes of the built project."""
    zmax = float(depths[-2])
    ap = load(cfg).get("geometry", {}).get("aperture_km")
    if ap is None:
        raise RuntimeError(".lant/auto_parameters.json has no geometry: run `lant run clean` first")
    cell = float("%.2g" % ((0.34 if cfg["project"]["array_type"] == "line" else 0.18) * ap))
    rev = [round(f * zmax, 3) for f in (0.21, 0.46, 0.71)]
    across = 0.0 if cfg["project"]["array_type"] == "line" else cell
    tests = [{"name": "cb_vertical", "cell_km": [cell, across], "reversal_depths_km": rev}]
    if cfg["project"]["array_type"] == "line":
        tests.insert(
            0, {"name": "cb_vertical_thick", "cell_km": [cell, 0.0], "reversal_depths_km": [round(0.46 * zmax, 3)]}
        )
    return tests
