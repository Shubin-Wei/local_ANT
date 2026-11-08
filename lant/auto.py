"""Parameters given as "auto": derived from the data by written rules after the clean step.

The rules (docs/parameters.md gives the reasons); the values and the numbers they came from are written to
<output>/records/auto_parameters.json and applied before every later step.  A value written in the configuration
always wins over the rule.

  project.array_type             line when the station layout is elongated (principal-axis ratio > 4), else surface
  data.day_start_utc_hour        the UTC hour at which most records longer than 12 h start (0 when they are spread)
  correlation.fs                 the largest of 50, 40, 25, 20, 10, 5 Hz that divides every sampling rate
  correlation.fmax_hz            min(12.5, 0.4 fs)
  dispersion.fmin_hz             the 90th-percentile distance is 2 wavelengths at 2.5 km/s, rounded down to the
                                 1/12-octave grid through 1 Hz (lower bound 0.05 Hz)
  dispersion.fmax_hz             0.84 x correlation.fmax_hz
  export.tmax_s                  1 / dispersion.fmin_hz
  correlation.maxlag_s           max distance / 1.2 km/s + 2 / dispersion.fmin_hz, rounded up to 10 s
  correlation.seg_len_s          at least 2 x maxlag (rounded up to 60 s)
  inversion.grid_spacing_deg     areal: half the median nearest-station distance, to a multiple of 0.0025 deg (at least
                                 0.0025), times grid_refine for the coarse grid; lines: 0.005 deg
  checks.checkerboards           cells of 0.12, 0.18, 0.25 x the aperture (lines: 0.165, 0.24, 0.34 x the line length,
                                 along the line), to 0.1 km
  checks.vertical                cell 0.18 (areal) / 0.34 (line) x aperture; reversals at 0.21, 0.46, 0.71 x the deepest
                                 updated node (lines: a second test with one reversal at 0.46); set by the inversion step
  plots.sections                 areal: two sections through the centre of the stations along the principal axes
"""
import json
import math
import os
from collections import Counter

import numpy as np

from . import config as C
from .signal import geodesic_km, local_xy_km

KEYS = [("project", "array_type"), ("data", "day_start_utc_hour"), ("correlation", "fs"), ("correlation", "fmax_hz"),
        ("correlation", "maxlag_s"), ("correlation", "seg_len_s"), ("dispersion", "fmin_hz"), ("dispersion", "fmax_hz"),
        ("export", "tmax_s"), ("inversion", "grid_spacing_deg"), ("checks", "checkerboards"), ("checks", "vertical"),
        ("plots", "sections")]


def path(cfg):
    return os.path.join(cfg["project"]["output_dir"], "records", "auto_parameters.json")


def geometry(stations):
    names = sorted(stations)
    lon = np.array([stations[n]["lon"] for n in names])
    lat = np.array([stations[n]["lat"] for n in names])
    iu, ju = np.triu_indices(len(names), 1)
    d, _ = geodesic_km(lon[iu], lat[iu], lon[ju], lat[ju])
    nn = np.array([d[(iu == i) | (ju == i)].min() for i in range(len(names))])
    x, y = local_xy_km(lon, lat)
    w, v = np.linalg.eigh(np.cov(np.c_[x, y].T))
    u_major = v[:, np.argmax(w)]
    s = x * u_major[0] + y * u_major[1]
    return dict(n=len(names), d=d, nn=nn, elongation=float(np.sqrt(w.max() / max(w.min(), 1e-12))),
                aperture=float(s.max() - s.min()), axes=v, lon=lon, lat=lat, x=x, y=y)


def user_value(cfg, keys):
    u = cfg["_user"]
    for k in keys:
        if not isinstance(u, dict) or k not in u:
            return None
        u = u[k]
    return u


def infer(cfg, stations, traces):
    """-> {dotted key: {value, reason}} for every key that is "auto" in the configuration."""
    g = geometry(stations)
    out = {}

    def put(keys, value, reason):
        uv = user_value(cfg, keys)
        if uv is not None and not C.is_auto(uv):
            out[".".join(keys)] = dict(value=uv, reason="given in the configuration")
        else:
            out[".".join(keys)] = dict(value=value, reason=reason)
        return out[".".join(keys)]["value"]

    at = put(("project", "array_type"), "line" if g["elongation"] > 4 else "surface",
             "principal-axis ratio of the station layout %.1f (line if > 4)" % g["elongation"])
    long_rec = [t for t in traces if t["end"] - t["start"] >= 12 * 3600]
    import obspy
    hours = Counter(obspy.UTCDateTime(t["start"] + 30).hour for t in long_rec)
    if hours and hours.most_common(1)[0][1] >= 0.6 * len(long_rec):
        h, n = hours.most_common(1)[0]
        put(("data", "day_start_utc_hour"), int(h), "%d of %d records longer than 12 h start at %02d h UTC" % (n, len(long_rec), h))
    else:
        put(("data", "day_start_utc_hour"), 0, "no common start hour of the long records (%s): UTC days" % dict(hours))
    rates = sorted({t["fs"] for t in traces})
    fs = next((c for c in (50, 40, 25, 20, 10, 5) if all(abs(r / c - round(r / c)) < 1e-6 and round(r / c) >= 1 for r in rates)), None)
    if fs is None:
        raise ValueError("no processing rate 5-50 Hz divides every sampling rate %s: set correlation.fs" % rates)
    fs = put(("correlation", "fs"), float(fs), "sampling rates %s Hz" % rates)
    fmax_c = put(("correlation", "fmax_hz"), float(min(12.5, 0.4 * fs)), "min(12.5, 0.4 x fs)")
    p90 = float(np.percentile(g["d"], 90))
    freq = 2 * 2.5 / p90
    fmin = max(0.05, 2 ** (math.floor(12 * math.log2(freq) + 1e-9) / 12))
    fmin = round(put(("dispersion", "fmin_hz"), float(round(fmin, 6)),
                     "90 %% of the distances <= %.2f km = 2 wavelengths at 2.5 km/s -> %.3f Hz, on the 1/12-octave grid" % (p90, freq)), 6)
    put(("dispersion", "fmax_hz"), float(round(0.84 * fmax_c, 4)), "0.84 x correlation.fmax_hz")
    put(("export", "tmax_s"), float(round(1 / fmin, 4)), "1 / dispersion.fmin_hz")
    ml = put(("correlation", "maxlag_s"), float(10 * math.ceil((g["d"].max() / 1.2 + 2 / fmin) / 10)),
             "max distance %.2f km / 1.2 km/s + 2 periods of %.3f Hz, up to 10 s" % (g["d"].max(), fmin))
    put(("correlation", "seg_len_s"), float(max(60, 60 * math.ceil(2 * ml / 60))), "60 s, at least 2 x maxlag")
    if at == "line":
        put(("inversion", "grid_spacing_deg"), [0.005, 0.005], "line arrays: 0.005 deg (default for lines)")
    else:
        r = cfg["inversion"]["grid_refine"]
        fine = max(0.0025, 0.0025 * round(np.median(g["nn"]) / 2 / 111.195 / 0.0025))
        put(("inversion", "grid_spacing_deg"), [round(fine * r, 6)] * 2,
            "median nearest-station distance %.3f km / 2 -> %.4f deg (x grid_refine %d)" % (np.median(g["nn"]), fine, r))
    ap = g["aperture"]
    frac = (0.165, 0.24, 0.34) if at == "line" else (0.12, 0.18, 0.25)
    cells = [round(f * ap, 1) for f in frac]
    put(("checks", "checkerboards"),
        [{"name": "cb_%gkm" % c, "cell_km": [c, 0.0 if at == "line" else c]} for c in cells],
        "%s x the aperture %.1f km" % (" / ".join("%g" % f for f in frac), ap))
    put(("checks", "vertical"), "auto", "set by the inversion step from the depth nodes")
    if at == "surface":
        k = np.radians(1.0) * 6371.0
        lat0, lon0 = float(np.mean(g["lat"])), float(np.mean(g["lon"]))
        secs = []
        for name, u in zip("AB", (g["axes"][:, 1], g["axes"][:, 0])):
            s = g["x"] * u[0] + g["y"] * u[1]
            ends = []
            for sv in (s.min(), s.max()):
                ends.append([round(lon0 + sv * u[0] / (k * math.cos(math.radians(lat0))), 5), round(lat0 + sv * u[1] / k, 5)])
            secs.append({"name": name, "start": ends[0], "end": ends[1]})
        put(("plots", "sections"), secs, "through the centre of the stations along the principal axes (A: long axis)")
    else:
        put(("plots", "sections"), [], "line arrays: section along the line")
    return out, dict(stations=g["n"], elongation=g["elongation"], aperture_km=ap, p90_distance_km=p90,
                     max_distance_km=float(g["d"].max()), median_nearest_km=float(np.median(g["nn"])))


def vertical_tests(cfg, depths):
    """checks.vertical "auto": from the depth nodes of the built project."""
    zmax = float(depths[-2])
    ap = None
    p = path(cfg)
    if os.path.exists(p):
        ap = json.load(open(p))["geometry"]["aperture_km"]
    cell = round((0.34 if cfg["project"]["array_type"] == "line" else 0.18) * (ap or 10.0), 1)
    rev = [round(f * zmax, 1) for f in (0.21, 0.46, 0.71)]
    across = 0.0 if cfg["project"]["array_type"] == "line" else cell
    tests = [{"name": "cb_vertical", "cell_km": [cell, across], "reversal_depths_km": rev}]
    if cfg["project"]["array_type"] == "line":
        tests.insert(0, {"name": "cb_vertical_thick", "cell_km": [cell, 0.0], "reversal_depths_km": [round(0.46 * zmax, 1)]})
    return tests


def write(cfg, decisions, geom):
    p = path(cfg)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    json.dump(dict(parameters=decisions, geometry=geom), open(p, "w"), indent=1, ensure_ascii=False, default=str)
    return p


def apply(cfg):
    """The configuration with the "auto" values resolved (needs records/auto_parameters.json from the clean step)."""
    p = path(cfg)
    if not os.path.exists(p):
        return cfg
    dec = json.load(open(p))["parameters"]
    at = dec.get("project.array_type", {}).get("value", "surface")
    if cfg["project"]["array_type"] == "auto":
        cfg = C.rebuild(cfg, at)
        cfg["project"]["array_type"] = at
    for key, d in dec.items():
        keys = key.split(".")
        cur = cfg
        for k in keys[:-1]:
            cur = cur[k]
        if C.is_auto(cur.get(keys[-1])) and not C.is_auto(d["value"]):
            cur[keys[-1]] = d["value"]
    C.validate(cfg)
    return cfg


def unresolved(cfg):
    """Keys still "auto" (other than those resolved by later steps)."""
    left = []
    for keys in KEYS:
        cur = cfg
        for k in keys:
            cur = cur[k]
        if C.is_auto(cur) and keys not in (("checks", "vertical"),):
            left.append(".".join(keys))
    return left
