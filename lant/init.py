"""lant init: a configuration file with the numbers suggested for this array.

The suggestions follow from the station geometry, the sampling rates of the records and the velocity range the user
gives; each one is written into the file with the rule that gave it, and every one may be changed:

  fmin  two-station measurements need the path to be at least disp.min_wavelengths (2) wavelengths long (Bensen
        et al. 2007; Luo et al. 2015): the longest pair r_max measures down to f = 2 v_fast / r_max; the correlation
        band starts lower by the half-gain width of the Gaussian filter, x (1 - sqrt(ln 2 / alpha));
  fmax  spatial sampling theorem: the stations sample a wave without aliasing when their spacing is at most half
        its wavelength (e.g. Rost & Thomas 2002), here the median distance to the nearest station dx:
        f = v_slow / (2 dx), widened by x (1 + sqrt(ln 2 / alpha)); time sampling: at most 0.4 x the lowest raw
        sampling rate (the anti-alias low-pass is at 0.8 x Nyquist);
  fs    the lowest whole-number rate reachable from every raw rate by an integer factor with fmax <= 0.4 fs (the
        anti-alias low-pass is at 0.8 x Nyquist);
  maxlag  r_max / (0.6 v_slow) + 4 / fmin: the group velocity of strongly dispersive layers can be well below the
        slowest phase velocity (0.6: an empirical margin), and four periods after the arrival give the noise window;
  segment  max(2.5 maxlag, 10 / fmin): longer than twice the lag window (circular correlation) and ten periods;
  substack_hours  24 (one sub-stack per day) when the records span at least 3 days, else the longest of 12, 8, 6,
        4, 3, 2, 1 h that gives 3 sub-stacks: the jackknife error of lant disp needs at least 3.
"""

import math

import numpy as np

from . import config as C
from . import geo, stations
from . import records as R

HALF = math.sqrt(math.log(2) / 5.0)  # sqrt(ln 2 / alpha) for the default alpha = 5


def _up(x, step):
    return step * math.ceil(x / step - 1e-9)


def suggest(table, rates, vrange, min_wavelengths=2.0):
    """-> {key: (value, rule)}."""
    names = sorted(table)
    lon = np.array([table[n]["lon"] for n in names])
    lat = np.array([table[n]["lat"] for n in names])
    iu, ju = np.triu_indices(len(names), 1)
    d = np.atleast_1d(geo.geodesic_km(lon[iu], lat[iu], lon[ju], lat[ju])[0])
    rmax, dnn = float(d.max()), float(np.median(geo.nearest_neighbour_km(lon, lat)))
    v0, v1 = float(vrange[0]), float(vrange[1])
    out = {}
    f_lo = min_wavelengths * v1 / rmax
    fmin = float("%.3g" % (f_lo * (1 - HALF)))
    out["cc.fmin_hz"] = (fmin, "%g x v_fast %g km/s / longest pair %.3f km = %.4g Hz, x (1 - %.3f)"
                         % (min_wavelengths, v1, rmax, f_lo, HALF))  # fmt: skip
    f_hi = v0 / (2 * dnn) if dnn > 0 else np.inf
    raw_cap = 0.4 * min(rates)
    fmax = float("%.3g" % min(f_hi * (1 + HALF), raw_cap))
    out["cc.fmax_hz"] = (fmax, "v_slow %g km/s / (2 x median station spacing %.4g km) = %.4g Hz, x (1 + %.3f); at "
                         "most 0.4 x the lowest raw rate (%g Hz)" % (v0, dnn, f_hi, HALF, raw_cap))  # fmt: skip
    if fmax <= fmin * 1.5:
        raise ValueError(
            "the suggested band %.3g-%.3g Hz is too narrow: the array is small for velocities up to %g km/s "
            "(the longest pair is %.3f km); check the velocity range" % (fmin, fmax, v1, rmax)
        )
    lo = min(rates)
    cand = [lo / q for q in range(1, int(lo) + 1) if abs(lo / q - round(lo / q)) < 1e-6
            and all(abs(r / (lo / q) - round(r / (lo / q))) < 1e-6 for r in rates)]  # fmt: skip
    ok = [f for f in cand if 0.4 * f >= fmax - 1e-9]
    fs = float(round(min(ok))) if ok else float(lo)
    out["clean.fs_hz"] = (fs, "lowest whole-number rate reachable from the raw rates %s with fmax <= 0.4 fs"
                       % sorted(set(rates)))  # fmt: skip
    lag = _up(rmax / (0.6 * v0) + 4 / fmin, 1.0 if rmax / v0 < 100 else 10.0)
    out["cc.maxlag_s"] = (lag, "longest pair %.3f km / (0.6 x %g km/s) + 4 / fmin" % (rmax, v0))
    seg = _up(max(2.5 * lag, 10 / fmin), 10.0)
    out["cc.segment_s"] = (seg, "max(2.5 x maxlag, 10 / fmin), up to 10 s")
    return out, dict(stations=len(names), longest_pair_km=rmax, median_spacing_km=dnn)


def write(path, values, notes):
    """The template with `values` ({dotted key: value}) filled in, the rule of each suggestion above its line."""
    out, section = [], None
    for ln in C.template_text().splitlines():
        s = ln.strip()
        if not s or s.startswith("#") or ":" not in s:
            out.append(ln)
            continue
        key, rest = s.split(":", 1)
        top = not ln.startswith(" ")
        if top:
            section = key if not rest.split("#", 1)[0].strip() else None
            full = key
        else:
            full = "%s.%s" % (section, key)
        if full not in values:
            out.append(ln)
            continue
        indent = ln[: len(ln) - len(ln.lstrip())]
        if full in notes:
            out.append("%s# lant init: %s" % (indent, notes[full]))
        body = "%s%s: %s" % (indent, key, _yaml(values[full]))
        comment = rest.split("#", 1)[1] if "#" in rest else None
        out.append("%-29s #%s" % (body, comment) if comment is not None else body)  # YAML needs a space before #
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


def _yaml(v):
    import yaml

    return yaml.safe_dump(v, default_flow_style=True).strip().replace("\n...", "")


def run(path, raw_dir, station_file, out_dir, vrange, channel="*Z", n_cores=1):
    import os

    table = stations.read_table(station_file)
    print("station table: %d stations" % len(table))
    rows, bad = R.build_index(raw_dir, os.path.join(out_dir, "1_clean", "record_index.csv"), n_cores)
    recs, notes = R.match(rows, table, channel)
    for n in notes:
        print("  " + n)
    if len(recs) < 2:
        raise R.DataError("fewer than two stations of the table have %s records under %s" % (channel, raw_dir))
    rates = sorted({r[0]["fs"] for r in recs.values()})
    sug, geom = suggest({n: table[n] for n in recs}, rates, vrange)
    hours = (max(r["end"] for v in recs.values() for r in v) - min(r["start"] for v in recs.values() for r in v)) / 3600
    sub = next((h for h in (24, 12, 8, 6, 4, 3, 2, 1) if hours / h >= 3), 1)
    sug["cc.substack_hours"] = (sub, "the records span %.1f h: the longest sub-stack giving at least 3 for the "
                                "jackknife" % hours)  # fmt: skip
    values = {"output_dir": os.path.abspath(out_dir), "data.raw_dir": os.path.abspath(raw_dir),
              "data.station_file": os.path.abspath(station_file), "data.channel": channel,
              "data.velocity_range_kms": [float(vrange[0]), float(vrange[1])]}  # fmt: skip
    values.update({k: v for k, (v, _) in sug.items()})
    write(path, values, {k: r for k, (_, r) in sug.items()})
    print("%d stations with records, longest pair %.3f km, median station spacing %.3f km"
          % (geom["stations"], geom["longest_pair_km"], geom["median_spacing_km"]))  # fmt: skip
    for k, (v, rule) in sug.items():
        print("  %s = %s   (%s)" % (k, v, rule))
    print("written: %s\nnext: check it, then `lant clean -c %s`" % (path, path))
