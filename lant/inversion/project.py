"""The input files of the 3-D inversion, in DSurfTomo's formats (Fang et al. 2015; DSurfTomo manual v1.4), so that
the same project can also be run with the original DSurfTomo:

  surfdata.dat   one block per (period, source station): "# lat lon period_index 2 0", then "lat lon c" of each
                 receiver (2 0: Rayleigh, phase velocity); period indices count the sorted periods (from 1)
  DSurfTomo.in   grid, regularization, velocity bounds, iterations, periods, synthetic flag, noise, rejection k
  MOD            the depth nodes, then for every depth ny lines of nx Vs (here a laterally uniform 1-D model)
  MOD.true       resolution test only: the true model, MOD without its depth line

Grid (x = latitude from north to south, y = longitude from west to east, as DSurfTomo): the stations with data plus
a margin of 3 cells on each side; DSurfTomo neither updates nor writes the outermost ring of nodes.
"""

import math
from collections import defaultdict

import numpy as np

KM_PER_DEG = 2 * np.pi * 6371.0 / 360.0  # on DSurfTomo's sphere
MARGIN_CELLS = 3
MAX_PERIODS = 60  # NP of DSurfTomo's CalSurfG


def on_grid(periods, per_octave):
    """Periods whose frequency lies on the 1/per_octave-octave grid through 1 Hz."""
    n = per_octave * np.log2(1.0 / np.asarray(periods, float))
    return np.abs(n - np.round(n)) < 1e-4


def records(table, per_octave):
    """Rows of 3_disp/dispersion.csv -> (records [dict(sta1=(lat, lon), sta2, values {period: c})], periods)."""
    by = defaultdict(dict)
    keep = on_grid(table["period_s"], per_octave)
    for i in np.where(keep)[0]:
        key = (table["lat1"][i], table["lon1"][i], table["lat2"][i], table["lon2"][i])
        by[key][round(float(table["period_s"][i]), 6)] = float(table["c_kms"][i])
    periods = sorted({t for v in by.values() for t in v})
    if len(periods) > MAX_PERIODS:
        raise RuntimeError("%d periods; DSurfTomo takes at most %d: lower inv.per_octave" % (len(periods), MAX_PERIODS))
    idx = {t: k + 1 for k, t in enumerate(periods)}
    recs = [dict(sta1=(k[0], k[1]), sta2=(k[2], k[3]), values=v, index={idx[t]: c for t, c in v.items()})
            for k, v in by.items()]  # fmt: skip
    return recs, periods


def write_surfdata(recs, path):
    grouped = defaultdict(list)
    stations = set()
    for r in recs:
        stations.update([r["sta1"], r["sta2"]])
        for k, v in r["index"].items():
            grouped[(k, r["sta1"])].append((r["sta2"], v))
    max_rcv, per_period = 0, defaultdict(set)
    with open(path, "w") as fh:
        for key in sorted(grouped, key=lambda it: (it[0], it[1][0], it[1][1])):
            k, (lat1, lon1) = key
            rcv = sorted(grouped[key])
            max_rcv = max(max_rcv, len(rcv))
            per_period[k].add((lat1, lon1))
            fh.write("# %.6f %.6f %d 2 0\n" % (lat1, lon1, k))
            for (lat2, lon2), v in rcv:
                fh.write("%.6f %.6f %.6f\n" % (lat2, lon2, v))
    nsrc = max(max_rcv, max(len(v) for v in per_period.values()), len(stations))
    return dict(nsrc=nsrc, travel_times=sum(len(r["index"]) for r in recs), stations=sorted(stations))


def grid(points, spacing_km):
    """points [(lat, lon)] -> dict(nx, ny, goxd, gozd, dvxd, dvzd): equal spacing in km in both directions."""
    lats, lons = [p[0] for p in points], [p[1] for p in points]
    lat0 = 0.5 * (min(lats) + max(lats))
    dvx = float("%.3g" % (spacing_km / KM_PER_DEG))
    dvz = float("%.3g" % (spacing_km / (KM_PER_DEG * math.cos(math.radians(lat0)))))
    nx = int(math.ceil((max(lats) - min(lats)) / dvx - 1e-9)) + 1 + 2 * MARGIN_CELLS
    ny = int(math.ceil((max(lons) - min(lons)) / dvz - 1e-9)) + 1 + 2 * MARGIN_CELLS
    return dict(nx=nx, ny=ny, goxd=max(lats) + MARGIN_CELLS * dvx, gozd=min(lons) - MARGIN_CELLS * dvz, dvxd=dvx,
                dvzd=dvz)  # fmt: skip


def node_coords(g):
    lat = g["goxd"] - np.arange(g["nx"]) * g["dvxd"]
    lon = g["gozd"] + np.arange(g["ny"]) * g["dvzd"]
    return lat, lon


def write_control(path, g, nz, nsrc, inv, periods, synthetic=0, noise=0.0):
    def line(v, c):
        return "%-33s%s\n" % (v, "c: " + c)

    from .start_model import MINTHK

    with open(path, "w") as fh:
        fh.write("c" * 69 + "\nc INPUT PARAMETERS (written by lant)\n" + "c" * 69 + "\n")
        fh.write(line("surfdata.dat", "data file"))
        fh.write(line("%d %d %d" % (g["nx"], g["ny"], nz), "nx ny nz (grid number in lat lon and depth direction)"))
        fh.write(line("%.6f %.6f" % (g["goxd"], g["gozd"]), "goxd gozd (upper left point,[lat,lon])"))
        fh.write(line("%.6f %.6f" % (g["dvxd"], g["dvzd"]), "dvxd dvzd (grid interval in lat and lon direction)"))
        fh.write(line("%d" % nsrc, "max(sources, receivers)"))
        fh.write(line("%.6g %.6g" % (inv["weight"], inv["damping"]), "weight damp (lant: of the model, as means)"))
        fh.write(line("%d" % MINTHK, "sublayers (for computing depth kernel, 2~5)"))
        fh.write(line("%.4f %.4f" % (inv["vs_min"], inv["vs_max"]), "minimum velocity, maximum velocity"))
        fh.write(line("%d" % inv["iterations"], "maximum iteration"))
        fh.write(line("0.100", "sparsity fraction (not used by lant)"))
        fh.write(line("%d" % len(periods), "kmaxRc (followed by periods)"))
        fh.write(" ".join("%.6f" % t for t in periods) + "\n")
        for k in ("kmaxRg", "kmaxLc", "kmaxLg"):
            fh.write(line("0", k))
        fh.write(line("%d" % synthetic, "synthetic flag(0:real data,1:synthetic)"))
        fh.write(line("%.3f" % noise, "noiselevel"))
        fh.write(line("%.3f" % inv["outlier_mad"], "threshold (lant: k of the rejection beyond k x 1.4826 MAD)"))


def write_model(path, z, vs, g, header=True, pattern=None, depth_sign=None):
    """1-D model on every node; pattern [ny, nx] (a fraction) applied as vs (1 + s pattern), s = depth_sign[k]."""
    with open(path, "w") as f:
        if header:
            f.write(" ".join("%7.3f" % d for d in z) + "\n")
        for k in range(len(z)):
            s = 1.0 if depth_sign is None else depth_sign[k]
            for j in range(g["ny"]):
                row = vs[k] * (1 + s * pattern[j]) if pattern is not None else np.full(g["nx"], vs[k])
                f.write(" ".join("%.5f" % v for v in row) + "\n")


def checker_pattern(g, cell_km, amplitude):
    """[ny, nx] checkerboard: amplitude x sin(pi x / cell) sin(pi y / cell), x, y in km from the grid corner."""
    lat, lon = node_coords(g)
    x = (g["goxd"] - lat) * KM_PER_DEG
    y = (lon - g["gozd"]) * KM_PER_DEG * math.cos(math.radians(np.mean(lat)))
    return amplitude * np.outer(np.sin(np.pi * y / cell_km), np.sin(np.pi * x / cell_km))


def read_inverted(path, depths):
    """DSurfTomo's model output (lon lat depth Vs of the updated nodes) -> rows [lon, lat, depth, vs]."""
    a = np.loadtxt(path, ndmin=2)
    return a[:, :4]
