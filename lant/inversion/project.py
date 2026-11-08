"""DSurfTomo project files (Fang et al. 2015; DSurfTomo manual v1.4).

  surfdata.dat   one block per (period, source station): "# lat lon period_index wave_type vel_type", then
                 "lat lon velocity" of every receiver; period indices count the sorted periods that carry data
  DSurfTomo.in   grid, regularisation, velocity bounds, iterations, sparsity, period list, synthetic flag, noise,
                 threshold (the layout of DSurfTomo's own GenerateDSurfTomoInputFile.py)
  MOD            depths, then for every depth ny lines of nx velocities (laterally uniform 1-D model)
  MOD.true       checkerboards only: the true model, same layout without the depth line

Grid: the extent follows the stations that carry data, with a margin of max(1 % of the range, 3 cells) on each side
(no station may fall into the boundary ring), spacing inversion.grid_spacing_deg (lat, lon); then every interval is
divided by grid_refine (same origin and extent).  Lines: across the line only line_strip_cells cells are kept on
each side of the stations (DSurfTomo computes travel times from the second node over ny-3 cells; a station on the
edge of that strip gets wrong travel times without an error).
"""
import json
import math
import os
import shutil
from collections import defaultdict

import numpy as np

KM_PER_DEG = 2 * np.pi * 6371.0 / 360.0
MAX_PERIODS = 80


def records_from_curves(curves):
    recs = []
    for c in curves:
        recs.append(dict(sta1=(c["lat1"], c["lon1"]), sta2=(c["lat2"], c["lon2"]),
                         values={round(float(t), 4): float(v) for t, v in zip(c["T"], c["v"])}, name=c["name"]))
    periods = sorted({t for r in recs for t in r["values"]})
    if len(periods) > MAX_PERIODS:
        raise RuntimeError("%d periods carry data; DSurfTomo (surfdisp96) accepts at most %d: narrow export.tmin_s / "
                           "tmax_s or use fewer frequencies per octave" % (len(periods), MAX_PERIODS))
    idx = {t: k + 1 for k, t in enumerate(periods)}
    for r in recs:
        r["index"] = {idx[t]: v for t, v in r["values"].items()}
    return recs, periods


def write_surfdata(recs, path, wave_type=2, vel_type=0):
    grouped = defaultdict(list)
    stations = set()
    for r in recs:
        stations.add(r["sta1"])
        stations.add(r["sta2"])
        for k, v in r["index"].items():
            grouped[(k, r["sta1"])].append((r["sta2"], v))
    n_blocks, max_rcv = 0, 0
    src_per_period = defaultdict(set)
    with open(path, "w") as fh:
        for key in sorted(grouped, key=lambda it: (it[0], it[1][0], it[1][1])):
            k, (lat1, lon1) = key
            rcv = sorted(grouped[key], key=lambda it: (it[0][0], it[0][1]))
            n_blocks += 1
            max_rcv = max(max_rcv, len(rcv))
            src_per_period[k].add((lat1, lon1))
            fh.write("# %.6f %.6f %d %d %d\n" % (lat1, lon1, k, wave_type, vel_type))
            for (lat2, lon2), v in rcv:
                fh.write("%.6f %.6f %.6f\n" % (lat2, lon2, v))
    nsrc = max(max_rcv, max(len(v) for v in src_per_period.values()), len(stations))
    return dict(blocks=n_blocks, nsrc=nsrc, travel_times=sum(len(r["index"]) for r in recs), stations=sorted(stations))


def grid(points, spacing, refine=1, strip_cells=None):
    """points [(lat, lon)] -> dict(nx, ny, goxd, gozd, dvxd, dvzd) (x = latitude, z = longitude, as in DSurfTomo)."""
    lats = [p[0] for p in points]
    lons = [p[1] for p in points]
    lat_min, lat_max, lon_min, lon_max = min(lats), max(lats), min(lons), max(lons)
    lat_rng, lon_rng = max(lat_max - lat_min, 0.01), max(lon_max - lon_min, 0.01)
    dvx, dvz = spacing
    eps_x = max(max(lat_rng * 0.01, 1e-4), 3 * dvx)
    eps_z = max(max(lon_rng * 0.01, 1e-4), 3 * dvz)
    goxd, gozd = lat_max + eps_x, lon_min - eps_z
    nx = max(int(((lat_rng + 2.0 * eps_x) / dvx) + 0.999999) + 3, 4)
    ny = max(int(((lon_rng + 2.0 * eps_z) / dvz) + 0.999999) + 3, 4)
    g = dict(nx=nx, ny=ny, goxd=goxd, gozd=gozd, dvxd=dvx, dvzd=dvz, coarse_dvxd=dvx, coarse_dvzd=dvz, refine=refine)
    if refine > 1:
        g.update(nx=(nx - 1) * refine + 1, ny=(ny - 1) * refine + 1, dvxd=dvx / refine, dvzd=dvz / refine)
    if strip_cells is not None:
        km_x = (lat_max - lat_min) * KM_PER_DEG
        km_z = (lon_max - lon_min) * KM_PER_DEG * math.cos(math.radians(0.5 * (lat_min + lat_max)))
        if km_x >= km_z:      # line runs north-south: narrow in longitude
            n2 = int(math.ceil((lon_max - lon_min) / g["dvzd"] + 2 * strip_cells - 1e-9)) + 3
            g.update(ny=n2, gozd=0.5 * (lon_min + lon_max) - 0.5 * (n2 - 3) * g["dvzd"], strip="longitude")
        else:                 # east-west: narrow in latitude
            n2 = int(math.ceil((lat_max - lat_min) / g["dvxd"] + 2 * strip_cells - 1e-9)) + 3
            g.update(nx=n2, goxd=0.5 * (lat_min + lat_max) + 0.5 * (n2 - 3) * g["dvxd"], strip="latitude")
    return g


def node_coords(g):
    lat = g["goxd"] - np.arange(g["nx"]) * g["dvxd"]
    lon = g["gozd"] + np.arange(g["ny"]) * g["dvzd"]
    return lat, lon


def write_in(path, g, nz, nsrc, inv, periods, synthetic=0):
    def line(v, c):
        return "%-33s%s\n" % (v, "c: " + c)
    with open(path, "w") as fh:
        fh.write("c" * 69 + "\nc INPUT PARAMETERS\n" + "c" * 69 + "\n")
        fh.write(line("surfdata.dat", "data file"))
        fh.write(line("%d %d %d" % (g["nx"], g["ny"], nz), "nx ny nz (grid number in lat lon and depth direction)"))
        fh.write(line("%.6f %.6f" % (g["goxd"], g["gozd"]), "goxd gozd (upper left point,[lat,lon])"))
        fh.write(line("%.6f %.6f" % (g["dvxd"], g["dvzd"]), "dvxd dvzd (grid interval in lat and lon direction)"))
        fh.write(line("%d" % nsrc, "max(sources, receivers)"))
        fh.write(line("%.3f %.3f" % (inv["weight"], inv["damp"]), "weight damp"))
        fh.write(line("%d" % inv["sublayers"], "sublayers (for computing depth kernel, 2~5)"))
        fh.write(line("%.3f %.3f" % (inv["vmin"], inv["vmax"]), "minimum velocity, maximum velocity (a priori information)"))
        fh.write(line("%d" % inv["iterations"], "maximum iteration"))
        fh.write(line("%.3f" % inv["sparsity"], "sparsity fraction"))
        fh.write(line("%d" % len(periods), "kmaxRc (followed by periods)"))
        fh.write(" ".join("%.4f" % t for t in periods) + "\n")
        fh.write(line("0", "kmaxRg"))
        fh.write(line("0", "kmaxLc"))
        fh.write(line("0", "kmaxLg"))
        fh.write(line("%d" % synthetic, "synthetic flag(0:real data,1:synthetic)"))
        fh.write(line("%.3f" % inv["noise"], "noiselevel"))
        fh.write(line("%.3f" % inv["threshold"], "threshold"))


def write_mod(path, z, vs, g, header=True, pattern=None, depth_factor=None):
    """1-D model on every node; with pattern [ny, nx] (km/s) added, scaled per depth by depth_factor."""
    with open(path, "w") as f:
        if header:
            f.write(" ".join("%7.3f" % d for d in z) + "\n")
        for k in range(len(z)):
            s = 1.0 if depth_factor is None else depth_factor[k]
            for j in range(g["ny"]):
                if pattern is None:
                    f.write("".join("%7.3f" % vs[k] for _ in range(g["nx"])) + "\n")
                else:
                    f.write("".join("%7.3f" % (vs[k] + s * pattern[j, i]) for i in range(g["nx"])) + "\n")


def checker_pattern(g, cell_km, amp, array_type, line_dir=None):
    """[ny, nx] pattern.  Areal arrays: sin(a_lat * i) sin(a_lon * j) on the node index of the coarse grid (cell size
    cell_km = [north-south, east-west] km; 0 = constant in that direction).  Lines: sin(pi s / cell) of the distance
    s along the line (cell_km[0]), constant across it unless cell_km[1] > 0."""
    nx, ny = g["nx"], g["ny"]
    r = g["refine"]
    xi = (np.arange(1, nx + 1) - 1) / r + 1          # index on the coarse grid
    yj = (np.arange(1, ny + 1) - 1) / r + 1
    latc = g["goxd"] - (nx - 1) / 2.0 * g["dvxd"]
    if array_type == "line":
        lat, lon = node_coords(g)
        ys = (g["goxd"] - lat) * KM_PER_DEG                                  # km south of the grid origin
        xs = (lon - g["gozd"]) * KM_PER_DEG * math.cos(math.radians(latc))   # km east
        u = np.asarray(line_dir, float) / np.hypot(*line_dir)               # (south, east) unit vector along the line
        s = u[0] * ys[None, :] + u[1] * xs[:, None]                         # [ny, nx]
        t = -u[1] * ys[None, :] + u[0] * xs[:, None]
        p = np.sin(np.pi * s / cell_km[0]) if cell_km[0] else np.ones_like(s)
        if cell_km[1]:
            p = p * np.sin(np.pi * t / cell_km[1])
        return amp * p
    a_lat = np.pi * g["coarse_dvxd"] * KM_PER_DEG / cell_km[0] if cell_km[0] else 0.0
    a_lon = np.pi * g["coarse_dvzd"] * KM_PER_DEG * math.cos(math.radians(latc)) / cell_km[1] if cell_km[1] else 0.0
    px = np.sin(a_lat * xi) if a_lat else np.ones(nx)
    py = np.sin(a_lon * yj) if a_lon else np.ones(ny)
    return amp * np.outer(py, px)


def depth_signs(z, seq, reversals):
    """+1 / -1 on the nodes of the standard sequence (sign changes at each reversal depth), linear on added nodes."""
    zs = np.array([v for v in seq if v <= max(z) + 1e-9], float)
    sg = np.array([(-1.0) ** int(np.sum(np.asarray(reversals) < v)) for v in zs])
    return np.interp(z, zs, sg)


def line_direction(points):
    """Principal direction of the station positions as a (south, east) unit vector."""
    lat = np.array([p[0] for p in points])
    lon = np.array([p[1] for p in points])
    y = -(lat - lat.mean()) * KM_PER_DEG
    x = (lon - lon.mean()) * KM_PER_DEG * math.cos(math.radians(lat.mean()))
    w, v = np.linalg.eigh(np.cov(np.c_[y, x].T))
    u = v[:, np.argmax(w)]
    return u if u[0] >= 0 else -u


def build_main(cfg, curves, z, vs, out_dir):
    from .start import check_layers
    inv = cfg["inversion"]
    check_layers(len(z), inv["sublayers"])
    recs, periods = records_from_curves(curves)
    os.makedirs(out_dir, exist_ok=True)
    info = write_surfdata(recs, os.path.join(out_dir, "surfdata.dat"))
    pts = info.pop("stations")
    g = grid(pts, inv["grid_spacing_deg"], inv["grid_refine"],
             inv["line_strip_cells"] if cfg["project"]["array_type"] == "line" else None)
    nparam = g["nx"] * g["ny"] * len(z)
    maxnar = int(inv["sparsity"] * info["travel_times"] * nparam)
    if 2 * maxnar + 1 > 2147483647:
        raise RuntimeError("DSurfTomo's sparse matrix would overflow (%d entries): lower inversion.sparsity or coarsen the "
                           "grid" % maxnar)
    write_in(os.path.join(out_dir, "DSurfTomo.in"), g, len(z), info["nsrc"], inv, periods)
    write_mod(os.path.join(out_dir, "MOD"), z, vs, g)
    np.savetxt(os.path.join(out_dir, "vs1d.txt"), np.c_[z, vs], fmt="%10.4f %10.4f")
    meta = dict(grid=g, periods=periods, curves=len(recs), travel_times=info["travel_times"], nsrc=info["nsrc"],
                depths=list(map(float, z)), weight=inv["weight"], damp=inv["damp"],
                line_direction=list(map(float, line_direction(pts))) if cfg["project"]["array_type"] == "line" else None,
                stations_with_data=len(pts))
    json.dump(meta, open(os.path.join(out_dir, "project.json"), "w"), indent=1)
    return meta


def copy_project(src, dst, files=("DSurfTomo.in", "MOD", "surfdata.dat", "vs1d.txt", "project.json")):
    if os.path.exists(os.path.join(dst, "DSurfTomo.in.log")):
        return False
    shutil.rmtree(dst, ignore_errors=True)
    os.makedirs(dst)
    for f in files:
        shutil.copyfile(os.path.join(src, f), os.path.join(dst, f))
    return True


def set_line(path, key, value):
    L = open(path).read().split("\n")
    k = next(i for i, l in enumerate(L) if "c: " + key in l)
    L[k] = "%-33s" % value + L[k][L[k].index("c: "):]
    open(path, "w").write("\n".join(L))


def build_checks(cfg, main_dir, out_root):
    """Checkerboards, depth-varying checkerboards and the control inversions, built from the main project."""
    meta = json.load(open(os.path.join(main_dir, "project.json")))
    g = meta["grid"]
    m = np.loadtxt(os.path.join(main_dir, "vs1d.txt"))
    z, vs = m[:, 0], m[:, 1]
    ck = cfg["checks"]
    at = cfg["project"]["array_type"]
    built = []
    vertical = ck["vertical"]
    if vertical == "auto":
        from ..auto import vertical_tests
        vertical = vertical_tests(cfg, z)
    specs = [(c["name"], c["cell_km"], None) for c in ck["checkerboards"]] + \
            [(c["name"], c["cell_km"], c["reversal_depths_km"]) for c in vertical]
    for name, cell, rev in specs:
        d = os.path.join(out_root, "checkerboards", name)
        if not copy_project(main_dir, d):
            built.append((name, "kept (already run)"))
            continue
        pat = checker_pattern(g, cell, ck["amplitude"], at, meta.get("line_direction"))
        fac = depth_signs(z, cfg["inversion"]["depth_sequence_km"], rev) if rev else None
        write_mod(os.path.join(d, "MOD.true"), z, vs, g, header=False, pattern=pat, depth_factor=fac)
        set_line(os.path.join(d, "DSurfTomo.in"), "synthetic flag", "1")
        built.append((name, "built"))
    for name in ck["controls"]:
        d = os.path.join(out_root, "controls", name)
        if not copy_project(main_dir, d):
            built.append((name, "kept (already run)"))
            continue
        if name.startswith("reg_w"):
            w, dp = name[5:].split("_d")
            set_line(os.path.join(d, "DSurfTomo.in"), "weight damp", "%.3f %.3f" % (float(w), float(dp)))
        elif name.startswith("start_x"):
            f = float(name[7:])
            M = open(os.path.join(main_dir, "MOD")).read().split("\n")
            out = [M[0]] + ["".join("%7.3f" % (f * float(v)) for v in l.split()) for l in M[1:] if l.strip()]
            open(os.path.join(d, "MOD"), "w").write("\n".join(out) + "\n")
        else:
            raise RuntimeError("checks.controls: unknown control %s (use reg_w<weight>_d<damp> or start_x<factor>)" % name)
        built.append((name, "built"))
    return dict(built=built, vertical=vertical)
