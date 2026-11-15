"""DSurfTomo project files (Fang et al. 2015; DSurfTomo manual v1.4).

  surfdata.dat   one block per (period, source station): "# lat lon period_index wave_type vel_type", then
                 "lat lon velocity" of every receiver; period indices count the sorted periods that carry data
  DSurfTomo.in   grid, regularisation, velocity bounds, iterations, period list, synthetic flag, noise, k of the
                 rejection (the layout of DSurfTomo's own GenerateDSurfTomoInputFile.py; the sparsity line is unused)
  MOD            depths, then for every depth ny lines of nx velocities (laterally uniform 1-D model)
  MOD.true       checkerboards only: the true model, same layout without the depth line

Grid: the extent follows the stations that carry data, with a margin of max(1 % of the range, 3 cells) on each side
(no station may fall into the boundary ring), spacing inversion.grid_spacing_deg (lat, lon; auto: grid_spacing, the
median measured wavelength / 3 in km, the same in both directions).  Lines: across the line only line_strip_cells cells are kept on
each side of the stations (DSurfTomo computes travel times from the second node over ny-3 cells; a station on the
edge of that strip gets wrong travel times without an error).
"""

import json
import math
import os
import shutil
from collections import defaultdict

import numpy as np

from .. import config as C

KM_PER_DEG = 2 * np.pi * 6371.0 / 360.0
SUBLAYERS = 3  # depth-kernel sublayers per node interval: DSurfTomo's manual default (its range 2-5)  # on DSurfTomo's sphere (grid and checkerboard cells)


def on_coarse_grid(periods, measured_per_octave, per_octave):
    """Which periods lie on the 1/per_octave-octave grid through 1 Hz, for periods measured on the
    1/measured_per_octave grid (every round(measured / per_octave)-th frequency; all when per_octave >= measured)."""
    step = max(1, int(round(measured_per_octave / per_octave)))
    n = np.rint(measured_per_octave * np.log2(1.0 / np.asarray(periods, float))).astype(int)
    return n % step == 0


def records_from_curves(curves, max_periods, per_octave=None):
    """per_octave: (measured, inversion) frequencies per octave; only the periods on the inversion's grid are kept
    (inversion.per_octave), curves left with none are dropped."""
    recs = []
    for c in curves:
        t, v = np.asarray(c["T"], float), np.asarray(c["v"], float)
        if per_octave is not None:
            keep = on_coarse_grid(t, *per_octave)
            t, v = t[keep], v[keep]
        if not len(t):
            continue
        recs.append(
            dict(
                sta1=(c["lat1"], c["lon1"]),
                sta2=(c["lat2"], c["lon2"]),
                values={round(float(a), 6): float(b) for a, b in zip(t, v)},
                name=c["name"],
            )
        )
    periods = sorted({t for r in recs for t in r["values"]})
    if len(periods) > max_periods:
        raise RuntimeError(
            "%d periods carry data (%.4g-%.4g s); the solver accepts at most %d (NP of DSurfTomo's CalSurfG): set "
            "inversion.per_octave to %d or fewer, or narrow export.tmin_s / export.tmax_s, then run the "
            "chain again from inversion"
            % (len(periods), periods[0], periods[-1], max_periods, _per_octave_for(periods, max_periods))
        )
    idx = {t: k + 1 for k, t in enumerate(periods)}
    for r in recs:
        r["index"] = {idx[t]: v for t, v in r["values"].items()}
    return recs, periods


def _per_octave_for(periods, cap):
    """The largest whole number of frequencies per octave that keeps the band of `periods` within `cap` periods."""
    octaves = np.log2(periods[-1] / periods[0])
    return max(1, int(np.floor((cap - 1) / max(octaves, 1e-9))))


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


def median_wavelength(recs):
    """Median of c T over every travel time of the project (km): the wavelength the data carry most."""
    return float(np.median([v * t for r in recs for t, v in r["values"].items()]))


def grid_spacing(recs, points):
    """inversion.grid_spacing_deg auto -> ([lat, lon] degrees, reason).  The cell is the median measured wavelength /
    3 in km in both directions (user's rule 2026-11-13; lateral resolution of surface-wave tomography is about a
    wavelength, three nodes sample it), on two significant digits."""
    lam = median_wavelength(recs)
    km = lam / 3.0
    lat0 = float(np.mean([p[0] for p in points]))
    dlat = float("%.2g" % (km / KM_PER_DEG))
    dlon = float("%.2g" % (km / (KM_PER_DEG * math.cos(math.radians(lat0)))))
    why = "median measured wavelength %.4g km / 3 = %.4g km -> %g deg (lat), %g deg (lon)" % (lam, km, dlat, dlon)
    return [dlat, dlon], why


def grid(points, spacing, refine=1, strip_cells=None):
    """points [(lat, lon)] -> dict(nx, ny, goxd, gozd, dvxd, dvzd) (x = latitude, z = longitude, as in DSurfTomo)."""
    lats = [p[0] for p in points]
    lons = [p[1] for p in points]
    lat_min, lat_max, lon_min, lon_max = min(lats), max(lats), min(lons), max(lons)
    dvx, dvz = spacing
    lat_rng, lon_rng = max(lat_max - lat_min, dvx), max(lon_max - lon_min, dvz)  # no fixed minimum extent
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
        if km_x >= km_z:  # line runs north-south: narrow in longitude
            n2 = int(math.ceil((lon_max - lon_min) / g["dvzd"] + 2 * strip_cells - 1e-9)) + 3
            g.update(ny=n2, gozd=0.5 * (lon_min + lon_max) - 0.5 * (n2 - 3) * g["dvzd"], strip="longitude")
        else:  # east-west: narrow in latitude
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
        fh.write(line("%.6g 0" % inv["weight"], "weight damp (lant: no damping)"))
        fh.write(line("%d" % SUBLAYERS, "sublayers (for computing depth kernel, 2~5)"))
        fh.write(
            line("%.3f %.3f" % (inv["vmin"], inv["vmax"]), "minimum velocity, maximum velocity (a priori information)")
        )
        fh.write(line("%d" % inv["iterations"], "maximum iteration"))
        fh.write(line("0.100", "sparsity fraction (not used by lant)"))
        fh.write(line("%d" % len(periods), "kmaxRc (followed by periods)"))
        fh.write(" ".join("%.6f" % t for t in periods) + "\n")
        fh.write(line("0", "kmaxRg"))
        fh.write(line("0", "kmaxLc"))
        fh.write(line("0", "kmaxLg"))
        fh.write(line("%d" % synthetic, "synthetic flag(0:real data,1:synthetic)"))
        fh.write(line("0.000", "noiselevel"))  # set for the checkerboards only (build_checks)
        fh.write(line("%.3f" % inv["outlier_mad"], "threshold (lant: k of the rejection beyond k x 1.4826 MAD)"))


def write_mod(path, z, vs, g, header=True, pattern=None, depth_factor=None, scale=None):
    """1-D model on every node; with pattern [ny, nx] (a fraction of the 1-D Vs) applied as vs (1 + f pattern), f the
    depth_factor of each depth (sign of the depth-varying checkerboards); scale [nz]: the nodes DSurfTomo updates
    (all depths but the last, the horizontal nodes inside the outer ring) multiplied by scale[k] (control
    inversions); the fixed nodes keep the 1-D model, as in the main project."""
    nx, ny = g["nx"], g["ny"]
    with open(path, "w") as f:
        if header:
            f.write(" ".join("%7.3f" % d for d in z) + "\n")
        for k in range(len(z)):
            s = 1.0 if depth_factor is None else depth_factor[k]
            for j in range(ny):
                row = []
                for i in range(nx):
                    v = vs[k]
                    if scale is not None and k < len(z) - 1 and 0 < j < ny - 1 and 0 < i < nx - 1:
                        v = v * scale[k]
                    if pattern is not None:
                        v = v * (1 + s * pattern[j, i])
                    row.append("%.5f" % v)  # 5 decimals: 0.1 km/s sediments keep 1e-4 relative precision
                f.write(" ".join(row) + "\n")


def control_scale(name, z):
    """Factors of the updatable depths of a starting-model control: start_x<f> multiplies them all by f (the
    residuals then all have one sign: a test of how far the data pull the model back); start_grad<f> multiplies them by 1 + f at the surface changing linearly to 1 - f at the
    deepest updatable node (a smooth change of the vertical gradient, on the scale of the depth kernels, so the data
    see it); the build step then centres it on the data (centre_on_data)."""
    n = len(z)
    if name.startswith("start_x"):
        f = float(name[7:])
        if not f > 0:
            raise RuntimeError("checks.controls: %s needs a positive factor" % name)
        return np.r_[np.full(n - 1, f), 1.0]
    if name.startswith("start_grad"):
        f = float(name[10:])
        if not 0 < abs(f) < 0.5:
            raise RuntimeError("checks.controls: %s needs 0 < |f| < 0.5" % name)
        zu = np.asarray(z[: n - 1], float)
        return np.r_[1 + f * (1 - 2 * zu / max(zu[-1], 1e-9)), 1.0]
    raise ValueError(name)


def centre_on_data(z, vs, scale, data, iters=6):
    """start_grad controls: one more factor u on the updatable depths, so that the median travel-time residual of
    the data (laterally uniform model, great-circle paths: t = distance / c_1D(T)) equals that of the main starting
    model.  The alternating perturbation is zero-mean in depth but not in the data (the shallow nodes govern the
    high frequencies; long paths weigh more in seconds); centred, the control differs from the main project in the
    shape of the model, not in an overall offset the data would first have to remove.  data: rows (lat1, lon1, lat2, lon2, period, velocity) of surfdata.dat.
    -> (scale x u on the updatable depths, info)."""
    from ..signal import sphere_km
    from .start import forward

    dist = sphere_km(data[:, 1], data[:, 0], data[:, 3], data[:, 2])  # as DSurfTomo turns velocities into times
    T = np.unique(data[:, 4])
    k = np.searchsorted(T, data[:, 4])
    t_obs = dist / data[:, 5]
    med = lambda v: float(np.median(t_obs - dist / forward(z, v, T)[k]))
    base = med(vs)
    upd = np.r_[np.ones(len(z) - 1, bool), False]
    lo, hi = 0.8, 1.25  # bisection on u: the median residual grows with u
    for _ in range(30):
        u = 0.5 * (lo + hi)
        if med(vs * np.where(upd, scale * u, scale)) > base:
            hi = u
        else:
            lo = u
    u = 0.5 * (lo + hi)
    sc = np.where(upd, scale * u, scale)
    return sc, dict(
        uniform_factor=round(u, 5), median_residual_main_ms=1e3 * base, median_residual_control_ms=1e3 * med(vs * sc)
    )


def checker_pattern(g, cell_km, amp, array_type, line_dir=None):
    """[ny, nx] pattern.  Areal arrays: sin(a_lat * i) sin(a_lon * j) on the node index of the coarse grid (cell size
    cell_km = [north-south, east-west] km; 0 = constant in that direction).  Lines: sin(pi (s - s_mid) / cell) of the
    distance s along the line from the middle of the grid (cell_km[0]), constant across it unless cell_km[1] > 0."""
    nx, ny = g["nx"], g["ny"]
    r = g["refine"]
    xi = (np.arange(1, nx + 1) - 1) / r + 1  # index on the coarse grid
    yj = (np.arange(1, ny + 1) - 1) / r + 1
    latc = g["goxd"] - (nx - 1) / 2.0 * g["dvxd"]
    if array_type == "line":
        lat, lon = node_coords(g)
        ys = (g["goxd"] - lat) * KM_PER_DEG  # km south of the grid origin
        xs = (lon - g["gozd"]) * KM_PER_DEG * math.cos(math.radians(latc))  # km east
        u = np.asarray(line_dir, float) / np.hypot(*line_dir)  # (south, east) unit vector along the line
        s = u[0] * ys[None, :] + u[1] * xs[:, None]  # [ny, nx]
        t = -u[1] * ys[None, :] + u[0] * xs[:, None]
        # odd about the middle of the grid along the line: the pattern averages to zero along the line, so the
        # residuals of the synthetic data do not lean to one side (a line holds only a few cells, and a pattern
        # starting at the grid edge can be mostly + or -)
        sc = 0.5 * (s.min() + s.max())
        p = np.sin(np.pi * (s - sc) / cell_km[0]) if cell_km[0] else np.ones_like(s)
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
    check_layers(len(z), SUBLAYERS)
    from .start import PROJECT_MAX_PERIODS

    recs, periods = records_from_curves(
        curves, PROJECT_MAX_PERIODS, (cfg["dispersion"]["per_octave"], inv["per_octave"])
    )
    os.makedirs(out_dir, exist_ok=True)
    info = write_surfdata(recs, os.path.join(out_dir, "surfdata.dat"))
    pts = info.pop("stations")
    spacing_reason = "inversion.grid_spacing_deg"
    if C.is_auto(inv["grid_spacing_deg"]):
        spacing, spacing_reason = grid_spacing(recs, pts)
        inv = dict(inv, grid_spacing_deg=spacing)
    g = grid(
        pts,
        inv["grid_spacing_deg"],
        strip_cells=inv["line_strip_cells"] if cfg["project"]["array_type"] == "line" else None,
    )
    write_in(os.path.join(out_dir, "DSurfTomo.in"), g, len(z), info["nsrc"], inv, periods)
    write_mod(os.path.join(out_dir, "MOD"), z, vs, g)
    np.savetxt(os.path.join(out_dir, "vs1d.txt"), np.c_[z, vs], fmt="%10.4f %10.4f")
    meta = dict(
        grid=dict(g, spacing_reason=spacing_reason),
        periods=periods,
        per_octave=inv["per_octave"],
        curves=len(recs),
        travel_times=info["travel_times"],
        nsrc=info["nsrc"],
        depths=list(map(float, z)),
        weight=inv["weight"],
        line_direction=list(map(float, line_direction(pts))) if cfg["project"]["array_type"] == "line" else None,
        stations_with_data=len(pts),
    )
    json.dump(meta, open(os.path.join(out_dir, "project.json"), "w"), indent=1)
    return meta


def copy_project(src, dst, files=("DSurfTomo.in", "MOD", "surfdata.dat", "vs1d.txt", "project.json")):
    """Copy the main project into a new directory dst (replaced if it exists)."""
    shutil.rmtree(dst, ignore_errors=True)
    os.makedirs(dst)
    for f in files:
        shutil.copyfile(os.path.join(src, f), os.path.join(dst, f))
    return True


def settle(new, dst):
    """Put the freshly built project `new` in place of dst.  A dst that has already run is kept when its inputs are
    the same, and is an error when they differ (its results belong to other settings: remove it to rebuild)."""
    import filecmp

    if os.path.exists(os.path.join(dst, "DSurfTomo.in.log")):
        names = sorted(os.listdir(new))
        same = all(
            os.path.exists(os.path.join(dst, f))
            and filecmp.cmp(os.path.join(new, f), os.path.join(dst, f), shallow=False)
            for f in names
        )
        shutil.rmtree(new)
        if same:
            return "kept (already run with the same inputs)"
        raise RuntimeError(
            "%s has run with other inputs than the present configuration gives: move it away (or delete "
            "it) to rebuild" % dst
        )
    shutil.rmtree(dst, ignore_errors=True)
    os.replace(new, dst)
    return "built"


def set_line(path, key, value):
    L = open(path).read().split("\n")
    k = next(i for i, l in enumerate(L) if "c: " + key in l)
    L[k] = "%-33s" % value + L[k][L[k].index("c: ") :]
    open(path, "w").write("\n".join(L))


def _number(name, text):
    try:
        x = float(text)
    except ValueError:
        x = float("nan")
    if not x > 0:
        raise RuntimeError("checks.controls: %s needs a positive number, not %r" % (name, text))
    return x


def set_regularization(d, weight, control):
    """Smoothing weight of a control project, in DSurfTomo.in and in its project.json (which was copied from the main
    project, so the weight there is replaced)."""
    set_line(os.path.join(d, "DSurfTomo.in"), "weight damp", "%.6g 0" % weight)
    path = os.path.join(d, "project.json")
    meta = json.load(open(path))
    meta.update(weight=weight, control=control)
    json.dump(meta, open(path, "w"), indent=1)


def build_checks(cfg, main_dir, out_root, noise=None, part="all"):
    """Checkerboards, depth-varying checkerboards and the control inversions, built from the main project.
    noise: a fraction of each synthetic travel time (DSurfTomo's noiselevel), or [[period_s, sigma_s], ...] written
    to dsurf_main.NOISE_FILE (lant's noise per period, what inversion.noise auto gives); None leaves the
    checkerboards to be built later (listed under 'deferred'), when inversion.noise is auto and the main project has
    to run first.  part: all | checkerboards."""
    meta = json.load(open(os.path.join(main_dir, "project.json")))
    g = meta["grid"]
    m = np.loadtxt(os.path.join(main_dir, "vs1d.txt"))
    z, vs = m[:, 0], m[:, 1]
    ck = cfg["checks"]
    at = cfg["project"]["array_type"]
    built = []
    info = {}
    horizontal, vertical, design = ck["checkerboards"], ck["vertical"], None
    if horizontal == "auto" or vertical == "auto":
        from .checker_design import design as rule

        h_auto, v_auto, design = rule(main_dir, z, at)
        horizontal = h_auto if horizontal == "auto" else horizontal
        vertical = v_auto if vertical == "auto" else vertical
    specs = [(c["name"], c["cell_km"], None) for c in horizontal] + [
        (c["name"], c["cell_km"], c["reversal_depths_km"]) for c in vertical
    ]
    if len({n for n, _, _ in specs}) != len(specs):
        raise RuntimeError("checkerboard names must be different: %s" % [n for n, _, _ in specs])
    deferred = []
    for name, cell, rev in specs:
        if noise is None:
            deferred.append(name)
            continue
        d = os.path.join(out_root, "checkerboards", name)
        new = d + ".new"
        copy_project(main_dir, new)
        pat = checker_pattern(g, cell, ck["amplitude"], at, meta.get("line_direction"))
        fac = depth_signs(z, cfg["inversion"]["depth_sequence_km"], rev) if rev else None
        write_mod(os.path.join(new, "MOD.true"), z, vs, g, header=False, pattern=pat, depth_factor=fac)
        set_line(os.path.join(new, "DSurfTomo.in"), "synthetic flag", "1")
        if isinstance(noise, (list, tuple)):
            from .dsurf_main import NOISE_FILE

            np.savetxt(os.path.join(new, NOISE_FILE), np.asarray(noise, float), fmt="%.6f",
                       header="period_s sigma_s: Gaussian noise of the synthetic travel times")  # fmt: skip
            set_line(os.path.join(new, "DSurfTomo.in"), "noiselevel", "0.0000")
        else:
            set_line(os.path.join(new, "DSurfTomo.in"), "noiselevel", "%.4f" % noise)
        built.append((name, settle(new, d)))
    for name in ck["controls"] if part == "all" else []:
        d = os.path.join(out_root, "controls", name)
        new = d + ".new"
        copy_project(main_dir, new)
        dd, d = d, new
        if name.startswith("reg_x"):  # the main project's smoothing weight times the factor
            factor = _number(name, name[5:])
            set_regularization(d, float("%.6g" % (factor * meta["weight"])), name)
        elif name.startswith("reg_w"):
            set_regularization(d, round(_number(name, name[5:]), 3), name)
        elif name.startswith(("start_x", "start_grad")):  # the updatable nodes only: fixed nodes as in main
            sc = control_scale(name, z)
            if name.startswith("start_grad"):
                from .run import read_data, read_in

                data = read_data(
                    os.path.join(main_dir, "surfdata.dat"), read_in(os.path.join(main_dir, "DSurfTomo.in"))[4]
                )
                sc, info[name] = centre_on_data(z, vs, sc, data)
            write_mod(os.path.join(d, "MOD"), z, vs, g, scale=sc)
        else:
            raise RuntimeError(
                "checks.controls: unknown control %s (use reg_x<factor>, reg_w<weight>, start_grad<f> or "
                "start_x<factor>)" % name
            )
        built.append((name, settle(new, dd)))
    return dict(
        built=built,
        deferred=deferred,
        horizontal=horizontal,
        vertical=vertical,
        design=design,
        controls=info,
        noise=noise,
    )
