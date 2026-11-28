"""Phase velocities of the paths of part A through a 3-D model (part B of the synthetic example).

    python make_3d_data.py A_CONFIG B_DIR [--model structures | spikes] [--cores 8]

A_CONFIG: the configuration of part A after `lant disp`; its accepted phase velocities give the station pairs and
periods.  The true model: the layered model of models.txt with a sedimentary basin, a dipping low-velocity fault
zone and a buried high-velocity body (structures), or with five isolated anomalies at three depths (spikes).  The travel times of every pair and period are computed by
ray tracing in the phase-velocity maps of the true model (lant's port of DSurfTomo's forward calculation: fast
marching and surfdisp96 on every column) on a grid twice as fine as the one lant inv derives, then each gets
Gaussian noise with the jackknife error of its measurement in part A (sigma_t = dist sigma_c / c^2), as in the
synthetic tests of Fang et al. (2015) and Rawlinson & Spakman (2016).

Writes B_DIR/config.yaml (part A's settings, output in B_DIR/output), B_DIR/output/3_disp/dispersion.csv (the
input of lant inv), B_DIR/output/1_clean/stations_used.csv (for the maps) and B_DIR/true_model.csv.  Then:
    lant inv -c B_DIR/config.yaml
    lant checker -c B_DIR/config.yaml --size 1
"""

import argparse
import os
import shutil
import sys

import modes as M
import numpy as np
import yaml


def background(z):
    """Vs of the layered model at depths z (km)."""
    lay = M.read_models()["layered"]
    tops = np.r_[0, np.cumsum(lay[:-1, 0])]
    return lay[np.searchsorted(tops, z, side="right") - 1, 2]


def _smooth_step(d, w):
    """1 inside (d < -w/2), 0 outside (d > w/2), a cosine in between."""
    return 0.5 * (1 - np.sin(np.pi * np.clip(d / w, -0.5, 0.5)))


def anomaly(x, y, z):
    """Relative Vs perturbation at local coordinates x (east), y (north), z (depth), km from the array centre; all
    within the depth range the data resolve (about the upper 0.6 km: half the longest wavelength)."""
    a = np.zeros(np.broadcast(x, y, z).shape)
    # sedimentary basin: a bowl 1.5 km in radius and 0.12 km deep at its centre, Vs 30 % lower inside
    r = np.hypot(x + 1.0, y + 0.9) / 1.5
    base = 0.12 * (1 - r**2)
    a -= 0.30 * _smooth_step(z - base, 0.02) * (r < 1)
    # fault zone: 0.3 km wide, striking N30E through (1.0, 0), dipping 60 degrees to the south-east, Vs 15 % lower,
    # down to 0.6 km
    s = np.radians(30.0)
    n = np.array([np.cos(s), -np.sin(s)])  # horizontal normal, pointing south-east
    dist = (x - 1.0) * n[0] + y * n[1] - z / np.tan(np.radians(60.0))  # horizontal offset from the dipping plane
    a -= 0.15 * _smooth_step(np.abs(dist) - 0.15, 0.08) * _smooth_step(z - 0.6, 0.05)
    # buried high-velocity body: an ellipsoid at 0.35 km depth (1.8 x 1.4 x 0.24 km), Vs 10 % higher
    e = np.sqrt(((x - 1.3) / 0.9) ** 2 + ((y - 1.4) / 0.7) ** 2 + ((z - 0.35) / 0.12) ** 2)
    a += 0.10 * _smooth_step(e - 1.0, 0.3)
    return a


SPIKES = [
    (-1.5, 1.2, 0.05, -0.1),
    (1.6, -1.4, 0.05, 0.1),
    (1.4, 1.3, 0.2, -0.1),
    (-1.4, -1.3, 0.2, 0.1),
    (0.0, 0.0, 0.4, 0.1),
]  # x, y, depth (km), relative perturbation


def spikes(x, y, z):
    """Sparse spike test (Rawlinson & Spakman 2016): isolated Gaussian anomalies of +-10 %, 0.4 km horizontal and
    0.03-0.1 km vertical half-width, at three depths."""
    a = np.zeros(np.broadcast(x, y, z).shape)
    for x0, y0, z0, amp in SPIKES:
        a += amp * np.exp(-(((x - x0) ** 2 + (y - y0) ** 2) / 0.4**2 + ((z - z0) / max(0.03, 0.25 * z0)) ** 2))
    return a


def true_model(lon, lat, z, lon0, lat0, kind="structures"):
    kx, ky = 111.195 * np.cos(np.radians(lat0)), 111.195
    x, y = (lon - lon0) * kx, (lat - lat0) * ky
    return background(z) * (1 + (anomaly if kind == "structures" else spikes)(x, y, z))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("config")
    ap.add_argument("out")
    ap.add_argument("--cores", type=int, default=os.cpu_count())
    ap.add_argument("--seed", type=int, default=3)
    ap.add_argument("--model", choices=("structures", "spikes"), default="structures")
    a = ap.parse_args()
    sys.path.insert(0, os.path.join(M.HERE, "..", ".."))
    from lant import config as C
    from lant.inversion import dsurf_input, dsurf_main
    from lant.inversion import project as P
    from lant.inversion.run import read_dispersion

    cfg = C.load(a.config)
    tab = read_dispersion(os.path.join(cfg["output_dir"], "3_disp", "dispersion.csv"))
    tab = tab[np.isfinite(tab["c_kms"])]
    po = int(cfg["inv"]["per_octave"])
    keep = P.on_grid(tab["period_s"], po)
    tab = tab[keep]
    recs, periods = P.records(tab, po, "phase")
    sta = np.genfromtxt(os.path.join(cfg["output_dir"], "1_clean", "stations_used.csv"), delimiter=",", names=True,
                        dtype=None, encoding="utf-8")  # fmt: skip
    lon0, lat0 = float(np.mean(sta["longitude"])), float(np.mean(sta["latitude"]))
    lam = float(np.median([v * t for r in recs for (k, t), v in r["values"].items()]))
    spacing = float("%.2g" % (lam / 3)) / 2  # half the grid lant inv derives (lant.inversion.run)
    zmax = 2.0
    z = np.round(np.r_[np.arange(0, 0.8, 0.02), np.arange(0.8, zmax + 1e-9, 0.1)], 4)

    work = os.path.join(a.out, "forward")
    if os.path.isdir(work):
        shutil.rmtree(work)
    os.makedirs(work)
    info = P.write_surfdata(recs, os.path.join(work, "surfdata.dat"))
    g = P.grid(info["stations"], spacing)
    inv = dict(cfg["inv"], vs_min=0.1, vs_max=10.0, weight=0.0, damping=0.0, iterations=0,
               outlier_mad=float(cfg["inv"]["outlier_mad"]))  # fmt: skip
    P.write_control(os.path.join(work, "DSurfTomo.in"), g, len(z), info["nsrc"], inv, periods, synthetic=1, noise=0.0)
    P.write_model(os.path.join(work, "MOD"), z, background(z), g)
    nlat, nlon = P.node_coords(g)
    LA, LO = np.meshgrid(nlat, nlon, indexing="xy")  # [ny (lon), nx (lat)]
    true = np.stack([true_model(LO, LA, zk, lon0, lat0, a.model) for zk in z])  # [nz, ny, nx]
    pr = dsurf_input.read_project(work)
    print("forward grid: %d x %d nodes every %.3g km, %d depth nodes; %d travel times, %d periods"
          % (g["nx"], g["ny"], spacing, len(z), len(pr.measurements.observed), len(periods["phase"])))  # fmt: skip
    obst = dsurf_main.synthetic_data(pr, true.astype(np.float32), a.cores, None)
    dist = pr.measurements.distance.astype(float)
    # rows of surfdata.dat in file order -> (pair, period) -> jackknife error of part A
    rows = []
    src = None
    for ln in open(os.path.join(work, "surfdata.dat")):
        f = ln.split()
        if ln.startswith("#"):
            src, k = (float(f[1]), float(f[2])), int(f[3])
        else:
            rows.append((src, (float(f[0]), float(f[1])), periods["phase"][k - 1]))
    sig = {}
    for r in tab:
        key = (round(r["lat1"], 6), round(r["lon1"], 6), round(r["lat2"], 6), round(r["lon2"], 6),
               round(float(r["period_s"]), 6))  # fmt: skip
        sig[key] = (float(r["sigma_kms"]), float(r["c_kms"]))
    rng = np.random.default_rng(a.seed)
    out_rows = []
    for (s, rcv, t), d, ob in zip(rows, dist, obst):
        key = (round(s[0], 6), round(s[1], 6), round(rcv[0], 6), round(rcv[1], 6), round(t, 6))
        sc, c0 = sig.get(key, sig.get((key[2], key[3], key[0], key[1], key[4]), (0.0, 1.0)))
        st = d * sc / c0**2
        tt = float(ob) + rng.normal() * st
        out_rows.append((s, rcv, d, t, d / tt, sc))
    os.makedirs(os.path.join(a.out, "output", "3_disp"), exist_ok=True)
    os.makedirs(os.path.join(a.out, "output", "1_clean"), exist_ok=True)
    shutil.copy(
        os.path.join(cfg["output_dir"], "1_clean", "stations_used.csv"), os.path.join(a.out, "output", "1_clean")
    )
    with open(os.path.join(a.out, "output", "3_disp", "dispersion.csv"), "w") as fh:
        fh.write("sta1,sta2,lon1,lat1,lon2,lat2,dist_km,frequency_hz,period_s,c_kms,sigma_kms,u_kms,u_sigma_kms,snr\n")
        for s, rcv, d, t, c, sc in out_rows:
            fh.write("-,-,%.6f,%.6f,%.6f,%.6f,%.4f,%.6f,%.6f,%.5f,%.5f,,,\n" % (s[1], s[0], rcv[1], rcv[0], d, 1 / t, t,
                                                                              c, sc))  # fmt: skip
    with open(os.path.join(a.out, "true_model.csv"), "w") as fh:
        fh.write("lon,lat,depth_km,vs_kms\n")
        for k, zk in enumerate(z):
            for j in range(len(nlon)):
                for i in range(len(nlat)):
                    fh.write("%.6f,%.6f,%.4f,%.5f\n" % (nlon[j], nlat[i], zk, true[k, j, i]))
    c2 = dict(cfg)
    c2.pop("_file", None)
    c2["output_dir"] = os.path.abspath(os.path.join(a.out, "output"))
    with open(os.path.join(a.out, "config.yaml"), "w") as fh:
        fh.write("# part B of the synthetic example: the settings of part A, output in this folder\n")
        yaml.safe_dump(c2, fh, sort_keys=False)
    print("written: %d phase velocities -> %s" % (len(out_rows), a.out))


if __name__ == "__main__":
    main()
