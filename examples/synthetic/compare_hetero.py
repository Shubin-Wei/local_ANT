#!/usr/bin/env python
"""Compare a run on a make_stack.py stack (laterally varying medium) with its truth file.

    python examples/synthetic/compare_hetero.py --output OUT --truth OUT/truth.npz

Prints numbers only:
  * phase velocities: accepted measurements against the path-average velocity of the truth, per frequency band
    (median and 95th percentile of |measured / true - 1|), and the correlation between the measured and the true
    path anomaly (velocity of a pair relative to the uniform background);
  * model: per depth node, at grid nodes within one cell of a station, correlation and regression slope of the
    recovered lateral Vs anomaly (relative to the mean of the depth) against the pattern put in.
"""

import argparse
import os
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))


def main():
    from make_stack import pattern

    from lant.inversion.run import read_model
    from lant.signal import local_xy_km

    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output", required=True)
    p.add_argument("--truth", required=True)
    p.add_argument("--model", default=None, help="another model file (default: OUT/06_inversion/main)")
    a = p.parse_args()
    T = np.load(a.truth)
    amp, cell = float(T["hetero"]), float(T["cell"])
    with h5py.File(os.path.join(a.output, "work", "05_dispersion", "disp.h5"), "r") as h:
        f, vc, ok = h["freqs"][:], h["vc"][:], h["ok_c"][:]
        pairs = h["pairs"][:].astype("U")
    order = {q: i for i, q in enumerate(T["pairs"].astype("U"))}
    idx = np.array([order[q] for q in pairs])
    lf = np.log(T["freqs"])
    true = np.array([np.interp(np.log(f), lf, T["c_path"][i]) for i in idx])
    c0 = np.interp(np.log(f), lf, T["c0"])
    rel = np.where(ok, vc / true - 1, np.nan)
    print("phase velocities against the truth (accepted samples)")
    print("  band (Hz)        samples   median(m/t-1) %   p95|m/t-1| %")
    edges = np.geomspace(f[0], f[-1] * 1.0001, 5)
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (f >= lo) & (f < hi)
        r = rel[:, m][ok[:, m]]
        if r.size:
            print(
                "  %5.2f - %5.2f   %7d   %+9.3f          %7.3f"
                % (lo, hi, r.size, 100 * np.median(r), 100 * np.percentile(np.abs(r), 95))
            )
    r = rel[ok]
    print(
        "  all             %7d   %+9.3f          %7.3f"
        % (r.size, 100 * np.median(r), 100 * np.percentile(np.abs(r), 95))
    )
    am, at = (vc / c0[None, :] - 1)[ok], (true / c0[None, :] - 1)[ok]
    print(
        "  path anomaly (relative to the uniform background): true sd %.3f %%, measured sd %.3f %%, correlation %.3f, "
        "slope %.3f"
        % (
            100 * at.std(),
            100 * am.std(),
            np.corrcoef(am, at)[0, 1],
            np.dot(am - am.mean(), at - at.mean()) / np.dot(at - at.mean(), at - at.mean()),
        )
    )
    lon, lat, dep, v = read_model(
        a.model or os.path.join(a.output, "work", "06_inversion", "main", "DSurfTomo.inMeasure.dat")
    )
    st = np.genfromtxt(
        os.path.join(a.output, "work", "01_clean", "stations.csv"), delimiter=",", names=True, dtype=None, encoding=None
    )
    G = np.meshgrid(lon, lat, indexing="ij")
    x, y = local_xy_km(G[0].ravel(), G[1].ravel(), float(T["lat0"]), float(T["lon0"]))
    sx, sy = local_xy_km(st["longitude"], st["latitude"], float(T["lat0"]), float(T["lon0"]))
    near = np.hypot(x[:, None] - sx[None, :], y[:, None] - sy[None, :]).min(axis=1) <= cell
    d = pattern(x, y, cell, amp)[near]
    print(
        "model against the pattern (amplitude %.1f %%, cell %.1f km; %d grid nodes within one cell of a station)"
        % (100 * amp, cell, near.sum())
    )
    print("  depth (km)   correlation   slope   recovered sd %")
    for k, z in enumerate(dep):
        vz = v[:, :, k].ravel()[near]
        an = vz / vz.mean() - 1
        if an.std() < 1e-6:
            continue
        print(
            "  %8.2f     %+7.3f     %+6.3f   %7.3f"
            % (z, np.corrcoef(an, d)[0, 1], np.dot(an, d) / np.dot(d, d), 100 * an.std())
        )


if __name__ == "__main__":
    main()
