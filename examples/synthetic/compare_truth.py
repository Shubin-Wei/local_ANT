#!/usr/bin/env python
"""Compare the results of the synthetic example with the known answer.

    python examples/synthetic/compare_truth.py      (after `lant run all -c synthetic.yaml` in examples/synthetic)
    python examples/synthetic/compare_truth.py --dir DIR [--output OUT]   (another example made with --out DIR)

Writes output/results/figures/truth_comparison.png and prints:
  * median measured phase velocity against the 1-D background (not the path-dependent truth with bodies);
  * the starting model and horizontal means of the inverted and true models, plus 3-D RMS error per depth.
true_checkerboard.json supplies the checkerboard truth (zero amplitude with --uniform).
Model summaries include all output nodes, including poorly constrained edges; they are not resolution estimates.
"""

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))


def _maps(fig, axs, out, description, lon, lat, dep, recovered, truth, tm):
    import h5py

    from lant import data
    from lant.plotting import style
    from lant.signal import local_xy_km

    kx, ky = description["km_per_degree"]
    centre_lon = description["origin_lon"] + description["size_km"] / (2 * kx)
    centre_lat = description["origin_lat"] + description["size_km"] / (2 * ky)
    stations = data.read_canonical(data.canonical_table(out))
    sx, sy = local_xy_km(
        [v["lon"] for v in stations.values()], [v["lat"] for v in stations.values()], centre_lat, centre_lon
    )
    gx, _ = local_xy_km(lon, np.full(len(lon), centre_lat), centre_lat, centre_lon)
    _, gy = local_xy_km(np.full(len(lat), centre_lon), lat, centre_lat, centre_lon)
    X, Y = np.meshgrid(gx, gy, indexing="ij")
    k = int(np.argmin(abs(dep - min(0.6, dep[-1] / 2))))
    base = np.interp(dep[k], tm[:, 0], tm[:, 1])
    limit = 100 * (description["amplitude"] or 0.06)
    for ax, field, title in (
        (axs[0], truth, "(a) Input checkerboard"),
        (axs[1], recovered, "(b) Recovered, %.2f km" % dep[k]),
    ):
        image = ax.pcolormesh(
            X, Y, 100 * (field[:, :, k] / base - 1), shading="gouraud", cmap=style.SEIS, vmin=-limit, vmax=limit
        )
        ax.plot(sx, sy, "^", ms=6, mfc="white", mec="black", mew=0.8, zorder=5)
        ax.set_aspect("equal")
        margin = 0.06 * description["size_km"]
        ax.set_xlim(sx.min() - margin, sx.max() + margin)
        ax.set_ylim(sy.min() - margin, sy.max() + margin)
        ax.set_xlabel("East (km)")
        ax.set_title(title)
    axs[0].set_ylabel("North (km)")
    with h5py.File(os.path.join(out, "work/03_stack/stack.h5"), "r") as h:
        px, py = local_xy_km(
            np.c_[h["lon1"][:], h["lon2"][:]], np.c_[h["lat1"][:], h["lat2"][:]], centre_lat, centre_lon
        )
    for x, y in zip(px, py):
        axs[0].plot(x, y, color="0.15", lw=0.25, alpha=0.25, zorder=2)
    cax = fig.add_axes([0.90, 0.56, 0.018, 0.32])
    fig.colorbar(image, cax=cax, label="Vs anomaly (%)")


def main():
    import argparse

    from lant.export import read_curves
    from lant.inversion.run import read_model

    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=HERE, help="the example directory (true_model.txt, true_dispersion.txt)")
    ap.add_argument("--output", default=None, help="the output directory (default DIR/output)")
    a = ap.parse_args()
    out = a.output or os.path.join(a.dir, "output")
    td = np.loadtxt(os.path.join(a.dir, "true_dispersion.txt"))
    tm = np.loadtxt(os.path.join(a.dir, "true_model.txt"))
    curves = read_curves(os.path.join(out, "work", "05_dispersion", "curves.csv"), "C")
    acc = {}
    for c in curves:
        for t, v in zip(c["T"], c["v"]):
            acc.setdefault(round(float(t), 4), []).append(v)
    T = np.array(sorted(acc))
    f = 1 / T
    med = np.array([np.median(acc[t]) for t in T])
    n = np.array([len(acc[t]) for t in T])
    true_c = np.interp(np.log(f), np.log(td[:, 0]), td[:, 1])
    print("phase velocity: median of %d curves against the 1-D background (not path truth with bodies)" % len(curves))
    print("   f (Hz)   pairs  median  background   difference")
    for k in range(0, len(f), max(1, len(f) // 12)):
        print("  %6.2f   %5d   %.3f  %.3f  %+5.1f %%" % (f[k], n[k], med[k], true_c[k], 100 * (med[k] / true_c[k] - 1)))
    rel = 100 * (med / true_c - 1)
    print(
        "  all frequencies: median difference %+.2f %%, max |difference| %.2f %%" % (np.median(rel), np.abs(rel).max())
    )
    model = os.path.join(out, "work", "06_inversion", "main", "DSurfTomo.inMeasure.dat")
    if not os.path.exists(model):
        print("\nno inverted model yet (%s): run the chain to the end (lant run all) for the Vs comparison" % model)
        return
    vs1d = np.loadtxt(os.path.join(out, "work", "06_inversion", "vs1d.txt"))
    lon, lat, dep, v = read_model(os.path.join(out, "work", "06_inversion", "main", "DSurfTomo.inMeasure.dat"))
    mean = np.nanmean(v, axis=(0, 1))
    std = np.nanstd(v, axis=(0, 1))
    checker_path = os.path.join(a.dir, "true_checkerboard.json")
    if not os.path.exists(checker_path):
        raise SystemExit("Missing checkerboard truth; generate a new example with python -m lant.synthetic")
    from lant import synthetic_checkerboard

    description = json.load(open(checker_path))
    truth = synthetic_checkerboard.model(description, lon, lat, dep, tm[:, 0], tm[:, 1])
    true_mean = truth.mean(axis=(0, 1))
    rms = np.sqrt(np.mean((v - truth) ** 2, axis=(0, 1)))
    print("\nVs (km/s) per depth node, all grid nodes: start, inverted (mean +- std), truth mean, 3-D RMS error")
    for k, z in enumerate(dep):
        print(
            "  %5.2f km   %.3f   %.3f +- %.3f   %.3f   %.4f"
            % (z, np.interp(z, vs1d[:, 0], vs1d[:, 1]), mean[k], std[k], true_mean[k], rms[k])
        )
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from lant.plotting import style

    style.apply()
    fig, panels = plt.subplots(2, 2, figsize=(10.4, 8.4))
    fig.subplots_adjust(wspace=0.30, hspace=0.34, bottom=0.09, right=0.88)
    _maps(fig, panels[0], out, description, lon, lat, dep, v, truth, tm)
    axs = panels[1]
    for c in curves:
        axs[0].semilogx(1 / c["T"], c["v"], "-", c="0.75", lw=0.5)
    axs[0].semilogx(f, med, "ko", ms=3, label="median of the measured curves")
    axs[0].semilogx(td[:, 0], td[:, 1], color=style.RED, lw=1.8, label="1-D background")
    axs[0].set_title("(c) Phase velocities", loc="left")
    axs[0].set_xlabel("Frequency (Hz)")
    axs[0].set_ylabel("Phase velocity (km/s)")
    axs[0].set_xlim(f.min() * 0.9, f.max() * 1.1)
    from matplotlib.ticker import FuncFormatter, NullFormatter

    axs[0].xaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    axs[0].xaxis.set_minor_formatter(NullFormatter())
    axs[0].legend(frameon=False)
    axs[1].plot(tm[:, 1], tm[:, 0], "--", color=style.RED, lw=1.5, label="True 1-D background")
    axs[1].plot(vs1d[:, 1], vs1d[:, 0], color="0.45", lw=1.2, label="1-D start")
    axs[1].fill_betweenx(dep, mean - std, mean + std, color=style.BLUE, alpha=0.2, lw=0)
    axs[1].plot(mean, dep, color=style.BLUE, lw=1.8, label="3-D mean ± std.")
    axs[1].invert_yaxis()
    axs[1].set_xlabel("Vs (km/s)")
    axs[1].set_ylabel("Depth (km)")
    axs[1].legend(frameon=False, loc="lower left")
    axs[1].set_title("(d) Vs profiles", loc="left")
    fig.text(
        0.06,
        0.02,
        "Straight-ray synthetic noise; bent-ray inversion. Shared colour limits; the map does not show resolution.",
        fontsize=9,
    )
    p = os.path.join(out, "results", "figures", "truth_comparison.png")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    from lant import artifacts, records
    from lant.provenance import Lock

    cfg = {"project": {"output_dir": os.path.abspath(out)}}
    with Lock(records.rec_dir(cfg)):
        before = artifacts.snapshot(cfg)
        fig.savefig(p, dpi=150)
        plt.close(fig)
        artifacts.register(cfg, before)
    print("\nwrote", p)


if __name__ == "__main__":
    main()
