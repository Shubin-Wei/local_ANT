#!/usr/bin/env python
"""Compare the results of the synthetic example with the known answer.

    python examples/synthetic/compare_truth.py      (after `lant run all -c synthetic.yaml` in examples/synthetic)
    python examples/synthetic/compare_truth.py --dir DIR [--output OUT]   (another example made with --out DIR)

Writes output/figures/truth_comparison.png and prints:
  * the median of the measured phase velocities against the true dispersion curve, per frequency;
  * the 1-D starting model and the horizontal mean of the inverted model against the true model, per depth node.
The synthetic medium is laterally uniform, so the inverted model should be close to the true one wherever the
data constrain it (down to about half the longest wavelength), and laterally almost uniform.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))


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
    print("phase velocity: median of %d curves against the true curve" % len(curves))
    print("   f (Hz)   pairs  median  true   difference")
    for k in range(0, len(f), max(1, len(f) // 12)):
        print("  %6.2f   %5d   %.3f  %.3f  %+5.1f %%" % (f[k], n[k], med[k], true_c[k], 100 * (med[k] / true_c[k] - 1)))
    rel = 100 * (med / true_c - 1)
    print("  all frequencies: median difference %+.2f %%, max |difference| %.2f %%" % (np.median(rel), np.abs(rel).max()))
    model = os.path.join(out, "work", "06_inversion", "main", "DSurfTomo.inMeasure.dat")
    if not os.path.exists(model):
        print("\nno inverted model yet (%s): run the chain to the end (lant run all) for the Vs comparison" % model)
        return
    vs1d = np.loadtxt(os.path.join(out, "work", "06_inversion", "vs1d.txt"))
    lon, lat, dep, v = read_model(os.path.join(out, "work", "06_inversion", "main", "DSurfTomo.inMeasure.dat"))
    mean = np.nanmean(v, axis=(0, 1))
    std = np.nanstd(v, axis=(0, 1))
    true_v = np.interp(dep, tm[:, 0], tm[:, 1])
    print("\nVs (km/s) per depth node: start model, inverted (horizontal mean +- std), true")
    for k, z in enumerate(dep):
        print("  %5.2f km   %.3f   %.3f +- %.3f   %.3f   %+5.1f %%" % (z, np.interp(z, vs1d[:, 0], vs1d[:, 1]), mean[k], std[k],
                                                                  true_v[k], 100 * (mean[k] / true_v[k] - 1)))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
    for c in curves:
        axs[0].semilogx(1 / c["T"], c["v"], "-", c="0.75", lw=0.5)
    axs[0].semilogx(f, med, "ko", ms=3, label="median of the measured curves")
    axs[0].semilogx(td[:, 0], td[:, 1], "r-", lw=1.8, label="true dispersion")
    axs[0].set_xlabel("Frequency (Hz)")
    axs[0].set_ylabel("Phase velocity (km/s)")
    axs[0].set_xlim(f.min() * 0.9, f.max() * 1.1)
    from matplotlib.ticker import FuncFormatter, NullFormatter
    axs[0].xaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    axs[0].xaxis.set_minor_formatter(NullFormatter())
    axs[0].legend(frameon=False)
    zt = np.r_[tm[:, 0], tm[-1, 0] * 1.2]
    axs[1].plot(np.r_[tm[:, 1], tm[-1, 1]], zt, "r-", lw=1.8, label="true model")
    axs[1].plot(vs1d[:, 1], vs1d[:, 0], "b--", lw=1.2, label="1-D starting model")
    axs[1].errorbar(mean, dep, xerr=std, fmt="ko", ms=3, label="inverted (horizontal mean, std)")
    axs[1].invert_yaxis()
    axs[1].set_xlabel("Vs (km/s)")
    axs[1].set_ylabel("Depth (km)")
    axs[1].legend(frameon=False, loc="lower left")
    fig.tight_layout()
    p = os.path.join(out, "results", "figures", "truth_comparison.png")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    fig.savefig(p, dpi=150)
    print("\nwrote", p)


if __name__ == "__main__":
    main()
