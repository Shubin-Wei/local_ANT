#!/usr/bin/env python
"""Redraw the published 3-km inversion benchmark; do not run or change either solver.

python tools/plot_validation_comparison.py --root VALIDATION/compare_v201 --out docs/assets
Reads the original inv_areal and inv_line cases, scores and raw model files. All paths in the output provenance
are relative to --root. The main map uses a common station-hull display mask; the full-domain companion retains
every finite model value, with outside nodes explicitly hatched. A station hull is not a resolution test.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import ConvexHull, Delaunay

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lant.plotting import style  # noqa: E402

TITLES = ("True", "lant", "DSurfTomo")
DEPTHS = (0.24, 0.48, 0.75)


def stations(path):
    points = []
    for line in path.read_text().splitlines():
        words = line.split()
        if words:
            points.append([float(words[2]), float(words[1])] if words[0] == "#" else
                          [float(words[1]), float(words[0])])  # fmt: skip
    return np.unique(points, axis=0)


def truth_anomaly(x, y, pattern):
    inside = (x >= 0) & (x <= pattern["size_km"]) & (y >= 0) & (y <= pattern["size_km"])
    return (
        100
        * pattern["amplitude"]
        * np.sin(np.pi * x / pattern["cell_km"])
        * np.sin(np.pi * y / pattern["cell_km"])
        * inside
    )


def areal(root, out):
    case = root / "cases_inv_areal/main"
    pattern = json.loads((root / "inv_areal/true_checkerboard.json").read_text())
    if pattern["cell_km"] != 3.0:
        raise ValueError("the published benchmark has 3-km cells; refusing a different synthetic case")
    background = np.loadtxt(root / "inv_areal/true_model.txt")
    scores = json.loads((root / "results_inv_areal.json").read_text())
    st = stations(case / "inputs/surfdata.dat")
    kx, ky = pattern["km_per_degree"]

    def xy(points):
        return (points - [pattern["origin_lon"], pattern["origin_lat"]]) * [kx, ky]

    sxy = xy(st)
    hull = ConvexHull(sxy)
    boundary = sxy[np.r_[hull.vertices, hull.vertices[0]]]
    triangulation = Delaunay(sxy)
    models = [np.loadtxt(case / backend / "DSurfTomo.inMeasure.dat") for backend in ("lant", "native")]
    if not np.array_equal(models[0][:, :3], models[1][:, :3]):
        raise ValueError("model coordinates differ between the two solvers")
    dws = np.loadtxt(case / "lant/DWS.dat")
    if not np.allclose(dws[:, :3], models[0][:, :3], atol=2e-5, rtol=0):
        raise ValueError("DWS coordinates do not match the model")
    panels, diagnostics = [], []
    for depth in DEPTHS:
        select = np.isclose(models[0][:, 2], depth, atol=1e-4)
        points = xy(models[0][select, :2])
        x, y = np.unique(points[:, 0]), np.unique(points[:, 1])
        xx, yy = np.meshgrid(x, y, indexing="ij")
        ij = np.searchsorted(x, points[:, 0]), np.searchsorted(y, points[:, 1])
        inside = triangulation.find_simplex(np.c_[xx.ravel(), yy.ravel()], tol=1e-9).reshape(xx.shape) >= 0
        base = np.interp(depth, background[:, 0], background[:, 1])
        fields = [truth_anomaly(xx, yy, pattern)]
        for model in models:
            field = np.full(xx.shape, np.nan)
            field[ij] = 100 * (model[select, 3] / base - 1)
            fields.append(field)
        outside = ~inside[ij]
        values = models[0][select, 3][outside] / base - 1
        diagnostics.append(dict(depth_km=depth, outside_nodes=int(outside.sum()),
                                outside_zero_dws=int((dws[select, 3][outside] == 0).sum()),
                                outside_anomaly_percent=[float(100 * values.min()), float(100 * values.max())]))  # fmt: skip
        # Display interpolation does not alter stored models or diagnostics on their original nodes.
        dx = np.arange(x.min(), x.max() + 0.0125, 0.025)
        dy = np.arange(y.min(), y.max() + 0.0125, 0.025)
        gx, gy = np.meshgrid(dx, dy, indexing="ij")
        points_fine = np.c_[gx.ravel(), gy.ravel()]
        shown_fields = [truth_anomaly(gx, gy, pattern)] + [
            RegularGridInterpolator((x, y), f, bounds_error=False)(points_fine).reshape(gx.shape) for f in fields[1:]
        ]
        mask = triangulation.find_simplex(points_fine, tol=1e-9).reshape(gx.shape) >= 0
        panels.append((depth, gx, gy, mask, shown_fields))
    amplitude = 100 * pattern["amplitude"]
    for full in (False, True):
        limit = amplitude if not full else np.ceil(max(np.nanmax(np.abs(f)) for p in panels for f in p[4]))
        fig, axes = style.plt.subplots(3, 3, figsize=(10.4, 10.8), sharex=True, sharey=True)
        fig.subplots_adjust(left=0.10, right=0.89, bottom=0.10, top=0.94, wspace=0.10, hspace=0.24)
        for row, (depth, xx, yy, inside, fields) in enumerate(panels):
            for col, field in enumerate(fields):
                ax = axes[row, col]
                ax.set_facecolor("0.91")
                shown = field if full else np.ma.masked_where(~inside, field)
                image = ax.pcolormesh(xx, yy, shown, shading="nearest", cmap=style.PERTURBATION, vmin=-limit, vmax=limit,
                                     rasterized=True)  # fmt: skip
                if full:
                    ax.contourf(xx, yy, (~inside).astype(float), levels=[0.5, 1.5], colors="none", hatches=["///"])
                ax.plot(boundary[:, 0], boundary[:, 1], "k--", lw=0.75)
                ax.plot(sxy[:, 0], sxy[:, 1], "^", ms=3.7, mfc="white", mec="black", mew=0.6)
                ax.set(xlim=(xx.min(), xx.max()) if full else (-0.5, 9.5),
                       ylim=(yy.min(), yy.max()) if full else (-0.5, 9.5),
                       aspect="equal", xticks=[0, 3, 6, 9], yticks=[0, 3, 6, 9])  # fmt: skip
                if row == 0:
                    ax.set_title(TITLES[col], pad=24, fontsize=16)
                ax.text(0, 1.025, "(%s)" % "abcdefghi"[row * 3 + col], transform=ax.transAxes, fontsize=11)
                if col:
                    stat = next(s for s in scores["main/" + ("lant" if col == 1 else "native")]["depths"]
                                if abs(s["z"] - depth) < 1e-4)  # fmt: skip
                    ax.text(1, 1.025, "RMS %.0f m/s" % (1000 * stat["rmse"]), transform=ax.transAxes,
                            ha="right", fontsize=10)  # fmt: skip
                else:
                    ax.set_ylabel("North (km)")
                    ax.text(-0.29, 0.5, "Depth %.2f km" % depth, transform=ax.transAxes, rotation=90,
                            va="center", ha="center", fontsize=12)  # fmt: skip
                if row == 2:
                    ax.set_xlabel("East (km)")
        cax = fig.add_axes([0.91, 0.22, 0.018, 0.56])
        fig.colorbar(image, cax=cax, label="Vs anomaly relative to true 1-D background (%)", extend="both")
        region_note = (
            "Hatched: outside station hull; extrapolation, not resolved structure."
            if full
            else ("Grey: outside station hull (not displayed); full-domain models are provided separately.")
        )
        fig.text(0.10, 0.025, "3-km cells; 64 stations; same original models. Dashed: station hull.\n" + region_note,
                 fontsize=10)  # fmt: skip
        name = "compare_model_full_domain.png" if full else "compare_model.png"
        fig.savefig(out / name, dpi=220)
        style.plt.close(fig)
    return dict(cell_km=pattern["cell_km"], amplitude=pattern["amplitude"], depths_km=DEPTHS,
                display_mask="same station convex hull on all main panels; full-domain companion is unmasked",
                resolution_note="The station hull and nonzero DWS alone do not establish resolution.",
                outside_lant_diagnostics=diagnostics)  # fmt: skip


def line(root, out):
    scores = json.loads((root / "results_inv_line.json").read_text())
    data = scores["main/lant"]
    x, z = np.asarray(scores["points"]["s"]), np.asarray(scores["points"]["depths"])
    st = json.loads((root / "line_stations.json").read_text())["s"]
    base = np.asarray(data["start"])
    fields = [data["truth"], data["section"], scores["main/native"]["section"]]
    fig, axes = style.plt.subplots(1, 3, figsize=(12, 4.8), sharey=True, sharex=True)
    fig.subplots_adjust(left=0.07, right=0.90, bottom=0.25, top=0.80, wspace=0.10)
    for ax, title, field in zip(axes, TITLES, fields):
        image = ax.pcolormesh(x, z, 100 * (np.asarray(field) / base - 1), shading="nearest",
                             cmap=style.PERTURBATION, vmin=-8, vmax=8, rasterized=True)  # fmt: skip
        ax.plot(st, np.full(len(st), -0.04), "v", ms=4, mfc="white", mec="black", mew=0.5, clip_on=False)
        ax.set(ylim=(z.max(), 0), xlim=(x.min(), x.max()), xlabel="Distance along line (km)")
        ax.set_title(title, pad=26)
    axes[0].set_ylabel("Depth (km)")
    cax = fig.add_axes([0.92, 0.27, 0.016, 0.55])
    fig.colorbar(image, cax=cax, label="Vs anomaly relative to starting model (%)")
    fig.text(
        0.07,
        0.08,
        "Same original 3-km synthetic model; 24 stations; no difference panels.\n"
        "Sections show the line between its end stations; benchmarked solver version: lant 2.0.1.",
        fontsize=10,
    )
    fig.savefig(out / "line_sections.png", dpi=220)
    style.plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    style.apply()
    report = areal(args.root, args.out)
    line(args.root, args.out)
    sources = [
        "inv_areal/true_checkerboard.json",
        "inv_areal/true_model.txt",
        "inv_areal/synthetic.json",
        "results_inv_areal.json",
        "results_inv_line.json",
        "line_stations.json",
    ]
    for kind in ("areal", "line"):
        sources += ["cases_inv_%s/main/inputs/%s" % (kind, name) for name in ("DSurfTomo.in", "MOD", "surfdata.dat")]
        sources += [
            "cases_inv_%s/main/%s/%s" % (kind, backend, name)
            for backend in ("lant", "native")
            for name in ("DSurfTomo.inMeasure.dat", "DSurfTomo.in", "MOD", "surfdata.dat")
        ]
    sources.append("cases_inv_areal/main/lant/DWS.dat")
    report.update(
        benchmarked_lant_version="2.0.1",
        titles=TITLES,
        difference_panels=False,
        calculation="existing inversion-level synthetic results; no solver rerun or model editing",
        source_sha256={p: hashlib.sha256((args.root / p).read_bytes()).hexdigest() for p in sources},
    )
    (args.out / "validation_figures.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "source_sha256"}, indent=2))


if __name__ == "__main__":
    main()
