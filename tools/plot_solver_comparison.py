#!/usr/bin/env python
"""Draw existing native/Python checkerboard or real-line comparisons in the README style; runs no solver.

python tools/plot_solver_comparison.py --case VALIDATION_DIR --out FIGURE_DIR
The case contains inputs/{DSurfTomo.in,MOD,surfdata.dat}, case.json, and native/ and python/ solver outputs.
Checkerboard truth may be stored separately when both solvers read shared measured-mode observations.
Writes checkerboard.png or section.png, timing.png and summary.json. Native timings use the external gdb log.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lant.inversion import dsurf_input  # noqa: E402
from lant.inversion.run import read_model  # noqa: E402
from lant.plotting import style  # noqa: E402
from lant.signal import local_xy_km  # noqa: E402


def trimmed(values):
    values = sorted(values)
    return float(np.mean(values[1:-1] if len(values) > 2 else values))


def station_pairs(path):
    pairs, source = [], None
    for line in path.read_text().splitlines():
        words = line.split()
        if not words:
            continue
        if words[0] == "#":
            source = [float(words[2]), float(words[1])]
        else:
            pairs.append([source, [float(words[1]), float(words[0])]])
    return np.unique(np.array(pairs), axis=0)


def line_figure(fields, axes, stations, out, iterations):
    """Initial model and the two solvers on the central longitude column; no difference panel."""
    lon, lat, dep = axes
    centre_lat, centre_lon = np.mean(stations[:, 1]), np.mean(stations[:, 0])
    j = int(np.argmin(abs(lon - centre_lon)))
    _, y = local_xy_km(np.full(len(lat), lon[j]), lat, centre_lat, centre_lon)
    _, sy = local_xy_km(stations[:, 0], stations[:, 1], centre_lat, centre_lon)
    origin = sy.min()
    y, sy = y - origin, sy - origin
    sections = [fields[i][j].T for i in (0, 2, 1)]
    low, high = min(float(v.min()) for v in sections), max(float(v.max()) for v in sections)
    fig, axs = style.plt.subplots(1, 3, figsize=(13, 4.2), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.07, right=0.89, bottom=0.23, top=0.88, wspace=0.12)
    for ax, title, values in zip(axs, ("Initial", "lant", "DSurfTomo"), sections):
        image = ax.pcolormesh(y, dep, values, shading="nearest", cmap=style.PERTURBATION, vmin=low, vmax=high)
        ax.plot(sy, np.zeros(len(sy)), "v", ms=4, mfc="white", mec="black", mew=0.5, zorder=5)
        ax.set(
            ylim=(dep[-1], -0.025), xlim=(y.min(), y.max()), title=title, xlabel="North from southernmost station (km)"
        )
    axs[0].set_ylabel("Depth (km)")
    cax = fig.add_axes([0.92, 0.25, 0.015, 0.58])
    fig.colorbar(image, cax=cax, label="Vs (km/s)")
    fig.text(
        0.07,
        0.04,
        "Common observations; %d stations; %d iterations. Central longitude %.5f°."
        % (len(stations), iterations, lon[j]),
        fontsize=9,
    )
    fig.savefig(out / "section.png", dpi=180)
    style.plt.close(fig)
    return dict(section_longitude=float(lon[j]), section="central longitude column; north-south strip")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--depth", type=float, default=0.6)
    args = parser.parse_args()
    case, out = Path(args.case), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    description = json.loads((case / "case.json").read_text())
    inputs = case / "inputs"
    control = dsurf_input.read_control(inputs / "DSurfTomo.in")
    shared_observations = not bool(control.ifsyn)
    if shared_observations:
        for backend in ("native", "python"):
            for name in ("DSurfTomo.in", "MOD", control.datafile):
                if (case / backend / name).read_bytes() != (inputs / name).read_bytes():
                    raise ValueError("Common-input comparison has different %s/%s" % (backend, name))
        observed = [
            np.loadtxt(case / backend / "residualFirst.dat")[:, 2].astype(np.float32)
            for backend in ("native", "python")
        ]
        if not np.array_equal(*observed):
            raise ValueError("The two solvers did not read identical observations")
    z, start = dsurf_input.read_model(inputs / "MOD", control.nx, control.ny, control.nz)
    filename = "DSurfTomo.inSyn.dat" if control.ifsyn else "DSurfTomo.inMeasure.dat"
    truth = case / "native/Vs_model.real" if control.ifsyn else case / "Vs_model.real"
    model_paths = ([truth] if truth.exists() else []) + [case / backend / filename for backend in ("native", "python")]
    models = [read_model(str(p)) for p in model_paths]
    lon, lat, dep = models[0][:3]
    for model in models[1:]:
        if not all(np.array_equal(a, b) for a, b in zip(model[:3], (lon, lat, dep))):
            raise ValueError("Model coordinates differ")
    fields = [m[3] for m in models]
    if not truth.exists():
        fields.insert(0, start[:-1, 1:-1, 1:-1].transpose(1, 2, 0))
    delta = fields[2] - fields[1]
    k = int(np.argmin(abs(dep - args.depth)))
    base = float(np.mean(start[int(np.argmin(abs(z - dep[k])))]))
    pairs = station_pairs(inputs / "surfdata.dat")
    stations = np.unique(pairs.reshape(-1, 2), axis=0)
    centre_lon = description.get("lon0", float(np.mean(stations[:, 0])))
    centre_lat = description.get("lat0", float(np.mean(stations[:, 1])))
    gx, _ = local_xy_km(lon, np.full(len(lon), centre_lat), centre_lat, centre_lon)
    _, gy = local_xy_km(np.full(len(lat), centre_lon), lat, centre_lat, centre_lon)
    X, Y = np.meshgrid(gx, gy, indexing="ij")
    sx, sy = local_xy_km(stations[:, 0], stations[:, 1], centre_lat, centre_lon)
    px, py = local_xy_km(pairs[:, :, 0], pairs[:, :, 1], centre_lat, centre_lon)
    style.apply()
    plot_metadata = {}
    if truth.exists():
        fig, axs = style.plt.subplots(1, 3, figsize=(13.2, 4.9))
        fig.subplots_adjust(left=0.07, right=0.89, bottom=0.22, top=0.88, wspace=0.17)
        amplitude = 100 * description["amplitude"]
        for i, (ax, title, index) in enumerate(zip(axs, ("True", "lant", "DSurfTomo"), (0, 2, 1))):
            values = 100 * (fields[index][:, :, k] / base - 1)
            image = ax.pcolormesh(
                X, Y, values, shading="gouraud", cmap=style.PERTURBATION, vmin=-amplitude, vmax=amplitude
            )
            ax.plot(sx, sy, "^", ms=4.5, mfc="white", mec="black", mew=0.6, zorder=5)
            ax.set(
                xlim=(sx.min() - 0.75, sx.max() + 0.75),
                ylim=(sy.min() - 0.75, sy.max() + 0.75),
                aspect="equal",
                xlabel="East (km)",
                title=title,
            )
            if i == 0:
                ax.set_ylabel("North (km)")
        cax = fig.add_axes([0.92, 0.27, 0.014, 0.53])
        fig.colorbar(image, cax=cax, label="Vs anomaly (%)")
        fig.text(
            0.07,
            0.06,
            "%.2f km depth; %d stations; %d iterations. Shared ±%g%% scale.\n"
            "Historical solver-port comparison; same saved models, no solver rerun or difference panels."
            % (dep[k], len(stations), control.maxiter, amplitude),
            fontsize=9,
        )
        fig.savefig(out / "checkerboard.png", dpi=180)
        style.plt.close(fig)

    else:
        plot_metadata = line_figure(fields, (lon, lat, dep), stations, out, control.maxiter)

    calls = [json.loads(line) for line in (case / "native/timing.jsonl").read_text().splitlines()]
    times = {
        name: [v["end"] - v["start"] for v in calls if v["name"] == name]
        for name in ("depthkernel", "calsurfg", "lsmr")
    }
    native = [times["depthkernel"], [c - d for c, d in zip(times["calsurfg"], times["depthkernel"])], times["lsmr"]]
    record = json.loads((case / "python/solver.json").read_text())
    py = [[row[name] for row in record["iterations"]] for name in ("senk", "traveltime", "lsmr")]
    values = np.array([[trimmed(v) for v in native], [trimmed(v) for v in py]])
    fig, axes = style.plt.subplots(1, 3, figsize=(10.4, 3.7))
    for j, (ax, name) in enumerate(zip(axes, ("Depth kernels", "Travel times", "LSMR"))):
        bars = ax.bar([0, 1], values[:, j], color=["0.65", style.BLUE], width=0.6)
        ax.bar_label(bars, fmt="%.2f", padding=3, fontsize=11)
        ax.set_xticks([0, 1], ["Native", "Python"])
        ax.set_ylim(0, max(values[:, j]) * 1.28)
        ax.set_title("(%s) %s" % ("abc"[j], name))
        if j == 0:
            ax.set_ylabel("Seconds per iteration")
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    fig.text(
        0.06,
        0.02,
        "%d threads; native travel times/LSMR are serial. Mean excludes the fastest and slowest iteration; gdb overhead is included."
        % record["threads"],
        fontsize=9,
    )
    fig.savefig(out / "timing.png", dpi=180)
    style.plt.close(fig)
    sources = model_paths + sorted(p for p in inputs.iterdir() if p.is_file()) + [case / "case.json"]
    if (case / "MOD.true").exists():
        sources.append(case / "MOD.true")
    summary = dict(
        depth_km=float(dep[k]) if truth.exists() else None,
        plot=plot_metadata,
        observation_protocol=description.get(
            "input_protocol", "Each implementation generates synthetic observations independently."
        ),
        shared_observations=shared_observations,
        final_max_abs_kms=float(abs(delta).max()),
        final_rms_kms=float(np.sqrt(np.mean(delta**2))),
        iterations=control.maxiter,
        threads=record["threads"],
        native_trimmed_seconds=values[0].tolist(),
        python_trimmed_seconds=values[1].tolist(),
        python_identity=record["identity"],
        sha256={str(p.relative_to(case)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
    )
    run_summary = case / "run_summary.json"
    if run_summary.exists():
        summary["comparison"] = json.loads(run_summary.read_text())
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k not in ("python_identity", "sha256")}, indent=2))


if __name__ == "__main__":
    main()
