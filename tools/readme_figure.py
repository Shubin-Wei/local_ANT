#!/usr/bin/env python
"""The overview figure of the README (docs/assets/overview.png), from a synthetic array with a known answer.

A random areal array records plane-wave noise of the 1-D "rock" model (lant.synthetic); its stack is then replaced
by correlations of a laterally varying medium, a +-6 % chequer pattern of 3-km cells with 2 % incoherent noise
(examples/synthetic/make_stack.py), and the chain runs from the reference step to the end.  The gates are confirmed
because the data are synthetic.  The figure shows the input pattern with the ray paths, the recovered model, the
phase-velocity measurements and the 1-D models.

    python tools/readme_figure.py --work DIR [--out docs/assets/overview.png]

Steps already done in DIR are not repeated; about 20 minutes on a shared server.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

STATIONS, HOURS, HETERO, CELL, NOISE = 36, 4, 0.06, 3.0, 0.02
DEPTH_KM = 0.6  # depth of the map panels
REASON = "Generated synthetic input with known truth (README figure)"


def lant(*args, check=True):
    return subprocess.run([sys.executable, "-m", "lant"] + list(args), cwd=str(ROOT), check=check).returncode


def run_with_gates(config, first=(), keep=()):
    """lant run all (with `first` only on the first call), confirming the gates of the synthetic data."""
    args = list(first) + list(keep)
    for _ in range(8):
        code = lant("run", "all", "-c", config, *args, check=False)
        if code == 0:
            return
        if code != 3:
            raise SystemExit(code)
        state = json.loads(
            subprocess.check_output([sys.executable, "-m", "lant", "status", "-c", config, "--json"], cwd=str(ROOT))
        )
        gate = next(step for step, value in state["steps"].items() if value["state"] == "to confirm")
        lant("confirm", gate, "-c", config, "--reason", REASON)
        args = list(keep)
    raise SystemExit("too many confirmation stops")


def make_data(work):
    from lant import synthetic

    config = str(work / "synthetic.yaml")
    if not (work / "synthetic.yaml").exists():
        synthetic.main(["--out", str(work), "--stations", str(STATIONS), "--days", "1", "--hours", str(HOURS)])
    stack = work / "output" / "work" / "03_stack" / "stack.h5"
    truth = work / "truth.npz"
    if not truth.exists():
        run_with_gates(config, keep=("--to", "stack"))
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "examples/synthetic/make_stack.py"),
                "--like",
                str(stack),
                "--out",
                str(stack),
                "--model",
                str(work / "true_model.txt"),
                "--hetero",
                str(HETERO),
                "--cell",
                str(CELL),
                "--noise",
                str(NOISE),
                "--truth",
                str(truth),
            ],
            check=True,
        )
        run_with_gates(config, first=("--from", "reference", "--yes"))
    else:
        run_with_gates(config)
    return work / "output", truth


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--work", required=True, type=Path, help="a directory for the synthetic project")
    parser.add_argument("--out", default=str(ROOT / "docs" / "assets" / "overview.png"))
    parser.add_argument("--data-only", action="store_true", help="make and run the project, do not draw")
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    output, truth = make_data(args.work.resolve())
    if not args.data_only:
        draw(output, truth, args.work.resolve(), args.out)


def _read_csv(path):
    import csv

    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def draw(output, truth, work, out):
    import matplotlib
    import matplotlib.ticker

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    from lant.plotting.style import SEIS, apply
    from lant.signal import local_xy_km

    t = np.load(str(truth))
    lon0, lat0, cell, amp = float(t["lon0"]), float(t["lat0"]), float(t["cell"]), float(t["hetero"])

    stations = _read_csv(output / "results/reproduce/stations.csv")
    sx, sy = local_xy_km(
        [float(r["longitude"]) for r in stations], [float(r["latitude"]) for r in stations], lat0, lon0
    )
    rows = [r for r in _read_csv(output / "results/dispersion/curves.csv") if r["used_in_inversion"] == "1"]
    pairs = {r["pair"]: tuple(float(r[k]) for k in ("lon1", "lat1", "lon2", "lat2")) for r in rows}

    model = np.array([[float(v) for v in r.values()] for r in _read_csv(output / "results/model/vs_model.csv")])
    depths = np.unique(model[:, 2])
    depth = depths[np.argmin(abs(depths - DEPTH_KM))]
    layer = model[model[:, 2] == depth]
    lons, lats = np.unique(layer[:, 0]), np.unique(layer[:, 1])
    grid = np.full((lats.size, lons.size), np.nan)
    grid[np.searchsorted(lats, layer[:, 1]), np.searchsorted(lons, layer[:, 0])] = layer[:, 3]
    gx, _ = local_xy_km(lons, np.full(lons.size, lat0), lat0, lon0)
    _, gy = local_xy_km(np.full(lats.size, lon0), lats, lat0, lon0)
    inside = (
        (gx[None, :] >= sx.min()) & (gx[None, :] <= sx.max()) & (gy[:, None] >= sy.min()) & (gy[:, None] <= sy.max())
    )
    recovered = 100 * (grid / np.nanmean(grid[inside]) - 1)
    X, Y = np.meshgrid(gx, gy)
    true = 100 * amp * np.sin(np.pi * X / cell) * np.sin(np.pi * Y / cell)

    apply(0.9)
    fig = plt.figure(figsize=(10.4, 8.4))
    grid_spec = fig.add_gridspec(2, 2, wspace=0.25, hspace=0.3)
    ax_a, ax_b = fig.add_subplot(grid_spec[0, 0]), fig.add_subplot(grid_spec[0, 1])
    ax_c, ax_d = fig.add_subplot(grid_spec[1, 0]), fig.add_subplot(grid_spec[1, 1])
    lim = 100 * amp
    for ax, field in ((ax_a, true), (ax_b, recovered)):
        image = ax.pcolormesh(X, Y, field, shading="gouraud", cmap=SEIS, vmin=-lim, vmax=lim, rasterized=True)
        ax.set_aspect("equal")
        ax.set_xlim(sx.min() - 1, sx.max() + 1)  # the array and a margin, not the padding of the model grid
        ax.set_ylim(sy.min() - 1, sy.max() + 1)
        ax.set_xlabel("East (km)")
    for lo1, la1, lo2, la2 in pairs.values():
        (x1, x2), (y1, y2) = local_xy_km([lo1, lo2], [la1, la2], lat0, lon0)
        ax_a.plot([x1, x2], [y1, y2], color="0.15", lw=0.25, alpha=0.35)
    for ax in (ax_a, ax_b):
        ax.plot(sx, sy, "^", ms=6.5, mfc="white", mec="black", mew=0.9, zorder=5)
    ax_a.set_ylabel("North (km)")
    ax_a.set_title("Input Vs anomaly, %d station pairs" % len(pairs), fontsize=11)
    ax_b.set_title("Recovered by lant, %.1f km depth" % depth, fontsize=11)
    bar = fig.colorbar(image, ax=[ax_a, ax_b], fraction=0.025, pad=0.02)
    bar.set_label("Vs anomaly (%)")

    by_pair = {}
    for r in rows:
        by_pair.setdefault(r["pair"], []).append((float(r["frequency_hz"]), float(r["inversion_velocity_kms"])))
    for values in by_pair.values():
        f, c = np.array(sorted(values)).T
        ax_c.plot(f, c, color="0.6", lw=0.5, alpha=0.6)
    ref = np.load(str(output / "work/04_reference/reference.npz"))
    valid = ref["ref_valid"].astype(bool)
    band = (min(min(v)[0] for v in by_pair.values()), max(max(v)[0] for v in by_pair.values()))
    keep = valid & (ref["freqs"] >= band[0]) & (ref["freqs"] <= band[1])  # the measured band only
    ax_c.plot(ref["freqs"][keep], ref["c_ref"][keep], color="black", lw=1.8, label="reference curve (F-J)")
    inband = (t["freqs"] >= band[0]) & (t["freqs"] <= band[1])
    ax_c.plot(t["freqs"][inband], t["c0"][inband], color=SEIS(0.08), lw=1.4, ls="--", label="true background")
    ax_c.plot([], [], color="0.6", lw=0.8, label="%d two-station curves" % len(by_pair))
    ax_c.set_xscale("log")
    ax_c.set_xticks([t_ for t_ in (0.5, 1, 2, 5, 10) if band[0] <= t_ <= band[1]])
    ax_c.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax_c.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax_c.set_xlabel("Frequency (Hz)")
    ax_c.set_ylabel("Phase velocity (km/s)")
    ax_c.legend(frameon=False, fontsize=9, loc="upper right")
    ax_c.set_title("Rayleigh-wave phase velocities", fontsize=11)

    start = np.loadtxt(str(output / "results/reproduce/vs1d.txt"))
    mean = np.array([np.nanmean(model[model[:, 2] == z, 3]) for z in depths])
    spread = np.array([np.nanstd(model[model[:, 2] == z, 3]) for z in depths])
    ax_d.plot(t["vs"], t["z"], color=SEIS(0.08), lw=1.4, ls="--", label="true background")
    ax_d.plot(start[:, 1], start[:, 0], color="0.45", lw=1.2, label="1-D starting model")
    ax_d.fill_betweenx(depths, mean - spread, mean + spread, color=SEIS(0.85), alpha=0.25, lw=0)
    ax_d.plot(mean, depths, color=SEIS(0.9), lw=1.8, label="3-D model: mean and spread")
    ax_d.set_ylim(depths.max(), 0)
    ax_d.set_xlabel("Vs (km/s)")
    ax_d.set_ylabel("Depth (km)")
    ax_d.legend(frameon=False, fontsize=9, loc="lower left")
    ax_d.set_title("1-D models", fontsize=11)

    for ax, label in zip((ax_a, ax_b, ax_c, ax_d), "abcd"):
        ax.text(-0.02, 1.02, "(%s)" % label, transform=ax.transAxes, ha="right", va="bottom", fontweight="bold")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
