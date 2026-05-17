#!/usr/bin/env python3

"""Plot the number of dispersion curves available at each period point."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Count DSurfTomo project measurements at each period from surfdata.dat "
            "and plot a bar chart."
        )
    )
    parser.add_argument(
        "--project-dir",
        type=Path,
        required=True,
        help="DSurfTomo project directory containing surfdata.dat and surfdata.json.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output figure path.",
    )
    parser.add_argument(
        "--csv-output",
        type=Path,
        default=None,
        help="Optional CSV path for exported counts.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Figure dpi.",
    )
    return parser.parse_args()


def build_period_maps(metadata):
    maps = {}
    for dataset in metadata.get("datasets", []):
        vel_type_name = dataset.get("vel_type_name")
        indices = dataset.get("used_period_indices", [])
        periods = dataset.get("used_periods", [])
        if len(indices) == len(periods) and indices:
            maps[vel_type_name] = {
                int(idx): round(float(period), 6) for idx, period in zip(indices, periods)
            }
    return maps


def read_project_counts(project_dir: Path):
    if not project_dir.is_dir():
        raise FileNotFoundError(f"Project directory not found: {project_dir}")

    surfdata_path = project_dir / "surfdata.dat"
    metadata_path = project_dir / "surfdata.json"
    if not surfdata_path.is_file():
        raise FileNotFoundError(f"surfdata.dat not found: {surfdata_path}")
    if not metadata_path.is_file():
        raise FileNotFoundError(f"surfdata.json not found: {metadata_path}")

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    period_maps = build_period_maps(metadata)
    default_min = round(float(metadata["min_period"]), 6)
    default_step = float(metadata["period_step"])

    counts = {"group": {}, "phase": {}}
    current = None
    with surfdata_path.open("r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                parts = line.split()
                current = {
                    "period_index": int(parts[3]),
                    "vel_type": int(parts[5]),
                }
                continue
            if current is None:
                continue
            vel_type_name = "group" if current["vel_type"] == 1 else "phase"
            period_map = period_maps.get(vel_type_name, {})
            period = period_map.get(
                current["period_index"],
                round(default_min + (current["period_index"] - 1) * default_step, 6),
            )
            counts[vel_type_name][period] = counts[vel_type_name].get(period, 0) + 1

    return counts["group"], counts["phase"], surfdata_path, metadata_path


def apply_publication_style():
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
            "font.size": 12,
            "axes.labelsize": 15,
            "axes.titlesize": 15,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "legend.fontsize": 12,
            "axes.linewidth": 1.1,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.major.size": 5.0,
            "ytick.major.size": 5.0,
        }
    )


def export_csv(path: Path, periods, group_counts, phase_counts):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["period", "group_count", "phase_count"])
        for period in periods:
            writer.writerow([f"{period:.6f}", group_counts.get(period, 0), phase_counts.get(period, 0)])


def plot_counts(periods, group_counts, phase_counts, output_path: Path, dpi: int):
    apply_publication_style()
    x = np.arange(len(periods), dtype=float)
    width = 0.72
    fig_width = max(9.0, min(14.0, 6.5 + 0.12 * len(periods)))
    fig, ax = plt.subplots(figsize=(fig_width, 4.8), constrained_layout=True)

    group_values = np.array([group_counts.get(period, 0) for period in periods], dtype=float)
    phase_values = np.array([phase_counts.get(period, 0) for period in periods], dtype=float)

    if phase_counts:
        ax.bar(
            x,
            phase_values,
            width=width,
            color="0.72",
            edgecolor="0.35",
            linewidth=0.45,
            label="Phase",
            zorder=2,
        )
    if group_counts:
        ax.bar(
            x,
            group_values,
            width=width * 0.62,
            color="0.30",
            edgecolor="0.10",
            linewidth=0.45,
            label="Group",
            zorder=3,
        )

    ax.set_xlabel("Period (s)")
    ax.set_ylabel("Counts")
    ax.set_xlim(-0.8, len(periods) - 0.2)
    tick_step = max(1, int(np.ceil(len(periods) / 12.0)))
    tick_idx = np.arange(0, len(periods), tick_step, dtype=int)
    if len(periods) - 1 not in tick_idx:
        tick_idx = np.append(tick_idx, len(periods) - 1)
    ax.set_xticks(tick_idx)
    ax.set_xticklabels([f"{periods[idx]:.2f}" for idx in tick_idx])
    ax.tick_params(direction="out", top=True, right=True, pad=4)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color="0.88", linewidth=0.8)
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)
    ax.legend(
        loc="upper right",
        frameon=True,
        fancybox=False,
        edgecolor="black",
        framealpha=1.0,
        borderpad=0.6,
        handlelength=1.8,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    group_counts, phase_counts, surfdata_path, metadata_path = read_project_counts(
        args.project_dir.resolve()
    )

    periods = sorted(set(group_counts) | set(phase_counts))
    if not periods:
        raise ValueError("No periods found in the current project.")

    output_path = args.output.resolve()
    print(f"project_dir: {args.project_dir.resolve()}")
    print(f"surfdata: {surfdata_path}")
    print(f"metadata: {metadata_path}")
    print(f"periods: {len(periods)}")
    print(f"period_range: {periods[0]:.3f} - {periods[-1]:.3f}")
    print(f"output: {output_path}")

    if args.csv_output is not None:
        csv_path = args.csv_output.resolve()
        export_csv(csv_path, periods, group_counts, phase_counts)
        print(f"csv_output: {csv_path}")

    plot_counts(periods, group_counts, phase_counts, output_path, args.dpi)
    print("done")


if __name__ == "__main__":
    main()
