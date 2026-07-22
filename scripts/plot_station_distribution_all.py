#!/usr/bin/env python3
"""Plot station distributions for an arbitrary set of line arrays."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


COLORS = ["#D62728", "#1F77B4", "#2CA02C", "#9467BD", "#FF7F0E", "#17BECF"]
MARKERS = ["^", "v", "o", "s", "D", "P"]


def read_station_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    return (
        np.array([float(row["longitude"]) for row in rows]),
        np.array([float(row["latitude"]) for row in rows]),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-root", default="LineData",
                   help="directory containing LineN/station.csv")
    p.add_argument("--lines", nargs="+", type=int, default=[1, 2, 3, 4, 5, 6])
    p.add_argument("--output", default="LineInversion/station_distribution.png")
    p.add_argument("--title", help="custom figure title")
    args = p.parse_args()

    datasets = []
    for line in args.lines:
        path = Path(args.data_root) / f"Line{line}" / "station.csv"
        lon, lat = read_station_csv(path)
        datasets.append((line, lon, lat))

    all_lon = np.concatenate([item[1] for item in datasets])
    all_lat = np.concatenate([item[2] for item in datasets])
    mean_lat = float(all_lat.mean())

    fig, ax = plt.subplots(figsize=(10.0, 9.0))
    for index, (line, lon, lat) in enumerate(datasets):
        ax.scatter(
            lon, lat, s=31, marker=MARKERS[index % len(MARKERS)],
            facecolor=COLORS[index % len(COLORS)], edgecolor="black",
            linewidth=0.4, zorder=2, label=f"Line{line} (n={len(lon)})")

    lat_pad = max(0.002, 0.025 * float(np.ptp(all_lat)))
    ymin = float(all_lat.min() - lat_pad)
    ymax = float(all_lat.max() + lat_pad)
    lat_span = ymax - ymin
    lon_center = float(0.5 * (all_lon.min() + all_lon.max()))
    geographic_half_width = 0.5 * lat_span / np.cos(np.radians(mean_lat))
    observed_half_width = 0.5 * float(np.ptp(all_lon)) * 1.10
    lon_half_width = max(geographic_half_width, observed_half_width)
    ax.set_xlim(lon_center - lon_half_width, lon_center + lon_half_width)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect(1.0 / np.cos(np.radians(mean_lat)), adjustable="box")

    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")
    title = args.title or f"Station Distribution: Line{args.lines[0]}–Line{args.lines[-1]}"
    ax.set_title(title, fontsize=13, pad=9)
    ax.grid(True, color="0.85", lw=0.7, ls="--")
    ax.legend(loc="upper right", ncol=2, frameon=True, framealpha=0.95, fontsize=8.5)
    ax.ticklabel_format(axis="both", style="plain", useOffset=False)
    fig.tight_layout()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    pdf = output.with_suffix(".pdf")
    fig.savefig(pdf, bbox_inches="tight")
    print(f"Saved {output}")
    print(f"Saved {pdf}")


if __name__ == "__main__":
    main()
