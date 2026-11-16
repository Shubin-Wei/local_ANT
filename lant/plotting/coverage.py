"""Observed source-receiver coverage, explicitly drawn as straight connections.

These are not bent rays from the final velocity model. They show which paths
remain after the final residual rejection at up to six representative periods.
"""

import numpy as np
from matplotlib.collections import LineCollection

from . import style
from .style import plt


def path_coverage(cfg, residuals, out):
    data = np.genfromtxt(residuals, delimiter=",", names=True, ndmin=1)
    if not len(data):
        return None
    periods = np.unique(data["period_s"])
    selected = periods[np.unique(np.round(np.linspace(0, len(periods) - 1, min(6, len(periods)))).astype(int))]
    style.apply(cfg["plots"]["font_scale"])
    cols = min(3, len(selected))
    rows = (len(selected) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4.0 * cols, 3.8 * rows), squeeze=False)
    for ax, period in zip(axes.flat, selected):
        d = data[data["period_s"] == period]
        segments = np.stack([np.c_[d["src_lon"], d["src_lat"]], np.c_[d["rcv_lon"], d["rcv_lat"]]], axis=1)
        kept = d["kept_by_rule"] > 0
        ax.add_collection(LineCollection(segments, colors="0.8", linewidths=0.45, alpha=0.5))
        ax.add_collection(LineCollection(segments[kept], colors=style.BLUE, linewidths=0.55, alpha=0.35))
        points = np.unique(segments.reshape(-1, 2), axis=0)
        ax.scatter(
            points[:, 0], points[:, 1], marker="^", s=18, c="white", edgecolors="black", linewidths=0.5, zorder=3
        )
        ax.autoscale()
        ax.margins(0.06)
        ax.set_aspect(1 / np.cos(np.deg2rad(points[:, 1].mean())))
        ax.ticklabel_format(useOffset=False)
        ax.set_xlabel("Longitude (deg)")
        ax.set_ylabel("Latitude (deg)")
        ax.set_title("T = %g s; retained %d / %d" % (period, kept.sum(), len(d)), fontsize=11)
    for ax in list(axes.flat)[len(selected) :]:
        ax.set_visible(False)
    fig.suptitle("Source-receiver coverage (straight connections)", fontsize=13)
    fig.text(
        0.5,
        0.01,
        "Blue: retained after residual rejection; grey: input paths. Not bent ray trajectories.",
        ha="center",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    return dict(
        files=style.save(fig, out, cfg),
        parameters=dict(
            periods_s=selected.tolist(), path_representation="straight source-receiver connections, not bent rays"
        ),
    )
