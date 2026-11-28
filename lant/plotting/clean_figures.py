"""Figures of step 1 (1_clean/figures):

stations.png      the stations with cleaned data
availability.png  hours of cleaned data per station and day: gaps, late starts, stations that stopped
"""

import csv

import numpy as np

from .. import stations as ST
from . import maps, style
from .style import plt


def make(cfg, step):
    style.apply()
    sta = ST.read_csv(step.path("stations_used.csv"))
    try:
        maps.station_map(cfg, sta, step.fig("stations"), step.say)
    except Exception as exc:  # a missing map must not stop the step: the numbers are written
        step.say("  stations map not drawn: %s" % exc)
    availability(step.path("availability.csv"), step.fig("availability"), cfg["plots"]["dpi"])


def availability(table, out, dpi=300):
    with open(table, newline="") as fh:
        rows = list(csv.reader(fh))
    if len(rows) < 2 or len(rows[0]) < 2:
        return
    days, names = rows[0][1:], [r[0] for r in rows[1:]]
    h = np.array([[float(v) for v in r[1:]] for r in rows[1:]])
    ns, nd = h.shape
    fig, ax = plt.subplots(figsize=(min(16, 4 + 0.18 * nd), min(14, 2.5 + 0.12 * ns)), constrained_layout=True)
    im = ax.imshow(h, aspect="auto", interpolation="nearest", cmap="viridis", vmin=0, vmax=24,
                   extent=(-0.5, nd - 0.5, ns - 0.5, -0.5))  # fmt: skip
    step = max(1, nd // 15)
    ax.set_xticks(range(0, nd, step))
    ax.set_xticklabels([days[k] for k in range(0, nd, step)], rotation=60, ha="right")
    if ns <= 80:
        ax.set_yticks(range(ns))
        ax.set_yticklabels(names, fontsize=6 if ns > 40 else 8)
    else:
        ax.set_ylabel("station (%d, sorted by name)" % ns)
    ax.set_xlabel("UTC day")
    fig.colorbar(im, ax=ax, label="hours of cleaned data", shrink=0.8)
    ax.set_title(
        "Cleaned data (stations: %d, days: %d, median %.1f h per station-day)" % (ns, nd, np.median(h)), fontsize=10
    )
    style.save(fig, out, dpi)
