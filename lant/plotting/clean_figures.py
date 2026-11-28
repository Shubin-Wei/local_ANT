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
    fig = plt.figure(figsize=(min(11.0, 5.0 + 0.16 * nd), min(9.0, 3.2 + 0.035 * ns)))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 4], width_ratios=[40, 1], hspace=0.06, wspace=0.03)
    ax0, ax, cax = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    ax0.bar(np.arange(nd), (h > 0).sum(axis=0), width=0.8, color="0.45")
    ax0.set_xlim(-0.5, nd - 0.5)
    ax0.set_ylabel("Stations")
    ax0.tick_params(labelbottom=False)
    im = ax.imshow(h, aspect="auto", interpolation="nearest", cmap="Blues", vmin=0, vmax=24,
                   extent=(-0.5, nd - 0.5, ns - 0.5, -0.5))  # fmt: skip
    every = max(1, int(np.ceil(nd / 12)))
    ax.set_xticks(range(0, nd, every))
    ax.set_xticklabels([days[k][5:] for k in range(0, nd, every)])
    ax.set_xlabel("UTC day (%s\u2013%s)" % (days[0], days[-1]))
    if ns <= 50:
        ax.set_yticks(range(ns))
        ax.set_yticklabels(names, fontsize=8)
    else:
        ax.set_ylabel("Station (%d, sorted by name)" % ns)
    fig.colorbar(im, cax=cax, label="Hours of cleaned data", ticks=[0, 6, 12, 18, 24])
    style.save(fig, out, dpi)
