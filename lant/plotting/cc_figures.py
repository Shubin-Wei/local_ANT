"""Figures of step 1 (1_cc/figures):

stations.png        the stations with records
ccf_section.png     the stacked correlations against distance, band-passed in the correlation band, each
                    normalised; lines of the slowest and fastest velocity of data.velocity_range_kms.  Surface waves
                    appear as a moveout between them, on both sides when the noise comes from all directions.
ccf_time_<pair>.png the sub-stacks of three pairs (short, middle, long) one above the other, in time order:
                    a stable arrival means a stable correlation; a sub-stack that differs points at a problem of
                    a station on that day (clock, instrument) or a strong local source.
"""

import os

import h5py
import numpy as np
from scipy.signal import butter, sosfiltfilt

from .. import stations as ST
from . import maps, style
from .style import plt

MAX_TRACES = 400


def _filter(x, fs, band):
    sos = butter(4, band, btype="band", fs=fs, output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def make(cfg, step):
    style.apply()
    sta = ST.read_csv(step.path("stations_used.csv"))
    try:
        maps.station_map(cfg, sta, step.fig("stations"), step.say)
    except Exception as exc:  # a missing map must not stop the step: the numbers are written
        step.say("  stations map not drawn: %s" % exc)
    with h5py.File(step.path("stack.h5"), "r") as h:
        fs, maxlag = float(h.attrs["fs"]), float(h.attrs["maxlag"])
        band = (float(h.attrs["fmin"]), min(float(h.attrs["fmax"]), 0.45 * fs))
        cnt = h["count"][:]
        dist = h["dist_km"][:]
        idx = np.where(cnt > 0)[0]
        if not idx.size:
            return
        by_dist = idx[np.argsort(dist[idx])]
        if by_dist.size > MAX_TRACES:  # evenly over distance
            by_dist = by_dist[np.linspace(0, by_dist.size - 1, MAX_TRACES).astype(int)]
        sel = np.sort(by_dist)  # h5py reads increasing indices
        x, d = h["ccf"][sel].astype(np.float64), dist[sel]
        order = np.argsort(d)
        x, d = x[order], d[order]
        names = ["%s - %s" % (a.decode(), b.decode()) for a, b in zip(h["sta1"][:], h["sta2"][:])]
    section(cfg, x, d, fs, maxlag, band, step.fig("ccf_section"))
    q = idx[np.argsort(dist[idx])]
    for p in sorted({q[int(f * (q.size - 1))] for f in (0.25, 0.5, 0.9)}):
        sub_stacks(cfg, step, p, names[p], dist[p], fs, maxlag, band)


def section(cfg, x, d, fs, maxlag, band, out):
    y = _filter(x, fs, band)
    y /= np.abs(y).max(axis=1, keepdims=True) + 1e-30
    lag = (np.arange(x.shape[1]) - x.shape[1] // 2) / fs
    span = (d.max() - d.min()) or 1.0
    amp = 1.5 * span / max(len(d), 1)
    fig, ax = plt.subplots(figsize=(7.5, 6))
    for k in range(len(d)):
        ax.plot(lag, d[k] + amp * y[k], color="k", lw=0.35)
    v = cfg["data"]["velocity_range_kms"]
    if v:
        dd = np.linspace(0, d.max(), 10)
        for vel, ls in ((v[0], "--"), (v[1], ":")):
            ax.plot(dd / vel, dd, color=style.RED, lw=1.0, ls=ls, label="%g km/s" % vel)
            ax.plot(-dd / vel, dd, color=style.RED, lw=1.0, ls=ls)
        ax.legend(loc="lower right", fontsize=9)
    ax.set_xlim(-maxlag, maxlag)
    ax.set_ylim(d.min() - 2 * amp, d.max() + 2 * amp)
    ax.set_xlabel("Lag time (s)")
    ax.set_ylabel("Distance (km)")
    ax.set_title("Stacked correlations, %.3g-%.3g Hz (%d pairs)" % (band[0], band[1], len(d)), fontsize=11)
    style.save(fig, out, int(cfg["plots"]["dpi"]))


def sub_stacks(cfg, step, p, name, dist, fs, maxlag, band):
    days = sorted(f for f in os.listdir(step.path("days")) if f.endswith(".h5"))
    rows, labels = [], []
    for f in days:
        with h5py.File(step.path("days", f), "r") as h:
            gp = h["pair"][:]
            k = np.searchsorted(gp, p)
            if k >= gp.size or gp[k] != p:
                continue
            for s in range(h["count"].shape[0]):
                if h["count"][s, k] > 0:
                    rows.append(h["ccf"][s, k].astype(np.float64))
                    labels.append(
                        f[:-3] if h["count"].shape[0] == 1 else "%s %02dh" % (f[:-3], s * 24 // h["count"].shape[0])
                    )
    if len(rows) < 2:
        return
    y = _filter(np.array(rows), fs, band)
    y /= np.abs(y).max(axis=1, keepdims=True) + 1e-30
    lag = (np.arange(y.shape[1]) - y.shape[1] // 2) / fs
    fig, ax = plt.subplots(figsize=(7.5, 0.25 * len(rows) + 2.0))
    ax.pcolormesh(lag, np.arange(len(rows)), y, cmap="RdBu_r", vmin=-1, vmax=1, shading="nearest", rasterized=True)
    ax.set_yticks(np.arange(len(rows))[:: max(1, len(rows) // 20)])
    ax.set_yticklabels(labels[:: max(1, len(rows) // 20)], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(-maxlag, maxlag)
    ax.set_xlabel("Lag time (s)")
    ax.set_title("%s, %.3f km: sub-stacks, %.3g-%.3g Hz" % (name, dist, band[0], band[1]), fontsize=11)
    style.save(fig, step.fig("ccf_time_" + name.replace(" - ", "_")), int(cfg["plots"]["dpi"]))
