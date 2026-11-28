"""Figures of step 4 (4_inv/figures) and of the resolution test (5_checker_*/figures):

start_model.png  the average curve with the fit of the starting model, and the starting model
residuals.png    travel-time residuals of the starting and the final model (kept travel times)
paths.png        the paths of the shortest, a middle and the longest period (GMT)
slices.png       Vs at the slice depths (plots.slice_depths_km) with the two section lines (GMT)
sections.png     Vs along the two principal directions of the station distribution through its centre, with
                 contours, the relief along the line and the stations near it (GMT)
checker.png      the input checkerboard and the recovered one at the slice depths (GMT)
"""

import os

import numpy as np

from .. import stations as ST
from . import maps, style
from .style import plt


def _stations(cfg):
    return ST.read_csv(os.path.join(cfg["output_dir"], "1_clean", "stations_used.csv"))


def grid_model(rows):
    """rows lon lat depth value -> (lon, lat, depth, V [nlon, nlat, nz])."""
    lon, lat, dep = (np.unique(np.round(rows[:, k], 6)) for k in range(3))
    V = np.full((len(lon), len(lat), len(dep)), np.nan)
    i = np.searchsorted(lon, np.round(rows[:, 0], 6))
    j = np.searchsorted(lat, np.round(rows[:, 1], 6))
    k = np.searchsorted(dep, np.round(rows[:, 2], 6))
    V[i, j, k] = rows[:, 3]
    return lon, lat, dep, V


def slice_depths(cfg, dep):
    want = cfg["plots"]["slice_depths_km"]
    if want:
        return [float(dep[np.argmin(np.abs(dep - z))]) for z in want]
    idx = np.unique(np.round(np.linspace(0, len(dep) - 1, min(8, len(dep)))).astype(int))
    return [float(dep[i]) for i in idx]


def make(cfg, step, pdir):
    style.apply()
    dpi = int(cfg["plots"]["dpi"])
    start_model(step, dpi)
    residuals(pdir, step.fig("residuals"), dpi)
    try:
        paths(cfg, pdir, _stations(cfg), step.fig("paths"), step.say)
    except Exception as exc:
        step.say("  paths not drawn: %s" % exc)
    a = np.loadtxt(step.path("vs_model.csv"), delimiter=",", skiprows=1)
    sta = _stations(cfg)
    lon, lat, dep, V = grid_model(a[:, :4])
    inside = _inside(sta, lon, lat, len(dep)) & np.isfinite(V)
    vr = cfg["plots"]["vs_range_kms"] or [float("%.2f" % np.nanpercentile(V[inside], q)) for q in (2, 98)]
    lines = maps.section_lines(sta)
    try:
        LO, LA = np.meshgrid(lon, lat, indexing="ij")
        zs = slice_depths(cfg, dep)
        layers = [(LO.ravel(), LA.ravel(), V[:, :, int(np.argmin(np.abs(dep - z)))].ravel()) for z in zs]
        maps.panels(cfg, layers, ["%g km" % z for z in zs], sta, step.fig("slices"), "seis", "%g/%g" % tuple(vr),
                    maps.colorbar_ticks(*vr), "Vs (km/s)", overlay=maps.line_overlay(lines), log=step.say)  # fmt: skip
    except Exception as exc:
        step.say("  slices not drawn: %s" % exc)
    try:
        maps.sections(cfg, (lon, lat, dep, V), sta, lines, vr, step.fig("sections"), zmax=dep[-2] if len(dep) > 2
                      else dep[-1], log=step.say)  # fmt: skip
    except Exception as exc:
        step.say("  sections not drawn: %s" % exc)


def _inside(sta, lon, lat, nz):
    """Nodes inside the convex hull of the stations [nlon, nlat, nz] (for the colour range)."""
    from matplotlib.path import Path

    ring = maps.hull_ring([s["lon"] for s in sta.values()], [s["lat"] for s in sta.values()])
    LO, LA = np.meshgrid(lon, lat, indexing="ij")
    m = Path(ring).contains_points(np.c_[LO.ravel(), LA.ravel()]).reshape(LO.shape)
    return np.repeat(m[:, :, None], nz, axis=2)


def start_model(step, dpi):
    z, vs = np.loadtxt(step.path("start_model.csv"), delimiter=",", skiprows=1, unpack=True)
    fit = step.path("start_fit.csv")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.0, 4.2))
    if os.path.exists(fit):
        a = np.genfromtxt(fit, delimiter=",", names=True, dtype=None, encoding="utf-8")
        for kind, mk in (("phase", "o"), ("group", "s")):
            m = a["velocity"] == kind
            if m.any():
                ax.plot(a["period_s"][m], a["v_average_kms"][m], "k" + mk, ms=4, mfc="none" if kind == "group" else
                        "k", label="Median of the pairs, %s" % kind)  # fmt: skip
                ax.plot(a["period_s"][m], a["v_start_model_kms"][m], color=style.RED, lw=1.5,
                        ls="-" if kind == "phase" else "--", label="Starting model, %s" % kind)  # fmt: skip
        ax.legend()
    style.log_axis(ax.xaxis)
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("Velocity (km/s)")
    bx.plot(vs, z, "k.-")
    bx.invert_yaxis()
    bx.set_xlabel("Vs (km/s)")
    bx.set_ylabel("Depth (km)")
    bx.text(0.96, 0.04, "last node: fixed half-space", transform=bx.transAxes, ha="right", va="bottom", fontsize=12)
    style.save(fig, step.fig("start_model"), dpi)


def residuals(pdir, out, dpi):
    first = np.loadtxt(os.path.join(pdir, "residualFirst.dat"), ndmin=2)
    last = np.loadtxt(os.path.join(pdir, "residualLast.dat"), ndmin=2)
    r0 = 1e3 * (first[:, 2] - first[:, 1])[first[:, 5] > 0]
    r1 = 1e3 * (last[:, 2] - last[:, 1])[last[:, 5] > 0]
    lim = np.percentile(np.abs(np.r_[r0, r1]), 99.5) if r0.size else 1.0
    bins = np.linspace(-lim, lim, 61)
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.hist(r0, bins, color="0.75", label="starting model, rms %.1f ms" % np.sqrt(np.mean(r0**2)))
    ax.hist(r1, bins, histtype="step", color="k", lw=1.5, label="final model, rms %.1f ms" % np.sqrt(np.mean(r1**2)))
    ax.set_xlabel("Travel-time residual, observed - predicted (ms)")
    ax.set_ylabel("Travel times")
    ax.legend()
    style.save(fig, out, dpi)


def paths(cfg, pdir, sta, out, log=print):
    """The paths of the shortest, a middle and the longest period (GMT)."""
    from ..inversion import dsurf_input

    c = dsurf_input.read_control(os.path.join(pdir, "DSurfTomo.in"))
    nrc = len(c.periods_rc)
    rows, per = [], None
    with open(os.path.join(pdir, "surfdata.dat")) as fh:
        for ln in fh:
            f = ln.split()
            if ln.startswith("#"):  # period block knum: phase blocks first, then group blocks
                src, per = (float(f[1]), float(f[2])), int(f[3]) + (nrc if int(f[5]) == 1 else 0)
            elif f:
                rows.append((per, src[1], src[0], float(f[1]), float(f[0])))
    a = np.array(rows)
    ks = np.unique(a[:, 0]).astype(int)
    pick = sorted({int(ks[int(q * (len(ks) - 1))]) for q in (0.0, 0.5, 1.0)})
    periods = c.periods_all
    labels = ["%.3g s %s, %d paths" % (periods[k - 1], "phase" if k <= nrc else "group", (a[:, 0] == k).sum())
              for k in pick]  # fmt: skip

    def overlay(i, tmp, P):
        m = a[a[:, 0] == pick[i]]
        with open("%s/paths%d.txt" % (tmp, i), "w") as fh:
            for r in m:
                fh.write(">\n%.6f %.6f\n%.6f %.6f\n" % tuple(r[1:]))
        pen, tr = (0.3, 50) if len(m) < 1000 else (0.15, 85)
        return ["gmt plot %s/paths%d.txt %s -W%gp,black@%d" % (tmp, i, P, pen, tr)]

    maps.panels(cfg, [None] * len(pick), labels, sta, out, ncol=3, show_stations=True, overlay=overlay, log=log)


def checker(cfg, step, pdir, z, vs1d):
    true = np.loadtxt(os.path.join(pdir, "Vs_model.real"), ndmin=2)
    syn = np.loadtxt(os.path.join(pdir, "DSurfTomo.inSyn.dat"), ndmin=2)
    dep = np.unique(np.round(syn[:, 2], 4))
    zs = slice_depths(cfg, dep)[:7]
    layers, labels = [], []

    def lay(a, zk):
        m = np.abs(a[:, 2] - zk) < 1e-4
        b = np.interp(zk, z, vs1d)
        return a[m, 0], a[m, 1], 100 * (a[m, 3] / b - 1)

    inputs, seen = [], []
    for zk in zs:  # one input panel per distinct pattern (two when the sign is reversed below a depth)
        v = lay(true, zk)[2]
        if not any(np.corrcoef(v, s0)[0, 1] > 0.9 for s0 in seen):
            seen.append(v)
            inputs.append(zk)
    zs = zs[: 4 * int(np.ceil((len(zs) + len(inputs)) / 4)) - len(inputs)] if len(zs) + len(inputs) > 4 else zs
    for i, zk in enumerate(inputs):
        layers.append(lay(true, zk))
        labels.append("Input" if len(inputs) == 1 else "Input, %g km" % zk)
    for zk in zs:
        layers.append(lay(syn, zk))
        labels.append("%g km" % zk)
    amp = 100 * float(cfg["checker"]["amplitude"])
    r = 5 * np.ceil(amp / 5)
    try:
        maps.panels(cfg, layers, labels, _stations(cfg), step.fig("checker"), "polar -I", "-%g/%g" % (r, r),
                    "a%gf%g" % (r / 2, r / 10), "@~D@~Vs (%)", ncol=4, show_stations=True, log=step.say)  # fmt: skip
    except Exception as exc:
        step.say("  checkerboard maps not drawn: %s" % exc)
