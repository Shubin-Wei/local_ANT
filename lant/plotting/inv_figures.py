"""Figures of step 4 (4_inv/figures) and of the resolution test (5_checker_*/figures):

start_model.png  the average curve with the fit of the starting model, and the starting model
residuals.png    travel-time residuals of the starting and the final model (kept travel times)
paths.png        the paths of the data at three periods
slices.png       Vs at the slice depths (GMT; plots.slice_depths_km)
profiles.png     Vs along the two principal axes of the station distribution through its centre, bilinear
checker.png      true and recovered checkerboard at the slice depths (GMT)
"""

import os

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .. import geo
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
    paths(pdir, step.fig("paths"), dpi)
    a = np.loadtxt(step.path("vs_model.csv"), delimiter=",", skiprows=1)
    sta = _stations(cfg)
    lon, lat, dep, V = grid_model(a[:, :4])
    inside = _inside(sta, lon, lat, len(dep)) & np.isfinite(V)
    vr = cfg["plots"]["vs_range_kms"] or [float("%.2f" % np.nanpercentile(V[inside], q)) for q in (2, 98)]
    try:
        LO, LA = np.meshgrid(lon, lat, indexing="ij")
        zs = slice_depths(cfg, dep)
        layers = [(LO.ravel(), LA.ravel(), V[:, :, int(np.argmin(np.abs(dep - z)))].ravel()) for z in zs]
        maps.panels(cfg, layers, ["%g km" % z for z in zs], sta, "seis", "%g/%g" % tuple(vr),
                    maps.colorbar_ticks(*vr), "Vs (km/s)", step.fig("slices"), log=step.say)  # fmt: skip
    except Exception as exc:
        step.say("  slices not drawn: %s" % exc)
    profiles(sta, lon, lat, dep, V, vr, step.fig("profiles"), dpi)


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
        ax.legend(fontsize=8)
    style.log_axis(ax.xaxis)
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("Velocity (km/s)")
    bx.plot(vs, z, "k.-")
    bx.invert_yaxis()
    bx.set_xlabel("Vs (km/s)")
    bx.set_ylabel("Depth (km)")
    bx.set_title("last node: fixed half-space", fontsize=9)
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
    ax.legend(fontsize=9)
    style.save(fig, out, dpi)


def paths(pdir, out, dpi):
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
    fig, axs = plt.subplots(1, len(pick), figsize=(4.2 * len(pick), 4.2), squeeze=False)
    for ax, k in zip(axs[0], pick):
        m = a[:, 0] == k
        for r in a[m][:: max(1, m.sum() // 3000)]:
            ax.plot([r[1], r[3]], [r[2], r[4]], color="k", lw=0.2, alpha=0.3)
        ax.set_title("%.3g s %s: %d paths" % (periods[k - 1], "phase" if k <= nrc else "group", m.sum()),
                     fontsize=10)  # fmt: skip
        ax.set_aspect(1 / np.cos(np.radians(np.mean(a[:, 2]))))
        ax.set_xlabel("Longitude")
    axs[0][0].set_ylabel("Latitude")
    style.save(fig, out, dpi)


def axes_of(sta):
    """Centre and the two principal directions (unit vectors in local km) of the station distribution."""
    lon = np.array([s["lon"] for s in sta.values()])
    lat = np.array([s["lat"] for s in sta.values()])
    x, y = geo.local_xy_km(lon, lat)
    w, v = np.linalg.eigh(np.cov(np.c_[x, y].T))
    return (float(lon.mean()), float(lat.mean())), [v[:, 1], v[:, 0]], x, y


def profiles(sta, lon, lat, dep, V, vr, out, dpi):
    (lon0, lat0), dirs, x, y = axes_of(sta)
    kx, ky = geo.km_per_degree(lat0)
    f = RegularGridInterpolator((lon, lat, dep), V, bounds_error=False, fill_value=np.nan)
    fig, axs = plt.subplots(2, 1, figsize=(8.0, 7.0), constrained_layout=True)
    for ax, u, nm in zip(axs, dirs, ("A", "B")):
        proj = x * u[0] + y * u[1]
        s = np.linspace(proj.min(), proj.max(), 300)
        zz = np.linspace(dep[0], dep[-1], 120)
        S, Z = np.meshgrid(s, zz)
        P = np.c_[(lon0 + S.ravel() * u[0] / kx), (lat0 + S.ravel() * u[1] / ky), Z.ravel()]
        img = f(P).reshape(S.shape)
        pm = ax.pcolormesh(s, zz, img, cmap=style.SEIS, vmin=vr[0], vmax=vr[1], shading="gouraud", rasterized=True)
        ax.plot(proj, np.zeros_like(proj), "kv", ms=4, clip_on=False)
        ax.invert_yaxis()
        ax.set_ylabel("Depth (km)")
        ax.set_xlabel("Distance along %s-%s' (km)" % (nm, nm))
        ax.set_title("%s-%s': azimuth %.0f deg through the array centre" % (
            nm, nm, np.degrees(np.arctan2(u[0], u[1])) % 180), fontsize=10)  # fmt: skip
        fig.colorbar(pm, ax=ax, label="Vs (km/s)")
    style.save(fig, out, dpi)


def checker(cfg, step, pdir, z, vs1d):
    true = np.loadtxt(os.path.join(pdir, "Vs_model.real"), ndmin=2)
    syn = np.loadtxt(os.path.join(pdir, "DSurfTomo.inSyn.dat"), ndmin=2)
    dep = np.unique(np.round(syn[:, 2], 4))
    zs = slice_depths(cfg, dep)[:6]
    layers, labels = [], []

    def lay(a, zk):
        m = np.abs(a[:, 2] - zk) < 1e-4
        b = np.interp(zk, z, vs1d)
        return a[m, 0], a[m, 1], 100 * (a[m, 3] / b - 1)

    for zk in zs:
        layers += [lay(true, zk), lay(syn, zk)]
        labels += ["Input %g km" % zk, "Output %g km" % zk]
    amp = 100 * float(cfg["checker"]["amplitude"])
    r = 5 * np.ceil(amp / 5)
    try:
        maps.panels(cfg, layers, labels, _stations(cfg), "polar -I", "-%g/%g" % (r, r), "a%gf%g" % (r / 2, r / 10),
                    "@~D@~Vs (%)", step.fig("checker"), ncol=4, show_stations=True, log=step.say)  # fmt: skip
    except Exception as exc:
        step.say("  checkerboard maps not drawn: %s" % exc)
