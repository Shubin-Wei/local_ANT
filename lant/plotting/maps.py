"""Map figures drawn with GMT: station map, depth slices of the model, checkerboard maps.

Style: light grey hill-shaded relief (plots.topography; plain background
without it), model values resampled bilinearly for display (node values unchanged) and drawn 35 % transparent,
clipped to the convex hull of the stations, GMT "seis" colour table red (slow) to blue (fast), stations as white
triangles, one plain scale bar, one colour bar for the figure.  No geological information is drawn.
"""
import os

import numpy as np
from scipy.spatial import ConvexHull

from . import gmt as G


def _region(lon, lat, pad=0.1, min_pad=0.005):
    W, E, S, N = min(lon), max(lon), min(lat), max(lat)
    px = max((E - W) * pad, min_pad)
    py = max((N - S) * pad, min_pad)
    return W - px, E + px, S - py, N + py


def _ticks(W, E, S, N):
    return G.nice_step(E - W, 3), G.nice_step(N - S, 4)


def station_map(cfg, stations, out, title=None):
    lon = [s["lon"] for s in stations.values()]
    lat = [s["lat"] for s in stations.values()]
    W, E, S, N = _region(lon, lat)
    tmp = G.workdir()
    np.savetxt(tmp + "/sta.xy", np.c_[lon, lat], fmt="%.6f")
    width = 14.0
    P = "-R%.5f/%.5f/%.5f/%.5f -JM%gc -Xa2.2c -Ya1.6c" % (W, E, S, N, width)
    pre, relief = G.relief_cmds(cfg, tmp, W, E, S, N, palette="elevation")
    fs = cfg["plots"]["font_scale"]
    cmd = []
    if relief:
        cmd.append("gmt grdimage %s/topo.nc -I%s/shade.nc -C%s/relief.cpt %s" % (tmp, tmp, tmp, P))
    cmd.append("gmt plot %s/sta.xy %s -St%.2fc -G#1f3864 -W0.5p,white" % (tmp, P, 0.36 * fs))
    dx, dy = _ticks(W, E, S, N)
    km, latc = G.scale_km(W, E, S, N)
    cmd += ["gmt basemap %s -BWSne -Bxa%gf%g -Bya%gf%g --FONT_ANNOT_PRIMARY=%gp,Helvetica" % (P, dx, dx / 4, dy, dy / 4, 12 * fs),
            "gmt basemap %s -LjBR+c%.4f+w%gk+f+l'km'+ar+o0.8c/0.6c --MAP_SCALE_HEIGHT=0.2c --FONT_ANNOT_PRIMARY=%gp,Helvetica "
            "--FONT_LABEL=%gp,Helvetica" % (P, latc, km, 11 * fs, 11 * fs)]
    open(tmp + "/legend.txt", "w").write("S 0.3c t %.2fc #1f3864 0.5p,white 0.8c %s (%d)\n"
                                         % (0.36 * fs, title or "Stations", len(lon)))
    cmd.append("gmt legend %s/legend.txt %s -DJTR+jTL+w6c+o0.3c/0c --FONT_ANNOT_PRIMARY=%gp,Helvetica" % (tmp, P, 12 * fs))
    if relief:
        cmd.append("gmt colorbar -C%s/relief.cpt %s -DJMR+jML+w5c/0.35c+o0.6c/-2c -Bxaf+l'Elevation (m)' "
                   "--FONT_ANNOT_PRIMARY=%gp,Helvetica --FONT_LABEL=%gp,Helvetica" % (tmp, P, 11 * fs, 12 * fs))
    return G.run(tmp, pre, cmd, out, cfg["plots"]["formats"], cfg["plots"]["dpi"])


def _hull(stations, tmp):
    xy = np.array([(s["lon"], s["lat"]) for s in stations.values()])
    h = ConvexHull(xy)
    ring = xy[np.r_[h.vertices, h.vertices[:1]]]
    np.savetxt(tmp + "/hull.xy", ring, fmt="%.6f %.6f")


def panels(cfg, tmp, layers, titles, stations, region, cpt, rng, annot, label, out, clip_hull=True, ncol=4,
           transp=35, panel_width=4.6):
    """Generic map panels: layers = list of (lon, lat, value) arrays on a regular grid."""
    W, E, S, N = region
    fs = cfg["plots"]["font_scale"]
    PH = G.mercator_height(panel_width, W, E, S, N)
    n = len(layers)
    ncol = min(ncol, n)
    nrow = (n + ncol - 1) // ncol
    gapx, gapy, x0, y0 = 0.25 * fs, 0.75 * fs, 1.6, 2.2
    pre = ["gmt makecpt -C%s -T%s -Z -D > %s/v.cpt" % (cpt, rng, tmp)]
    r_pre, relief = G.relief_cmds(cfg, tmp, W, E, S, N)
    pre += r_pre
    np.savetxt(tmp + "/sta.xy", np.array([(s["lon"], s["lat"]) for s in stations.values()]), fmt="%.6f")
    if clip_hull:
        _hull(stations, tmp)
    R = "-R%.5f/%.5f/%.5f/%.5f" % region
    dx, dy = _ticks(W, E, S, N)
    cmd = []
    for k, (lo, la, v) in enumerate(layers):
        ulo, ula = np.unique(lo), np.unique(la)
        ddx, ddy = np.diff(ulo).min(), np.diff(ula).min()
        np.savetxt("%s/val_%d.xyz" % (tmp, k), np.c_[lo, la, v], fmt="%.6f %.6f %.5f")
        pre += ["gmt xyz2grd %s/val_%d.xyz -R%.6f/%.6f/%.6f/%.6f -I%.7f/%.7f -G%s/a%d.nc" % (tmp, k, ulo[0], ulo[-1], ula[0], ula[-1], ddx, ddy, tmp, k),
                "gmt grdsample %s/a%d.nc -I%.7f -nl -G%s/s%d.nc" % (tmp, k, min(ddx, ddy) / 5, tmp, k)]
        col, row = k % ncol, k // ncol
        P = "%s -JM%gc -Xa%.3fc -Ya%.3fc" % (R, panel_width, x0 + col * (panel_width + gapx), y0 + (nrow - 1 - row) * (PH + gapy))
        if relief:
            cmd.append("gmt grdimage %s/topo.nc -I%s/shade.nc -C%s/relief.cpt %s" % (tmp, tmp, tmp, P))
        if clip_hull:
            cmd.append("gmt clip %s/hull.xy %s" % (tmp, P))
        cmd.append("gmt grdimage %s/s%d.nc -C%s/v.cpt -Q %s %s" % (tmp, k, tmp, "-t%d" % transp if relief else "", P))
        if clip_hull:
            cmd.append("gmt clip -C %s" % P)
        cmd += ["gmt plot %s/sta.xy -St%.2fc -Gwhite -W0.2p,black %s" % (tmp, 0.1 * fs, P),
                "gmt basemap %s -B%s%sne -Bxa%gf%g -Bya%gf%g" % (P, "W" if col == 0 else "w", "S" if row == nrow - 1 else "s",
                                                             dx, dx / 4, dy, dy / 4),
                "gmt text %s/title_%d.txt %s -F+cTC+f%gp,Helvetica+jBC -D0/0.12c -N" % (tmp, k, P, 10 * fs)]
        open("%s/title_%d.txt" % (tmp, k), "w").write(titles[k] + "\n")
        if k == n - 1:
            km, latc = G.scale_km(W, E, S, N)
            cmd.append("gmt basemap %s -LjBR+c%.4f+w%gk+o0.3c/0.35c+l'%g km' --MAP_SCALE_HEIGHT=0.1c --FONT_LABEL=%gp,Helvetica"
                       % (P, latc, km, km, 8 * fs))
    total = ncol * (panel_width + gapx) - gapx
    cmd.append("gmt colorbar -C%s/v.cpt -Dx%.3fc/%.3fc+w%.2fc/0.3c+h+jTC -Bx%s+l'%s' --FONT_ANNOT_PRIMARY=%gp,Helvetica "
               "--FONT_LABEL=%gp,Helvetica" % (tmp, x0 + total / 2, y0 - 0.9, min(9.0, total * 0.7), annot, label, 9 * fs, 10 * fs))
    return G.run(tmp, pre, cmd, out, cfg["plots"]["formats"], cfg["plots"]["dpi"], font=(7 * fs, 8 * fs))


def depth_slices(cfg, model, stations, out, depths=None):
    from ..inversion.run import read_model
    lon, lat, dep, v = read_model(model)
    depths = depths or cfg["plots"]["slice_depths_km"]
    use = [z for z in depths if np.min(np.abs(dep - z)) < 1e-3]
    if not use:
        raise RuntimeError("none of plots.slice_depths_km %s is a model depth node %s" % (depths, list(dep)))
    LO, LA = np.meshgrid(lon, lat, indexing="ij")
    layers = [(LO.ravel(), LA.ravel(), v[:, :, int(np.argmin(np.abs(dep - z)))].ravel()) for z in use]
    lo_, hi_ = cfg["plots"]["vs_range"]
    tmp = G.workdir()
    return panels(cfg, tmp, layers, ["%g km" % z for z in use], stations, (lon.min(), lon.max(), lat.min(), lat.max()),
                  cfg["plots"]["vs_cpt"], "%g/%g" % (lo_, hi_), "a0.2f0.1" if hi_ - lo_ <= 2 else "a0.5f0.1", "Vs (km/s)", out)


def checkerboard(cfg, project, stations, out, depths=None, vertical=False):
    """Input pattern (first panel) and the recovered model at several depths, in per cent of the 1-D background.
    vertical: the input changes with depth, so the top row shows the input and the bottom row the recovery at
    every depth."""
    syn = np.loadtxt(os.path.join(project, "DSurfTomo.inSyn.dat"))
    tru = np.loadtxt(os.path.join(project, "Vs_model.real"))
    M = open(os.path.join(project, "MOD")).read().split("\n")
    zmod = [round(float(v), 3) for v in M[0].split()]
    nrow = len([l for l in M[1:] if l.strip()]) // len(zmod)
    bg = {z: float(M[1 + k * nrow].split()[0]) for k, z in enumerate(zmod)}
    dep = np.unique(np.round(syn[:, 2], 3))
    depths = depths or [z for z in cfg["plots"]["slice_depths_km"] if np.min(np.abs(dep - z)) < 1e-3][:5]
    layers, titles = [], []

    def lay(a, z):
        m = np.round(a[:, 2], 3) == z
        return (a[m, 0], a[m, 1], 100 * (a[m, 3] - bg[z]) / bg[z])
    if vertical:
        for z in depths:
            layers.append(lay(tru, z))
            titles.append("Input %g km" % z)
    else:
        layers.append(lay(tru, depths[0]))
        titles.append("Input")
    for z in depths:
        layers.append(lay(syn, z))
        titles.append("Recovered %g km" % z)
    amp = 100 * cfg["checks"]["amplitude"] / np.mean(list(bg.values()))
    r = 5 * np.ceil(amp / 5)
    tmp = G.workdir()
    return panels(cfg, tmp, layers, titles, stations,
                  (syn[:, 0].min(), syn[:, 0].max(), syn[:, 1].min(), syn[:, 1].max()), G.CHECK_CPT, "-%g/%g" % (r, r),
                  "a%gf%g" % (r / 2, r / 10), "@~D@~Vs (%)", out, ncol=len(depths) if vertical else 3)
