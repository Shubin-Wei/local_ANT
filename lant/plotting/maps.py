"""Maps drawn with GMT: the station map and panels of model slices (Vs, checkerboards).

Model values are resampled bilinearly for display (the node values are not changed), clipped to the convex hull of
the stations, over a pale hill-shaded relief when one is available; GMT "seis" (red slow, blue fast) for Vs, red -
white - blue for perturbations.
"""

import numpy as np
from scipy.spatial import ConvexHull

from . import gmt as G


def _write(path, text):
    with open(path, "w") as f:
        f.write(text)


def station_map(cfg, stations, out, log=print):
    """stations {name: dict(lon, lat)}."""
    lon = [s["lon"] for s in stations.values()]
    lat = [s["lat"] for s in stations.values()]
    W, E, S, N = G.region(lon, lat)
    width = min(12.0, 14.0 / G.mercator_height(1.0, W, E, S, N))
    P = "-R%.5f/%.5f/%.5f/%.5f -JM%gc" % (W, E, S, N, width)
    dx, dy = G.nice_step(E - W, 3), G.nice_step(N - S, 4)
    km = G.scale_km(W, E, S, N)
    sym = 0.17 if len(stations) < 120 else 0.12

    def make(tmp, relief):
        np.savetxt(tmp + "/sta.xy", np.c_[lon, lat], fmt="%.6f")
        pre = G.relief_cmds(relief, tmp, (W, E, S, N), "geo") if relief else []
        cmd = ["gmt grdimage %s/topo.nc -I%s/shade.nc -C%s/relief.cpt %s" % (tmp, tmp, tmp, P)] if relief else []
        cmd += [
            "gmt plot %s/sta.xy %s -St%.2fc -Gwhite -W0.6p,black" % (tmp, P, sym),
            "gmt basemap %s -BWSne -Bxa%gf%g -Bya%gf%g %s" % (P, dx, dx / 4, dy, dy / 4, G.geo_format(dx, dy)),
            "gmt basemap %s -LjBL+c%.4f+w%gk+f+l'%g km'+o0.4c/0.4c" % (P, 0.5 * (S + N), km, km),
        ]
        _write(tmp + "/legend.txt", "S 0.14c t %.2fc white 0.6p,black 0.38c Stations (%d)\n" % (sym, len(stations)))
        cmd.append("gmt legend %s/legend.txt %s -DjTR+o0.15c -F+p0.3p,gray60+gwhite@10" % (tmp, P))
        if relief:
            cmd.append("gmt colorbar -C%s/relief.cpt %s -DJBC+w5c/0.22c+h+o0c/0.9c -Bxaf+l'Elevation (m)'" % (tmp, P))
        return pre, cmd

    return G.draw(cfg, make, out, log)


def hull_ring(lon, lat):
    xy = np.c_[lon, lat]
    h = ConvexHull(xy)
    return xy[np.r_[h.vertices, h.vertices[:1]]]


def panels(cfg, layers, labels, stations, cpt, rng, annot, label, out, ncol=4, show_stations=False, log=print):
    """Map panels of layers [(lon, lat, value) on a regular grid], clipped to the station hull, one colour bar."""
    slon = np.array([s["lon"] for s in stations.values()])
    slat = np.array([s["lat"] for s in stations.values()])
    W, E, S, N = G.region(slon, slat)
    pw = 4.6
    PH = G.mercator_height(pw, W, E, S, N)
    n = len(layers)
    ncol = min(ncol, n)
    nrow = (n + ncol - 1) // ncol
    gx, gy, x0, y0 = 0.8, 0.75, 1.6, 2.4
    R = "-R%.5f/%.5f/%.5f/%.5f" % (W, E, S, N)
    dx, dy = G.nice_step(E - W, 3), G.nice_step(N - S, 4)

    def make(tmp, relief):
        pre = ["gmt makecpt -C%s -T%s -Z -D > %s/v.cpt" % (cpt, rng, tmp)]
        pre += G.relief_cmds(relief, tmp, (W, E, S, N)) if relief else []
        np.savetxt(tmp + "/sta.xy", np.c_[slon, slat], fmt="%.6f")
        np.savetxt(tmp + "/hull.xy", hull_ring(slon, slat), fmt="%.6f")
        cmd = []
        for k, (lo, la, v) in enumerate(layers):
            # node spacing: the median step (model files keep 5 decimals, so single steps jitter by 1e-5)
            ulo, ula = np.unique(np.round(lo, 5)), np.unique(np.round(la, 5))
            ddx, ddy = float(np.median(np.diff(ulo))), float(np.median(np.diff(ula)))
            np.savetxt("%s/v%d.xyz" % (tmp, k), np.c_[lo, la, v], fmt="%.6f %.6f %.5f")
            pre += [
                "gmt xyz2grd %s/v%d.xyz -R%.6f/%.6f/%.6f/%.6f -I%.7f/%.7f -G%s/a%d.nc"
                % (tmp, k, ulo[0], ulo[-1], ula[0], ula[-1], ddx, ddy, tmp, k),
                "gmt grdsample %s/a%d.nc -I%.7f -nl -G%s/s%d.nc" % (tmp, k, min(ddx, ddy) / 5, tmp, k),
            ]
            col, row = k % ncol, k // ncol
            P = "%s -JM%gc -Xa%.3fc -Ya%.3fc" % (R, pw, x0 + col * (pw + gx), y0 + (nrow - 1 - row) * (PH + gy))
            if relief:
                cmd.append("gmt grdimage %s/topo.nc -I%s/shade.nc -C%s/relief.cpt %s" % (tmp, tmp, tmp, P))
            cmd += ["gmt clip %s/hull.xy %s" % (tmp, P),
                    "gmt grdimage %s/s%d.nc -C%s/v.cpt -Q %s %s" % (tmp, k, tmp, "-t30" if relief else "", P),
                    "gmt clip -C %s" % P]  # fmt: skip
            if show_stations:
                cmd.append("gmt plot %s/sta.xy -St0.12c -W0.3p,black %s" % (tmp, P))
            cmd.append("gmt basemap %s -B%s%sne -Bxa%gf%g -Bya%gf%g %s" % (
                P, "W" if col == 0 else "w", "S" if row == nrow - 1 else "s", dx, dx / 4, dy, dy / 4,
                G.geo_format(dx, dy)))  # fmt: skip
            _write("%s/label%d.txt" % (tmp, k), labels[k] + "\n")
            cmd.append("gmt text %s/label%d.txt %s -F+cTL+f9p,Helvetica-Bold+jTL -D0.08c/-0.08c -Gwhite -W0.4p "
                       "-C0.05c/0.04c" % (tmp, k, P))  # fmt: skip
            if k == n - 1:
                km = G.scale_km(W, E, S, N)
                cmd.append("gmt basemap %s -LjBR+c%.4f+w%gk+o0.3c/0.35c+l'%g km' --MAP_SCALE_HEIGHT=0.1c"
                           % (P, 0.5 * (S + N), km, km))  # fmt: skip
        total = ncol * (pw + gx) - gx
        cmd.append("gmt colorbar -C%s/v.cpt -Dx%.3fc/%.3fc+w%.2fc/0.3c+h+jTC+e -Bx%s+l'%s'"
                   % (tmp, x0 + total / 2, y0 - 0.9, min(9.0, 0.7 * total), annot, label))  # fmt: skip
        return pre, cmd

    return G.draw(cfg, make, out, log)


def colorbar_ticks(lo, hi):
    from ..geo import nice_number

    a = nice_number((hi - lo) / 5, "nearest")
    return "a%gf%g" % (a, a / 5)
