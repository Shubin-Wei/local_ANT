"""Maps and sections drawn with GMT: the station map, map panels (Vs slices, checkerboards, ray paths) and vertical
sections of the model.

Model values are resampled bilinearly for display (the node values are not changed) and clipped to the convex hull
of the stations, over a pale hill-shaded relief when one is available; GMT "seis" (red slow, blue fast) for Vs, red -
white - blue for perturbations.
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import ConvexHull

from .. import geo
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
    sym = 0.2 if len(stations) < 120 else 0.15

    def make(tmp, relief):
        np.savetxt(tmp + "/sta.xy", np.c_[lon, lat], fmt="%.6f")
        pre = G.relief_cmds(relief, tmp, (W, E, S, N), "geo") if relief else []
        cmd = ["gmt grdimage %s/topo.nc -I%s/shade.nc -C%s/relief.cpt %s" % (tmp, tmp, tmp, P)] if relief else []
        cmd += [
            "gmt plot %s/sta.xy %s -St%.2fc -Gwhite -W0.6p,black" % (tmp, P, sym),
            "gmt basemap %s -BWSne -Bxa%gf%g -Bya%gf%g %s" % (P, dx, dx / 4, dy, dy / 4, G.geo_format(dx, dy)),
            "gmt basemap %s -LjBL+c%.4f+w%gk+f+l'%g km'+o0.4c/0.4c -F+gwhite@30+c0.1c" % (P, 0.5 * (S + N), km, km),
        ]
        _write(tmp + "/legend.txt", "S 0.16c t %.2fc white 0.6p,black 0.42c Stations (%d)\n" % (sym, len(stations)))
        cmd.append("gmt legend %s/legend.txt %s -DjTR+o0.15c -F+p0.3p,gray60+gwhite@10" % (tmp, P))
        if relief:
            cmd.append("gmt colorbar -C%s/relief.cpt %s -DJBC+w6c/0.25c+h+o0c/1.0c -Bxaf+l'Elevation (m)'" % (tmp, P))
        return pre, cmd

    return G.draw(cfg, make, out, log)


def hull_ring(lon, lat):
    xy = np.c_[lon, lat]
    h = ConvexHull(xy)
    return xy[np.r_[h.vertices, h.vertices[:1]]]


def panels(cfg, layers, labels, stations, out, cpt=None, rng=None, annot=None, label=None, ncol=4,
           show_stations=False, overlay=None, log=print):  # fmt: skip
    """Map panels with the label above each; layers [(lon, lat, value) on a regular grid, or None for no image],
    images clipped to the station hull, one colour bar.  overlay(k, tmp, P) -> GMT commands drawn on panel k."""
    slon = np.array([s["lon"] for s in stations.values()])
    slat = np.array([s["lat"] for s in stations.values()])
    W, E, S, N = G.region(slon, slat)
    pw = 5.0
    PH = G.mercator_height(pw, W, E, S, N)
    n = len(labels)
    ncol = min(ncol, n)
    nrow = (n + ncol - 1) // ncol
    gx, gy, x0 = 0.45, 1.0, 1.8
    y0 = 2.6 if cpt else 1.2
    R = "-R%.5f/%.5f/%.5f/%.5f" % (W, E, S, N)
    dx, dy = G.nice_step(E - W, 3), G.nice_step(N - S, 4)

    def make(tmp, relief):
        pre = ["gmt makecpt -C%s -T%s -Z -D > %s/v.cpt" % (cpt, rng, tmp)] if cpt else []
        pre += G.relief_cmds(relief, tmp, (W, E, S, N)) if relief else []
        np.savetxt(tmp + "/sta.xy", np.c_[slon, slat], fmt="%.6f")
        np.savetxt(tmp + "/hull.xy", hull_ring(slon, slat), fmt="%.6f")
        cmd = []
        for k in range(n):
            col, row = k % ncol, k // ncol
            P = "%s -JM%gc -Xa%.3fc -Ya%.3fc" % (R, pw, x0 + col * (pw + gx), y0 + (nrow - 1 - row) * (PH + gy))
            if relief:
                cmd.append("gmt grdimage %s/topo.nc -I%s/shade.nc -C%s/relief.cpt %s" % (tmp, tmp, tmp, P))
            if layers[k] is not None:
                lo, la, v = layers[k]
                # node spacing: the median step (model files keep 5 decimals, so single steps jitter by 1e-5)
                ulo, ula = np.unique(np.round(lo, 5)), np.unique(np.round(la, 5))
                ddx, ddy = float(np.median(np.diff(ulo))), float(np.median(np.diff(ula)))
                np.savetxt("%s/v%d.xyz" % (tmp, k), np.c_[lo, la, v], fmt="%.6f %.6f %.5f")
                pre += [
                    "gmt xyz2grd %s/v%d.xyz -R%.6f/%.6f/%.6f/%.6f -I%.7f/%.7f -G%s/a%d.nc"
                    % (tmp, k, ulo[0], ulo[-1], ula[0], ula[-1], ddx, ddy, tmp, k),
                    "gmt grdsample %s/a%d.nc -I%.7f -nl -G%s/s%d.nc" % (tmp, k, min(ddx, ddy) / 5, tmp, k),
                ]
                cmd += ["gmt clip %s/hull.xy %s" % (tmp, P),
                        "gmt grdimage %s/s%d.nc -C%s/v.cpt -Q %s %s" % (tmp, k, tmp, "-t25" if relief else "", P),
                        "gmt clip -C %s" % P]  # fmt: skip
            if overlay is not None:
                cmd += overlay(k, tmp, P)
            if show_stations:
                cmd.append("gmt plot %s/sta.xy -St0.13c -W0.35p,black %s" % (tmp, P))
            cmd.append("gmt basemap %s -B%s%sne -Bxa%gf%g -Bya%gf%g %s" % (
                P, "W" if col == 0 else "w", "S" if row == nrow - 1 else "s", dx, dx / 4, dy, dy / 4,
                G.geo_format(dx, dy)))  # fmt: skip
            _write("%s/label%d.txt" % (tmp, k), labels[k] + "\n")
            cmd.append("gmt text %s/label%d.txt %s -F+cTC+f11p,Helvetica+jBC -D0/0.15c -N" % (tmp, k, P))
        total = ncol * (pw + gx) - gx
        km = G.scale_km(W, E, S, N)  # scale bar below the panels, at the right
        Ps = "%s -JM%gc -Xa%.3fc -Ya%.3fc" % (R, pw, x0 + total - pw, y0 - PH - (1.0 if cpt else 0.6))
        cmd.append("gmt basemap %s -LjTR+c%.4f+w%gk+l'%g km' --MAP_SCALE_HEIGHT=0.12c" % (Ps, 0.5 * (S + N), km, km))
        if cpt:
            cmd.append("gmt colorbar -C%s/v.cpt -Dx%.3fc/%.3fc+w%.2fc/0.32c+h+jTC+e -Bx%s+l'%s' "
                       "--FONT_ANNOT_PRIMARY=10p --FONT_LABEL=11p" % (tmp, x0 + total / 2, y0 - 1.0, min(10.0, 0.6 * total),
                                                                      annot, label))  # fmt: skip
        return pre, cmd

    return G.draw(cfg, make, out, log)


def line_overlay(lines):
    """Overlay for panels: the section lines [(name, (lon0, lat0), (lon1, lat1))] dashed, ends labelled."""

    def cmds(k, tmp, P):
        out = []
        for i, (nm, a, b) in enumerate(lines):
            f = "%s/line%d.xy" % (tmp, i)
            np.savetxt(f, np.array([a, b]), fmt="%.6f")
            out.append("gmt plot %s %s -W1.0p,black,-" % (f, P))
            _write("%s/end%d.txt" % (tmp, i), "%.6f %.6f %s\n%.6f %.6f %s'\n" % (a[0], a[1], nm, b[0], b[1], nm))
            out.append("gmt text %s/end%d.txt %s -F+f10p,Helvetica-Bold,black=2p,white+jCM -N" % (tmp, i, P))
            out.append("gmt text %s/end%d.txt %s -F+f10p,Helvetica-Bold,black+jCM -N" % (tmp, i, P))
        return out

    return cmds


def section_lines(stations):
    """Two sections through the centre of the stations along the principal directions of their distribution, from
    the first to the last station projected on each -> [(name, (lon, lat) start, (lon, lat) end)]."""
    lon = np.array([s["lon"] for s in stations.values()])
    lat = np.array([s["lat"] for s in stations.values()])
    lon0, lat0 = float(lon.mean()), float(lat.mean())
    kx, ky = geo.km_per_degree(lat0)
    x, y = (lon - lon0) * kx, (lat - lat0) * ky
    w, v = np.linalg.eigh(np.cov(np.c_[x, y].T))
    out = []
    for name, u in (("A", v[:, 1]), ("B", v[:, 0])):
        if u[0] < 0:  # west to east (a north-south line south to north)
            u = -u
        if abs(u[0]) < 1e-6 and u[1] < 0:
            u = -u
        p = x * u[0] + y * u[1]
        s0, s1 = float(p.min()), float(p.max())
        out.append(
            (name, (lon0 + s0 * u[0] / kx, lat0 + s0 * u[1] / ky), (lon0 + s1 * u[0] / kx, lat0 + s1 * u[1] / ky))
        )
    return out


def sections(cfg, model, stations, lines, vr, out, zmax=None, cpt="seis", label="Vs (km/s)", contours=True,
             log=print):  # fmt: skip
    """Vertical sections of model (lon, lat, depth, V [nlon, nlat, nz]) along lines (section_lines), Vs in colour with
    contours, the relief along the line above each when available, the stations within half a station spacing of the
    line as triangles; no vertical exaggeration unless the section would be flatter than 1:6."""
    lon, lat, dep, V = model
    f = RegularGridInterpolator((lon, lat, dep), V, bounds_error=False, fill_value=np.nan)
    zmax = float(dep[-1]) if zmax is None else float(zmax)
    slon = np.array([s["lon"] for s in stations.values()])
    slat = np.array([s["lat"] for s in stations.values()])
    lat0 = float(np.mean(slat))
    kx, ky = geo.km_per_degree(lat0)
    spacing = float(np.median(geo.nearest_neighbour_km(slon, slat)))
    lens = [float(np.hypot((b[0] - a[0]) * kx, (b[1] - a[1]) * ky)) for _, a, b in lines]
    cm = 16.0 / max(lens)
    ve = max(1.0, float(np.ceil(max(lens) / (6 * zmax))))
    topo_h, gap = 1.3, 2.3
    from ..geo import nice_number

    cstep = nice_number((vr[1] - vr[0]) / 8, "nearest")
    zstep = nice_number(zmax / 4, "nearest")

    def make(tmp, relief):
        pre = ["gmt makecpt -C%s -T%g/%g -Z -D > %s/vs.cpt" % (cpt, vr[0], vr[1], tmp)]
        if relief:
            W, E, S, N = G.region(slon, slat)
            pre += G.relief_cmds(relief, tmp, (W, E, S, N), "geo")[:1]
        cmd = ["gmt set FONT_ANNOT_PRIMARY 10p,Helvetica FONT_LABEL 11p,Helvetica"]
        y = 1.6
        heights = []
        for i, ((nm, a, b), L) in enumerate(reversed(list(zip(lines, lens)))):
            w, hgt = L * cm, zmax * cm * ve
            s = np.linspace(0, L, max(60, int(L / 0.02)))
            t = s / L
            lo, la = a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])
            z = np.linspace(0, zmax, max(40, int(zmax / 0.01)))
            S_, Z_ = np.meshgrid(s, z)
            v = f(np.c_[np.tile(lo, z.size), np.tile(la, z.size), Z_.ravel()]).reshape(Z_.shape)
            tag = "%s/s%d" % (tmp, i)
            np.c_[S_.ravel(), Z_.ravel(), v.ravel()].astype("<f4").tofile(tag + ".bin")
            P = "-JX%.3fc/-%.3fc -R0/%.4f/0/%g -Xa%.3fc -Ya%.3fc" % (w, hgt, L, zmax, 2.0, y)
            pre.append("gmt xyz2grd %s.bin -bi3f -R0/%.4f/0/%g -I%d+n/%d+n -G%s.nc" % (tag, s[-1], zmax, s.size, z.size,
                                                                                       tag))  # fmt: skip
            cmd.append("gmt grdimage %s.nc %s -C%s/vs.cpt -Q" % (tag, P, tmp))
            if contours:
                cmd.append("gmt grdcontour %s.nc %s -C%g -A%g+f9p,Helvetica,gray10 -Gd6c -W0.4p,gray25 "
                           "--FORMAT_FLOAT_MAP=%%.4g" % (tag, P, cstep, 2 * cstep))  # fmt: skip
            # stations near the line
            x, yy = (slon - a[0]) * kx, (slat - a[1]) * ky
            ux, uy = (b[0] - a[0]) * kx / L, (b[1] - a[1]) * ky / L
            along, across = x * ux + yy * uy, np.abs(-x * uy + yy * ux)
            near = (across <= 0.5 * spacing) & (along >= 0) & (along <= L)
            if near.any():
                np.savetxt(tag + "_sta.txt", np.c_[along[near], np.zeros(near.sum())], fmt="%.4f")
                cmd.append("gmt plot %s_sta.txt %s -Si0.25c -Gwhite -W0.5p,black -D0/0.13c -N" % (tag, P))
            xl = "+l'Distance along the section (km)'" if i == 0 else ""
            cmd.append("gmt basemap %s -BWSen -Bxaf%s -Bya%gf%g+l'Depth (km)'" % (P, xl, zstep, zstep / 2))
            top = y + hgt
            if relief:
                np.savetxt(tag + "_ll.txt", np.c_[lo, la], fmt="%.6f")
                pre.append("gmt grdtrack %s_ll.txt -G%s/topo.nc > %s_el.txt" % (tag, tmp, tag))
                Pt = "-JX%.3fc/%.3fc -R0/%.4f/${ZL%d}/${ZH%d} -Xa%.3fc -Ya%.3fc" % (w, topo_h, L, i, i, 2.0, top + 0.12)
                pre.append("read ZL%d ZH%d <<< \"$(awk '{print $3}' %s_el.txt | sort -g | awk 'NR==1{a=$1} {b=$1} END"
                           "{d=b-a; if (d<10) d=10; printf \"%%d %%d\", a-0.15*d-1, b+0.15*d+1}')\"" % (i, i, tag))  # fmt: skip
                pre.append("awk -v s=%.6f -v n=%d '{printf \"%%.4f %%s\\n\", (NR-1)*s/(n-1), $3}' %s_el.txt > %s_p.txt"
                           % (L, s.size, tag, tag))  # fmt: skip
                pre.append('(cat %s_p.txt; echo "%.4f ${ZL%d}"; echo "0 ${ZL%d}") > %s_f.txt' % (tag, L, i, i, tag))
                cmd += ["gmt plot %s_f.txt %s -Ggray85" % (tag, Pt), "gmt plot %s_p.txt %s -W0.8p,gray30" % (tag, Pt),
                        "gmt basemap %s -BWsen -Bxf -Byaf+l'm'" % Pt]  # fmt: skip
                top += topo_h + 0.12
            _write(tag + "_ends.txt", "0 0 %s\n%.4f 0 %s'\n" % (nm, L, nm))
            Pn = "-JX%.3fc/1c -R0/%.4f/0/1 -Xa%.3fc -Ya%.3fc" % (w, L, 2.0, top)
            cmd.append("gmt text %s_ends.txt %s -F+f13p,Helvetica-Bold+jBC -D0/0.12c -N" % (tag, Pn))
            heights.append((y, hgt, w))
            y = top + gap
        yb, hb, wb = heights[-1]
        cmd.append("gmt colorbar -C%s/vs.cpt -Dx%.3fc/%.3fc+w%.2fc/0.35c+jBL+e -Bxa%gf%g -By+l'%s' "
                   "--FONT_ANNOT_PRIMARY=10p --FONT_LABEL=11p"
                   % (tmp, 2.0 + max(w for _, _, w in heights) + 0.5, heights[0][0], min(6.0, sum(h for _, h, _ in
                      heights) + 1.0), 2 * cstep, cstep, label))  # fmt: skip
        if ve > 1:
            _write(tmp + "/ve.txt", "V.E. %g:1\n" % ve)
            cmd.append("gmt text %s/ve.txt -JX%gc/1c -R0/1/0/1 -Xa%.3fc -Ya%.3fc -F+cBR+f10p,Helvetica+jBR -N"
                       % (tmp, heights[0][2], 2.0, 0.3))  # fmt: skip
        return pre, cmd

    return G.draw(cfg, make, out, log)


def colorbar_ticks(lo, hi):
    from ..geo import nice_number

    a = nice_number((hi - lo) / 5, "nearest")
    return "a%gf%g" % (a, a / 5)
