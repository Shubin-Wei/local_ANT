"""Vertical sections of the model drawn with GMT, the usual form in ambient-noise papers: bilinear-resampled image with
a continuous colour table, thin labelled isolines (plots.contour_step), stations as inverted triangles on top, an
elevation profile above the section when the station table has elevations.  No masks (resolution is shown by the
checkerboards).

Areal arrays: plots.sections = [{name, start: [lon, lat], end: [lon, lat]}]; the model is interpolated
(trilinear) every 0.1 km along the line; stations within 0.5 km of the line are projected on it.
Lines: one section along the principal direction of the stations (from the south end, or the west end for an
east-west line), the velocity beneath each point of the line axis (Xie et al. 2024).
"""
import os

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from . import gmt as G
from ..signal import local_xy_km

CM_PER_KM = 1.7
TOPO_H = 1.1


def model_interp(model):
    from ..inversion.run import read_model
    lon, lat, dep, v = read_model(model)
    return lon, lat, dep, v, RegularGridInterpolator((lon, lat, dep), v, bounds_error=False, fill_value=np.nan)


def cut(f, dep, p0, p1, step=0.1, zmax=None):
    """Section from p0 to p1 (lon, lat): x (km), z, V[z, x], and the lon / lat of the samples."""
    x0, y0 = local_xy_km([p0[0], p1[0]], [p0[1], p1[1]])
    length = float(np.hypot(x0[1] - x0[0], y0[1] - y0[0]))
    x = np.arange(0, length + 1e-9, step)
    lon = p0[0] + (p1[0] - p0[0]) * x / max(length, 1e-9)
    lat = p0[1] + (p1[1] - p0[1]) * x / max(length, 1e-9)
    z = dep if zmax is None else dep[dep <= zmax + 1e-9]
    pts = np.array([(lo, la, zz) for zz in z for lo, la in zip(lon, lat)])
    return x, z, f(pts).reshape(len(z), len(x)), lon, lat


def line_axis(stations):
    lon = np.array([s["lon"] for s in stations.values()])
    lat = np.array([s["lat"] for s in stations.values()])
    x, y = local_xy_km(lon, lat)
    w, v = np.linalg.eigh(np.cov(np.c_[x, y].T))
    u = v[:, np.argmax(w)]
    if (abs(u[1]) >= abs(u[0]) and u[1] < 0) or (abs(u[1]) < abs(u[0]) and u[0] < 0):
        u = -u                                   # from the south end (north-south lines) or the west end
    s = x * u[0] + y * u[1]
    lat0, lon0 = float(np.mean(lat)), float(np.mean(lon))
    k = np.radians(1.0) * 6371.0
    def ll(sv):
        return lon0 + sv * u[0] / (k * np.cos(np.radians(lat0))), lat0 + sv * u[1] / k
    start, end = ll(s.min()), ll(s.max())
    return start, end, u, s - s.min(), (lon, lat)


def project_stations(stations, p0, p1, width_km):
    lon = np.array([s["lon"] for s in stations.values()])
    lat = np.array([s["lat"] for s in stations.values()])
    el = np.array([s["elev"] for s in stations.values()])
    lat0, lon0 = 0.5 * (p0[1] + p1[1]), 0.5 * (p0[0] + p1[0])
    xs, ys = local_xy_km(lon, lat, lat0, lon0)
    xa, ya = local_xy_km([p0[0], p1[0]], [p0[1], p1[1]], lat0, lon0)
    u = np.array([xa[1] - xa[0], ya[1] - ya[0]])
    L = np.hypot(*u)
    u = u / L
    s = (xs - xa[0]) * u[0] + (ys - ya[0]) * u[1]
    d = np.abs(-(xs - xa[0]) * u[1] + (ys - ya[0]) * u[0])
    m = (d <= width_km) & (s >= -0.05) & (s <= L + 0.05)
    return s[m], el[m]


def draw(cfg, rows, out, xlabel="Distance (km)", width=12.0):
    """rows: list of dict(x, z, v, sta, elev, label, zmax).  One section per row, top to bottom."""
    fs = cfg["plots"]["font_scale"]
    lo, hi = cfg["plots"]["vs_range"]
    step = cfg["plots"]["contour_step"]
    tmp = G.workdir()
    pre = ["gmt makecpt -C%s -T%g/%g -Z -D > %s/vs.cpt" % (cfg["plots"]["vs_cpt"], lo, hi, tmp)]
    cmd = []
    y = 2.2
    pos_rows = []
    for r in rows[::-1]:
        h = CM_PER_KM * r["zmax"]
        has_topo = r.get("elev") is not None and np.isfinite(r["elev"]).sum() >= 2
        pos_rows.append((r, y, h, has_topo))
        y += h + (TOPO_H + 0.45 if has_topo else 0) + 1.1 * fs
    for i, (r, y0, h, has_topo) in enumerate(pos_rows[::-1]):
        tag = "%s/r%d" % (tmp, i)
        x, z, v = np.asarray(r["x"], float), np.asarray(r["z"], float), np.asarray(r["v"], float)
        keep = z <= r["zmax"] + 1e-9
        z, v = z[keep], v[keep]
        V = v.copy()
        for k in range(V.shape[0]):                      # fill isolated gaps along x for display
            m = np.isfinite(V[k])
            if m.any() and not m.all():
                V[k, ~m] = np.interp(np.flatnonzero(~m), np.flatnonzero(m), V[k, m])
        xf = np.arange(x[0], x[-1] + 1e-9, 0.02)
        zf = np.arange(z[0], z[-1] + 1e-9, 0.01)
        fi = RegularGridInterpolator((z, x), V, bounds_error=False, fill_value=None)
        Z, X = np.meshgrid(zf, xf, indexing="ij")
        vf = fi(np.c_[Z.ravel(), X.ravel()])
        np.c_[X.ravel(), Z.ravel(), vf].astype("<f4").tofile(tag + ".bin")
        R = "-R%.4f/%.4f/%.4f/%.4f" % (xf[0], xf[-1], zf[0], zf[-1])
        J = "-JX%.3fc/-%.3fc -R%.4f/%.4f/0/%.4f" % (width, h, x[0], x[-1], r["zmax"])
        pos = "-Xa2.0c -Ya%.3fc" % y0
        cmd += ["gmt xyz2grd %s.bin -bi3f %s -I%d+n/%d+n -G%s.nc" % (tag, R, len(xf), len(zf), tag),
                "gmt grdimage %s.nc %s %s -C%s/vs.cpt" % (tag, J, pos, tmp)]
        if step:
            cmd.append("gmt grdcontour %s.nc %s %s -C%g -A%g+f%gp,Helvetica,gray10 -Gd3.5c -W0.2p,gray10 -S4 --FORMAT_FLOAT_MAP=%%.1f"
                       % (tag, J, pos, step, step, 5 * fs))
        if r.get("sta") is not None and len(r["sta"]):
            np.savetxt(tag + "_sta.txt", np.c_[r["sta"], np.zeros(len(r["sta"]))], fmt="%.4f")
            cmd.append("gmt plot %s_sta.txt %s %s -Si%.2fc -Gblack -D0/0.07c -N" % (tag, J, pos, 0.13 * fs))
        dx = G.nice_step(x[-1] - x[0], 5)
        dz = G.nice_step(r["zmax"], 3)
        cmd.append("gmt basemap %s %s -BWSen -Bxa%gf%g+l'%s' -Bya%gf%g+l'Depth (km)'" % (J, pos, dx, dx / 2, xlabel, dz, dz / 5))
        if has_topo:
            se, el = np.asarray(r["sta"]), np.asarray(r["elev"])
            ok = np.isfinite(el)
            o = np.argsort(se[ok])
            xs, e = se[ok][o], el[ok][o]
            elo, ehi = 50 * np.floor(e.min() / 50) - 50, 50 * np.ceil(e.max() / 50) + 50
            np.savetxt(tag + "_e.txt", np.c_[xs, e], fmt="%.4f")
            np.savetxt(tag + "_ef.txt", np.r_[np.c_[xs, e], [[xs[-1], elo], [xs[0], elo]]], fmt="%.4f")
            JT = "-JX%.3fc/%.3fc -R%.4f/%.4f/%g/%g" % (width, TOPO_H, x[0], x[-1], elo, ehi)
            pt = "-Xa2.0c -Ya%.3fc" % (y0 + h + 0.35)
            cmd += ["gmt plot %s_ef.txt %s %s -Ggray85" % (tag, JT, pt), "gmt plot %s_e.txt %s %s -W0.5p,gray30" % (tag, JT, pt),
                    "gmt basemap %s %s -BWsen -Bxf1 -Bya%gf%g+l'Elev. (m)'" % (JT, pt, G.nice_step(ehi - elo, 2), G.nice_step(ehi - elo, 2) / 2)]
            ylab = y0 + h + 0.35 + TOPO_H
        else:
            ylab = y0 + h + 0.25
        if r.get("label"):
            with open(tag + "_label.txt", "w") as fh:
                fh.write(r["label"] + "\n")
            cmd.append("gmt text %s_label.txt -R0/1/0/1 -JX%.3fc/1c -Xa2.0c -Ya%.3fc -F+cTL+f%gp,Helvetica-Bold+jBL -N"
                       % (tag, width, ylab - 0.75, 9 * fs))
    cmd.append("gmt colorbar -C%s/vs.cpt -Dx%.3fc/0.9c+w%.2fc/0.28c+h+jTC -Bxa0.2f0.1+l'Vs (km/s)' --FONT_ANNOT_PRIMARY=%gp,Helvetica "
               "--FONT_LABEL=%gp,Helvetica" % (tmp, 2.0 + width / 2, min(8.0, 0.7 * width), 8 * fs, 9 * fs))
    return G.run(tmp, pre, cmd, out, cfg["plots"]["formats"], cfg["plots"]["dpi"], font=(7 * fs, 8 * fs))


def model_zmax(cfg, project):
    import json
    meta = json.load(open(os.path.join(project, "project.json")))
    zs = meta["depths"]
    return cfg["plots"]["section_zmax_km"] or zs[-2]          # the deepest updated node


def surface_sections(cfg, project, stations, out):
    secs = cfg["plots"]["sections"]
    if not secs:
        return None
    lon, lat, dep, v, f = model_interp(os.path.join(project, "DSurfTomo.inMeasure.dat"))
    zmax = model_zmax(cfg, project)
    rows = []
    for s in secs:
        x, z, V, _, _ = cut(f, dep, s["start"], s["end"], zmax=zmax)
        sta, el = project_stations(stations, s["start"], s["end"], width_km=0.5)
        rows.append(dict(x=x, z=z, v=V, sta=sta, elev=el, label="%s-%s'" % (s["name"], s["name"]) if len(s["name"]) == 1 else s["name"], zmax=zmax))
    return draw(cfg, rows, out)


def line_section(cfg, project, stations, out, label=None):
    lon, lat, dep, v, f = model_interp(os.path.join(project, "DSurfTomo.inMeasure.dat"))
    zmax = model_zmax(cfg, project)
    p0, p1, u, s, _ = line_axis(stations)
    x, z, V, _, _ = cut(f, dep, p0, p1, zmax=zmax)
    el = np.array([st["elev"] for st in stations.values()])
    north_south = abs(u[1]) >= abs(u[0])
    rows = [dict(x=x, z=z, v=V, sta=s, elev=el, label=label or cfg["project"]["name"], zmax=zmax)]
    return draw(cfg, rows, out, xlabel="Distance along the line from the %s end (km)" % ("south" if north_south else "west"))
