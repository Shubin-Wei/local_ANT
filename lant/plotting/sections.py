"""Vertical sections of the model drawn with GMT, the usual form in ambient-noise papers: bilinear-resampled image with
a continuous colour table, thin labelled isolines (plots.contour_step), stations as inverted triangles on top, an
elevation profile above the section when the station table has elevations (areal arrays: interpolated linearly
between all stations, inside their convex hull; lines: the station elevations).  No masks (resolution is shown by the
checkerboards).

Areal arrays: plots.sections = [{name, start: [lon, lat], end: [lon, lat]}]; the model is interpolated
(trilinear) every 0.1 km along the line; stations within 0.5 km of the line are projected on it.
Lines: one section along the principal direction of the stations (from the south end, or the west end for an
east-west line), the velocity beneath each point of the line axis; the checkerboards of a line are drawn the same way
(input above, recovery below, in per cent of the 1-D background), since the model of a line is a narrow strip.
"""

import os

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from . import gmt as G
from ..signal import km_per_degree, local_xy_km

H_MIN, H_MAX = 2.5, 8.0  # height (cm) of a section panel: true scale inside these limits


def panel_height(length_km, zmax_km, width_cm):
    """Height (cm) of a section of length_km x zmax_km drawn width_cm wide, and its vertical exaggeration: true
    scale (1) when that gives a panel between H_MIN and H_MAX, otherwise the nearest of the two limits."""
    true = width_cm * zmax_km / max(length_km, 1e-9)
    h = min(max(true, H_MIN), H_MAX)
    return h, h / true


TOPO_H = 1.1


def model_interp(model):
    from ..inversion.run import read_model

    lon, lat, dep, v = read_model(model)
    return lon, lat, dep, v, RegularGridInterpolator((lon, lat, dep), v, bounds_error=False, fill_value=np.nan)


def cut(f, dep, p0, p1, step=None, zmax=None):
    """Section from p0 to p1 (lon, lat): x (km), z, V[z, x], and the lon / lat of the samples (300 samples along
    the section unless a step in km is given)."""
    x0, y0 = local_xy_km([p0[0], p1[0]], [p0[1], p1[1]])
    length = float(np.hypot(x0[1] - x0[0], y0[1] - y0[0]))
    x = np.arange(0, length + 1e-9, step or length / 300.0)
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
        u = -u  # from the south end (north-south lines) or the west end
    s = x * u[0] + y * u[1]
    lat0, lon0 = float(np.mean(lat)), float(np.mean(lon))
    kx, ky = km_per_degree(lat0)

    def ll(sv):
        return lon0 + sv * u[0] / kx, lat0 + sv * u[1] / ky

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
    width_km = 0.05 * L if width_km is None else width_km  # stations within 5 % of the section length
    m = (d <= width_km) & (s >= -0.005 * L) & (s <= 1.005 * L)
    return s[m], el[m]


def draw(cfg, rows, out, xlabel="Distance (km)", width=12.0, scale=None):
    """rows: list of dict(x, z, v, sta, elev, label, zmax[, topo]).  One section per row, top to bottom.  The elevation
    profile above a section is topo = (x, elevation) when given, otherwise the elevations of the stations drawn on it.
    scale: dict(cpt, lo, hi, label) for values other than Vs (no isolines then); default the Vs colour table."""
    fs = cfg["plots"]["font_scale"]
    if scale:
        cpt, (lo, hi), cblabel, step = scale["cpt"], (scale["lo"], scale["hi"]), scale["label"], None
    else:
        cpt, (lo, hi), cblabel, step = (
            cfg["plots"]["vs_cpt"],
            cfg["plots"]["vs_range"],
            "Vs (km/s)",
            cfg["plots"]["contour_step"],
        )
    tmp = G.workdir()
    pre = ["gmt makecpt -C%s -T%g/%g -Z -D > %s/vs.cpt" % (cpt, lo, hi, tmp)]
    cmd = []
    y = 2.2
    pos_rows = []
    for r in rows[::-1]:
        h, ve = panel_height(float(r["x"][-1] - r["x"][0]), float(r["zmax"]), width)
        r["ve"] = ve
        if r.get("topo") is None and r.get("elev") is not None:
            r["topo"] = (np.asarray(r["sta"], float), np.asarray(r["elev"], float))
        has_topo = r.get("topo") is not None and np.isfinite(r["topo"][1]).sum() >= 2
        pos_rows.append((r, y, h, has_topo))
        y += h + (TOPO_H + 0.45 if has_topo else 0) + 2.1 * fs  # x label of this row, head line of the next
    for i, (r, y0, h, has_topo) in enumerate(pos_rows[::-1]):
        tag = "%s/r%d" % (tmp, i)
        x, z, v = np.asarray(r["x"], float), np.asarray(r["z"], float), np.asarray(r["v"], float)
        keep = z <= r["zmax"] + 1e-9
        z, v = z[keep], v[keep]
        V = v.copy()
        for k in range(V.shape[0]):  # fill isolated gaps along x for display
            m = np.isfinite(V[k])
            if m.any() and not m.all():
                V[k, ~m] = np.interp(np.flatnonzero(~m), np.flatnonzero(m), V[k, m])
        xf = np.linspace(x[0], x[-1], 601)  # smooth display at any scale of the model
        zf = np.linspace(z[0], z[-1], 201)
        fi = RegularGridInterpolator((z, x), V, bounds_error=False, fill_value=None)
        Z, X = np.meshgrid(zf, xf, indexing="ij")
        vf = fi(np.c_[Z.ravel(), X.ravel()])
        np.c_[X.ravel(), Z.ravel(), vf].astype("<f4").tofile(tag + ".bin")
        R = "-R%.4f/%.4f/%.4f/%.4f" % (xf[0], xf[-1], zf[0], zf[-1])
        J = "-JX%.3fc/-%.3fc -R%.4f/%.4f/0/%.4f" % (width, h, x[0], x[-1], r["zmax"])
        pos = "-Xa2.0c -Ya%.3fc" % y0
        cmd += [
            "gmt xyz2grd %s.bin -bi3f %s -I%d+n/%d+n -G%s.nc" % (tag, R, len(xf), len(zf), tag),
            "gmt grdimage %s.nc %s %s -C%s/vs.cpt" % (tag, J, pos, tmp),
        ]
        if step:
            nd = max(1, len(("%.6f" % step).rstrip("0").split(".")[1]))  # 0.25 km/s: labels 1.25, 1.50
            cmd.append(
                "gmt grdcontour %s.nc %s %s -C%g -A%g+f%gp,Helvetica,gray10 -Gd3.5c -W0.2p,gray10 -S4 --FORMAT_FLOAT_MAP=%%.%df"
                % (tag, J, pos, step, step, 5 * fs, nd)
            )
        if r.get("sta") is not None and len(r["sta"]):
            np.savetxt(tag + "_sta.txt", np.c_[r["sta"], np.zeros(len(r["sta"]))], fmt="%.4f")
            cmd.append("gmt plot %s_sta.txt %s %s -Si%.2fc -Gblack -D0/0.07c -N" % (tag, J, pos, 0.13 * fs))
        dx = G.nice_step(x[-1] - x[0], 5)
        dz = G.nice_step(r["zmax"], 3)
        cmd.append(
            "gmt basemap %s %s -BWSen -Bxa%gf%g+l'%s' -Bya%gf%g+l'Depth (km)'"
            % (J, pos, dx, dx / 2, xlabel, dz, dz / 5)
        )
        if has_topo:
            se, el = np.asarray(r["topo"][0], float), np.asarray(r["topo"][1], float)
            ok = np.isfinite(el)
            o = np.argsort(se[ok])
            xs, e = se[ok][o], el[ok][o]
            elo, ehi = 50 * np.floor(e.min() / 50) - 50, 50 * np.ceil(e.max() / 50) + 50
            np.savetxt(tag + "_e.txt", np.c_[xs, e], fmt="%.4f")
            np.savetxt(tag + "_ef.txt", np.r_[np.c_[xs, e], [[xs[-1], elo], [xs[0], elo]]], fmt="%.4f")
            JT = "-JX%.3fc/%.3fc -R%.4f/%.4f/%g/%g" % (width, TOPO_H, x[0], x[-1], elo, ehi)
            pt = "-Xa2.0c -Ya%.3fc" % (y0 + h + 0.35)
            cmd += [
                "gmt plot %s_ef.txt %s %s -Ggray85" % (tag, JT, pt),
                "gmt plot %s_e.txt %s %s -W0.5p,gray30" % (tag, JT, pt),
                "gmt basemap %s %s -BWsen -Bxf1 -Bya%gf%g+l'Elev. (m)'"
                % (JT, pt, G.nice_step(ehi - elo, 3), G.nice_step(ehi - elo, 3) / 2),
            ]
            ylab = y0 + h + 0.35 + TOPO_H
        else:
            ylab = y0 + h + 0.25
        # one line above the section (and its elevation panel): the ends (letter and direction) at the two sides,
        # a name and the vertical exaggeration in the middle
        ve = r["ve"]
        mid = " ".join(
            t for t in (r.get("label") or "", "" if abs(ve - 1) < 0.05 else "(vertical exaggeration %.2g)" % ve) if t
        )
        with open(tag + "_head.txt", "w") as fh:
            if r.get("ends"):
                fh.write(
                    "0 0 %gp,Helvetica-Bold BL %s\n1 0 %gp,Helvetica-Bold BR %s\n"
                    % (9 * fs, r["ends"][0], 9 * fs, r["ends"][1])
                )
            if mid:
                fh.write("0.5 0 %gp,Helvetica BC %s\n" % (8 * fs, mid))
        cmd.append("gmt text %s_head.txt -R0/1/0/1 -JX%.3fc/1c -Xa2.0c -Ya%.3fc -F+f+j -N" % (tag, width, ylab + 0.12))
    cmd.append(
        "gmt colorbar -C%s/vs.cpt -Dx%.3fc/0.9c+w%.2fc/0.28c+h+jTC -Bxaf+l'%s' --FONT_ANNOT_PRIMARY=%gp,Helvetica "
        "--FONT_LABEL=%gp,Helvetica" % (tmp, 2.0 + width / 2, min(8.0, 0.7 * width), cblabel, 8 * fs, 9 * fs)
    )
    return G.run(tmp, pre, cmd, out, cfg["plots"]["formats"], cfg["plots"]["dpi"], font=(7 * fs, 8 * fs))


COMPASS = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]


def compass(p0, p1):
    """Compass names of the start and the end of a section from p0 to p1 (lon, lat): the direction the start lies
    in, seen from the middle, and the direction of the end."""
    x, y = local_xy_km([p0[0], p1[0]], [p0[1], p1[1]])
    az = np.degrees(np.arctan2(x[1] - x[0], y[1] - y[0])) % 360
    name = lambda a: COMPASS[int(((a + 11.25) % 360) // 22.5)]
    return name((az + 180) % 360), name(az)


def relief_along(cfg, lon, lat):
    """Elevation (m) of the relief grid (plots.topography) at the points of a section, or None without relief."""
    import subprocess
    import tempfile
    from . import relief

    W, E, S, N = min(lon) - 0.01, max(lon) + 0.01, min(lat) - 0.01, max(lat) + 0.01
    tmp = tempfile.mkdtemp(prefix="lant_track_")
    try:
        grid, _ = relief.grid(cfg, W, E, S, N, tmp)
        if not grid:
            return None
        pts = os.path.join(tmp, "p.xy")
        np.savetxt(pts, np.c_[lon, lat], fmt="%.6f")
        r = subprocess.run(
            ["gmt", "grdtrack", pts, "-G" + grid], capture_output=True, text=True, cwd=tmp
        )  # gmt.history there
        z = np.array([float(l.split()[2]) for l in r.stdout.split("\n") if len(l.split()) >= 3])
        return z if len(z) == len(lon) else None
    finally:
        import shutil

        shutil.rmtree(tmp, ignore_errors=True)


def model_zmax(cfg, project):
    import json

    meta = json.load(open(os.path.join(project, "project.json")))
    zs = meta["depths"]
    return cfg["plots"]["section_zmax_km"] or zs[-2]  # the deepest updated node


def elevation_along(stations, lon, lat):
    """Elevation (m) at the points of a section, interpolated linearly between all stations (Delaunay triangles of
    the station positions); NaN outside the convex hull of the stations."""
    from scipy.interpolate import LinearNDInterpolator

    slon = np.array([s["lon"] for s in stations.values()])
    slat = np.array([s["lat"] for s in stations.values()])
    el = np.array([s["elev"] for s in stations.values()], float)
    ok = np.isfinite(el)
    if ok.sum() < 3:
        return np.full(len(lon), np.nan)
    lat0, lon0 = float(np.mean(slat)), float(np.mean(slon))
    xs, ys = local_xy_km(slon[ok], slat[ok], lat0, lon0)
    xp, yp = local_xy_km(np.asarray(lon), np.asarray(lat), lat0, lon0)
    try:
        return LinearNDInterpolator(np.c_[xs, ys], el[ok])(np.c_[xp, yp])
    except Exception:  # stations on one line: no triangles
        return np.full(len(lon), np.nan)


def surface_sections(cfg, project, stations, out):
    secs = cfg["plots"]["sections"]
    if not secs:
        return None
    lon, lat, dep, v, f = model_interp(os.path.join(project, "DSurfTomo.inMeasure.dat"))
    zmax = model_zmax(cfg, project)
    rows = []
    for s in secs:
        x, z, V, slon, slat = cut(f, dep, s["start"], s["end"], zmax=zmax)
        sta, el = project_stations(stations, s["start"], s["end"], width_km=None)
        topo = relief_along(cfg, slon, slat)
        if topo is None:
            topo = elevation_along(stations, slon, slat)
        inside = np.where(np.isfinite(V).any(axis=0))[0]  # a section drawn beyond the model is cut to it
        if not inside.size:
            print("[lant] section %s lies outside the model: not drawn" % s["name"], flush=True)
            continue
        a, b = inside[0], inside[-1] + 1
        if a > 0 or b < len(x):
            print(
                "[lant] section %s: cut to the model, %.2f km left out at its start and %.2f km at its end"
                % (s["name"], x[a] - x[0], x[-1] - x[b - 1]),
                flush=True,
            )
            keep = (sta >= x[a]) & (sta <= x[b - 1])
            sta, el = sta[keep] - x[a], el[keep]
            x, V, topo = x[a:b] - x[a], V[:, a:b], topo[a:b]
        c0, c1 = compass(s["start"], s["end"])
        one = len(s["name"]) == 1
        rows.append(
            dict(
                x=x,
                z=z,
                v=V,
                sta=sta,
                elev=el,
                topo=(x, topo),
                zmax=zmax,
                label="" if one else s["name"],
                ends=("%s  %s" % (s["name"], c0), "%s  %s'" % (c1, s["name"])) if one else (c0, c1),
            )
        )
    return draw(cfg, rows, out) if rows else None


def line_section(cfg, project, stations, out, label=None):
    lon, lat, dep, v, f = model_interp(os.path.join(project, "DSurfTomo.inMeasure.dat"))
    zmax = model_zmax(cfg, project)
    p0, p1, u, s, _ = line_axis(stations)
    x, z, V, slon, slat = cut(f, dep, p0, p1, zmax=zmax)
    el = np.array([st["elev"] for st in stations.values()])
    north_south = abs(u[1]) >= abs(u[0])
    topo = relief_along(cfg, slon, slat)
    rows = [
        dict(
            x=x,
            z=z,
            v=V,
            sta=s,
            elev=el,
            label=label or cfg["project"]["name"],
            zmax=zmax,
            ends=compass(p0, p1),
            topo=(x, topo) if topo is not None else None,
        )
    ]
    return draw(
        cfg, rows, out, xlabel="Distance along the line from the %s end (km)" % ("south" if north_south else "west")
    )


def line_checkerboard(cfg, project, stations, out, label=None):
    """Checkerboard of a line: input and recovered anomaly (per cent of the 1-D background of each depth) as two
    sections along the line axis."""
    from ..inversion.run import read_model

    M = open(os.path.join(project, "MOD")).read().split("\n")
    zmod = [round(float(v), 3) for v in M[0].split()]
    nrow = len([l for l in M[1:] if l.strip()]) // len(zmod)
    bg = np.array([float(M[1 + k * nrow].split()[0]) for k in range(len(zmod))])
    zmax = model_zmax(cfg, project)
    p0, p1, u, s, _ = line_axis(stations)
    rows = []
    for fname, title in (("Vs_model.real", "Input"), ("DSurfTomo.inSyn.dat", "Recovered")):
        lon, lat, dep, v = read_model(os.path.join(project, fname))
        b = np.interp(dep, zmod, bg)
        f = RegularGridInterpolator((lon, lat, dep), 100 * (v - b) / b, bounds_error=False, fill_value=np.nan)
        x, z, V, _, _ = cut(f, dep, p0, p1, zmax=zmax)
        rows.append(
            dict(
                x=x, z=z, v=V, sta=s, elev=None, label="%s: %s" % (label or os.path.basename(project), title), zmax=zmax
            )
        )
    r = 5 * np.ceil(100 * cfg["checks"]["amplitude"] / 5)
    north_south = abs(u[1]) >= abs(u[0])
    return draw(
        cfg,
        rows,
        out,
        xlabel="Distance along the line from the %s end (km)" % ("south" if north_south else "west"),
        scale=dict(cpt=G.CHECK_CPT, lo=-r, hi=r, label="@~D@~Vs (%)"),
    )
