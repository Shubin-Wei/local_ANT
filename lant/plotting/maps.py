"""Map figures drawn with GMT: station map, depth slices of the model, checkerboard maps, 3-D view of the slices.

Style (as in ambient-noise papers): hill-shaded relief (plots.topography; coloured under the station map, light grey
under the model maps, plain background without relief), model values resampled bilinearly for display (node values
unchanged), clipped to the convex hull of the stations and drawn 35 % transparent over the relief, GMT "seis" colour
table red (slow) to blue (fast) for Vs, red - white - blue for the checkerboards, one scale bar, one colour bar per
figure.  No geological information is drawn.
"""

import os

import numpy as np
from scipy.spatial import ConvexHull

from . import gmt as G


def _region(lon, lat, pad=0.1, min_pad=0.005, min_aspect=0.4):
    """Map region around the stations: pad on each side; the shorter side widened to at least min_aspect x the
    longer one (in km), so that a line array is not drawn as a thin strip (display choice)."""
    W, E, S, N = min(lon), max(lon), min(lat), max(lat)
    px = max((E - W) * pad, min_pad)
    py = max((N - S) * pad, min_pad)
    W, E, S, N = W - px, E + px, S - py, N + py
    kx = np.cos(np.radians(0.5 * (S + N)))
    wx, wy = (E - W) * kx, N - S
    if wx < min_aspect * wy:
        d = 0.5 * (min_aspect * wy - wx) / kx
        W, E = W - d, E + d
    elif wy < min_aspect * wx:
        d = 0.5 * (min_aspect * wx - wy)
        S, N = S - d, N + d
    return W, E, S, N


def _model_region(stations, lon, lat, traces=None):
    """Region of the model maps: the stations (and the ends of the section lines) with the padding of the station
    map, inside the model grid (the grid's outer ring of nodes lies far outside the array and is not drawn)."""
    plon = [s["lon"] for s in stations.values()] + [p[0] for t in (traces or []) for p in t[1:]]
    plat = [s["lat"] for s in stations.values()] + [p[1] for t in (traces or []) for p in t[1:]]
    W, E, S, N = _region(plon, plat)
    return (
        max(W, float(np.min(lon))),
        min(E, float(np.max(lon))),
        max(S, float(np.min(lat))),
        min(N, float(np.max(lat))),
    )


def _ticks(W, E, S, N):
    return G.nice_step(E - W, 3), G.nice_step(N - S, 4)


def _geo_format(*steps):
    """--FORMAT_GEO_MAP with as many decimals as the annotation steps need (0.025 degrees: three), at least two."""
    nd = max([2] + [len(("%.6f" % s).rstrip("0").split(".")[1]) for s in steps])
    return "--FORMAT_GEO_MAP=ddd.%sF" % ("x" * nd)


def _write(path, text):
    with open(path, "w") as f:
        f.write(text)


def station_map(cfg, stations, out, excluded=None, title=None):
    """Stations over coloured, hill-shaded relief; excluded stations (data.exclude_stations) as red crosses; legend
    outside the map, elevation colour bar, scale bar."""
    excluded = excluded or {}
    allst = list(stations.values()) + list(excluded.values())
    W, E, S, N = _region([s["lon"] for s in allst], [s["lat"] for s in allst])
    tmp = G.workdir()
    np.savetxt(tmp + "/sta.xy", np.array([(s["lon"], s["lat"]) for s in stations.values()]), fmt="%.6f")
    if excluded:
        np.savetxt(tmp + "/exc.xy", np.array([(s["lon"], s["lat"]) for s in excluded.values()]), fmt="%.6f")
    fs = cfg["plots"]["font_scale"]
    width = 14.0
    H = G.mercator_height(width, W, E, S, N)
    P = "-R%.5f/%.5f/%.5f/%.5f -JM%gc -Xa2.2c -Ya1.6c" % (W, E, S, N, width)
    pre, relief = G.relief_cmds(cfg, tmp, W, E, S, N, palette="elevation")
    dx, dy = _ticks(W, E, S, N)
    km, latc = G.scale_km(W, E, S, N)
    sym = 0.32 * fs * (1.0 if len(stations) < 120 else 0.8)
    cmd = []
    if relief:
        cmd.append("gmt grdimage %s/topo.nc -I%s/shade.nc -C%s/relief.cpt %s" % (tmp, tmp, tmp, P))
    cmd.append("gmt plot %s/sta.xy %s -St%.2fc -Gwhite -W0.6p,black" % (tmp, P, sym))
    if excluded:
        cmd.append("gmt plot %s/exc.xy %s -Sx%.2fc -W1.6p,red" % (tmp, P, 0.8 * sym))
    _write(tmp + "/title.txt", (title or cfg["project"]["name"]) + "\n")
    cmd += [
        "gmt basemap %s -BWSne -Bxa%gf%g -Bya%gf%g --FONT_ANNOT_PRIMARY=%gp,Helvetica %s"
        % (P, dx, dx / 4, dy, dy / 4, 11 * fs, _geo_format(dx, dy)),
        "gmt text %s/title.txt -R0/1/0/1 -JX%gc/%gc -Xa2.2c -Ya1.6c -F+cTC+f%gp,Helvetica+jBC -D0/0.35c -N"
        % (tmp, width, H, 14 * fs),
        "gmt basemap %s -LjBL+c%.4f+w%gk+f+l'%g km'+o0.6c/0.6c --MAP_SCALE_HEIGHT=0.15c --FONT_ANNOT_PRIMARY=%gp,Helvetica "
        "--FONT_LABEL=%gp,Helvetica" % (P, latc, km, km, 10 * fs, 10 * fs),
    ]
    leg = "S 0.3c t %.2fc white 0.6p,black 0.75c Stations used (%d)\n" % (sym, len(stations))
    if excluded:
        leg += "S 0.3c x %.2fc - 1.6p,red 0.75c Stations excluded (%d)\n" % (0.8 * sym, len(excluded))
    _write(tmp + "/legend.txt", leg)
    cmd.append(
        "gmt legend %s/legend.txt %s -DJTR+jTL+w6.2c+o0.4c/0c -F+p0.5p+gwhite --FONT_ANNOT_PRIMARY=%gp,Helvetica"
        % (tmp, P, 11 * fs)
    )
    if relief:
        cmd.append(
            "gmt colorbar -C%s/relief.cpt %s -DJBR+jBL+w%.2fc/0.35c+o0.6c/0c -Bxaf+l'Elevation (m)' "
            "--FONT_ANNOT_PRIMARY=%gp,Helvetica --FONT_LABEL=%gp,Helvetica"
            % (tmp, P, min(5.0, 0.45 * H), 10 * fs, 11 * fs)
        )
    return G.run(tmp, pre, cmd, out, cfg["plots"]["formats"], cfg["plots"]["dpi"])


def _hull(stations, tmp):
    xy = np.array([(s["lon"], s["lat"]) for s in stations.values()])
    h = ConvexHull(xy)
    ring = xy[np.r_[h.vertices, h.vertices[:1]]]
    np.savetxt(tmp + "/hull.xy", ring, fmt="%.6f %.6f")
    return ring


def panels(
    cfg,
    tmp,
    layers,
    labels,
    stations,
    region,
    cpt,
    rng,
    annot,
    label,
    out,
    ncol=4,
    transp=35,
    panel_width=4.6,
    show_stations=True,
    boxed=False,
    extend=False,
    traces=None,
):
    """Map panels: layers = list of (lon, lat, value) arrays on a regular grid, clipped to the station hull.
    boxed: the label of each panel in a white box at its top-left corner (otherwise a title above it).
    traces: section lines [(name, (lon, lat), (lon, lat))], dashed, end letters on the first panel."""
    W, E, S, N = region
    fs = cfg["plots"]["font_scale"]
    PH = G.mercator_height(panel_width, W, E, S, N)
    n = len(layers)
    ncol = min(ncol, n)
    nrow = (n + ncol - 1) // ncol
    gapx, gapy, x0, y0 = 0.25 * fs, (0.3 if boxed else 0.75) * fs, 1.6, 2.4
    pre = ["gmt makecpt -C%s -T%s -Z %s > %s/v.cpt" % (cpt, rng, "-D" if extend else "-D", tmp)]
    r_pre, relief = G.relief_cmds(cfg, tmp, W, E, S, N)
    pre += r_pre
    np.savetxt(tmp + "/sta.xy", np.array([(s["lon"], s["lat"]) for s in stations.values()]), fmt="%.6f")
    _hull(stations, tmp)
    if traces:
        with open(tmp + "/traces.gmt", "w") as fo, open(tmp + "/ends.txt", "w") as fe:
            for name, p0, p1 in traces:
                fo.write(">\n%.6f %.6f\n%.6f %.6f\n" % (p0[0], p0[1], p1[0], p1[1]))
                fe.write("%.6f %.6f %s\n%.6f %.6f %s'\n" % (p0[0], p0[1], name, p1[0], p1[1], name))
    R = "-R%.5f/%.5f/%.5f/%.5f" % region
    dx, dy = _ticks(W, E, S, N)
    cmd = []
    for k, (lo, la, v) in enumerate(layers):
        ulo, ula = np.unique(lo), np.unique(la)
        ddx, ddy = np.diff(ulo).min(), np.diff(ula).min()
        np.savetxt("%s/val_%d.xyz" % (tmp, k), np.c_[lo, la, v], fmt="%.6f %.6f %.5f")
        pre += [
            "gmt xyz2grd %s/val_%d.xyz -R%.6f/%.6f/%.6f/%.6f -I%.7f/%.7f -G%s/a%d.nc"
            % (tmp, k, ulo[0], ulo[-1], ula[0], ula[-1], ddx, ddy, tmp, k),
            "gmt grdsample %s/a%d.nc -I%.7f -nl -G%s/s%d.nc" % (tmp, k, min(ddx, ddy) / 5, tmp, k),
        ]
        col, row = k % ncol, k // ncol
        P = "%s -JM%gc -Xa%.3fc -Ya%.3fc" % (
            R,
            panel_width,
            x0 + col * (panel_width + gapx),
            y0 + (nrow - 1 - row) * (PH + gapy),
        )
        if relief:
            cmd.append("gmt grdimage %s/topo.nc -I%s/shade.nc -C%s/relief.cpt %s" % (tmp, tmp, tmp, P))
        cmd += [
            "gmt clip %s/hull.xy %s" % (tmp, P),
            "gmt grdimage %s/s%d.nc -C%s/v.cpt -Q %s %s" % (tmp, k, tmp, "-t%d" % transp if relief else "", P),
            "gmt clip -C %s" % P,
        ]
        if traces:
            cmd.append("gmt plot %s/traces.gmt -W%.2fp,black,4_2:0 %s" % (tmp, 0.7 * fs, P))
            if k == 0:
                cmd += [
                    "gmt text %s/ends.txt -F+f%gp,Helvetica-Bold,white=1.2p,white+jCM %s" % (tmp, 9 * fs, P),
                    "gmt text %s/ends.txt -F+f%gp,Helvetica-Bold+jCM %s" % (tmp, 9 * fs, P),
                ]
        if show_stations:
            cmd.append("gmt plot %s/sta.xy -St%.2fc -W0.3p,black %s" % (tmp, 0.12 * fs, P))
        cmd.append(
            "gmt basemap %s -B%s%sne -Bxa%gf%g -Bya%gf%g %s"
            % (
                P,
                "W" if col == 0 else "w",
                "S" if row == nrow - 1 else "s",
                dx,
                dx / 4,
                dy,
                dy / 4,
                _geo_format(dx, dy),
            )
        )
        _write("%s/label_%d.txt" % (tmp, k), labels[k] + "\n")
        if boxed:
            cmd.append(
                "gmt text %s/label_%d.txt %s -F+cTL+f%gp,Helvetica-Bold+jTL -D0.08c/-0.08c -Gwhite -W0.4p -C0.05c/0.04c"
                % (tmp, k, P, 9 * fs)
            )
        else:
            cmd.append("gmt text %s/label_%d.txt %s -F+cTC+f%gp,Helvetica+jBC -D0/0.12c -N" % (tmp, k, P, 10 * fs))
        if k == n - 1:
            km, latc = G.scale_km(W, E, S, N)
            cmd.append(
                "gmt basemap %s -LjBR+c%.4f+w%gk+o0.3c/0.35c+l'%g km' --MAP_SCALE_HEIGHT=0.1c --FONT_LABEL=%gp,Helvetica"
                % (P, latc, km, km, 8 * fs)
            )
    total = ncol * (panel_width + gapx) - gapx
    cmd.append(
        "gmt colorbar -C%s/v.cpt -Dx%.3fc/%.3fc+w%.2fc/0.3c+h+jTC%s -Bx%s+l'%s' --FONT_ANNOT_PRIMARY=%gp,Helvetica "
        "--FONT_LABEL=%gp,Helvetica"
        % (tmp, x0 + total / 2, y0 - 0.9, min(9.0, total * 0.7), "+e" if extend else "", annot, label, 9 * fs, 10 * fs)
    )
    return G.run(tmp, pre, cmd, out, cfg["plots"]["formats"], cfg["plots"]["dpi"], font=(7 * fs, 8 * fs))


def depth_slices(cfg, model, stations, out, depths=None, traces=None):
    """Vs at the slice depths (plots.slice_depths_km), over the relief; the traces of the sections, no stations."""
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
    return panels(
        cfg,
        tmp,
        layers,
        ["%g km" % z for z in use],
        stations,
        _model_region(stations, lon, lat, traces),
        cfg["plots"]["vs_cpt"],
        "%g/%g" % (lo_, hi_),
        colorbar_ticks(lo_, hi_),
        "Vs (km/s)",
        out,
        show_stations=False,
        traces=traces,
    )


def colorbar_ticks(lo, hi):
    """GMT colour-bar annotation: about five labelled steps of the 1-2-2.5-5 sequence, minor ticks at a fifth."""
    from ..signal import nice_number

    a = nice_number((hi - lo) / 5, "nearest")
    return "a%gf%g" % (a, a / 5)


def checker_layers(project):
    """(true, recovered) arrays lon lat depth v, and the 1-D background of each depth (MOD)."""
    syn = np.loadtxt(os.path.join(project, "DSurfTomo.inSyn.dat"))
    tru = np.loadtxt(os.path.join(project, "Vs_model.real"))
    M = open(os.path.join(project, "MOD")).read().split("\n")
    zmod = [round(float(v), 3) for v in M[0].split()]
    nrow = len([l for l in M[1:] if l.strip()]) // len(zmod)
    bg = {z: float(M[1 + k * nrow].split()[0]) for k, z in enumerate(zmod)}
    return tru, syn, bg


def checkerboard(cfg, project, stations, out, depths=None, vertical=False):
    """Input pattern and the recovered model at up to seven depths, in per cent of the 1-D background, red (slow)
    - white - blue (fast); labels boxed in the panels.  vertical: the input changes with depth, so every recovered
    depth is shown next to its input (input and recovery in alternate panels)."""
    tru, syn, bg = checker_layers(project)
    dep = np.unique(np.round(syn[:, 2], 3))
    avail = [z for z in cfg["plots"]["slice_depths_km"] if np.min(np.abs(dep - z)) < 1e-3]
    depths = depths or (avail[:4] if vertical else avail[:7])
    layers, labels = [], []

    def lay(a, z):
        m = np.round(a[:, 2], 3) == z
        return (a[m, 0], a[m, 1], 100 * (a[m, 3] - bg[z]) / bg[z])

    if vertical:
        for z in depths:
            layers += [lay(tru, z), lay(syn, z)]
            labels += ["Input %g km" % z, "%g km" % z]
    else:
        layers.append(lay(tru, depths[0]))
        labels.append("Input")
        for z in depths:
            layers.append(lay(syn, z))
            labels.append("%g km" % z)
    amp = 100 * cfg["checks"]["amplitude"]
    r = 5 * np.ceil(amp / 5)
    tmp = G.workdir()
    return panels(
        cfg,
        tmp,
        layers,
        labels,
        stations,
        _model_region(stations, syn[:, 0], syn[:, 1]),
        G.CHECK_CPT,
        "-%g/%g" % (r, r),
        "a%gf%g" % (r / 2, r / 10),
        "@~D@~Vs (%)",
        out,
        ncol=4,
        boxed=True,
        extend=True,
    )


def slices_3d(cfg, model, stations, out, depths=None):
    """Areal arrays: Vs slices at five depths stacked in a perspective view on a depth axis drawn to scale, the
    hill-shaded relief floating above them; each slice clipped to the station hull."""
    from ..inversion.run import read_model

    lon, lat, dep, v = read_model(model)
    avail = [z for z in cfg["plots"]["slice_depths_km"] if np.min(np.abs(dep - z)) < 1e-3]
    if not depths:
        idx = np.unique(np.round(np.linspace(0, len(avail) - 1, min(5, len(avail)))).astype(int))
        depths = [avail[i] for i in idx]
    fs = cfg["plots"]["font_scale"]
    ring = np.array([(s["lon"], s["lat"]) for s in stations.values()])
    ring = ring[np.r_[ConvexHull(ring).vertices, ConvexHull(ring).vertices[:1]]]
    W, E, S, N = ring[:, 0].min(), ring[:, 0].max(), ring[:, 1].min(), ring[:, 1].max()
    pw = 0.03 * max(E - W, N - S)
    W, E, S, N = W - pw, E + pw, S - pw, N + pw
    tmp = G.workdir()
    PW, VIEW = 10.0, (155, 28)
    zmax = max(depths)
    TOP, RELIEF_Z = 0.75 * zmax, 0.12 * zmax  # base of the relief above depth 0 and its height (depth units)
    ZD = [-z for z in depths]
    ZLO, ZHI = -zmax - 0.02 * zmax, TOP + RELIEF_Z + 0.02 * zmax
    HZ = 17.0  # height of the depth axis (cm)

    def P(z=None, zhi=ZHI):
        return "-R%.5f/%.5f/%.5f/%.5f/%g/%g -JM%gc -JZ%gc -p%g/%g%s -Xa2.6c -Ya3.4c" % (
            W,
            E,
            S,
            N,
            ZLO,
            zhi,
            PW,
            HZ * (zhi - ZLO) / (ZHI - ZLO),
            VIEW[0],
            VIEW[1],
            "" if z is None else "/%g" % z,
        )

    np.savetxt(tmp + "/hull.xy", ring, fmt="%.6f")
    LO, LA = np.meshgrid(lon, lat, indexing="ij")
    lo_, hi_ = cfg["plots"]["vs_range"]
    pre = ["gmt makecpt -C%s -T%g/%g -Z -D > %s/vs.cpt" % (cfg["plots"]["vs_cpt"], lo_, hi_, tmp)]
    ddx, ddy = lon[1] - lon[0], lat[1] - lat[0]
    for k, z in enumerate(depths):
        i = int(np.argmin(np.abs(dep - z)))
        np.savetxt("%s/v%d.xyz" % (tmp, k), np.c_[LO.ravel(), LA.ravel(), v[:, :, i].ravel()], fmt="%.6f %.6f %.4f")
        np.savetxt("%s/hull%d.xyz" % (tmp, k), np.c_[ring, np.full(len(ring), ZD[k])], fmt="%.6f")
        pre += [
            "gmt xyz2grd %s/v%d.xyz -R%.6f/%.6f/%.6f/%.6f -I%.7f/%.7f -G%s/a%d.nc"
            % (tmp, k, lon[0], lon[-1], lat[0], lat[-1], ddx, ddy, tmp, k),
            "gmt grdsample %s/a%d.nc -I%.7f -nl -G%s/s%d.nc" % (tmp, k, min(ddx, ddy) / 5, tmp, k),
        ]
    r_pre, relief = G.relief_cmds(cfg, tmp, W, E, S, N, palette="elevation")
    pre += r_pre
    if relief:  # elevations mapped onto [TOP, TOP + RELIEF_Z] of the depth axis (relief exaggerated)
        pre += [
            "read ZMIN ZMAX <<< \"$(gmt grdinfo -C %s/topo.nc | awk '{print $6, $7}')\"" % tmp,
            "gmt grdsample %s/topo.nc -R%.5f/%.5f/%.5f/%.5f -I%.7f -G%s/topo_s.nc"
            % (tmp, W, E, S, N, max(E - W, N - S) / 400, tmp),
            "gmt grdgradient %s/topo_s.nc -A315 -Nt0.45 -G%s/shade_s.nc" % (tmp, tmp),
            "gmt grdmath %s/topo_s.nc $ZMIN SUB $ZMAX $ZMIN SUB 1 MAX DIV %g MUL %g ADD = %s/topoz.nc"
            % (tmp, RELIEF_Z, TOP, tmp),
            "gmt makecpt -C%s -T%g/%g -Z > %s/reliefz.cpt" % (G.ELEV_CPT, TOP, TOP + RELIEF_Z, tmp),
        ]
    zlab, last = [], None  # labels closer than 6 % of the axis to the last one: tick only
    for z in [0] + list(depths):
        lab = last is None or z - last >= 0.06 * zmax or z == depths[-1]
        if lab and zlab and z == depths[-1] and z - last < 0.06 * zmax:
            zlab[-1] = (zlab[-1][0], False)  # the deepest slice keeps its label, the one above gives way
        zlab.append((z, lab))
        last = z if lab else last
    _write(
        tmp + "/zann.txt", "".join("%g %s %s\n" % (-z, "a" if lab else "f", "%g" % z if lab else "") for z, lab in zlab)
    )
    dx, dy = _ticks(W, E, S, N)
    for nm, a, b, step, fmt in (("x", W, E, dx, "%.4f"), ("y", S, N, dy, "%.4f")):  # inner annotations only: the
        t = np.arange(np.ceil(a / step) * step, b, step / 4)  # corners carry the depth axis
        inner = lambda x: a + 0.15 * (b - a) < x < b - 0.15 * (b - a)
        _write(
            "%s/%sann.txt" % (tmp, nm),
            "".join(
                (fmt + " %s\n") % (x, "af" if abs(x / step - round(x / step)) < 1e-6 and inner(x) else "f") for x in t
            ),
        )
    cmd = ["gmt set MAP_ANNOT_ORTHO wez FONT_ANNOT_PRIMARY %gp,Helvetica FONT_LABEL %gp,Helvetica" % (10 * fs, 11 * fs)]
    for k in reversed(range(len(depths))):
        z = ZD[k]
        cmd += [
            "gmt clip %s/hull.xy %s" % (tmp, P(z)),
            "gmt grdimage %s/s%d.nc -C%s/vs.cpt -Q %s" % (tmp, k, tmp, P(z)),
            "gmt clip -C %s" % P(z),
            "gmt plot3d %s/hull%d.xyz -W0.4p,gray25 %s" % (tmp, k, P()),
        ]
        if k == len(depths) - 1:
            cmd.append(
                "gmt basemap %s -BSE -Bxc%s/xann.txt -Byc%s/yann.txt --MAP_FRAME_TYPE=plain %s"
                % (P(z), tmp, tmp, _geo_format(dx, dy))
            )
    cmd.append("gmt basemap %s -Bzc%s/zann.txt+l'Depth (km)' -BZ --MAP_LABEL_OFFSET=%gc" % (P(zhi=0), tmp, 0.9 * fs))
    if relief:
        cmd.append("gmt grdview %s/topoz.nc -I%s/shade_s.nc -C%s/reliefz.cpt -Qi300 %s" % (tmp, tmp, tmp, P()))
    R2 = "-R%.5f/%.5f/%.5f/%.5f -JM%gc -Xa2.6c -Ya3.4c" % (W, E, S, N, PW)
    cmd.append(
        "gmt colorbar -C%s/vs.cpt -DjBC+w6c/0.3c+h+o%gc/-2.6c -B%s+l'Vs (km/s)' %s"
        % (tmp, -2.2 if relief else 0, colorbar_ticks(lo_, hi_), R2)
    )
    if relief:
        pre.append(
            "gmt makecpt -C%s -T$ZMIN/$(awk -v a=$ZMIN -v b=$ZMAX 'BEGIN{print (b>a)?b:a+1}') -Z > %s/elev.cpt"
            % (G.ELEV_CPT, tmp)
        )
        pre.append(
            "ESTEP=$(awk -v a=$ZMIN -v b=$ZMAX 'BEGIN{r=(b-a)/2; if(r<=0)r=1; p=10^int(log(r)/log(10)); if(p>r)p/=10; "
            "m=r/p; s=(m<1.5)?1:(m<3.5)?2:(m<7.5)?5:10; print s*p}')"
        )
        cmd.append(
            "gmt colorbar -C%s/elev.cpt -DjBC+w4c/0.3c+h+o4.6c/-2.6c -Bxa${ESTEP}f$(awk -v s=$ESTEP 'BEGIN{print s/5}')"
            "+l'Elevation (m)' %s" % (tmp, R2)
        )
    return G.run(tmp, pre, cmd, out, cfg["plots"]["formats"], cfg["plots"]["dpi"], font=(9 * fs, 10 * fs))
