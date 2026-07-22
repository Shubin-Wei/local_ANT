#!/usr/bin/env python3
"""
Stitch per-pair 1D Vs models (from lant_invert.py) into a pseudo-2D section.

Each model is assigned to its pair midpoint. Midpoints are projected onto the
principal direction of the array (PCA) to give an along-line distance, then the
1D Vs(depth) columns are binned/averaged along distance onto a regular grid and
plotted as a depth section.

All parameters (paths, misfit cutoff, bins, smoothing, contours, color scale)
are read from the YAML config; section.png/.npz are written to output_dir.

Usage:
    python lant_stitch.py --config my_line.yaml   # same config as lant_invert.py
"""

import os
import glob
import csv
import argparse

import yaml
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def fill_nan_rows(section, x):
    """Fill empty (NaN) distance bins by linear interpolation along each depth row."""
    out = section.copy()
    for i in range(out.shape[0]):
        row = out[i]
        ok = ~np.isnan(row)
        if ok.sum() >= 2:
            out[i] = np.interp(x, x[ok], row[ok])
        elif ok.sum() == 1:
            out[i] = row[ok][0]
    return out


def apply_depth_smoothing(section, depth, default_sigma=0.0, bands=None):
    """Smooth shallow/mid/deep rows with different Gaussian sigmas."""
    if not bands:
        return gaussian_filter(section, sigma=default_sigma) if default_sigma > 0 else section

    out = np.full_like(section, np.nan)
    start_depth = -np.inf
    for band in bands:
        max_depth = band.get("max_depth")
        sigma = float(band.get("sigma", default_sigma))
        if max_depth is None:
            row_mask = depth >= start_depth
        else:
            max_depth = float(max_depth)
            row_mask = (depth >= start_depth) & (depth < max_depth)
        smoothed = gaussian_filter(section, sigma=sigma) if sigma > 0 else section
        out[row_mask] = smoothed[row_mask]
        if max_depth is not None:
            start_depth = max_depth

    # If the configured bands do not cover all depths, fall back to the default.
    uncovered = ~np.isfinite(out).all(axis=1)
    if np.any(uncovered):
        fallback = gaussian_filter(section, sigma=default_sigma) if default_sigma > 0 else section
        out[uncovered] = fallback[uncovered]
    return out


def describe_smoothing(default_sigma, bands):
    if not bands:
        return f"uniform sigma={default_sigma:g}"
    start = 0.0
    parts = []
    for band in bands:
        max_depth = band.get("max_depth")
        sigma = float(band.get("sigma", default_sigma))
        if max_depth is None:
            parts.append(f">={start:g} km: sigma={sigma:g}")
        else:
            parts.append(f"{start:g}-{float(max_depth):g} km: sigma={sigma:g}")
            start = float(max_depth)
    return "; ".join(parts)


def load_models(models_dir, misfit_max):
    """Load every models/<pair>.txt -> (midlon, midlat, misfit, depth, vs)."""
    # accept either the project dir (contains models/) or the models/ dir itself
    paths = sorted(glob.glob(os.path.join(models_dir, "models", "*.txt")))
    if not paths:
        paths = sorted(glob.glob(os.path.join(models_dir, "*.txt")))
    recs = []
    depth_ref = None
    for path in paths:
        with open(path) as f:
            header = f.readline().split()
        meta = {header[i]: header[i + 1] for i in range(1, len(header) - 1, 2)}
        misfit = float(meta.get("misfit", "nan"))
        if misfit_max is not None and misfit > misfit_max:
            continue
        d = np.loadtxt(path)
        depth, vs = d[:, 0], d[:, 1]
        if depth_ref is None:
            depth_ref = depth
        recs.append((float(meta["midlon"]), float(meta["midlat"]), misfit, vs))
    return recs, depth_ref


def project_to_line(lons, lats, transform=None):
    """Project coordinates onto the model-midpoint principal axis in km.

    With ``transform=None`` the line is fitted to the supplied coordinates and
    both the distances and fitted transform are returned.  Reusing that
    transform for station coordinates keeps the topography and Vs section on
    exactly the same along-line axis.
    """
    lons = np.asarray(lons, dtype=float)
    lats = np.asarray(lats, dtype=float)
    if transform is not None:
        lon0, lat0, direction, shift = transform
        x = (lons - lon0) * 111.195 * np.cos(np.radians(lat0))
        y = (lats - lat0) * 111.195
        return np.column_stack([x, y]) @ direction - shift

    lon0 = float(np.mean(lons))
    lat0 = float(np.mean(lats))
    x = (lons - lon0) * 111.195 * np.cos(np.radians(lat0))
    y = (lats - lat0) * 111.195
    pts = np.column_stack([x, y])
    # principal direction = first right-singular vector
    _, _, vt = np.linalg.svd(pts - pts.mean(0), full_matrices=False)
    direction = vt[0]
    # Fix PCA's otherwise arbitrary sign so profile distance always increases
    # from south to north (left=S, right=N in the output section).
    if direction[1] < 0:
        direction = -direction
    raw_dist = pts @ direction
    shift = float(raw_dist.min())
    return raw_dist - shift, (lon0, lat0, direction, shift)


def load_topography(station_csv, transform, elevation_min=None):
    """Read station lon/lat/elevation and project them onto the section axis."""
    aliases = {
        "longitude": ("longitude", "lon", "lng"),
        "latitude": ("latitude", "lat"),
        "elevation": ("elevation", "elev", "altitude", "height"),
    }
    with open(station_csv, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        fields = {name.strip().lower(): name for name in (reader.fieldnames or [])}

        columns = {}
        for canonical, candidates in aliases.items():
            match = next((fields[name] for name in candidates if name in fields), None)
            if match is None:
                raise ValueError(
                    f"{station_csv}: missing {canonical} column; found {reader.fieldnames}")
            columns[canonical] = match

        rows = []
        for row in reader:
            try:
                lon = float(row[columns["longitude"]])
                lat = float(row[columns["latitude"]])
                elevation = float(row[columns["elevation"]])
            except (TypeError, ValueError):
                continue
            if (np.all(np.isfinite([lon, lat, elevation])) and
                    (elevation_min is None or elevation >= elevation_min)):
                rows.append((lon, lat, elevation))

    if len(rows) < 2:
        raise ValueError(f"{station_csv}: need at least two valid station elevations.")
    rows = np.asarray(rows)
    distance = project_to_line(rows[:, 0], rows[:, 1], transform=transform)
    order = np.argsort(distance)
    return distance[order], rows[order, 2]


def profile_end_labels(line_transform):
    direction = line_transform[2]
    north_at_left = direction[1] < 0
    return ("N", "S") if north_at_left else ("S", "N")


def relative_velocity(velocity):
    """Return lateral Vs perturbation (%) and the depth-wise 1-D reference."""
    reference = np.nanmean(velocity, axis=1, keepdims=True)
    return 100.0 * (velocity - reference) / reference, reference[:, 0]


def section_geometry(x, depth, topo_dist=None, topo_elevation=None):
    """Build pcolormesh coordinates for depth or topography-following plots."""
    if topo_dist is not None and topo_elevation is not None:
        surface = np.interp(x, topo_dist, topo_elevation) / 1000.0
        plot_x = np.broadcast_to(x, (len(depth), len(x)))
        plot_y = surface[np.newaxis, :] - depth[:, np.newaxis]
    else:
        surface = None
        plot_x, plot_y = x, depth
    return plot_x, plot_y, surface


def taper_vertices(x, depth, taper_slope, surface=None):
    dmax = float(np.max(x))
    zmax = float(np.max(depth))
    inset = min(float(taper_slope or 0) * zmax, 0.45 * dmax)
    if surface is not None:
        top_x = np.linspace(0.0, dmax, max(200, len(x) * 2))
        top_y = np.interp(top_x, x, surface)
        bottom_x = np.linspace(dmax - inset, inset, max(100, len(x)))
        bottom_y = np.interp(bottom_x, x, surface) - zmax
        return list(zip(top_x, top_y)) + list(zip(bottom_x, bottom_y))
    return [(0.0, 0.0), (dmax, 0.0), (dmax - inset, zmax), (inset, zmax)]


def apply_taper(ax, artists, vertices):
    from matplotlib.patches import Polygon
    clip = Polygon(vertices, closed=True, transform=ax.transData)
    for artist in artists:
        if artist is None:
            continue
        try:
            artist.set_clip_path(clip)
        except Exception:
            for item in getattr(artist, "collections", []):
                item.set_clip_path(clip)
    ax.plot([v[0] for v in vertices] + [vertices[0][0]],
            [v[1] for v in vertices] + [vertices[0][1]], color="k", lw=0.8)


def finish_section_axes(ax, x, depth, surface, line_transform, title):
    dmax = float(np.max(x))
    zmax = float(np.max(depth))
    ax.set_xlim(0, dmax)
    if surface is not None:
        vertical_margin = max(0.015, 0.04 * zmax)
        ax.set_ylim(float((surface - zmax).min()) - vertical_margin,
                    float(surface.max()) + vertical_margin)
        ax.set_ylabel("Elevation (km)")
        if line_transform is not None:
            left_label, right_label = profile_end_labels(line_transform)
            ax.text(0.018, 0.978, left_label, transform=ax.transAxes,
                    ha="left", va="top", fontsize=9, fontweight="bold", zorder=10)
            ax.text(0.982, 0.978, right_label, transform=ax.transAxes,
                    ha="right", va="top", fontsize=9, fontweight="bold", zorder=10)
    else:
        ax.invert_yaxis()
        ax.set_ylabel("Depth (km)")
    ax.set_xlabel("Distance along profile (km)")
    ax.set_title(title, fontsize=13, pad=9)
    ax.set_aspect("auto")


def plot_relative_section(x, depth, velocity, config, line_transform=None,
                          topo_dist=None, topo_elevation=None, output_path=None,
                          limit=None, save_pdf=None):
    """Plot and save the lateral relative Vs section."""
    st = config["stitch"]
    output_dir = config["paths"]["output_dir"]
    relative, reference = relative_velocity(velocity)
    plot_x, plot_y, surface = section_geometry(x, depth, topo_dist, topo_elevation)
    dmax = float(np.max(x))
    zmax = float(np.max(depth))
    taper_slope = float(st.get("taper_slope") or 0)
    relative_limit = float(limit if limit is not None else st.get("relative_limit", 30.0))

    fig, ax = plt.subplots(figsize=(12, 6.75))
    pcm = ax.pcolormesh(plot_x, plot_y, relative, cmap=st.get("relative_cmap", "RdBu"),
                        vmin=-relative_limit, vmax=relative_limit, shading="gouraud")

    if taper_slope > 0:
        vertices = taper_vertices(x, depth, taper_slope, surface)
        apply_taper(ax, [pcm], vertices)

    if surface is not None:
        ax.plot(x, surface, color="k", lw=1.2, zorder=5)
        if bool(st.get("show_stations", True)) and topo_dist is not None and topo_elevation is not None:
            visible = (topo_dist >= 0) & (topo_dist <= dmax)
            ax.scatter(topo_dist[visible], topo_elevation[visible] / 1000.0,
                       marker="v", s=10, facecolor="k", edgecolor="none",
                       zorder=6, clip_on=True)

    absolute_title = st.get("title") or "S-wave Velocity Section"
    title = absolute_title.replace("S-wave Velocity Section",
                                   "Relative S-wave Velocity Section")
    finish_section_axes(ax, x, depth, surface, line_transform, title)

    cb = fig.colorbar(pcm, ax=ax, pad=0.025, extend="neither")
    cb.set_label("Relative Vs (%)")
    cb.ax.tick_params(labelsize=8)
    fig.tight_layout()

    out = output_path or os.path.join(output_dir, "relative_velocity.png")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    dpi = int(st.get("dpi", 300))
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"Saved {out}")

    if save_pdf is None:
        save_pdf = bool(st.get("save_pdf", False))
    if save_pdf and os.path.splitext(out)[1].lower() != ".pdf":
        pdf_path = os.path.splitext(out)[0] + ".pdf"
        fig.savefig(pdf_path, bbox_inches="tight")
        print(f"Saved {pdf_path}")

    npz = os.path.splitext(out)[0] + ".npz"
    saved = dict(distance=x, depth=depth, relative_vs_percent=relative,
                 reference_vs=reference)
    if topo_dist is not None:
        saved.update(topography_distance=topo_dist,
                     topography_elevation=topo_elevation)
    if line_transform is not None:
        saved.update(axis_labels=np.asarray(profile_end_labels(line_transform)))
    np.savez(npz, **saved)
    print(f"Saved {npz}")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(
        description="Stitch 1D Vs models into a pseudo-2D section. Defaults come "
                    "from --config; any flag below overrides the config value.")
    p.add_argument("--config", required=True,
                   help="path to YAML config (template: configs/lant.example.yaml).")
    # optional overrides of the config's 'stitch' section (handy for quick tuning)
    p.add_argument("--misfit-max", type=float, help="override stitch.misfit_max")
    p.add_argument("--nbins", type=int, help="override stitch.nbins")
    p.add_argument("--smooth", type=float, help="override stitch.smooth")
    p.add_argument("--contour-step", type=float, help="override stitch.contour_step")
    p.add_argument("--taper-slope", type=float, help="override stitch.taper_slope")
    p.add_argument("--cmap", help="override stitch.cmap")
    p.add_argument("--vmin", type=float, help="override stitch.vmin")
    p.add_argument("--vmax", type=float, help="override stitch.vmax")
    p.add_argument("--station-csv",
                   help="override paths.station_csv used for the topography profile")
    p.add_argument("--output-dir",
                   help="override paths.output_dir, used to load models and write outputs")
    p.add_argument("--topography", dest="topography", action="store_true",
                   help="show the station-elevation profile above the section")
    p.add_argument("--no-topography", dest="topography", action="store_false",
                   help="hide the station-elevation profile")
    p.set_defaults(topography=None)
    p.add_argument("--output", help="output path (default: <output_dir>/section.png)")
    p.add_argument("--relative-output",
                   help="relative section output path (default: <output_dir>/relative_velocity.png)")
    p.add_argument("--relative-limit", type=float,
                   help="symmetric relative-velocity color limit in percent")
    args = p.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    output_dir = args.output_dir or config["paths"]["output_dir"]
    config["paths"]["output_dir"] = output_dir
    st = config["stitch"]

    def pick(cli, key):
        return cli if cli is not None else st.get(key)

    misfit_max = pick(args.misfit_max, "misfit_max")
    zmax = float(config["inversion"].get("zmax", np.inf))
    nbins = int(pick(args.nbins, "nbins"))
    cmap = pick(args.cmap, "cmap")
    vmin_cfg = pick(args.vmin, "vmin")
    vmax_cfg = pick(args.vmax, "vmax")
    smooth = float(pick(args.smooth, "smooth"))
    smooth_bands = None if args.smooth is not None else st.get("depth_smoothing")
    contour_step = float(pick(args.contour_step, "contour_step"))
    taper_slope = float(pick(args.taper_slope, "taper_slope") or 0)
    show_stations = bool(st.get("show_stations", True))
    title = st.get("title") or "Pseudo-2D S-wave velocity section"
    dpi = int(st.get("dpi", 300))
    save_pdf = bool(st.get("save_pdf", False))
    station_csv = args.station_csv or config["paths"].get("station_csv")
    topography = (args.topography if args.topography is not None
                  else bool(st.get("topography", station_csv is not None)))
    elevation_min = st.get("topography_elevation_min")
    if elevation_min is not None:
        elevation_min = float(elevation_min)
    if topography and not station_csv:
        raise SystemExit(
            "Topography is enabled, but paths.station_csv is not configured.")
    if topography and not os.path.isfile(station_csv):
        raise SystemExit(f"Station CSV not found: {station_csv}")
    out_path = args.output if args.output else os.path.join(output_dir, "section.png")
    if os.path.isdir(out_path) or out_path.endswith(os.sep):
        out_path = os.path.join(out_path, "section.png")

    recs, depth = load_models(output_dir, misfit_max)
    if not recs:
        raise SystemExit("No models loaded (check paths.output_dir / stitch.misfit_max).")
    if np.isfinite(zmax):
        depth_mask = depth <= zmax + 1e-9
        if not np.any(depth_mask):
            raise SystemExit(f"No model samples within configured zmax={zmax} km.")
        depth = depth[depth_mask]
        recs = [(lon, lat, misfit, vs[depth_mask]) for lon, lat, misfit, vs in recs]
    lons = [r[0] for r in recs]
    lats = [r[1] for r in recs]
    vs_all = np.array([r[3] for r in recs])           # (nmodels, ndepth)
    dist, line_transform = project_to_line(lons, lats)
    print(f"Loaded {len(recs)} models | distance 0-{dist.max():.2f} km | depth 0-{depth.max():.2f} km")

    topo_dist = topo_elevation = None
    if topography:
        try:
            topo_dist, topo_elevation = load_topography(
                station_csv, line_transform, elevation_min=elevation_min)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc
        print(f"Loaded {len(topo_dist)} station elevations | "
              f"range {topo_elevation.min():.0f}-{topo_elevation.max():.0f} m")

    # bin along distance and average the Vs columns in each bin
    edges = np.linspace(0, dist.max(), nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    section = np.full((len(depth), nbins), np.nan)
    counts = np.zeros(nbins, dtype=int)
    idx = np.clip(np.digitize(dist, edges) - 1, 0, nbins - 1)
    for b in range(nbins):
        sel = idx == b
        if np.any(sel):
            section[:, b] = vs_all[sel].mean(axis=0)
            counts[b] = sel.sum()
    filled = int((counts > 0).sum())
    print(f"Bins filled: {filled}/{nbins}")

    vmin = vmin_cfg if vmin_cfg is not None else np.nanpercentile(section, 2)
    vmax = vmax_cfg if vmax_cfg is not None else np.nanpercentile(section, 98)

    # fill empty bins, then smooth for a continuous-looking section
    filled = fill_nan_rows(section, centers)
    smoothed = apply_depth_smoothing(
        filled, depth, default_sigma=smooth, bands=smooth_bands)
    print(f"Smoothing: {describe_smoothing(smooth, smooth_bands)}")

    fig, ax = plt.subplots(figsize=(12, 6.75))
    Dmax = float(centers.max())
    if topography:
        # Each model depth is relative to the local ground surface.  Shift every
        # section column by its interpolated station elevation so the upper edge
        # of the Vs image is the actual terrain profile.
        surface = np.interp(centers, topo_dist, topo_elevation) / 1000.0
        plot_x, plot_y, surface = section_geometry(
            centers, depth, topo_dist, topo_elevation)
    else:
        plot_x, plot_y, surface = section_geometry(centers, depth)

    pcm = ax.pcolormesh(plot_x, plot_y, smoothed, cmap=cmap,
                        vmin=vmin, vmax=vmax, shading="gouraud")
    # velocity contour lines with labels
    if contour_step > 0:
        lo = np.floor(vmin / contour_step) * contour_step
        hi = np.ceil(vmax / contour_step) * contour_step
        levels = np.arange(lo, hi + 1e-9, contour_step)
        cs = ax.contour(plot_x, plot_y, smoothed, levels=levels,
                        colors="k", linewidths=0.5, alpha=0.7)
        clabels = ax.clabel(cs, fmt="%.2f", fontsize=6, inline=True)
    # mask into an inverted trapezoid (wide top, narrow bottom): edges become
    # less reliable with depth. Clip the image and contours to the trapezoid.
    if taper_slope > 0:
        verts = taper_vertices(centers, depth, taper_slope, surface)
        apply_taper(ax, [pcm, cs if contour_step > 0 else None], verts)
        if contour_step > 0:
            # clabel makes separate Text objects that set_clip_path doesn't touch;
            # hide any contour label whose anchor falls outside the trapezoid.
            from matplotlib.path import Path as MplPath
            poly = MplPath(verts)
            for t in clabels:
                if not poly.contains_point(t.get_position()):
                    t.set_visible(False)

    if topography:
        ax.plot(centers, surface, color="k", lw=1.2, zorder=5)
        if show_stations:
            station_view = (topo_dist >= 0) & (topo_dist <= Dmax)
            ax.scatter(topo_dist[station_view], topo_elevation[station_view] / 1000.0,
                       marker="v", s=10, facecolor="k", edgecolor="none",
                       zorder=6, clip_on=True)

    finish_section_axes(ax, centers, depth, surface, line_transform, title)
    data_min, data_max = float(np.nanmin(smoothed)), float(np.nanmax(smoothed))
    extend = ("both" if data_min < vmin and data_max > vmax else
              "min" if data_min < vmin else "max" if data_max > vmax else "neither")
    cb = fig.colorbar(pcm, ax=ax, pad=0.025, extend=extend)
    cb.set_label(r"$V_S$ (km s$^{-1}$)")
    cb.ax.tick_params(labelsize=8)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    print(f"Saved {out_path}")
    if save_pdf and os.path.splitext(out_path)[1].lower() != ".pdf":
        pdf_path = os.path.splitext(out_path)[0] + ".pdf"
        fig.savefig(pdf_path, bbox_inches="tight")
        print(f"Saved {pdf_path}")

    # also save the gridded section for further use
    npz = os.path.splitext(out_path)[0] + ".npz"
    saved = dict(distance=centers, depth=depth, vs=section, counts=counts)
    if topo_dist is not None:
        saved.update(topography_distance=topo_dist,
                     topography_elevation=topo_elevation)
    saved.update(axis_labels=np.asarray(profile_end_labels(line_transform)))
    np.savez(npz, **saved)
    print(f"Saved {npz}")
    plt.close(fig)

    rel_out = args.relative_output
    if rel_out is None and args.output:
        root, ext = os.path.splitext(out_path)
        rel_out = f"{root}_relative{ext or '.png'}"
    plot_relative_section(centers, depth, smoothed, config,
                          line_transform=line_transform,
                          topo_dist=topo_dist,
                          topo_elevation=topo_elevation,
                          output_path=rel_out,
                          limit=args.relative_limit,
                          save_pdf=save_pdf)


if __name__ == "__main__":
    main()
