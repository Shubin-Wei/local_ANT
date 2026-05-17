#!/usr/bin/env python3

"""Shared utilities for DSurfTomo scripts."""

import ast
import math
import os


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BOUNDARY_EPS_FRACTION = 0.01


# ---------------------------------------------------------------------------
# Interactive prompts
# ---------------------------------------------------------------------------

def prompt(label, default=None):
    if default is None:
        raw = input("%s: " % label).strip()
    else:
        raw = input("%s [%s]: " % (label, default)).strip()
    if raw:
        return raw
    if default is None:
        raise ValueError("Input required: %s" % label)
    return str(default)


def prompt_choice(label, choices, default=None):
    normalized = [item.lower() for item in choices]
    alias_map = {}
    if len(choices) == 2:
        for item in choices:
            alias_map[item[0].lower()] = item.lower()

    while True:
        value = prompt(label, default).strip().lower()
        value = alias_map.get(value, value)
        if value in normalized:
            return value
        print("Invalid choice. Allowed values: %s" % ", ".join(choices))


def prompt_float(label, default=None):
    while True:
        value = prompt(label, default)
        try:
            return float(value)
        except ValueError:
            print("Please enter a number.")


def prompt_int(label, default=None):
    while True:
        value = prompt(label, default)
        try:
            return int(value)
        except ValueError:
            print("Please enter an integer.")


def _normalize_brackets(text):
    return (
        text.replace("【", "[")
        .replace("】", "]")
        .replace("（", "(")
        .replace("）", ")")
        .replace("，", ",")
    )


def prompt_depths(label, default=None):
    while True:
        text = _normalize_brackets(prompt(label, default))
        try:
            value = ast.literal_eval(text)
        except Exception:
            print("Please enter a Python-style list, for example [0,0.2,0.4].")
            continue
        if not isinstance(value, (list, tuple)) or not value:
            print("Please enter a non-empty list.")
            continue
        try:
            values = [float(item) for item in value]
        except ValueError:
            print("All depth values must be numeric.")
            continue
        if any(values[i] >= values[i + 1] for i in range(len(values) - 1)):
            print("Depth values must be strictly increasing.")
            continue
        return values


def prompt_depths_with_arithmetic(label, default=None):
    mode = prompt_choice("depth input mode (list/arithmetic)", ["list", "arithmetic"], "list")
    if mode == "list":
        return prompt_depths(label, default)

    start = prompt_float("depth start (km)", "0")
    end = prompt_float("depth end (km)")
    interval = prompt_float("depth interval (km)")
    if interval <= 0:
        raise ValueError("depth interval must be positive.")
    if end < start:
        raise ValueError("depth end must be greater than or equal to depth start.")

    values = []
    current = start
    eps = interval * 1e-6
    while current <= end + eps:
        values.append(round(current, 10))
        current += interval
    if not values:
        raise ValueError("Failed to generate depth array.")
    return values


# ---------------------------------------------------------------------------
# File I/O
# ---------------------------------------------------------------------------

def load_station_file(path):
    stations = {}
    with open(path, "r") as fobj:
        for raw_line in fobj:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            name = parts[0]
            first = float(parts[1])
            second = float(parts[2])
            if abs(first) > 90 and abs(second) <= 90:
                lon = first
                lat = second
            elif abs(second) > 90 and abs(first) <= 90:
                lat = first
                lon = second
            else:
                lat = first
                lon = second
            stations[name] = (lat, lon)
    if not stations:
        raise ValueError("No stations found in %s" % path)
    return stations


def metadata_path_for(data_path):
    base, _ = os.path.splitext(data_path)
    return base + ".json"


# ---------------------------------------------------------------------------
# Grid / region computation
# ---------------------------------------------------------------------------

def compute_min_distance(points):
    if len(points) < 2:
        return 0.01
    min_dist = None
    for i in range(len(points)):
        lat1, lon1 = points[i]
        for j in range(i + 1, len(points)):
            lat2, lon2 = points[j]
            dist = ((lat1 - lat2) ** 2 + (lon1 - lon2) ** 2) ** 0.5
            if min_dist is None or dist < min_dist:
                min_dist = dist
    if min_dist is None or min_dist <= 0:
        return 0.01
    return min_dist


def infer_grid(points):
    lats = [pt[0] for pt in points]
    lons = [pt[1] for pt in points]
    lat_min, lat_max = min(lats), max(lats)
    lon_min, lon_max = min(lons), max(lons)
    lat_range = max(lat_max - lat_min, 0.01)
    lon_range = max(lon_max - lon_min, 0.01)
    min_dist = compute_min_distance(points)
    default_grid = math.ceil((min_dist / 3.0) * 1000.0) / 1000.0
    if default_grid <= 0:
        default_grid = 0.01
    return {
        "lat_min": lat_min,
        "lat_max": lat_max,
        "lon_min": lon_min,
        "lon_max": lon_max,
        "lat_range": lat_range,
        "lon_range": lon_range,
        "grid": default_grid,
    }


def compute_region_origin(grid_info):
    lat_eps = max(grid_info["lat_range"] * BOUNDARY_EPS_FRACTION, 1e-4)
    lon_eps = max(grid_info["lon_range"] * BOUNDARY_EPS_FRACTION, 1e-4)
    goxd = grid_info["lat_max"] + lat_eps
    gozd = grid_info["lon_min"] - lon_eps
    return lat_eps, lon_eps, goxd, gozd


def compute_region_defaults(grid_info, nx, ny):
    lat_eps, lon_eps, goxd, gozd = compute_region_origin(grid_info)
    dvxd = (grid_info["lat_range"] + 2.0 * lat_eps) / max(nx - 3, 1)
    dvzd = (grid_info["lon_range"] + 2.0 * lon_eps) / max(ny - 3, 1)
    return dvxd, dvzd, goxd, gozd


def compute_grid_counts_from_spacing(grid_info, dvxd, dvzd):
    if dvxd <= 0 or dvzd <= 0:
        raise ValueError("dvxd and dvzd must be positive.")
    _, _, goxd, gozd = compute_region_origin(grid_info)
    lat_eps = max(grid_info["lat_range"] * BOUNDARY_EPS_FRACTION, 1e-4)
    lon_eps = max(grid_info["lon_range"] * BOUNDARY_EPS_FRACTION, 1e-4)
    nx = int(((grid_info["lat_range"] + 2.0 * lat_eps) / dvxd) + 0.999999) + 3
    ny = int(((grid_info["lon_range"] + 2.0 * lon_eps) / dvzd) + 0.999999) + 3
    return max(nx, 4), max(ny, 4), goxd, gozd


def infer_default_nx_ny(points):
    grid_info = infer_grid(points)
    lat_span = grid_info["lat_range"]
    lon_span = grid_info["lon_range"]
    target_lat_grid = lat_span / 15.0
    target_lon_grid = lon_span / 15.0
    default_grid = max(grid_info["grid"], target_lat_grid, target_lon_grid, 0.005)
    nx = int((lat_span / default_grid) + 0.999999) + 3
    ny = int((lon_span / default_grid) + 0.999999) + 3
    return max(nx, 4), max(ny, 4), grid_info


# ---------------------------------------------------------------------------
# Velocity bounds
# ---------------------------------------------------------------------------

def infer_velocity_bounds(values):
    if not values:
        return 0.5, 3.0
    min_v = min(values)
    max_v = max(values)
    default_min = max(0.1, math.floor((min_v - 0.2) * 10.0) / 10.0)
    default_max = math.ceil((max_v + 0.2) * 10.0) / 10.0
    if default_max <= default_min:
        default_max = default_min + 1.0
    return default_min, default_max


def infer_velocity_bounds_from_mod(path):
    if not path or not os.path.exists(path):
        return None
    with open(path, "r") as fobj:
        lines = [line.strip() for line in fobj if line.strip()]
    if len(lines) < 2:
        return None
    values = []
    for line in lines[1:]:
        values.extend(float(x) for x in line.split())
    if not values:
        return None
    return infer_velocity_bounds(values)


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------

def print_banner(title, lines=None):
    width = max(72, len(title) + 8)
    border = "#" * width
    print("")
    print(border)
    print("## %s" % title)
    if lines:
        for line in lines:
            print("## %s" % line)
    print(border)


def print_section(title):
    width = max(72, len(title) + 8)
    border = "-" * width
    print("")
    print(border)
    print("[ %s ]" % title)
    print(border)


def print_kv(label, value):
    print("  - %s: %s" % (label, value))
