#!/usr/bin/env python3

"""Generate an initial DSurfTomo MOD file from surfdata.dat."""

import json
import os

import numpy as np
from dsurftomo_utils import (
    compute_grid_counts_from_spacing,
    compute_region_defaults,
    infer_default_nx_ny,
    metadata_path_for,
    prompt,
    prompt_choice,
    prompt_depths_with_arithmetic,
    prompt_float,
    prompt_int,
)


def load_surfdata(path):
    records = []
    current = None
    with open(path, "r") as fobj:
        for raw_line in fobj:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                parts = line.split()
                if len(parts) < 6:
                    raise ValueError("Invalid header line: %s" % raw_line.rstrip())
                current = {
                    "sta1_lat": float(parts[1]),
                    "sta1_lon": float(parts[2]),
                    "period_index": int(parts[3]),
                    "wave_type": int(parts[4]),
                    "vel_type": int(parts[5]),
                }
                continue
            if current is None:
                raise ValueError("Receiver line encountered before header.")
            parts = line.split()
            if len(parts) < 3:
                raise ValueError("Invalid receiver line: %s" % raw_line.rstrip())
            record = current.copy()
            record.update(
                {
                    "sta2_lat": float(parts[0]),
                    "sta2_lon": float(parts[1]),
                    "velocity": float(parts[2]),
                }
            )
            records.append(record)
    return records


def load_surfdata_metadata(surfdata_path):
    path = metadata_path_for(surfdata_path)
    if not os.path.exists(path):
        return None
    with open(path, "r") as fobj:
        return json.load(fobj)


def build_period_map(mode, records):
    unique_indices = sorted({rec["period_index"] for rec in records})
    if not unique_indices:
        raise ValueError("No period indices found in surfdata.")

    if mode == "range":
        min_period = prompt_float("minimum period in surfdata (s)")
        max_period = prompt_float("maximum period in surfdata (s)")
        period_step = prompt_float("period step in surfdata (s)")
        period_map = {}
        for idx in unique_indices:
            period_map[idx] = min_period + (idx - 1) * period_step
        expected_max = period_map[unique_indices[-1]]
        print(
            "Mapped period indices using min_period=%.3f, step=%.3f, last_period=%.3f"
            % (min_period, period_step, expected_max)
        )
        if abs(expected_max - max_period) > max(period_step * 0.51, 1e-6):
            print(
                "Warning: provided max period %.3f does not match mapped last period %.3f"
                % (max_period, expected_max)
            )
        return period_map

    periods = []
    for idx in unique_indices:
        periods.append(prompt_float("period for index %d (s)" % idx))
    return dict(zip(unique_indices, periods))


def build_period_map_from_metadata(records, metadata):
    unique_indices = sorted({rec["period_index"] for rec in records})
    for dataset in metadata.get("datasets", []):
        used_indices = dataset.get("used_period_indices", [])
        used_periods = dataset.get("used_periods", [])
        if len(used_indices) == len(used_periods) and used_indices:
            mapping = {int(idx): float(period) for idx, period in zip(used_indices, used_periods)}
            if all(idx in mapping for idx in unique_indices):
                print(
                    "Using exact period mapping from %s: first_period=%.3f s, last_period=%.3f s"
                    % (
                        os.path.basename(metadata_path_for(metadata["output_surfdata"])),
                        mapping[unique_indices[0]],
                        mapping[unique_indices[-1]],
                    )
                )
                return mapping

    min_period = float(metadata["min_period"])
    period_step = float(metadata["period_step"])
    period_map = {}
    for idx in unique_indices:
        period_map[idx] = min_period + (idx - 1) * period_step
    max_period = min_period + (unique_indices[-1] - 1) * period_step
    print(
        "Using metadata range from %s: min_period=%.3f s, step=%.3f s, last_period=%.3f s"
        % (os.path.basename(metadata_path_for(metadata["output_surfdata"])), min_period, period_step, max_period)
    )
    return period_map


def mean_velocity_at_period(records, period_map, target_period, wave_type, vel_type, prefer=None):
    grouped = {}
    for rec in records:
        if rec["wave_type"] != wave_type or rec["vel_type"] != vel_type:
            continue
        period = period_map.get(rec["period_index"])
        if period is None:
            continue
        grouped.setdefault(period, []).append(rec["velocity"])
    if not grouped:
        raise ValueError(
            "No measurements found for wave_type=%d, vel_type=%d"
            % (wave_type, vel_type)
        )
    available_periods = sorted(grouped)
    if prefer == "lower":
        candidates = [value for value in available_periods if value <= target_period + 1e-9]
        if candidates:
            chosen_period = max(candidates)
        else:
            chosen_period = min(available_periods, key=lambda value: abs(value - target_period))
    elif prefer == "upper":
        candidates = [value for value in available_periods if value >= target_period - 1e-9]
        if candidates:
            chosen_period = min(candidates)
        else:
            chosen_period = min(available_periods, key=lambda value: abs(value - target_period))
    else:
        chosen_period = min(available_periods, key=lambda value: abs(value - target_period))
    values = grouped[chosen_period]
    return float(np.mean(values)), len(values), chosen_period


def generate_mod(nx, ny, dep1, minvel, velgrad, output_path):
    nz = len(dep1)
    mod = np.zeros((nz * ny, nx))
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                mod[k * ny + j, i] = minvel + dep1[k] * velgrad

    with open(output_path, "w") as fp:
        for depth in dep1:
            fp.write("%5.1f" % depth)
        fp.write("\n")
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    fp.write("%7.3f" % mod[k * ny + j, i])
                fp.write("\n")
    return mod


def main():
    print("GenerateIniMOD.py")
    print("Generate initial MOD from surfdata.dat")
    print("")

    surfdata_path = prompt("surfdata path", default="surfdata.dat")
    if not os.path.exists(surfdata_path):
        raise ValueError("File not found: %s" % surfdata_path)

    records = load_surfdata(surfdata_path)
    metadata = load_surfdata_metadata(surfdata_path)

    if metadata is not None:
        period_map = build_period_map_from_metadata(records, metadata)
    else:
        mapping_mode = prompt_choice(
            "period mapping mode (range/manual)", ["range", "manual"], default="range"
        )
        period_map = build_period_map(mapping_mode, records)

    if metadata is not None:
        default_wave_type = str(metadata.get("wave_type", 2))
    else:
        default_wave_type = "2"
    print("Wave type in surfdata: Love=1, Rayleigh=2")
    wave_type = prompt_int("wave type used to estimate initial model", default=default_wave_type)

    if metadata is not None:
        curve_types = metadata.get("curve_types", [])
        default_vel_selection = "phase" if "phase" in curve_types else "group"
    else:
        default_vel_selection = "phase"
    print("Velocity type in surfdata: phase=0, group=1")
    selection = prompt_choice(
        "velocity type used to estimate initial model (phase/group)",
        ["phase", "group"],
        default=default_vel_selection,
    )
    vel_type = 0 if selection == "phase" else 1

    period1 = prompt_float("short period T1 (s)", default="1.0")
    period2 = prompt_float("long period T2 (s)", default="3.0")
    if period2 <= period1:
        raise ValueError("T2 must be greater than T1.")

    vs_1, n1, used_period1 = mean_velocity_at_period(
        records, period_map, period1, wave_type, vel_type, prefer="lower"
    )
    vs_2, n2, used_period2 = mean_velocity_at_period(
        records, period_map, period2, wave_type, vel_type, prefer="upper"
    )

    d1 = vs_1 * 0.92 * used_period1 / 3.0
    d2 = vs_2 * 0.92 * used_period2 / 3.0
    if d2 <= d1:
        raise ValueError(
            "Computed depth bounds are invalid: d1=%.3f, d2=%.3f "
            "(used T1=%.3f s, T2=%.3f s). Try a smaller T1, a larger T2, or manual mode."
            % (d1, d2, used_period1, used_period2)
        )

    velgrad = (vs_2 - vs_1) / (d2 - d1)
    minvel = vs_1 - velgrad * d1

    print("")
    print(
        "Average velocity near T1=%.3f s uses T=%.3f s: %.4f km/s from %d records"
        % (period1, used_period1, vs_1, n1)
    )
    print(
        "Average velocity near T2=%.3f s uses T=%.3f s: %.4f km/s from %d records"
        % (period2, used_period2, vs_2, n2)
    )
    if n1 < 10:
        print("Warning: T1 estimate uses only %d records, which is not very stable." % n1)
    if n2 < 10:
        print("Warning: T2 estimate uses only %d records, which is not very stable." % n2)
    print("Depth bounds from empirical rule: d1=%.4f km, d2=%.4f km" % (d1, d2))
    print("Recommended minvel=%.4f km/s, velgrad=%.4f km/s/km" % (minvel, velgrad))
    if velgrad <= 0:
        print("Warning: recommended velgrad is non-positive. Check period selection or use manual override.")
    print("")

    points = sorted(
        {
            (rec["sta1_lat"], rec["sta1_lon"])
            for rec in records
        }
        | {
            (rec["sta2_lat"], rec["sta2_lon"])
            for rec in records
        }
    )
    nx_default, ny_default, grid_info = infer_default_nx_ny(points)

    grid_mode = prompt_choice(
        "grid parameter mode (count/spacing)", ["count", "spacing"], default="count"
    )
    if grid_mode == "count":
        print("Grid mode: user specifies nx and ny; dvxd and dvzd are computed.")
        nx = prompt_int("nx", str(nx_default))
        ny = prompt_int("ny", str(ny_default))
        dvxd, dvzd, _, _ = compute_region_defaults(grid_info, nx, ny)
    else:
        dvxd_default, dvzd_default, _, _ = compute_region_defaults(grid_info, nx_default, ny_default)
        print("Grid mode: user specifies dvxd and dvzd; nx and ny are computed.")
        dvxd = prompt_float("dvxd", "%.5f" % dvxd_default)
        dvzd = prompt_float("dvzd", "%.5f" % dvzd_default)
        nx, ny, _, _ = compute_grid_counts_from_spacing(grid_info, dvxd, dvzd)

    print("Grid parameters:")
    print("  nx=%d ny=%d" % (nx, ny))
    print("  dvxd=%.5f dvzd=%.5f" % (dvxd, dvzd))
    dep1 = np.array(
        prompt_depths_with_arithmetic(
            "depth array dep1 in km", "[0,0.2,0.4,0.6,0.8,1.1,1.4,1.8,2.5]"
        ),
        dtype=float,
    )

    override = prompt_choice(
        "use recommended minvel/velgrad? (yes/no)", ["yes", "no"], default="yes"
    )
    if override == "no":
        minvel = prompt_float("minvel (km/s)", default="%.4f" % minvel)
        velgrad = prompt_float("velgrad (km/s/km)", default="%.4f" % velgrad)

    output_path = prompt("output MOD path", default="MOD")
    mod = generate_mod(nx, ny, dep1, minvel, velgrad, output_path)

    print("")
    print("Wrote %s" % output_path)
    print("nx=%d ny=%d nz=%d" % (nx, ny, len(dep1)))
    print("minvel=%.4f velgrad=%.4f" % (minvel, velgrad))
    print("top velocity=%.4f km/s" % mod[0, 0])
    print("bottom velocity=%.4f km/s" % mod[(len(dep1) - 1) * ny, 0])


if __name__ == "__main__":
    main()
