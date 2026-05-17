#!/usr/bin/env python3

"""Reformat picked dispersion curves for DSurfTomo.

Supported inputs
1. DisperPicker post-processed results:
   - phase: <input_dir>/phase_velocity/STA1.STA2.txt
   - group: <input_dir>/group_velocity/STA1.STA2.txt
   Requires a station file because the curve files only contain T-V pairs.

2. Manual picks from EGFAnalysisTimeFreq:
   - phase: CDisp.T.STA1_STA2.dat
   - group: GDisp.STA1_STA2.dat
   These files contain station coordinates in the first two rows.

Output
  - surfdata.dat in DSurfTomo format
  - surfdata.json metadata for standalone downstream scripts

Period sampling
  - Sampling is strict.
  - A curve sample is kept only when its period lies on the
    min_period + n * period_step grid within a small floating-point tolerance.
  - Off-grid samples are skipped instead of being snapped to the nearest grid point.
"""

import os
from collections import OrderedDict, defaultdict
import json

import numpy as np
from dsurftomo_utils import (
    load_station_file,
    metadata_path_for,
    prompt,
    prompt_choice,
    prompt_float,
    prompt_int,
)


class Args(object):
    pass


STRICT_PERIOD_TOL = 1e-6


def collect_inputs():
    args = Args()
    print("DSurfTomo extractnew.py")
    print("Convert picked dispersion curves to surfdata.dat")
    print("")
    args.input_dir = prompt("input directory")
    selected_input_type = prompt_choice(
        "input source (auto/DisperPicker/manual)",
        ["auto", "disperpicker", "manual"],
        default="auto",
    )
    if selected_input_type == "auto":
        args.input_type = infer_input_type(args.input_dir)
        print("Detected input source: %s" % args.input_type)
    else:
        args.input_type = selected_input_type
    args.curve_type = prompt_choice(
        "curve type (phase/group/both)", ["phase", "group", "both"], default="both"
    )
    args.station_file = None
    if args.input_type == "disperpicker":
        args.station_file = prompt("station coordinate file for DisperPicker station pairs", default="")
    args.period_step = prompt_float("strict period step (s, exact grid only)")
    args.min_period = None
    args.max_period = None
    args.wave_type = prompt_int("wave type wavetp (1=Love, 2=Rayleigh)", default="2")
    args.vel_type = None if args.curve_type == "both" else (0 if args.curve_type == "phase" else 1)
    args.output = prompt("output surfdata path", default="surfdata.dat")
    return args


def list_input_files(input_dir, curve_type):
    input_type = infer_input_type(input_dir)
    return list_input_files_by_type(input_dir, input_type, curve_type)


def infer_input_type(input_dir):
    phase_dir = os.path.join(input_dir, "phase_velocity")
    group_dir = os.path.join(input_dir, "group_velocity")
    has_subdirs = os.path.isdir(phase_dir) or os.path.isdir(group_dir)

    top_level_names = os.listdir(input_dir)
    has_manual_top = any(
        name.startswith("CDisp.T.") or name.startswith("GDisp.") for name in top_level_names
    )

    def subdir_has_manual_dat(search_dir):
        if not os.path.isdir(search_dir):
            return False
        for name in os.listdir(search_dir):
            if not name.endswith(".dat"):
                continue
            if name.startswith("CDisp.T.") or name.startswith("GDisp."):
                return True
        return False

    has_manual_subdir = subdir_has_manual_dat(phase_dir) or subdir_has_manual_dat(group_dir)
    has_manual = has_manual_top or has_manual_subdir

    if has_manual:
        return "manual"
    if has_subdirs:
        return "disperpicker"
    raise ValueError("Cannot infer input source from directory contents.")


def list_input_files_by_type(input_dir, input_type, curve_type):
    if input_type == "manual":
        prefix = "CDisp.T." if curve_type == "phase" else "GDisp."
        subdir = "phase_velocity" if curve_type == "phase" else "group_velocity"
        if os.path.basename(os.path.normpath(input_dir)) == subdir:
            search_dir = input_dir
        elif os.path.isdir(os.path.join(input_dir, subdir)):
            search_dir = os.path.join(input_dir, subdir)
        else:
            search_dir = input_dir
        return sorted(
            os.path.join(search_dir, name)
            for name in os.listdir(search_dir)
            if name.startswith(prefix) and name.endswith(".dat")
        )
    subdir = "phase_velocity" if curve_type == "phase" else "group_velocity"
    if os.path.basename(os.path.normpath(input_dir)) == subdir:
        search_dir = input_dir
    else:
        search_dir = os.path.join(input_dir, subdir)
    if not os.path.isdir(search_dir):
        raise ValueError("Missing directory: %s" % search_dir)
    return sorted(
        os.path.join(search_dir, name)
        for name in os.listdir(search_dir)
        if name.endswith(".txt")
    )


def pair_from_ml_filename(path):
    name = os.path.splitext(os.path.basename(path))[0]
    parts = name.split(".")
    if len(parts) != 2:
        raise ValueError("Unexpected ml file name: %s" % path)
    return parts[0], parts[1]


def pair_from_manual_filename(path, curve_type):
    name = os.path.basename(path)
    prefix = "CDisp.T." if curve_type == "phase" else "GDisp."
    pair = os.path.splitext(name[len(prefix) :])[0]
    parts = pair.split("_")
    if len(parts) != 2:
        raise ValueError("Unexpected manual file name: %s" % path)
    return parts[0], parts[1]


def period_to_index(period, min_period, period_step, nf):
    offset = (period - min_period) / period_step
    nearest = round(offset)
    if abs(offset - nearest) > STRICT_PERIOD_TOL:
        return None
    index = int(nearest) + 1
    if index < 1 or index > nf:
        return None
    expected = min_period + (index - 1) * period_step
    if abs(period - expected) > STRICT_PERIOD_TOL:
        return None
    return index


def detect_period_bounds_by_type(input_dir, input_type, curve_type):
    input_files = list_input_files_by_type(input_dir, input_type, curve_type)
    periods = []
    for path in input_files:
        if input_type == "manual":
            with open(path, "r") as fobj:
                lines = [line.strip() for line in fobj if line.strip()]
            for line in lines[2:]:
                parts = line.split()
                if len(parts) < 2:
                    continue
                period = float(parts[0])
                velocity = float(parts[1])
                if len(parts) >= 4:
                    try:
                        flag = int(round(float(parts[3])))
                    except ValueError:
                        flag = 1
                else:
                    flag = 1
                if flag == 1 and velocity > 0:
                    periods.append(period)
        else:
            data = np.loadtxt(path)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            for row in data:
                period = float(row[0])
                velocity = float(row[1])
                if velocity > 0:
                    periods.append(period)
    if not periods:
        raise ValueError("No valid periods found in %s (%s)" % (input_dir, curve_type))
    return round(min(periods), 10), round(max(periods), 10)


def read_ml_curve(path, stations, min_period, max_period, period_step, nf):
    sta1, sta2 = pair_from_ml_filename(path)
    if sta1 not in stations or sta2 not in stations:
        raise KeyError("Station not found in station file: %s or %s" % (sta1, sta2))

    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise ValueError("Expected at least 2 columns in %s" % path)

    values = OrderedDict()
    for row in data:
        period = float(row[0])
        velocity = float(row[1])
        if velocity <= 0 or period < min_period or period > max_period:
            continue
        period_index = period_to_index(period, min_period, period_step, nf)
        if period_index is None:
            continue
        snapped_period = round(min_period + (period_index - 1) * period_step, 10)
        values[snapped_period] = velocity

    lat1, lon1 = stations[sta1]
    lat2, lon2 = stations[sta2]
    return values, (lat1, lon1), (lat2, lon2)


def read_manual_curve(path, min_period, max_period, period_step, nf):
    with open(path, "r") as fobj:
        lines = [line.strip() for line in fobj if line.strip()]

    if len(lines) < 3:
        raise ValueError("Unexpected manual dispersion file format: %s" % path)

    sta1_parts = lines[0].split()
    sta2_parts = lines[1].split()
    if len(sta1_parts) < 2 or len(sta2_parts) < 2:
        raise ValueError("Unexpected coordinate rows in %s" % path)

    lon1, lat1 = float(sta1_parts[0]), float(sta1_parts[1])
    lon2, lat2 = float(sta2_parts[0]), float(sta2_parts[1])

    values = OrderedDict()
    for line in lines[2:]:
        parts = line.split()
        if len(parts) < 2:
            continue
        period = float(parts[0])
        velocity = float(parts[1])
        if len(parts) >= 4:
            try:
                flag = int(round(float(parts[3])))
            except ValueError:
                flag = 1
        else:
            flag = 1
        if flag != 1 or velocity <= 0 or period < min_period or period > max_period:
            continue
        period_index = period_to_index(period, min_period, period_step, nf)
        if period_index is None:
            continue
        snapped_period = round(min_period + (period_index - 1) * period_step, 10)
        values[snapped_period] = velocity

    return values, (lat1, lon1), (lat2, lon2)


def remap_period_indices(records):
    all_periods = sorted(
        {period for rec in records for period in rec["values"].keys()}
    )
    period_to_index = {period: idx + 1 for idx, period in enumerate(all_periods)}
    for rec in records:
        rec["values"] = OrderedDict(
            (period_to_index[period], velocity)
            for period, velocity in rec["values"].items()
        )
    return all_periods


def build_records(args):
    input_type = args.input_type

    if args.vel_type is None:
        vel_type = 0 if args.curve_type == "phase" else 1
    else:
        vel_type = args.vel_type

    if args.min_period is None or args.max_period is None:
        args.min_period, args.max_period = detect_period_bounds_by_type(
            args.input_dir, input_type, args.curve_type
        )
    if args.period_step <= 0:
        raise ValueError("period_step must be positive.")
    nf = int(round((args.max_period - args.min_period) / args.period_step)) + 1
    input_files = list_input_files_by_type(args.input_dir, input_type, args.curve_type)
    if not input_files:
        raise ValueError("No input files found.")

    stations = None
    if input_type == "disperpicker":
        if not args.station_file:
            args.station_file = prompt("station coordinate file for DisperPicker station pairs")
        if not args.station_file:
            raise ValueError("station coordinate file is required for DisperPicker input.")
        stations = load_station_file(args.station_file)

    records = []
    for path in input_files:
        if input_type == "disperpicker":
            values, sta1, sta2 = read_ml_curve(
                path,
                stations,
                args.min_period,
                args.max_period,
                args.period_step,
                nf,
            )
        else:
            values, sta1, sta2 = read_manual_curve(
                path,
                args.min_period,
                args.max_period,
                args.period_step,
                nf,
            )
        if not values:
            continue
        records.append(
            {
                "sta1": sta1,
                "sta2": sta2,
                "values": values,
                "path": path,
            }
        )
    used_periods = remap_period_indices(records)
    return records, len(used_periods), vel_type, input_type, used_periods


def write_outputs(records, wave_type, vel_type, output_path):
    grouped = defaultdict(list)
    all_stations = set()
    for rec in records:
        sta1 = rec["sta1"]
        sta2 = rec["sta2"]
        all_stations.add(sta1)
        all_stations.add(sta2)
        for period_index, velocity in rec["values"].items():
            grouped[(period_index, sta1)].append((sta2, velocity))

    nsrcout = 0
    max_receivers = 0
    sources_per_period = defaultdict(set)
    with open(output_path, "w") as fout:
        for key in sorted(grouped, key=lambda item: (item[0], item[1][0], item[1][1])):
            period_index, (lat1, lon1) = key
            receivers = grouped[key]
            receivers.sort(key=lambda item: (item[0][0], item[0][1]))
            nsrcout += 1
            max_receivers = max(max_receivers, len(receivers))
            sources_per_period[period_index].add((lat1, lon1))
            fout.write(
                "# %.6f %.6f %d %d %d\n"
                % (lat1, lon1, period_index, wave_type, vel_type)
            )
            for (lat2, lon2), velocity in receivers:
                fout.write("%.6f %.6f %.6f\n" % (lat2, lon2, velocity))

    max_sources_per_period = max([len(v) for v in sources_per_period.values()] or [0])
    recommended_nsrc = max(max_receivers, max_sources_per_period, len(all_stations))
    return nsrcout, recommended_nsrc


def write_outputs_multi(runs, output_path):
    total_nsrcout = 0
    global_max_receivers = 0
    global_max_sources_per_period = 0
    all_stations = set()
    with open(output_path, "w") as fout:
        for run in runs:
            grouped = defaultdict(list)
            for rec in run["records"]:
                sta1 = rec["sta1"]
                sta2 = rec["sta2"]
                all_stations.add(sta1)
                all_stations.add(sta2)
                for period_index, velocity in rec["values"].items():
                    grouped[(period_index, sta1)].append((sta2, velocity))
            sources_per_period = defaultdict(set)
            for key in sorted(grouped, key=lambda item: (item[0], item[1][0], item[1][1])):
                period_index, (lat1, lon1) = key
                receivers = grouped[key]
                receivers.sort(key=lambda item: (item[0][0], item[0][1]))
                total_nsrcout += 1
                global_max_receivers = max(global_max_receivers, len(receivers))
                sources_per_period[period_index].add((lat1, lon1))
                fout.write(
                    "# %.6f %.6f %d %d %d\n"
                    % (lat1, lon1, period_index, run["wave_type"], run["vel_type"])
                )
                for (lat2, lon2), velocity in receivers:
                    fout.write("%.6f %.6f %.6f\n" % (lat2, lon2, velocity))
            run_max_sources = max([len(v) for v in sources_per_period.values()] or [0])
            global_max_sources_per_period = max(global_max_sources_per_period, run_max_sources)

    recommended_nsrc = max(global_max_receivers, global_max_sources_per_period, len(all_stations))
    return total_nsrcout, recommended_nsrc


def build_metadata(args, runs, recommended_nsrc):
    metadata = OrderedDict()
    metadata["input_dir"] = args.input_dir
    metadata["input_type"] = runs[0]["input_type"] if runs else args.input_type
    metadata["station_file"] = args.station_file
    metadata["wave_type"] = args.wave_type
    metadata["wave_type_name"] = "Rayleigh" if args.wave_type == 2 else "Love"
    metadata["curve_type"] = args.curve_type
    metadata["curve_types"] = [run["curve_type"] for run in runs]
    all_used_periods = [period for run in runs for period in run["used_periods"]]
    if all_used_periods:
        metadata["min_period"] = min(all_used_periods)
        metadata["max_period"] = max(all_used_periods)
    else:
        metadata["min_period"] = args.min_period
        metadata["max_period"] = args.max_period
    metadata["period_step"] = args.period_step
    metadata["period_count"] = len(sorted(set(all_used_periods))) if all_used_periods else 0
    metadata["output_surfdata"] = args.output
    metadata["output_metadata"] = metadata_path_for(args.output)
    metadata["recommended_nsrc"] = recommended_nsrc
    metadata["datasets"] = []
    for run in runs:
        item = OrderedDict()
        item["curve_type"] = run["curve_type"]
        item["vel_type"] = run["vel_type"]
        item["vel_type_name"] = "phase" if run["vel_type"] == 0 else "group"
        item["input_file_count"] = len(run["records"])
        item["used_period_indices"] = list(range(1, len(run["used_periods"]) + 1))
        item["used_periods"] = [round(period, 10) for period in run["used_periods"]]
        metadata["datasets"].append(item)
    return metadata


def write_metadata(metadata):
    path = metadata["output_metadata"]
    with open(path, "w") as fobj:
        json.dump(metadata, fobj, indent=2, ensure_ascii=True)
    return path


def main():
    args = collect_inputs()
    print("Strict sampling enabled: only exact period-step grid points are kept.")
    curve_types = ["phase", "group"] if args.curve_type == "both" else [args.curve_type]
    runs = []

    for curve_type in curve_types:
        run_args = Args()
        run_args.input_dir = args.input_dir
        run_args.input_type = args.input_type
        run_args.curve_type = curve_type
        run_args.station_file = args.station_file
        run_args.min_period = args.min_period
        run_args.max_period = args.max_period
        run_args.period_step = args.period_step
        run_args.wave_type = args.wave_type
        run_args.vel_type = None
        run_args.output = args.output

        records, nf, vel_type, input_type, used_periods = build_records(run_args)
        runs.append(
            {
                "curve_type": curve_type,
                "records": records,
                "nf": nf,
                "min_period": run_args.min_period,
                "max_period": run_args.max_period,
                "used_periods": used_periods,
                "vel_type": vel_type,
                "input_type": input_type,
                "wave_type": run_args.wave_type,
            }
        )

    if args.curve_type == "both":
        total_nsrcout, recommended_nsrc = write_outputs_multi(
            runs=runs, output_path=args.output
        )
        metadata_file = write_metadata(build_metadata(args, runs, recommended_nsrc))
        for run in runs:
            print(
                "Converted %d files (%s, %s). nf=%d, period range=%.3f-%.3f s"
                % (
                    len(run["records"]),
                    run["input_type"],
                    run["curve_type"],
                    run["nf"],
                    run["min_period"],
                    run["max_period"],
                )
            )
        print("Output: %s" % args.output)
        print("Recommended line 8 max(sources, receivers): %d" % recommended_nsrc)
        print("Wrote %s" % metadata_file)
    else:
        run = runs[0]
        _, recommended_nsrc = write_outputs(
            records=run["records"],
            wave_type=run["wave_type"],
            vel_type=run["vel_type"],
            output_path=args.output,
        )
        metadata_file = write_metadata(build_metadata(args, runs, recommended_nsrc))
        print(
            "Converted %d files (%s, %s). nf=%d, period range=%.3f-%.3f s"
            % (
                len(run["records"]),
                run["input_type"],
                run["curve_type"],
                run["nf"],
                run["min_period"],
                run["max_period"],
            )
        )
        print("Output: %s" % args.output)
        print("Recommended line 8 max(sources, receivers): %d" % recommended_nsrc)
        print("Wrote %s" % metadata_file)


if __name__ == "__main__":
    main()
