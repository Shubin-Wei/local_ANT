#!/usr/bin/env python3

"""Interactive generator for DSurfTomo.in from existing surfdata.dat and MOD."""

import json
import os

from dsurftomo_utils import (
    compute_grid_counts_from_spacing,
    compute_region_defaults,
    infer_velocity_bounds,
    infer_velocity_bounds_from_mod,
    load_station_file,
    metadata_path_for,
    prompt,
    prompt_choice,
    prompt_float,
    prompt_int,
)


def load_json_if_exists(path):
    if not path or not os.path.exists(path):
        return None
    with open(path, "r") as fobj:
        return json.load(fobj)


def load_surfdata(path):
    records = []
    all_points = set()
    source_blocks = 0
    current = None
    with open(path, "r") as fobj:
        for raw_line in fobj:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                parts = line.split()
                current = {
                    "sta1": (float(parts[1]), float(parts[2])),
                    "period_index": int(parts[3]),
                    "wave_type": int(parts[4]),
                    "vel_type": int(parts[5]),
                }
                all_points.add(current["sta1"])
                source_blocks += 1
                continue
            parts = line.split()
            sta2 = (float(parts[0]), float(parts[1]))
            all_points.add(sta2)
            records.append(
                {
                    "sta1": current["sta1"],
                    "sta2": sta2,
                    "velocity": float(parts[2]),
                    "period_index": current["period_index"],
                    "wave_type": current["wave_type"],
                    "vel_type": current["vel_type"],
                }
            )
    return records, sorted(all_points), source_blocks


def build_period_map_from_surfdata(records, metadata=None):
    unique_indices = sorted({rec["period_index"] for rec in records})
    if metadata is not None and "min_period" in metadata and "period_step" in metadata:
        for dataset in metadata.get("datasets", []):
            used_indices = dataset.get("used_period_indices", [])
            used_periods = dataset.get("used_periods", [])
            if len(used_indices) == len(used_periods) and used_indices:
                mapping = {int(idx): float(period) for idx, period in zip(used_indices, used_periods)}
                if all(idx in mapping for idx in unique_indices):
                    return mapping
        min_period = float(metadata["min_period"])
        step = float(metadata["period_step"])
        mapping = {}
        for idx in unique_indices:
            mapping[idx] = min_period + (idx - 1) * step
        return mapping
    print("surfdata.dat stores period indices, not actual periods.")
    mode = prompt_choice("period mapping mode (range/manual)", ["range", "manual"], default="range")
    mapping = {}
    if mode == "range":
        min_period = prompt_float("minimum period in surfdata (s)")
        max_period = prompt_float("maximum period in surfdata (s)")
        step = prompt_float("period step in surfdata (s)")
        for idx in unique_indices:
            mapping[idx] = min_period + (idx - 1) * step
        predicted_last = mapping[unique_indices[-1]]
        if abs(predicted_last - max_period) > max(step * 0.51, 1e-6):
            print(
                "Warning: last mapped period is %.3f s, different from supplied max %.3f s"
                % (predicted_last, max_period)
            )
    else:
        for idx in unique_indices:
            mapping[idx] = prompt_float("period for index %d (s)" % idx)
    return mapping


def infer_model_shape_from_mod(path):
    if not path or not os.path.exists(path):
        return None
    with open(path, "r") as fobj:
        lines = [line.strip() for line in fobj if line.strip()]
    if not lines:
        return None
    nz = len(lines[0].split())
    if nz <= 0:
        return None
    body = lines[1:]
    if not body:
        return None
    nx = len(body[0].split())
    if nx <= 0:
        return None
    if len(body) % nz != 0:
        return None
    ny = len(body) // nz
    return nx, ny, nz


def write_input_file(
    output_path,
    data_file_name,
    nx,
    ny,
    nz,
    goxd,
    gozd,
    dvxd,
    dvzd,
    nsrcmax,
    weight,
    damp,
    sublayers,
    minvel,
    maxvel,
    maxiteration,
    sparsityfraction,
    periods_rc,
    periods_rg,
    periods_lc,
    periods_lg,
    synthetic_flag,
    noiselevel,
    threshold,
):
    # Fortran list-directed reads (read(10,*) ...) take only the leading
    # value(s) on each line and ignore the rest, so trailing "c: ..." comments
    # are safe. The period lines carry no comment (read(10,*)(tRc(i),...) would
    # otherwise stop mid-line anyway), matching the manual's example .in.
    dummy = "c" * 69

    def line(value_str, comment):
        return "%-33s%s\n" % (value_str, "c: " + comment)

    with open(output_path, "w") as inputfile:
        inputfile.write(dummy + "\n")
        inputfile.write("c INPUT PARAMETERS\n")
        inputfile.write(dummy + "\n")
        inputfile.write(line("%s" % data_file_name, "data file"))
        inputfile.write(line("%d %d %d" % (nx, ny, nz),
                             "nx ny nz (grid number in lat lon and depth direction)"))
        inputfile.write(line("%.6f %.6f" % (goxd, gozd),
                             "goxd gozd (upper left point,[lat,lon])"))
        inputfile.write(line("%.6f %.6f" % (dvxd, dvzd),
                             "dvxd dvzd (grid interval in lat and lon direction)"))
        inputfile.write(line("%d" % nsrcmax, "max(sources, receivers)"))
        inputfile.write(line("%.3f %.3f" % (weight, damp), "weight damp"))
        inputfile.write(line("%d" % sublayers,
                             "sublayers (for computing depth kernel, 2~5)"))
        inputfile.write(line("%.3f %.3f" % (minvel, maxvel),
                             "minimum velocity, maximum velocity (a priori information)"))
        inputfile.write(line("%d" % maxiteration, "maximum iteration"))
        inputfile.write(line("%.3f" % sparsityfraction, "sparsity fraction"))
        kmax_comments = ["kmaxRc (followed by periods)", "kmaxRg", "kmaxLc", "kmaxLg"]
        for periods, comment in zip(
            [periods_rc, periods_rg, periods_lc, periods_lg], kmax_comments
        ):
            inputfile.write(line("%d" % len(periods), comment))
            if periods:
                inputfile.write(" ".join("%.3f" % value for value in periods) + "\n")
        inputfile.write(line("%d" % synthetic_flag,
                             "synthetic flag(0:real data,1:synthetic)"))
        inputfile.write(line("%.3f" % noiselevel, "noiselevel"))
        inputfile.write(line("%.3f" % threshold, "threshold"))


def main():
    from dsurftomo_utils import infer_grid

    print("GenerateDSurfTomoInputFile.py")
    print("Generate DSurfTomo.in from existing surfdata.dat, MOD and station.txt")
    print("")

    surfdata_path = prompt("existing surfdata path", "surfdata.dat")
    if not os.path.exists(surfdata_path):
        raise ValueError("File not found: %s" % surfdata_path)
    mod_path = prompt("MOD file path", "MOD")
    if not os.path.exists(mod_path):
        raise ValueError("File not found: %s" % mod_path)
    surfdata_metadata = load_json_if_exists(metadata_path_for(surfdata_path))
    station_file = None
    if surfdata_metadata is not None:
        station_file = surfdata_metadata.get("station_file")
    if station_file and os.path.exists(station_file):
        print("Using station coordinate file: %s" % station_file)
    else:
        fallback_station = "station.txt"
        if os.path.exists(fallback_station):
            station_file = fallback_station
            print("Using station coordinate file: %s" % station_file)
        else:
            station_file = prompt("station coordinate file")

    surf_records, points_from_surfdata, source_blocks = load_surfdata(surfdata_path)
    period_mapping = build_period_map_from_surfdata(surf_records, surfdata_metadata)
    period_lists = {
        ("rayleigh", "phase"): sorted({period_mapping[rec["period_index"]] for rec in surf_records if rec["wave_type"] == 2 and rec["vel_type"] == 0}),
        ("rayleigh", "group"): sorted({period_mapping[rec["period_index"]] for rec in surf_records if rec["wave_type"] == 2 and rec["vel_type"] == 1}),
        ("love", "phase"): sorted({period_mapping[rec["period_index"]] for rec in surf_records if rec["wave_type"] == 1 and rec["vel_type"] == 0}),
        ("love", "group"): sorted({period_mapping[rec["period_index"]] for rec in surf_records if rec["wave_type"] == 1 and rec["vel_type"] == 1}),
    }
    data_file_name = os.path.basename(surfdata_path)
    print("Using data file name in DSurfTomo.in: %s" % data_file_name)

    if station_file and os.path.exists(station_file):
        points = sorted(load_station_file(station_file).values())
    else:
        if station_file:
            print("Warning: station file not found, use coordinates from surfdata.dat instead.")
        points = points_from_surfdata

    grid_info = infer_grid(points)
    print("Region bounds: tight envelope from station coverage with small safety margin")
    model_shape = infer_model_shape_from_mod(mod_path)
    if model_shape is None:
        default_nx, default_ny, default_nz = 18, 18, 9
    else:
        default_nx, default_ny, default_nz = model_shape

    print("")
    print("Detected model shape from MOD:")
    print("  nx=%d ny=%d nz=%d" % (default_nx, default_ny, default_nz))
    use_model_shape = prompt_choice(
        "use detected model shape? (yes/no)", ["yes", "no"], default="yes"
    )
    if use_model_shape == "yes":
        nx, ny, nz = default_nx, default_ny, default_nz
    else:
        nx = prompt_int("nx", str(default_nx))
        ny = prompt_int("ny", str(default_ny))
        nz = prompt_int("nz", str(default_nz))

    default_dvxd, default_dvzd, default_goxd, default_gozd = compute_region_defaults(grid_info, nx, ny)

    print("")
    grid_mode = prompt_choice(
        "grid parameter mode (count/spacing)", ["count", "spacing"], default="count"
    )
    if grid_mode == "count":
        print("Grid mode: user specifies nx and ny; dvxd and dvzd are computed.")
        nx = prompt_int("nx", str(nx))
        ny = prompt_int("ny", str(ny))
        dvxd, dvzd, goxd, gozd = compute_region_defaults(grid_info, nx, ny)
    else:
        print("Grid mode: user specifies dvxd and dvzd; nx and ny are computed.")
        dvxd = prompt_float("dvxd", "%.5f" % default_dvxd)
        dvzd = prompt_float("dvzd", "%.5f" % default_dvzd)
        nx, ny, goxd, gozd = compute_grid_counts_from_spacing(grid_info, dvxd, dvzd)

    print("")
    print("Computed region/grid parameters:")
    print("  grid mode=%s" % grid_mode)
    print("  nx=%d ny=%d nz=%d" % (nx, ny, nz))
    print("  dvxd=%.5f dvzd=%.5f" % (dvxd, dvzd))
    print("  goxd=%.5f gozd=%.5f" % (goxd, gozd))

    default_nsrcmax = str(source_blocks if source_blocks is not None else 1)
    if surfdata_metadata is not None:
        if surfdata_metadata.get("recommended_nsrc"):
            default_nsrcmax = str(surfdata_metadata["recommended_nsrc"])
        elif surfdata_metadata.get("source_blocks"):
            default_nsrcmax = str(surfdata_metadata["source_blocks"])
    print("")
    print("Computed line-8 parameter:")
    print("  max(sources, receivers)=%s" % default_nsrcmax)
    use_nsrcmax = prompt_choice(
        "use computed line-8 value? (yes/no)", ["yes", "no"], default="yes"
    )
    if use_nsrcmax == "yes":
        nsrcmax = int(default_nsrcmax)
    else:
        nsrcmax = prompt_int("max(sources, receivers)", default_nsrcmax)

    weight = prompt_float("weight", "4.0")
    damp = prompt_float("damp", "1.0")
    sublayers = prompt_int("sublayers", "3")
    mod_bounds = infer_velocity_bounds_from_mod(mod_path)
    if mod_bounds is not None:
        default_minvel, default_maxvel = mod_bounds
    else:
        default_minvel, default_maxvel = infer_velocity_bounds(
            [rec["velocity"] for rec in surf_records]
        )
    print("")
    print("Computed velocity bounds:")
    print("  minimum velocity=%.3f maximum velocity=%.3f" % (default_minvel, default_maxvel))
    use_velocity_bounds = prompt_choice(
        "use computed velocity bounds? (yes/no)", ["yes", "no"], "yes"
    )
    if use_velocity_bounds == "yes":
        minvel, maxvel = default_minvel, default_maxvel
    else:
        minvel = prompt_float("minimum velocity", "%.3f" % default_minvel)
        maxvel = prompt_float("maximum velocity", "%.3f" % default_maxvel)
    maxiteration = prompt_int("maximum iteration", "10")
    sparsityfraction = prompt_float("sparsity fraction", "0.2")
    synthetic_flag = prompt_int("synthetic flag (0=real, 1=synthetic)", "0")
    noiselevel = prompt_float("noise level", "0.02")
    threshold = prompt_float("threshold", "0.5")

    input_output = prompt("output DSurfTomo.in path", "DSurfTomo.in")
    write_input_file(
        output_path=input_output,
        data_file_name=data_file_name,
        nx=nx,
        ny=ny,
        nz=nz,
        goxd=goxd,
        gozd=gozd,
        dvxd=dvxd,
        dvzd=dvzd,
        nsrcmax=nsrcmax,
        weight=weight,
        damp=damp,
        sublayers=sublayers,
        minvel=minvel,
        maxvel=maxvel,
        maxiteration=maxiteration,
        sparsityfraction=sparsityfraction,
        periods_rc=period_lists[("rayleigh", "phase")],
        periods_rg=period_lists[("rayleigh", "group")],
        periods_lc=period_lists[("love", "phase")],
        periods_lg=period_lists[("love", "group")],
        synthetic_flag=synthetic_flag,
        noiselevel=noiselevel,
        threshold=threshold,
    )

    print("")
    print("surfdata input: %s" % surfdata_path)
    print("MOD input: %s" % mod_path)
    if station_file:
        print("station file: %s" % station_file)
    print("DSurfTomo.in output: %s" % input_output)
    print("period counts: Rc=%d Rg=%d Lc=%d Lg=%d" % (
        len(period_lists[("rayleigh", "phase")]),
        len(period_lists[("rayleigh", "group")]),
        len(period_lists[("love", "phase")]),
        len(period_lists[("love", "group")]),
    ))
    print("grid: nx=%d ny=%d nz=%d" % (nx, ny, nz))


if __name__ == "__main__":
    main()
