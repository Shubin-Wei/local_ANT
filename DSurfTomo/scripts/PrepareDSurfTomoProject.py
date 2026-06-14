#!/usr/bin/env python3

"""Run extractnew, GenerateIniMOD and GenerateDSurfTomoInputFile in sequence.

Two ways to run it:
  * interactive  : `python3 PrepareDSurfTomoProject.py`  (asks questions)
  * command-line : `python3 PrepareDSurfTomoProject.py --input-dir ... --project-dir ...`
    Any argument switches it into the non-interactive `build_project()` path,
    which builds surfdata.dat, surfdata.json, MOD, (optional) MOD.true and
    DSurfTomo.in in one shot. See `--help` for all knobs.
"""

import argparse
import os
import sys

import numpy as np

import extractnew as ex
import GenerateDSurfTomoInputFile as gdi
import GenerateIniMOD as gim
import GenerateTrueMOD as gtm
from dsurftomo_utils import (
    check_sparse_matrix_size,
    compute_grid_counts_from_spacing,
    compute_region_defaults,
    infer_default_nx_ny,
    infer_grid,
    infer_velocity_bounds,
    infer_velocity_bounds_from_mod,
    print_banner,
    print_kv,
    print_section,
    prompt,
    prompt_choice,
    prompt_depths_with_arithmetic,
    prompt_float,
    prompt_int,
)


def infer_curve_type(input_dir, input_type):
    if input_type != "egfpy":
        raise ValueError("Unsupported input type: %s" % input_type)
    search_dir = ex.egfpy_search_dir(input_dir)
    names = os.listdir(search_dir)
    has_phase = any(name.startswith("CDisp.") and name.endswith(".txt") for name in names)
    has_group = any(name.startswith("GDisp.") and name.endswith(".txt") for name in names)

    if has_phase and has_group:
        return "both"
    if has_phase:
        return "phase"
    if has_group:
        return "group"
    raise ValueError("Cannot infer curve type from directory contents.")


def build_runs(args):
    curve_types = ["phase", "group"] if args.curve_type == "both" else [args.curve_type]
    runs = []
    for curve_type in curve_types:
        run_args = ex.Args()
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
        records, nf, vel_type, input_type, used_periods = ex.build_records(run_args)
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
                "wave_type": args.wave_type,
            }
        )
    return runs


def write_surfdata(args, runs):
    if args.curve_type == "both":
        _, recommended_nsrc = ex.write_outputs_multi(runs=runs, output_path=args.output)
    else:
        run = runs[0]
        _, recommended_nsrc = ex.write_outputs(
            records=run["records"],
            wave_type=run["wave_type"],
            vel_type=run["vel_type"],
            output_path=args.output,
        )
    metadata_path = ex.write_metadata(ex.build_metadata(args, runs, recommended_nsrc))
    return recommended_nsrc, metadata_path


def choose_station_points(station_file, surfdata_path):
    if station_file and os.path.exists(station_file):
        return sorted(gdi.load_station_file(station_file).values())
    surf_records, points_from_surfdata, _ = gdi.load_surfdata(surfdata_path)
    del surf_records
    return points_from_surfdata


def build_period_lists(surfdata_path):
    surf_records, _, _ = gdi.load_surfdata(surfdata_path)
    surfdata_metadata = gdi.load_json_if_exists(gdi.metadata_path_for(surfdata_path))
    period_mapping = gdi.build_period_map_from_surfdata(surf_records, surfdata_metadata)
    period_lists = {
        ("rayleigh", "phase"): sorted(
            {
                period_mapping[rec["period_index"]]
                for rec in surf_records
                if rec["wave_type"] == 2 and rec["vel_type"] == 0
            }
        ),
        ("rayleigh", "group"): sorted(
            {
                period_mapping[rec["period_index"]]
                for rec in surf_records
                if rec["wave_type"] == 2 and rec["vel_type"] == 1
            }
        ),
        ("love", "phase"): sorted(
            {
                period_mapping[rec["period_index"]]
                for rec in surf_records
                if rec["wave_type"] == 1 and rec["vel_type"] == 0
            }
        ),
        ("love", "group"): sorted(
            {
                period_mapping[rec["period_index"]]
                for rec in surf_records
                if rec["wave_type"] == 1 and rec["vel_type"] == 1
            }
        ),
    }
    return surf_records, period_lists


def _parse_depths(text):
    cleaned = text.strip().strip("[]").replace(",", " ")
    return np.array([float(item) for item in cleaned.split() if item], dtype=float)


def build_project(opts):
    """Non-interactive equivalent of main(): build a full DSurfTomo project."""
    input_dir = opts.input_dir
    input_type = ex.infer_input_type(input_dir)
    curve_type = opts.curve_type or infer_curve_type(input_dir, input_type)
    project_dir = opts.project_dir
    os.makedirs(project_dir, exist_ok=True)
    surfdata_path = os.path.join(project_dir, "surfdata.dat")

    # ---- Step 1: surfdata.dat (+ .json) ----
    args = ex.Args()
    args.input_dir = input_dir
    args.input_type = input_type
    args.curve_type = curve_type
    args.station_file = None
    args.min_period = opts.min_period
    args.max_period = opts.max_period
    args.period_step = opts.period_step
    args.wave_type = opts.wave_type
    args.vel_type = None
    args.output = surfdata_path
    runs = build_runs(args)
    recommended_nsrc, metadata_path = write_surfdata(args, runs)

    # ---- Step 2: initial MOD ----
    records = gim.load_surfdata(surfdata_path)
    metadata = gim.load_surfdata_metadata(surfdata_path)
    period_map = gim.build_period_map_from_metadata(records, metadata)
    points = choose_station_points(None, surfdata_path)
    grid_info = infer_grid(points)
    if opts.grid_spacing is not None:
        dvxd, dvzd = opts.grid_spacing
        nx, ny, goxd, gozd = compute_grid_counts_from_spacing(grid_info, dvxd, dvzd)
    elif opts.grid_count is not None:
        nx, ny = opts.grid_count
        dvxd, dvzd, goxd, gozd = compute_region_defaults(grid_info, nx, ny)
    else:
        nx, ny, _ = infer_default_nx_ny(points)
        dvxd, dvzd, goxd, gozd = compute_region_defaults(grid_info, nx, ny)

    dep1 = _parse_depths(opts.depths)

    # Fail fast if DSurfTomo's 32-bit sparse indexing would overflow for this
    # region/data volume (scales with dall*nx*ny*nz*sparsity).
    check_sparse_matrix_size(len(records), nx, ny, len(dep1), opts.sparsity)

    if opts.minvel is not None and opts.velgrad is not None:
        minvel, velgrad = opts.minvel, opts.velgrad
    else:
        model_vel_type = 0 if opts.model_vel == "phase" else 1
        vs1, _, ut1 = gim.mean_velocity_at_period(
            records, period_map, opts.t1, opts.wave_type, model_vel_type, prefer="lower")
        vs2, _, ut2 = gim.mean_velocity_at_period(
            records, period_map, opts.t2, opts.wave_type, model_vel_type, prefer="upper")
        minvel, velgrad, fellback, _, _ = gim.estimate_minvel_velgrad(vs1, ut1, vs2, ut2)
        if fellback:
            print("Warning: flat/decreasing dispersion; using fallback initial model "
                  "(minvel=%.4f, velgrad=%.4f). Pass --minvel/--velgrad to override."
                  % (minvel, velgrad))

    mod_path = os.path.join(project_dir, "MOD")
    gim.generate_mod(nx, ny, dep1, minvel, velgrad, mod_path)

    # ---- Step 2b: optional checkerboard MOD.true ----
    true_mod_path = None
    if opts.checkerboard:
        _, background_mod = gtm.load_mod(mod_path)
        true_mod = gtm.generate_true_mod(background_mod, anosize=opts.anosize, amplitude=opts.amplitude)
        true_mod_path = os.path.join(project_dir, "MOD.true")
        gtm.write_true_mod(true_mod_path, true_mod)

    # ---- Step 3: DSurfTomo.in ----
    surf_records, period_lists = build_period_lists(surfdata_path)
    nsrcmax = opts.nsrcmax if opts.nsrcmax else recommended_nsrc
    if opts.inv_minvel is not None and opts.inv_maxvel is not None:
        inv_minvel, inv_maxvel = opts.inv_minvel, opts.inv_maxvel
    else:
        bounds = infer_velocity_bounds_from_mod(mod_path)
        if bounds is not None:
            inv_minvel, inv_maxvel = bounds
        else:
            inv_minvel, inv_maxvel = infer_velocity_bounds([r["velocity"] for r in surf_records])
    synthetic_flag = opts.synthetic
    if synthetic_flag is None:
        synthetic_flag = 1 if true_mod_path is not None else 0

    input_path = os.path.join(project_dir, "DSurfTomo.in")
    gdi.write_input_file(
        output_path=input_path,
        data_file_name=os.path.basename(surfdata_path),
        nx=nx, ny=ny, nz=len(dep1), goxd=goxd, gozd=gozd, dvxd=dvxd, dvzd=dvzd,
        nsrcmax=nsrcmax, weight=opts.weight, damp=opts.damp, sublayers=opts.sublayers,
        minvel=inv_minvel, maxvel=inv_maxvel, maxiteration=opts.max_iter,
        sparsityfraction=opts.sparsity,
        periods_rc=period_lists[("rayleigh", "phase")],
        periods_rg=period_lists[("rayleigh", "group")],
        periods_lc=period_lists[("love", "phase")],
        periods_lg=period_lists[("love", "group")],
        synthetic_flag=synthetic_flag, noiselevel=opts.noise, threshold=opts.threshold,
    )

    print_section("Project Ready (non-interactive)")
    print_kv("project dir", project_dir)
    print_kv("curve type", curve_type)
    for run in runs:
        print_kv("%s periods" % run["curve_type"],
                 "%.2f-%.2f s (%d used)" % (run["min_period"], run["max_period"], run["nf"]))
    print_kv("grid nx ny nz", "%d %d %d" % (nx, ny, len(dep1)))
    print_kv("dvxd dvzd", "%.5f %.5f" % (dvxd, dvzd))
    print_kv("goxd gozd", "%.5f %.5f" % (goxd, gozd))
    print_kv("kmaxRc kmaxRg",
             "%d %d" % (len(period_lists[("rayleigh", "phase")]),
                        len(period_lists[("rayleigh", "group")])))
    print_kv("init minvel velgrad", "%.4f %.4f" % (minvel, velgrad))
    print_kv("inv vel bounds", "%.3f %.3f" % (inv_minvel, inv_maxvel))
    print_kv("nsrcmax (line 8)", nsrcmax)
    print_kv("weight damp sublayers", "%.2f %.2f %d" % (opts.weight, opts.damp, opts.sublayers))
    print_kv("max_iter sparsity", "%d %.2f" % (opts.max_iter, opts.sparsity))
    print_kv("synthetic noise threshold",
             "%d %.3f %.3f" % (synthetic_flag, opts.noise, opts.threshold))
    if true_mod_path is not None:
        print_kv("MOD.true", "anosize=%.3f amplitude=%.3f km/s" % (opts.anosize, opts.amplitude))
    print("Files written: surfdata.dat, surfdata.json, MOD%s, DSurfTomo.in"
          % (", MOD.true" if true_mod_path is not None else ""))
    return project_dir


def parse_cli(argv):
    p = argparse.ArgumentParser(
        description="Build a DSurfTomo project (surfdata/MOD/MOD.true/DSurfTomo.in) non-interactively.")
    p.add_argument("--input-dir", required=True,
                   help="EGFAnalysisPy dispersion dir (qc_result or its parent).")
    p.add_argument("--project-dir", required=True, help="output project directory.")
    p.add_argument("--curve-type", choices=["both", "phase", "group"], default=None,
                   help="default: auto-detect from files present.")
    p.add_argument("--period-step", type=float, required=True, help="strict period grid step (s).")
    p.add_argument("--min-period", type=float, default=None, help="default: auto-detect.")
    p.add_argument("--max-period", type=float, default=None, help="default: auto-detect.")
    p.add_argument("--wave-type", type=int, default=2, help="1=Love, 2=Rayleigh (default 2).")
    grid = p.add_mutually_exclusive_group()
    grid.add_argument("--grid-spacing", type=float, nargs=2, metavar=("DVXD", "DVZD"),
                      default=None, help="grid spacing in deg (lat lon); nx,ny computed.")
    grid.add_argument("--grid-count", type=int, nargs=2, metavar=("NX", "NY"),
                      default=None, help="grid counts; dvxd,dvzd computed.")
    p.add_argument("--depths", default="0,0.2,0.4,0.6,0.8,1.0,1.3,1.6,2.0,2.5",
                   help="comma/space list of depth nodes (km).")
    p.add_argument("--minvel", type=float, default=None, help="manual initial-model surface vel.")
    p.add_argument("--velgrad", type=float, default=None, help="manual initial-model vel gradient.")
    p.add_argument("--t1", type=float, default=0.5, help="short period for model estimate (s).")
    p.add_argument("--t2", type=float, default=1.8, help="long period for model estimate (s).")
    p.add_argument("--model-vel", choices=["phase", "group"], default="phase",
                   help="vel type used to estimate initial model.")
    p.add_argument("--nsrcmax", type=int, default=None, help="line 8; default: recommended.")
    p.add_argument("--weight", type=float, default=5.0, help="smoothing weight (default 5).")
    p.add_argument("--damp", type=float, default=1.0, help="damping (default 1).")
    p.add_argument("--sublayers", type=int, default=3)
    p.add_argument("--inv-minvel", type=float, default=None, help="a priori min vel; default from MOD.")
    p.add_argument("--inv-maxvel", type=float, default=None, help="a priori max vel; default from MOD.")
    p.add_argument("--max-iter", type=int, default=8)
    p.add_argument("--sparsity", type=float, default=0.2)
    p.add_argument("--checkerboard", action="store_true", help="also write a checkerboard MOD.true.")
    p.add_argument("--anosize", type=float, default=1.0, help="checkerboard angular size.")
    p.add_argument("--amplitude", type=float, default=0.4, help="checkerboard amplitude (km/s).")
    p.add_argument("--synthetic", type=int, default=None, choices=[0, 1],
                   help="DSurfTomo synthetic flag; default 1 if --checkerboard else 0.")
    p.add_argument("--noise", type=float, default=0.02)
    p.add_argument("--threshold", type=float, default=3.0)
    return p.parse_args(argv)


def main():
    print_banner(
        "PrepareDSurfTomoProject.py",
        [
            "Build a DSurfTomo project from dispersion data in one interactive run.",
            "Equivalent to running extractnew.py, GenerateIniMOD.py and GenerateDSurfTomoInputFile.py in sequence.",
        ],
    )
    print_section("Step 1: Build surfdata")

    input_dir = prompt("input dispersion directory")
    input_type = ex.infer_input_type(input_dir)
    curve_type = infer_curve_type(input_dir, input_type)
    print_kv("Detected input source", input_type)
    print_kv("Detected curve type", curve_type)

    # EGFAnalysisPy curves carry station coordinates in their headers, so no
    # separate station file is needed; station points are derived from surfdata.
    station_file = None

    project_dir = prompt("output project directory", "DSurfTomo_project")
    os.makedirs(project_dir, exist_ok=True)

    period_step = prompt_float("strict period step (s, exact grid only)")
    if period_step <= 0:
        raise ValueError("period step must be positive.")
    wave_type = prompt_int("wave type (1=Love, 2=Rayleigh)", "2")

    surfdata_path = os.path.join(project_dir, "surfdata.dat")
    args = ex.Args()
    args.input_dir = input_dir
    args.input_type = input_type
    args.curve_type = curve_type
    args.station_file = station_file
    args.min_period = None
    args.max_period = None
    args.period_step = period_step
    args.wave_type = wave_type
    args.vel_type = None
    args.output = surfdata_path

    runs = build_runs(args)
    recommended_nsrc, metadata_path = write_surfdata(args, runs)

    print_section("Step 1 Complete")
    print_kv("surfdata", surfdata_path)
    print_kv("surfdata metadata", metadata_path)
    for run in runs:
        print_kv(
            "%s periods auto-detected" % run["curve_type"],
            "%.3f-%.3f s (%d used)"
            % (
                run["min_period"],
                run["max_period"],
                run["nf"],
            ),
        )
    print_kv("recommended line 8 max(sources, receivers)", recommended_nsrc)

    records = gim.load_surfdata(surfdata_path)
    metadata = gim.load_surfdata_metadata(surfdata_path)
    if metadata is not None:
        period_map = gim.build_period_map_from_metadata(records, metadata)
        default_wave_type = str(metadata.get("wave_type", 2))
        curve_types = metadata.get("curve_types", [])
        default_vel_selection = "phase" if "phase" in curve_types else "group"
    else:
        mapping_mode = prompt_choice(
            "period mapping mode (range/manual)", ["range", "manual"], default="range"
        )
        period_map = gim.build_period_map(mapping_mode, records)
        default_wave_type = "2"
        default_vel_selection = "phase"

    print_section("Step 2: Build MOD")
    points = choose_station_points(station_file, surfdata_path)
    nx_default, ny_default, grid_info = infer_default_nx_ny(points)
    dvxd_default, dvzd_default, goxd, gozd = compute_region_defaults(grid_info, nx_default, ny_default)

    print_section("Grid Parameter Mode")
    grid_mode = prompt_choice(
        "grid parameter mode (count/spacing)", ["count", "spacing"], default="count"
    )
    if grid_mode == "count":
        print_kv("Mode", "user specifies nx and ny; dvxd and dvzd are computed")
        nx = prompt_int("nx", str(nx_default))
        ny = prompt_int("ny", str(ny_default))
        dvxd, dvzd, goxd, gozd = compute_region_defaults(grid_info, nx, ny)
    else:
        print_kv("Mode", "user specifies dvxd and dvzd; nx and ny are computed")
        dvxd = prompt_float("dvxd", "%.5f" % dvxd_default)
        dvzd = prompt_float("dvzd", "%.5f" % dvzd_default)
        nx, ny, goxd, gozd = compute_grid_counts_from_spacing(grid_info, dvxd, dvzd)

    print_section("Grid Parameters")
    print_kv("grid mode", grid_mode)
    print_kv("nx ny", "%d %d" % (nx, ny))
    print_kv("dvxd dvzd", "%.5f %.5f" % (dvxd, dvzd))
    print_kv("goxd gozd", "%.5f %.5f" % (goxd, gozd))
    dep1 = np.array(
        prompt_depths_with_arithmetic(
            "depth array dep1 in km", "[0,0.2,0.4,0.6,0.8,1.1,1.4,1.8,2.5]"
        ),
        dtype=float,
    )

    minvel = 2.0
    velgrad = 0.3
    mod_mode = prompt_choice(
        "initial model parameter mode (estimate/manual)",
        ["estimate", "manual"],
        default="estimate",
    )
    if mod_mode == "estimate":
        print("Wave type in surfdata: Love=1, Rayleigh=2")
        model_wave_type = prompt_int(
            "wave type used to estimate initial model", default=default_wave_type
        )
        print("Velocity type in surfdata: phase=0, group=1")
        selection = prompt_choice(
            "velocity type used to estimate initial model (phase/group)",
            ["phase", "group"],
            default=default_vel_selection,
        )
        model_vel_type = 0 if selection == "phase" else 1

        period1 = prompt_float("short period T1 (s)", "1.0")
        period2 = prompt_float("long period T2 (s)", "3.0")
        if period2 <= period1:
            raise ValueError("T2 must be greater than T1.")

        vs_1, n1, used_period1 = gim.mean_velocity_at_period(
            records, period_map, period1, model_wave_type, model_vel_type, prefer="lower"
        )
        vs_2, n2, used_period2 = gim.mean_velocity_at_period(
            records, period_map, period2, model_wave_type, model_vel_type, prefer="upper"
        )
        minvel, velgrad, fellback, d1, d2 = gim.estimate_minvel_velgrad(
            vs_1, used_period1, vs_2, used_period2
        )

        print_section("Step 2 Summary")
        print_kv(
            "Average velocity near T1=%.3f s" % period1,
            "uses T=%.3f s: %.4f km/s from %d records" % (used_period1, vs_1, n1),
        )
        print_kv(
            "Average velocity near T2=%.3f s" % period2,
            "uses T=%.3f s: %.4f km/s from %d records" % (used_period2, vs_2, n2),
        )
        if n1 < 10:
            print("  ! Warning: T1 estimate uses only %d records, which is not very stable." % n1)
        if n2 < 10:
            print("  ! Warning: T2 estimate uses only %d records, which is not very stable." % n2)
        print_kv("Depth bounds from empirical rule", "d1=%.4f km, d2=%.4f km" % (d1, d2))
        print_kv("Recommended model parameters", "minvel=%.4f km/s, velgrad=%.4f km/s/km" % (minvel, velgrad))
        if fellback:
            print("  ! Warning: flat/decreasing dispersion (d1=%.4f, d2=%.4f); using a "
                  "fallback gentle-gradient model. Consider manual minvel/velgrad." % (d1, d2))

        override = prompt_choice(
            "use recommended minvel/velgrad? (yes/no)", ["yes", "no"], default="yes"
        )
        if override == "no":
            minvel = prompt_float("minvel (km/s)", "%.4f" % minvel)
            velgrad = prompt_float("velgrad (km/s/km)", "%.4f" % velgrad)
    else:
        minvel = prompt_float("minvel (km/s)", "%.4f" % minvel)
        velgrad = prompt_float("velgrad (km/s/km)", "%.4f" % velgrad)

    mod_path = os.path.join(project_dir, "MOD")
    gim.generate_mod(nx, ny, dep1, minvel, velgrad, mod_path)

    print_section("Step 2 Complete")
    print_kv("MOD", mod_path)

    print_section("Step 2b: Optional MOD.true Checkerboard")
    generate_true_mod = prompt_choice(
        "generate checkerboard MOD.true? (yes/no)",
        ["yes", "no"],
        default="no",
    )
    true_mod_path = None
    if generate_true_mod == "yes":
        checkerboard_mode = prompt_choice(
            "checkerboard parameter mode (default/manual)",
            ["default", "manual"],
            default="default",
        )
        anosize = 0.5
        amplitude = 0.4
        if checkerboard_mode == "manual":
            anosize = prompt_float("checkerboard angular size", "%.3f" % anosize)
            amplitude = prompt_float("checkerboard amplitude (km/s)", "%.3f" % amplitude)
        else:
            print_kv("checkerboard angular size", "%.3f" % anosize)
            print_kv("checkerboard amplitude", "%.3f km/s" % amplitude)
        true_mod_path = os.path.join(project_dir, "MOD.true")
        _, background_mod = gtm.load_mod(mod_path)
        true_mod = gtm.generate_true_mod(background_mod, anosize=anosize, amplitude=amplitude)
        gtm.write_true_mod(true_mod_path, true_mod)

        print_section("Step 2b Complete")
        print_kv("MOD.true", true_mod_path)
        print_kv(
            "checkerboard velocity range",
            "%.3f-%.3f km/s" % (float(np.min(true_mod)), float(np.max(true_mod))),
        )
    else:
        print_kv("MOD.true", "not generated")

    surf_records, period_lists = build_period_lists(surfdata_path)
    data_file_name = os.path.basename(surfdata_path)
    print_section("Step 3: Build DSurfTomo.in")
    if station_file and os.path.exists(station_file):
        print_kv("Using station coordinate file", station_file)
    else:
        print_kv("Coordinate source", "surfdata.dat")
    print_kv("Using data file name in DSurfTomo.in", data_file_name)
    print_kv("Region bounds", "tight envelope from station coverage with small safety margin")

    print_section("Grid Parameters")
    print_kv("grid mode", grid_mode)
    print_kv("nx ny", "%d %d" % (nx, ny))
    print_kv("dvxd dvzd", "%.5f %.5f" % (dvxd, dvzd))
    print_kv("goxd gozd", "%.5f %.5f" % (goxd, gozd))

    print_section("Line 8 Parameter")
    print_kv("max(sources, receivers)", recommended_nsrc)
    use_line8 = prompt_choice(
        "use computed line-8 value? (yes/no)", ["yes", "no"], "yes"
    )
    if use_line8 == "yes":
        nsrcmax = recommended_nsrc
    else:
        nsrcmax = prompt_int("max(sources, receivers)", str(recommended_nsrc))

    mod_bounds = infer_velocity_bounds_from_mod(mod_path)
    if mod_bounds is not None:
        default_minvel, default_maxvel = mod_bounds
    else:
        from dsurftomo_utils import infer_velocity_bounds
        default_minvel, default_maxvel = infer_velocity_bounds(
            [rec["velocity"] for rec in surf_records]
        )

    print_section("Inversion Parameters")
    weight = prompt_float("weight", "4.0")
    damp = prompt_float("damp", "1.0")
    sublayers = prompt_int("sublayers", "3")

    print_section("Velocity Bounds")
    print_kv(
        "Computed velocity bounds",
        "minimum velocity=%.3f maximum velocity=%.3f" % (default_minvel, default_maxvel),
    )
    use_velocity_bounds = prompt_choice(
        "use computed velocity bounds? (yes/no)", ["yes", "no"], "yes"
    )
    if use_velocity_bounds == "yes":
        minvel_input, maxvel_input = default_minvel, default_maxvel
    else:
        minvel_input = prompt_float("minimum velocity", "%.3f" % default_minvel)
        maxvel_input = prompt_float("maximum velocity", "%.3f" % default_maxvel)

    print_section("Stopping / Noise Parameters")
    maxiteration = prompt_int("maximum iteration", "10")
    sparsityfraction = prompt_float("sparsity fraction", "0.2")
    check_sparse_matrix_size(len(records), nx, ny, len(dep1), sparsityfraction)
    synthetic_flag_default = "1" if true_mod_path is not None else "0"
    synthetic_flag = prompt_int(
        "synthetic flag (0=real, 1=synthetic)", synthetic_flag_default
    )
    noiselevel = prompt_float("noise level", "0.02")
    threshold = prompt_float("threshold", "0.5")

    input_path = os.path.join(project_dir, "DSurfTomo.in")
    gdi.write_input_file(
        output_path=input_path,
        data_file_name=data_file_name,
        nx=nx,
        ny=ny,
        nz=len(dep1),
        goxd=goxd,
        gozd=gozd,
        dvxd=dvxd,
        dvzd=dvzd,
        nsrcmax=nsrcmax,
        weight=weight,
        damp=damp,
        sublayers=sublayers,
        minvel=minvel_input,
        maxvel=maxvel_input,
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

    print_section("Step 3 Complete")
    print_kv("DSurfTomo.in", input_path)
    print_section("Project Ready")
    print_kv("Project directory", project_dir)
    print("  - surfdata.dat")
    print("  - surfdata.json")
    print("  - MOD")
    if true_mod_path is not None:
        print("  - MOD.true")
    print("  - DSurfTomo.in")


if __name__ == "__main__":
    # Any command-line argument switches to the non-interactive build_project path.
    if len(sys.argv) > 1:
        build_project(parse_cli(sys.argv[1:]))
    else:
        main()
