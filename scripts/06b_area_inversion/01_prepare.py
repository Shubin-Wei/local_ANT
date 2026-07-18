#!/usr/bin/env python3
"""
Step 06b-1 — build a DSurfTomo project from step-05 dispersion results.

Thin wrapper: translates the `prepare:` section of the YAML config into the
command line of DSurfTomo/scripts/PrepareDSurfTomoProject.py (the actual
implementation), which writes surfdata.dat, MOD and DSurfTomo.in into the
project directory.

    python scripts/06b_area_inversion/01_prepare.py \
        --config configs/06b_area_inversion/dsurftomo.yaml
"""

import argparse
import os
import subprocess
import sys

import _common


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    _common.add_common_args(parser)
    parser.add_argument("--overwrite", action="store_true",
                        help="allow writing into a non-empty project directory")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the underlying command without running it")
    args = parser.parse_args()

    cfg = _common.load_yaml(args.config)
    input_dir = (cfg.get("paths") or {}).get("input_dir")
    if not input_dir:
        sys.exit("ERROR: paths.input_dir missing in the stage config")
    input_dir = os.path.abspath(input_dir)
    if not os.path.isdir(input_dir):
        sys.exit(f"ERROR: input_dir does not exist: {input_dir}\n"
                 f"       (expected the qc_result directory from step 05)")

    project_dir = _common.get_project_dir(cfg)
    if os.path.isdir(project_dir) and os.listdir(project_dir) and not args.overwrite:
        sys.exit(f"ERROR: {project_dir} is not empty; use --overwrite to continue.")

    prep = cfg.get("prepare") or {}
    spacing, count = prep.get("grid_spacing"), prep.get("grid_count")
    if bool(spacing) == bool(count):
        sys.exit("ERROR: set exactly one of prepare.grid_spacing / prepare.grid_count")
    if prep.get("period_step") is None:
        sys.exit("ERROR: prepare.period_step is required")

    tool = os.path.join(_common.REPO_ROOT, "DSurfTomo", "scripts",
                        "PrepareDSurfTomoProject.py")
    cmd = [sys.executable, tool,
           "--input-dir", input_dir, "--project-dir", project_dir,
           "--period-step", str(prep["period_step"])]
    if spacing:
        cmd += ["--grid-spacing", str(spacing[0]), str(spacing[1])]
    else:
        cmd += ["--grid-count", str(count[0]), str(count[1])]
    for key, flag in [("curve_type", "--curve-type"), ("min_period", "--min-period"),
                      ("max_period", "--max-period"), ("wave_type", "--wave-type"),
                      ("depths", "--depths"), ("minvel", "--minvel"),
                      ("velgrad", "--velgrad"), ("t1", "--t1"), ("t2", "--t2"),
                      ("model_vel", "--model-vel"), ("nsrcmax", "--nsrcmax"),
                      ("weight", "--weight"), ("damp", "--damp"),
                      ("sublayers", "--sublayers")]:
        if prep.get(key) is not None:
            cmd += [flag, str(prep[key])]

    print("command:", " ".join(cmd))
    if args.dry_run:
        return
    os.makedirs(project_dir, exist_ok=True)
    ret = subprocess.call(cmd)
    if ret != 0:
        sys.exit(f"ERROR: PrepareDSurfTomoProject.py exited with status {ret}")
    print(f"done. DSurfTomo project prepared in {project_dir}")


if __name__ == "__main__":
    main()
