#!/usr/bin/env python3
"""
Step 06b-3 — plot DSurfTomo results with the bundled GMT scripts.

Unified plotting entry: every script listed under `plot.scripts` in the
YAML config (from DSurfTomo/scripts/, e.g. plotdepth_subplot.gmt,
plotcross_subplot.gmt, plotray_subplot.gmt, plotchecker_subplot.gmt) is
run with --project-dir/--output-dir, plus any per-script extra arguments.
Requires GMT 6 (docs/INSTALL.md §5). Figures go to <project_dir>/figs by
default.

Config forms accepted under plot.scripts:

    plot:
      scripts:
        - plotdepth_subplot.gmt                # plain name
        - name: plotcross_subplot.gmt          # name + extra arguments
          args: ["--profiles", "lon1/lat1/lon2/lat2/label"]

plotdepth_subplot.gmt automatically receives --depths from prepare.depths.

    python scripts/06b_area_inversion/03_plot.py \
        --config configs/06b_area_inversion/dsurftomo.yaml
"""

import argparse
import os
import subprocess
import sys

import _common


def normalize_entries(cfg):
    entries = (cfg.get("plot") or {}).get("scripts") or ["plotdepth_subplot.gmt"]
    normalized = []
    for entry in entries:
        if isinstance(entry, str):
            normalized.append((entry, []))
        elif isinstance(entry, dict) and entry.get("name"):
            args = entry.get("args") or []
            normalized.append((entry["name"], [str(a) for a in args]))
        else:
            sys.exit(f"ERROR: invalid plot.scripts entry: {entry!r} "
                     f"(use a script name or {{name: ..., args: [...]}})")
    return normalized


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    _common.add_common_args(parser)
    parser.add_argument("--output-dir", default=None,
                        help="figure output directory (default: <project_dir>/figs)")
    args = parser.parse_args()

    cfg = _common.load_yaml(args.config)
    project_dir = _common.get_project_dir(cfg)
    if not os.path.isfile(os.path.join(project_dir, "DSurfTomo.in")):
        sys.exit(f"ERROR: {project_dir} does not look like a DSurfTomo project "
                 f"(no DSurfTomo.in). Run 01_prepare.py / 02_invert.py first.")

    output_dir = os.path.abspath(args.output_dir or os.path.join(project_dir, "figs"))
    os.makedirs(output_dir, exist_ok=True)
    depths = (cfg.get("prepare") or {}).get("depths")

    for name, extra in normalize_entries(cfg):
        script = os.path.join(_common.REPO_ROOT, "DSurfTomo", "scripts", name)
        if not os.path.isfile(script):
            sys.exit(f"ERROR: plot script not found: {script}")
        cmd = ["bash", script, "--project-dir", project_dir,
               "--output-dir", output_dir]
        if name.startswith("plotdepth") and "--depths" not in extra:
            if not depths:
                sys.exit("ERROR: prepare.depths missing in the stage config "
                         "(plotdepth_subplot.gmt needs the depth list)")
            cmd += ["--depths", str(depths)]
        cmd += extra
        print("command:", " ".join(cmd))
        ret = subprocess.call(cmd)
        if ret != 0:
            sys.exit(f"ERROR: {name} exited with status {ret}")
    print(f"done. figures in {output_dir}")


if __name__ == "__main__":
    main()
