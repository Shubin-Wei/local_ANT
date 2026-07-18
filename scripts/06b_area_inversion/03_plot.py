#!/usr/bin/env python3
"""
Step 06b-3 — plot DSurfTomo results with the bundled GMT scripts.

Thin wrapper around the scripts listed in `plot.scripts` (from
DSurfTomo/scripts/, e.g. plotdepth_subplot.gmt). Requires GMT 6
(docs/INSTALL.md §5). Figures go to <project_dir>/figs by default.

    python scripts/06b_area_inversion/03_plot.py \
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
    parser.add_argument("--output-dir", default=None,
                        help="figure output directory (default: <project_dir>/figs)")
    args = parser.parse_args()

    cfg = _common.load_yaml(args.config)
    project_dir = _common.get_project_dir(cfg)
    if not os.path.isfile(os.path.join(project_dir, "DSurfTomo.in")):
        sys.exit(f"ERROR: {project_dir} does not look like a DSurfTomo project "
                 f"(no DSurfTomo.in). Run 01_prepare.py / 02_invert.py first.")

    depths = (cfg.get("prepare") or {}).get("depths")
    if not depths:
        sys.exit("ERROR: prepare.depths missing in the stage config "
                 "(plotdepth_subplot.gmt needs the depth list)")
    output_dir = os.path.abspath(args.output_dir or os.path.join(project_dir, "figs"))
    os.makedirs(output_dir, exist_ok=True)

    scripts = (cfg.get("plot") or {}).get("scripts") or ["plotdepth_subplot.gmt"]
    for name in scripts:
        script = os.path.join(_common.REPO_ROOT, "DSurfTomo", "scripts", name)
        if not os.path.isfile(script):
            sys.exit(f"ERROR: plot script not found: {script}")
        cmd = ["bash", script, "--project-dir", project_dir,
               "--depths", str(depths), "--output-dir", output_dir]
        print("command:", " ".join(cmd))
        ret = subprocess.call(cmd)
        if ret != 0:
            sys.exit(f"ERROR: {name} exited with status {ret}")
    print(f"done. figures in {output_dir}")


if __name__ == "__main__":
    main()
