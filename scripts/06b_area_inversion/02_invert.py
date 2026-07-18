#!/usr/bin/env python3
"""
Step 06b-2 — run the compiled DSurfTomo executable on a prepared project.

The executable reads DSurfTomo.in, surfdata.dat and MOD from the project
directory (built by 01_prepare.py), so it is executed with that directory as
the working directory. Console output is also saved to dsurftomo_run.log.

    python scripts/06b_area_inversion/02_invert.py \
        --config configs/06b_area_inversion/dsurftomo.yaml
"""

import argparse
import os
import subprocess
import sys

import _common

REQUIRED_INPUTS = ["DSurfTomo.in", "surfdata.dat", "MOD"]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    _common.add_common_args(parser)
    args = parser.parse_args()

    cfg = _common.load_yaml(args.config)
    project = _common.load_yaml(args.project_config, required=False)
    project_dir = _common.get_project_dir(cfg)
    binary = _common.get_dsurftomo_bin(cfg, project)

    missing = [f for f in REQUIRED_INPUTS
               if not os.path.isfile(os.path.join(project_dir, f))]
    if missing:
        sys.exit(f"ERROR: missing input file(s) in {project_dir}: "
                 f"{', '.join(missing)}\n       Run 01_prepare.py first.")

    log_path = os.path.join(project_dir, "dsurftomo_run.log")
    print(f"running {binary}\n  cwd: {project_dir}\n  log: {log_path}")
    with open(log_path, "w", encoding="utf-8") as log:
        proc = subprocess.Popen([binary], cwd=project_dir,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True)
        # The executable asks for the input file name on stdin
        # (interactive prompt in main.f90); answer with the default.
        proc.stdin.write("DSurfTomo.in\n")
        proc.stdin.close()
        for line in proc.stdout:
            print(line, end="")
            log.write(line)
        ret = proc.wait()
    if ret != 0:
        sys.exit(f"ERROR: DSurfTomo exited with status {ret} (see {log_path})")
    print("done. inversion outputs are in the project directory.")


if __name__ == "__main__":
    main()
