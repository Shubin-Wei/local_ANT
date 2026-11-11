#!/usr/bin/env python
"""Run a known-truth synthetic example without hiding unexpected failures."""

import argparse
import json
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dir", required=True, type=Path)
    parser.add_argument("--full", action="store_true", help="include solver, evaluation and figures")
    args = parser.parse_args()
    if not all((args.dir / name).is_file() for name in ("true_model.txt", "true_dispersion.txt", "synthetic.yaml")):
        parser.error("this runner is only for generated synthetic examples with known truth")
    config = str((args.dir / "synthetic.yaml").resolve())
    base = [sys.executable, "-m", "lant"]
    target = [] if args.full else ["--to", "inversion"]
    for _ in range(6):
        result = subprocess.run(base + ["run", "all", "-c", config] + target)
        if result.returncode == 0:
            return
        if result.returncode != 3:
            raise SystemExit(result.returncode)
        state = json.loads(subprocess.check_output(base + ["status", "-c", config, "--json"]))
        gate = next((step for step, value in state["steps"].items() if value["state"] == "to confirm"), None)
        if gate is None:
            raise SystemExit("stopped without a scientific gate")
        subprocess.run(
            base + ["confirm", gate, "-c", config, "--reason", "Generated synthetic input with known truth"], check=True
        )
    raise SystemExit("too many synthetic confirmation stops")


if __name__ == "__main__":
    main()
