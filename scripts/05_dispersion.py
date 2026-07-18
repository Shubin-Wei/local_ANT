#!/usr/bin/env python3
"""
Step 05 — dispersion analysis & curve picking (EGFAnalysisPy + DisperPicker).

Reads configs/05_dispersion.yaml, renders it into the runtime config module
EGFAnalysisPy expects (written into the result directory, so the bundled
EGFAnalysisPy/config.py is never modified), then invokes
EGFAnalysisPy/run_parallel.py with that module. Serial runs simply use one
worker.

Run from the repo root, inside the *egfanalysispy* conda environment
(Python 3.6 + TensorFlow 1.13 — see docs/INSTALL.md §3):

    python scripts/05_dispersion.py --config configs/05_dispersion.yaml
"""

import argparse
import os
import subprocess
import sys
from datetime import datetime

import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# YAML section/key -> Config attribute (all optional; bundled config.py
# supplies every default). Paths are handled separately.
_PARAM_MAP = {
    "signal": ["isEGF", "WinAlpha", "NoiseTime", "MinSNR", "WinPeriodNum",
               "WinMinTime", "FilterKaiserPara", "MaxFilterLengthLog"],
    "dispersion": ["StartV", "EndV", "StartT", "EndT", "DeltaT", "DeltaV",
                   "MinStaDist", "MinWavelengthRatio"],
    "picker": ["confidence_G", "mean_confidence_C", "confidence_C",
               "random_plot_ratio", "ref_T", "ref_T2"],
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--config", required=True,
                   help="stage config, e.g. configs/05_dispersion.yaml")
    p.add_argument("--project-config",
                   default=os.path.join(REPO_ROOT, "configs", "00_project.yaml"),
                   help="project config supplying egfanalysispy_dir "
                        "(optional; default: configs/00_project.yaml)")
    p.add_argument("--workers", type=int, default=None,
                   help="override run.workers from the config")
    p.add_argument("--overwrite", action="store_true",
                   help="allow writing into a non-empty result directory")
    p.add_argument("--dry-run", action="store_true",
                   help="generate the runtime config and print the command, "
                        "but do not run EGFAnalysisPy")
    return p.parse_args()


def load_yaml(path, required=True):
    if not os.path.isfile(path):
        if required:
            sys.exit(f"ERROR: config file not found: {path}")
        return {}
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def render_runtime_config(cfg, egf_dir, result_path, source_path):
    """Write the runtime config module and return its path."""
    overrides = {}
    for section, keys in _PARAM_MAP.items():
        for key in keys:
            value = (cfg.get(section) or {}).get(key)
            if value is not None:
                overrides[key] = value
    green = (cfg.get("signal") or {}).get("GreenFcnObjects")

    lines = [
        "# AUTO-GENERATED runtime config — do not edit.",
        f"# source : {os.path.abspath(source_path)}",
        f"# created: {datetime.now().isoformat(timespec='seconds')}",
        "# Loaded as module 'config' by EGFAnalysisPy/run_parallel.py; the",
        "# bundled EGFAnalysisPy/config.py stays untouched.",
        "import importlib.util as _ilu",
        "import os as _os",
        f"_EGF_DIR = {os.path.abspath(egf_dir)!r}",
        "_spec = _ilu.spec_from_file_location('_egf_base_config',",
        "                                     _os.path.join(_EGF_DIR, 'config.py'))",
        "_base = _ilu.module_from_spec(_spec)",
        "_spec.loader.exec_module(_base)",
        "GreenFcnObjectsType = _base.GreenFcnObjectsType",
        "",
        "",
        "class Config(_base.Config):",
        "    def __init__(self):",
        "        super().__init__()",
        f"        self.waveform_path = {os.path.abspath(cfg['paths']['waveform_path'])!r}",
        f"        self.result_path = {os.path.abspath(result_path)!r}",
    ]
    for key, value in overrides.items():
        lines.append(f"        self.{key} = {value!r}")
    if green is not None:
        lines.append(f"        self.GreenFcnObjects = GreenFcnObjectsType.{green}")
    lines += [
        "        # Recompute quantities derived in the base class (same formulas",
        "        # as EGFAnalysisPy/config.py) now that overrides are applied.",
        "        self.dT = self.DeltaT",
        "        self.dV = self.DeltaV",
        "        self.range_T = [self.StartT, self.EndT,",
        "                        round((self.EndT - self.StartT)/self.dT) + 1]",
        "        self.range_V = [self.StartV, self.EndV,",
        "                        round((self.EndV - self.StartV)/self.dV) + 1]",
        "        self.input_size = [self.range_V[2], self.range_T[2], 2]",
        "        self.min_len = round(self.range_T[2]/8)",
        "",
    ]
    out = os.path.join(result_path, "runtime_config.py")
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    return out


def main():
    args = parse_args()
    cfg = load_yaml(args.config)
    project = load_yaml(args.project_config, required=False)

    paths = cfg.get("paths") or {}
    for key in ("waveform_path", "result_path"):
        if not paths.get(key):
            sys.exit(f"ERROR: paths.{key} missing in {args.config}")

    egf_dir = project.get("egfanalysispy_dir", "EGFAnalysisPy")
    if not os.path.isabs(egf_dir):
        egf_dir = os.path.join(REPO_ROOT, egf_dir)
    runner = os.path.join(egf_dir, "run_parallel.py")
    if not os.path.isfile(runner):
        sys.exit(f"ERROR: run_parallel.py not found under {egf_dir}")

    waveform_path = os.path.abspath(paths["waveform_path"])
    if not os.path.isdir(waveform_path) or not os.listdir(waveform_path):
        sys.exit(f"ERROR: waveform_path missing or empty: {waveform_path}\n"
                 f"       (expected output of step 04)")

    result_path = os.path.abspath(paths["result_path"])
    qc_dir = os.path.join(result_path, "qc_result")
    if os.path.isdir(qc_dir) and os.listdir(qc_dir) and not args.overwrite:
        sys.exit(f"ERROR: {qc_dir} already has results; "
                 f"use --overwrite to continue anyway.")
    os.makedirs(result_path, exist_ok=True)

    runtime_cfg = render_runtime_config(cfg, egf_dir, result_path, args.config)
    print(f"runtime config written: {runtime_cfg}")

    run_cfg = cfg.get("run") or {}
    workers = args.workers or (run_cfg.get("workers") if run_cfg.get("parallel") else 1) or 1
    cmd = [sys.executable, "run_parallel.py",
           "--config", runtime_cfg, "-w", str(workers)]
    print(f"command: {' '.join(cmd)}  (cwd: {egf_dir})")
    if args.dry_run:
        print("dry-run: stopping before EGFAnalysisPy is invoked.")
        return

    try:
        import tensorflow  # noqa: F401  (fail early with a clear message)
    except ImportError:
        sys.exit("ERROR: TensorFlow is not importable. Step 05 must run inside "
                 "the 'egfanalysispy' conda environment (docs/INSTALL.md §3).")

    ret = subprocess.call(cmd, cwd=egf_dir)
    if ret != 0:
        sys.exit(f"ERROR: EGFAnalysisPy exited with status {ret}")
    print(f"done. results in {result_path} (picked curves in qc_result/)")


if __name__ == "__main__":
    main()
