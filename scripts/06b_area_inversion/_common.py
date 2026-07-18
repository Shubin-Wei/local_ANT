"""Shared helpers for the 06b areal-inversion entry scripts."""

import os
import sys

import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_PROJECT_CONFIG = os.path.join(REPO_ROOT, "configs", "00_project.yaml")


def load_yaml(path, required=True):
    if not os.path.isfile(path):
        if required:
            sys.exit(f"ERROR: config file not found: {path}")
        return {}
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def add_common_args(parser):
    parser.add_argument("--config", required=True,
                        help="stage config, e.g. configs/06b_area_inversion/dsurftomo.yaml")
    parser.add_argument("--project-config", default=DEFAULT_PROJECT_CONFIG,
                        help="project config (optional; default: configs/00_project.yaml)")


def resolve_in_repo(path):
    return path if os.path.isabs(path) else os.path.join(REPO_ROOT, path)


def get_project_dir(cfg):
    project_dir = (cfg.get("paths") or {}).get("project_dir")
    if not project_dir:
        sys.exit("ERROR: paths.project_dir missing in the stage config")
    return os.path.abspath(project_dir)


def get_dsurftomo_bin(cfg, project):
    """Stage config wins over 00_project.yaml (config precedence rule)."""
    binary = (cfg.get("paths") or {}).get("dsurftomo_bin") \
        or project.get("dsurftomo_bin", "DSurfTomo/bin/DSurfTomo")
    binary = resolve_in_repo(binary)
    if not (os.path.isfile(binary) and os.access(binary, os.X_OK)):
        sys.exit(f"ERROR: DSurfTomo executable not found or not executable: {binary}\n"
                 f"       Compile it first (docs/INSTALL.md §4), then set "
                 f"dsurftomo_bin in configs/00_project.yaml.")
    return binary
