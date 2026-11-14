"""The rules file for AI assistants in a project: one file, for the assistant chosen in project.ai_assistant, with
the full text from one template; AGENTS.md (the open convention read by Codex and many other coding agents) when
none is chosen.  An existing file is never overwritten (it may hold the user's own rules)."""

import os
from pathlib import Path

from .delivery import atomic_text

RULES_FILES = {None: "AGENTS.md", "codex": "AGENTS.md", "claude": "CLAUDE.md", "gemini": "GEMINI.md"}


def rules_file(cfg):
    return RULES_FILES[cfg["project"]["ai_assistant"]]


def write(cfg, root):
    """Write the rules file into the project directory `root` unless it exists -> its path."""
    root = Path(root)
    config = cfg.get("_config_file") or "CONFIG"
    if config != "CONFIG":
        rel = os.path.relpath(os.path.abspath(config), os.path.abspath(str(root)))
        config = rel if not rel.startswith("..") else os.path.abspath(config)
    rules = root / rules_file(cfg)
    if not rules.exists():
        atomic_text(rules, PROJECT_AGENTS.replace("CONFIG", config))
    return rules


PROJECT_AGENTS = """# Rules for AI assistants in this lant project

This directory is written by lant.  Start with REPORT.md and `lant status -c CONFIG --json`.

## Files

- Create files only in two ways: notes with `lant note "text" -c CONFIG`, and temporary scripts, extracts or
  figures under `scratch/` (create it when needed; delete what you no longer need).  Nowhere else in this
  directory, and not in the raw-data directory or the lant repository.
- Do not edit or delete anything under `work/`, `results/` or `.lant/`, nor REPORT.md, README.md or this file;
  lant writes them.  Change the configuration only when the user asks, and record the reason with `lant note`.
- `lant status -c CONFIG --json` lists under `artifacts.unregistered` every file lant did not write; it should hold
  only what the user put here.

## Decisions

- When `next.requires_user` is true, show the user the items and figures of `.lant/steps/<step>.md` and wait.  Record
  the user's decision with `lant confirm STEP -c CONFIG --reason "the user's decision and evidence"`.
- Do not loosen quality control or change scientific priors to get more points or nicer figures.
- Report numbers and figures, and say what was not checked.

## Running

- `lant run all -c CONFIG` continues from the first step that does not stand; after an interruption run it again.
- While `writer` is set in the status, a step is running: wait, do not start another lant command that writes.
- The number of cores is the program's (resources.cores: auto, 80 % of the idle cores, taken at every step).  Do not
  set or change it unless the user asks; report the cores and memory each step used (its COMPUTING RESOURCES).
"""
