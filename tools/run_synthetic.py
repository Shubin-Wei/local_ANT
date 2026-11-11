#!/usr/bin/env python
"""Run a known-truth synthetic example without hiding unexpected failures (lant.demo.run_known_truth)."""

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from lant.demo import run_known_truth  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dir", required=True, type=Path)
    parser.add_argument("--full", action="store_true", help="include solver, evaluation and figures")
    args = parser.parse_args()
    try:
        run_known_truth(args.dir, args.full)
    except RuntimeError as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
