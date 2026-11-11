#!/usr/bin/env python
"""Synthetic test data with a known answer, written next to this script by default.

The generator is lant.synthetic (`python -m lant.synthetic --help`); `lant demo` uses it to run a whole example.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from lant.synthetic import main  # noqa: E402

if __name__ == "__main__":
    main(default_out=HERE)
