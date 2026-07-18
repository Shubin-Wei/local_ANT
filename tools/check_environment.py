#!/usr/bin/env python3
"""
Environment checker — reports what is usable, installs or changes nothing.

Run it inside whichever conda environment you want to inspect:

    conda activate noisepy       && python tools/check_environment.py
    conda activate egfanalysispy && python tools/check_environment.py

It prints PASS/FAIL per item; items irrelevant to the current environment
simply show as FAIL there (e.g. TensorFlow inside the noisepy env is fine
to fail). See docs/INSTALL.md for how to fix failures.
"""

import importlib
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CHECKS = [
    # (label, modules) — grouped by which environment normally provides them
    ("noisepy env / steps 01-04", [
        "noisepy", "obspy", "numpy", "scipy", "pandas", "matplotlib",
        "h5py", "yaml", "datetimerange",
    ]),
    ("linear inversion / step 06a", [
        "evodcinv", "disba", "tqdm",
    ]),
    ("egfanalysispy env / step 05", [
        "tensorflow", "geopy", "seaborn",
    ]),
]

FILES = [
    ("DisperPicker pretrained model",
     os.path.join(REPO_ROOT, "EGFAnalysisPy", "DisperPicker", "saver", "checkpoint")),
    ("DSurfTomo executable (step 06b)",
     os.path.join(REPO_ROOT, "DSurfTomo", "bin", "DSurfTomo")),
]


def try_import(name):
    try:
        mod = importlib.import_module(name)
        version = getattr(mod, "__version__", "")
        return True, version
    except Exception as exc:  # ImportError or any init failure
        return False, str(exc).splitlines()[0]


def main():
    print("Python {}  ({})\n".format(sys.version.split()[0], sys.executable))
    failures = 0
    for group, modules in CHECKS:
        print("[{}]".format(group))
        for name in modules:
            ok, info = try_import(name)
            mark = "PASS" if ok else "FAIL"
            if not ok:
                failures += 1
            print("  {:4s}  {:15s} {}".format(mark, name, info))
        print()

    print("[file checks]")
    for label, path in FILES:
        if os.path.exists(path):
            executable = os.access(path, os.X_OK)
            note = "" if (not path.endswith("DSurfTomo") or executable) \
                else "  (exists but NOT executable)"
            print("  PASS  {}{}".format(label, note))
        else:
            failures += 1
            print("  FAIL  {}: {} missing".format(label, path))

    print("\n{} item(s) failed in this environment. "
          "See docs/INSTALL.md".format(failures))


if __name__ == "__main__":
    main()
