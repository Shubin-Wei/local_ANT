"""Command line: lant init | clean | cc | disp | inv | checker | export-sac | check.

Run the steps one at a time and look at each step's log and figures before the next one (docs/tutorial.md).
"""

import argparse
import sys

from . import __version__

STEPS = """steps (run them in this order, checking each one's figures before the next):
  lant init -o CONFIG --raw DIR --stations FILE --out DIR --velocity SLOW FAST
  lant clean -c CONFIG       station-days cleaned and decimated     -> OUT/1_clean
  lant cc -c CONFIG          cross-correlation and stack            -> OUT/2_cc
  lant disp -c CONFIG        reference curve and phase velocities   -> OUT/3_disp
  lant inv -c CONFIG         starting model and 3-D inversion       -> OUT/4_inv
  lant checker -c CONFIG --size KM [--flip-depth KM]   resolution test -> OUT/5_checker_*
"""


def main(argv=None):
    p = argparse.ArgumentParser(prog="lant", description="Ambient-noise surface-wave tomography for dense arrays.",
                                epilog=STEPS, formatter_class=argparse.RawDescriptionHelpFormatter)  # fmt: skip
    p.add_argument("--version", action="version", version="lant %s" % __version__)
    sub = p.add_subparsers(dest="cmd", metavar="COMMAND")
    i = sub.add_parser("init", help="write a configuration with the values suggested for an array")
    i.add_argument("-o", "--config", required=True, help="configuration file to write")
    i.add_argument("--raw", required=True, help="folder of the continuous records")
    i.add_argument("--stations", required=True, help="station table")
    i.add_argument("--out", required=True, help="output folder")
    i.add_argument("--velocity", nargs=2, type=float, required=True, metavar=("SLOW", "FAST"),
                   help="slowest and fastest Rayleigh phase velocity expected (km/s)")  # fmt: skip
    i.add_argument("--channel", default="*Z", help="channel codes of the vertical component (default *Z)")
    for name, text in (("clean", "step 1: clean the records"), ("cc", "step 2: cross-correlation"),
                       ("disp", "step 3: phase velocities"), ("inv", "step 4: 3-D inversion")):  # fmt: skip
        s = sub.add_parser(name, help=text)
        s.add_argument("-c", "--config", required=True)
    k = sub.add_parser("checker", help="checkerboard resolution test of the inversion")
    k.add_argument("-c", "--config", required=True)
    k.add_argument("--size", type=float, required=True, help="cell size (km)")
    k.add_argument("--flip-depth", type=float, default=None, help="the pattern changes sign below this depth (km)")
    e = sub.add_parser("export-sac", help="write the stacked correlations as SAC files")
    e.add_argument("-c", "--config", required=True)
    e.add_argument("-o", "--out", required=True, help="folder for the SAC files")
    sub.add_parser("check", help="check that the programs lant needs are installed")
    a = p.parse_args(argv)
    if a.cmd is None:
        p.print_help()
        return 0
    try:
        return _run(a)
    except KeyboardInterrupt:
        print("interrupted", file=sys.stderr)
        return 130
    except (ValueError, FileNotFoundError, RuntimeError) as exc:
        print("lant %s: %s" % (a.cmd, exc), file=sys.stderr)
        return 1


def _run(a):
    if a.cmd == "check":
        return check()
    if a.cmd == "init":
        from . import init, resources

        init.run(a.config, a.raw, a.stations, a.out, a.velocity, a.channel, resources.cores()[0])
        return 0
    from . import config

    cfg = config.load(a.config)
    if a.cmd == "clean":
        from . import clean

        clean.run(cfg)
    elif a.cmd == "cc":
        from . import correlate

        correlate.run(cfg)
    elif a.cmd == "disp":
        from . import disp_run

        disp_run.run(cfg)
    elif a.cmd == "inv":
        from .inversion import run

        run.run(cfg)
    elif a.cmd == "checker":
        from .inversion import checker

        checker.run(cfg, a.size, a.flip_depth)
    elif a.cmd == "export-sac":
        from . import export

        export.run(cfg, a.out)
    return 0


def check():
    import importlib

    ok = True
    for mod in ("numpy", "scipy", "numba", "h5py", "obspy", "matplotlib", "yaml"):
        try:
            m = importlib.import_module(mod)
            print("  %-10s %s" % (mod, getattr(m, "__version__", "")))
        except ImportError:
            print("  %-10s MISSING" % mod)
            ok = False
    from .plotting import gmt

    v = gmt.available()
    print("  %-10s %s" % ("gmt", v or "MISSING (maps are not drawn; conda install -c conda-forge gmt)"))
    return 0 if ok else 1
