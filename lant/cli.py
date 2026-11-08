"""Command line: `lant <command> [-c config.yaml] ...`  (`python -m lant` works as well).

  lant init surface|line -o my_project.yaml   write an example configuration to edit
  lant check                                  check the installation (Python packages, GMT, DSurfTomo, CC-FJpy)
  lant validate -c CFG                        read and check the configuration and the station table
  lant run STEP [STEP ...] -c CFG             run steps (clean cc stack notch reference disp1 disp2 inversion run evaluate figures)
  lant run all -c CFG [--from STEP] [--to STEP]
  lant run run -c CFG --what main|checks|all  which DSurfTomo projects to run (default all)
  lant status -c CFG                          the state of every step
  lant inspect TOPIC -c CFG                   key numbers of a step (TOPIC: data or a step name), as JSON
  lant plot [NAME ...] -c CFG                 draw the figures (default all; names: see lant.plotting.figures)
"""
import argparse
import json
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _cfg(path):
    from . import config
    if not path:
        sys.exit("give the project configuration with -c CONFIG.yaml (create one with `lant init surface -o CONFIG.yaml`)")
    try:
        return config.load(path)
    except (config.ConfigError, FileNotFoundError) as exc:
        sys.exit("configuration error: %s" % exc)


def cmd_init(a):
    src = os.path.join(ROOT, "config", "%s.example.yaml" % a.kind)
    if os.path.exists(a.output) and not a.force:
        sys.exit("%s exists (use --force to overwrite)" % a.output)
    shutil.copyfile(src, a.output)
    print("wrote %s: edit project.output_dir, data.raw_dir and data.station_file, then run `lant validate -c %s`"
          % (a.output, a.output))


def cmd_check(a):
    from .tools import external_path
    ok = True
    rows = []
    for mod in ("numpy", "scipy", "h5py", "obspy", "matplotlib", "yaml"):
        try:
            m = __import__(mod)
            rows.append(("python: %s" % mod, getattr(m, "__version__", "ok")))
        except ImportError:
            rows.append(("python: %s" % mod, "MISSING (pip install -e . or conda env create -f environment.yml)"))
            ok = False
    from .plotting.gmt import available
    v = available()
    rows.append(("GMT 6 (maps, sections)", v or "MISSING (conda install -c conda-forge gmt); figures of maps/sections cannot be drawn"))
    for name, rel in (("DSurfTomo", "bin/DSurfTomo"), ("fwd1d (1-D starting model)", "bin/fwd1d"),
                      ("pDSurfTomo (optional)", "bin/pDSurfTomo")):
        p = external_path(rel)
        rows.append((name, p if os.path.exists(p) else "missing: bash tools/install_external.sh"))
        if not os.path.exists(p) and "optional" not in name:
            ok = False
    try:
        from .reference import import_ccfj
        import_ccfj()
        rows.append(("CC-FJpy (F-J reference)", "ok"))
    except Exception:
        rows.append(("CC-FJpy (F-J reference)", "missing: bash tools/install_external.sh ccfj (or borrow a reference curve)"))
    w = max(len(r[0]) for r in rows)
    for k, v in rows:
        print("%-*s  %s" % (w, k, v))
    if not ok:
        sys.exit(1)


def cmd_validate(a):
    from . import data as D, stations as ST
    cfg = _cfg(a.config)
    print("configuration %s: ok (array type %s, output %s)" % (cfg["_config_file"], cfg["project"]["array_type"],
                                                              cfg["project"]["output_dir"]))
    try:
        rows, how = ST.read_table(cfg["data"]["station_file"])
        print("station table: %d stations read (%s)" % (len(rows), how))
    except (ST.StationError, TypeError, FileNotFoundError) as exc:
        print("station table: PROBLEM: %s" % exc)
    files = D.find_files(cfg["data"]["raw_dir"], cfg["data"]["file_patterns"], cfg["data"]["ignore_patterns"])
    print("raw directory %s: %d candidate files; `lant run clean` reads their headers, matches the names and derives "
          "the 'auto' parameters" % (cfg["data"]["raw_dir"], len(files)))


def cmd_run(a):
    from . import pipeline as PL
    cfg = _cfg(a.config)
    steps = a.steps
    if steps == ["all"]:
        i0 = PL.STEPS.index(a.from_step) if a.from_step else 0
        i1 = PL.STEPS.index(a.to_step) + 1 if a.to_step else len(PL.STEPS)
        steps = PL.STEPS[i0:i1]
    for s in steps:
        if s not in PL.STEPS:
            sys.exit("unknown step %s (steps: %s)" % (s, " ".join(PL.STEPS)))
    for s in steps:
        kw = {"what": a.what} if s == "run" else {}
        try:
            PL.run_step(cfg, s, **kw)
        except Exception as exc:
            if a.debug:
                raise
            sys.exit("step %s failed: %s\n(details: %s/records/run_log.jsonl; rerun with --debug for the traceback)"
                     % (s, exc, cfg["project"]["output_dir"]))


def cmd_status(a):
    from .inspect import status
    cfg = _cfg(a.config)
    for k, v in status(cfg).items():
        print("%-10s %-8s %s %s" % (k, v["state"], v.get("time", ""), v.get("error", "")))


def cmd_inspect(a):
    from .inspect import inspect
    cfg = _cfg(a.config)
    print(json.dumps(inspect(cfg, a.topic), indent=1, default=str, ensure_ascii=False))


def cmd_plot(a):
    from . import pipeline as PL
    cfg = _cfg(a.config)
    try:
        s = PL.run_step(cfg, "figures", which=a.names or None)
    except Exception as exc:
        sys.exit("figures failed: %s" % exc)
    for k, v in s["not_drawn"].items():
        print("not drawn: %s (%s)" % (k, v))


def main(argv=None):
    p = argparse.ArgumentParser(prog="lant", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    q = sub.add_parser("init")
    q.add_argument("kind", choices=("surface", "line"))
    q.add_argument("-o", "--output", required=True)
    q.add_argument("--force", action="store_true")
    q.set_defaults(fn=cmd_init)
    sub.add_parser("check").set_defaults(fn=cmd_check)
    for name, fn in (("validate", cmd_validate), ("status", cmd_status)):
        q = sub.add_parser(name)
        q.add_argument("-c", "--config")
        q.set_defaults(fn=fn)
    q = sub.add_parser("run")
    q.add_argument("steps", nargs="+")
    q.add_argument("-c", "--config")
    q.add_argument("--from", dest="from_step")
    q.add_argument("--to", dest="to_step")
    q.add_argument("--what", default="all", choices=("main", "checks", "all"))
    q.add_argument("--debug", action="store_true")
    q.set_defaults(fn=cmd_run)
    q = sub.add_parser("inspect")
    q.add_argument("topic")
    q.add_argument("-c", "--config")
    q.set_defaults(fn=cmd_inspect)
    q = sub.add_parser("plot")
    q.add_argument("names", nargs="*")
    q.add_argument("-c", "--config")
    q.set_defaults(fn=cmd_plot)
    a = p.parse_args(argv)
    a.fn(a)
