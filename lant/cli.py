"""Command line: `lant <command> [-c config.yaml] ...`  (`python -m lant` works as well).

  lant check                                  check the installation (Python packages and versions, GMT, DSurfTomo)
  lant run all --raw DIR --stations FILE --out DIR
                                              the shortest way: writes DIR(out)/config.yaml and runs the chain
  lant init -o CFG [--raw DIR --stations FILE --out DIR]
                                              write a configuration: the three required paths, every other key
                                              with its default and what it rests on
  lant validate -c CFG                        read and check the configuration and the station table
  lant plan -c CFG                            what every step will do, with which values and where they come from
  lant run all -c CFG                         run the chain from the first step that does not stand; stops where a
                                              judgement is yours (--from STEP: run again from there; --to STEP;
                                              --rebuild: move the old results of the steps run to *_old first)
  lant run STEP [STEP ...] -c CFG             run single steps (clean cc stack reference disp inversion run evaluate figures)
  lant run run -c CFG --what main|checks|all  which DSurfTomo projects to run (default all)
  lant confirm STEP -c CFG                    confirm what a step left for you to judge (stack: --exclude-candidates)
  lant status -c CFG [--json]                 the state of every step, why, and the next command
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
        sys.exit("give the project configuration with -c CONFIG.yaml (create one with `lant init -o CONFIG.yaml`), or --raw DIR --stations FILE --out DIR")
    try:
        return config.load(path)
    except (config.ConfigError, FileNotFoundError) as exc:
        sys.exit("configuration error: %s" % exc)


def cmd_init(a):
    from . import params
    if os.path.exists(a.output) and not a.force:
        sys.exit("%s exists (use --force to overwrite)" % a.output)
    given = {k: os.path.abspath(v) for k, v in (("raw_dir", a.raw), ("station_file", a.stations), ("output_dir", a.out)) if v}
    with open(a.output, "w") as fh:
        fh.write(params.template(full=a.full, **given))
    print("wrote %s\n  1. %s\n  2. lant validate -c %s\n  3. lant plan -c %s\n  4. lant run all -c %s"
          % (a.output, "the three paths are filled in" if len(given) == 3 else "fill in the three paths at its top",
             a.output, a.output, a.output))


def _quick_config(a):
    """`lant run ... --raw --stations --out` without -c: the configuration is <out>/config.yaml, written on the
    first call from the three paths and used as it is afterwards."""
    from . import params
    if not (a.raw and a.stations and a.out):
        sys.exit("give -c CONFIG, or all of --raw DIR --stations FILE --out DIR")
    from . import stations as ST
    out = os.path.abspath(a.out)
    try:
        lines = sorted({r["line"] for r in ST.read_table(os.path.abspath(a.stations))[0] if r.get("line")})
    except Exception:
        lines = []
    paths_ = []
    for line in (lines if len(lines) > 1 else [None]):       # a table with several lines: one project per line
        d = os.path.join(out, "line_%s" % line) if line else out
        path = os.path.join(d, "config.yaml")
        if os.path.exists(path):
            print("[lant] using the configuration %s (written earlier; edit it to change parameters)" % path)
        else:
            os.makedirs(d, exist_ok=True)
            with open(path, "w") as fh:
                fh.write(params.template(os.path.abspath(a.raw), os.path.abspath(a.stations), ".", line))
            print("[lant] wrote the configuration %s" % path)
        paths_.append(path)
    return paths_


MIN_VERSIONS = {"numpy": "1.19", "scipy": "1.5", "h5py": "2.10", "obspy": "1.2", "matplotlib": "3.5", "yaml": "6.0",
                "openpyxl": "3.1"}                   # as in pyproject.toml
MIN_PYTHON = (3, 8)


def _vtuple(v):
    import re
    return tuple(int(x) for x in re.findall(r"\d+", str(v))[:3])


def cmd_check(a):
    from .tools import external_path
    ok = True
    rows = []
    py = sys.version_info[:3]
    rows.append(("python", "%d.%d.%d%s" % (py + ("" if py >= MIN_PYTHON else "  TOO OLD: needs >= %d.%d" % MIN_PYTHON,))))
    ok &= py >= MIN_PYTHON
    for mod, need in MIN_VERSIONS.items():
        try:
            m = __import__(mod)
            v = getattr(m, "__version__", None)
            if v is not None and _vtuple(v) < _vtuple(need):
                rows.append(("python: %s" % mod, "%s  TOO OLD: needs >= %s (pip install -U %s)"
                             % (v, need, "pyyaml" if mod == "yaml" else mod)))
                ok = False
            else:
                rows.append(("python: %s" % mod, v or "ok"))
        except ImportError:
            rows.append(("python: %s" % mod, "MISSING: pip install -e .   (or: conda env create -f environment.yml)"))
            ok = False
    from .plotting.gmt import available
    v = available()
    rows.append(("GMT 6 (maps, sections)", v or "MISSING: conda install -c conda-forge gmt   (maps and sections "
                                                "cannot be drawn without it)"))
    for name, rel, what in (("DSurfTomo", "bin/DSurfTomo", "dsurftomo"), ("fwd1d (1-D starting model)", "bin/fwd1d", "dsurftomo"),
                            ("pDSurfTomo (optional)", "bin/pDSurfTomo", "pdsurftomo")):
        p = external_path(rel)
        good = os.path.isfile(p) and os.access(p, os.X_OK)
        rows.append((name, p if good else ("not executable: chmod +x %s" % p if os.path.exists(p)
                                           else "MISSING: lant install-external %s" % what)))
        if not good and "optional" not in name:
            ok = False
    w = max(len(r[0]) for r in rows)
    for k, v in rows:
        print("%-*s  %s" % (w, k, v))
    if not ok:
        sys.exit("\nthe installation is incomplete: run the commands printed after MISSING / TOO OLD, then `lant check` again")
    print("\nthe installation is complete; next: lant run all --raw RAW_DIR --stations STATION_TABLE --out OUTPUT_DIR")


def cmd_validate(a):
    from . import data as D, stations as ST
    cfg = _cfg(a.config)
    print("configuration %s: ok (array type %s, output %s)" % (cfg["_config_file"], cfg["project"]["array_type"],
                                                              cfg["project"]["output_dir"]))
    problems = 0
    try:
        rows, how = ST.read_table(cfg["data"]["station_file"])
        print("station table: %d stations read (%s)" % (len(rows), how))
    except (ST.StationError, TypeError, FileNotFoundError) as exc:
        print("station table: PROBLEM: %s" % exc)
        problems += 1
    try:
        files = D.find_files(cfg["data"]["raw_dir"], cfg["data"]["file_patterns"], cfg["data"]["ignore_patterns"])
        print("raw directory %s: %d candidate files; `lant run clean` reads their headers, matches the names and derives "
              "the 'auto' parameters" % (cfg["data"]["raw_dir"], len(files)))
        if not files:
            print("raw directory: PROBLEM: no candidate files")
            problems += 1
    except D.DataError as exc:
        print("raw directory: PROBLEM: %s" % exc)
        problems += 1
    if problems:
        sys.exit(1)


STEP_DIRS = {"clean": ["01_clean"], "cc": ["02_ccf"], "stack": ["03_stack"], "reference": ["04_reference"],
             "disp": ["05_dispersion"], "inversion": ["06_inversion"], "run": [], "evaluate": [], "figures": ["figures"]}


def _rebuild(cfg, steps):
    """Keep prior directories in a unique backup; never replace an older backup."""
    import datetime
    from . import artifacts, paths, records
    from .provenance import Lock
    artifacts.check_roots(cfg)
    with Lock(records.rec_dir(cfg, create=True)):
        before = artifacts.snapshot(cfg)
        backup = os.path.join(records.rec_dir(cfg), "backups", datetime.datetime.now().strftime("%Y%m%dT%H%M%S.%f"))
        try:
            for step in steps:
                for directory in STEP_DIRS[step]:
                    source = paths.P(cfg, directory)
                    if os.path.isdir(source):
                        os.makedirs(backup, exist_ok=True)
                        os.rename(source, os.path.join(backup, directory))
                        print("[lant] backup: %s" % backup, file=sys.stderr)
        finally:
            artifacts.register(cfg, before)


def _accept_impact(message, yes):
    print(message, file=sys.stderr)
    if yes:
        return
    if not sys.stdin.isatty():
        raise RuntimeError("noninteractive execution requires --yes after reviewing the impact")
    if input("Continue? Type yes: ").strip() != "yes":
        raise RuntimeError("cancelled; no changes made")


def cmd_run(a):
    from . import manifest, report
    from . import pipeline as PL
    from .confirm import NeedsConfirmation
    if not a.config and (a.raw or a.stations or a.out):
        cfgs = _quick_config(a)
        if len(cfgs) > 1:                           # one project per line, one after the other
            import subprocess
            codes = []
            for c in cfgs:
                print("\n[lant] ===== project %s" % c, flush=True)
                codes.append(subprocess.call([sys.executable, "-m", "lant", "run"] + a.steps + ["-c", c]))
            sys.exit(max(codes))
        a.config = cfgs[0]
    cfg = _cfg(a.config)
    cfg["_verbose"] = a.verbose
    chain = a.steps == ["all"]
    steps = a.steps
    if chain:
        for s in (a.from_step, a.to_step):
            if s and s not in PL.STEPS:
                sys.exit("unknown step %s (steps: %s)" % (s, " ".join(PL.STEPS)))
        i1 = PL.STEPS.index(a.to_step) + 1 if a.to_step else len(PL.STEPS)
        steps = PL.STEPS[:i1]
    for s in steps:
        if s not in PL.STEPS:
            sys.exit("unknown step %s (steps: %s)" % (s, " ".join(PL.STEPS)))
    forced = set(steps)                             # single steps named by the user always run
    if chain:                                       # the chain leaves steps that stand alone, unless told otherwise
        forced = set(PL.STEPS[PL.STEPS.index(a.from_step):]) if a.from_step else (set(steps) if a.rebuild else set())
    if a.rebuild or a.from_step:
        _accept_impact("Recompute steps: " + ", ".join(s for s in steps if s in forced) +
                       ". Retained results will be refreshed; --rebuild keeps a unique working backup.", a.yes)
    if a.rebuild:
        _rebuild(cfg, [s for s in steps if s in forced])
    for s in steps:
        kw = {"what": a.what} if s == "run" else {}
        try:
            if chain and s not in forced:
                st = manifest.status(report.resolved(cfg))[s]
                if st["state"] == "done":
                    print("[lant] step %s stands (run %s): not run again" % (s, st.get("time", "")), flush=True)
                    continue
                if st["state"] == "to confirm":
                    manifest.gate(report.resolved(cfg), PL.STEPS[PL.STEPS.index(s) + 1])
            PL.run_step(cfg, s, **kw)
        except NeedsConfirmation as exc:
            print("\nSTOPPED for your judgement: %s" % exc)
            sys.exit(3)
        except Exception as exc:
            if a.debug:
                raise
            sys.exit("step %s failed: %s\n(details: %s/.lant/run_log.jsonl; the traceback: add --debug)"
                     % (s, exc, cfg["project"]["output_dir"]))
    if chain:
        rc = report.resolved(cfg)
        cmd, why = manifest.next_action(rc, manifest.status(rc))
        print("\n[lant] %s%s" % (why, "\n  next: " + cmd if cmd else ""))
        if cmd and cmd.startswith("lant confirm"):
            sys.exit(3)


def cmd_confirm(a):
    from . import auto, confirm, manifest, report
    cfg = _cfg(a.config)
    cfg = auto.apply(cfg)
    try:
        entry, rerun = confirm.confirm(cfg, a.step, exclude_candidates=a.exclude_candidates,
                                        reason=a.reason, invocation="terminal" if sys.stdin.isatty() else "noninteractive")
    except (RuntimeError, ValueError) as exc:
        sys.exit("confirm failed: %s" % exc)
    print("confirmed %s: written into %s (block 'confirmed:')" % (a.step, cfg["_config_file"]))
    for t in entry.get("_judged", []):
        print("  judged: %s" % t)
    for k, v in entry.items():
        if not k.startswith("_"):
            print("  %s = %s" % (k, json.dumps(v, default=str)))
    rc = report.resolved(cfg)
    cmd, why = manifest.next_action(rc, manifest.status(rc))
    print("%s%s" % (why, "\n  next: " + cmd if cmd else ""))


def cmd_status(a):
    from .inspect import status
    cfg = _cfg(a.config)
    st = status(cfg)
    if a.json:
        print(json.dumps(st, indent=1, default=str, ensure_ascii=False))
        return
    for k, v in st["steps"].items():
        print("%-10s %-12s %-20s %s" % (k, v["state"].upper() if v["state"] != "done" else "done", v.get("time", ""),
                                        v.get("reason", "")))
    for key in ("unregistered", "missing", "modified"):
        files = st["artifacts"][key]
        if files:
            print("%s files: %d (%s)" % (key, len(files), ", ".join(files[:5])))
    print("\n%s%s" % (st["next"]["why"], "\n  next: " + st["next"]["command"] if st["next"]["command"] else ""))


def cmd_plan(a):
    from . import auto, report
    cfg = auto.apply(_cfg(a.config))
    print(report.plan(cfg))


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


def _external_program(value):
    if value not in ("dsurftomo", "pdsurftomo"):
        raise argparse.ArgumentTypeError("choose dsurftomo or pdsurftomo")
    return value


def cmd_install_external(a):
    import subprocess
    from .tools import external_path
    script = os.path.join(os.path.dirname(__file__), "_external", "install.sh")
    for exe in ("git", "make", "gfortran", "bash"):
        if not shutil.which(exe):
            raise RuntimeError("external installation needs %s on PATH" % exe)
    env = dict(os.environ, LANT_EXTERNAL=external_path())
    subprocess.run(["bash", script] + a.programs, env=env, check=True)


def cmd_note(a):
    from . import artifacts, records
    from .delivery import atomic_text, initialize
    from .provenance import Lock
    cfg = _cfg(a.config)
    artifacts.check_roots(cfg)
    with Lock(records.rec_dir(cfg, create=True)):
        before = artifacts.snapshot(cfg)
        try:
            initialize(cfg)
            file = os.path.join(cfg["project"]["output_dir"], "notes.md")
            prior = open(file).read() if os.path.exists(file) else "# Notes and decisions\n"
            atomic_text(file, prior + "\n- %s: %s\n" % (records._now(), a.text))
            records.log(cfg, event="note", text=a.text)
        finally:
            artifacts.register(cfg, before)
    print(file)


def cmd_clean_work(a):
    from .artifacts import clean
    cfg = _cfg(a.config)
    selection = a.cache or ["archive", "cc"]
    plan = clean(cfg, selection)
    print(json.dumps(plan, indent=1))
    if a.execute:
        _accept_impact("Delete %d registered cache files (%.3f GB). %s" %
                       (len(plan["files"]), plan["bytes"] / 1e9, plan["impact"]), a.yes)
        clean(cfg, selection, execute=True)


def main(argv=None):
    p = argparse.ArgumentParser(prog="lant", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    q = sub.add_parser("init")
    q.add_argument("-o", "--output", required=True)
    for f, h in (("--raw", "raw data directory"), ("--stations", "station table"), ("--out", "output directory")):
        q.add_argument(f, help=h)
    q.add_argument("--force", action="store_true")
    q.add_argument("--full", action="store_true", help="write all keys with their defaults and sources")
    q.set_defaults(fn=cmd_init)
    sub.add_parser("check").set_defaults(fn=cmd_check)
    for name, fn in (("validate", cmd_validate), ("status", cmd_status), ("plan", cmd_plan)):
        q = sub.add_parser(name)
        q.add_argument("-c", "--config")
        if name == "status":
            q.add_argument("--json", action="store_true", help="machine-readable: steps, reasons, next command")
        q.set_defaults(fn=fn)
    q = sub.add_parser("run")
    q.add_argument("steps", nargs="+")
    q.add_argument("-c", "--config")
    for f, h in (("--raw", "raw data directory"), ("--stations", "station table"), ("--out", "output directory")):
        q.add_argument(f, help=h + " (instead of -c: the configuration is OUT/config.yaml)")
    q.add_argument("--from", dest="from_step")
    q.add_argument("--to", dest="to_step")
    q.add_argument("--what", default="all", choices=("main", "checks", "all"))
    q.add_argument("-v", "--verbose", action="store_true", help="show detailed step output")
    q.add_argument("--debug", action="store_true")
    q.add_argument("--yes", action="store_true", help="acknowledge the displayed rebuild impact")
    q.add_argument("--rebuild", action="store_true", help="move the old results of the steps run to *_old first")
    q.set_defaults(fn=cmd_run)
    q = sub.add_parser("confirm")
    q.add_argument("step", choices=("clean", "stack", "reference", "inversion"))
    q.add_argument("-c", "--config")
    q.add_argument("--reason", required=True, help="user decision and evidence (terminal use is not identity verification)")
    q.add_argument("--exclude-candidates", action="store_true", help="stack: remove the candidate station-days/-hours")
    q.set_defaults(fn=cmd_confirm)
    q = sub.add_parser("inspect")
    q.add_argument("topic")
    q.add_argument("-c", "--config")
    q.set_defaults(fn=cmd_inspect)
    q = sub.add_parser("plot")
    q.add_argument("names", nargs="*")
    q.add_argument("-c", "--config")
    q.set_defaults(fn=cmd_plot)
    q = sub.add_parser("install-external", help="build fixed third-party programs in a writable location")
    q.add_argument("programs", nargs="*", type=_external_program, metavar="PROGRAM")
    q.set_defaults(fn=cmd_install_external)
    q = sub.add_parser("note", help="append a project note")
    q.add_argument("text")
    q.add_argument("-c", "--config")
    q.set_defaults(fn=cmd_note)
    q = sub.add_parser("clean-work", help="preview registered cache cleanup")
    q.add_argument("-c", "--config")
    q.add_argument("--cache", action="append", choices=("archive", "cc"))
    q.add_argument("--execute", action="store_true", help="perform cleanup (default is preview only)")
    q.add_argument("--yes", action="store_true")
    q.set_defaults(fn=cmd_clean_work)
    a = p.parse_args(argv)
    try:
        a.fn(a)
    except KeyboardInterrupt:
        sys.exit("interrupted; state and logs are retained; inspect lant status before continuing")
    except (RuntimeError, ValueError, OSError) as exc:
        sys.exit(str(exc))
