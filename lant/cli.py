"""lant: ambient-noise surface-wave tomography for dense arrays, from raw records to a 3-D Vs model.

Getting started
  lant check                                check the installation (Python packages, GMT, DSurfTomo)
  lant install-external [PROGRAM ...]       build DSurfTomo and fwd1d (dsurftomo) and/or pDSurfTomo (pdsurftomo)
  lant demo [--out DIR]                     run a small synthetic array with a known answer (a few minutes)
  lant run all --raw DIR --stations FILE --out DIR
                                            your data: writes OUT/config.yaml and runs the whole chain

Working with a project (CFG: its configuration file)
  lant init -o CFG [--raw DIR --stations FILE --out DIR] [--full]
                                            write a configuration (compact; --full: every key with its basis)
  lant validate -c CFG                      check the configuration, the station table and the raw directory
  lant plan -c CFG                          what every step will do, with which values and why
  lant run all -c CFG                       run from the first step that does not stand; stops where you judge
                                            (--from STEP, --to STEP, --rebuild; -v: detailed output)
  lant run STEP [STEP ...] -c CFG           single steps: clean cc stack reference disp inversion run evaluate figures
  lant status -c CFG [--json]               the state of every step, why, and the next command
  lant confirm STEP -c CFG --reason TEXT    record your decision at a judgement point
  lant note "TEXT" -c CFG                   append a note to OUT/notes.md (also while a step runs)
  lant inspect TOPIC -c CFG                 key numbers of a step (TOPIC: data or a step name), as JSON
  lant plot [NAME ...] -c CFG               draw the figures again
  lant clean-work -c CFG [--execute --yes]  preview (or remove) rebuildable caches

Results: open OUT/REPORT.md.  `lant COMMAND --help` lists the options of a command.
Documentation: docs/quickstart.md and docs/tutorial.md (https://github.com/Shubin-Wei/local_ANT).
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
        sys.exit(
            "give the project configuration with -c CONFIG.yaml (create one with `lant init -o CONFIG.yaml`), "
            "or --raw DIR --stations FILE --out DIR"
        )
    if not os.path.isfile(path):
        sys.exit("configuration file not found: %s\n  create one with: lant init -o %s" % (path, path))
    try:
        return config.load(path)
    except config.ConfigError as exc:
        sys.exit("configuration error in %s: %s" % (path, exc))


def cmd_init(a):
    from . import params

    if os.path.exists(a.output) and not a.force:
        sys.exit("%s exists (use --force to overwrite)" % a.output)
    given = {
        k: os.path.abspath(v) for k, v in (("raw_dir", a.raw), ("station_file", a.stations), ("output_dir", a.out)) if v
    }
    with open(a.output, "w") as fh:
        fh.write(params.template(full=a.full, **given))
    print(
        "wrote %s\n  1. %s\n  2. lant validate -c %s\n  3. lant plan -c %s\n  4. lant run all -c %s"
        % (
            a.output,
            "the three paths are filled in" if len(given) == 3 else "fill in the three paths at its top",
            a.output,
            a.output,
            a.output,
        )
    )


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
    for line in lines if len(lines) > 1 else [None]:  # a table with several lines: one project per line
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


MIN_VERSIONS = {
    "numpy": "1.19",
    "scipy": "1.5",
    "h5py": "2.10",
    "obspy": "1.2",
    "matplotlib": "3.5",
    "yaml": "6.0",
    "openpyxl": "3.1",
}  # as in pyproject.toml
MIN_PYTHON = (3, 8)


def _vtuple(v):
    import re

    return tuple(int(x) for x in re.findall(r"\d+", str(v))[:3])


def cmd_check(a):
    from .tools import external_path

    ok = True
    rows = []
    py = sys.version_info[:3]
    rows.append(
        ("python", "%d.%d.%d%s" % (py + ("" if py >= MIN_PYTHON else "  TOO OLD: needs >= %d.%d" % MIN_PYTHON,)))
    )
    ok &= py >= MIN_PYTHON
    for mod, need in MIN_VERSIONS.items():
        try:
            m = __import__(mod)
            v = getattr(m, "__version__", None)
            if v is not None and _vtuple(v) < _vtuple(need):
                rows.append(
                    (
                        "python: %s" % mod,
                        "%s  TOO OLD: needs >= %s (pip install -U %s)" % (v, need, "pyyaml" if mod == "yaml" else mod),
                    )
                )
                ok = False
            else:
                rows.append(("python: %s" % mod, v or "ok"))
        except ImportError:
            rows.append(("python: %s" % mod, "MISSING: pip install -e .   (or: conda env create -f environment.yml)"))
            ok = False
    from .plotting.gmt import available

    v = available()
    rows.append(
        (
            "GMT 6 (maps, sections)",
            v or "MISSING: conda install -c conda-forge gmt   (maps and sections cannot be drawn without it)",
        )
    )
    for name, rel, what in (
        ("DSurfTomo", "bin/DSurfTomo", "dsurftomo"),
        ("fwd1d (1-D starting model)", "bin/fwd1d", "dsurftomo"),
        ("pDSurfTomo (optional)", "bin/pDSurfTomo", "pdsurftomo"),
    ):
        p = external_path(rel)
        good = os.path.isfile(p) and os.access(p, os.X_OK)
        rows.append(
            (
                name,
                p
                if good
                else (
                    "not executable: chmod +x %s" % p
                    if os.path.exists(p)
                    else "MISSING: lant install-external %s" % what
                ),
            )
        )
        if not good and "optional" not in name:
            ok = False
    w = max(len(r[0]) for r in rows)
    for k, v in rows:
        print("%-*s  %s" % (w, k, v))
    if not ok:
        sys.exit(
            "\nthe installation is incomplete: run the commands printed after MISSING / TOO OLD, then `lant check` again"
        )
    print("\nthe installation is complete; next: lant run all --raw RAW_DIR --stations STATION_TABLE --out OUTPUT_DIR")


def cmd_validate(a):
    from . import data as D
    from . import stations as ST

    cfg = _cfg(a.config)
    placeholders = [
        k
        for k, v in (
            ("project.output_dir", cfg["project"]["output_dir"]),
            ("data.raw_dir", cfg["data"]["raw_dir"]),
            ("data.station_file", cfg["data"]["station_file"]),
        )
        if "/path/to/" in str(v)
    ]
    if placeholders:
        sys.exit(
            "%s still has the template paths of %s: write your own paths at the top of the file"
            % (cfg["_config_file"], ", ".join(placeholders))
        )
    print(
        "configuration %s: read (array type %s, output %s)"
        % (cfg["_config_file"], cfg["project"]["array_type"], cfg["project"]["output_dir"])
    )
    problems = 0
    try:
        rows, how = ST.read_table(cfg["data"]["station_file"])
        print("station table: %d stations read (%s)" % (len(rows), how))
    except (ST.StationError, TypeError, FileNotFoundError) as exc:
        print("station table: PROBLEM: %s" % exc)
        problems += 1
    try:
        files = D.find_files(cfg["data"]["raw_dir"], cfg["data"]["file_patterns"], cfg["data"]["ignore_patterns"])
        print(
            "raw directory %s: %d candidate files; `lant run clean` reads their headers, matches the names and derives "
            "the 'auto' parameters" % (cfg["data"]["raw_dir"], len(files))
        )
        if not files:
            print("raw directory: PROBLEM: no candidate files")
            problems += 1
    except D.DataError as exc:
        print("raw directory: PROBLEM: %s" % exc)
        problems += 1
    if problems:
        sys.exit("%d problem(s) to fix before running" % problems)
    print("ok: next `lant plan -c %s` (what will be done) or `lant run all -c %s`" % (a.config, a.config))


STEP_DIRS = {
    "clean": ["01_clean"],
    "cc": ["02_ccf"],
    "stack": ["03_stack"],
    "reference": ["04_reference"],
    "disp": ["05_dispersion"],
    "inversion": ["06_inversion"],
    "run": [],
    "evaluate": [],
    "figures": ["figures"],
}


def _rebuild(cfg, steps):
    """Keep prior directories in a unique backup; never replace an older backup."""
    import datetime

    from . import artifacts, paths, records

    with artifacts.writing(cfg, "rebuild"):
        backup = os.path.join(records.rec_dir(cfg), "backups", datetime.datetime.now().strftime("%Y%m%dT%H%M%S.%f"))
        for step in steps:
            for directory in STEP_DIRS[step]:
                source = paths.P(cfg, directory)
                if os.path.isdir(source):
                    os.makedirs(backup, exist_ok=True)
                    os.rename(source, os.path.join(backup, directory))
                    print("[lant] backup: %s" % backup, file=sys.stderr)


def _accept_impact(message, yes):
    print(message, file=sys.stderr)
    if yes:
        return
    if not sys.stdin.isatty():
        raise RuntimeError("noninteractive execution requires --yes after reviewing the impact")
    sys.stderr.write("Continue? Type yes: ")  # stdout may carry JSON
    sys.stderr.flush()
    if sys.stdin.readline().strip() != "yes":
        raise RuntimeError("cancelled; no changes made")


def cmd_run(a):
    from . import manifest, report
    from . import pipeline as PL
    from .confirm import NeedsConfirmation

    if not a.config and (a.raw or a.stations or a.out):
        cfgs = _quick_config(a)
        if len(cfgs) > 1:  # one project per line, one after the other
            import subprocess

            codes = []
            for c in cfgs:
                print("\n[lant] ===== project %s" % c, flush=True)
                flags = ["--what", a.what]
                for option, value in (("--from", a.from_step), ("--to", a.to_step)):
                    if value:
                        flags += [option, value]
                for option, enabled in (
                    ("--debug", a.debug),
                    ("--verbose", a.verbose),
                    ("--rebuild", a.rebuild),
                    ("--yes", a.yes),
                ):
                    if enabled:
                        flags.append(option)
                codes.append(subprocess.call([sys.executable, "-m", "lant", "run"] + a.steps + ["-c", c] + flags))
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
    forced = set(steps)  # single steps named by the user always run
    if chain:  # the chain leaves steps that stand alone, unless told otherwise
        forced = set(PL.STEPS[PL.STEPS.index(a.from_step) :]) if a.from_step else (set(steps) if a.rebuild else set())
    if a.rebuild or a.from_step:
        _accept_impact(
            "Recompute steps: "
            + ", ".join(s for s in steps if s in forced)
            + ". Retained results will be refreshed; --rebuild keeps a unique working backup.",
            a.yes,
        )
    if a.rebuild:
        _rebuild(cfg, [s for s in steps if s in forced])
    for k, s in enumerate(steps):
        kw = {"what": a.what} if s == "run" else {}
        cfg["_position"] = "%d/%d" % (k + 1, len(steps)) if len(steps) > 1 else ""
        try:
            if chain and s not in forced:
                st = manifest.status(report.resolved(cfg))[s]
                if st["state"] == "done":
                    print(
                        "[lant] [%s] step %s stands (run %s): not run again"
                        % (cfg["_position"], s, st.get("time", "")),
                        flush=True,
                    )
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
            sys.exit(
                "step %s failed: %s\n  kept: the results of the earlier steps and the record of this one "
                "(%s/.lant/lant.log, run_log.jsonl)\n  next: fix the cause, then `lant run all -c %s` "
                "(the traceback: add --debug)" % (s, exc, cfg["project"]["output_dir"], a.config)
            )
    if chain:
        rc = report.resolved(cfg)
        cmd, why = manifest.next_action(rc, manifest.status(rc))
        print("\n[lant] %s%s" % (why, "\n  next: " + cmd if cmd else ""))
        if not cmd:
            from .delivery import finish_lines

            print("\n".join(finish_lines(rc)))
        if cmd and cmd.startswith("lant confirm"):
            sys.exit(3)


def cmd_confirm(a):
    from . import auto, confirm, manifest, report

    cfg = _cfg(a.config)
    cfg = auto.apply(cfg)
    try:
        entry, rerun = confirm.confirm(
            cfg,
            a.step,
            exclude_candidates=a.exclude_candidates,
            reason=a.reason,
            invocation="terminal" if sys.stdin.isatty() else "noninteractive",
        )
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
        print(
            "%-10s %-12s %-20s %s"
            % (k, v["state"].upper() if v["state"] != "done" else "done", v.get("time", ""), v.get("reason", ""))
        )
    for key in ("unregistered", "missing", "modified"):
        files = st["artifacts"][key]
        if files:
            print("%s files: %d (%s)" % (key, len(files), ", ".join(files[:5])))
    scratch = st["artifacts"].get("scratch") or {}
    if scratch.get("files"):
        print(
            "scratch/: %d temporary files (%.1f MB); delete them when no longer needed"
            % (scratch["files"], scratch["bytes"] / 1e6)
        )
    if st["artifacts"].get("deferred"):
        print("file audit: %s" % st["artifacts"]["deferred"])
    nxt = st["next"]
    print(
        "\n%s%s"
        % (
            nxt["why"],
            "\n  next%s: %s" % (" (the user's decision)" if nxt["requires_user"] else "", nxt["command"])
            if nxt["command"]
            else "",
        )
    )


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
    result = subprocess.run(["bash", script] + a.programs, env=env)
    if result.returncode:
        raise RuntimeError(
            "external compilation failed (exit %d); inspect the compiler output and verify "
            "git, make and gfortran before retrying" % result.returncode
        )


def cmd_note(a):
    from . import artifacts, records
    from .delivery import append_note, initialize
    from .provenance import active_writer

    cfg = _cfg(a.config)
    line = "%s: %s" % (records._now(), a.text)
    if os.path.isdir(records.rec_dir(cfg)) and active_writer(records.rec_dir(cfg)):
        print(append_note(cfg, line))  # a step is running: the note does not wait for it
        return
    with artifacts.writing(cfg, "note"):
        initialize(cfg)
        file = append_note(cfg, line)
        records.log(cfg, event="note", text=a.text)
    print(file)


def cmd_clean_work(a):
    from .artifacts import clean

    cfg = _cfg(a.config)
    selection = a.cache or ["archive", "cc"]
    plan = clean(cfg, selection)
    print(json.dumps(plan, indent=1))
    if a.execute:
        _accept_impact(
            "Delete %d registered cache files (%.3f GB). %s"
            % (len(plan["files"]), plan["bytes"] / 1e9, plan["impact"]),
            a.yes,
        )
        clean(cfg, selection, execute=True)


def _interrupt_on_termination():
    """SIGTERM and SIGHUP (kill, job schedulers, a closed terminal) take the same path as Ctrl-C: solver process
    groups are stopped, the step is recorded as failed and the lock released.  A signal already ignored (nohup)
    stays ignored; forked workers keep the default action, so that Pool.terminate still ends them."""
    import signal

    owner = os.getpid()

    def interrupt(signum, frame):
        if os.getpid() != owner:
            signal.signal(signum, signal.SIG_DFL)
            os.kill(os.getpid(), signum)
            return
        raise KeyboardInterrupt("signal %d" % signum)

    for sig in (signal.SIGTERM, signal.SIGHUP):
        if signal.getsignal(sig) != signal.SIG_IGN:
            signal.signal(sig, interrupt)


def _reason(text):
    """A confirmation records the user's decision; an empty or copied placeholder reason is not one."""
    if not text.strip() or "DECISION_AND_EVIDENCE" in text.upper():
        raise argparse.ArgumentTypeError("give the user's decision and its evidence in your own words")
    return text


def cmd_demo(a):
    from . import demo

    demo.run(a.out)


def main(argv=None):
    from . import __version__

    p = argparse.ArgumentParser(
        prog="lant",
        usage="lant [--version] COMMAND [options]",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--version", action="version", version="lant %s" % __version__)
    sub = p.add_subparsers(dest="cmd", metavar="COMMAND", help=argparse.SUPPRESS)
    config = dict(dest="config", metavar="CFG", help="the configuration file of the project")

    def command(name, fn, text):
        q = sub.add_parser(name, prog="lant " + name, help=text, description=text)
        q.set_defaults(fn=fn)
        return q

    q = command("init", cmd_init, "write a configuration file")
    q.add_argument("-o", "--output", required=True, metavar="CFG", help="the file to write")
    for f, h in (("--raw", "raw data directory"), ("--stations", "station table"), ("--out", "output directory")):
        q.add_argument(f, help=h)
    q.add_argument("--force", action="store_true", help="overwrite an existing file")
    q.add_argument("--full", action="store_true", help="write all keys with their defaults and sources")
    command("check", cmd_check, "check the installation: Python packages, GMT, DSurfTomo")
    q = command("demo", cmd_demo, "run a small synthetic array with a known answer through the whole chain")
    q.add_argument(
        "--out", default="lant-demo", metavar="DIR", help="a new directory for the demo (default: %(default)s)"
    )
    for name, fn, text in (
        ("validate", cmd_validate, "check the configuration, the station table and the raw directory"),
        ("status", cmd_status, "the state of every step, why, and the next command"),
        ("plan", cmd_plan, "what every step will do, with which values and why"),
    ):
        q = command(name, fn, text)
        q.add_argument("-c", "--config", **config)
        if name == "status":
            q.add_argument("--json", action="store_true", help="machine-readable: steps, reasons, next command")
    q = command("run", cmd_run, "run steps of the chain (all: from the first step that does not stand)")
    q.add_argument(
        "steps",
        nargs="+",
        metavar="STEP",
        help="all, or step names: clean cc stack reference disp inversion run evaluate figures",
    )
    q.add_argument("-c", "--config", **config)
    for f, h in (("--raw", "raw data directory"), ("--stations", "station table"), ("--out", "output directory")):
        q.add_argument(f, help=h + " (instead of -c: the configuration is OUT/config.yaml)")
    q.add_argument("--from", dest="from_step", metavar="STEP", help="run all: run again from this step on")
    q.add_argument("--to", dest="to_step", metavar="STEP", help="run all: stop after this step")
    q.add_argument(
        "--what",
        default="all",
        choices=("main", "checks", "all"),
        help="step run: which DSurfTomo projects (default: all)",
    )
    q.add_argument("-v", "--verbose", action="store_true", help="show detailed step output")
    q.add_argument("--debug", action="store_true", help="show the Python traceback of a failure")
    q.add_argument("--yes", action="store_true", help="acknowledge the displayed rebuild impact")
    q.add_argument("--rebuild", action="store_true", help="keep prior step outputs in unique .lant/backups directories")
    q = command("confirm", cmd_confirm, "record your decision at a judgement point")
    q.add_argument("step", choices=("clean", "stack", "reference", "inversion"), help="the step that waits for you")
    q.add_argument("-c", "--config", **config)
    q.add_argument(
        "--reason",
        required=True,
        type=_reason,
        help="your decision and its evidence (terminal use is not identity verification)",
    )
    q.add_argument("--exclude-candidates", action="store_true", help="stack: remove the candidate station-days/-hours")
    q = command("inspect", cmd_inspect, "key numbers of a step, as JSON")
    q.add_argument("topic", help="data, or a step name")
    q.add_argument("-c", "--config", **config)
    q = command("plot", cmd_plot, "draw the figures again")
    q.add_argument("names", nargs="*", metavar="NAME", help="figures to draw (default: all)")
    q.add_argument("-c", "--config", **config)
    q = command("install-external", cmd_install_external, "build fixed third-party programs in a writable location")
    q.add_argument(
        "programs",
        nargs="*",
        type=_external_program,
        metavar="PROGRAM",
        help="dsurftomo and/or pdsurftomo (default: both)",
    )
    q = command("note", cmd_note, "append a note to OUT/notes.md (also while a step runs)")
    q.add_argument("text", help="the note, in quotes")
    q.add_argument("-c", "--config", **config)
    q = command("clean-work", cmd_clean_work, "preview, or remove, the rebuildable caches of a finished project")
    q.add_argument("-c", "--config", **config)
    q.add_argument("--cache", action="append", choices=("archive", "cc"), help="which caches (default: both)")
    q.add_argument("--execute", action="store_true", help="remove them (default: preview only)")
    q.add_argument("--yes", action="store_true", help="acknowledge the displayed impact")
    a = p.parse_args(argv)
    if not a.cmd:
        p.print_help()
        print("\nFirst time: lant check, then lant demo.")
        return
    _interrupt_on_termination()
    try:
        a.fn(a)
    except KeyboardInterrupt:
        sys.exit("interrupted; state and logs are retained; inspect lant status before continuing")
    except (RuntimeError, ValueError, OSError) as exc:
        sys.exit(str(exc))
