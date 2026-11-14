"""Running solver projects (the Python DSurfTomo, lant.inversion.solver), and the forward evaluation of a model on a
fixed data set.

A project counts as finished only if its log ends with "Program finishes successfully", its model file is complete
(one finite row per output node), its solver.json names the current implementation and none of its inputs is newer
than the log.  The output of an interrupted run is never overwritten: it is moved to .lant/backups (reported) and
the project runs again from its inputs.  A project whose inputs are newer than its log, or that an earlier
implementation computed (stale), is an error.
"""

import json
import os
import subprocess
import sys
import threading

import numpy as np

from .. import records, resources
from ..processes import run as run_process
from ..processes import workers
from ..progress import track
from .model_evaluation import check_node_order as check_node_order
from .model_evaluation import checker_recovery as checker_recovery
from .model_evaluation import evaluate as evaluate
from .model_evaluation import full_grid as full_grid
from .model_evaluation import logged_residuals as logged_residuals
from .model_evaluation import precheck as precheck
from .model_evaluation import read_data as read_data
from .model_evaluation import read_in as read_in
from .model_evaluation import read_mod as read_mod
from .model_evaluation import read_model as read_model
from .model_evaluation import residual_stats as residual_stats
from .model_evaluation import residuals as residuals

INPUTS = ("DSurfTomo.in", "MOD", "MOD.true", "surfdata.dat")
BUILT = INPUTS + ("vs1d.txt", "project.json")


def model_complete(d):
    """The model file exists, has one row per output node ((nx - 2) x (ny - 2) x (nz - 1): DSurfTomo writes neither
    the outer ring of the grid nor the fixed bottom node) and only finite, positive velocities."""
    import json

    synthetic = os.path.exists(os.path.join(d, "MOD.true"))  # checkerboards write the recovered model here
    p = os.path.join(d, "DSurfTomo.inSyn.dat" if synthetic else "DSurfTomo.inMeasure.dat")
    meta = os.path.join(d, "project.json")
    if not os.path.exists(p) or not os.path.getsize(p):
        return False
    try:
        m = np.loadtxt(p, ndmin=2)
    except ValueError:
        return False
    if os.path.exists(meta):
        g = json.load(open(meta))
        try:  # DSurfTomo drops the rim and the fixed bottom node
            check_node_order(m, g["grid"]["nx"], g["grid"]["ny"], len(g["depths"]), p)
        except RuntimeError:
            return False
    return m.shape[1] >= 4 and bool(np.all(np.isfinite(m[:, 3]))) and bool(np.all(m[:, 3] > 0))


def model_says_done(d):
    log = os.path.join(d, "DSurfTomo.in.log")
    return os.path.exists(log) and "Program finishes successfully" in open(log, errors="ignore").read()[-2000:]


def status(d):
    """missing | new | stale (an input newer than the log) | done (the log ends with the success line and the model
    file is complete) | unfinished."""
    log = os.path.join(d, "DSurfTomo.in.log")
    if not os.path.exists(os.path.join(d, "DSurfTomo.in")):
        return "missing"
    if not os.path.exists(log) and not os.path.exists(os.path.join(d, "run.log")):
        return "new"
    for f in INPUTS:
        p = os.path.join(d, f)
        if os.path.exists(p) and os.path.exists(log) and os.path.getmtime(p) > os.path.getmtime(log):
            return "stale"
    if os.path.exists(log):
        tail = open(log, errors="ignore").read()[-2000:].strip().splitlines()[-5:]
        record = os.path.join(d, "solver.json")
        if any("Program finishes successfully" in ln for ln in tail) and model_complete(d) and os.path.exists(record):
            from .solver import identity

            try:
                if json.load(open(record)).get("identity") != identity():
                    return "stale"
            except (ValueError, OSError):
                return "stale"
            return "done"
    return "unfinished"


def observed_iteration(project):
    """Read only a bounded log tail; buffered logs may have no iteration information yet."""
    import re

    path = os.path.join(project, "DSurfTomo.in.log")
    try:
        with open(path, "rb") as fh:
            fh.seek(max(0, os.path.getsize(path) - 65536))
            text = fh.read().decode("utf-8", errors="replace")
    except FileNotFoundError:
        return None
    matches = re.findall(r"(\d+)(?:st|nd|rd|th)\s+iteration", text, re.IGNORECASE)
    return int(matches[-1]) if matches else None


def max_iterations(project):
    """The iteration count written into DSurfTomo.in ("c: maximum iteration"), or None."""
    try:
        for line in open(os.path.join(project, "DSurfTomo.in")):
            if "c: maximum iteration" in line:
                return int(line.split()[0])
    except (OSError, ValueError):
        pass
    return None


def set_aside_unfinished(cfg, dirs):
    """Move the solver output of interrupted projects to one .lant/backups directory, so that they run again from
    their inputs (nothing is deleted) -> {project: backup}."""
    import datetime

    from .. import records

    root = cfg["project"]["output_dir"]
    backup = os.path.join(
        records.rec_dir(cfg), "backups", datetime.datetime.now().strftime("%Y%m%dT%H%M%S.%f") + "-interrupted"
    )
    moved = {}
    for d in dirs:
        if status(d) != "unfinished":
            continue
        target = os.path.join(backup, os.path.relpath(d, root))
        os.makedirs(target)
        for name in sorted(os.listdir(d)):
            if name not in BUILT:
                os.rename(os.path.join(d, name), os.path.join(target, name))
        moved[d] = target
        message = "%s: the output of an interrupted run is kept in %s; running it again" % (d, target)
        print("  " + message, flush=True)
        sys.stderr.write("[lant] %s\n" % message)  # the user should see it without -v
    return moved


_running = {}  # project name -> observed iteration, for the progress line
_running_guard = threading.Lock()


def _show_running(name=None, value=None):
    from ..progress import detail

    with _running_guard:
        if value is None:
            _running.pop(name, None)
        else:
            _running[name] = value
        text = ", ".join("%s %s" % kv for kv in sorted(_running.items()))
    detail("iterations: " + text if text else "")


def solver_command(threads, estimate):
    """The solver process of a project (run in its directory)."""
    return [sys.executable, "-m", "lant.inversion.solver", "--threads", str(threads), "--estimate", "%.0f" % estimate]


def run_one(d, threads):
    """Run the solver in project d with `threads` threads -> (d, message)."""
    st = status(d)
    if st == "done":
        return d, "skipped (finished, inputs unchanged)"
    if st == "stale":
        return (
            d,
            "FAILED: its inputs or numerical implementation changed (stale); rebuild with lant run all --from inversion",
        )
    if st == "unfinished":  # lant.pipeline sets interrupted runs aside before this
        return d, "FAILED: output of an earlier run is present (unfinished); move it away or rerun the run step"
    from ..runtime_identity import assert_unchanged, worker_environment
    from .solver import memory_estimate

    # Numba gets the threads; BLAS and OpenMP inside it stay single-threaded so that they do not multiply them.
    env = dict(os.environ, **resources.thread_environment(1), **worker_environment())
    env["NUMBA_NUM_THREADS"] = str(threads)
    progress_errors = []
    name, total = os.path.basename(d), max_iterations(d)
    _show_running(name, "starting")

    def show_iteration():
        try:
            iteration = observed_iteration(d)
        except OSError as exc:
            if not progress_errors:
                print("  iteration progress unavailable: %s" % exc, flush=True)
                progress_errors.append(str(exc))
            return
        if iteration is not None:
            _show_running(name, "%d/%s" % (iteration, total or "?"))

    command = solver_command(threads, memory_estimate(d))
    try:
        with open(os.path.join(d, "run.log"), "w") as fh:
            rc = run_process(
                command, cwd=d, env=env, stdout=fh, stderr=subprocess.STDOUT, on_poll=show_iteration
            ).returncode
    finally:
        _show_running(name)
    assert_unchanged()
    if rc != 0:
        return d, "FAILED: exit status %d (see %s/run.log)" % (rc, d)
    st = status(d)
    if st != "done":
        if st == "stale":
            return d, "FAILED: solver.json names another implementation than this lant (or an input changed)"
        return d, (
            "FAILED: the model file is incomplete (rows, NaN or non-positive velocities)"
            if model_says_done(d)
            else "FAILED: no 'Program finishes successfully' at the end of the log"
        )
    return d, "ok"


@records.timed_phase("solver_batch")
def run_projects(cfg, dirs, on_done=None):
    """The projects, as many at a time as the cores and the memory allow (lant.resources.plan); each takes its
    threads when it starts: its share of the cores of the step, fewer when the machine has become busier since."""
    from .solver import memory_estimate

    todo = [d for d in dirs if status(d) != "done"]
    peak = max((memory_estimate(d) for d in todo), default=0)
    groups = sorted({"main" if os.path.basename(d) == "main" else os.path.basename(os.path.dirname(d)) for d in todo})
    plan = resources.plan(cfg, max(1, len(todo)), peak or None, "solver projects (%s)" % ", ".join(groups))
    own = resources.OwnLoad()
    print(
        "  cores %d (%s): %d project(s) at a time, up to %d threads each"
        % (plan["cores"], plan["cores_rule"], plan["parallel"], plan["threads_each"]),
        flush=True,
    )

    def timed_project(directory):
        with records.phase(cfg, "solver_project", project=str(directory)) as timing:
            now, how = resources.cores(cfg["resources"]["cores"], own.estimate())
            threads = max(1, min(plan["threads_each"], now // plan["parallel"]))
            timing.update(threads=threads, cores_rule=how)
            job = own.start(threads)
            try:
                result = run_one(directory, threads)
            finally:
                own.stop(job)
            timing["result"] = result[1]
            if result[1].startswith("FAILED"):
                timing["status"] = "failed"
            elif result[1].startswith("skipped"):
                timing["status"] = "skipped"
            else:
                record = json.load(open(os.path.join(directory, "solver.json")))
                resources.note_solver(
                    dict(project=os.path.relpath(directory, cfg["project"]["output_dir"]), threads=threads,
                         wall_seconds=record["wall_seconds"], peak_memory_gb=record["peak_memory_gb"],
                         memory_estimate_gb=record["memory_estimate_gb"])
                )  # fmt: skip
            return result

    with workers(max_workers=plan["parallel"]) as ex:
        from concurrent.futures import as_completed

        futures = [ex.submit(timed_project, d) for d in dirs]
        kinds = sorted(
            {"main" if os.path.basename(d) == "main" else os.path.basename(os.path.dirname(d)) for d in dirs}
        )
        unit = "solver projects (%s)" % ", ".join(kinds)
        done = []
        for future in track(as_completed(futures), len(futures), unit, eta=False):
            result = future.result()
            done.append(result)
            if on_done is not None:
                on_done(*result)
    res = sorted(done, key=lambda r: dirs.index(r[0]))  # summaries in project order, not in order of completion
    for d, msg in res:
        print("  %s %s" % (d, msg), flush=True)
    failed = [d for d, m in res if m.startswith("FAILED")]
    return dict(results=dict(res), failed=failed)


# ------------------------------------------------------------------------------------------------ models
