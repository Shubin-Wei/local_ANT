"""Child processes: solver process groups reaped when parallel work is interrupted, and results of
multiprocessing pools that fail instead of waiting forever when a worker is killed."""

import contextlib
import os
import signal
import subprocess
import threading
from concurrent.futures import ThreadPoolExecutor

_active = set()
_guard = threading.Lock()
_cancelled = threading.Event()


def _terminate(process):
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()
    except ProcessLookupError:
        process.wait()


def run(args, on_poll=None, **kwargs):
    """The small subset of subprocess.run needed by scientific solver calls."""
    if _cancelled.is_set():
        raise RuntimeError("solver batch interrupted")
    process = subprocess.Popen(args, start_new_session=True, **kwargs)
    with _guard:
        _active.add(process)
    try:
        if _cancelled.is_set():
            _terminate(process)
            raise RuntimeError("solver batch interrupted")
        if on_poll:
            if kwargs.get("stdout") == subprocess.PIPE or kwargs.get("stderr") == subprocess.PIPE:
                raise ValueError("polling requires file-backed output to avoid pipe deadlocks")
            while True:
                try:
                    process.wait(timeout=1)
                    break
                except subprocess.TimeoutExpired:
                    on_poll()
        stdout, stderr = process.communicate()
        return subprocess.CompletedProcess(args, process.returncode, stdout, stderr)
    except BaseException:
        _terminate(process)
        raise
    finally:
        with _guard:
            _active.discard(process)


@contextlib.contextmanager
def workers(max_workers):
    executor = ThreadPoolExecutor(max_workers=max_workers)
    try:
        yield executor
    except BaseException:
        _cancelled.set()
        with _guard:
            active = list(_active)
        for process in active:
            _terminate(process)
        raise
    finally:
        executor.shutdown(wait=True)
        _cancelled.clear()


def _run_chunk(fn, chunk):
    return [fn(task) for task in chunk]


def pool_results(pool, fn, tasks, chunksize=1, ordered=False, poll=5.0):
    """pool.imap_unordered (or pool.imap when ordered) whose iteration fails when a worker is killed.

    When the kernel kills a worker (typically for lack of memory), multiprocessing.Pool starts a replacement but the
    task the worker held never returns, so iterating over the results waits forever.  Workers are replaced only
    then (the pools here have no maxtasksperchild), so a changed set of worker processes means a lost task.  Pool
    keeps its workers in pool._pool (CPython 3.8 to 3.13).  Tasks are sent in chunks of `chunksize` as Pool itself
    does (in order within a chunk), because Pool's own chunked iterator cannot wait with a timeout."""
    import functools
    import multiprocessing

    workers = {w.pid for w in pool._pool}
    tasks = list(tasks)
    chunks = [tasks[i : i + chunksize] for i in range(0, len(tasks), chunksize)]
    results = (pool.imap if ordered else pool.imap_unordered)(functools.partial(_run_chunk, fn), chunks)
    while True:
        try:
            chunk = results.next(timeout=poll)
        except StopIteration:
            return
        except multiprocessing.TimeoutError:
            if {w.pid for w in pool._pool} != workers:
                raise RuntimeError(
                    "a worker process ended unexpectedly, most often killed for lack of memory (see `dmesg`); "
                    "its task is lost, so the step stops instead of waiting. Free memory or use fewer processes "
                    "(the nproc of this step), then run again"
                ) from None
            continue
        yield from chunk
