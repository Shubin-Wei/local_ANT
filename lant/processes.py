"""Reap solver process groups when parallel work is interrupted."""

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


def run(args, **kwargs):
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
