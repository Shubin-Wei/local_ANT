"""An interrupted parallel batch reaps its solver before releasing its caller."""

import sys
import threading
import time

import pytest

from lant import processes


def test_interrupt_terminates_running_solver(tmp_path):
    started = threading.Event()
    result = []

    def child():
        started.set()
        result.append(processes.run([sys.executable, "-c", "import time; time.sleep(60)"]))

    begin = time.monotonic()
    with pytest.raises(KeyboardInterrupt):
        with processes.workers(1) as executor:
            future = executor.submit(child)
            assert started.wait(2)
            # Wait for process registration, with a bounded wait rather than relying on a sleep race.
            for _ in range(200):
                with processes._guard:
                    if processes._active:
                        break
                time.sleep(0.01)
            else:
                pytest.fail("solver did not start")
            raise KeyboardInterrupt
    assert time.monotonic() - begin < 10
    assert future.done() and not processes._active
    assert (result and result[0].returncode != 0) or isinstance(future.exception(), RuntimeError)
    assert processes.run([sys.executable, "-c", "pass"]).returncode == 0
