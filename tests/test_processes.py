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


def test_file_backed_process_can_report_progress(tmp_path):
    polled = []
    with (tmp_path / "run.log").open("w") as output:
        result = processes.run(
            [sys.executable, "-c", "import time; time.sleep(1.2)"], stdout=output, on_poll=lambda: polled.append(True)
        )
    assert result.returncode == 0 and polled


def test_termination_signal_stops_the_solver_group(tmp_path):
    import os
    import signal
    import subprocess

    marker = "sleep %s" % (60 + os.getpid() % 1000 / 1000)  # unique, to find exactly this process
    child = (
        "import sys\n"
        "from lant import cli, processes\n"
        "cli._interrupt_on_termination()\n"
        "try:\n"
        "    with processes.workers(1) as ex:\n"
        "        ex.submit(processes.run, %r.split()).result()\n"
        "except KeyboardInterrupt:\n"
        "    sys.exit(130)\n" % marker
    )
    lant_process = subprocess.Popen(
        [sys.executable, "-c", child],
        cwd=str(tmp_path.parents[0]),
        env=dict(os.environ, PYTHONPATH=os.path.dirname(os.path.dirname(__file__))),
    )

    def solver_running():
        out = subprocess.run(["ps", "-eo", "args"], capture_output=True, text=True).stdout
        return any(line.strip() == marker for line in out.splitlines())

    for _ in range(100):
        if solver_running():
            break
        time.sleep(0.05)
    else:
        lant_process.kill()
        pytest.fail("solver did not start")
    lant_process.send_signal(signal.SIGTERM)
    assert lant_process.wait(timeout=10) == 130  # the interrupt path ran, not the default kill
    assert not solver_running()  # no orphaned solver keeps writing
