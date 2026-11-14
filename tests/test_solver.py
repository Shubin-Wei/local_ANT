"""lant's use of the Python DSurfTomo: the solver process, its record, DWS, and the in-process forward evaluation."""

import json
import os

import numpy as np
import pytest

pytest.importorskip("numba")
from lant.inversion import run, solver  # noqa: E402
from tests.util import solver_project  # noqa: E402


def test_solver_run_writes_record_dws_and_counts_as_done(tmp_path):
    d = solver_project(tmp_path / "p", iterations=1)
    record = solver.run(d, 1, estimate=solver.memory_estimate(str(d)))
    assert record["identity"] == solver.identity() and record["threads"] == 1
    assert record["peak_memory_gb"] > 0 and record["memory_estimate_gb"] > 0 and len(record["iterations"]) == 1
    assert json.loads((d / "solver.json").read_text()) == json.loads(json.dumps(record))
    dws = np.loadtxt(d / solver.DWS_FILE)
    model = np.loadtxt(d / "DSurfTomo.inMeasure.dat")
    assert np.array_equal(dws[:, :3], model[:, :3]) and (dws[:, 3] > 0).any() and (dws[:, 3] >= 0).all()
    assert run.status(str(d)) == "done"
    (d / "solver.json").unlink()  # killed between the model and the record: run again
    assert run.status(str(d)) == "unfinished"


def test_evaluation_matches_the_first_iteration_and_not_the_threads(tmp_path):
    d = solver_project(tmp_path / "p", iterations=1)
    nx, ny, nz = run.read_in(os.path.join(d, "DSurfTomo.in"))[1:4]
    depz, vs = run.full_grid(str(d), "start", nx, ny, nz)
    one = solver.evaluate(d, depz, vs, threads=1)
    two = solver.evaluate(d, depz, vs, threads=2)
    assert one.tobytes() == two.tobytes()
    solver.run(d, 1)
    first = np.loadtxt(d / "residualFirst.dat")
    np.testing.assert_allclose(one, first, rtol=1e-6, atol=1e-7)


def test_synthetic_evaluation_uses_the_true_model(tmp_path):
    d = solver_project(tmp_path / "p", iterations=1)
    nx, ny, nz = run.read_in(os.path.join(d, "DSurfTomo.in"))[1:4]
    depz, vs = run.full_grid(str(d), "start", nx, ny, nz)
    true = vs.copy()
    true[:-1] *= 1.05
    (d / "MOD.true").write_text("\n".join(" ".join("%.5f" % v for v in row) for row in true.reshape(-1, nx)) + "\n")
    measured, synthetic = solver.evaluate(d, depz, vs), solver.evaluate(d, depz, vs, synthetic=True)
    assert np.array_equal(measured[:, 1], synthetic[:, 1])  # same prediction
    assert np.all(synthetic[:, 2] < measured[:, 1])  # 5 % faster truth: earlier synthetic arrivals
    assert not (d / "velmap2dRc.dat").exists()


def test_memory_estimate_grows_with_the_data(tmp_path):
    small = solver.memory_estimate(str(solver_project(tmp_path / "a", layout="line")))
    large = solver.memory_estimate(str(solver_project(tmp_path / "b")))
    assert solver.BASE_BYTES < small and solver.BASE_BYTES < large


def test_identity_is_the_same_in_the_solver_process(tmp_path):
    import subprocess
    import sys

    d = solver_project(tmp_path / "p", iterations=1)
    env = dict(os.environ, PYTHONPATH=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    subprocess.run([sys.executable, "-m", "lant.inversion.solver", "--threads", "1"], cwd=d, env=env, check=True,
                   capture_output=True)  # fmt: skip
    assert json.loads((d / "solver.json").read_text())["identity"] == solver.identity()
    assert run.status(str(d)) == "done"
