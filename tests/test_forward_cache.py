"""Exact reuse across lifecycle stages, with dependency and failure regressions."""

import json
import os

import numpy as np
import pytest

from lant import config, pipeline, records
from lant.inversion import forward_cache as cache
from lant.inversion import model_evaluation as run
from tests.util import solver_project


def fixture(tmp_path, monkeypatch):
    output = tmp_path / "output"
    project = solver_project(output / "area", velocity=[2.5, 2.5, 2.5])
    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(output)
    identity = {"solver": "fixed"}
    monkeypatch.setattr(cache, "solver_identity", lambda: dict(identity))
    calls = []

    def forward(cfg, project, depths, model, **kwargs):
        _, _, _, _, periods = run.read_in(os.path.join(project, "DSurfTomo.in"))
        data = run.read_data(os.path.join(project, "surfdata.dat"), periods)
        distance = np.arange(1, len(data) + 1, dtype=float)
        observed = distance / data[:, 5]
        predicted = observed + 0.01 * np.sin(distance)
        calls.append(model.copy())
        return np.c_[distance, predicted, observed, predicted, observed, np.ones(len(data))]

    monkeypatch.setattr(run, "_forward_once", forward)
    return cfg, project, calls, identity


def events(cfg):
    return [
        item
        for line in (cache.Path(records.rec_dir(cfg)) / "run_log.jsonl").read_text().splitlines()
        if (item := json.loads(line))["event"] == "forward_evaluation"
    ]


def test_precheck_noise_and_residuals_share_exact_predictions(tmp_path, monkeypatch):
    cfg, project, calls, _ = fixture(tmp_path, monkeypatch)
    run.precheck(cfg, str(project))
    run.evaluate(cfg, str(project), "start")
    assert len(calls) == 1
    # A final model and its last iteration copy are the same full-grid request.
    from lant.inversion import dsurf_input, dsurf_main

    p = dsurf_input.read_project(project)
    p.vs[:-1, 1:-1, 1:-1] += np.float32(0.01)
    rows = "".join(line + "\n" for line in dsurf_main.model_lines(p.control, p.depths, p.vs))
    for name in ("DSurfTomo.inMeasure.dat", "DSurfTomo.inMeasure.dat.iter002"):
        (project / name).write_text(rows)
    monkeypatch.setattr(pipeline, "inv_dir", lambda cfg, *args: str(project.joinpath(*args)))
    noise = pipeline.checkerboard_noise(cfg, str(project))
    # Even a legacy mtime-only noise file must not bypass checked predictions.
    (project / "checkerboard_noise.json").write_text(
        json.dumps({"model_mtime": (project / "DSurfTomo.inMeasure.dat").stat().st_mtime, "value": 999})
    )
    assert pipeline.checkerboard_noise(cfg, str(project))["value"] == noise["value"]
    run.residuals(cfg, str(project), str(project / "residuals_final"))
    assert len(calls) == 2  # initial and final: no repeated native computation
    assert [item["cache"] for item in events(cfg)] == ["miss", "hit", "miss", "hit", "hit", "hit", "hit"]
    assert sum(item["compute_seconds"] == 0 for item in events(cfg)) == 5


def test_same_mtime_edits_and_observation_row_order_invalidate(tmp_path, monkeypatch):
    cfg, project, calls, _ = fixture(tmp_path, monkeypatch)
    run.evaluate(cfg, str(project), "start")
    data = project / "surfdata.dat"
    stamp = data.stat().st_mtime_ns
    lines = data.read_text().splitlines()
    indices = [i for i, line in enumerate(lines) if not line.startswith("#")]
    # Exchange receiver rows in one source block: count and observation values stay valid.
    i, j = indices[:2]
    lines[i], lines[j] = lines[j], lines[i]
    data.write_text("\n".join(lines) + "\n")
    os.utime(data, ns=(stamp, stamp))
    run.evaluate(cfg, str(project), "start")
    assert len(calls) == 2
    model = project / "MOD"
    stamp = model.stat().st_mtime_ns
    model.write_text(model.read_text().replace("2.00000", "2.01000"))
    os.utime(model, ns=(stamp, stamp))
    run.evaluate(cfg, str(project), "start")
    assert len(calls) == 3


def test_implementation_and_input_settings_invalidate_but_threads_do_not(tmp_path, monkeypatch):
    cfg, project, calls, identity = fixture(tmp_path, monkeypatch)
    run.evaluate(cfg, str(project), "start")
    identity["solver"] = "source or dependency changed"
    run.evaluate(cfg, str(project), "start")
    run.evaluate(cfg, str(project), "start", threads=2)  # the result does not depend on the threads
    settings = project / "DSurfTomo.in"
    settings.write_text(settings.read_text() + "\n")
    run.evaluate(cfg, str(project), "start", threads=2)
    assert len(calls) == 3


@pytest.mark.parametrize("failure", [RuntimeError("solver failed"), KeyboardInterrupt()])
def test_failed_or_interrupted_forward_is_not_published(tmp_path, monkeypatch, failure):
    cfg, project, calls, _ = fixture(tmp_path, monkeypatch)
    forward = run._forward_once

    def fail(*args, **kwargs):
        raise failure

    monkeypatch.setattr(run, "_forward_once", fail)
    with pytest.raises(type(failure)):
        run.evaluate(cfg, str(project), "start")
    assert not list((cache.Path(records.rec_dir(cfg)) / "forward_cache").glob("*.npz"))
    assert events(cfg)[-1]["status"] == "failed"
    monkeypatch.setattr(run, "_forward_once", forward)
    run.evaluate(cfg, str(project), "start")
    assert len(calls) == 1


def test_corrupt_or_partial_cache_is_recomputed(tmp_path, monkeypatch):
    cfg, project, calls, _ = fixture(tmp_path, monkeypatch)
    run.evaluate(cfg, str(project), "start")
    entry = next((cache.Path(records.rec_dir(cfg)) / "forward_cache").glob("*.npz"))
    with np.load(entry, allow_pickle=False) as saved:
        result, metadata = saved["residuals"], saved["metadata"]
    result[0, 1] += 1
    np.savez(entry, residuals=result, metadata=metadata)
    run.evaluate(cfg, str(project), "start")
    assert len(calls) == 2 and events(cfg)[-1]["cache"] == "rejected"
    entry.write_bytes(b"partial")
    run.evaluate(cfg, str(project), "start")
    assert len(calls) == 3


def test_midflight_input_or_solver_changes_never_publish(tmp_path, monkeypatch):
    cfg, project, _, identity = fixture(tmp_path, monkeypatch)
    forward = run._forward_once

    def changed(*args, **kwargs):
        result = forward(*args, **kwargs)
        identity["solver"] = "replaced while running"
        return result

    monkeypatch.setattr(run, "_forward_once", changed)
    with pytest.raises(RuntimeError, match="dependencies changed"):
        run.evaluate(cfg, str(project), "start")
    assert not list((cache.Path(records.rec_dir(cfg)) / "forward_cache").glob("*.npz"))

    def changed_input(*args, **kwargs):
        result = forward(*args, **kwargs)
        data = project / "surfdata.dat"
        data.write_text(data.read_text() + "\n")
        return result

    monkeypatch.setattr(run, "_forward_once", changed_input)
    with pytest.raises(RuntimeError, match="input files changed"):
        run.evaluate(cfg, str(project), "start")
    assert not list((cache.Path(records.rec_dir(cfg)) / "forward_cache").glob("*.npz"))


def test_stochastic_synthetic_precheck_is_never_reused(tmp_path, monkeypatch):
    cfg, project, calls, _ = fixture(tmp_path, monkeypatch)
    (project / "MOD.true").write_text("truth")
    run.precheck(cfg, str(project))
    run.precheck(cfg, str(project))
    assert len(calls) == 2 and all(item["cache"] == "bypass_synthetic" for item in events(cfg))


def test_solver_identity_follows_the_source_of_the_port(monkeypatch):
    from lant.inversion import least_squares, solver

    original = cache.solver_identity()
    assert original == solver.identity() and "lant.inversion.least_squares" in original["sources"]
    read = cache.Path.read_bytes
    monkeypatch.setattr(
        cache.Path,
        "read_bytes",
        lambda p: read(p) + (b"#" if p.name == cache.Path(least_squares.__file__).name else b""),
    )
    assert cache.solver_identity() != original


def test_numerical_environment_changes_invalidate(tmp_path, monkeypatch):
    cfg, project, calls, _ = fixture(tmp_path, monkeypatch)
    run.evaluate(cfg, str(project), "start")
    monkeypatch.setenv("OPENBLAS_NUM_THREADS", "17")
    run.evaluate(cfg, str(project), "start")
    assert len(calls) == 2


def test_invalid_results_never_become_a_cache_entry(tmp_path, monkeypatch):
    cfg, project, _, _ = fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(run, "_forward_once", lambda *a, **k: np.full((45, 6), np.nan))
    with pytest.raises(RuntimeError, match="invalid travel times"):
        run.evaluate(cfg, str(project), "start")
    assert not list((cache.Path(records.rec_dir(cfg)) / "forward_cache").glob("*.npz"))


def test_interrupted_atomic_publication_removes_temporary_file(tmp_path, monkeypatch):
    path = tmp_path / "cache/result.npz"
    monkeypatch.setattr(cache.os, "replace", lambda *args: (_ for _ in ()).throw(KeyboardInterrupt()))
    with pytest.raises(KeyboardInterrupt):
        cache._write(path, {"basis": "checked"}, np.ones((4, 6)))
    assert list(path.parent.iterdir()) == []
