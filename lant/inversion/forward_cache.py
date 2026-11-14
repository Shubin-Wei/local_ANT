"""Exact, checksum-verified forward evaluations, shared by precheck/noise/residuals.

Cache only measured-data evaluations; synthetic prechecks always run (cheap, and their
noise belongs to the run). Never infer scientific identity from mtime. The number of
threads is not part of the key: the solver's results do not depend on it.
The project's existing write lock owns publication; readers see one atomic NPZ.
"""

import hashlib
import json
import os
import platform
import sys
import tempfile
import time
import zipfile
from pathlib import Path

import numpy as np

from .. import records, runtime_identity


def file_digest(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def array_digest(array):
    array = np.ascontiguousarray(array)
    return dict(shape=list(array.shape), dtype=array.dtype.str, sha256=hashlib.sha256(array.tobytes()).hexdigest())


def input_identity(project, model="start", synthetic=False):
    """Also guard model parsing against input edits between hashing and evaluation."""
    names = ["DSurfTomo.in", "MOD", "surfdata.dat"] + (["MOD.true"] if synthetic else [])
    result = {name: file_digest(Path(project) / name) for name in names}
    if model != "start":
        result["evaluation_model"] = file_digest(model)
    return result


def solver_identity():
    from .solver import identity

    return identity()


def basis(cfg, project, depths, model):
    runtime_identity.assert_unchanged()
    prefixes = ("OMP_", "OPENBLAS_", "MKL_", "NUMBA_", "NUMPY_", "GFORTRAN_", "LD_", "DYLD_", "VECLIB_", "BLIS_")
    environment = {
        k: v
        for k, v in os.environ.items()
        if k.startswith(prefixes) or k in ("LANG", "LC_ALL", "LC_NUMERIC", "PYTHONHASHSEED")
    }
    return dict(
        schema=2,
        inputs=input_identity(project),
        depths=array_digest(depths),
        model=array_digest(model),
        implementation=solver_identity(),
        runtime_source=runtime_identity.SOURCE_SHA256,
        python=sys.version,
        numpy=np.__version__,
        host=platform.node(),
        machine=platform.machine(),
        environment_sha256=hashlib.sha256(json.dumps(environment, sort_keys=True).encode()).hexdigest(),
    )


def _read(path, expected, validate):
    loaded = np.load(path, allow_pickle=False)
    if not isinstance(loaded, np.lib.npyio.NpzFile):
        raise ValueError("forward cache is not an NPZ entry")
    with loaded as entry:
        metadata = json.loads(str(entry["metadata"].item()))
        result = entry["residuals"]
    if not isinstance(metadata, dict):
        raise ValueError("forward cache metadata is not a mapping")
    if metadata.get("basis") != expected or metadata.get("result") != array_digest(result):
        raise ValueError("forward cache failed its dependency or result checksum")
    validate(result)
    return result


def _write(path, identity, result):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".tmp", delete=False) as output:
            temporary = output.name
            np.savez_compressed(
                output,
                residuals=result,
                metadata=json.dumps(dict(basis=identity, result=array_digest(result)), sort_keys=True),
            )
        os.replace(temporary, path)
    finally:
        if temporary and os.path.exists(temporary):
            os.remove(temporary)


def prediction(cfg, project, depths, model, compute, validate, unchanged_inputs, synthetic=False):
    """Reuse only identical successful requests; failed/interrupted work never becomes a cache entry."""
    started = time.perf_counter()
    state, key, compute_seconds, identity_seconds = "miss", None, 0.0, 0.0
    try:
        runtime_identity.assert_unchanged()
        unchanged_inputs()
        if synthetic:
            state, identity, path = "bypass_synthetic", None, None
        else:
            identity_start = time.perf_counter()
            identity = basis(cfg, project, depths, model)
            identity_seconds = time.perf_counter() - identity_start
            key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
            path = Path(records.rec_dir(cfg)) / "forward_cache" / (key + ".npz")
        result = None
        if path is not None and path.exists():
            try:
                result = _read(path, identity, validate)
                state = "hit"
            except (OSError, ValueError, KeyError, RuntimeError, EOFError, zipfile.BadZipFile) as exc:
                # A rebuildable cache cannot stand in for a corrupt scientific product.
                print("  forward cache rejected for %s: %s; recomputing" % (project, exc), flush=True)
                state = "rejected"
        if result is None:
            compute_start = time.perf_counter()
            try:
                result = compute()
            finally:
                compute_seconds = time.perf_counter() - compute_start
            validate(result)
        unchanged_inputs()
        runtime_identity.assert_unchanged()
        # Check executable/runtime identity again before publishing or returning a hit.
        # This is independent of file timestamps, including a same-size edited solver.
        if identity is not None:
            identity_start = time.perf_counter()
            current = basis(cfg, project, depths, model)
            identity_seconds += time.perf_counter() - identity_start
            if current != identity:
                raise RuntimeError("forward dependencies changed during evaluation; no cached result published")
        if path is not None and state != "hit":
            _write(path, identity, result)
        return result
    except BaseException as exc:
        records.log(
            cfg,
            event="forward_evaluation",
            project=str(project),
            cache=state,
            key=key,
            status="failed",
            seconds=time.perf_counter() - started,
            compute_seconds=compute_seconds,
            identity_seconds=identity_seconds,
            error="%s: %s" % (type(exc).__name__, exc),
        )
        raise
    finally:
        if sys.exc_info()[0] is None:
            records.log(
                cfg,
                event="forward_evaluation",
                project=str(project),
                cache=state,
                key=key,
                status="done",
                seconds=time.perf_counter() - started,
                compute_seconds=compute_seconds,
                identity_seconds=identity_seconds,
            )
