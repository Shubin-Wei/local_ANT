"""Freeze package source bytes at process import and reject mixed-source long runs.

This detects changes to an editable checkout; it cannot make a mutable environment
immutable. Production runs should still use a fixed wheel or dedicated checkout.
"""

import datetime
import hashlib
import json
import os
from pathlib import Path

PACKAGE = Path(__file__).resolve().parent


def source_files():
    return {
        str(p.relative_to(PACKAGE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(PACKAGE.rglob("*.py"))
    }


def source_digest(files):
    return hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()


SOURCES = source_files()
SOURCE_SHA256 = source_digest(SOURCES)
STARTED = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
expected = os.environ.get("LANT_EXPECTED_SOURCE_SHA256")
if expected and expected != SOURCE_SHA256:
    raise RuntimeError("worker source differs from the launching lant process; restart from a fixed installation")


def changed_sources():
    current = source_files()
    return sorted(name for name in SOURCES.keys() | current.keys() if SOURCES.get(name) != current.get(name))


def assert_unchanged():
    changed = changed_sources()
    if changed:
        raise RuntimeError(
            "lant source changed after process startup (%s); stop development in this checkout and restart "
            "the command from a fixed wheel or separate checkout" % ", ".join(changed[:8])
        )


def identity():
    return dict(source_sha256=SOURCE_SHA256, sources=SOURCES, root=str(PACKAGE), started=STARTED, pid=os.getpid())


def worker_environment():
    assert_unchanged()
    return {"LANT_EXPECTED_SOURCE_SHA256": SOURCE_SHA256}


# The modules whose code computes each processing step's products (lant/<name>.py).  Their code identity is part of
# the step record (lant.manifest.numerical_dependencies) and of the caches the steps keep (cleaned station-days, cc day
# files, the stack of every hour), so that a change of the computation makes them out of date.
STEP_SOURCES = {
    "clean": ("clean", "clean_archive", "data", "pulse", "response", "record_index", "station_match", "stations"),
    "cc": ("cc", "cc_kernel", "data", "pulse", "signal"),
    "stack": ("stack", "station_qc", "cc", "signal"),
    "reference": ("reference", "fj", "disp", "signal"),
    "disp": ("disp", "qc", "export", "signal"),
}
_CODE = {}


def _code_only(source):
    """The syntax tree without docstrings, comments, formatting or line numbers."""
    import ast

    tree = ast.parse(source)
    for node in ast.walk(tree):
        body = getattr(node, "body", None)
        if (
            isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
            and body
            and isinstance(body[0], ast.Expr)
            and isinstance(body[0].value, ast.Constant)
            and isinstance(body[0].value.value, str)
        ):
            node.body = body[1:] or [ast.Pass()]
    return ast.dump(tree, include_attributes=False)


def code_digest(step):
    """sha256 (16 hex digits) of the code of the modules of STEP_SOURCES[step]: changed by any change of the code,
    not by docstrings, comments or formatting."""
    if step not in _CODE:
        h = hashlib.sha256()
        for name in STEP_SOURCES[step]:
            h.update(name.encode())
            h.update(_code_only((PACKAGE / (name + ".py")).read_text()).encode())
        _CODE[step] = h.hexdigest()[:16]
    return _CODE[step]
