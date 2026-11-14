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
