"""Locations of the external programs built by lant install-external (user cache or source checkout)."""

import os
import shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def external_path(*parts):
    # A source checkout keeps its existing development installation. Wheels use a writable user cache.
    checkout = os.path.join(ROOT, "external")
    default = (
        checkout
        if os.path.isfile(os.path.join(ROOT, "pyproject.toml"))
        else os.path.join(os.environ.get("XDG_CACHE_HOME", os.path.expanduser("~/.cache")), "lant", "external")
    )
    return os.path.join(os.environ.get("LANT_EXTERNAL", default), *parts)


# Periods a project may carry: the NP of the solver's CalSurfG (DSurfTomo src/CalSurfG.f90 NP=60; pDSurfTomo
# src_pDSurfTomo/CalSurfG_parallel.f90 NP=80).  surfdisp96 itself (and fortran/fwd1d.f90) takes 80.
SOLVER_MAX_PERIODS = {"dsurftomo": 60, "pdsurftomo": 80}


def max_periods(solver):
    return SOLVER_MAX_PERIODS[solver]


def solver_binary(name):
    """dsurftomo -> external/bin/DSurfTomo, pdsurftomo -> external/bin/pDSurfTomo (or the same name on PATH)."""
    exe = {"dsurftomo": "DSurfTomo", "pdsurftomo": "pDSurfTomo"}[name]
    p = external_path("bin", exe)
    if os.path.exists(p):
        return p
    q = shutil.which(exe)
    if q:
        return q
    raise RuntimeError(
        "%s not found in %s or on PATH: run `lant install-external %s`" % (exe, external_path("bin"), name)
    )
