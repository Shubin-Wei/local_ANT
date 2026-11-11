"""Locations of the external programs built by tools/install_external.sh (inside the repository, git-ignored)."""
import os
import shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def external_path(*parts):
    return os.path.join(os.environ.get("LANT_EXTERNAL", os.path.join(ROOT, "external")), *parts)


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
    raise RuntimeError("%s not found in %s or on PATH: run `bash tools/install_external.sh %s`" % (exe, external_path("bin"), name))
