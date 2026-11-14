"""Bounded, spawned main-figure workers; Matplotlib and GMT state stay private.

Each task owns one named figure group and publishes through _Figures.attempt.
Only the parent removes obsolete variants, after merging successful and failed
tasks. No worker changes measurements, models, quality decisions or plot ranges.
"""

import os
from multiprocessing import get_context

from .. import paths, resources, runtime_identity
from ..processes import pool_results
from ..progress import track


def memory_estimate(cfg):
    """Scheduling estimate, not an RSS bound: empirical 600 MB interpreter/rendering budget plus work arrays.

    Allow eight float64 buffers for the complete stack and one pair's hourly
    traces, and reuse the existing dispersion work estimate for pair images.
    The resource planner retains its available-memory headroom.
    """
    import h5py

    estimate = 600e6
    stack = paths.stack_product(cfg)
    if os.path.exists(stack):
        with h5py.File(stack, "r") as h:
            pairs, samples = h["stack"].shape
            blocks = h["daily"].shape[0]
        estimate += 64 * samples * (pairs + 24 * blocks)
        disp = paths.P(cfg, paths.DISP, "disp.h5")
        if os.path.exists(disp):
            from ..disp import memory_estimate as pair_memory

            with h5py.File(disp, "r") as h:
                nf, nv = len(h["freqs"]), len(h["vel_axis"])
            estimate = max(estimate, pair_memory(stack, 1, nf, nv, 0) + 200e6)
    model = paths.P(cfg, paths.INV, "main", "DSurfTomo.inMeasure.dat")
    if os.path.exists(model):
        estimate += 16 * os.path.getsize(model)  # parsed rows, grids and interpolation/contour work
    return estimate


def _draw(task):
    cfg, name, expected_source = task
    if runtime_identity.SOURCE_SHA256 != expected_source:
        raise RuntimeError("figure worker source differs from its parent; restart from a fixed checkout")
    from . import gmt as G
    from .figures import GROUPS, _Figures

    G.RELIEF_NOTES.clear()
    figure = _Figures(cfg, [name])
    for group in GROUPS:
        group(figure)
    runtime_identity.assert_unchanged()
    return dict(
        drawn=figure.done,
        not_drawn=figure.skipped,
        failed=sorted(figure.failed),
        figure_seconds=figure.seconds,
        relief_notes=list(G.RELIEF_NOTES),
    )


def render(cfg, which):
    """Ordered worker results, or None to use the existing serial path when one process suffices."""
    which = list(dict.fromkeys(which))
    plan = resources.plan(cfg, len(which), memory_estimate(cfg), "main figure groups")
    plan["threads_each"] = 1
    if plan["parallel"] == 1:
        return None
    runtime_identity.assert_unchanged()
    tasks = [(cfg, name, runtime_identity.SOURCE_SHA256) for name in which]
    # Spawn avoids forking a parent that already used Numba, BLAS, HDF5 or Matplotlib.
    with get_context("spawn").Pool(plan["parallel"], initializer=resources.limit_threads, initargs=(1,)) as pool:
        results = list(track(pool_results(pool, _draw, tasks, ordered=True), len(tasks), "main figure groups"))
        resources.note_workers(w.pid for w in pool._pool)
    runtime_identity.assert_unchanged()
    return results
