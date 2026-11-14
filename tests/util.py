"""Helpers of the tests: configuration files, temporary directories, a small set of raw records."""

import textwrap

import numpy as np

REQ = "velocity_range_kms: [0.5, 3.0]"


def write_cfg(tmp_path, text):
    import pathlib

    tmp_path = pathlib.Path(tmp_path)
    p = tmp_path / "cfg.yaml"
    p.write_text(textwrap.dedent(text))
    return str(p)


def tmp_dir():
    import tempfile

    return tempfile.mkdtemp()


def make_raw(tmp_path):
    import obspy

    rng = np.random.default_rng(0)
    t0 = obspy.UTCDateTime(2024, 1, 1)
    raw = tmp_path / "raw"
    for i, name in enumerate(["A01", "A02", "A03"]):
        tr = obspy.Trace(
            rng.standard_normal(int(3 * 3600 * 100)).astype(np.float32),
            header=dict(
                network="",
                station=name,
                channel="HHZ",
                sampling_rate=100.0,
                starttime=t0 + (0.004 if i == 2 else 0.0003),
            ),
        )
        d = raw / ("deep" if i else "flat") / "x"
        d.mkdir(parents=True, exist_ok=True)
        if i == 0:
            tr.stats.network = "ZZ"
            tr.write(str(d / "rec1.mseed"), format="MSEED", encoding="FLOAT32")
            tr.copy().write(str(d / "dup.mseed"), format="MSEED", encoding="FLOAT32")  # exact duplicate
        else:
            tr.write(str(d / ("%s.sac" % name)), format="SAC")
    (raw / "notes.txt").write_text("not data")
    (tmp_path / "sta.csv").write_text(
        "network,station,latitude,longitude,elevation,instrument\n"
        "ZZ,A01,30.0,100.0,10,N1\nZZ,A02,30.01,100.0,20,N1\nZZ,A03,30.0,100.01,,N1\n"
    )
    return raw


def solver_project(folder, layout="area", velocity=None, iterations=2):
    """A tiny solver project (DSurfTomo.in, MOD, surfdata.dat) of 6 stations 0.01 deg apart: pairs measured at
    0.2, 0.5 and 1 s with a phase velocity `velocity` (3 values; default 2.4, 2.6, 2.9 km/s) x (1 + 0.01 sin k)."""
    import copy
    from itertools import combinations
    from pathlib import Path

    from lant import config
    from lant.inversion import project

    folder = Path(folder)
    z = np.array([0, 0.1, 0.3, 1, 3])
    vs = np.array([2, 2.2, 2.6, 3.1, 3.5])
    periods = np.array([0.2, 0.5, 1.0])
    velocity = np.array([2.4, 2.6, 2.9]) if velocity is None else np.asarray(velocity, float)
    points = (
        [(30 + 0.01 * i, 100 + 0.01 * j) for i in range(2) for j in range(3)]
        if layout == "area"
        else [(30, 100 + 0.01 * j) for j in range(6)]
    )
    curves = [
        dict(name="pair%d" % k, lat1=a[0], lon1=a[1], lat2=b[0], lon2=b[1], T=periods,
             v=velocity * (1 + 0.01 * np.sin(k)))
        for k, (a, b) in enumerate(combinations(points, 2))
    ]  # fmt: skip
    folder.mkdir(parents=True, exist_ok=True)
    recs, ts = project.records_from_curves(curves, 60)
    info = project.write_surfdata(recs, folder / "surfdata.dat")
    grid = project.grid(points, (0.01, 0.01), strip_cells=2 if layout == "line" else None)
    inv = copy.deepcopy(config.DEFAULTS["inversion"])
    need = project.sparsity_needed(recs, grid, len(z))
    inv.update(iterations=iterations, vmin=1.2, vmax=5.0, weight=1.0, damp=1.0,
               sparsity=float(min(1.0, max(0.05, round(2 * need + 0.005, 2)))))  # fmt: skip
    project.write_in(folder / "DSurfTomo.in", grid, len(z), info["nsrc"], inv, ts)
    project.write_mod(folder / "MOD", z, vs, grid)
    return folder
