"""Export of the curves and the DSurfTomo projects: one segment per curve, project files, period limit, node order,
depth nodes, forward modelling."""
import os
from pathlib import Path

import numpy as np
import pytest

from lant import config, export
from lant.inversion import project


def test_one_segment():
    """export.longest_run keeps one run per curve: the longest, on a tie the one with the higher mean SNR."""
    ok = np.array([0, 1, 1, 0, 1, 1, 1, 0], bool)
    snr = np.array([0, 50, 50, 0, 5, 5, 5, 0], float)
    assert export.longest_run(ok, snr).tolist() == [0, 0, 0, 0, 1, 1, 1, 0]
    ok = np.array([1, 1, 0, 1, 1], bool)
    snr = np.array([5, 5, 0, 9, 9], float)
    assert export.longest_run(ok, snr).tolist() == [0, 0, 0, 1, 1]


def test_project_files(tmp_path):
    """surfdata.dat from the curves, the inversion grid with a ring outside the stations, the depth signs of the
    vertical checks."""
    curves = [dict(name="A_B", lat1=30.0, lon1=100.0, lat2=30.05, lon2=100.02, T=np.array([0.5, 1.0]), v=np.array([2.5, 2.7])),
              dict(name="A_C", lat1=30.0, lon1=100.0, lat2=30.02, lon2=100.06, T=np.array([1.0]), v=np.array([2.8]))]
    recs, periods = project.records_from_curves(curves, 60)
    assert periods == [0.5, 1.0]
    info = project.write_surfdata(recs, str(tmp_path / "surfdata.dat"))
    lines = (tmp_path / "surfdata.dat").read_text().splitlines()
    assert lines[0] == "# 30.000000 100.000000 1 2 0" and info["travel_times"] == 3
    g = project.grid(info["stations"], (0.005, 0.005), refine=2)
    assert g["dvxd"] == 0.0025 and (g["nx"] - 1) % 2 == 0
    lat, lon = project.node_coords(g)
    assert lat[1] > 30.05 and lon[1] < 100.0              # stations inside the boundary ring
    z = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    s = project.depth_signs(z, [0, 0.2, 0.4, 0.6], [0.3])
    assert np.allclose(s, [1, 1, 1, 0, -1, -1, -1])


def test_period_cap():
    """DSurfTomo's CalSurfG holds 60 periods (pDSurfTomo 80): 61 periods are refused before a project is written."""
    from lant.tools import max_periods
    T = 1 / np.geomspace(0.5, 20, 61)
    curves = [dict(name="A_B", lat1=30.0, lon1=100.0, lat2=30.05, lon2=100.02, T=T, v=np.full(61, 2.5))]
    assert max_periods("dsurftomo") == 60 and max_periods("pdsurftomo") == 80
    project.records_from_curves(curves[:1], max_periods("pdsurftomo"), "pdsurftomo")
    with pytest.raises(RuntimeError, match="at most 60"):
        project.records_from_curves(curves, max_periods("dsurftomo"), "dsurftomo")
    curves[0]["T"], curves[0]["v"] = T[:60], np.full(60, 2.5)
    project.records_from_curves(curves, 60)


def test_model_node_order(tmp_path):
    """A model file with a repeated node in place of a missing one is not taken as complete."""
    from lant.inversion.run import check_node_order
    nx, ny, nz = 4, 5, 3
    rows = [(100 + 0.01 * j, 30 - 0.01 * i, z, 2.0) for z in (0.0, 0.5) for j in range(1, ny - 1) for i in range(1, nx - 1)]
    m = np.array(rows)
    check_node_order(m, nx, ny, nz)
    m[3] = m[2]
    with pytest.raises(RuntimeError):
        check_node_order(m, nx, ny, nz)


def test_depth_nodes_layers_curvature():
    """Depth nodes from the wavelengths (the last but one at the deepest sensitivity), the layer limit of DSurfTomo,
    and the 1-D curvature operator of the starting model."""
    from lant.inversion import start as St
    z = St.auto_sequence(np.c_[1 / np.array([10.0, 1.0, 0.5]), [2.0, 3.0, 3.3]])
    assert z[0] == 0 and abs(z[1] - 0.05) < 1e-9 and z[-2] <= 3.3 < z[-1]
    St.check_layers(50, 3)
    with pytest.raises(RuntimeError):
        St.check_layers(51, 3)
    z = St.depth_nodes(np.array([[0.1, 2.0], [2.0, 3.0]]), [0, 0.001, 0.002])[0]
    assert z[-2] <= 3.0 < z[-1]
    # curvature: zero for a straight line on any spacing; about the same measure for 6 or 24 steps
    for z in (np.array([0, 0.1, 0.3, 0.7, 1.5, 3.0]), np.linspace(0, 3, 25)):
        assert np.allclose(St.curvature_operator(z) @ (2 + 0.3 * z), 0)
    q = [np.sum((St.curvature_operator(z) @ (z / 3) ** 3) ** 2) for z in (np.linspace(0, 3, 7), np.linspace(0, 3, 25))]
    assert all(abs(x / 12.0 - 1) < 0.25 for x in q)                 # integral of (6 zeta)^2 over 0..1 = 12


def test_forward_period_order():
    """surfdisp96 needs increasing periods: forward() sorts them (strongly layered model, decreasing periods)."""
    from lant.inversion.start import forward
    from lant.tools import external_path
    if not os.path.exists(external_path("bin", "fwd1d")):
        pytest.skip("fwd1d not built")
    z = np.array([0.0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2.0])
    vs = np.array([0.25, 0.30, 0.40, 0.55, 0.80, 1.00, 1.40, 1.80, 2.10, 2.50])
    T = 1 / np.geomspace(0.5, 15, 12)                  # decreasing periods
    c = forward(z, vs, T)
    assert np.all(np.diff(c) < 0)                      # phase velocity decreases with frequency, no jump
    assert np.allclose(c, forward(z, vs, T[::-1])[::-1])


def test_observed_iteration_handles_missing_and_buffered_logs(tmp_path):
    from lant.inversion.run import observed_iteration
    assert observed_iteration(str(tmp_path)) is None
    log = tmp_path / "DSurfTomo.in.log"
    log.write_text("initializing\n")
    assert observed_iteration(str(tmp_path)) is None
    log.write_text("1st iteration...\nresiduals\n 2th iteration...\n")
    assert observed_iteration(str(tmp_path)) == 2


def test_solver_run_captures_output_and_reports_observed_iteration(tmp_path, monkeypatch):
    import sys
    from lant import progress
    from lant.inversion.run import run_one

    shown = []
    monkeypatch.setattr(progress, "detail", shown.append)
    (tmp_path / "DSurfTomo.in").write_text("4                                c: maximum iteration\n")
    binary = tmp_path / "fake_solver"
    binary.write_text("#!%s\n" % sys.executable +
                      "from pathlib import Path\nimport os, time\n"
                      "Path('env.txt').write_text(os.environ.get('GFORTRAN_UNBUFFERED_ALL', ''))\n"
                      "Path('DSurfTomo.in.log').write_text('1st iteration...\\n')\n"
                      "print('solver output', flush=True)\ntime.sleep(1.2)\n"
                      "Path('DSurfTomo.inMeasure.dat').write_text('100 30 0 2.5\\n')\n"
                      "with open('DSurfTomo.in.log', 'a') as f: f.write('Program finishes successfully\\n')\n")
    binary.chmod(0o755)
    _, result = run_one(str(tmp_path), str(binary), 1)
    assert result == "ok" and "solver output" in (tmp_path / "run.log").read_text()
    assert "iterations: %s 1/4" % tmp_path.name in shown and shown[-1] == ""   # cleared when the project ends
    assert (tmp_path / "env.txt").read_text() == "y"                           # the log reaches the file at once


def test_interrupted_project_is_set_aside_and_runs_again(tmp_path):
    from lant.inversion import run as RUN

    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path)
    project = tmp_path / "work/06_inversion/main"
    project.mkdir(parents=True)
    for name in RUN.BUILT:
        (project / name).write_text("built\n")
    (project / "DSurfTomo.in.log").write_text("1st iteration...\n")        # killed: no success line
    (project / "DSurfTomo.inMeasure.dat.iter001").write_text("partial\n")
    assert RUN.status(str(project)) == "unfinished"
    moved = RUN.set_aside_unfinished(cfg, [str(project)])
    assert RUN.status(str(project)) == "new"
    assert sorted(p.name for p in project.iterdir()) == sorted(RUN.BUILT)
    backup = Path(moved[str(project)])
    assert str(backup).startswith(str(tmp_path / ".lant/backups"))
    assert sorted(p.name for p in backup.iterdir()) == ["DSurfTomo.in.log", "DSurfTomo.inMeasure.dat.iter001"]
