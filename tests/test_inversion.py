"""Export of the curves and the DSurfTomo projects: one segment per curve, project files, period limit, node order,
depth nodes, forward modelling."""

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
    curves = [
        dict(
            name="A_B", lat1=30.0, lon1=100.0, lat2=30.05, lon2=100.02, T=np.array([0.5, 1.0]), v=np.array([2.5, 2.7])
        ),
        dict(name="A_C", lat1=30.0, lon1=100.0, lat2=30.02, lon2=100.06, T=np.array([1.0]), v=np.array([2.8])),
    ]
    recs, periods = project.records_from_curves(curves, 60)
    assert periods == [0.5, 1.0]
    info = project.write_surfdata(recs, str(tmp_path / "surfdata.dat"))
    lines = (tmp_path / "surfdata.dat").read_text().splitlines()
    assert lines[0] == "# 30.000000 100.000000 1 2 0" and info["travel_times"] == 3
    g = project.grid(info["stations"], (0.005, 0.005), refine=2)
    assert g["dvxd"] == 0.0025 and (g["nx"] - 1) % 2 == 0
    lat, lon = project.node_coords(g)
    assert lat[1] > 30.05 and lon[1] < 100.0  # stations inside the boundary ring
    z = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    s = project.depth_signs(z, [0, 0.2, 0.4, 0.6], [0.3])
    assert np.allclose(s, [1, 1, 1, 0, -1, -1, -1])


def test_period_cap():
    """DSurfTomo's CalSurfG holds 60 periods: more are refused before a project is written (the 1-D forward
    calculation of the starting model, surfdisp96, holds 80)."""
    from lant.inversion.start import MAX_PERIODS, PROJECT_MAX_PERIODS

    T = 1 / np.geomspace(0.5, 20, 61)
    curves = [dict(name="A_B", lat1=30.0, lon1=100.0, lat2=30.05, lon2=100.02, T=T, v=np.full(61, 2.5))]
    assert (MAX_PERIODS, PROJECT_MAX_PERIODS) == (80, 60)
    with pytest.raises(RuntimeError, match="at most 60"):
        project.records_from_curves(curves, PROJECT_MAX_PERIODS)
    curves[0]["T"], curves[0]["v"] = T[:60], np.full(60, 2.5)
    project.records_from_curves(curves, PROJECT_MAX_PERIODS)


def test_model_node_order(tmp_path):
    """A model file with a repeated node in place of a missing one is not taken as complete."""
    from lant.inversion.run import check_node_order

    nx, ny, nz = 4, 5, 3
    rows = [
        (100 + 0.01 * j, 30 - 0.01 * i, z, 2.0) for z in (0.0, 0.5) for j in range(1, ny - 1) for i in range(1, nx - 1)
    ]
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
    assert all(abs(x / 12.0 - 1) < 0.25 for x in q)  # integral of (6 zeta)^2 over 0..1 = 12


def test_forward_period_order():
    """surfdisp96 needs increasing periods: forward() sorts them (strongly layered model, decreasing periods)."""
    from lant.inversion.start import forward

    z = np.array([0.0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2.0])
    vs = np.array([0.25, 0.30, 0.40, 0.55, 0.80, 1.00, 1.40, 1.80, 2.10, 2.50])
    T = 1 / np.geomspace(0.5, 15, 12)  # decreasing periods
    c = forward(z, vs, T)
    assert np.all(np.diff(c) < 0)  # phase velocity decreases with frequency, no jump
    assert np.allclose(c, forward(z, vs, T[::-1])[::-1])


def test_forward_is_dsurftomo_caldespersion_and_fits_a_start_model():
    """The 1-D forward is DSurfTomo's own (the native depthkernel's phase velocity of the same column, bit for bit);
    the starting-model fit runs on it; a Rayleigh half-space gives the Rayleigh-polynomial speed."""
    from pathlib import Path

    from lant.inversion import start

    ref = np.load(Path(__file__).parent / "data" / "native_depthkernel.npz")
    for j, i in ((0, 0), (1, 2)):
        c = start.forward(ref["depths"], ref["vs"][:, j, i], ref["periods"], minthk=int(ref["minthk"]))
        np.testing.assert_array_equal(c, ref["phase"][:, j, i])
    t, z = np.array([0.1, 0.3, 1]), np.array([0, 0.1, 0.5, 2])
    measured = start.forward(z, [2, 2.2, 2.6, 3.3], t)
    vs, predicted, history = start.invert1d(t, measured, z, iters=2, log=lambda _: None)
    assert len(history) == 2 and np.all(vs > 0) and np.all(np.isfinite(predicted))
    v = np.float32(3.0)
    vp = 0.9409 + 2.0947 * v - 0.8206 * v**2 + 0.2683 * v**3 - 0.0251 * v**4
    ratio = (v / vp) ** 2
    roots = np.roots([1, -8, 24 - 16 * ratio, -16 * (1 - ratio)])
    analytic = v * np.sqrt([x.real for x in roots if abs(x.imag) < 1e-10 and 0 < x.real < 1][0])
    values = start.forward([0, 1], [3, 3], [0.02, 0.01])
    assert np.array_equal(values, start.forward([0, 1], [3, 3], [0.01, 0.02])[::-1])  # solved in increasing order
    assert np.max(np.abs(values - analytic)) < 2e-4  # flat half-space, short periods
    with pytest.raises(NotImplementedError):
        start.forward(z, [2, 2.2, 2.6, 3.3], t, igr=1)


def test_numerical_dependency_change_invalidates_selected_fit(tmp_path, monkeypatch):
    from lant import config, manifest, report

    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path)
    monkeypatch.setattr("lant.artifacts.required", lambda *args: [])
    identity = {"numba": "old"}
    monkeypatch.setattr(
        manifest, "numerical_dependencies", lambda cfg, step: dict(identity) if step == "inversion" else {}
    )
    manifest.record(cfg, "inversion", "done", params=report.used_params(cfg, "inversion"))
    identity["numba"] = "new"
    state = manifest.status(cfg, check_inputs=False)["inversion"]
    assert state["state"] == "out of date" and "dependency versions" in state["reason"]


def test_observed_iteration_handles_missing_and_buffered_logs(tmp_path):
    from lant.inversion.run import observed_iteration

    assert observed_iteration(str(tmp_path)) is None
    log = tmp_path / "DSurfTomo.in.log"
    log.write_text("initializing\n")
    assert observed_iteration(str(tmp_path)) is None
    log.write_text("1st iteration...\nresiduals\n 2th iteration...\n")
    assert observed_iteration(str(tmp_path)) == 2


def test_solver_run_captures_output_and_reports_observed_iteration(tmp_path, monkeypatch):
    import json
    import sys

    from lant import progress
    from lant.inversion.run import run_one

    shown = []
    monkeypatch.setattr(progress, "detail", shown.append)
    (tmp_path / "DSurfTomo.in").write_text("4                                c: maximum iteration\n")
    from lant.inversion import run, solver

    record = json.dumps({"identity": solver.identity()})
    fake = tmp_path / "fake_solver.py"
    fake.write_text(
        "from pathlib import Path\nimport os, time\n"
        "Path('env.txt').write_text(os.environ.get('NUMBA_NUM_THREADS', '') + ' ' + os.environ['OMP_NUM_THREADS'])\n"
        "Path('DSurfTomo.in.log').write_text('1st iteration...\\n')\n"
        "print('solver output', flush=True)\ntime.sleep(1.2)\n"
        "Path('DSurfTomo.inMeasure.dat').write_text('100 30 0 2.5\\n')\n"
        "Path('solver.json').write_text(%r)\n"
        "with open('DSurfTomo.in.log', 'a') as f: f.write('Program finishes successfully\\n')\n" % record
    )
    monkeypatch.setattr(run, "solver_command", lambda threads, estimate: [sys.executable, str(fake)])
    monkeypatch.setattr(solver, "memory_estimate", lambda d: 1e6)
    _, result = run_one(str(tmp_path), 3)
    assert result == "ok" and "solver output" in (tmp_path / "run.log").read_text()
    assert "iterations: %s 1/4" % tmp_path.name in shown and shown[-1] == ""  # cleared when the project ends
    assert (tmp_path / "env.txt").read_text() == "3 1"  # Numba has the threads, BLAS/OpenMP one each


def test_interrupted_project_is_set_aside_and_runs_again(tmp_path):
    from lant.inversion import run as RUN

    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path)
    project = tmp_path / "work/06_inversion/main"
    project.mkdir(parents=True)
    for name in RUN.BUILT:
        (project / name).write_text("built\n")
    (project / "DSurfTomo.in.log").write_text("1st iteration...\n")  # killed: no success line
    (project / "DSurfTomo.inMeasure.dat.iter001").write_text("partial\n")
    assert RUN.status(str(project)) == "unfinished"
    moved = RUN.set_aside_unfinished(cfg, [str(project)])
    assert RUN.status(str(project)) == "new"
    assert sorted(p.name for p in project.iterdir()) == sorted(RUN.BUILT)
    backup = Path(moved[str(project)])
    assert str(backup).startswith(str(tmp_path / ".lant/backups"))
    assert sorted(p.name for p in backup.iterdir()) == ["DSurfTomo.in.log", "DSurfTomo.inMeasure.dat.iter001"]


def _design_project(tmp_path, spacing_deg, periods, velocity):
    """A main project with stations on a 4 x 4 grid, every pair measured at `periods` with one phase velocity."""
    pts = [(30 + spacing_deg * i, 100 + spacing_deg * j) for i in range(4) for j in range(4)]
    curves = [
        dict(name="%d_%d" % (a, b), lat1=pts[a][0], lon1=pts[a][1], lat2=pts[b][0], lon2=pts[b][1], T=np.array(periods),
             v=np.full(len(periods), velocity))
        for a in range(len(pts)) for b in range(a + 1, len(pts))
    ]  # fmt: skip
    recs, T = project.records_from_curves(curves, 60)
    info = project.write_surfdata(recs, str(tmp_path / "surfdata.dat"))
    g = project.grid(info["stations"], (spacing_deg / 2, spacing_deg / 2))
    inv = dict(weight=1, damp=1, sublayers=3, vmin=1, vmax=4, iterations=1, sparsity=0.1, threshold=3)
    project.write_in(str(tmp_path / "DSurfTomo.in"), g, 6, info["nsrc"], inv, T)
    return str(tmp_path)


def test_checker_design_from_wavelength(tmp_path):
    """Three horizontal tests of 1, 2 and 3 median wavelengths and one depth-varying test of 2, reversal at the median
    wavelength / 3; notes when the stations are farther apart than the smallest cell or the aperture is short."""
    from lant.inversion.checker_design import design

    z = [0, 0.1, 0.2, 0.4, 0.6, 0.8]
    main = _design_project(tmp_path, 0.01, [0.2, 0.4, 0.6], 2.0)  # wavelengths 0.4, 0.8, 1.2 km; spacing ~1 km
    h, v, info = design(main, z, "surface")
    assert info["median_wavelength_km"] == pytest.approx(0.8)
    assert 0.9 < info["median_nearest_station_km"] < 1.2
    assert [c["cell_km"] for c in h] == [[0.8, 0.8], [1.6, 1.6], [2.4, 2.4]]
    assert [c["name"] for c in h] == ["cb_0.8km", "cb_1.6km", "cb_2.4km"]
    assert "more coarsely" in info["note"]
    assert len(v) == 1 and v[0]["reversal_depths_km"] == [0.267] and v[0]["cell_km"] == [1.6, 1.6]
    h, v, info = design(main, z, "line")
    assert all(c["cell_km"][1] == 0.0 for c in h + v)
    (tmp_path / "w").mkdir()
    main = _design_project(tmp_path / "w", 0.002, [0.5, 1.0, 1.5], 2.0)  # wavelengths 1-3 km, spacing ~0.2 km
    h, v, info = design(main, z, "surface")
    assert h[0]["name"] == "cb_2km" and h[0]["cell_km"] == [2.0, 2.0]
    assert "fewer than two cells" in info["note"] and "more coarsely" not in info["note"]
    h, v, info = design(main, [0, 0.1, 0.2, 0.3], "surface")  # reversal 0.667 km below every updated node
    assert v == [] and "no depth-varying test" in info["vertical_skipped"]


def test_checker_recovery_counts_only_nodes_with_rays(tmp_path):
    from lant.inversion.run import checker_recovery

    lon, lat = np.meshgrid([100.0, 100.01, 100.02], [30.0, 30.01], indexing="ij")
    nodes = np.c_[lon.ravel(), lat.ravel()]
    true = 2.0 + 0.2 * np.array([1, -1, 1, -1, 1, -1])
    recovered = 2.0 + 0.1 * np.array([1, -1, 1, -1, -9, 9])  # the last two nodes: no ray, wrong sign
    dws = np.array([1.0, 2.0, 0.5, 3.0, 0.0, 0.0])
    for name, v in (("Vs_model.real", true), ("DSurfTomo.inSyn.dat", recovered), ("DWS.dat", dws)):
        np.savetxt(tmp_path / name, np.c_[nodes, np.zeros(6), v])
    (tmp_path / "MOD").write_text("0.0 1.0\n" + "2.0 2.0 2.0\n" * 4 + "3.0 3.0 3.0\n" * 4)
    (row,) = checker_recovery(str(tmp_path))
    assert row["n"] == 4 and row["n_all"] == 6 and row["r"] == pytest.approx(1.0) and row["slope"] == pytest.approx(0.5)


def test_regularization_controls_carry_their_weight(tmp_path):
    """reg_x<f> scales the main project's smoothing weight, reg_w<w>_d<d> sets both; DSurfTomo.in and project.json
    of the control agree; a factor that is not a positive number is refused."""
    import json

    from tests.util import solver_project

    main = solver_project(tmp_path / "main")
    z = np.array([0, 0.1, 0.3, 1, 3])
    np.savetxt(main / "vs1d.txt", np.c_[z, [2, 2.2, 2.6, 3.1, 3.5]])
    meta = dict(grid=project.grid([(30, 100), (30.01, 100.02)], (0.01, 0.01)), weight=1.4, damp=1.0)
    (main / "project.json").write_text(json.dumps(meta))
    cfg = config.defaults()
    cfg["checks"].update(checkerboards=[], vertical=[], controls=["reg_x0.5", "reg_w2_d0.5"])
    built = project.build_checks(cfg, str(main), str(tmp_path / "checks"))
    assert [n for n, _ in built["built"]] == ["reg_x0.5", "reg_w2_d0.5"]
    for name, w, d in (("reg_x0.5", 0.7, 1.0), ("reg_w2_d0.5", 2.0, 0.5)):
        c = tmp_path / "checks" / "controls" / name
        line = next(ln for ln in (c / "DSurfTomo.in").read_text().splitlines() if "c: weight damp" in ln)
        assert [float(x) for x in line.split()[:2]] == [w, d]
        m = json.loads((c / "project.json").read_text())
        assert (m["weight"], m["damp"], m["control"]) == (w, d, name)
    for bad in ("reg_x0", "reg_xabc", "reg_w1"):
        cfg["checks"]["controls"] = [bad]
        with pytest.raises(RuntimeError, match="needs a positive number"):
            project.build_checks(cfg, str(main), str(tmp_path / "bad"))


def test_grid_spacing_is_a_third_of_the_median_wavelength():
    """inversion.grid_spacing_deg auto: the final cell is the median c T / 3 km in both directions (longitude
    widened by 1 / cos(latitude)), two significant digits; the coarse spacing is grid_refine times it."""
    import math

    from lant.inversion.project import KM_PER_DEG, grid_spacing, median_wavelength

    recs = [dict(values={0.5: 2.0, 1.0: 2.4}), dict(values={0.5: 1.8, 2.0: 3.0})]  # c T: 1.0, 2.4, 0.9, 6.0
    assert median_wavelength(recs) == pytest.approx(1.7)
    pts = [(60.0, 10.0), (60.0, 10.1)]
    (dlat, dlon), why = grid_spacing(recs, pts, 2)
    assert dlat == pytest.approx(2 * float("%.2g" % (1.7 / 3 / KM_PER_DEG)))
    assert dlon == pytest.approx(2 * float("%.2g" % (1.7 / 3 / KM_PER_DEG / math.cos(math.radians(60)))))
    assert "/ 3" in why
