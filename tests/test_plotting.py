"""Choices the figures make from the data: band edges and bands of the moveout figures, their distance and velocity
window, the example pairs, the relief tiles; and the checkerboard noise from the fit of the main model."""

import os

import h5py
import numpy as np
import pytest

from lant import config
from lant.plotting import ccf, dispersion, relief


def test_start_model_draws_parameter_nodes_without_inventing_interfaces(tmp_path, monkeypatch):
    from lant.plotting import inversion, style

    model = np.array([[0, 1.2], [0.1, 1.7], [0.4, 2.8]])
    np.savetxt(tmp_path / "vs1d.txt", model)
    np.savetxt(tmp_path / "vs1d_fit.txt", [[0.1, 1.2, 1.21], [0.5, 2.0, 1.99]])
    monkeypatch.setattr(style, "save", lambda figure, *args: figure)
    fig = inversion.start_model(config.defaults(), str(tmp_path), "unused")
    try:
        model_line = fig.axes[1].lines[0]
        assert model_line.get_drawstyle() == "default"
        assert np.array_equal(model_line.get_xdata(), model[:, 1])
        assert np.array_equal(model_line.get_ydata(), model[:, 0])
        assert fig.axes[1].yaxis_inverted()
        assert "not a detected interface" in fig.texts[0].get_text()
    finally:
        style.plt.close(fig)


def test_residual_histogram_keeps_outliers_and_separates_mean_std_and_rms(tmp_path, monkeypatch):
    from lant.plotting import inversion, style

    # Include a tail beyond the previous 99.5-percentile window, with a nonzero mean.
    observed = np.ones(201)
    residual = np.full(201, 0.01)
    residual[-1] = 0.8
    csv = tmp_path / "residuals_final.csv"
    np.savetxt(
        csv,
        np.c_[np.ones(201), observed, observed - residual, np.r_[np.ones(200), 0]],
        delimiter=",",
        header="period_s,t_obs_s,t_pred_s,kept_by_rule",
        comments="",
    )
    monkeypatch.setattr(style, "save", lambda figure, *args: figure)
    fig = inversion.residuals(config.defaults(), str(csv), {"residuals": {"iterations": []}}, "unused")
    try:
        assert "(ms)" in fig.axes[0].get_xlabel()
        histogram = fig.axes[0].patches[0].get_xy()
        assert histogram[:, 0].max() >= 800  # full tail is visible
        table = next(iter(fig.axes[-1].tables))
        assert table[(1, 4)].get_text().get_text() == "200 / 201"
        assert table[(1, 1)].get_text().get_text() == "%.1f" % np.mean(1000 * residual)
        assert table[(1, 2)].get_text().get_text() == "%.1f" % np.std(1000 * residual)
        assert table[(1, 3)].get_text().get_text() == "%.1f" % np.sqrt(np.mean((1000 * residual) ** 2))
    finally:
        style.plt.close(fig)


def test_band_edges_and_bands():
    """Band edges at whole Hz from 2 Hz, 0.5 Hz from 0.5 Hz, 0.1 Hz below; the bands join, are at least 2/3 octave
    wide, at most four, and cover the rounded measurement band."""
    assert [ccf.nice_hz(f) for f in (13.6, 2.4, 1.26, 0.63, 0.31, 0.04)] == [14.0, 2.0, 1.5, 0.5, 0.3, 0.1]
    for lo, hi in ((0.63, 14.3), (1.0, 1.3), (2.2, 47.0), (0.12, 3.1)):
        b = ccf.figure_bands(lo, hi)
        assert b[0][0] == ccf.nice_hz(lo) and b[-1][1] == ccf.nice_hz(hi) and 1 <= len(b) <= 4
        assert all(b[k][1] == b[k + 1][0] for k in range(len(b) - 1))
        assert len(b) == 1 or all(np.log2(h / l) >= 2 / 3 - 1e-9 for l, h in b)
    assert ccf.figure_bands(0.63, 14.3) == [(0.5, 1.0), (1.0, 2.0), (2.0, 5.0), (5.0, 14.0)]


def _disp(path, n=40):
    rng = np.random.default_rng(1)
    freqs = np.geomspace(1.0, 10.0, 30)
    dist = np.sort(rng.uniform(0.5, 10.0, n))
    ok = np.zeros((n, 30), bool)
    for j, d in enumerate(dist):
        ok[j, (freqs * d / 2.0 >= 1.0) & (np.arange(30) < 10 + j % 15)] = True  # at least one wavelength
    vc = np.tile(2.0 - 0.05 * np.log(freqs), (n, 1)).astype(np.float32)
    with h5py.File(path, "w") as h:
        h["freqs"], h["dist_km"], h["ok_c"], h["vc"], h["vg"] = freqs, dist, ok, vc, (0.8 * vc).astype(np.float32)
        h["pairs"] = np.array([("P%02d_Q%02d" % (j, j)).encode() for j in range(n)])
    return dist, ok, freqs


def test_moveout_window_and_example_pairs(tmp_path):
    """The moveout window of a band from the measurements in it (distance, slow and fast velocity on 0.1 km/s); the
    example pairs at the 20 / 50 / 80 % distance quantiles, ordered by distance, with rounded bands."""
    p = str(tmp_path / "disp.h5")
    dist, ok, freqs = _disp(p)
    w = ccf.band_window(p, (2.0, 5.0))
    m = (freqs >= 2.0) & (freqs <= 5.0)
    assert abs(w["dmax"] - np.percentile(dist[ok[:, m].any(1)], 99)) < 1e-9
    assert w["vslow"] <= 0.8 * 1.95 and w["vfast"] >= 2.0 - 0.05 * np.log(2.0) - 1e-6
    assert ccf.band_window(p, (20.0, 30.0)) is None  # nothing measured in the band
    pairs = dispersion.choose_pairs(p)
    d = {("P%02d_Q%02d" % (j, j)): dist[j] for j in range(len(dist))}
    assert len(pairs) == 3 and [d[n] for n, _ in pairs] == sorted(d[n] for n, _ in pairs)
    assert all(ccf.nice_hz(lo) == lo and ccf.nice_hz(hi) == hi and hi > lo for _, (lo, hi) in pairs)


def test_pair_figure_separates_reference_and_measured_wavelength_gates(tmp_path, monkeypatch):
    from lant.plotting import style

    path = tmp_path / "disp.h5"
    f = np.geomspace(0.5, 10, 20)
    ref = np.full(len(f), 2.0)
    accepted = f >= 4
    with h5py.File(path, "w") as h:
        h["pairs"], h["freqs"], h["c_ref"] = [b"A_B"], f, ref
    image = dict(freqs=f, vel=np.linspace(1, 3, 11), okc=accepted, vc=ref, vg=ref * 0.8, dist=1.0,
                 branches={}, cref=ref, gimg=np.ones((len(f), 11)), tt=np.linspace(0, 2, 20), sym=np.zeros(20))  # fmt: skip
    monkeypatch.setattr(dispersion, "_images", lambda *args: image)
    monkeypatch.setattr(style, "save", lambda fig, *args: fig)
    cfg = config.defaults()  # initial: d/lambda_ref >= 1; final: d/lambda_measured >= 2
    fig = dispersion.pair_example(cfg, str(path), "unused", "A_B", "unused", vwin=[3, 1])
    try:
        ax = fig.axes[1]
        limit = next(line for line in ax.lines if line.get_label().startswith("Final limit:"))
        assert np.allclose(limit.get_ydata(), limit.get_xdata() / 2)
        shade = next(c for c in ax.collections if c.get_label().startswith("Below"))
        assert "= 1$" in shade.get_label()
        assert any(line.get_label() == "Rejected candidate" for line in ax.lines)
    finally:
        style.plt.close(fig)


def test_relief_tiles_and_limits(tmp_path):
    """SRTM tile names on both sides of the equator and of Greenwich; no automatic relief for wide or polar maps,
    none when switched off, a missing grid file named."""
    assert relief.tile_names(99.5, 100.5, 29.5, 30.2) == ["N29E099", "N29E100", "N30E099", "N30E100"]
    assert relief.tile_names(-0.5, 0.5, -0.5, 0.5) == ["S01W001", "S01E000", "N00W001", "N00E000"]
    assert relief.tile_names(100.0, 101.0, 30.0, 31.0) == ["N30E100"]  # edges on whole degrees
    cfg = config.defaults()
    cfg["plots"]["topography"] = "auto"
    assert relief.grid(cfg, 100, 104, 30, 31, str(tmp_path))[0] is None
    assert "60 N" in relief.grid(cfg, 10, 11, 65, 66, str(tmp_path))[1]
    cfg["plots"]["topography"] = None
    assert relief.grid(cfg, 100, 101, 30, 31, str(tmp_path)) == (None, "plots.topography is none")
    cfg["plots"]["topography"] = str(tmp_path / "none.nc")
    assert "not found" in relief.grid(cfg, 100, 101, 30, 31, str(tmp_path))[1]


def test_checkerboard_noise_from_the_main_fit(tmp_path, monkeypatch):
    """inversion.noise auto: per period, 1.4826 MAD of the residuals of the main model in seconds; computed again
    through the dependency-checked forward evaluator (never by model mtime)."""
    from types import SimpleNamespace

    from lant import pipeline
    from lant.inversion import dsurf_input

    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path)
    main = tmp_path / "work" / "06_inversion" / "main"
    main.mkdir(parents=True)
    (main / "DSurfTomo.inMeasure.dat").write_text("model")
    rng = np.random.default_rng(0)
    res = np.r_[rng.normal(0, 0.01, 40), rng.normal(0, 0.04, 40)]  # the longer period twice as noisy and more
    pred = np.full(80, 2.0)
    R = np.c_[np.zeros(80), pred, pred + res, np.zeros(80), np.zeros(80), np.ones(80)]
    periods = np.r_[np.full(40, 1), np.full(40, 2)]
    monkeypatch.setattr(dsurf_input, "read_control", lambda p: SimpleNamespace(periods_rc=[0.5, 1.0], datafile="d"))
    monkeypatch.setattr(dsurf_input, "read_measurements", lambda p, c: SimpleNamespace(row_periods=lambda: periods))
    calls = []
    monkeypatch.setattr(pipeline.RUN, "evaluate", lambda cfg, d: calls.append(d) or (None, R, None))
    nz = pipeline.checkerboard_noise(cfg, str(main))
    (t1, s1), (t2, s2) = nz["value"]
    assert (t1, t2) == (0.5, 1.0) and nz["n"] == 80
    for s, r in ((s1, res[:40]), (s2, res[40:])):
        assert abs(s - 1.4826 * np.median(np.abs(r - np.median(r)))) < 1e-6
    config.validate(dict(cfg, inversion=dict(cfg["inversion"], noise=nz["value"])))  # a valid inversion.noise
    assert pipeline.checkerboard_noise(cfg, str(main))["value"] == nz["value"] and len(calls) == 2
    os.utime(main / "DSurfTomo.inMeasure.dat", (1e9, 1e9))
    pipeline.checkerboard_noise(cfg, str(main))
    assert len(calls) == 3  # mtime alone no longer decides whether the forward result stands


def test_failed_figure_refresh_keeps_previous_complete_image(tmp_path):
    from lant.plotting.figures import _Figures

    figure = _Figures.__new__(_Figures)
    figure.out, figure.done, figure.skipped, figure.failed = str(tmp_path), {}, {}, set()
    old = tmp_path / "model.png"
    old.write_bytes(b"previous complete image")

    def fail():
        with open(figure.path("model.png"), "wb") as fh:
            fh.write(b"partial")
        raise RuntimeError("renderer failed")

    figure.attempt("model", fail)
    assert old.read_bytes() == b"previous complete image" and "model" in figure.skipped

    def success():
        path = figure.path("model.png")
        with open(path, "wb") as fh:
            fh.write(b"new complete image")
        return [path]

    figure.attempt("replacement", success)
    assert old.read_bytes() == b"new complete image"
    assert figure.done["replacement"] == [str(old)]
    assert sorted(p.name for p in tmp_path.iterdir()) == ["model.png"]


def test_parallel_failure_keeps_previous_group_and_other_completed_outputs(tmp_path, monkeypatch):
    from lant import config
    from lant.plotting import figure_batch, figures

    figure = figures._Figures.__new__(figures._Figures)
    figure.out, figure.done, figure.skipped, figure.failed = str(tmp_path), {}, {}, set()
    figure.todo, figure.seconds, figure.model = ["ccf_moveout", "stations"], {}, "unused"
    figure.stations = {}
    old = tmp_path / "2_ccf_moveout_1-2Hz.png"
    old.write_bytes(b"previous complete moveout")
    complete = tmp_path / "1_stations.png"
    complete.write_bytes(b"completed worker output")
    worker_results = [
        dict(
            drawn={},
            not_drawn={"2_ccf_moveout_1-2Hz": "renderer failed"},
            failed=["2_ccf_moveout_1-2Hz"],
            figure_seconds={},
            relief_notes=[],
        ),
        dict(
            drawn={"stations": [str(complete)]},
            not_drawn={},
            failed=[],
            figure_seconds={"stations": 1.0},
            relief_notes=[],
        ),
    ]
    monkeypatch.setattr(figures, "_Figures", lambda cfg, which: figure)
    monkeypatch.setattr(figures, "resolve", lambda *args: {})
    monkeypatch.setattr(figure_batch, "render", lambda cfg, which: worker_results)
    result = figures.make_all(config.defaults())
    assert old.read_bytes() == b"previous complete moveout"
    assert result["drawn"]["stations"] == [str(complete)]
    assert result["not_drawn"]["2_ccf_moveout_1-2Hz"] == "renderer failed"
    assert result["failed"] == ["2_ccf_moveout_1-2Hz: renderer failed"]  # the step is then incomplete
    assert result["figure_seconds"] == {"stations": 1.0}


def test_nothing_to_draw_is_not_a_failure_and_obsolete_variants_go(tmp_path):
    from lant.plotting.figures import _Figures

    figure = _Figures.__new__(_Figures)
    figure.out, figure.done, figure.skipped, figure.failed = str(tmp_path), {}, {}, set()
    figure.todo = ["ccf_moveout", "recovery"]
    figure.attempt("recovery", lambda: None)  # no checkerboards: nothing to draw
    assert figure.skipped["recovery"] == "nothing to draw" and not figure.failed
    stale = tmp_path / "2_ccf_moveout_7-10Hz.png"  # a band without measurements in the present data
    stale.write_bytes(b"figure of an earlier run")
    figure.skipped["2_ccf_moveout_7-10Hz"] = "no phase velocity was measured in this band"
    figure.remove_old()
    assert not stale.exists()
    kept = tmp_path / "2_ccf_moveout_1-2Hz.png"
    kept.write_bytes(b"previous complete image")
    figure.attempt("2_ccf_moveout_1-2Hz", lambda: (_ for _ in ()).throw(RuntimeError("GMT failed")))
    figure.remove_old()
    assert kept.exists()  # a failed group keeps its previous figures


def test_figure_with_numbers_keeps_them_and_its_staged_files(tmp_path):
    from lant.plotting.figures import _Figures

    figure = _Figures.__new__(_Figures)
    figure.out, figure.done, figure.skipped, figure.failed = str(tmp_path), {}, {}, set()

    def stack_time():  # like ccf.stack_time: files and key numbers
        path = figure.path("2_ccf_stack_time_1.png")
        open(path, "wb").write(b"image")
        return dict(files=[path], snr_final=12.5)

    figure.attempt("ccf_stack_time_1", stack_time)
    assert figure.done["ccf_stack_time_1"] == dict(files=[str(tmp_path / "2_ccf_stack_time_1.png")], snr_final=12.5)
    figure.todo = ["ccf_stack_time"]
    figure.remove_old()
    assert (tmp_path / "2_ccf_stack_time_1.png").exists()


def test_recovery_combines_heatmaps_and_depth_profiles_with_readable_labels(monkeypatch):
    from lant.plotting import inversion, style

    monkeypatch.setattr(style, "save", lambda figure, *args: figure)
    evaluation = {
        "checkerboards": {
            "cb_1_2.5km": [{"depth_km": 0.0, "slope": 0.32, "r": 0.61}, {"depth_km": 0.2, "slope": 0.14, "r": 0.52}],
            "cb_vertical_thick": [
                {"depth_km": 0.0, "slope": 0.21, "r": 0.48},
                {"depth_km": 0.2, "slope": -0.02, "r": -0.15},
            ],
        }
    }
    fig = inversion.recovery(config.defaults(), evaluation, "unused")
    try:
        titles = [ax.get_title(loc="left") for ax in fig.axes if ax.get_title(loc="left")]
        assert "(a) Amplitude recovery" in titles and "(b) Pattern recovery" in titles
        assert any("100% is ideal" in ax.get_ylabel() for ax in fig.axes)
        assert any("1.0 is ideal" in ax.get_ylabel() for ax in fig.axes)
        ticklabels = [t.get_text() for t in fig.axes[0].get_yticklabels()]
        assert ticklabels == ["2.5 km cells", "Thick vertical"]
        assert any("damped amplitude" in t.get_text() for t in fig.texts)
        assert any("DWS > 0" in t.get_text() for t in fig.texts)
    finally:
        style.plt.close(fig)


def test_station_map_and_line_section_layout_semantics(tmp_path, monkeypatch):
    from lant.plotting import gmt, maps, sections

    captured = []

    def fake_run(tmp, pre, cmd, out, *args, **kwargs):
        captured.append((tmp, list(pre), list(cmd)))
        return [out + ".png"]

    monkeypatch.setattr(gmt, "run", fake_run)
    cfg = config.defaults()
    cfg["plots"]["topography"] = None
    stations = {
        "S1": {"lon": 100.01, "lat": 30.00, "elev": 950.0},
        "S2": {"lon": 100.01, "lat": 30.07, "elev": 980.0},
        "S3": {"lon": 100.01, "lat": 30.14, "elev": 1020.0},
    }
    maps.station_map(cfg, stations, str(tmp_path / "1_stations"), inversion_count=2)
    tmp, _, cmd = captured[-1]
    legend = open(tmp + "/legend.txt").read()
    assert "Recorded (3)" in legend and "Inversion (2)" in legend
    assert any("-DjT" in c for c in cmd)

    proj = tmp_path / "work" / "06_inversion" / "main"
    proj.mkdir(parents=True)
    np.savetxt(proj.parent / "vs1d.txt", [[0.0, 1.5], [0.2, 2.0], [0.5, 2.8]])
    import json

    (proj / "project.json").write_text(json.dumps({"depths": [0.0, 0.2, 0.5]}))
    rows = [
        [lon, lat, dep, 1.5 + dep + 0.05 * (lat - 30.00)]
        for dep in (0.0, 0.2)
        for lat in (29.99, 30.07, 30.15)
        for lon in (100.00, 100.02)
    ]
    np.savetxt(proj / "DSurfTomo.inMeasure.dat", rows)
    cfg["plots"]["vs_range"] = [1.2, 2.5]
    cfg["plots"]["contour_step"] = 0.2
    sections.line_section(cfg, str(proj), stations, str(tmp_path / "5_sections"), label="line_1")
    tmp_sec, pre_sec, cmd_sec = captured[-1]
    head0 = open(tmp_sec + "/r0_head.txt").read()
    head1 = open(tmp_sec + "/r1_head.txt").read()
    note = open(tmp_sec + "/note.txt").read()
    assert "line_1 | Vs model" in head0 and "vertical exaggeration" not in head0
    assert "Anomaly relative to 1-D reference" in head1
    assert "V.E. =" in note and "resolution is not masked" in note
    assert len(pre_sec) == 2


def test_recovery_single_depth_and_full_negative_correlation(monkeypatch):
    from lant.plotting import inversion, style

    monkeypatch.setattr(style, "save", lambda figure, *args: figure)
    evaluation = {"checkerboards": {"one": [{"depth_km": 0.1, "slope": 1.7, "r": -0.8}]}}
    fig = inversion.recovery(config.defaults(), evaluation, "unused")
    try:
        assert fig.axes[0].collections[0].get_array().size == 1
        assert fig.axes[0].collections[0].norm.vmax >= 170
        assert fig.axes[2].get_xlim()[1] > 1.7
        assert fig.axes[3].get_xlim()[0] < -0.8
    finally:
        style.plt.close(fig)
    evaluation["checkerboards"]["cb_vertical"] = [{"depth_km": 0.2, "slope": 0.2, "r": 0.3}]
    fig = inversion.recovery(config.defaults(), evaluation, "unused")
    try:  # one horizontal and one depth-varying test on different depths: blank cells, no error
        values = fig.axes[0].collections[0].get_array()
        assert values.size == 4 and int(np.ma.count_masked(values)) == 2
        assert [t.get_text() for t in fig.axes[0].get_yticklabels()] == ["one", "Depth-varying"]
    finally:
        style.plt.close(fig)
    evaluation["checkerboards"]["bad"] = [{"depth_km": 0.2, "slope": 0.2, "r": 0.3}, {"depth_km": 0.1}]
    with pytest.raises(ValueError, match="strictly increasing"):
        inversion.recovery(config.defaults(), evaluation, "unused")


def test_checker_labels():
    from lant.plotting.inversion import _checker_label

    assert _checker_label("cb_1.2km") == "1.2 km cells"
    assert _checker_label("cb_2_3.7km") == "3.7 km cells"
    assert _checker_label("cb_vertical") == "Depth-varying"


def test_checkerboard_figures_follow_the_present_checks(tmp_path):
    """A checkerboard directory left from an earlier configuration is not drawn (nor evaluated)."""
    import json
    import types

    from lant import paths
    from lant.plotting import figures

    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path)
    root = tmp_path / "work" / paths.INV
    for name in ("cb_1km", "cb_old"):
        (root / "checkerboards" / name).mkdir(parents=True)
        (root / "checkerboards" / name / "DSurfTomo.inSyn.dat").write_text("")
    (root / "checks.json").write_text(json.dumps(dict(projects=[["checkerboards", "cb_1km"]])))
    drawn = []
    F = types.SimpleNamespace(
        cfg=cfg, todo=["checkerboards"], line=True, stations={}, path=str, attempt=lambda n, fn: drawn.append(n)
    )
    figures._checkerboards(F)
    assert drawn == ["checkerboard_cb_1km"]
