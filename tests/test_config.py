"""Configuration: defaults, errors, integer keys, the parameter documentation and the template."""

import pytest
from util import REQ, tmp_dir, write_cfg

from lant import config


@pytest.mark.parametrize("host", [True, 1, "Gemini", "unknown", ["gemini"]])
def test_invalid_assistant_host_is_reported(tmp_path, host):
    cfg = config.defaults()
    cfg["project"]["ai_assistant"] = host
    with pytest.raises(config.ConfigError, match="project.ai_assistant"):
        config.validate(cfg)


def test_config_defaults_and_errors(tmp_path):
    """Defaults by array type; removed, unknown and inconsistent keys are refused with a message."""
    p = write_cfg(
        tmp_path,
        """
        project: {output_dir: out, array_type: surface}
        data: {raw_dir: raw, station_file: s.csv, %s}
        """
        % REQ,
    )
    cfg = config.load(p)
    assert cfg["correlation"]["fs"] == "auto" and cfg["correlation"]["seg_step_s"] == "auto"
    assert cfg["project"]["output_dir"] == str(tmp_path / "out")
    assert "notch" not in cfg and "pass2" not in cfg["reference"] and not cfg["stack"]["exclude_candidates"]
    p_old = write_cfg(
        tmp_path,
        """
        project: {output_dir: out, array_type: surface}
        data: {raw_dir: raw, station_file: s.csv, %s}
        notch: {enabled: true}
        """
        % REQ,
    )
    with pytest.raises(config.ConfigError, match="removed in 1.2"):
        config.load(p_old)
    p = write_cfg(
        tmp_path,
        """
        project: {output_dir: out, array_type: line}
        data: {raw_dir: raw, station_file: s.csv, %s}
        """
        % REQ,
    )
    cfg = config.load(p)
    assert (
        cfg["reference"]["fj_mode"] == "track_fundamental"
        and cfg["inversion"]["weight"] == config.DEFAULTS["inversion"]["weight"]
    )
    assert cfg["checks"]["controls"] == ["start_grad0.05"]
    assert cfg["inversion"]["noise"] == "auto"  # checkerboard noise from the fit of the main model
    cfg = config.load(write_cfg(tmp_path, "project: {output_dir: out}\ndata: {raw_dir: raw, station_file: s.csv}"))
    assert cfg["data"]["velocity_range_kms"] == [0.1, 5.0] and cfg["data"]["day_label"] == "utc_start"
    assert cfg["project"]["array_type"] == "surface" and not cfg["_array_type"]["given"]  # no table: areal
    for missing, msg in (
        ("project: {output_dir: out}\ndata: {raw_dir: raw}", "station_file"),
        ("project: {output_dir: out}\ndata: {raw_dir: raw, station_file: s.csv, instrument: x}", "removed in 1.3"),
        (
            "project: {output_dir: out, array_type: surface}\ndata: {raw_dir: raw, station_file: s.csv, "
            "velocity_range_kms: [2, 1]}",
            "slow < fast",
        ),
        (
            "project: {output_dir: out, array_type: surface}\ndata: {raw_dir: raw, station_file: s.csv}\n"
            "inversion: {noise: -0.01}",
            "inversion.noise",
        ),
    ):
        with pytest.raises(config.ConfigError, match=msg):
            config.load(write_cfg(tmp_path, missing))
    p = write_cfg(
        tmp_path,
        """
        project: {output_dir: out, array_type: surface}
        data: {raw_dir: raw, station_file: s.csv, chanel: "*Z", %s}
        """
        % REQ,
    )
    with pytest.raises(config.ConfigError, match="unknown configuration key 'data.chanel'"):
        config.load(p)
    p = write_cfg(
        tmp_path,
        """
        project: {output_dir: out, array_type: surface}
        data: {raw_dir: raw, station_file: s.csv, %s}
        correlation: {seg_len_s: 30, maxlag_s: 20}
        """
        % REQ,
    )
    with pytest.raises(config.ConfigError, match="seg_len_s"):
        config.load(p)


def test_config_integers():
    """Fractional hours or thread counts are refused when the configuration is read, not deep inside a step."""
    p = write_cfg(
        tmp_dir(),
        "project: {output_dir: o, array_type: surface}\ndata: {raw_dir: r, station_file: s, %s}\n"
        "stack: {block_hours: 1.5}\n" % REQ,
    )
    with pytest.raises(config.ConfigError, match="block_hours"):
        config.load(p)
    p = write_cfg(
        tmp_dir(),
        "project: {output_dir: o, array_type: surface}\ndata: {raw_dir: r, station_file: s, %s}\n"
        "dispersion: {lam_min: 2.5}\n" % REQ,
    )
    with pytest.raises(config.ConfigError, match="anchor_max_lambda"):
        config.load(p)


def test_config_rename():
    """data.rename keeps the dotted names as written."""
    p = write_cfg(
        tmp_dir(),
        "project: {output_dir: o, array_type: surface}\ndata: {raw_dir: r, station_file: s, "
        "rename: {'XX.A': 'YY.A'}, %s}\n" % REQ,
    )
    assert config.load(p)["data"]["rename"] == {"XX.A": "YY.A"}


def test_params_documented_and_template(tmp_path):
    """Every configuration key has a description and a basis; `lant init` writes a file that loads once the
    required keys are filled in, and lists every key."""
    import yaml

    from lant import params

    keys = [k for k, _ in params.flat(config.DEFAULTS)]
    assert sorted(keys) == sorted(params.DOC) and all(v[0] in params.KINDS and v[1] for v in params.DOC.values())
    t = params.template(full=True)
    assert all(("%s:" % k.split(".")[-1]) in t for k in keys)
    (tmp_path / "c.yaml").write_text(t)
    config.load(str(tmp_path / "c.yaml"))
    assert set(yaml.safe_load(t)) == {"project", "data"} and [".".join(k) for k in config.REQUIRED] == [
        "project.output_dir",
        "data.raw_dir",
        "data.station_file",
    ]
    # array type from the geometry: a 5 x 5 grid is areal, 12 stations on a line are a line (line defaults apply)
    for name, pts, kind in (
        ("g", [(i * 0.01, j * 0.01) for i in range(5) for j in range(5)], "surface"),
        ("l", [(i * 0.01, 0.0002 * (i % 2)) for i in range(12)], "line"),
    ):
        (tmp_path / (name + ".csv")).write_text(
            "station,lon,lat\n" + "".join("S%02d,%.5f,%.5f\n" % (k, 100 + x, 30 + y) for k, (x, y) in enumerate(pts))
        )
        c = config.load(
            write_cfg(tmp_path, "project: {output_dir: o}\ndata: {raw_dir: r, station_file: %s.csv}" % name)
        )
        assert c["project"]["array_type"] == kind and c["_array_type"]["clear"]
        assert (c["reference"]["fj_mode"] == "track_fundamental") == (kind == "line")


def test_compact_and_full_templates_have_identical_effective_defaults(tmp_path):
    from lant import params

    compact = tmp_path / "compact.yaml"
    full = tmp_path / "full.yaml"
    compact.write_text(params.template())
    full.write_text(params.template(full=True))
    a, b = config.load(str(compact)), config.load(str(full))
    assert {k: v for k, v in a.items() if not k.startswith("_")} == {
        k: v for k, v in b.items() if not k.startswith("_")
    }
    assert len(compact.read_text().splitlines()) < 25
    assert a["plots"]["formats"] == ["png"]


@pytest.mark.parametrize("formats", ["[pdf]", "[png, pdf]", "[]", "png"])
def test_workflow_figures_are_png_only(tmp_path, formats):
    text = "project: {output_dir: out}\ndata: {raw_dir: raw, station_file: s.csv}\nplots: {formats: %s}" % formats
    with pytest.raises(config.ConfigError, match="plots.formats must be"):
        config.load(write_cfg(tmp_path, text))


@pytest.mark.parametrize(
    "text, replaced",
    [
        ("clean: {nproc: 4}", "resources.cores"),
        ("correlation: {threads_per_process: 2}", "resources.cores"),
        ("inversion: {njobs: 2}", "resources.cores"),
        ("inversion: {solver: dsurftomo}", "Python port of DSurfTomo"),
        ("inversion: {start_model: {forward: fwd1d}}", "always DSurfTomo's surfdisp96"),
    ],
)
def test_removed_keys_say_what_replaced_them(tmp_path, text, replaced):
    p = write_cfg(
        tmp_path, "project: {output_dir: out}\ndata: {raw_dir: raw, station_file: sta.csv, %s}\n%s" % (REQ, text)
    )
    with pytest.raises(config.ConfigError, match=replaced):
        config.load(p)


def test_resources_cores_auto_or_whole_number(tmp_path):
    base = "project: {output_dir: out}\ndata: {raw_dir: raw, station_file: sta.csv, %s}\n" % REQ
    assert config.load(write_cfg(tmp_path, base))["resources"]["cores"] == "auto"
    assert config.load(write_cfg(tmp_path, base + "resources: {cores: 3}"))["resources"]["cores"] == 3
    with pytest.raises(config.ConfigError, match="resources.cores"):
        config.load(write_cfg(tmp_path, base + "resources: {cores: 0}"))


@pytest.mark.parametrize("layout", [None, 1])
def test_project_of_an_earlier_version_is_not_continued(tmp_path, layout):
    import json

    out = tmp_path / "out" / ".lant"
    out.mkdir(parents=True)
    if layout is None:
        (out / "manifest.json").write_text("{}")
    else:
        (out / "layout.json").write_text(json.dumps({"version": layout}))
    p = write_cfg(tmp_path, "project: {output_dir: out}\ndata: {raw_dir: raw, station_file: sta.csv, %s}\n" % REQ)
    with pytest.raises(config.ConfigError, match="earlier version of lant .* run the chain again"):
        config.load(p)


def test_controls_are_checked_before_anything_runs():
    cfg = config.defaults()
    cfg["data"].update(raw_dir="raw", station_file="stations.csv")
    cfg["project"]["output_dir"] = "out"
    cfg["checks"]["controls"] = ["reg_x0.5", "reg_w2", "start_grad-0.05", "start_x1.1", "reg_x1e-1"]
    config.validate(cfg)
    for bad in (
        ["reg_x0"],
        ["reg_xabc"],
        ["reg_w0"],
        ["reg_w2_d0.5"],
        ["start_grad0.6"],
        ["start_x-1"],
        ["other"],
        ["reg_x2", "reg_x2"],
    ):
        cfg["checks"]["controls"] = bad
        with pytest.raises(config.ConfigError, match="checks.controls"):
            config.validate(cfg)


@pytest.mark.parametrize("key", ["weight"])
def test_automatic_regularization_was_removed(tmp_path, key):
    """inversion.weight: auto is refused with a message (a fixed calibrated default)."""
    p = write_cfg(
        tmp_path,
        """
        project: {output_dir: out, array_type: surface}
        data: {raw_dir: raw, station_file: s.csv, %s}
        inversion: {%s: auto}
        """
        % (REQ, key),
    )
    with pytest.raises(config.ConfigError, match="auto was removed"):
        config.load(p)


@pytest.mark.parametrize("key", ["grid_refine", "sublayers", "sparsity", "threshold", "damp"])
def test_removed_inversion_keys_are_refused(tmp_path, key):
    """Keys that lant no longer needs (2.0.1) are refused with what to do instead."""
    p = write_cfg(
        tmp_path,
        """
        project: {output_dir: out, array_type: surface}
        data: {raw_dir: raw, station_file: s.csv, %s}
        inversion: {%s: 3}
        """
        % (REQ, key),
    )
    with pytest.raises(config.ConfigError, match="2.0.1"):
        config.load(p)


def test_band_edge_was_removed(tmp_path):
    """dispersion.qc.band_edge (no literature basis) is refused with a message since 2.0.1."""
    p = write_cfg(
        tmp_path,
        """
        project: {output_dir: out, array_type: surface}
        data: {raw_dir: raw, station_file: s.csv, %s}
        dispersion: {qc: {band_edge: 0.5}}
        """
        % REQ,
    )
    with pytest.raises(config.ConfigError, match="removed in 2.0.1"):
        config.load(p)


@pytest.mark.parametrize("out", ["raw", "raw/output"])
def test_output_inside_raw_is_refused(tmp_path, out):
    """The raw directory is searched recursively: an output directory in it would feed lant its own products."""
    p = write_cfg(
        tmp_path,
        """
        project: {output_dir: %s, array_type: surface}
        data: {raw_dir: raw, station_file: s.csv, %s}
        """
        % (out, REQ),
    )
    with pytest.raises(config.ConfigError, match="must not be data.raw_dir or inside it"):
        config.load(p)


@pytest.mark.parametrize(
    "change, msg",
    [
        (("data", "start_date", "2024-02-31"), "not a date"),
        (("data", "min_trace_length_s", -1.0), "min_trace_length_s"),
        (("reference", "zone_km", -2.0), "zone_km"),
        (("inversion", "depth_sequence_km", []), "at least two depths"),
    ],
)
def test_values_that_cannot_work_are_refused(change, msg):
    cfg = config.defaults()
    *keys, value = change
    config._set(cfg, keys, value)
    with pytest.raises(config.ConfigError, match=msg):
        config.validate(cfg)


def test_dates_in_order_and_steps_of_whole_samples():
    cfg = config.defaults()
    cfg["data"].update(start_date="2024-03-02", end_date="2024-03-01")
    with pytest.raises(config.ConfigError, match="start_date must not be after"):
        config.validate(cfg)
    cfg = config.defaults()
    cfg["correlation"].update(fs=50.0, seg_step_s=0.001)  # round(0.05) = 0 samples
    with pytest.raises(config.ConfigError, match="less than one sample"):
        config.validate(cfg)
