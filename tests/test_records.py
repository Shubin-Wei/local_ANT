"""Records: the manifest and its gates, confirmations, derived parameters, the lock."""

import os

import pytest
from util import REQ, make_raw, write_cfg

from lant import config


def test_manifest_gate_confirm(tmp_path):
    """One manifest decides what stands: a step with nothing to judge does not stop the chain; a judgement stops it
    until confirmed; changed parameters, changed raw files or a rerun upstream make results out of date."""
    from lant import auto, manifest, pipeline, report
    from lant import confirm as CF

    make_raw(tmp_path)
    p = write_cfg(
        tmp_path,
        """
        project: {output_dir: out, array_type: surface}
        data: {raw_dir: raw, station_file: sta.csv, %s}
        clean: {nproc: 2}
        correlation: {fs: auto}
        """
        % REQ,
    )
    cfg = config.load(p)
    assert auto.user_value(cfg, ["correlation", "fs"]) is None  # "auto" written out = the key left out
    with pytest.raises(manifest.NotCurrent, match="lant run all"):
        pipeline.run_step(cfg, "cc")  # clean has not run
    pipeline.run_step(cfg, "clean")
    rc = report.resolved(cfg)
    st = manifest.status(rc)
    assert st["clean"]["state"] == "done" and st["cc"]["state"] == "not run"  # names matched: nothing to judge
    manifest.gate(rc, "cc")
    txt = open(tmp_path / "out" / ".lant" / "steps" / "clean.md").read()
    assert "WHAT THIS STEP DID" in txt and "1 interpolated" in txt and "derived: " in txt and "[rule]" in txt
    assert "nothing; the chain goes on" in txt
    with pytest.raises(RuntimeError, match="nothing to confirm"):
        CF.confirm(rc, "clean")
    # a judgement: the gate stops until it is confirmed, for exactly these results
    r = manifest.load(cfg)["steps"]["clean"]
    manifest.record(
        cfg,
        "clean",
        "done",
        params=r["params"],
        inputs=r["inputs"],
        to_confirm=["is A03 right?"],
        basis=CF.basis(rc, "clean", ["is A03 right?"]),
    )
    assert manifest.status(rc)["clean"]["state"] == "to confirm"
    with pytest.raises(CF.NeedsConfirmation, match="is A03 right"):
        manifest.gate(rc, "cc")
    entry, _ = CF.confirm(rc, "clean")
    assert entry["_judged"] == ["is A03 right?"]
    rc = report.resolved(config.load(p))
    assert manifest.status(rc)["clean"]["state"] == "done"
    text = open(p).read()
    assert text.count(CF.BEGIN) == 1 and "clean: {nproc: 2}" in text  # the user's part is untouched
    # computing resources do not change results: the step and its confirmation stand
    rc["clean"]["nproc"] = 7
    assert manifest.status(rc)["clean"]["state"] == "done"
    # a changed parameter, then a changed raw file
    rc["data"]["min_trace_length_s"] = 5.0
    assert "min_trace_length_s" in manifest.status(rc)["clean"]["reason"]
    rc = report.resolved(config.load(p))
    f = tmp_path / "raw" / "flat" / "x" / "rec1.mseed"
    os.utime(f, (1e9, 1e9))
    e = manifest.status(rc)["clean"]
    assert e["state"] == "out of date" and "raw files" in e["reason"]
    cmd, why = manifest.next_action(rc, manifest.status(rc))
    assert cmd.startswith("lant run all -c") and "out of date" in why


def test_noise_derived_by_run_keeps_inversion_current():
    """inversion.noise auto is derived by the run step; the inversion step records the setting, so deriving the value
    does not make the inversion step (and every later one) out of date."""
    from lant import auto, report

    cfg = config.defaults("surface")
    cfg["_user"] = {}
    assert auto.noise_is_auto(cfg) and report.used_params(cfg, "inversion")["inversion.noise"] == "auto"
    cfg["inversion"]["noise"] = 0.0011  # as auto.apply resolves it after the run step
    assert report.used_params(cfg, "inversion")["inversion.noise"] == "auto"
    assert report.used_params(cfg, "run")["inversion.noise"] == 0.0011
    cfg["_user"] = {"inversion": {"noise": 0.03}}  # a value written by the user: read by both steps
    cfg["inversion"]["noise"] = 0.03
    assert not auto.noise_is_auto(cfg) and report.used_params(cfg, "inversion")["inversion.noise"] == 0.03


def test_lock_is_exclusive(tmp_path):
    """.lant/lock: a second process cannot take it while the first holds it; nested steps of one process can."""
    import subprocess
    import sys

    from lant.provenance import Lock

    with Lock(str(tmp_path)):
        with Lock(str(tmp_path)):
            pass
        r = subprocess.run(
            [sys.executable, "-c", "from lant.provenance import Lock\nwith Lock(%r): pass" % str(tmp_path)],
            capture_output=True,
            text=True,
            cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        )
        assert r.returncode != 0 and "another lant process" in r.stderr
    r = subprocess.run(
        [sys.executable, "-c", "from lant.provenance import Lock\nwith Lock(%r): pass" % str(tmp_path)],
        capture_output=True,
        text=True,
        cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    )
    assert r.returncode == 0


def test_missing_cache_is_stale_until_explicitly_marked_cleaned(tmp_path):
    from lant import manifest, report

    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path)
    clean = tmp_path / "work/01_clean"
    (clean / "archive").mkdir(parents=True)
    for name in ("archive_index.json", "stations.csv", "archive/day.mseed"):
        (clean / name).write_text("test")
    manifest.record(cfg, "clean", "done", params=report.used_params(cfg, "clean"))
    (clean / "archive/day.mseed").unlink()
    assert manifest.status(cfg, check_inputs=False)["clean"]["state"] == "out of date"
    m = manifest.load(cfg)
    m["steps"]["clean"]["cache_cleaned"] = True
    manifest._save(cfg, m)
    state = manifest.status(cfg, check_inputs=False)["clean"]
    assert state["state"] == "done" and state["cache_cleaned"]


def test_status_reports_active_writer_without_trusting_stale_pid(tmp_path):
    from lant import inspect, records
    from lant.provenance import Lock

    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path)
    rec = records.rec_dir(cfg, create=True)
    with Lock(rec):
        state = inspect.status(cfg)
        assert state["writer"]["pid"] == os.getpid()
        assert state["next"]["command"] is None
        with open(os.path.join(rec, "lock"), "w") as fh:
            fh.write("{")
        assert inspect.status(cfg)["writer"] == {"state": "starting"}
    assert inspect.status(cfg)["writer"] is None


def test_step_progress_names_its_place_in_the_chain(tmp_path, capsys):
    from lant import records

    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path)
    cfg["_position"] = "6/9"
    with records.Step(cfg, "run"):
        pass
    assert "[lant] [6/9] run " in capsys.readouterr().err
