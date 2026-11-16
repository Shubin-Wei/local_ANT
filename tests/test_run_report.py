"""REPORT.md: the run report (lant.run_report)."""

from pathlib import Path

from lant import config, run_report


def test_step_times_come_from_the_latest_completed_run(tmp_path):
    """The steps table of REPORT.md shows the wall time of the latest completed run of each step."""
    import json

    from lant import records

    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path)
    log = Path(records.rec_dir(cfg, create=True)) / "run_log.jsonl"
    with open(log, "w") as fh:
        for e in (
            dict(step="cc", event="done", seconds=50.0),
            dict(step="cc", event="failed", seconds=3.0),
            dict(step="cc", event="done", seconds=136.0),
            dict(step="run", event="done", seconds=1256.2),
        ):
            fh.write(json.dumps(e) + "\n")
        fh.write("not json\n")
    assert run_report.step_seconds(cfg) == {"cc": 136.0, "run": 1256.2}
    assert [run_report._duration(s) for s in (None, 52.0, 1256.2, 9000.0)] == ["", "52 s", "21 min", "2.5 h"]


def test_report_shows_current_step_duration_before_its_done_log_is_written(tmp_path):
    from lant import manifest

    cfg = config.defaults()
    cfg["project"]["output_dir"] = str(tmp_path)
    manifest.record(cfg, "clean", "done", seconds=17.0)
    run_report.write_report(cfg)
    assert "| clean |" in (tmp_path / "REPORT.md").read_text()
    assert "17 s" in (tmp_path / "REPORT.md").read_text()
