"""CLI boundaries: machine output, read-only status and explicit destructive intent."""

import io
import json

import pytest

from lant import cli


def test_note_and_json_status_have_one_output_contract(tmp_path, capsys):
    cfg, out = tmp_path / "config.yaml", tmp_path / "out"
    cli.main(
        [
            "init",
            "-o",
            str(cfg),
            "--out",
            str(out),
            "--raw",
            str(tmp_path / "raw"),
            "--stations",
            str(tmp_path / "stations.csv"),
        ]
    )
    capsys.readouterr()
    cli.main(["note", "user supplied a velocity prior", "-c", str(cfg)])
    capsys.readouterr()
    assert "user supplied a velocity prior" in (out / "notes.md").read_text()
    unexpected = out / "results" / "extra.txt"
    unexpected.parent.mkdir(exist_ok=True)
    unexpected.write_text("preserve")
    before = {p: p.stat().st_mtime_ns for p in out.rglob("*")}
    cli.main(["status", "-c", str(cfg), "--json"])
    status = json.loads(capsys.readouterr().out)
    assert status["artifacts"]["unregistered"] == ["results/extra.txt"]
    assert before == {p: p.stat().st_mtime_ns for p in out.rglob("*")}


@pytest.mark.parametrize("command", [["run", "all", "--rebuild"], ["clean-work", "--execute"]])
def test_noninteractive_rebuild_and_cleanup_require_yes(tmp_path, monkeypatch, command):
    cfg = tmp_path / "config.yaml"
    cfg.write_text("project: {output_dir: out}\ndata: {raw_dir: raw, station_file: stations.csv}\n")
    monkeypatch.setattr(cli.sys, "stdin", io.StringIO())
    with pytest.raises(SystemExit, match="requires --yes"):
        cli.main(command + ["-c", str(cfg)])
    assert not (tmp_path / "out").exists()


def test_install_external_defaults_to_both_programs(monkeypatch):
    called = []
    monkeypatch.setattr(cli, "cmd_install_external", lambda args: called.append(args.programs))
    cli.main(["install-external"])
    assert called == [[]]


def test_multiline_dispatch_preserves_execution_options(monkeypatch):
    import subprocess

    called = []
    monkeypatch.setattr(cli, "_quick_config", lambda args: ["line1.yaml", "line2.yaml"])
    monkeypatch.setattr(subprocess, "call", lambda command: called.append(command) or 0)
    with pytest.raises(SystemExit) as result:
        cli.main(
            [
                "run",
                "all",
                "--raw",
                "raw",
                "--stations",
                "stations.csv",
                "--out",
                "output",
                "--from",
                "reference",
                "--to",
                "inversion",
                "--what",
                "main",
                "--verbose",
                "--rebuild",
                "--yes",
            ]
        )
    assert result.value.code == 0 and len(called) == 2
    for command in called:
        for option, value in (("--from", "reference"), ("--to", "inversion"), ("--what", "main")):
            assert command[command.index(option) + 1] == value
        assert all(option in command for option in ("--verbose", "--rebuild", "--yes"))


@pytest.mark.parametrize("reason", ["USER_DECISION_AND_EVIDENCE", "  "])
def test_confirmation_needs_a_reason_of_the_user(tmp_path, reason, capsys):
    with pytest.raises(SystemExit) as stop:
        cli.main(["confirm", "clean", "-c", str(tmp_path / "config.yaml"), "--reason", reason])
    assert stop.value.code == 2 and "own words" in capsys.readouterr().err


def test_a_confirmation_is_marked_as_the_users_decision(monkeypatch):
    from lant import inspect, manifest, report

    monkeypatch.setattr(report, "resolved", lambda cfg: cfg)
    monkeypatch.setattr(manifest, "status", lambda cfg: {"stack": dict(state="to confirm")})
    monkeypatch.setattr(inspect, "active_writer", lambda rec: None, raising=False)
    import lant.artifacts
    import lant.provenance

    monkeypatch.setattr(lant.provenance, "active_writer", lambda rec: None)
    monkeypatch.setattr(lant.artifacts, "audit", lambda cfg: dict(unregistered=[], missing=[], modified=[]))
    st = inspect.status({"project": {"output_dir": "/nonexistent"}, "_config_file": "c.yaml"})
    assert st["next"]["command"].startswith("lant confirm stack") and st["next"]["requires_user"]


def test_note_does_not_wait_for_a_running_step(tmp_path, capsys):
    from lant.provenance import Lock

    cfg, out = tmp_path / "config.yaml", tmp_path / "out"
    cfg.write_text("project: {output_dir: out}\ndata: {raw_dir: raw, station_file: stations.csv}\n")
    (out / ".lant").mkdir(parents=True)
    with Lock(str(out / ".lant"), "cc"):  # a step is writing
        cli.main(["note", "noise source seen near station A01", "-c", str(cfg)])
        cli.main(["note", "second note", "-c", str(cfg)])
    text = (out / "notes.md").read_text()
    assert text.count("# Notes and decisions") == 1
    assert text.index("noise source seen") < text.index("second note")
    assert not (out / ".lant/run_log.jsonl").exists()  # managed records wait for the writer


def test_finished_chain_says_where_to_look(tmp_path):
    import json

    from lant import delivery

    rec = tmp_path / ".lant/summary"
    rec.mkdir(parents=True)
    summaries = {
        "clean": {"stations": 10, "days": 1},
        "disp": {"export": {"phase_curves": 35, "phase_points": 840}},
        "evaluate": {
            "residuals": {"final": {"rms_ms": 8.14, "n": 810}, "iterations": [{"model": "start", "rms_ms": 20.68}]},
            "checkerboards": {"cb_1": {}},
        },
        "figures": {"drawn": {"stations": []}},
    }
    for step, data in summaries.items():
        (rec / ("%s.json" % step)).write_text(json.dumps(data))
    lines = delivery.finish_lines({"project": {"output_dir": str(tmp_path)}})
    assert lines[0].endswith("REPORT.md") and "results" in lines[1]
    assert lines[2] == (
        "  10 stations, 1 processing day; 35 phase-velocity curves (840 points); residual rms of all 810 "
        "travel times 20.7 -> 8.1 ms; 1 checkerboard test; 1 figure"
    )


def test_version_and_guidance_without_a_command(capsys):
    from lant import __version__

    with pytest.raises(SystemExit) as stop:
        cli.main(["--version"])
    assert stop.value.code == 0 and capsys.readouterr().out.strip() == "lant %s" % __version__
    cli.main([])
    out = capsys.readouterr().out
    assert "lant demo" in out and "First time: lant check, then lant demo." in out


def test_messages_say_what_to_do(tmp_path, capsys):
    with pytest.raises(SystemExit, match="lant init -o"):
        cli.main(["status", "-c", str(tmp_path / "missing.yaml")])
    typo = tmp_path / "typo.yaml"
    typo.write_text("project: {output_dir: out}\ndata: {raw_dir: r, station_file: s, velocty_range_kms: [1, 3]}\n")
    with pytest.raises(SystemExit, match="did you mean 'data.velocity_range_kms'"):
        cli.main(["validate", "-c", str(typo)])
    template = tmp_path / "template.yaml"
    cli.main(["init", "-o", str(template)])
    capsys.readouterr()
    with pytest.raises(SystemExit, match="template paths"):
        cli.main(["validate", "-c", str(template)])


def test_demo_never_works_in_a_directory_it_did_not_make(tmp_path):
    from lant import demo

    (tmp_path / "real.mseed").write_bytes(b"data of a user")
    with pytest.raises(RuntimeError, match="not made by lant demo"):
        demo.run(tmp_path)
    with pytest.raises(RuntimeError, match="not a generated synthetic example"):
        demo.run_known_truth(tmp_path)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["real.mseed"]


def test_init_creates_the_directory_of_the_configuration(tmp_path, capsys):
    config = tmp_path / "result" / "config.yaml"  # the output directory does not exist yet
    cli.main(["init", "-o", str(config), "--out", str(tmp_path / "result"), "--raw", "raw", "--stations", "s.csv"])
    assert config.is_file() and "output_dir" in config.read_text()
