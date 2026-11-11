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
