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
