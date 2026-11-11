"""Progress stays on its chosen stream and never invents totals."""

import io

from lant.progress import Progress, track


def test_nonterminal_progress_and_counts(monkeypatch):
    monkeypatch.delenv("LANT_PROGRESS", raising=False)
    output = io.StringIO()
    with Progress("cc", stream=output, interval=100) as p:
        assert list(track([1, 2], 2, "days")) == [1, 2]
        assert p.done == 2
    text = output.getvalue()
    assert "running" in text and "2/2 days" in text
    assert "ETA" not in text and "\r" not in text and "\033" not in text


def test_progress_can_be_disabled(monkeypatch):
    monkeypatch.setenv("LANT_PROGRESS", "off")
    output = io.StringIO()
    with Progress("unknown", stream=output, interval=100):
        pass
    assert output.getvalue() == ""
