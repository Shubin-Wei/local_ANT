"""Step reports (lant.report)."""

import os

from lant import report


def test_key_numbers_show_project_paths_relative_to_the_root(tmp_path):
    """Summary keys that hold a path inside the project are shown relative to its root, not as absolute paths."""
    root = str(tmp_path / "out")
    main = os.path.join(root, "work", "06_inversion", "main")
    nums = dict(report._numbers({"results": {main: "ok"}, "noise": 0.05}, root=root))
    assert nums == {"results.work/06_inversion/main": "ok", "noise": 0.05}
    assert dict(report._numbers({"results": {main: "ok"}}))["results." + main] == "ok"  # no root: unchanged
