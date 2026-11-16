"""Synthetic figure labels and cell size must stay tied to the published benchmark."""

import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from lant import __version__
from lant.plotting import style
from tools import plot_validation_comparison as comparison


def test_release_versions_match_package_and_citation():
    import tomllib

    root = Path(__file__).parents[1]
    metadata = tomllib.loads((root / "pyproject.toml").read_text())
    citation = yaml.safe_load((root / "CITATION.cff").read_text())
    assert metadata["project"]["version"] == citation["version"] == __version__


def test_synthetic_pattern_keeps_the_original_three_km_cells():
    pattern = dict(size_km=9, cell_km=3, amplitude=0.06)
    values = comparison.truth_anomaly(np.array([1.5, 4.5, 7.5, 10]), np.full(4, 1.5), pattern)
    assert np.allclose(values, [6, -6, 6, 0])
    assert comparison.TITLES == ("True", "lant", "DSurfTomo")
    assert np.allclose(style.PERTURBATION(0.5)[:3], [1, 1, 1])


def test_redraw_refuses_to_substitute_a_denser_synthetic_case(tmp_path):
    source = tmp_path / "inv_areal"
    source.mkdir()
    (source / "true_checkerboard.json").write_text(json.dumps(dict(cell_km=1.5)))
    with pytest.raises(ValueError, match="3-km cells"):
        comparison.areal(tmp_path, tmp_path)
