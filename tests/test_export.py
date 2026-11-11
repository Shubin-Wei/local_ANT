"""Selection and numeric semantics of the solver-facing dispersion table."""
import csv

import h5py
import numpy as np

from lant import config, export


def test_table_retains_qc_points_but_only_exports_selected_segment(tmp_path):
    cfg = config.defaults()
    cfg["export"].update(tmin_s=0.1, tmax_s=2, one_segment=True)
    with h5py.File(str(tmp_path / "disp.h5"), "w") as h:
        h["freqs"] = [1., 2., 3., 4., 5.]
        h["pairs"] = np.array([b"A_B"])
        h["vc"] = [[1.234567, 1.234567, 1.234567, 1.234567, 1.234567]]
        h["snr"] = [[10.12345, 10.12345, 1., 20.12345, 20.12345]]
        h["jk_c"] = [[0.01234567] * 5]
        h["ok_c"] = [[True, True, False, True, True]]
    with h5py.File(str(tmp_path / "stack.h5"), "w") as h:
        for k, value in zip(("lon1", "lat1", "lon2", "lat2"), (100.123456789, 30., 100.2, 30.1)):
            h[k] = [value]
    result = export.export(cfg, str(tmp_path / "disp.h5"), str(tmp_path / "stack.h5"), str(tmp_path / "out"))
    table = tmp_path / "out/curves.csv"
    rows = list(csv.DictReader(table.open()))
    assert len(rows) == 5 and sum(int(r["accepted_qc"]) for r in rows) == 4
    assert result["phase_points"] == 2 and result["points_dropped_by_one_segment"] == 2
    curves = export.read_curves(table)
    assert np.array_equal(curves[0]["T"], [0.2, 0.25])
    assert np.array_equal(curves[0]["v"], [1.23457, 1.23457])
    assert curves[0]["lon1"] == 100.12345679
    assert float(rows[0]["velocity_kms"]) == 1.234567
    assert sorted(p.name for p in (tmp_path / "out").iterdir()) == ["curves.csv"]
