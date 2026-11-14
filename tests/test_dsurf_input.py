"""Reading DSurfTomo projects with main.f90 semantics (single precision, grouping rules)."""

import math

import numpy as np
import pytest

from lant.inversion import dsurf_input as di

CONTROL = """c
c INPUT PARAMETERS
c
surfdata.dat                     c: data file
5 6 3                            c: nx ny nz
30.020000 99.970000              c: goxd gozd
0.010000 0.010000                c: dvxd dvzd
3                                c: max(sources, receivers)
1.0 1.0                          c: weight damp
3                                c: sublayers
1.2 5.0                          c: velocity bounds
2                                c: iterations
0.3                              c: sparsity fraction
3                                c: kmaxRc
0.2 0.5
1.0
0                                c: kmaxRg
0                                c: kmaxLc
0                                c: kmaxLg
0                                c: synthetic flag
0.0                              c: noiselevel
3.0                              c: threshold
"""


def test_control_reads_list_directed_records(tmp_path):
    (tmp_path / "DSurfTomo.in").write_text(CONTROL)
    c = di.read_control(tmp_path / "DSurfTomo.in")
    assert (c.nx, c.ny, c.nz, c.nsrc, c.kmax) == (5, 6, 3, 3, 3)
    # Periods may continue on the next line, as a Fortran list-directed READ does.
    assert c.periods_rc.tolist() == [0.2, 0.5, 1.0]
    assert c.goxd.dtype == np.float32 and c.goxd == np.float32(30.02)


def test_measurements_follow_main_grouping_and_overwrite(tmp_path):
    (tmp_path / "DSurfTomo.in").write_text(CONTROL)
    (tmp_path / "surfdata.dat").write_text(
        "# 30.0 100.0 1 2 0\n30.0 100.01 2.0\n30.0 100.02 2.0\n"
        "# 30.0 100.01 1 2 0\n30.0 100.02 2.1\n"
        "# 30.0 100.0 2 2 0\n30.0 100.03 2.2\n"
        # Period 1 again after period 2: main.f90 restarts its source counter.
        "# 30.0 100.03 1 2 0\n30.0 100.0 2.3\n"
    )
    c = di.read_control(tmp_path / "DSurfTomo.in")
    m = di.read_measurements(tmp_path / "surfdata.dat", c)
    assert m.nsrc1.tolist() == [1, 1, 0]
    assert m.nrc1[0, :2].tolist() == [1, 1]  # the first source and its receiver count were overwritten
    assert m.calsurfg_rows().tolist() == [[1, 0, 0], [2, 0, 0]]
    # The overwritten slot points at the last file row; its stale second receiver stays unused.
    assert m.file_row[0, 0, :2].tolist() == [4, 1]
    assert len(m.observed) == 5 and m.observed.dtype == np.float32


def test_single_precision_coordinates_and_haversine():
    x = di.colatitude_radians(30.0)
    assert x.dtype == np.float32
    assert x == (np.float32(90) - np.float32(30)) * di.PI32 / np.float32(180)
    a = (di.colatitude_radians(30.0), di.longitude_radians(100.0))
    b = (di.colatitude_radians(30.0), di.longitude_radians(100.01))
    d = di.delsph(*a, *b)
    expected = 6371 * math.radians(0.01) * math.cos(math.radians(30))
    # REAL radians of 100.01 deg carry ~1e-7 rad error, i.e. ~0.6 m here: the reference's own precision.
    assert d.dtype == np.float32 and abs(float(d) - expected) < 2e-3


def test_model_records_may_span_lines(tmp_path):
    values = np.arange(1, 1 + 3 * 6 * 5) / 10 + 1
    text = "0 0.5 2\n" + "\n".join(" ".join("%g" % v for v in values[i : i + 5]) for i in range(0, values.size, 5))
    (tmp_path / "MOD").write_text(text + "\n")
    z, vs = di.read_model(tmp_path / "MOD", 5, 6, 3)
    assert z.tolist() == [0, 0.5, 2] and vs.shape == (3, 6, 5)
    assert vs[1, 2, 3] == np.float32(values[1 * 30 + 2 * 5 + 3])
    with pytest.raises(ValueError, match="ended"):
        di.read_model(tmp_path / "MOD", 5, 7, 3)
