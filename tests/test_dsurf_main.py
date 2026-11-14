"""DSurfTomo's main program in Python: outputs and the synthetic branch (the identity with the native binary was
shown with comparison tools kept outside the repository, docs/python_backends.md)."""

import io

import numpy as np
import pytest

pytest.importorskip("numba")
from lant.inversion import dsurf_input, dsurf_main  # noqa: E402
from tests.test_dsurf_forward import uniform_kernels, write_project  # noqa: E402


def project(tmp_path, iterations=2, synthetic=False):
    folder = write_project(tmp_path / "p")
    text = (
        (folder / "DSurfTomo.in")
        .read_text()
        .replace("1                                c: iterations", "%d c" % iterations)
    )
    if synthetic:
        text = text.replace("0                                c: synthetic flag", "1 c")
        model = (folder / "MOD").read_text().splitlines()
        (folder / "MOD.true").write_text("\n".join(model[1:]) + "\n")
    (folder / "DSurfTomo.in").write_text(text)
    return folder


def fake_kernels(folder):
    p = dsurf_input.read_project(folder)

    def kernels(depths, vs):
        k = uniform_kernels(p, phase=2.0)
        k.phase = k.phase + 0.1 * (vs[0] - 2.0)[None]  # depends on the model, so iterations change
        return k

    return kernels


def test_writes_native_files_and_log_lines(tmp_path):
    folder = project(tmp_path)
    screen = io.StringIO()
    timings = dsurf_main.run(folder, "DSurfTomo.in", kernels=fake_kernels(folder), screen=screen)
    assert len(timings) == 2 and {"senk", "traveltime", "lsmr"} <= set(timings[0])
    names = {p.name for p in folder.iterdir()}
    assert {"DSurfTomo.inMeasure.dat", "DSurfTomo.inMeasure.dat.iter001", "DSurfTomo.inMeasure.dat.iter002"} <= names
    assert {"residualFirst.dat", "residualLast.dat", "DSurfTomo.in.log"} <= names
    rows = (folder / "DSurfTomo.inMeasure.dat").read_text().splitlines()
    assert len(rows) == 5 * 6 * 3 and len(rows[0]) == 40
    log = (folder / "DSurfTomo.in.log").read_text()
    assert log.count("mean,std_devs and rms of residual:") == 2 and "Maximum and Average DWS values:" in log
    assert " mean,std_devs and rms of residual after weighting: " in screen.getvalue()
    residual = np.loadtxt(folder / "residualFirst.dat")
    assert residual.shape == (12, 6) and set(residual[:, 5]) <= {0.0, 1.0}


def test_synthetic_branch_writes_true_and_inverted_models(tmp_path):
    folder = project(tmp_path, iterations=1, synthetic=True)
    dsurf_main.run(folder, "DSurfTomo.in", kernels=fake_kernels(folder), screen=io.StringIO())
    names = {p.name for p in folder.iterdir()}
    assert {"Vs_model.real", "DSurfTomo.inSyn.dat", "velmap2dRc.dat"} <= names
    assert "DSurfTomo.inMeasure.dat" not in names
    # velmap2dRc.dat: one '(5f8.4)' line per interior node and period.
    assert len((folder / "velmap2dRc.dat").read_text().splitlines()) == 2 * 5 * 6


def test_stops_like_native_when_the_sparsity_fraction_is_too_small(tmp_path):
    folder = project(tmp_path)
    text = (folder / "DSurfTomo.in").read_text().replace("0.5                              c: sparsity", "0.001 c")
    (folder / "DSurfTomo.in").write_text(text)
    with pytest.raises(SystemExit, match="increase sparsity fraction"):
        dsurf_main.run(folder, "DSurfTomo.in", kernels=fake_kernels(folder), screen=io.StringIO())
