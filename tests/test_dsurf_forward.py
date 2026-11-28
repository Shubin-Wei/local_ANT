"""DSurfTomo-compatible forward: end-to-end rows (the identity with native CalSurfG was shown with comparison tools
kept outside the repository, docs/python_backends.md)."""

import numpy as np
import pytest

pytest.importorskip("numba")
from lant.inversion import dsurf_forward as forward  # noqa: E402
from lant.inversion import dsurf_input  # noqa: E402

CONTROL = """c
c INPUT PARAMETERS
c
surfdata.dat                     c: data file
7 8 4                            c: nx ny nz
30.060000 99.970000              c: goxd gozd
0.010000 0.010000                c: dvxd dvzd
4                                c: max(sources, receivers)
1.0 1.0                          c: weight damp
3                                c: sublayers
1.0 5.0                          c: velocity bounds
1                                c: iterations
0.5                              c: sparsity fraction
2                                c: kmaxRc
0.3 0.8
0                                c: kmaxRg
0                                c: kmaxLc
0                                c: kmaxLg
0                                c: synthetic flag
0.0                              c: noiselevel
3.0                              c: threshold
"""
STATIONS = [(30.0513, 99.9843), (30.0337, 100.0062), (30.0471, 100.0128), (30.0259, 99.9911)]


def write_project(folder):
    folder.mkdir(exist_ok=True)
    (folder / "DSurfTomo.in").write_text(CONTROL)
    lines = []
    for period in (1, 2):
        for i, (a, b) in enumerate(STATIONS):
            lines.append("# %.4f %.4f %d 2 0" % (a, b, period))
            lines.extend("%.4f %.4f %.3f" % (c, d, 2.0 + 0.1 * period) for c, d in STATIONS[i + 1 :])
            if i == len(STATIONS) - 2:
                break
    (folder / "surfdata.dat").write_text("\n".join(lines) + "\n")
    z = [0.0, 0.2, 0.6, 1.5]
    rng = np.random.default_rng(3)
    vs = 1.6 + 0.4 * np.arange(4)[:, None, None] + 0.15 * rng.random((4, 8, 7))
    text = " ".join("%g" % v for v in z) + "\n"
    text += "\n".join(" ".join("%.4f" % v for v in vs[k, j]) for k in range(4) for j in range(8))
    (folder / "MOD").write_text(text + "\n")
    return folder


def uniform_kernels(project, phase=2.0):
    c = project.control
    shape = (len(c.periods_rc), c.ny, c.nx)
    sens = np.full(shape + (c.nz,), 0.2)
    return forward.Kernels(np.full(shape, phase), sens, 0.5 * sens, 0.25 * sens)


@pytest.mark.parametrize("real", [np.float32, np.float64])
def test_rows_follow_calsurfg_order_and_uniform_times(tmp_path, real):
    project = dsurf_input.read_project(write_project(tmp_path / "p"))
    result = forward.forward(project, real=real, kernels=uniform_kernels(project))
    m = project.measurements
    assert result.rows.tolist()[:4] == [[1, 0, 0], [1, 0, 1], [1, 0, 2], [1, 1, 0]]
    assert result.file_rows.tolist() == list(range(len(m.observed)))
    assert result.times.dtype == real and result.jacobian.shape == (12, 5 * 6 * 3)
    # Uniform 2 km/s: the reference's FMM on this 140 m mesh agrees with distance / 2 to about 1 %.
    assert np.allclose(result.times, m.distance / 2.0, rtol=1e-2)
    assert np.all(result.jacobian.toarray() <= 0)


def test_rejects_other_wave_types(tmp_path):
    folder = write_project(tmp_path / "p")
    text = (folder / "surfdata.dat").read_text().replace("1 2 0", "1 1 0", 1)  # Love
    (folder / "surfdata.dat").write_text(text)
    with pytest.raises(NotImplementedError, match="Rayleigh phase"):
        dsurf_input.read_project(folder)


def write_group_project(folder):
    """The same pairs with the second period as Rayleigh group velocity (kind 1, period index 1 of kmaxRg)."""
    write_project(folder)
    text = (folder / "surfdata.dat").read_text().replace(" 2 2 0\n", " 1 2 1\n")
    (folder / "surfdata.dat").write_text(text)
    ctl = (folder / "DSurfTomo.in").read_text()
    ctl = ctl.replace(
        "2                                c: kmaxRc\n0.3 0.8\n0                                c: kmaxRg",
        "1                                c: kmaxRc\n0.3\n1                                c: kmaxRg\n0.8",
    )
    (folder / "DSurfTomo.in").write_text(ctl)
    return folder


def test_group_blocks_follow_the_phase_blocks(tmp_path):
    project = dsurf_input.read_project(write_group_project(tmp_path / "p"))
    c, m = project.control, project.measurements
    assert list(c.periods_rc) == [0.3] and list(c.periods_rg) == [0.8] and list(c.periods_all) == [0.3, 0.8]
    assert m.nsrc1.tolist() == [3, 3]  # knum 1: phase, knum 2: group


def test_group_rows_are_the_derivative_of_the_group_times(tmp_path):
    """G of the group travel times is d(time)/d(Vs) of the times the forward computes (all in the group-velocity
    map): a 1 % change of one node's Vs changes the times as G predicts, to the second-order ray-path effect."""
    from lant.inversion import dsurf_kernels

    project = dsurf_input.read_project(write_group_project(tmp_path / "p"))
    c = project.control
    k0 = dsurf_kernels.depth_kernels(project.depths, project.vs, c.periods_all, c.minthk, nrc=len(c.periods_rc))
    base = forward.forward(project, real=np.float64, kernels=k0)
    G = base.jacobian.toarray()
    group = base.rows[:, 0] == 2
    col = int(np.argmax(np.abs(G[group]).sum(axis=0)))  # the parameter the group rays sense most
    nvx, nvz = c.nx - 2, c.ny - 2
    kz, rest = divmod(col, nvz * nvx)
    j, i = divmod(rest, nvx)
    vs = project.vs.copy()
    dv = np.float32(0.01) * vs[kz, j + 1, i + 1]
    vs[kz, j + 1, i + 1] += dv
    moved = dsurf_input.Project(c, project.measurements, project.depths, vs)
    k1 = dsurf_kernels.depth_kernels(moved.depths, moved.vs, c.periods_all, c.minthk, nrc=len(c.periods_rc))
    after = forward.forward(moved, real=np.float64, kernels=k1, derivatives=False)
    predicted = G[group, col] * float(dv)
    actual = (after.times - base.times)[group]
    assert np.abs(predicted).max() > 0
    assert np.allclose(actual, predicted, rtol=0.05, atol=0.02 * np.abs(predicted).max())


def test_partials_only_where_rays_go_give_the_same_matrix(tmp_path, monkeypatch):
    from lant.inversion import dsurf_kernels

    project = dsurf_input.read_project(write_project(tmp_path / "p"))
    c = project.control
    full = dsurf_kernels.depth_kernels(project.depths, project.vs, c.periods_rc, c.minthk)
    lazy = dsurf_kernels.LazyKernels(project.depths, project.vs, c.periods_rc, c.minthk)
    monkeypatch.setattr(forward, "FIELD_BLOCK_BYTES", 1.0)  # one source per block: partials asked for block by block
    a = forward.forward(project, kernels=full)
    b = forward.forward(project, kernels=lazy)
    assert np.array_equal(a.times, b.times) and np.array_equal(full.phase, lazy.phase)
    assert (a.jacobian != b.jacobian).nnz == 0 and a.jacobian.nnz == b.jacobian.nnz
    interior = np.zeros(lazy.done.shape, bool)
    interior[1:-1, 1:-1] = True
    assert lazy.done.any() and not lazy.done[~interior].any()
    assert np.array_equal(lazy.vs[:, lazy.done], full.vs[:, lazy.done])
