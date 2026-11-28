"""surfdisp96 of DSurfTomo, phase (igr = 0) and group (igr = 1) velocities, bit for bit against the original Fortran.

tests/data/native_surfdisp96.npz: 40 random layer models (3-59 layers, some with low-velocity layers, Brocher Vp and
density) and 35 periods of 0.05-3 s, run through DSurfTomo's surfdisp96.f compiled with its Makefile flags
(gfortran -O -ffixed-line-length-none -ffloat-store -fbounds-check -m64 -mcmodel=medium) by a driver that calls
surfdisp96(thkm, vpm, vsm, rhom, nlayer, 1, 2, 1, igr, kmax, t, cg): phase, group with its root-search step
ddc0 = 0.005; phase_fine, group_fine with ddc0 = 0.0005, the step lant takes where Vs decreases with depth.
"""

import os

import numpy as np

from lant.inversion import dsurf_disp96 as D

DATA = os.path.join(os.path.dirname(__file__), "data", "native_surfdisp96.npz")


def test_phase_and_group_velocities_equal_the_fortran():
    z = np.load(DATA)
    t = z["periods"]
    for i, n in enumerate(z["nlayer"]):
        vs = z["vs"][i, :n]
        step = D.root_step(vs)
        fine = "_fine" if step != D._DC else ""
        assert bool(fine) == bool(np.any(np.diff(vs) < 0))
        for igr, key in ((0, "phase"), (1, "group")):
            cg = np.zeros(len(t))
            D.surfdisp96(z["thk"][i, :n], z["vp"][i, :n], vs, z["rho"][i, :n], int(n), t, cg, igr, step)
            assert np.array_equal(cg, z[key + fine][i]), (i, key + fine)


def test_the_fine_step_finds_the_fundamental_mode_under_a_low_velocity_layer():
    """A column of a real inversion: with DSurfTomo's step the search jumps to a higher mode at 0.105 s and finds no
    root from 0.25 s on; with the fine step the curve is the smooth fundamental mode."""
    z = np.array([0, 0.054, 0.11, 0.16, 0.21, 0.26, 0.33, 0.41, 0.51, 0.64, 0.8, 1.0, 1.2], np.float32)
    vs = np.array(
        [2.803, 2.705, 2.607, 2.538, 2.494, 2.475, 2.507, 2.565, 2.614, 2.576, 2.42, 2.249, 2.597], np.float32
    )
    t = np.array([0.08839, 0.105112, 0.125, 0.148651, 0.176777, 0.210224, 0.25, 0.297302, 0.353553, 0.420448, 0.5])
    c = D.column_phase(z, vs, np.float32(3), t)
    assert np.all(c > 0) and np.all(np.abs(np.diff(c)) < 0.03)
