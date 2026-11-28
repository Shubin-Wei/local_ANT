"""surfdisp96 of DSurfTomo, phase (igr = 0) and group (igr = 1) velocities, bit for bit against the original Fortran.

tests/data/native_surfdisp96.npz: 40 random layer models (3-59 layers, some with low-velocity layers, Brocher Vp and
density) and 35 periods of 0.05-3 s, run through DSurfTomo's surfdisp96.f compiled with its Makefile flags
(gfortran -O -ffixed-line-length-none -ffloat-store -fbounds-check -m64 -mcmodel=medium) by a driver that calls
surfdisp96(thkm, vpm, vsm, rhom, nlayer, 1, 2, 1, igr, kmax, t, cg).
"""

import os

import numpy as np

from lant.inversion import dsurf_disp96 as D

DATA = os.path.join(os.path.dirname(__file__), "data", "native_surfdisp96.npz")


def test_phase_and_group_velocities_equal_the_fortran():
    z = np.load(DATA)
    t = z["periods"]
    for i, n in enumerate(z["nlayer"]):
        for igr, key in ((0, "phase"), (1, "group")):
            cg = np.zeros(len(t))
            D.surfdisp96(z["thk"][i, :n], z["vp"][i, :n], z["vs"][i, :n], z["rho"][i, :n], int(n), t, cg, igr)
            assert np.array_equal(cg, z[key][i]), (i, key)
