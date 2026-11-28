"""Rayleigh modes of the 1-D models of models.txt, computed with Computer Programs in Seismology (Herrmann 2013:
sprep96, sdisp96, sregn96, sdpder96) and kept as tables (modes_<model>.csv) so that the waveforms can be made
without CPS.

    python modes.py [CPS bin folder]     -> modes_layered.csv, modes_lvl.csv

Columns: mode (0 = fundamental), frequency_hz, c_kms (phase velocity), u_kms (group velocity), ar (the amplitude
factor of the mode for a source and a receiver at the surface: CPS normalises the vertical eigenfunction to 1 at
the surface and puts the excitation into AR; Herrmann's far-field vertical displacement of a vertical point force is
AR / sqrt(8 pi k r) exp(-i (k r - pi / 4)), i.e. AR / 4 times the Hankel function H0^(2)(k r)).  Elastic
(sregn96 -NOQ): attenuation is applied when the waveforms are made.
"""

import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FMIN, FMAX, NFREQ, NMOD = 0.15, 25.0, 400, 20


def read_models(path=os.path.join(HERE, "models.txt")):
    """-> {name: array [layers, 4] of thickness, Vp, Vs, density}."""
    out, name = {}, None
    for line in open(path):
        s = line.split("#")[0].strip()
        if not s:
            continue
        if s.startswith("["):
            name = s.strip("[]")
            out[name] = []
        else:
            out[name].append([float(v) for v in s.split()])
    return {k: np.array(v) for k, v in out.items()}


def _cps_model(layers, path):
    with open(path, "w") as f:
        f.write("MODEL.01\nsynthetic\nISOTROPIC\nKGS\nFLAT EARTH\n1-D\nCONSTANT VELOCITY\nLINE08\nLINE09\nLINE10\n"
                "LINE11\n H VP VS RHO QP QS ETAP ETAS FREFP FREFS\n")  # fmt: skip
        for h, a, b, r in layers:
            f.write("%10.4f %10.4f %10.4f %10.4f 0 0 0 0 1 1\n" % (h, a, b, r))


def cps_table(layers, cps_bin):
    """-> rows (mode, f, c, u, ar) from CPS."""
    freqs = FMIN * (FMAX / FMIN) ** (np.arange(NFREQ) / (NFREQ - 1))
    with tempfile.TemporaryDirectory() as tmp:
        _cps_model(layers, os.path.join(tmp, "m.mod"))
        np.savetxt(os.path.join(tmp, "f.txt"), freqs, fmt="%.6f")
        run = lambda *a: subprocess.run([os.path.join(cps_bin, a[0]), *a[1:]], cwd=tmp, check=True,  # noqa: E731
                                        capture_output=True)  # fmt: skip
        run("sprep96", "-M", "m.mod", "-HS", "0", "-HR", "0", "-R", "-NMOD", str(NMOD), "-FARR", "f.txt")
        run("sdisp96")
        run("sregn96", "-NOQ", "-DER")
        run("sdpder96", "-R", "-TXT")
        text = open(os.path.join(tmp, "SRDER.TXT")).read().split("RAYLEIGH WAVE")[1:]
    rows = []
    for block in text:
        w = block.replace("=", " = ").split()
        mode = int(w[w.index("#") + 1])
        t, c, u, ar = (float(w[w.index(k) + 2]) for k in ("T", "C", "U", "AR"))
        rows.append((mode, 1 / t, c, u, ar))
    return np.array(sorted(rows))


def read_table(path):
    """-> {mode: (f, c, u, ar)} in increasing frequency."""
    a = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    out = {}
    for m in np.unique(a[:, 0]).astype(int):
        b = a[a[:, 0] == m]
        b = b[np.argsort(b[:, 1])]
        out[m] = tuple(b[:, k] for k in range(1, 5))
    return out


def main():
    cps_bin = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(subprocess.run(
        ["bash", "-c", "command -v sprep96"], capture_output=True, text=True).stdout.strip())  # fmt: skip
    for name, layers in read_models().items():
        rows = cps_table(layers, cps_bin)
        path = os.path.join(HERE, "modes_%s.csv" % name)
        np.savetxt(path, rows, fmt=("%d", "%.6f", "%.6f", "%.6f", "%.6e"), delimiter=",",
                   header="mode,frequency_hz,c_kms,u_kms,ar", comments="")  # fmt: skip
        print("%s: %d modes, %d rows -> %s" % (name, len(np.unique(rows[:, 0])), len(rows), path))


if __name__ == "__main__":
    main()
