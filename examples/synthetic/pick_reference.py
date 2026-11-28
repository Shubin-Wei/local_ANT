"""A reference curve picked along one ridge of the F-J spectrum of lant disp (disp.reference: file).

    python pick_reference.py OUTPUT/3_disp/fj.h5 reference.txt  f1,c1 f2,c2 ...

f,c: points read by eye along the ridge on figures/fj_reference.png (Hz, km/s), from low to high frequency.  At every
frequency between the first and the last point the curve takes the strongest local maximum of the F-J spectrum
within 8 % of the line through the points; frequencies without one are left out.  This is how
reference_layered.txt and reference_lvl.txt were made, with the points given in README.md.
"""

import sys

import h5py
import numpy as np


def main():
    with h5py.File(sys.argv[1], "r") as h:
        f, c, img = h["freqs"][:], h["c"][:], h["image"][:]
    guide = np.array([[float(x) for x in p.split(",")] for p in sys.argv[3:]])
    rows = []
    for i, fi in enumerate(f):
        if not guide[0, 0] <= fi <= guide[-1, 0]:
            continue
        g = np.exp(np.interp(np.log(fi), np.log(guide[:, 0]), np.log(guide[:, 1])))
        y = img[i] / img[i].max()
        peaks = [j for j in range(1, len(y) - 1) if y[j] >= y[j - 1] and y[j] >= y[j + 1] and abs(c[j] / g - 1) <= 0.08]
        if peaks:
            rows.append((fi, c[max(peaks, key=lambda j: y[j])]))
    np.savetxt(sys.argv[2], rows, fmt="%.4f %.3f",
               header="frequency_Hz phase_velocity_km/s (F-J ridge picked along points read by eye)")  # fmt: skip
    print("%d points -> %s" % (len(rows), sys.argv[2]))


if __name__ == "__main__":
    main()
