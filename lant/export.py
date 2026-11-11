"""A single dispersion table with explicit QC/selection flags and solver rounding.

The HDF5 measurement remains the full audit product. The CSV retains measured values
and separate inversion columns matching the previous text export exactly. Jackknife
uncertainty is conditional, uncalibrated and not used as a weight by DSurfTomo.
"""
import csv
import os

import h5py
import numpy as np

from .signal import runs_of

FIELDS = ("pair", "lon1", "lat1", "lon2", "lat2", "frequency_hz", "period_s", "velocity_kms", "snr",
          "uncertainty_kms", "accepted_qc", "used_in_inversion", "inversion_period_s", "inversion_velocity_kms",
          "inversion_snr", "inversion_uncertainty_kms")


def longest_run(ok, snr):
    rr = runs_of(ok)
    if not rr:
        return ok
    k0, k1 = max(rr, key=lambda t: (t[1] - t[0], np.nanmean(snr[t[0]:t[1]])))
    out = np.zeros_like(ok)
    out[k0:k1] = True
    return out


def export(cfg, disp_path, stack_path, out_dir):
    e = cfg["export"]
    with h5py.File(disp_path, "r") as h:
        freqs = h["freqs"][:]
        pairs = [p.decode() for p in h["pairs"][:]]
        vc, snr, jc, okc = (h[k][:] for k in ("vc", "snr", "jk_c", "ok_c"))
    with h5py.File(stack_path, "r") as h:
        coords = [h[k][:] for k in ("lon1", "lat1", "lon2", "lat2")]
    inside = (1 / freqs >= e["tmin_s"] - 1e-9) & (1 / freqs <= e["tmax_s"] + 1e-9)
    ncurves, npoints, dropped, several = 0, 0, 0, 0
    os.makedirs(out_dir, exist_ok=True)
    target = os.path.join(out_dir, "curves.csv")
    with open(target + ".tmp", "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(FIELDS)
        for i, pair in enumerate(pairs):
            selected = okc[i] & inside
            if e["one_segment"]:
                keep = longest_run(selected, snr[i])
                if len(runs_of(selected)) > 1:
                    several += 1
                    dropped += int(selected.sum() - keep.sum())
                selected = keep
            ncurves += int(selected.any())
            npoints += int(selected.sum())
            geometry = ["%.8f" % a[i] for a in coords]
            for k, f in enumerate(freqs):
                used = bool(selected[k])
                writer.writerow([pair] + geometry + [repr(float(f)), repr(float(1 / f)), repr(float(vc[i, k])),
                    repr(float(snr[i, k])), repr(float(jc[i, k])), int(okc[i, k]), int(used),
                    "%.6f" % (1 / f), "%.5f" % vc[i, k] if used else "0.000",
                    "%.3f" % snr[i, k] if used else "0.000", "%.5f" % jc[i, k] if used else "0.0000"])
    os.replace(target + ".tmp", target)
    return dict(phase_curves=ncurves, phase_points=npoints, curves_with_several_segments=several,
                points_dropped_by_one_segment=dropped)


def read_curves(table, kind="C", tmin=0.0, tmax=99.0):
    """Return solver-ready curves from selected rows, preserving legacy rounding and order."""
    if kind != "C":
        raise ValueError("only phase-velocity curves (C) are supported")
    grouped = {}
    with open(table, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["used_in_inversion"] != "1":
                continue
            t = float(row["inversion_period_s"])
            if not (tmin - 1e-9 <= t <= tmax + 1e-9):
                continue
            name = row["pair"]
            curve = grouped.setdefault(name, dict(name=name, **{k: float(row[k]) for k in ("lon1", "lat1", "lon2", "lat2")},
                                                  values=[]))
            curve["values"].append((t, float(row["inversion_velocity_kms"]), float(row["inversion_snr"]),
                                    float(row["inversion_uncertainty_kms"])))
    curves = []
    for name in sorted(grouped):
        curve = grouped[name]
        values = np.array(sorted(curve.pop("values")))
        curve.update(zip(("T", "v", "snr", "unc"), values.T))
        curves.append(curve)
    return curves
