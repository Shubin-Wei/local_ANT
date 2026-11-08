"""Export of the QC'd dispersion curves.

  06_dispersion/logf/CDisp|GDisp.PAIR.txt     every accepted sample on the measurement (log-frequency) grid:
                                               freq_Hz period_s velocity_km/s snr uncertainty_km/s
  06_dispersion/curves/CDisp|GDisp.PAIR.txt   the curves used by the inversion, on the MEASURED periods (1/f of the
                                               measurement grid between export.tmin_s and tmax_s; nothing interpolated):
                                               two header lines "lon lat" of the two stations, then
                                               period_s velocity_km/s snr uncertainty_km/s   (velocity 0 = no value)
With export.one_segment, each phase-velocity curve keeps only its longest contiguous run on the measurement grid
(a tie goes to the higher mean SNR; EGFAnalysisTimeFreq manual, Yao 2015, p. 9: "select only the most continuous
part with higher SNR"); the samples left out are written as unmeasured rows.  Group curves are exported as they are.
Uncertainty = max(jackknife std, redundant-path robust std of the band x velocity).
"""
import os

import h5py
import numpy as np

from .qc import BANDS
from .signal import runs_of


def longest_run(ok, snr):
    rr = runs_of(ok)
    if not rr:
        return ok
    k0, k1 = max(rr, key=lambda t: (t[1] - t[0], np.nanmean(snr[t[0]:t[1]])))
    out = np.zeros_like(ok)
    out[k0:k1] = True
    return out


def export(cfg, disp_path, stack_path, out_dir, precision=None):
    e = cfg["export"]
    with h5py.File(disp_path, "r") as h:
        freqs = h["freqs"][:]
        pairs = [p.decode() for p in h["pairs"][:]]
        vg, vc, snr, snr_c = h["vg"][:], h["vc"][:], h["snr"][:], h["snr_c"][:]
        jg, jc = h["jk_g"][:], h["jk_c"][:]
        okg, okc = h["ok_g"][:], h["ok_c"][:]
    with h5py.File(stack_path, "r") as h:
        lon1, lat1, lon2, lat2 = (h[k][:] for k in ("lon1", "lat1", "lon2", "lat2"))
    inside = (1 / freqs >= e["tmin_s"] - 1e-9) & (1 / freqs <= e["tmax_s"] + 1e-9)
    order = np.argsort(1 / freqs)
    order = order[inside[order]]
    floor = {"G": np.zeros(len(freqs)), "C": np.zeros(len(freqs))}
    if precision:
        for kind, key in (("G", "group"), ("C", "phase")):
            for (lo, hi) in BANDS:
                v = precision[key].get("%.1f-%.1f" % (lo, hi), {}).get("robust_std")
                m = (freqs >= lo) & ((freqs < hi) | (hi >= 10.4))
                floor[kind][m] = v if v is not None and np.isfinite(v) else 0.0
    for sub in ("logf", "curves"):
        d = os.path.join(out_dir, sub)
        os.makedirs(d, exist_ok=True)
        for old in os.listdir(d):
            os.remove(os.path.join(d, old))
    n = {"G": 0, "C": 0, "C_points": 0, "C_points_dropped": 0, "C_two_or_more_segments": 0}
    for i, p in enumerate(pairs):
        for kind, v, ok, jk0, sn in (("G", vg, okg, jg, snr), ("C", vc, okc, jc, snr_c)):
            m = ok[i] & inside
            if not m.any():
                continue
            if kind == "C" and e["one_segment"]:
                keep = longest_run(m, sn[i])
                if len(runs_of(m)) > 1:
                    n["C_two_or_more_segments"] += 1
                    n["C_points_dropped"] += int(m.sum() - keep.sum())
                m = keep
            jk = np.maximum(jk0[i], floor[kind] * v[i])
            head = "%.8f    %.8f\n%.8f    %.8f\n" % (lon1[i], lat1[i], lon2[i], lat2[i])
            with open(os.path.join(out_dir, "logf", "%sDisp.%s.txt" % (kind, p)), "w") as fh:
                fh.write(head)
                for k in np.where(m)[0]:
                    fh.write("%.4f %.4f %.4f %.2f %.4f\n" % (freqs[k], 1 / freqs[k], v[i, k], sn[i, k], jk[k]))
            with open(os.path.join(out_dir, "curves", "%sDisp.%s.txt" % (kind, p)), "w") as fh:
                fh.write(head)
                for k in order:
                    if m[k]:
                        fh.write("%.4f  %.3f  %.3f  %.4f\n" % (1 / freqs[k], v[i, k], sn[i, k], jk[k]))
                    else:
                        fh.write("%.4f  0.000  0.000  0.0000\n" % (1 / freqs[k]))
            n[kind] += 1
            if kind == "C":
                n["C_points"] += int(m.sum())
    return dict(phase_curves=n["C"], phase_points=n["C_points"], group_curves=n["G"],
                curves_with_several_segments=n["C_two_or_more_segments"], points_dropped_by_one_segment=n["C_points_dropped"])


def read_curves(curve_dir, kind="C", tmin=0.0, tmax=99.0):
    """-> list of dict(name, lon1, lat1, lon2, lat2, T [n], v [n]) for the exported curves (velocity > 0 only)."""
    import glob
    out = []
    for f in sorted(glob.glob(os.path.join(curve_dir, "%sDisp.*.txt" % kind))):
        L = [l.split() for l in open(f) if l.strip()]
        a = np.array(L[2:], float) if len(L) > 2 else np.zeros((0, 4))
        m = (a[:, 1] > 0) & (a[:, 0] >= tmin - 1e-9) & (a[:, 0] <= tmax + 1e-9) if len(a) else np.zeros(0, bool)
        if m.any():
            out.append(dict(name=os.path.basename(f)[6:-4], lon1=float(L[0][0]), lat1=float(L[0][1]),
                            lon2=float(L[1][0]), lat2=float(L[1][1]), T=a[m, 0], v=a[m, 1], snr=a[m, 2], unc=a[m, 3]))
    return out
