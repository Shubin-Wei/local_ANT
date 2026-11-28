"""Step 3 (lant disp): reference curve and phase velocities of every pair.

Results in OUTPUT/3_disp:
  reference.csv            the reference curve used, its group velocity and the F-J picks
  fj.h5                    the F-J spectrum (frequencies, velocities, |I|)
  measurements.h5          every pair and frequency: c, group velocity u, snr, sigma (jackknife) and the reasons of
                           every rejection (bit mask, lant.dispersion.REASONS)
  dispersion.csv           the accepted values, one row per pair and frequency: the input of `lant inv`
  median_vs_reference.csv  per frequency: accepted pairs, their median c and the reference
"""

import glob
import math
import os

import h5py
import numpy as np

from . import config as C
from . import dispersion as D
from . import reference, resources
from .signal import log_freqs
from .step import Step, need


def measurement_band(cfg, step):
    """disp.fmin_hz / fmax_hz, or the correlation band narrowed so that the Gaussian filters of the edge
    frequencies reach half gain inside it, f (1 -+ sqrt(ln 2 / alpha)); on the 1/per_octave grid through 1 Hz."""
    d, c = cfg["disp"], cfg["cc"]
    po = int(d["per_octave"])
    half = math.sqrt(math.log(2) / float(d["alpha"]))
    if d["fmin_hz"] is not None:
        lo = float(d["fmin_hz"])
    else:
        lo = step.derive(
            "disp.fmin_hz",
            round(float(c["fmin_hz"]) / (1 - half), 4),
            "cc.fmin_hz / (1 - sqrt(ln 2 / alpha)): the lowest filter keeps half gain in the band",
        )
    if d["fmax_hz"] is not None:
        hi = float(d["fmax_hz"])
    else:
        hi = step.derive(
            "disp.fmax_hz",
            round(float(c["fmax_hz"]) / (1 + half), 4),
            "cc.fmax_hz / (1 + sqrt(ln 2 / alpha)): the highest filter keeps half gain in the band",
        )
    k0, k1 = math.ceil(po * math.log2(lo) - 1e-9), math.floor(po * math.log2(hi) + 1e-9)
    if k1 - k0 < 2:
        raise C.ConfigError("measurement band %.4g-%.4g Hz holds fewer than 3 frequencies" % (lo, hi))
    return log_freqs(2.0 ** (k0 / po), 2.0 ** (k1 / po), po)


def settings(cfg):
    d = cfg["disp"]
    return dict(alpha=float(d["alpha"]), group_window=float(d["group_window"]), snr_min=float(d["snr_min"]),
                min_wavelengths=float(d["min_wavelengths"]), reference_error=float(d["reference_error"]),
                max_jackknife_rel=None if d["max_jackknife_rel"] is None else float(d["max_jackknife_rel"]),
                slip_cycles=float(d["slip_cycles"]), min_points=int(d["min_points"]))  # fmt: skip


def run(cfg):
    C.check_disp(cfg)
    stack = need(os.path.join(cfg["output_dir"], "2_cc", "stack.h5"), "cc")
    day_files = sorted(glob.glob(os.path.join(cfg["output_dir"], "2_cc", "days", "*.h5")))
    with Step(cfg, "3_disp", ["data", "disp", "cores"]) as step:
        n_cores, how = resources.cores(cfg["cores"])
        step.say("cores: %d (%s)" % (n_cores, how))
        freqs = measurement_band(cfg, step)
        step.say("measurement frequencies: %d, %.4g-%.4g Hz (%d per octave)" % (len(freqs), freqs[0], freqs[-1],
                                                                               cfg["disp"]["per_octave"]))  # fmt: skip
        with h5py.File(stack, "r") as h:
            cnt = h["count"][:]
            use = np.where(cnt > 0)[0]
            ccf, dist, fs = h["ccf"][use], h["dist_km"][use], float(h.attrs["fs"])
            sta1, sta2 = [x.decode() for x in h["sta1"][:]], [x.decode() for x in h["sta2"][:]]
            coords = {k: h[k][:] for k in ("lon1", "lat1", "lon2", "lat2")}
            dist_all = h["dist_km"][:]
        ref = reference.build(cfg, ccf, dist, fs, freqs, step.say)
        del ccf
        reference.write_csv(step.path("reference.csv"), ref)
        with h5py.File(step.path("fj.h5"), "w") as h:
            h["freqs"], h["c"], h["image"], h["picks"] = freqs, ref["fj_c"], ref["fj_image"], ref["fj_picks"]
        defined = np.isfinite(ref["c"])
        step.say("reference: %s; defined at %d of %d frequencies (%.4g-%.4g Hz)" % (
            ref["source"], defined.sum(), len(freqs), freqs[defined].min(), freqs[defined].max()))  # fmt: skip
        step.say("CHECK figures/fj_reference.png: the reference must follow the fundamental-mode ridge (the slowest "
                 "ridge at each frequency). If it does not, set disp.fj_band_hz to where it does, or give a model or a "
                 "curve (disp.reference), and run `lant disp` again.")  # fmt: skip
        nsub = sum(_n_substacks(f) for f in day_files)
        if nsub < 3 and cfg["disp"]["max_jackknife_rel"] is not None:
            raise C.ConfigError("%d sub-stacks: the jackknife needs at least 3; use shorter cc.substack_hours (and run "
                                "lant cc again) or set disp.max_jackknife_rel to null" % nsub)  # fmt: skip
        step.say("measuring %d pairs (%d jackknife groups: the sub-stacks)" % (len(use), nsub))
        p = settings(cfg)
        m = D.measure_all(stack, day_files, freqs, ref["c"], ref["u"], p, n_cores, step.say)
        with h5py.File(step.path("measurements.h5.tmp"), "w") as h:
            for k, v in m.items():
                h[k] = v.astype(np.float32) if k != "reason" else v
            h["freqs"], h["c_ref"], h["u_ref"] = freqs, ref["c"], ref["u"]
            h["sta1"], h["sta2"], h["dist_km"] = np.array(sta1, "S"), np.array(sta2, "S"), dist_all
            h.attrs.update({k: ("none" if v is None else v) for k, v in p.items()})
            h.attrs["reasons"] = str(D.REASONS)
        os.replace(step.path("measurements.h5.tmp"), step.path("measurements.h5"))
        ok = D.accepted(m["reason"], m["c"])
        write_table(step.path("dispersion.csv"), ok, m, freqs, sta1, sta2, coords, dist_all)
        med = write_medians(step.path("median_vs_reference.csv"), ok, m["c"], freqs, ref["c"])
        counts = D.reason_counts(m["reason"][cnt > 0])
        step.say("rejections (values with each reason; a value can have several): %s" % counts)
        step.number("pairs measured", int(len(use)))
        step.number("pairs with a curve", int(ok.any(axis=1).sum()))
        step.number("accepted values", int(ok.sum()))
        step.number("rejections by reason", counts)
        if med is not None:
            step.number("median of (median c / reference - 1), %", round(med, 2))
        from .plotting import disp_figures

        disp_figures.make(cfg, step, ref, m, ok, freqs, dist_all, sta1, sta2)


def _n_substacks(path):
    with h5py.File(path, "r") as h:
        return h["count"].shape[0]


def write_table(path, ok, m, freqs, sta1, sta2, coords, dist):
    with open(path, "w") as fh:
        fh.write("sta1,sta2,lon1,lat1,lon2,lat2,dist_km,frequency_hz,period_s,c_kms,sigma_kms,u_kms,snr\n")
        for q, k in zip(*np.where(ok)):
            sg = m["sigma"][q, k]
            fh.write("%s,%s,%.6f,%.6f,%.6f,%.6f,%.4f,%.6f,%.6f,%.5f,%s,%.5f,%.2f\n" % (
                sta1[q], sta2[q], coords["lon1"][q], coords["lat1"][q], coords["lon2"][q], coords["lat2"][q], dist[q],
                freqs[k], 1 / freqs[k], m["c"][q, k], "" if not np.isfinite(sg) else "%.5f" % sg, m["u"][q, k],
                m["snr"][q, k]))  # fmt: skip


def write_medians(path, ok, c, freqs, c_ref):
    """Per frequency the median of the accepted c against the reference; -> median of the deviations (%) over the
    frequencies with at least 10 pairs (a systematic offset: the reference, and so the cycle count, may be off)."""
    devs = []
    with open(path, "w") as fh:
        fh.write("frequency_hz,pairs,median_c_kms,c_ref_kms,deviation_percent\n")
        for k, f in enumerate(freqs):
            v = c[ok[:, k], k]
            med = float(np.median(v)) if v.size else np.nan
            dev = 100 * (med / c_ref[k] - 1) if v.size and np.isfinite(c_ref[k]) else np.nan
            if v.size >= 10 and np.isfinite(dev):
                devs.append(dev)
            fh.write("%.6f,%d,%s,%s,%s\n" % (f, v.size, reference._v(med), reference._v(c_ref[k]),
                                              "" if not np.isfinite(dev) else "%.2f" % dev))  # fmt: skip
    return float(np.median(devs)) if devs else None
