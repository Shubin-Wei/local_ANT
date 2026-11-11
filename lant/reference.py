"""Step 04: reference phase-velocity curve, and the measurement band and velocity axis that follow from it.

The reference guides the group ridge, fixes the cycle count of the phase and centres the deviation screens.
  reference.source = fj (default): F-J spectrum of this array and its fundamental-mode ridge (lant/fj.py);
  reference.source = file: a curve of another array (npz written by lant, or text "freq_Hz c_km/s"), interpolated
      in log f onto the frequencies only inside its own valid segments (borrowed_curve): no value across a stretch
      the source marks invalid and none outside (no extrapolation).
From the curve (lant.auto.from_reference): dispersion.fmin_hz (the longest pair = qc.lam_c wavelengths),
dispersion.vmin / vmax / dv, export.tmax_s; the reference is then cut to the measurement band.  The step stops
for `lant confirm reference`: are the picks on the fundamental-mode ridge, are the band and the velocities right?
"""

import os

import h5py
import numpy as np

from . import fj as FJ
from . import paths
from .signal import log_freqs


def freq_grid(cfg):
    """The measurement frequencies (dispersion.fmin_hz - fmax_hz on the 1/per_octave grid through 1 Hz)."""
    d = cfg["dispersion"]
    f = log_freqs(d["fmin_hz"], d["fmax_hz"], d["per_octave"])
    return f[f <= d["fmax_hz"] + 1e-9]


def borrowed_curve(path, freqs):
    """A borrowed reference curve on `freqs`, keeping the validity its source gives.

    npz written by lant: the points marked valid (ref_valid) form contiguous segments in the source's own frequency
    order; a target frequency gets a value only inside a segment (interpolated in log f between the segment's
    points), never across an invalid stretch and never outside (no extrapolation).  Text "freq_Hz c_km/s": the file
    says nothing about validity, so its listed points are taken as one segment and labelled "validity unknown".
    Support (the curve is backed by measurements there; only supported frequencies anchor the cycle count and get the
    deviation screens): an npz written by lant carries ref_support, and a target frequency is supported only when the
    two source points around it are both supported (a held / extrapolated stretch of the source stays unsupported);
    a text curve is taken as supported where it has values (it is given as a measured or model curve).
    -> c (NaN where there is no value), valid mask, support mask, origin labels per frequency, note."""
    freqs = np.asarray(freqs, float)
    if path.endswith(".npz"):
        z = np.load(path)
        f, c = np.asarray(z["freqs"], float), np.asarray(z["c_ref"], float)
        v = np.asarray(z["ref_valid"], bool) if "ref_valid" in z.files else np.ones(len(f), bool)
        sp = np.asarray(z["ref_support"], bool) if "ref_support" in z.files else v.copy()
        known = "ref_valid" in z.files
    else:
        a = np.loadtxt(path, ndmin=2)
        f, c = a[:, 0], a[:, 1]
        v = np.ones(len(f), bool)
        sp = v.copy()
        known = False
    o = np.argsort(f)
    f, c, v, sp = f[o], c[o], v[o], sp[o]
    v &= np.isfinite(c) & (c > 0) & (f > 0)
    if v.sum() < 2:
        raise RuntimeError("%s: fewer than 2 valid points (first column frequency in Hz, second c in km/s)" % path)
    if np.any(np.diff(f) <= 0):
        raise RuntimeError("%s: repeated frequencies" % path)
    if c[v].min() < 0.05 or c[v].max() > 10:
        raise RuntimeError("%s: velocities %.3f-%.3f look wrong (km/s expected)" % (path, c[v].min(), c[v].max()))
    out = np.full(len(freqs), np.nan)
    sup = np.zeros(len(freqs), bool)
    origin = np.array(["none"] * len(freqs), dtype=object)
    k = 0
    while k < len(f):
        if not v[k]:
            k += 1
            continue
        j = k
        while j + 1 < len(f) and v[j + 1]:
            j += 1
        if j > k:
            m = (freqs >= f[k] * (1 - 1e-6)) & (freqs <= f[j] * (1 + 1e-6))
            out[m] = np.exp(np.interp(np.log(freqs[m]), np.log(f[k : j + 1]), np.log(c[k : j + 1])))
            origin[m] = "borrowed" if known else "borrowed, validity unknown"
            for i in np.where(m)[0]:  # supported only between two supported source points
                b = int(np.clip(np.searchsorted(f[k : j + 1], freqs[i] * (1 - 1e-9)), 0, j - k))
                a = max(0, b - 1)
                at = abs(np.log(freqs[i] / f[k + b])) < 1e-6
                sup[i] = sp[k + b] if at else (sp[k + a] and sp[k + b])
        k = j + 1
    note = (
        ("valid segments of %s" % path)
        if known
        else (
            "%s has no validity information: its points are taken as one "
            "segment; gaps of the source cannot be recognised" % path
        )
    )
    valid = np.isfinite(out)
    sup &= valid
    for i in np.where(sup)[0]:
        origin[i] = origin[i] + ", measured"
    return out, valid, sup, origin, note


def from_file(cfg, path, freqs=None):
    freqs = FJ.scan_grid(cfg) if freqs is None else freqs
    c, valid, support, origin, note = borrowed_curve(path, freqs)
    return freqs, c, valid, support, origin, note


def _band_cut(cfg, freqs, *arrays):
    """The arrays on the measurement grid (dispersion.fmin_hz - fmax_hz, a subset of the scan grid)."""
    f = freq_grid(cfg)
    idx = np.array([int(np.argmin(np.abs(np.log(freqs / x)))) for x in f])
    if not np.allclose(freqs[idx], f, rtol=1e-6):
        raise RuntimeError(
            "dispersion.fmin_hz / fmax_hz %g-%g are not on the 1/%d-octave grid through 1 Hz"
            % (f[0], f[-1], cfg["dispersion"]["per_octave"])
        )
    return (f,) + tuple(np.asarray(a)[idx] for a in arrays)


def zone_references(cfg, stack, d, freqs, c, curve, support):
    """reference.zone_km: one F-J curve per square zone of that size, from the pairs whose path midpoint lies in
    the zone (local references where the structure changes across the array; a choice, not a literature value); every
    pair is then measured against the curve of its zone.  A zone whose image does not give a curve (too few
    pairs, no clear ridge) uses the curve of the whole array.  Writes work/04_reference/zones.npz, zones.txt, zones.png."""
    from .signal import local_xy_km

    r = cfg["reference"]
    size = float(r["zone_km"])
    with h5py.File(stack, "r") as h:
        mlon, mlat = 0.5 * (h["lon1"][:] + h["lon2"][:]), 0.5 * (h["lat1"][:] + h["lat2"][:])
        has = h["nseg"][:] > 0
    x, y = local_xy_km(mlon, mlat)
    ix, iy = np.floor((x - x[has].min()) / size).astype(int), np.floor((y - y[has].min()) / size).astype(int)
    keys = sorted({(a, b) for a, b in zip(ix[has], iy[has])})
    zone = np.full(len(has), -1, int)
    curves, sups, L = (
        [],
        [],
        ["zone references (reference.zone_km = %g): zone, centre x y (km), pairs, picks used, curve" % size],
    )
    for a, b in keys:
        sel = has & (ix == a) & (iy == b)
        cx, cy = float(x[sel].mean()), float(y[sel].mean())
        try:
            f2, c2, img, npair = FJ.fj_image(cfg, stack, freqs, select=sel, c_limit=(freqs, curve))
            if r["fj_mode"] == "max_peak":
                _, used, cz, _ = FJ.pick_max(cfg, f2, c2, img, npair)
            else:
                _, used, cz, _ = FJ.pick_track(cfg, f2, c2, img, start_band=FJ.track_start(cfg, f2, c2, img, npair))
            fu = freqs[used]
            sz = (freqs >= fu.min() * (1 - 1e-9)) & (freqs <= fu.max() * (1 + 1e-9))
            zone[sel] = len(curves)
            curves.append(cz)
            sups.append(sz)
            L.append(
                "  %d  %+.2f %+.2f  %d pairs  %d picks  own curve (%.3f-%.3f km/s)"
                % (len(curves) - 1, cx, cy, sel.sum(), used.sum(), cz.min(), cz.max())
            )
        except RuntimeError as exc:
            L.append("  -  %+.2f %+.2f  %d pairs  whole-array curve (%s)" % (cx, cy, sel.sum(), str(exc)[:70]))
    FJ.fj_image(cfg, stack, freqs[:1])  # restores the pair threshold of the whole array
    f, *cut = _band_cut(cfg, freqs, *curves, *sups)
    n = len(curves)
    np.savez(
        os.path.join(d, "zones.npz"),
        freqs=f,
        zone=zone,
        c_ref=np.array(cut[:n]).reshape(n, len(f)),
        ref_support=np.array(cut[n:]).reshape(n, len(f)).astype(bool),
    )
    open(os.path.join(d, "zones.txt"), "w").write("\n".join(L) + "\n")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for k, cz in enumerate(curves):
        ax.plot(freqs, cz, lw=1.0, label="zone %d" % k if n <= 12 else None)
    ax.plot(freqs, curve, "k-", lw=2.2, label="whole array")
    ax.set_xscale("log")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Phase velocity (km/s)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(d, "zones.png"), dpi=150)
    plt.close(fig)
    return dict(
        zone_km=size, zones=len(keys), with_own_curve=n, pairs_on_whole_array_curve=int((has & (zone < 0)).sum())
    )


def run(cfg):
    from . import auto, disp

    d = paths.mk(paths.P(cfg, paths.REFERENCE))
    r = cfg["reference"]
    geom = auto.load(cfg).get("geometry")
    if not geom:
        raise RuntimeError(".lant/auto_parameters.json has no geometry: run `lant run clean` first")
    auto.reset(cfg, auto.REFERENCE_KEYS)
    cfg = auto.restore(cfg, auto.REFERENCE_KEYS)
    if r["source"] == "file":
        freqs, c, valid, support, origin, note = from_file(cfg, r["file"])
        if valid.sum() < 2:
            raise RuntimeError(
                "the borrowed curve %s covers fewer than 2 of the frequencies %.3g-%.3g Hz"
                % (r["file"], freqs[0], freqs[-1])
            )
        idx = np.where(valid)[0]
        fill = c[
            idx[np.argmin(np.abs(np.arange(len(freqs))[:, None] - idx[None, :]), axis=1)]
        ]  # held for the band rule
        _, u, _, _ = disp.reference_work(freqs, c, valid.copy())
        u = np.where(np.isfinite(u), u, fill)
        picked = (float(freqs[valid].min()), float(freqs[valid].max()))
        dec, warn = auto.from_reference(cfg, freqs, fill, u, picked, geom)
        auto.write(cfg, dec, None)
        cfg = auto.apply(cfg)
        f, c2, v2, s2, o2 = _band_cut(cfg, freqs, c, valid, support, origin)
        if not v2.any():
            raise RuntimeError("the borrowed curve %s covers none of the measurement frequencies" % r["file"])
        np.savez(
            os.path.join(d, "reference.npz"),
            freqs=f,
            c_ref=c2,
            ref_valid=v2,
            ref_support=s2,
            origin=o2.astype(str),
            c_ref_source=np.array("borrowed: %s" % r["file"]),
        )
        return dict(
            source="file",
            file=r["file"],
            note=note,
            defined_hz=[float(f[v2].min()), float(f[v2].max())],
            frequencies_without_reference=int((~v2).sum()),
            warnings=warn,
            band_hz=[cfg["dispersion"]["fmin_hz"], cfg["dispersion"]["fmax_hz"]],
            velocity_axis_kms=[cfg["dispersion"]["vmin"], cfg["dispersion"]["vmax"], cfg["dispersion"]["dv"]],
        )
    stack = paths.stack_product(cfg)
    lam_alias = FJ.alias_wavelength_km(cfg, geom)

    def pick(c_limit):
        freqs, c, img, npair = FJ.fj_image(cfg, stack, c_limit=c_limit)
        start = None
        if r["fj_mode"] == "max_peak":
            picks, used, curve, flat_note = FJ.pick_max(cfg, freqs, c, img, npair)
        else:
            start = FJ.track_start(cfg, freqs, c, img, npair, lam_min_km=lam_alias)
            picks, used, curve, flat_note = FJ.pick_track(cfg, freqs, c, img, start_band=start, lam_min_km=lam_alias)
        return freqs, c, img, npair, start, picks, used, curve, flat_note

    c_limit, limit_text = None, "%.2f km/s" % float(r["fj_c_nominal"])
    if auto.user_value(cfg, ("reference", "fj_c_nominal")) is None:
        # the velocity is not known before the spectrum: a first pass with every pair gives a curve, and the
        # second pass limits the pairs of each frequency to fj_rmax_lambda wavelengths of that curve
        freqs, curve = pick((np.array([1e-6, 1e6]), np.array([1e9, 1e9])))[0::7][:2]
        c_limit, limit_text = (freqs, curve), "the curve of a first pass with every pair"
    freqs, c, img, npair, start, picks, used, curve, flat_note = pick(c_limit)
    np.savez(os.path.join(d, "fj.npz"), freqs=freqs, c=c, img=img, npair=npair, picks=picks, used=used)
    fu = freqs[used]
    picked = (float(fu.min()), float(fu.max()))
    u = disp.group_reference(freqs, curve)
    u = np.where(np.isfinite(u) & (u > 0), u, curve)
    dec, warn = auto.from_reference(cfg, freqs, curve, u, picked, geom)
    auto.write(cfg, dec, None, extra={"reference_picks": dict(picked_hz=list(picked), track_start_band_hz=start)})
    cfg = auto.apply(cfg)
    support = (freqs >= picked[0] * (1 - 1e-9)) & (freqs <= picked[1] * (1 + 1e-9))  # the picks back these; held beyond
    f, cb, sb = _band_cut(cfg, freqs, curve, support)
    np.savez(
        os.path.join(d, "reference.npz"),
        freqs=f,
        c_ref=cb,
        ref_valid=np.ones(len(f), bool),
        ref_support=sb,
        c_ref_source=np.array("F-J (%s), rmax %.1f lambda" % (r["fj_mode"], r["fj_rmax_lambda"])),
    )
    L = [
        "F-J reference (%s), pairs within %.1f wavelengths of %s; picks used %.3g-%.3g Hz; measurement band "
        "%.3g-%.3g Hz"
        % (
            r["fj_mode"],
            r["fj_rmax_lambda"],
            limit_text,
            picked[0],
            picked[1],
            cfg["dispersion"]["fmin_hz"],
            cfg["dispersion"]["fmax_hz"],
        )
        + (
            "; picks only where the wavelength >= %.3f km (twice the mean station spacing: spatial Nyquist "
            "limit of the line)" % lam_alias
            if lam_alias
            else ""
        ),
        "  f(Hz)  pairs   pick    used  curve",
    ]
    for i in range(len(freqs)):
        L.append(
            "  %5.2f  %5d  %s  %s  %6.3f"
            % (
                freqs[i],
                npair[i],
                "%6.3f" % picks[i] if np.isfinite(picks[i]) else "   -  ",
                "y" if used[i] else "-",
                curve[i],
            )
        )
    L += [""] + ["  WARNING " + w for w in warn]
    open(os.path.join(d, "fj.txt"), "w").write("\n".join(L) + "\n")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.5))
    Y = img / (img.max(axis=1, keepdims=True) + 1e-30)
    ax.pcolormesh(freqs, c, Y.T, cmap="jet", vmin=0, vmax=1, shading="nearest")
    ax.plot(freqs[used], picks[used], "o", ms=3.5, mfc="w", mec="k", mew=0.6, label="picks used")
    ax.plot(freqs, curve, "w-", lw=1.8, label="reference (held outside the picks)")
    for fb in (cfg["dispersion"]["fmin_hz"], cfg["dispersion"]["fmax_hz"]):
        ax.axvline(fb, color="m", lw=1.0)
    for v in cfg["data"]["velocity_range_kms"]:
        ax.axhline(v, color="w", lw=0.8, ls=":")
    if lam_alias:
        ax.plot(freqs, lam_alias * freqs, "w--", lw=1.0, label="wavelength = 2 x mean station spacing (picks above)")
        ax.set_ylim(c[0], c[-1])
    ax.set_xscale("log")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Phase velocity (km/s)")
    ax.set_title(
        "F-J spectrum, normalised per frequency; magenta: measurement band; dotted: your velocity range", fontsize=9
    )
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.join(d, "fj.png"), dpi=150)
    plt.close(fig)
    zp = os.path.join(d, "zones.npz")
    if os.path.exists(zp):
        os.remove(zp)
    zsum = zone_references(cfg, stack, d, freqs, c, curve, support) if r["zone_km"] else None
    out = dict(
        source="fj",
        zones=zsum,
        mode=r["fj_mode"],
        picks_used=int(used.sum()),
        picked_hz=list(picked),
        band_hz=[cfg["dispersion"]["fmin_hz"], cfg["dispersion"]["fmax_hz"]],
        velocity_axis_kms=[cfg["dispersion"]["vmin"], cfg["dispersion"]["vmax"], cfg["dispersion"]["dv"]],
        pairs_per_frequency=[int(n) for n in npair],
        warnings=warn,
        alias_wavelength_km=lam_alias,
    )
    for w in warn:
        print("[lant] reference WARNING: " + w, flush=True)
    if flat_note:
        out["held_flat"] = flat_note
        print("[lant] reference: " + flat_note, flush=True)
    return out
