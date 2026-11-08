"""Step 04: notch table by rule and the final cross-correlation product.

A frequency is notched when it is
  (a) a persistent narrow line of the RAW waveforms: in 1-h Welch spectra (100 s segments, 0.01 Hz) of
      notch.raw_nsta random stations on notch.raw_ndays days, bins more than 6 dB above a 1-Hz running median
      in more than notch.raw_min_frac of the station-hours; or
  (b) a narrow peak of the pair-averaged amplitude spectrum of the un-notched stack, at least notch.db above its
      1-Hz running median (lines too weak in single records that build up in the stack); applied only where the
      running median has its full window inside the correlation band.
Half-width = half the observed extent + 0.03 Hz, at least 0.05 Hz; entries closer than 0.06 Hz are merged.
The notch is a raised-cosine gain per line (lant.signal.notch_lines), applied once to the total and to every
sub-stack; the dispersion step does not notch again.  notch.table_file replaces the rule by a given table.
"""
import os
from multiprocessing import Pool

import h5py
import numpy as np
from scipy.ndimage import median_filter

from . import data as D
from . import paths
from .signal import notch_table


def _raw_psd(args):
    rows, t0, fmin, fmax = args
    import obspy
    from scipy.signal import welch
    from scipy.ndimage import median_filter as mf
    st = obspy.Stream()
    for path in sorted({r["file"] for r in rows}):
        try:
            st += obspy.read(path, starttime=obspy.UTCDateTime(t0), endtime=obspy.UTCDateTime(t0 + 86400))
        except Exception:
            continue
    st = st.select(network=rows[0]["net"], station=rows[0]["sta"], location=rows[0]["location"],
                   channel=rows[0]["channel"])
    if not len(st):
        return None
    tr = st.merge(fill_value=0)[0]
    fs = tr.stats.sampling_rate
    x = tr.data.astype(np.float64)
    nh = int(3600 * fs)
    out, f = [], None
    for k in range(len(x) // nh):
        y = x[k * nh:(k + 1) * nh]
        if np.std(y) == 0:
            continue
        f, P = welch(y - y.mean(), fs=fs, nperseg=int(100 * fs))
        m = (f >= fmin - 0.5) & (f <= fmax + 0.5)
        lp = np.log10(P[m] + 1e-30)
        out.append(lp - mf(lp, size=101))
        fm = f[m]
    return (fm, np.array(out)) if out else None


def raw_lines(cfg, out_txt):
    n = cfg["notch"]
    fmin, fmax = n["band_hz"]
    stations = D.used_stations(cfg)
    rows, _ = D.station_rows(D.processing_rows(cfg), stations)
    days = D.day_windows([r for v in rows.values() for r in v], cfg)
    labs = sorted(days)
    if n["raw_days"]:
        use_days = [d for d in labs if d in set(map(str, n["raw_days"]))]
    else:
        k = min(n["raw_ndays"], len(labs))
        use_days = [labs[int(round(i))] for i in np.linspace(0, len(labs) - 1, k + 2)[1:-1]] if len(labs) > k else labs
    stas = sorted(rows)
    rng = np.random.default_rng(n["seed"])
    stas = sorted(rng.choice(stas, min(n["raw_nsta"], len(stas)), replace=False).tolist())
    jobs = []
    for s in stas:
        for d in use_days:
            t0 = days[d][0]
            r = [x for x in rows[s] if x["end"] > t0 and x["start"] < t0 + 86400]
            if r:
                jobs.append((r, t0, fmin, fmax))
    with Pool(int(cfg["correlation"]["nproc"])) as p:
        res = [r for r in p.map(_raw_psd, jobs) if r is not None]
    if not res:
        raise RuntimeError("no raw data could be read for the spectral-line search")
    nb = min(len(r[0]) for r in res)
    f = res[0][0][:nb]
    E = np.concatenate([r[1][:, :nb] for r in res])
    inb = (f >= fmin) & (f <= fmax)
    frac = (E > 0.6).mean(0)
    on = (frac > n["raw_min_frac"]) & inb
    L = ["persistent lines on raw data: %d stations x %d days = %d station-hours" % (len(stas), len(use_days), len(E)),
         "stations: " + " ".join(stas), "days: " + " ".join(use_days),
         "bin 'on' = > 6 dB above 1-Hz running median of the 1-h PSD (0.01 Hz bins); line = run of bins on in > %.0f %% "
         "of station-hours" % (100 * n["raw_min_frac"]), "", "%9s %17s %10s" % ("centre", "extent (Hz)", "occupancy")]
    lines = []
    i = 0
    while i < len(f):
        if on[i]:
            j = i
            while j + 1 < len(f) and on[j + 1]:
                j += 1
            k = i + int(np.argmax(frac[i:j + 1]))
            L.append("%9.3f %8.3f-%8.3f %9.0f %%" % (f[k], f[i], f[j], 100 * frac[k]))
            lines.append((f[i], f[j], f[k], 100 * frac[k]))
            i = j + 1
        else:
            i += 1
    open(out_txt, "w").write("\n".join(L) + "\n")
    return lines, dict(stations=stas, days=use_days, station_hours=len(E))


def mean_spectrum(stack_path, chunk=1000):
    """Pair-averaged amplitude spectrum of a stack in dB above its 1-Hz running median."""
    with h5py.File(stack_path, "r") as h:
        fs = float(h.attrs["fs"])
        n = h["stack"].shape[1]
        npair = h["stack"].shape[0]
        nfft = 4 * n
        acc = np.zeros(nfft // 2 + 1)
        for i in range(0, npair, chunk):
            x = h["stack"][i:i + chunk].astype(np.float64)
            x -= x.mean(axis=1, keepdims=True)
            x *= np.hanning(n)
            acc += (np.abs(np.fft.rfft(x, n=nfft, axis=1)) ** 2).sum(axis=0)
    f = np.arange(nfft // 2 + 1) * fs / nfft
    amp = np.sqrt(acc / npair)
    bg = median_filter(amp, size=int(round(1.0 / f[1])) | 1, mode="nearest")
    return f, 20 * np.log10(amp / bg)


def notch_rule(cfg, raw, stack_path, out_table, report):
    n = cfg["notch"]
    fmin, fmax = n["band_hz"]
    ent = [(lo, hi, c, "raw", "%.0f %% of station-hours" % occ) for lo, hi, c, occ in raw if fmin <= c <= fmax]
    f, db = mean_spectrum(stack_path)
    b0, b1 = cfg["correlation"]["fmin_hz"], cfg["correlation"]["fmax_hz"]
    on = (db >= n["db"]) & (f >= max(fmin, b0 + 0.5)) & (f <= min(fmax, b1 - 0.5))
    i = 0
    while i < len(f):
        if on[i]:
            j = i
            while j + 1 < len(f) and on[j + 1]:
                j += 1
            k = i + int(np.argmax(db[i:j + 1]))
            ent.append((f[i], f[j], f[k], "ccf", "%.1f dB" % db[k]))
            i = j + 1
        else:
            i += 1
    ent.sort()
    merged = []
    for e in ent:
        if merged and e[0] - merged[-1][1] <= 0.06:
            lo, hi, c, s, v = merged[-1]
            merged[-1] = (lo, max(hi, e[1]), None, s + "+" + e[3], v + "; " + e[4])
        else:
            merged.append(e)
    rows = []
    for lo, hi, c, s, v in merged:
        centre = 0.5 * (lo + hi) if c is None else c
        hw = max(0.05, max(centre - lo, hi - centre) + 0.03)
        rows.append((centre, np.ceil(hw * 100 - 1e-9) / 100, lo, hi, s, v))
    with open(out_table, "w") as fh:
        fh.write("# Spectral lines to notch (lant notch rule): persistent raw-waveform line, or narrow peak >= %.0f dB of the\n"
                 "# pair-averaged spectrum of the un-notched stack; half-width = half the extent + 0.03 Hz (min 0.05).\n"
                 "# freq_Hz  half_width_Hz\n" % n["db"])
        for r in rows:
            fh.write("%.3f   %.2f\n" % (r[0], r[1]))
    L = ["notch table by rule: %d entries in %.1f-%.1f Hz" % (len(rows), fmin, fmax), "",
         "  freq_Hz  half_width  extent_Hz        source   evidence"]
    for r in rows:
        L.append("  %7.3f   %5.2f     %.3f-%.3f   %-8s %s" % r)
    open(report + ".txt", "w").write("\n".join(L) + "\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(14, 4.8))
    m = (f >= fmin) & (f <= fmax)
    ax.plot(f[m], db[m], "k", lw=0.6)
    ax.axhline(n["db"], color="r", lw=0.8, ls="--")
    for r in rows:
        ax.axvspan(r[0] - r[1], r[0] + r[1], color="C1" if "raw" in r[4] else "C0", alpha=0.35)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("dB above the 1-Hz running median")
    ax.set_title("Pair-averaged spectrum of the un-notched stack and the notches (orange: raw-waveform line; blue: stack only)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(report + ".png", dpi=140)
    plt.close(fig)
    return [(r[0], r[1]) for r in rows]


def apply(stack_path, table, out):
    with h5py.File(stack_path, "r") as h, h5py.File(out + ".tmp", "w") as o:
        fs = float(h.attrs["fs"])
        o["stack"] = notch_table(h["stack"][:], fs, table).astype(np.float32) if len(table) else h["stack"][:]
        nd, npair, nlag = h["daily"].shape
        dd = o.create_dataset("daily", (nd, npair, nlag), np.float32, chunks=(1, npair, nlag))
        for k in range(nd):
            dd[k] = notch_table(h["daily"][k], fs, table).astype(np.float32) if len(table) else h["daily"][k]
        for k in h:
            if k not in ("stack", "daily"):
                o[k] = h[k][:]
        for k, v in h.attrs.items():
            o.attrs[k] = v
        o.attrs["notch_table"] = np.array(table, float).reshape(-1, 2)
    os.replace(out + ".tmp", out)


def run(cfg):
    d = paths.mk(paths.P(cfg, paths.NOTCH))
    clean = paths.P(cfg, paths.STACK, "stack_clean.h5")
    if not os.path.exists(clean):
        raise RuntimeError("%s not found: run `lant run stack` first" % clean)
    n = cfg["notch"]
    table_path = os.path.join(d, "notch_table.txt")
    summary = {}
    if not n["enabled"]:
        table = []
        open(table_path, "w").write("# notch.enabled is false: no line notched\n")
    elif n["table_file"]:
        table = np.atleast_2d(np.loadtxt(n["table_file"])).tolist()
        np.savetxt(table_path, np.array(table), fmt="%.3f   %.2f", header="freq_Hz half_width_Hz (from %s)" % n["table_file"])
    else:
        raw, summary["raw_search"] = raw_lines(cfg, os.path.join(d, "raw_lines.txt"))
        table = notch_rule(cfg, raw, clean, table_path, os.path.join(d, "notch_report"))
    apply(clean, table, paths.stack_product(cfg))
    summary["notch_table"] = table
    return summary


def read_table(cfg):
    p = paths.P(cfg, paths.NOTCH, "notch_table.txt")
    a = np.loadtxt(p, ndmin=2) if os.path.exists(p) and os.path.getsize(p) else np.zeros((0, 2))
    return a.reshape(-1, 2)
