"""Step 02: hourly cross-coherence of every station pair, one HDF5 file per processing day.

Per station-day (lant.data.read_day): exact UTC sample grid, zero-phase anti-alias low-pass and decimation to
correlation.fs.  Then
  1. optional hourly removal of a UTC-locked 1-s periodic signal (template from folding the trace modulo 1 s;
     the GPS timing pulse of nodal instruments, Takagi et al. 2015 / Wang et al. 2017);
  2. segments of seg_len_s with seg_step_s step; a segment is dropped if it contains a gap or if its peak
     |amplitude| exceeds max_over_std times the standard deviation of the whole station-day (NoisePy
     max_over_std); linear detrend and a 5 % cosine taper per segment;
  3. cross-coherence  C_ij = conj(F_i) F_j / (|F_i| |F_j|)  in fmin-fmax with 0.1 Hz cosine edges
     (Xu et al. 2021, GJI 227:875), no spectral smoothing by default (smooth_bins = 1);
  4. hourly frequency-domain averaging, back to lag time +-maxlag (circular correlation, seg_len >= 2 maxlag).
Sign convention: positive lag = energy travelling from station i to j, (i, j) in sorted order of the names.

Each day file records its settings; a file written with other settings is never reused silently.
"""
import os
import time
from multiprocessing import Pool

import h5py
import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import detrend

from . import data as D
from . import paths

SETTINGS = ("fs", "maxlag", "seg_len", "seg_step", "fmin", "fmax", "decomb", "smooth_bins", "max_over_std", "detrend")


def remove_comb(y, ok, fs):
    """Subtract the UTC-locked 1-s periodic component hour by hour (hours with < 50 % valid samples untouched)."""
    p = int(round(fs))
    ph_all = np.arange(y.size) % p
    for h in range(24):
        a, b = int(h * 3600 * fs), int((h + 1) * 3600 * fs)
        m = ok[a:b]
        if m.mean() < 0.5:
            continue
        seg = y[a:b].astype(np.float64)
        ph = ph_all[a:b]
        fold = np.bincount(ph[m], weights=seg[m], minlength=p)
        cnt = np.bincount(ph[m], minlength=p)
        tmpl = fold / np.maximum(cnt, 1)
        tmpl -= tmpl.mean()
        y[a:b] = (seg - tmpl[ph]).astype(np.float32)
    return y


def settings_of(opt):
    return dict(fs=opt["fs"], maxlag=opt["maxlag"], seg_len=opt["seg_len"], seg_step=opt["seg_step"],
                fmin=opt["fmin"], fmax=opt["fmax"], decomb=bool(opt["decomb"]), smooth_bins=opt["smooth_bins"],
                max_over_std=opt["max_over_std"], detrend="linear")


def process_day(task):
    day, t0, by_station, opt = task
    out_path = os.path.join(opt["out_dir"], "day_%s.h5" % day)
    if os.path.exists(out_path):
        want = settings_of(opt)
        with h5py.File(out_path, "r") as h:
            have = {k: h.attrs.get(k, None) for k in want}
            have_st = sorted(x.decode() for x in h["stations"][:])
        diff = {k: (have[k], v) for k, v in want.items() if have[k] is None or have[k] != v}
        if diff:
            return day, "STALE %s was written with other settings (file, requested): %s" % (out_path, diff)
        return day, "exists (%d stations)" % len(have_st)
    t_start = time.time()
    fs = opt["fs"]
    names, data, oks = [], [], []
    for name in sorted(by_station):
        try:
            y, ok = D.read_day(by_station[name], t0, fs)
        except Exception as exc:
            print("[%s] read failed %s: %s" % (day, name, exc), flush=True)
            continue
        if y is None or ok.mean() < 0.05:
            continue
        if opt["decomb"]:
            y = remove_comb(y, ok, fs)
        names.append(name)
        data.append(y)
        oks.append(ok)
    ns = len(names)
    if ns < 2:
        return day, "too few stations (%d)" % ns
    X = np.stack(data)
    OK = np.stack(oks)
    del data, oks

    n = int(round(opt["seg_len"] * fs))
    step = int(round(opt["seg_step"] * fs))
    L = int(round(opt["maxlag"] * fs))
    nfft = n
    nfreq = nfft // 2 + 1
    freqs = np.arange(nfreq) * fs / nfft
    band = np.where((freqs >= opt["fmin"]) & (freqs <= opt["fmax"]))[0]
    b0, b1 = band[0], band[-1] + 1
    ntap = max(3, int(round(0.1 / (fs / nfft))))  # 0.1 Hz cosine edges
    spec_taper = np.ones(b1 - b0)
    edge = 0.5 * (1 - np.cos(np.pi * np.arange(ntap) / ntap))
    spec_taper[:ntap] = edge
    spec_taper[-ntap:] = edge[::-1]
    win = np.ones(n)
    tl = int(0.05 * n)
    w = 0.5 * (1 - np.cos(np.pi * np.arange(tl) / tl))
    win[:tl] = w
    win[-tl:] = w[::-1]

    iu, ju = np.triu_indices(ns, 1)
    npair = iu.size
    hourly = np.zeros((24, npair, 2 * L + 1), np.float32)
    count = np.zeros((24, npair), np.int32)

    starts = np.arange(0, X.shape[1] - n + 1, step)
    seg_std = np.full((ns, starts.size), np.nan)
    seg_peak = np.full((ns, starts.size), np.nan)
    for k, s0 in enumerate(starts):
        good = OK[:, s0:s0 + n].all(axis=1)
        if good.any():
            seg_std[good, k] = X[good, s0:s0 + n].std(axis=1)
            seg_peak[good, k] = np.abs(X[good, s0:s0 + n]).max(axis=1)
    day_std = np.array([X[i, OK[i]].std() if OK[i].any() else np.nan for i in range(ns)])
    n_gapfree = np.zeros(ns, np.int64)
    n_used = np.zeros(ns, np.int64)

    acc = np.zeros((npair, b1 - b0), np.complex128)
    acc_n = np.zeros(npair, np.int32)
    cur_hour = 0

    def flush(hour):
        if not acc_n.any():
            return
        full = np.zeros((npair, nfreq), np.complex128)
        full[:, b0:b1] = acc / np.maximum(acc_n, 1)[:, None]
        cc = np.fft.irfft(full, n=nfft, axis=1)
        hourly[hour] = np.concatenate([cc[:, -L:], cc[:, :L + 1]], axis=1).astype(np.float32)
        count[hour] = acc_n
        acc[:] = 0
        acc_n[:] = 0

    chunk = opt["pair_chunk"]
    for k, s0 in enumerate(starts):
        hour = int((s0 + n // 2) / fs // 3600)
        if hour != cur_hour:
            flush(cur_hour)
            cur_hour = hour
        good = OK[:, s0:s0 + n].all(axis=1) & (seg_std[:, k] > 0)
        n_gapfree += good
        if opt["max_over_std"] > 0:
            good &= seg_peak[:, k] < opt["max_over_std"] * day_std
        n_used += good
        if good.sum() < 2:
            continue
        seg = detrend(X[:, s0:s0 + n].astype(np.float64), axis=1, type="linear")
        F = np.fft.rfft(seg * win, n=nfft, axis=1)[:, b0:b1]
        A = np.abs(F)
        if opt["smooth_bins"] > 1:
            A = uniform_filter1d(A, size=opt["smooth_bins"], axis=1, mode="nearest")
        F = (F / (A + 1e-30) * spec_taper).astype(np.complex64)
        F[~good] = 0
        pm = good[iu] & good[ju]
        for c0 in range(0, npair, chunk):
            sel = slice(c0, min(c0 + chunk, npair))
            acc[sel] += np.conj(F[iu[sel]]) * F[ju[sel]]
        acc_n += pm
    flush(cur_hour)

    pair_names = ["%s_%s" % (names[i], names[j]) for i, j in zip(iu, ju)]
    tmp = out_path + ".tmp"
    with h5py.File(tmp, "w") as h:
        h.create_dataset("hourly", data=hourly)
        h.create_dataset("count", data=count)
        h.create_dataset("pairs", data=paths.ascii_names(pair_names))
        h.create_dataset("stations", data=paths.ascii_names(names))
        h.create_dataset("segments_gapfree", data=n_gapfree)
        h.create_dataset("segments_used", data=n_used)
        h.attrs.update(settings_of(opt))
        h.attrs.update(dict(t0_utc=t0, method="cross-coherence"))
    os.replace(tmp, out_path)
    return day, "ok %d stations %d pairs %.0f s" % (ns, npair, time.time() - t_start)


def run(cfg):
    c = cfg["correlation"]
    out_dir = paths.mk(paths.P(cfg, paths.CCF))
    stations = D.used_stations(cfg)
    rows, _ = D.station_rows(D.processing_rows(cfg), stations)
    days = D.day_windows([r for v in rows.values() for r in v], cfg)
    if not days:
        raise D.DataError("no processing day holds data of the stations in the table")
    opt = dict(fs=float(c["fs"]), seg_len=float(c["seg_len_s"]), seg_step=float(c["seg_step_s"]),
               maxlag=float(c["maxlag_s"]), fmin=float(c["fmin_hz"]), fmax=float(c["fmax_hz"]),
               smooth_bins=int(c["smooth_bins"]), max_over_std=float(c["max_over_std"]),
               decomb=bool(c["remove_1s_comb"]), out_dir=out_dir, pair_chunk=3000)
    tasks = []
    for lab, (t0, trs) in days.items():
        by = {}
        for r in trs:
            by.setdefault(r["station"], []).append(r)
        tasks.append((lab, t0, by, opt))
    tasks.sort(key=lambda t: -len(t[2]))          # biggest days first so memory peaks early
    os.environ.setdefault("OMP_NUM_THREADS", str(c["threads_per_process"]))
    stale, msgs = [], {}
    with Pool(int(c["nproc"])) as pool:
        for day, msg in pool.imap_unordered(process_day, tasks):
            print("  %s %s" % (day, msg), flush=True)
            msgs[day] = msg
            if msg.startswith("STALE"):
                stale.append(day)
    if stale:
        raise RuntimeError("day files written with other settings in %s (days %s): delete them or use another "
                           "output directory" % (out_dir, ",".join(sorted(stale))))
    return dict(days=len(tasks), messages=msgs)
