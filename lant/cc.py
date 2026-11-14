"""Step 02: hourly cross-coherence of every station pair, one HDF5 file per processing day.

Per station-day (lant.data.read_day): exact UTC sample grid, zero-phase anti-alias low-pass and decimation to
correlation.fs.  Then
  1. hourly removal of a UTC-locked 1-s periodic signal (template from folding the trace modulo 1 s; the timing
     pulse of some loggers, Takagi et al. 2015) for the stations where the clean step detected it
     (correlation.remove_1s_comb = auto; true: every station; false: none);
  2. segments of seg_len_s with seg_step_s step; a segment is dropped if it contains a gap or if its peak
     |amplitude| exceeds max_over_std times the standard deviation of the whole station-day (NoisePy
     max_over_std); linear detrend and a 5 % cosine taper per segment;
  3. cross-coherence  C_ij = conj(F_i) F_j / (|F_i| |F_j|)  in fmin-fmax with quarter-octave cosine edges
     (Xu et al. 2021, GJI 227:875), no spectral smoothing by default (smooth_bins = 1);
  4. hourly frequency-domain averaging, back to lag time +-maxlag (circular correlation, seg_len >= 2 maxlag).
Sign convention: positive lag = energy travelling from station i to j, (i, j) in sorted order of the names.

Each day file records its settings and a fingerprint of its input files (path, size, modification time and time
span of every record of the day).  An existing day file is reused only when both agree with the present run;
otherwise the day is computed again and the file replaced.  (A file whose content changed while its size and
modification time stayed the same is not noticed.)
"""

import os
import time
from multiprocessing import Pool

import h5py
import numpy as np
from scipy.ndimage import uniform_filter1d

from . import cc_kernel as K
from . import data as D
from . import paths
from . import pulse as PU
from .processes import pool_results
from .progress import track

SETTINGS = (
    "fs",
    "maxlag",
    "seg_len",
    "seg_step",
    "fmin",
    "fmax",
    "max_dist",
    "decomb",
    "smooth_bins",
    "max_over_std",
    "detrend",
    "band_taper",
)
READ_THREADS = 4  # station-days read and decimated at the same time in one day process (memory: resources())


def settings_of(opt):
    return dict(
        fs=opt["fs"],
        maxlag=opt["maxlag"],
        seg_len=opt["seg_len"],
        seg_step=opt["seg_step"],
        fmin=opt["fmin"],
        fmax=opt["fmax"],
        max_dist=float(opt.get("max_dist") or 0),
        decomb="all" if opt["decomb"] is True else ",".join(sorted(opt["decomb"] or [])) or "none",
        smooth_bins=opt["smooth_bins"],
        max_over_std=opt["max_over_std"],
        detrend="linear",
        band_taper="quarter-octave cosine",
    )


def band_taper(freqs, b0, b1):
    """Cosine edges over a quarter octave inside each end of the band (relative, so the same for every band)."""
    f = freqs[b0:b1]
    lo, hi = f[0], f[-1]
    t = np.ones(f.size)
    r = np.log2(f / lo) / 0.25
    t = np.where(r < 1, 0.5 * (1 - np.cos(np.pi * np.clip(r, 0, 1))), t)
    r = np.log2(hi / f) / 0.25
    t = np.where(r < 1, t * 0.5 * (1 - np.cos(np.pi * np.clip(r, 0, 1))), t)
    return t


def day_ccf(by_station, t0, opt, nbin=24, failures=None, out=None):
    """Cross-coherence of every pair of one processing day, averaged in nbin equal time bins (24 = hours).
    out(names, shape): where the bins are written as they are done (e.g. an HDF5 dataset; default an array in memory)
    -> (names, ccf [nbin, npair, 2L+1], count [nbin, npair], segments gap-free, segments used) or a message."""
    fs = opt["fs"]
    sel = opt.get("stations")
    todo = [name for name in sorted(by_station) if sel is None or name in sel]
    nout = int(round(86400 * fs))
    # filled in station order as the readers finish (no list of days to stack: one copy of the day in memory)
    X = np.empty((len(todo), nout), np.float32)
    OK = np.empty((len(todo), nout), bool)
    names = []
    for name, got in zip(todo, read_days(by_station, todo, t0, opt)):
        if isinstance(got, Exception):
            print("[%s] read failed %s: %s" % (D.utc_window(t0), name, got), flush=True)
            if failures is not None:
                failures.append("%s: %s" % (name, got))
            continue
        if got is None:
            continue
        X[len(names)], OK[len(names)] = got
        names.append(name)
    ns = len(names)
    if ns < 2:
        return "too few stations (%d)" % ns
    X, OK = X[:ns], OK[:ns]

    n = int(round(opt["seg_len"] * fs))
    step = int(round(opt["seg_step"] * fs))
    L = int(round(opt["maxlag"] * fs))
    nfft = n
    nfreq = nfft // 2 + 1
    freqs = np.arange(nfreq) * fs / nfft
    band = np.where((freqs >= opt["fmin"]) & (freqs <= opt["fmax"]))[0]
    b0, b1 = band[0], band[-1] + 1
    spec_taper = band_taper(freqs, b0, b1)
    win = np.ones(n)
    tl = int(0.05 * n)
    w = 0.5 * (1 - np.cos(np.pi * np.arange(tl) / tl))
    win[:tl] = w
    win[-tl:] = w[::-1]

    iu, ju = pair_indices(names, opt.get("coords"), opt.get("max_dist"))
    npair = iu.size
    if not npair:
        return "too few stations (no pair within correlation.max_distance_km)"
    nlag = 2 * L + 1
    hourly = out(names, (nbin, npair, nlag)) if out is not None else np.zeros((nbin, npair, nlag), np.float32)
    count = np.zeros((nbin, npair), np.int32)
    bin_len = 86400.0 / nbin

    starts = np.arange(0, X.shape[1] - n + 1, step)
    seg_std = np.full((ns, starts.size), np.nan)
    seg_peak = np.full((ns, starts.size), np.nan)
    for k, s0 in enumerate(starts):
        good = OK[:, s0 : s0 + n].all(axis=1)
        if good.any():
            seg_std[good, k] = X[good, s0 : s0 + n].std(axis=1)
            seg_peak[good, k] = np.abs(X[good, s0 : s0 + n]).max(axis=1)
    day_std = np.array([X[i, OK[i]].std() if OK[i].any() else np.nan for i in range(ns)])
    n_gapfree = np.zeros(ns, np.int64)
    n_used = np.zeros(ns, np.int64)

    # the spectra of the segments of one time bin, then their pair products summed in segment order (cc_kernel)
    seg_bin = [min(nbin - 1, int((s0 + n // 2) / fs // bin_len)) for s0 in starts]
    nb = b1 - b0
    spec = np.zeros((max(np.bincount(seg_bin)) if seg_bin else 1, ns, nb), np.complex64)
    spec_good = np.zeros((spec.shape[0], ns), bool)
    acc = np.empty((npair, nb), np.complex128)
    acc_n = np.zeros(npair, np.int32)
    kept = 0
    cur_hour = 0
    chunk = opt["pair_chunk"]
    threads = max(1, int(opt.get("threads") or 1))

    def flush(hour):
        nonlocal kept
        if acc_n.any():
            K.pair_sums(spec, spec_good, kept, iu, ju, acc, threads)
            for c0 in range(0, npair, chunk):  # back to lag time, a block of pairs at a time
                c1 = min(c0 + chunk, npair)
                full = np.zeros((c1 - c0, nfreq), np.complex128)
                full[:, b0:b1] = acc[c0:c1] / np.maximum(acc_n[c0:c1], 1)[:, None]
                cc = np.fft.irfft(full, n=nfft, axis=1)
                hourly[hour, c0:c1] = np.concatenate([cc[:, -L:], cc[:, : L + 1]], axis=1).astype(np.float32)
            count[hour] = acc_n
        kept = 0
        acc_n[:] = 0

    def segment_good(k):
        good = OK[:, starts[k] : starts[k] + n].all(axis=1) & (seg_std[:, k] > 0)
        gapfree = good.copy()
        if opt["max_over_std"] > 0:
            good &= seg_peak[:, k] < opt["max_over_std"] * day_std
        return gapfree, good

    def segment_spectra(k):
        """Whitened band spectra of segment k (every station; the bad ones are zeroed by the caller)."""
        s0 = starts[k]
        seg = K.detrend_linear(X[:, s0 : s0 + n].astype(np.float64))
        F = np.fft.rfft(seg * win, n=nfft, axis=1)[:, b0:b1]
        A = np.abs(F)
        if opt["smooth_bins"] > 1:
            A = uniform_filter1d(A, size=opt["smooth_bins"], axis=1, mode="nearest")
        return (F / (A + 1e-30) * spec_taper).astype(np.complex64)

    # The segments of one time bin are independent until their spectra are summed: their spectra are computed in
    # parallel (threads), then taken in segment order exactly as before.
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(threads) as pool:
        k = 0
        while k < len(starts):
            hour = seg_bin[k]
            if hour != cur_hour:
                flush(cur_hour)
                cur_hour = hour
            ks = [j for j in range(k, len(starts)) if seg_bin[j] == hour]
            k = ks[-1] + 1
            goods = [segment_good(j) for j in ks]
            use = [j for j, (_, g) in zip(ks, goods) if g.sum() >= 2]
            spectra = dict(zip(use, pool.map(segment_spectra, use)))
            for j, (gapfree, good) in zip(ks, goods):
                n_gapfree += gapfree
                n_used += good
                if good.sum() < 2:
                    continue
                F = spectra.pop(j)
                F[~good] = 0
                spec[kept] = F
                spec_good[kept] = good
                kept += 1
                acc_n += good[iu] & good[ju]
    flush(cur_hour)
    return names, hourly, count, n_gapfree, n_used


def read_days(by_station, todo, t0, opt):
    """The decimated day of every station in `todo` (lant.data.read_day, then the 1-s comb removal where asked), read
    by up to READ_THREADS threads (the decoding and filtering release the interpreter lock) -> in the order of
    `todo`: (data, valid mask), None (no data) or the exception raised."""
    from concurrent.futures import ThreadPoolExecutor

    fs = opt["fs"]

    def one(name):
        try:
            y, ok = D.read_day(by_station[name], t0, fs)
            if y is None or not ok.any():
                return None
            if opt["decomb"] is True or (opt["decomb"] and name in opt["decomb"]):
                y = PU.remove_comb(y, ok, fs)
            return y, ok
        except Exception as exc:
            return exc

    n = max(1, min(READ_THREADS, int(opt.get("threads") or 1)))
    if n == 1:
        yield from map(one, todo)
        return
    with ThreadPoolExecutor(n) as ex:
        pending = [ex.submit(one, name) for name in todo[:n]]
        for k in range(len(todo)):  # at most n days in flight: memory does not grow with the stations
            got = pending[k].result()
            pending[k] = None
            if k + n < len(todo):
                pending.append(ex.submit(one, todo[k + n]))
            yield got


def pair_indices(names, coords=None, max_dist=None):
    """Index arrays (i < j) of the pairs of `names` (sorted); with max_dist only those not farther apart (km)."""
    iu, ju = np.triu_indices(len(names), 1)
    if max_dist and coords:
        from .signal import geodesic_km

        lon = np.array([coords[n][0] for n in names])
        lat = np.array([coords[n][1] for n in names])
        keep = geodesic_km(lon[iu], lat[iu], lon[ju], lat[ju])[0] <= float(max_dist)
        iu, ju = iu[keep], ju[keep]
    return iu, ju


def input_fingerprint(by_station, sel=None):
    """sha1 over the records of one day: station, file, size, modification time, time span, polarity."""
    import hashlib

    h = hashlib.sha1()
    for name in sorted(by_station):
        if sel is not None and name not in sel:
            continue
        for r in sorted(by_station[name], key=lambda r: (r["file"], r["start"])):
            try:
                st = os.stat(r["file"])
                sig = (st.st_size, st.st_mtime_ns)
            except OSError:
                sig = ("missing",)
            h.update(repr((name, r["file"], sig, r["start"], r["end"], r.get("fs"), r.get("polarity", 1))).encode())
    return h.hexdigest()


def process_day(task):
    day, t0, by_station, opt = task
    out_path = os.path.join(opt["out_dir"], "day_%s.h5" % day)
    inputs = input_fingerprint(by_station, opt.get("stations"))
    redo = ""
    if os.path.exists(out_path):
        want = settings_of(opt)
        try:
            with h5py.File(out_path, "r") as h:
                have = {k: h.attrs.get(k, None) for k in want}
                nst = len(h["stations"])
                same = (
                    float(h.attrs.get("t0_utc", np.nan)) == t0
                    and h.attrs.get("inputs", "") == inputs
                    and all(have[k] is not None and have[k] == v for k, v in want.items())
                )
        except OSError:
            same = False  # an unreadable (half-written) file is computed again
        if same:
            return day, "exists (%d stations)" % nst, []
        redo = "recomputed (settings or input files changed): "
    t_start = time.time()
    failures = []
    tmp = out_path + ".tmp"
    files = []

    def out(names, shape):  # the hourly array goes to the file hour by hour instead of staying in memory
        files.append(h5py.File(tmp, "w"))
        return files[0].create_dataset("hourly", shape, np.float32)

    try:
        res = day_ccf(by_station, t0, opt, 24, failures, out=out)
        if isinstance(res, str):
            if files:
                files[0].close()
                os.remove(tmp)
            if os.path.exists(out_path):  # a file of an earlier run no longer stands for this day
                os.remove(out_path)
            return day, res, failures
        names, _, count, n_gapfree, n_used = res
        ns = len(names)
        iu, ju = pair_indices(names, opt.get("coords"), opt.get("max_dist"))
        npair = iu.size
        pair_names = ["%s_%s" % (names[i], names[j]) for i, j in zip(iu, ju)]
        h = files[0]
        h.create_dataset("count", data=count)
        h.create_dataset("pairs", data=paths.ascii_names(pair_names))
        h.create_dataset("stations", data=paths.ascii_names(names))
        h.create_dataset("segments_gapfree", data=n_gapfree)
        h.create_dataset("segments_used", data=n_used)
        h.attrs.update(settings_of(opt))
        h.attrs.update(dict(t0_utc=t0, method="cross-coherence", inputs=inputs))
    except BaseException:
        for h in files:
            h.close()
        if files and os.path.exists(tmp):
            os.remove(tmp)
        raise
    finally:
        for h in files:
            h.close()
    os.replace(tmp, out_path)
    return day, "%sok %d stations %d pairs %.0f s" % (redo, ns, npair, time.time() - t_start), failures


def resources(cfg, nsta, npair, ndays, raw_fs=None):
    """What the cc step needs -> dict(memory_per_process_gb, disk_gb, ram_gb, lags, segment).  Memory of one process
    (one day): the day of every station (samples in single precision and their mask), and the larger of the reading
    (READ_THREADS station-days at the raw rate, about 16 bytes per raw sample while decoded and filtered) and the
    correlation (the pair sums of one hour in double precision, the spectra of its segments in single precision and
    the transform of a block of pairs); the hourly array is written to the file hour by hour.  Disk: the hourly files
    of every day."""
    c = cfg["correlation"]
    fs = float(c["fs"])
    nlag = 2 * float(c["maxlag_s"]) * fs + 1
    hourly = 24.0 * npair * nlag * 4
    day = nsta * 86400.0 * fs * 5
    nseg = float(c["seg_len_s"]) * fs
    nband = nseg / 2 + 1  # at most: the band is narrower
    per_hour = 3600.0 / float(c["seg_step_s"]) + 2
    reading = min(READ_THREADS, nsta) * 16 * 86400.0 * float(raw_fs or fs)
    correlating = npair * nband * 16 + per_hour * nsta * nband * 8 + 3 * min(3000, npair) * nseg * 8
    from .resources import available_memory

    memory = available_memory()
    ram = memory / 1e9 if memory is not None else None
    return dict(
        memory_per_process_gb=(day + max(reading, correlating)) / 1e9,
        disk_gb=ndays * hourly / 1e9,
        ram_gb=ram,
        maxlag_s=float(c["maxlag_s"]),
        seg_len_s=float(c["seg_len_s"]),
        pairs=int(npair),
    )


def resources_text(cfg, r):
    v = cfg["data"]["velocity_range_kms"]
    return (
        "correlations: %d pairs, lags +-%g s, segments %g s: about %.1f GB of memory per process "
        "(one day%s) and %.0f GB of disk for the hourly files.  Both grow with the lag window, "
        "which is set by the slowest velocity of data.velocity_range_kms (%g km/s now): a range that fits the "
        "area shortens it"
        % (
            r["pairs"],
            r["maxlag_s"],
            r["seg_len_s"],
            r["memory_per_process_gb"],
            "; %.0f GB effectively available" % r["ram_gb"] if r["ram_gb"] else "",
            r["disk_gb"],
            v[0],
        )
    )


def run(cfg):
    c = cfg["correlation"]
    out_dir = paths.mk(paths.P(cfg, paths.CCF))
    stations = D.used_stations(cfg)
    rows, _ = D.station_rows(D.processing_rows(cfg), stations)
    days = D.day_windows([r for v in rows.values() for r in v], cfg)
    if not days:
        raise D.DataError("no processing day holds data of the stations in the table")
    decomb = c["remove_1s_comb"]
    if decomb == "auto":  # the stations where the clean step detected the signal
        import json

        idx = json.load(open(paths.P(cfg, "01_clean", "archive_index.json")))
        decomb = sorted(set(idx.get("pulse_1s") or []) & set(stations))
    if decomb and abs(float(c["fs"]) - round(float(c["fs"]))) > 1e-6:
        raise D.DataError(
            "correlation.remove_1s_comb needs a whole number of samples per second; correlation.fs is "
            "%g Hz: set remove_1s_comb to false or choose a whole-number rate" % float(c["fs"])
        )
    opt = dict(
        fs=float(c["fs"]),
        seg_len=float(c["seg_len_s"]),
        seg_step=float(c["seg_step_s"]),
        maxlag=float(c["maxlag_s"]),
        fmin=float(c["fmin_hz"]),
        fmax=float(c["fmax_hz"]),
        smooth_bins=int(c["smooth_bins"]),
        max_over_std=float(c["max_over_std"]),
        decomb=decomb,
        out_dir=out_dir,
        pair_chunk=3000,
        max_dist=c["max_distance_km"],
        coords={k: (v["lon"], v["lat"]) for k, v in stations.items()},
    )
    npair = pair_indices(sorted(stations), opt["coords"], opt["max_dist"])[0].size
    raw_fs = max((r["fs"] for v in rows.values() for r in v), default=None)
    res = resources(cfg, len(stations), npair, len(days), raw_fs)
    from .resources import preflight

    preflight(cfg, res, out_dir)
    tasks = []
    edge = 250.0 / opt["fs"] + 5.0  # read_day looks this far beyond the day (filter run-in, guard)
    for lab, (t0, trs) in days.items():
        by = {}
        for r in trs:
            by.setdefault(r["station"], []).append(r)
        for s, rr in by.items():  # the records next to the day (the day files of the archive end
            have = {id(r) for r in rr}  # at its limits): without them its first and last seconds are
            rr += [
                r
                for r in rows[s]
                if id(r) not in have  # the edge of a record and their segments are lost
                and r["end"] > t0 - edge
                and r["start"] < t0 + 86400 + edge
            ]
        tasks.append((lab, t0, by, opt))
    tasks.sort(key=lambda t: -len(t[2]))  # biggest days first so memory peaks early
    from . import resources as R

    plan = R.plan(cfg, len(tasks), res["memory_per_process_gb"] * 1e9, "days")
    opt["threads"] = plan["threads_each"]  # readers and pair sums of a day (cc_kernel) share the cores of its process
    if plan["limited_by_memory"]:
        print("  %s\n  -> %d days at a time, so that they fit into the memory" % (resources_text(cfg, res),
              plan["parallel"]), flush=True)  # fmt: skip
    msgs, failed, skipped = {}, [], []
    with Pool(plan["parallel"], initializer=R.limit_threads, initargs=(plan["threads_each"],)) as pool:
        for day, msg, fails in track(pool_results(pool, process_day, tasks), len(tasks), "days"):
            print("  %s %s" % (day, msg), flush=True)
            msgs[day] = msg
            failed += ["%s %s" % (day, f) for f in fails]
            if msg.startswith("too few stations"):
                skipped.append("%s: %s" % (day, msg))
        R.note_workers(w.pid for w in pool._pool)
    keep = {"day_%s.h5" % t[0] for t in tasks}
    for f in os.listdir(out_dir):  # day files of days no longer processed (other dates, other labels)
        if f.startswith("day_") and f.endswith(".h5") and f not in keep:
            os.remove(os.path.join(out_dir, f))
    # failed: station-days that could not be read (the step is then incomplete); skipped: days with fewer than two
    # stations (nothing to correlate, not an error)
    n = lambda w: sum(1 for m in msgs.values() if m.startswith(w))
    return dict(
        pairs=int(npair),
        max_distance_km=c["max_distance_km"],
        pulse_removed="all stations" if decomb is True else list(decomb or []),
        days=len(tasks),
        computed=n("ok"),
        recomputed=n("recomputed"),
        reused=n("exists"),
        messages=msgs,
        failed=failed,
        skipped_days=skipped,
    )
