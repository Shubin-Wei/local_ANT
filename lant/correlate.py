"""Step 2 (lant cc): cross-correlation of every station pair, day by day, and the stack.

Input: the cleaned days of lant clean (1_clean/days), without the stations of data.exclude_stations and the days
outside data.start_date / end_date.  For every day:
  1. segments of cc.segment_s, overlapping by cc.overlap (Welch averaging, Seats et al. 2012); a segment is not used
     when it has a gap or when its largest |amplitude| exceeds cc.max_over_std x the standard deviation of the
     station-day (earthquakes, glitches; NoisePy's max_over_std, Jiang & Denolle 2020); each segment: linear trend
     removed, 5 % cosine taper;
  2. cross-coherence of the two spectra, C_ij(f) = F_i*(f) F_j(f) / (|F_i(f)| |F_j(f)|): spectral whitening of both
     records (Bensen et al. 2007), the form compared and recommended by Xu et al. (2021) for high-frequency dense
     arrays; kept between cc.fmin_hz and cc.fmax_hz with quarter-octave cosine edges;
  3. the segments of each sub-stack (cc.substack_hours) averaged in the frequency domain and transformed back to lags
     -maxlag..+maxlag.  Positive lag: energy travelling from station 1 to station 2 of the pair (names sorted).
The stack is the average of all segments (linear stack, Bensen et al. 2007).

Results in OUTPUT/2_cc:
  days/YYYY-MM-DD.h5   sub-stacks of one day: ccf [sub-stack, pair, lag], count [sub-stack, pair] (segments),
                       pair (index into stack.h5), stations; a day computed with the same settings from the same
                       cleaned day is reused
  stack.h5             ccf [pair, lag], count [pair], sta1, sta2, lon1 lat1 lon2 lat2, dist_km (WGS84), az_deg
"""

import glob
import os
import time
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Pool

import h5py
import numpy as np

from . import cc_kernel as K
from . import clean, geo, resources, stations
from . import config as C
from . import records as R
from .step import Step, need

PAIR_CHUNK = 2000  # pairs transformed back to lags at a time
SEGMENT_BLOCK = 64  # segment spectra held in memory at a time


def band_taper(freqs, b0, b1):
    """Cosine edges over a quarter octave inside each end of the band."""
    f = freqs[b0:b1]
    lo, hi = f[0], f[-1]
    t = np.ones(f.size)
    r = np.log2(f / lo) / 0.25
    t = np.where(r < 1, 0.5 * (1 - np.cos(np.pi * np.clip(r, 0, 1))), t)
    r = np.log2(hi / f) / 0.25
    t = np.where(r < 1, t * 0.5 * (1 - np.cos(np.pi * np.clip(r, 0, 1))), t)
    return t


def segment_starts(n_samples, seg, overlap):
    step = max(1, int(round(seg * (1 - overlap))))
    return np.arange(0, n_samples - seg + 1, step)


def day_ccf(X, OK, iu, ju, opt, nsub, threads=1):
    """Cross-coherence of the pairs (iu, ju) of one day.

    X [nsta, 86400 fs] float32 data, OK the valid mask.  -> (ccf [nsub, npair, 2L+1] float32, count [nsub, npair],
    used [nsta] segments used per station, gapfree [nsta] segments without a gap)."""
    fs = opt["fs"]
    n = int(round(opt["segment"] * fs))
    L = int(round(opt["maxlag"] * fs))
    freqs = np.arange(n // 2 + 1) * fs / n
    band = np.where((freqs >= opt["fmin"]) & (freqs <= opt["fmax"]))[0]
    b0, b1 = band[0], band[-1] + 1
    taper = band_taper(freqs, b0, b1)
    win = np.ones(n)
    tl = max(1, int(0.05 * n))
    w = 0.5 * (1 - np.cos(np.pi * np.arange(tl) / tl))
    win[:tl], win[-tl:] = w, w[::-1]
    ns, npair = X.shape[0], iu.size
    starts = segment_starts(X.shape[1], n, opt["overlap"])
    day_std = np.array([X[i, OK[i]].std() if OK[i].any() else np.nan for i in range(ns)])
    sub_len = X.shape[1] / nsub
    sub_of = np.minimum(nsub - 1, ((starts + n // 2) // sub_len).astype(int))
    ccf = np.zeros((nsub, npair, 2 * L + 1), np.float32)
    count = np.zeros((nsub, npair), np.int32)
    used = np.zeros(ns, np.int64)
    gapfree = np.zeros(ns, np.int64)

    def good_of(s0):
        seg = X[:, s0 : s0 + n]
        g = OK[:, s0 : s0 + n].all(axis=1)
        gf = g.copy()
        if opt["max_over_std"] > 0:
            with np.errstate(invalid="ignore"):
                g &= np.abs(seg).max(axis=1) <= opt["max_over_std"] * day_std
        g &= seg.std(axis=1) > 0
        return gf, g

    def spectra(s0):
        seg = K.detrend_linear(X[:, s0 : s0 + n].astype(np.float64))
        F = np.fft.rfft(seg * win, axis=1)[:, b0:b1]
        return (F / (np.abs(F) + 1e-30) * taper).astype(np.complex64)

    with ThreadPoolExecutor(threads) as pool:
        for s in range(nsub):
            ks = starts[sub_of == s]
            if not ks.size:
                continue
            goods = [good_of(s0) for s0 in ks]
            keep = [(s0, g) for s0, (gf, g) in zip(ks, goods) if g.sum() >= 2]
            for gf, g in goods:
                gapfree += gf
                used += g
            if not keep:
                continue
            acc = np.zeros((npair, b1 - b0), np.complex128)
            part = np.empty_like(acc)
            cnt = np.zeros(npair, np.int32)
            for k0 in range(0, len(keep), SEGMENT_BLOCK):  # a block of segments at a time: bounded memory
                block = keep[k0 : k0 + SEGMENT_BLOCK]
                spec = np.zeros((len(block), ns, b1 - b0), np.complex64)
                good = np.zeros((len(block), ns), bool)
                for k, F in enumerate(pool.map(spectra, [s0 for s0, _ in block])):
                    g = block[k][1]
                    F[~g] = 0
                    spec[k], good[k] = F, g
                K.pair_sums(spec, good, len(block), iu, ju, part, threads)
                acc += part
                cnt += (good[:, iu] & good[:, ju]).sum(axis=0).astype(np.int32)
            for c0 in range(0, npair, PAIR_CHUNK):
                c1 = min(npair, c0 + PAIR_CHUNK)
                full = np.zeros((c1 - c0, n // 2 + 1), np.complex128)
                full[:, b0:b1] = acc[c0:c1] / np.maximum(cnt[c0:c1], 1)[:, None]
                x = np.fft.irfft(full, n=n, axis=1)
                ccf[s, c0:c1] = np.concatenate([x[:, -L:], x[:, : L + 1]], axis=1)
            count[s] = cnt
    return ccf, count, used, gapfree


SETTINGS = ("fs", "fmin", "fmax", "maxlag", "segment", "overlap", "max_over_std", "substack_hours")


def _settings(c, fs):
    return dict(fs=fs, fmin=float(c["fmin_hz"]), fmax=float(c["fmax_hz"]), maxlag=float(c["maxlag_s"]),
                segment=float(c["segment_s"]), overlap=float(c["overlap"]), max_over_std=float(c["max_over_std"]),
                substack_hours=int(c["substack_hours"]))  # fmt: skip


def _reusable(path, opt, names, made):
    if not os.path.exists(path):
        return False
    try:
        with h5py.File(path, "r") as h:
            same = all(h.attrs.get(k) == opt[k] for k in SETTINGS)
            return same and list(h.attrs.get("requested", [])) == names and h.attrs.get("cleaned") == made
    except OSError:
        return False


G = {}  # per worker process: what every day needs (set by _init)


def _init(shared, threads):
    import numba

    G.update(shared)
    G["threads"] = threads
    numba.set_num_threads(max(1, min(threads, numba.config.NUMBA_NUM_THREADS)))
    resources.limit_threads(1)


def _one_day(src):
    opt = G["opt"]
    name = os.path.basename(src)[:-3]
    path = os.path.join(G["out_dir"], name + ".h5")
    with h5py.File(src, "r") as h:
        want = sorted(set(s.decode() for s in h["stations"][:]) & G["names"])
        made = h.attrs["made"]
    if len(want) < 2:
        return name, "fewer than two stations: skipped", 0
    if _reusable(path, opt, want, made):
        return name, "reused (same settings and cleaned day)", 0
    start = time.time()
    have, X, OK, _ = clean.load_day(src, want)
    keep = OK.any(axis=1)
    have, X, OK = [s for s, k in zip(have, keep) if k], X[keep], OK[keep]
    if len(have) < 2:
        return name, "fewer than two stations with data: skipped", 0
    pi = G["pair_index"]
    idx = np.array([G["station_pos"][s] for s in have])
    iu, ju = np.triu_indices(len(have), 1)
    gp = np.array([pi.get((idx[a], idx[b]), -1) for a, b in zip(iu, ju)], np.int64)
    sel = gp >= 0
    iu, ju, gp = iu[sel], ju[sel], gp[sel]
    nsub = 24 // opt["substack_hours"]
    ccf, count, used, gapfree = day_ccf(X, OK, iu, ju, opt, nsub, G["threads"])
    with h5py.File(path + ".tmp", "w") as h:
        h["ccf"] = ccf
        h["count"] = count
        h["pair"] = gp
        h["stations"] = np.array(have, "S")
        h["segments_used"] = used
        h["segments_gapfree"] = gapfree
        h.attrs.update(opt)
        h.attrs["requested"] = want
        h.attrs["cleaned"] = made
    os.replace(path + ".tmp", path)
    return name, "%d stations, %d pairs, %.0f s" % (len(have), gp.size, time.time() - start), int(used.sum())


def stack_days(day_files, npair, nlag):
    """Linear stack: the average of all segments, sum_k count_k ccf_k / sum_k count_k."""
    tot = np.zeros((npair, nlag))
    cnt = np.zeros(npair, np.int64)
    for f in day_files:
        with h5py.File(f, "r") as h:
            gp = h["pair"][:]
            c = h["count"][:].astype(np.float64)
            for s in range(c.shape[0]):
                tot[gp] += c[s][:, None] * h["ccf"][s].astype(np.float64)
                cnt[gp] += c[s].astype(np.int64)
    with np.errstate(invalid="ignore"):
        return (tot / np.maximum(cnt, 1)[:, None]).astype(np.float32), cnt


def pairs_of(names, coords, max_dist=None):
    """Pairs (i < j) of the sorted names with distance (km, WGS84) and azimuth from station 1 to 2."""
    iu, ju = np.triu_indices(len(names), 1)
    lon = np.array([coords[n]["lon"] for n in names])
    lat = np.array([coords[n]["lat"] for n in names])
    d, az = geo.geodesic_km(lon[iu], lat[iu], lon[ju], lat[ju])
    d, az = np.atleast_1d(d), np.atleast_1d(az)
    keep = np.ones(iu.size, bool) if not max_dist else d <= float(max_dist)
    return iu[keep], ju[keep], d[keep], az[keep]


def _cleaned(cfg):
    """The cleaned days and stations of lant clean within the selection of data.* -> (files, names, coords, fs)."""
    folder = os.path.join(cfg["output_dir"], "1_clean")
    table = need(os.path.join(folder, "stations_used.csv"), "clean")
    coords = stations.read_csv(table)
    names = sorted(set(coords) - set(cfg["data"]["exclude_stations"]))
    lo, hi = cfg["data"]["start_date"], cfg["data"]["end_date"]
    files = []
    for f in sorted(glob.glob(os.path.join(folder, "days", "*.h5"))):
        day = os.path.basename(f)[:-3]
        if (lo is None or day >= str(lo)) and (hi is None or day <= str(hi)):
            files.append(f)
    if not files:
        raise R.DataError("no cleaned day in %s within data.start_date / end_date: run `lant clean`" % folder)
    rates = set()
    for f in files:
        with h5py.File(f, "r") as h:
            rates.add(float(h.attrs["fs"]))
    if len(rates) > 1:
        raise R.DataError("the cleaned days have different sampling rates %s: run `lant clean` again" % sorted(rates))
    fs = rates.pop()
    if cfg["clean"]["fs_hz"] is not None and float(cfg["clean"]["fs_hz"]) != fs:
        raise R.DataError("the cleaned days are at %g Hz but clean.fs_hz is %g: run `lant clean` again"
                          % (fs, float(cfg["clean"]["fs_hz"])))  # fmt: skip
    return files, names, {n: coords[n] for n in names}, fs


def run(cfg):
    C.check_cc(cfg)
    with Step(cfg, "2_cc", ["data", "cc", "cores"]) as step:
        n_cores, how = resources.cores(cfg["cores"])
        step.say("cores: %d (%s)" % (n_cores, how))
        sources, names, coords, fs = _cleaned(cfg)
        step.derive("fs_hz", fs, "sampling rate of the cleaned days (1_clean/days)")
        C.check_fs(fs, float(cfg["cc"]["fmax_hz"]), "the sampling rate of the cleaned days")
        if cfg["data"]["response"]:
            with h5py.File(sources[0], "r") as h:
                pre = list(h.attrs["prefilter_hz"])
            if not pre or pre[0] > cfg["cc"]["fmin_hz"] + 1e-9 or pre[1] < cfg["cc"]["fmax_hz"] - 1e-9:
                raise R.DataError("the responses were removed with a pre-filter flat over %s Hz, not the whole "
                                  "correlation band: run `lant clean` again" % pre)  # fmt: skip
        if len(names) < 2:
            raise R.DataError("fewer than two cleaned stations (after data.exclude_stations)")
        iu, ju, dist, az = pairs_of(names, coords, cfg["cc"]["max_distance_km"])
        step.say("%d stations, %d pairs, %d cleaned days (%s to %s) at %g Hz"
                 % (len(names), iu.size, len(sources), os.path.basename(sources[0])[:-3],
                    os.path.basename(sources[-1])[:-3], fs))  # fmt: skip
        opt = _settings(cfg["cc"], fs)
        nlag = 2 * int(round(opt["maxlag"] * opt["fs"])) + 1
        day_bytes = len(names) * R.DAY * opt["fs"] * 5 + iu.size * nlag * 4 * (24 // opt["substack_hours"] + 2)
        nproc = resources.parallel(len(sources), n_cores, day_bytes * 1.5)
        threads = max(1, n_cores // nproc)
        step.say("days in parallel: %d, threads each: %d (about %.1f GB per day)" % (nproc, threads, day_bytes / 1e9))
        out_dir = step.path("days")
        os.makedirs(out_dir, exist_ok=True)
        shared = dict(opt=opt, out_dir=out_dir, names=set(names),
                      station_pos={n: k for k, n in enumerate(names)},
                      pair_index={(int(a), int(b)): k for k, (a, b) in enumerate(zip(iu, ju))})  # fmt: skip
        with Pool(nproc, initializer=_init, initargs=(shared, threads)) as pool:
            for k, (day, msg, _) in enumerate(pool.imap_unordered(_one_day, sources)):
                step.say("  day %s: %s   [%d/%d]" % (day, msg, k + 1, len(sources)))
        files = [os.path.join(out_dir, os.path.basename(f)) for f in sources]
        files = [f for f in files if os.path.exists(f)]
        stale = set(os.listdir(out_dir)) - {os.path.basename(f) for f in files}
        for f in sorted(stale):  # days outside the present selection would otherwise enter the stack
            os.remove(os.path.join(out_dir, f))
            step.say("  removed %s (not a day of this run)" % f)
        if not files:
            raise R.DataError("no day has two stations with data")
        step.say("stacking %d days" % len(files))
        ccf, cnt = stack_days(files, iu.size, nlag)
        with h5py.File(step.path("stack.h5.tmp"), "w") as h:
            h["ccf"] = ccf
            h["count"] = cnt
            h["sta1"] = np.array([names[i] for i in iu], "S")
            h["sta2"] = np.array([names[j] for j in ju], "S")
            for key, idx in (("1", iu), ("2", ju)):
                h["lon" + key] = np.array([coords[names[i]]["lon"] for i in idx])
                h["lat" + key] = np.array([coords[names[i]]["lat"] for i in idx])
            h["dist_km"] = dist
            h["az_deg"] = az
            h["days"] = np.array([os.path.basename(f) for f in files], "S")
            h.attrs.update(opt)
        os.replace(step.path("stack.h5.tmp"), step.path("stack.h5"))
        with_data = int((cnt > 0).sum())
        step.number("stations", len(names))
        step.number("pairs", int(iu.size))
        step.number("pairs with data", with_data)
        step.number("days", len(files))
        step.number("median segments per pair", int(np.median(cnt[cnt > 0])) if with_data else 0)
        from .plotting import cc_figures

        cc_figures.make(cfg, step)
