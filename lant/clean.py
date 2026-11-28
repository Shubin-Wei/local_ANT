"""Step 1 (lant clean): the raw records cut into UTC days and cleaned station by station, one file per day.

The records of the station table are found (lant.records: index and matching), and every station-day is read and
cleaned (lant.records.read_day): the samples of the UTC day 00:00-24:00, gaps kept as gaps, median removed,
instrument response removed when data.response is given, anti-alias low-pass and decimation to clean.fs_hz.
Nothing is discarded here for its amplitude: lant cc leaves out the segments with gaps or large amplitudes.

Results in OUTPUT/1_clean:
  record_index.csv     every raw record found (headers)
  stations_used.csv    stations with records: coordinates, channel used, hours of cleaned data
  availability.csv     hours of cleaned data per station (rows) and day (columns)
  days/YYYY-MM-DD.h5   data [station, sample] float32 at fs, NaN where there is no valid sample; stations (names);
                       attributes t0 (UTC timestamp of 00:00), fs, response, prefilter_hz, made (when written).
                       A day written with the same settings from the same raw files is reused.
"""

import hashlib
import os
import time
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Pool

import h5py
import numpy as np

from . import config as C
from . import records as R
from . import resources, stations
from .step import Step

READ_THREADS = 4
SETTINGS = ("fs", "response", "prefilter_hz")


def save_day(path, names, data, ok, attrs):
    """One cleaned day: data [station, sample] with NaN where ok is False; hours [station] of valid data."""
    x = np.where(ok, data, np.nan).astype(np.float32)
    with h5py.File(path + ".tmp", "w") as h:
        h["data"] = x
        h["stations"] = np.array(names, "S")
        h["hours"] = ok.sum(axis=1) / float(attrs["fs"]) / 3600.0
        h.attrs.update(attrs)
        h.attrs["made"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    os.replace(path + ".tmp", path)


def load_day(path, names=None):
    """-> (names, data [station, sample] float32 with 0 in gaps, valid mask, attrs); names: only these stations."""
    with h5py.File(path, "r") as h:
        have = [s.decode() for s in h["stations"][:]]
        rows = list(range(len(have))) if names is None else [k for k, s in enumerate(have) if s in set(names)]
        x = h["data"][rows] if rows else np.zeros((0, h["data"].shape[1]), np.float32)
        attrs = dict(h.attrs)
    ok = np.isfinite(x)
    x[~ok] = 0.0
    return [have[k] for k in rows], x, ok, attrs


def _sources(recs, names, t0):
    """Fingerprint of the raw files of a day (path, size, modification time): a changed file is cleaned again."""
    rows = sorted({(r["file"], r["size"], r["mtime"]) for s in names for r in recs[s]
                   if r["end"] > t0 and r["start"] < t0 + R.DAY})  # fmt: skip
    return hashlib.sha1(repr(rows).encode()).hexdigest()


def _reusable(path, opt, src):
    if not os.path.exists(path):
        return False
    try:
        with h5py.File(path, "r") as h:
            same = all(k in h.attrs and np.array_equal(h.attrs[k], opt[k]) for k in SETTINGS)
            return same and h.attrs.get("sources") == src
    except OSError:
        return False


G = {}  # per worker process (set by _init)


def _init(shared, threads):
    G.update(shared)
    G["threads"] = threads
    resources.limit_threads(1)


def _one_day(t0):
    opt, recs = G["opt"], G["records"]
    name = R.day_name(t0)
    path = os.path.join(G["out_dir"], name + ".h5")
    want = sorted(s for s, rows in recs.items() if any(r["end"] > t0 and r["start"] < t0 + R.DAY for r in rows))
    src = _sources(recs, want, t0)
    if _reusable(path, opt, src):
        with h5py.File(path, "r") as h:
            names = [s.decode() for s in h["stations"][:]]
            hours = dict(zip(names, h["hours"][:].tolist()))
        return name, "reused (same settings and raw files)", hours
    start = time.time()
    band = tuple(opt["prefilter_hz"]) if G.get("inv") is not None else None
    nout = int(round(R.DAY * opt["fs"]))
    X = np.zeros((len(want), nout), np.float32)
    OK = np.zeros((len(want), nout), bool)
    have, failed = [], []

    def read(s):
        try:
            return s, R.read_day(recs[s], t0, opt["fs"], G.get("inv"), band)
        except Exception as exc:
            return s, exc

    with ThreadPoolExecutor(min(READ_THREADS, G["threads"])) as ex:
        for s, got in ex.map(read, want):
            if isinstance(got, Exception):
                failed.append("%s: %s" % (s, got))
            elif got[0] is not None and got[1].any():
                X[len(have)], OK[len(have)] = got
                have.append(s)
    if failed:  # a station that cannot be read stops the step: the user decides (exclude it or fix the file)
        raise R.DataError("%s: could not read %s" % (name, "; ".join(failed)))
    if not have:
        if os.path.exists(path):
            os.remove(path)
        return name, "no valid data: no file", {}
    save_day(path, have, X[: len(have)], OK[: len(have)], dict(opt, t0=t0, sources=src))
    hours = dict(zip(have, (OK[: len(have)].sum(axis=1) / opt["fs"] / 3600).tolist()))
    return name, "%d stations, %.0f s" % (len(have), time.time() - start), hours


def settings(cfg):
    c = cfg["clean"]
    resp = cfg["data"]["response"]
    pre = [float(cfg["cc"]["fmin_hz"]), float(cfg["cc"]["fmax_hz"])] if resp else []
    return dict(fs=float(c["fs_hz"]), response=resp or "none", prefilter_hz=pre)


def _write_availability(path, hours, names, day_names):
    with open(path, "w") as fh:
        fh.write(",".join(["station"] + day_names) + "\n")
        for s in names:
            fh.write(",".join([s] + ["%.2f" % hours[d].get(s, 0.0) for d in day_names]) + "\n")


def run(cfg):
    C.check_clean(cfg)
    with Step(cfg, "1_clean", ["data", "clean", "cores"]) as step:
        n_cores, how = resources.cores(cfg["cores"])
        step.say("cores: %d (%s)" % (n_cores, how))
        table = stations.read_table(cfg["data"]["station_file"])
        step.say("station table: %d stations" % len(table))
        rows, bad = R.build_index(cfg["data"]["raw_dir"], step.path("record_index.csv"), n_cores, step.say)
        for f, why in bad:
            step.say("  not a seismic record (skipped): %s (%s)" % (f, why))
        recs, notes = R.match(rows, table, cfg["data"]["channel"], set(cfg["data"]["exclude_stations"]))
        for n in notes:
            step.say("  " + n)
        if len(recs) < 2:
            raise R.DataError("fewer than two stations of the table have %s records" % cfg["data"]["channel"])
        days = R.days(recs, cfg["data"]["start_date"], cfg["data"]["end_date"])
        if not days:
            raise R.DataError("no UTC day with records within data.start_date / end_date")
        names = sorted(recs)
        step.say("%d stations, %d days (%s to %s)" % (len(names), len(days), R.day_name(days[0]),
                                                      R.day_name(days[-1])))  # fmt: skip
        opt = settings(cfg)
        inv = R.load_responses(cfg["data"]["response"]) if cfg["data"]["response"] else None
        step.say(
            "instrument responses: %s"
            % (
                "removed (%s), pre-filter flat over cc.fmin_hz-cc.fmax_hz %s Hz" % (cfg["data"]["response"],
                                                                                    opt["prefilter_hz"])
                if inv
                else "not removed (data.response is empty: same instrument everywhere)"
            )
        )  # fmt: skip
        raw_fs = max(r[0]["fs"] for r in recs.values())  # the cleaned day (4 + 1 bytes) and, per reading thread,
        day_bytes = len(names) * R.DAY * opt["fs"] * 5 + READ_THREADS * R.DAY * raw_fs * 24  # raw float64 copies
        nproc = resources.parallel(len(days), n_cores, day_bytes * 1.5)
        threads = max(1, n_cores // nproc)
        step.say("days in parallel: %d, threads each: %d (about %.1f GB per day)" % (nproc, threads, day_bytes / 1e9))
        out_dir = step.path("days")
        os.makedirs(out_dir, exist_ok=True)
        shared = dict(opt=opt, records=recs, out_dir=out_dir, inv=inv)
        hours = {}
        with Pool(nproc, initializer=_init, initargs=(shared, threads)) as pool:
            for k, (day, msg, h) in enumerate(pool.imap_unordered(_one_day, days)):
                step.say("  day %s: %s   [%d/%d]" % (day, msg, k + 1, len(days)))
                hours[day] = h
        keep = {R.day_name(t) + ".h5" for t in days if hours[R.day_name(t)]}
        for f in sorted(set(os.listdir(out_dir)) - keep):  # days outside the present selection
            os.remove(os.path.join(out_dir, f))
            step.say("  removed %s (not a day of this run)" % f)
        day_names = [R.day_name(t) for t in days if hours[R.day_name(t)]]
        total = {s: sum(hours[d].get(s, 0.0) for d in day_names) for s in names}
        used = [s for s in names if total[s] > 0]
        for s in names:
            if total[s] == 0:
                step.say("  %s: no valid cleaned data (records too short or all gaps)" % s)
        extra = {s: dict(hours=round(total[s], 1), channel="%s.%s %g Hz" % (recs[s][0]["location"],
                         recs[s][0]["channel"], recs[s][0]["fs"])) for s in used}  # fmt: skip
        stations.write_csv(step.path("stations_used.csv"), {s: table[s] for s in used}, extra)
        _write_availability(step.path("availability.csv"), hours, used, day_names)
        step.number("stations with data", len(used))
        step.number("days", len(day_names))
        step.number("hours of data (all stations)", round(sum(total.values()), 1))
        step.number("median hours per station", round(float(np.median([total[s] for s in used])), 1) if used else 0)
        from .plotting import clean_figures

        clean_figures.make(cfg, step)
