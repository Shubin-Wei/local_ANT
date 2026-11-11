"""`lant inspect <topic>`: the numbers needed to judge a step and to choose the parameters of the next one.

Topics: data (station geometry, record length, sampling), auto (the parameters derived from the data and their
reasons), and every step name (its summary
plus derived quantities).  The output is JSON (printed, and written to .lant/summary/inspect_<topic>.json), with
the files and figures of the step to look at.  Nothing here changes a result.
"""
import json
import os

import numpy as np

from . import paths, records
from .signal import geodesic_km, nearest_neighbour_km


FIGURES = {
    "clean": ["01_clean/report.txt", "01_clean/coverage.png", "01_clean/station_days.csv", "01_clean/processing_days.txt"],
    "stack": ["03_stack/candidates.txt", "03_stack/station_qc/station_qc.png", "03_stack/station_qc/station_qc.txt",
              "03_stack/removed.txt"],
    "reference": ["04_reference/fj.png", "04_reference/fj.txt"],
    "disp": ["05_dispersion/median_vs_reference.txt", "05_dispersion/pair_figures/", "figures/3_dispersion_curves.png"],
    "inversion": ["06_inversion/average_curve.txt", "06_inversion/vs1d.txt", "06_inversion/main/DSurfTomo.in",
                  "figures/4_start_model.png"],
    "evaluate": ["06_inversion/evaluation/evaluation.json", "figures/6_residuals.png", "figures/6_recovery.png"],
    "figures": ["figures/"],
}


def geometry(stations):
    names = sorted(stations)
    lon = np.array([stations[n]["lon"] for n in names])
    lat = np.array([stations[n]["lat"] for n in names])
    iu, ju = np.triu_indices(len(names), 1)
    d, _ = geodesic_km(lon[iu], lat[iu], lon[ju], lat[ju])
    nn = nearest_neighbour_km(lon, lat)
    from .signal import local_xy_km
    x, y = local_xy_km(lon, lat)
    w = np.linalg.eigvalsh(np.cov(np.c_[x, y].T)) if len(names) > 2 else np.array([0, 1])
    elong = float(np.sqrt(w.max() / max(w.min(), 1e-12)))
    return dict(stations=len(names), pairs=int(len(d)),
                interstation_km=dict(min=float(d.min()), p10=float(np.percentile(d, 10)), median=float(np.median(d)),
                                     p90=float(np.percentile(d, 90)), max=float(d.max())),
                nearest_neighbour_km=dict(median=float(np.nanmedian(nn)), max=float(np.nanmax(nn))),
                extent_km=dict(east_west=float(x.max() - x.min()), north_south=float(y.max() - y.min())),
                elongation=elong, looks_like="line" if elong > 4 else "areal",
                note="the longest usable period is about distance / (2 c): the reference step measures c and derives the band")


def data(cfg):
    from . import data as D
    out = {"configuration": cfg.get("_config_file")}
    try:
        st = D.used_stations(cfg)
        out["geometry"] = geometry(st)
    except D.DataError as exc:
        out["station_table_problem"] = str(exc)
    s = records.read_summary(cfg, "clean")
    if s:
        out["clean"] = s
    else:
        out["clean"] = "not run yet: `lant run clean` scans the raw files"
    out["correlation_band"] = dict(fs=cfg["correlation"]["fs"], fmin=cfg["correlation"]["fmin_hz"], fmax=cfg["correlation"]["fmax_hz"])
    return out


def stack(cfg, s):
    import h5py
    p = paths.P(cfg, paths.STACK, "stack.h5")
    if os.path.exists(p) and not isinstance(cfg["correlation"]["seg_step_s"], str):
        with h5py.File(p, "r") as h:
            nseg = h["nseg"][:]
            dist = h["dist_km"][:]
            hours = nseg * (cfg["correlation"]["seg_step_s"] / 3600)
        ok = nseg > 0
        s["pairs_with_data"] = int(ok.sum())
        s["pair_hours"] = dict(p10=float(np.percentile(hours[ok], 10)), median=float(np.median(hours[ok])), max=float(hours.max()))
        s["distance_km_of_pairs_with_data"] = dict(min=float(dist[ok].min()), median=float(np.median(dist[ok])), max=float(dist[ok].max()))
    return s


def disp(cfg, s):
    import h5py
    from .qc import reason_counts
    p = paths.P(cfg, paths.DISP, "disp.h5")
    if os.path.exists(p):
        with h5py.File(p, "r") as h:
            f = h["freqs"][:]
            okc = h["ok_c"][:]
            pre = h["pre_c"][:]
        s["phase_points_per_frequency"] = {"%.3f" % x: int(n) for x, n in zip(f, okc.sum(0))}
        s["candidate_points_before_phase_qc"] = int(pre.sum())
        s["rejected_samples_by_reason"] = reason_counts(p)
    return s


def inspect(cfg, topic):
    from . import auto
    if topic == "auto":
        p = auto.path(cfg)
        out = json.load(open(p)) if os.path.exists(p) else {"status": "run `lant run clean` first"}
    elif topic == "data":
        cfg = auto.apply(cfg)
        out = data(cfg)
    else:
        cfg = auto.apply(cfg)
        s = records.read_summary(cfg, topic)
        if s is None:
            out = {"status": "step %s has no summary yet (not run, or it failed: see .lant/run_log.jsonl)" % topic}
        else:
            out = {"stack": stack, "disp": disp}.get(topic, lambda c, x: x)(cfg, s)
    files = [paths.P(cfg, *f.rstrip("/").split("/")) for f in FIGURES.get(topic, [])]
    out = {"topic": topic, "result": out, "look_at": [f for f in files if os.path.exists(f.rstrip("/"))]}
    return out


def status(cfg):
    """The state of every step (lant.manifest.status) and the next thing to do."""
    from . import manifest, report
    rcfg = report.resolved(cfg)
    st = manifest.status(rcfg)
    cmd, why = manifest.next_action(rcfg, st)
    from .artifacts import audit
    from .delivery import RULES
    return dict(steps=st, next=dict(command=cmd, why=why), artifacts=audit(cfg), rules=RULES)
