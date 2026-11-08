"""`lant inspect <topic>`: the numbers needed to judge a step and to choose the parameters of the next one.

Topics: data (station geometry, record length, sampling), auto (the parameters derived from the data and their
reasons), and every step name (its summary
plus derived quantities).  The output is JSON (printed, and written to records/summary/inspect_<topic>.json), with
the files and figures of the step to look at.  Nothing here changes a result.
"""
import json
import os

import numpy as np

from . import paths, records
from .signal import geodesic_km


FIGURES = {
    "clean": ["01_clean/report.txt", "01_clean/coverage.png", "01_clean/station_days.csv"],
    "stack": ["03_stack/station_day_qc/station_day_qc.png", "03_stack/station_day_qc/station_day_qc.txt",
              "03_stack/station_hour_qc/station_hour_qc.txt", "03_stack/exclude_station_days.txt"],
    "notch": ["04_notch/notch_report.png", "04_notch/notch_table.txt", "04_notch/raw_lines.txt"],
    "reference": ["05_reference/fj.png", "05_reference/fj.txt"],
    "disp1": ["05_reference/reference_pass2.png", "05_reference/reference_pass2.txt"],
    "disp2": ["06_dispersion/precision.txt", "06_dispersion/pair_figures/", "figures/3_dispersion_curves.png"],
    "inversion": ["07_inversion/average_curve.txt", "07_inversion/vs1d.txt", "07_inversion/main/DSurfTomo.in",
                  "figures/4_start_model.png"],
    "evaluate": ["07_inversion/evaluation/evaluation.json", "figures/6_residuals.png", "figures/6_recovery.png"],
    "figures": ["figures/"],
}


def geometry(stations):
    names = sorted(stations)
    lon = np.array([stations[n]["lon"] for n in names])
    lat = np.array([stations[n]["lat"] for n in names])
    iu, ju = np.triu_indices(len(names), 1)
    d, _ = geodesic_km(lon[iu], lat[iu], lon[ju], lat[ju])
    nn = []
    for i in range(len(names)):
        m = (iu == i) | (ju == i)
        nn.append(d[m].min() if m.any() else np.nan)
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
                note="longest usable period ~ (distance / n_lambda) / c; e.g. with 2 wavelengths and c = 2.5 km/s the "
                     "p90 distance allows T up to %.2f s; shortest periods need distances >= 1-2 wavelengths"
                     % (np.percentile(d, 90) / 2 / 2.5))


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
    p = paths.P(cfg, paths.STACK, "stack_clean.h5")
    if os.path.exists(p):
        with h5py.File(p, "r") as h:
            nseg = h["nseg"][:]
            dist = h["dist_km"][:]
            hours = nseg * (cfg["correlation"]["seg_step_s"] / 3600)
        ok = nseg > 0
        s["pairs_with_data"] = int(ok.sum())
        s["pair_hours"] = dict(p10=float(np.percentile(hours[ok], 10)), median=float(np.median(hours[ok])), max=float(hours.max()))
        s["distance_km_of_pairs_with_data"] = dict(min=float(dist[ok].min()), median=float(np.median(dist[ok])), max=float(dist[ok].max()))
    return s


def disp2(cfg, s):
    import h5py
    p = paths.P(cfg, paths.DISP, "disp_final.h5")
    if os.path.exists(p):
        with h5py.File(p, "r") as h:
            f = h["freqs"][:]
            okc = h["ok_c"][:]
            okg = h["ok_g"][:]
            pre = h["pre_c"][:]
        s["phase_points_per_frequency"] = {"%.3f" % x: int(n) for x, n in zip(f, okc.sum(0))}
        s["group_points_per_frequency"] = {"%.3f" % x: int(n) for x, n in zip(f, okg.sum(0))}
        s["candidate_points_before_phase_qc"] = int(pre.sum())
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
        s = records.read_summary(cfg, topic)
        if s is None:
            out = {"status": "step %s has no summary yet (not run, or it failed: see records/run_log.jsonl)" % topic}
        else:
            out = {"stack": stack, "disp2": disp2}.get(topic, lambda c, x: x)(cfg, s)
    files = [os.path.join(cfg["project"]["output_dir"], f) for f in FIGURES.get(topic, [])]
    out = {"topic": topic, "result": out, "look_at": [f for f in files if os.path.exists(f.rstrip("/"))]}
    p = os.path.join(records.rec_dir(cfg), "summary", "inspect_%s.json" % topic)
    with open(p, "w") as fh:
        json.dump(out, fh, indent=1, default=records._json_default, ensure_ascii=False)
    return out


def status(cfg):
    """Which steps have run, with their last state, from records/run_log.jsonl."""
    from .pipeline import STEPS
    log = os.path.join(cfg["project"]["output_dir"], "records", "run_log.jsonl")
    last = {}
    if os.path.exists(log):
        for line in open(log):
            e = json.loads(line)
            last[e["step"]] = e
    return {s: ({"state": last[s]["event"], "time": last[s]["time"], **({"error": last[s]["error"]} if "error" in last[s] else {})}
                if s in last else {"state": "not run"}) for s in STEPS}
