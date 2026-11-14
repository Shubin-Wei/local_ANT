"""Model I/O, exact forward evaluation and recovery statistics; scheduling stays in lant.inversion.run."""

import glob
import json
import os

import numpy as np

from .. import resources


def read_in(path):
    L = open(path).read().split("\n")
    nx, ny, nz = (int(v) for v in L[4].split()[:3])
    k = next(i for i, l in enumerate(L) if "c: kmaxRc" in l)
    n = int(L[k].split()[0])
    periods = np.array([float(v) for v in L[k + 1].split()[:n]])
    return L, nx, ny, nz, periods


def read_mod(path, nx, ny, nz):
    v = np.array(open(path).read().split(), float)
    if v.size != nz + nx * ny * nz:
        raise RuntimeError("%s: %d numbers, expected %d" % (path, v.size, nz + nx * ny * nz))
    return v[:nz], v[nz:].reshape(nz, ny, nx)


def check_node_order(m, nx, ny, nz, path="model"):
    """The rows of a DSurfTomo model file ("lon lat depth vs", depth slowest, then longitude, latitude fastest) must
    be the (nx - 2) x (ny - 2) x (nz - 1) inner nodes, each once, in that order; raises otherwise."""
    want = (nx - 2) * (ny - 2) * (nz - 1)
    if m.ndim != 2 or m.shape[0] != want or m.shape[1] < 4:
        raise RuntimeError("%s: %d rows, expected %d" % (path, m.shape[0] if m.ndim == 2 else 0, want))
    M = m[:, :3].reshape(nz - 1, ny - 2, nx - 2, 3)
    tol = 1e-4
    ok = (
        np.ptp(M[..., 2], axis=(1, 2)).max() <= tol  # depth constant within a depth block
        and np.ptp(M[..., 0], axis=(0, 2)).max() <= tol  # longitude constant along depth and latitude
        and np.ptp(M[..., 1], axis=(0, 1)).max() <= tol
    )  # latitude constant along depth and longitude
    if ok:
        ok = (
            len(np.unique(np.round(M[:, 0, 0, 2], 4))) == nz - 1
            and len(np.unique(np.round(M[0, :, 0, 0], 5))) == ny - 2
            and len(np.unique(np.round(M[0, 0, :, 1], 5))) == nx - 2
        )
    if not ok:
        raise RuntimeError(
            "%s: the rows are not the inner grid nodes, each once, in DSurfTomo's order (repeated, "
            "missing or reordered nodes)" % path
        )


def full_grid(project, model, nx, ny, nz):
    """Starting MOD with the inner nodes (not boundary, not the fixed bottom node) taken from `model`."""
    depz, vs = read_mod(os.path.join(project, "MOD"), nx, ny, nz)
    if model != "start":
        m = np.loadtxt(model)
        if len(m) != (nx - 2) * (ny - 2) * (nz - 1):
            raise RuntimeError("%s: %d rows, expected %d" % (model, len(m), (nx - 2) * (ny - 2) * (nz - 1)))
        check_node_order(m, nx, ny, nz, model)
        vs = vs.copy()
        vs[: nz - 1, 1 : ny - 1, 1 : nx - 1] = m[:, 3].reshape(nz - 1, ny - 2, nx - 2)
    return depz, vs


def read_data(path, periods):
    rows, src = [], None
    for line in open(path):
        p = line.split()
        if not p:
            continue
        if p[0] == "#":
            src = (float(p[1]), float(p[2]), int(p[3]))
        else:
            rows.append((src[0], src[1], float(p[0]), float(p[1]), periods[src[2] - 1], float(p[2])))
    return np.array(rows)


def _forward_once(cfg, project, depz, vs, synthetic=False, threads=None):
    """DSurfTomo's first iteration on (depz, vs) without the inversion -> rows of residualFirst.dat (distance,
    predicted time, observed time, weighted predicted, weighted observed, weight by the rejection rule).  With
    synthetic=True the project's MOD.true gives the observed times, as in its own run."""
    from .solver import evaluate

    R = evaluate(project, depz, vs, synthetic, threads or resources.cores(cfg["resources"]["cores"])[0])
    if not np.all(np.isfinite(R[:, :3])) or (R[:, 1] <= 0).any():
        raise RuntimeError("forward run of %s: non-finite or non-positive predicted travel times" % project)
    return R


def logged_residuals(project):
    """Standard deviation (ms) of the weighted residuals DSurfTomo printed for every iteration (run.log)."""
    import re

    out = []
    try:
        for line in open(os.path.join(project, "run.log"), errors="replace"):
            m = re.search(r"residual after weighting:\s*(-?[\d.]+)\s*ms\s+(-?[\d.]+)\s*ms", line)
            if m:
                out.append(float(m.group(2)))
            m = re.search(r"weighted residual standard deviation:\s*([\d.eE+-]+)\s*ms", line)
            if m:
                out.append(float(m.group(1)))
    except OSError:
        pass
    return out


def evaluate(cfg, project, model="", threads=None):
    """Predicted travel times of a model for every row of the project's data file, by the solver's own forward
    operator (its first iteration from the model; the columns of residualFirst.dat: distance, predicted and
    observed time).  DSurfTomo's residualLast.dat belongs to the model of the previous iteration and its log mixes
    rejected rows in, so neither is used."""
    return _prediction(cfg, project, model or os.path.join(project, "DSurfTomo.inMeasure.dat"), threads=threads)


def _prediction(cfg, project, model="start", synthetic=False, threads=None):
    from . import forward_cache

    snapshot = forward_cache.input_identity(project, model, synthetic)
    _, nx, ny, nz, periods = read_in(os.path.join(project, "DSurfTomo.in"))
    depz, vs = full_grid(project, model, nx, ny, nz)
    data = read_data(os.path.join(project, "surfdata.dat"), periods)

    def unchanged_inputs():
        if forward_cache.input_identity(project, model, synthetic) != snapshot:
            raise RuntimeError("forward input files changed during evaluation; no cached result published")

    def validate(result):
        if result.shape != (len(data), 6) or not np.all(np.isfinite(result)) or np.any(result[:, 1:3] <= 0):
            raise RuntimeError("forward run of %s: incomplete or invalid travel times/weights" % project)
        if not synthetic and not np.allclose(result[:, 0] / data[:, 5], result[:, 2], rtol=2e-3):
            raise RuntimeError("row order check failed (distance / velocity differs from the program's observed time)")

    result = forward_cache.prediction(
        cfg,
        project,
        depz,
        vs,
        lambda: _forward_once(cfg, project, depz, vs, synthetic=synthetic, threads=threads),
        validate,
        unchanged_inputs,
        synthetic=synthetic,
    )
    return data, result, vs


def precheck(cfg, project, threads=None):
    """Before a project runs: all its data against its starting model, by DSurfTomo's own first iteration
    (checkerboards: against their synthetic data).  DSurfTomo's rejection rule (residual < threshold x q25 or >
    threshold x q75, on the signed quartiles; not changed) drops nearly every datum when the residuals all have
    one sign, e.g. a starting model that is too fast or too slow overall; such a project would 'finish' without
    using any data.  It is invalid when fewer than checks.precheck.min_kept of the data survive the rule, or when
    the median residual is further from 0 than checks.precheck.max_median_iqr x the interquartile range.
    Writes <project>/precheck.json; -> dict(valid, reasons, stats)."""

    pc = cfg["checks"]["precheck"]
    _, R, _ = _prediction(cfg, project, synthetic=os.path.exists(os.path.join(project, "MOD.true")), threads=threads)
    st = residual_stats(R)
    iqr = st["q75_ms"] - st["q25_ms"]
    reasons = []
    if st["kept_fraction"] < pc["min_kept"]:
        reasons.append(
            "only %.1f %% of the %d travel times survive DSurfTomo's rejection rule (checks.precheck.min_kept "
            "%.0f %%)" % (100 * st["kept_fraction"], st["n"], 100 * pc["min_kept"])
        )
    if abs(st["median_ms"]) > pc["max_median_iqr"] * iqr:
        reasons.append(
            "median residual %.2f ms is beyond %.1f x the interquartile range %.2f ms: the starting model is "
            "too fast or too slow overall" % (st["median_ms"], pc["max_median_iqr"], iqr)
        )
    out = dict(valid=not reasons, reasons=reasons, stats=st)
    json.dump(out, open(os.path.join(project, "precheck.json"), "w"), indent=1)
    return out


def residual_stats(R):
    """Residual (observed - predicted travel time) of all data, and how many DSurfTomo's rejection rule drops.
    The rule (residual < q25 x threshold or > q75 x threshold, on the signed quartiles) is DSurfTomo's own and is
    not changed; it is centred on the median only when the residuals are: median and quartiles are given so that a
    model that is too fast or too slow overall (all residuals of one sign) shows."""
    r = R[:, 2] - R[:, 1]
    q25, q50, q75 = np.percentile(r, [25, 50, 75])
    kept = R[:, 5] != 0
    return dict(
        n=len(r),
        rms_ms=1e3 * float(np.sqrt(np.mean(r**2))),
        rel_rms_percent=100 * float(np.sqrt(np.mean((r / R[:, 2]) ** 2))),
        mean_ms=1e3 * float(r.mean()),
        median_ms=1e3 * float(q50),
        q25_ms=1e3 * float(q25),
        q75_ms=1e3 * float(q75),
        rejected_by_rule=int((~kept).sum()),
        kept_fraction=float(kept.mean()),
        kept_median_ms=1e3 * float(np.median(r[kept])) if kept.any() else None,
        zero_predicted=int((R[:, 1] <= 0).sum()),
    )


def residuals(cfg, project, out_prefix, iterations=True):
    """Final-model residuals (csv + stats) and the residual of every iteration's model on the same data."""
    D, R, _ = evaluate(cfg, project)
    with open(out_prefix + ".csv", "w") as f:
        f.write("src_lat,src_lon,rcv_lat,rcv_lon,period_s,dist_km,t_obs_s,t_pred_s,kept_by_rule\n")
        for d, r in zip(D, R):
            f.write(
                "%.6f,%.6f,%.6f,%.6f,%.4f,%.5f,%.5f,%.5f,%d\n" % (d[0], d[1], d[2], d[3], d[4], r[0], r[2], r[1], r[5])
            )
    final = residual_stats(R)
    iters = []
    prev = None
    models = (
        ["start"] + sorted(glob.glob(os.path.join(project, "DSurfTomo.inMeasure.dat.iter???")))
        if iterations
        else ["start"]
    )
    for m in models:
        Dm, Rm, vs = evaluate(cfg, project, m)
        if m == "start":  # the residuals of the starting model, for the residual figure
            with open(out_prefix.replace("residuals_final", "residuals_start") + ".csv", "w") as f:
                f.write("src_lat,src_lon,rcv_lat,rcv_lon,period_s,dist_km,t_obs_s,t_pred_s,kept_by_rule\n")
                for d, r in zip(Dm, Rm):
                    f.write(
                        "%.6f,%.6f,%.6f,%.6f,%.4f,%.5f,%.5f,%.5f,%d\n"
                        % (d[0], d[1], d[2], d[3], d[4], r[0], r[2], r[1], r[5])
                    )
        s = residual_stats(Rm)
        s["model"] = "start" if m == "start" else m[-7:]
        if prev is not None:
            s["dm_rms"] = float(np.sqrt(np.mean((vs - prev) ** 2)))
            s["dm_max"] = float(np.abs(vs - prev).max())
        prev = vs
        iters.append(s)
        print(
            "    %-8s rms %.2f ms (%.3f %%), median %.2f ms, kept by the rejection rule %.1f %%"
            % (s["model"], s["rms_ms"], s["rel_rms_percent"], s["median_ms"], 100 * s["kept_fraction"]),
            flush=True,
        )
    return dict(final=final, iterations=iters)


def read_model(path):
    """DSurfTomo "lon lat depth vs" rows -> lon, lat, depth axes and v[lon, lat, depth]."""
    m = np.loadtxt(path)
    lon, lat, dep = np.unique(m[:, 0]), np.unique(m[:, 1]), np.unique(np.round(m[:, 2], 4))
    v = np.full((lon.size, lat.size, dep.size), np.nan)
    if len(m) != lon.size * lat.size * dep.size:
        raise RuntimeError(
            "%s: %d rows for %d x %d x %d nodes (missing or repeated nodes)"
            % (path, len(m), lon.size, lat.size, dep.size)
        )
    v[np.searchsorted(lon, m[:, 0]), np.searchsorted(lat, m[:, 1]), np.searchsorted(dep, np.round(m[:, 2], 4))] = m[
        :, 3
    ]
    if np.isnan(v).any():
        raise RuntimeError("%s: %d nodes have no row (repeated rows in their place)" % (path, int(np.isnan(v).sum())))
    return lon, lat, dep, v


def checker_recovery(project):
    """r, slope, sd of the recovered against the true anomaly per depth, at the nodes with DWS > 0 in the last
    iteration (nonzero data-Jacobian column support; DWS.dat of the solver).  The maps are drawn without this selection."""
    from .solver import DWS_FILE

    syn = np.loadtxt(os.path.join(project, "DSurfTomo.inSyn.dat"))
    tru = np.loadtxt(os.path.join(project, "Vs_model.real"))
    dws = np.loadtxt(os.path.join(project, DWS_FILE), ndmin=2)
    if dws.shape[0] != syn.shape[0] or not np.allclose(dws[:, :3], syn[:, :3], atol=1e-4):
        raise RuntimeError("%s: %s does not hold the nodes of the recovered model" % (project, DWS_FILE))
    M = open(os.path.join(project, "MOD")).read().split("\n")
    zmod = [round(float(v), 3) for v in M[0].split()]
    nrow = len([l for l in M[1:] if l.strip()]) // len(zmod)
    bg = {z: float(M[1 + k * nrow].split()[0]) for k, z in enumerate(zmod)}
    out = []
    for z in sorted(set(np.round(syn[:, 2], 3))):
        ms, mt = np.round(syn[:, 2], 3) == z, np.round(tru[:, 2], 3) == z
        keep = dws[ms, 3] > 0
        ar, at = syn[ms, 3][keep] - bg[z], tru[mt, 3][keep] - bg[z]
        if not keep.any() or np.std(at) < 1e-6:
            continue
        out.append(
            dict(
                depth_km=float(z),
                r=float(np.corrcoef(ar, at)[0, 1]),
                slope=float(np.dot(ar, at) / np.dot(at, at)),
                sd=float(np.std(ar) / np.std(at)),
                n=int(keep.sum()),
                n_all=int(ms.sum()),
            )
        )
    return out
