"""Running DSurfTomo / pDSurfTomo projects, and the forward evaluation of a model on a fixed data set.

A project counts as finished only if its log ends with "Program finishes successfully", its model file is complete
(one finite row per output node) and none of its inputs is newer than the log; the output of an interrupted run is never overwritten silently (status FAILED, move it away).
pDSurfTomo (Zhu et al. 2026) runs in its Mode 2 (parallel travel times and kernels, the original LSMR), which gives
the same model as DSurfTomo.
"""
import glob
import os
import shutil
import subprocess
import tempfile
from ..processes import run as run_process, workers
from ..progress import track

import numpy as np

from ..tools import solver_binary

INPUTS = ("DSurfTomo.in", "MOD", "MOD.true", "surfdata.dat")


def model_complete(d):
    """The model file exists, has one row per output node ((nx - 2) x (ny - 2) x (nz - 1): DSurfTomo writes neither
    the outer ring of the grid nor the fixed bottom node) and only finite, positive velocities."""
    import json
    synthetic = os.path.exists(os.path.join(d, "MOD.true"))          # checkerboards write the recovered model here
    p = os.path.join(d, "DSurfTomo.inSyn.dat" if synthetic else "DSurfTomo.inMeasure.dat")
    meta = os.path.join(d, "project.json")
    if not os.path.exists(p) or not os.path.getsize(p):
        return False
    try:
        m = np.loadtxt(p, ndmin=2)
    except ValueError:
        return False
    if os.path.exists(meta):
        g = json.load(open(meta))
        try:                                       # DSurfTomo drops the rim and the fixed bottom node
            check_node_order(m, g["grid"]["nx"], g["grid"]["ny"], len(g["depths"]), p)
        except RuntimeError:
            return False
    return m.shape[1] >= 4 and bool(np.all(np.isfinite(m[:, 3]))) and bool(np.all(m[:, 3] > 0))


def model_says_done(d):
    log = os.path.join(d, "DSurfTomo.in.log")
    return os.path.exists(log) and "Program finishes successfully" in open(log, errors="ignore").read()[-2000:]


def status(d):
    """missing | new | stale (an input newer than the log) | done (the log ends with the success line and the model
    file is complete) | unfinished."""
    log = os.path.join(d, "DSurfTomo.in.log")
    if not os.path.exists(os.path.join(d, "DSurfTomo.in")):
        return "missing"
    if not os.path.exists(log) and not os.path.exists(os.path.join(d, "run.log")):
        return "new"
    for f in INPUTS:
        p = os.path.join(d, f)
        if os.path.exists(p) and os.path.exists(log) and os.path.getmtime(p) > os.path.getmtime(log):
            return "stale"
    if os.path.exists(log):
        tail = open(log, errors="ignore").read()[-2000:].strip().splitlines()[-5:]
        if any("Program finishes successfully" in ln for ln in tail) and model_complete(d):
            return "done"
    return "unfinished"


def observed_iteration(project):
    """Read only a bounded log tail; buffered logs may have no iteration information yet."""
    import re
    path = os.path.join(project, "DSurfTomo.in.log")
    try:
        with open(path, "rb") as fh:
            fh.seek(max(0, os.path.getsize(path) - 65536))
            text = fh.read().decode("utf-8", errors="replace")
    except FileNotFoundError:
        return None
    matches = re.findall(r"(\d+)(?:st|nd|rd|th)\s+iteration", text, re.IGNORECASE)
    return int(matches[-1]) if matches else None


def run_one(d, binary, threads):
    st = status(d)
    if st == "done":
        return d, "skipped (finished, inputs unchanged)"
    if st in ("stale", "unfinished"):
        return d, "FAILED: output of an earlier or interrupted run is present (%s); delete the outputs or the directory" % st
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    if os.path.basename(binary) == "pDSurfTomo":
        open(os.path.join(d, "ParallelConfig.in"), "w").write(
            "tt_senK_lsmr = parallel, parallel, default\nThreadNum = %d\nInvResultDir = InvResult\n" % threads)
    progress_errors = []

    def show_iteration():
        from ..progress import detail
        try:
            iteration = observed_iteration(d)
        except OSError as exc:
            if not progress_errors:
                print("  iteration progress unavailable: %s" % exc, flush=True)
                progress_errors.append(str(exc))
            return
        if iteration is not None:
            detail("%s: observed iteration %d" % (os.path.basename(d), iteration))

    with open(os.path.join(d, "run.log"), "w") as fh:
        rc = run_process([binary, "DSurfTomo.in"], cwd=d, env=env, stdout=fh,
                         stderr=subprocess.STDOUT, on_poll=show_iteration).returncode
    if rc != 0:
        return d, "FAILED: exit status %d (see %s/run.log)" % (rc, d)
    if status(d) != "done":
        return d, ("FAILED: the model file is incomplete (rows, NaN or non-positive velocities)" if model_says_done(d)
                   else "FAILED: no 'Program finishes successfully' at the end of the log")
    return d, "ok"


def run_projects(cfg, dirs):
    inv = cfg["inversion"]
    binary = solver_binary(inv["solver"])
    with workers(max_workers=int(inv["njobs"])) as ex:
        from concurrent.futures import as_completed
        futures = [ex.submit(run_one, d, binary, int(inv["threads"])) for d in dirs]
        done = [future.result() for future in track(as_completed(futures), len(futures), "solver projects", eta=False)]
    res = sorted(done, key=lambda r: dirs.index(r[0]))  # summaries in project order, not in order of completion
    for d, msg in res:
        print("  %s %s" % (d, msg), flush=True)
    failed = [d for d, m in res if m.startswith("FAILED")]
    return dict(results=dict(res), failed=failed)


# ------------------------------------------------------------------------------------------------ forward eval
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
    ok = (np.ptp(M[..., 2], axis=(1, 2)).max() <= tol and          # depth constant within a depth block
          np.ptp(M[..., 0], axis=(0, 2)).max() <= tol and          # longitude constant along depth and latitude
          np.ptp(M[..., 1], axis=(0, 1)).max() <= tol)             # latitude constant along depth and longitude
    if ok:
        ok = (len(np.unique(np.round(M[:, 0, 0, 2], 4))) == nz - 1 and len(np.unique(np.round(M[0, :, 0, 0], 5))) == ny - 2
              and len(np.unique(np.round(M[0, 0, :, 1], 5))) == nx - 2)
    if not ok:
        raise RuntimeError("%s: the rows are not the inner grid nodes, each once, in DSurfTomo's order (repeated, "
                           "missing or reordered nodes)" % path)


def full_grid(project, model, nx, ny, nz):
    """Starting MOD with the inner nodes (not boundary, not the fixed bottom node) taken from `model`."""
    depz, vs = read_mod(os.path.join(project, "MOD"), nx, ny, nz)
    if model != "start":
        m = np.loadtxt(model)
        if len(m) != (nx - 2) * (ny - 2) * (nz - 1):
            raise RuntimeError("%s: %d rows, expected %d" % (model, len(m), (nx - 2) * (ny - 2) * (nz - 1)))
        check_node_order(m, nx, ny, nz, model)
        vs = vs.copy()
        vs[:nz - 1, 1:ny - 1, 1:nx - 1] = m[:, 3].reshape(nz - 1, ny - 2, nx - 2)
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
    """One DSurfTomo iteration on (depz, vs) in a scratch directory -> rows of residualFirst.dat (distance,
    predicted time, observed time, weighted predicted, weighted observed, weight by the rejection rule).  With
    synthetic=True the project's MOD.true gives the observed times, as in its own run."""
    binary = solver_binary(cfg["inversion"]["solver"])
    threads = threads or int(cfg["inversion"]["threads"])
    L = open(os.path.join(project, "DSurfTomo.in")).read().split("\n")
    ny, nz = vs.shape[1], vs.shape[0]
    tmp = tempfile.mkdtemp(prefix="lant_fwd_")
    try:
        with open(os.path.join(tmp, "MOD"), "w") as f:
            f.write(" ".join("%9.4f" % d for d in depz) + "\n")
            for k in range(nz):
                for j in range(ny):
                    f.write(" ".join("%9.5f" % v for v in vs[k, j]) + "\n")
        shutil.copyfile(os.path.join(project, "surfdata.dat"), os.path.join(tmp, "surfdata.dat"))
        if synthetic:
            shutil.copyfile(os.path.join(project, "MOD.true"), os.path.join(tmp, "MOD.true"))
        L2 = list(L)
        for key, val in (("maximum iteration", "1"), ("synthetic flag", "1" if synthetic else "0")):
            k = next(i for i, l in enumerate(L2) if "c: " + key in l)
            L2[k] = "%-33s" % val + L2[k][L2[k].index("c: "):]
        open(os.path.join(tmp, "DSurfTomo.in"), "w").write("\n".join(L2))
        env = dict(os.environ, OMP_NUM_THREADS=str(threads))
        if os.path.basename(binary) == "pDSurfTomo":
            open(os.path.join(tmp, "ParallelConfig.in"), "w").write(
                "tt_senK_lsmr = parallel, parallel, default\nThreadNum = %d\nInvResultDir = InvResult\n" % threads)
        r = run_process([binary, "DSurfTomo.in"], cwd=tmp, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        res = os.path.join(tmp, "residualFirst.dat")
        if r.returncode != 0 or not os.path.exists(res):
            raise RuntimeError("forward run failed (exit status %d):\n" % r.returncode + r.stdout.decode(errors="ignore")[-2000:])
        R = np.loadtxt(res, ndmin=2)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
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
    except OSError:
        pass
    return out


def evaluate(cfg, project, model="", threads=None):
    """Predicted travel times of a model for every row of the project's data file, by the program's own forward
    operator (one iteration from the model in a scratch directory; residualFirst.dat holds distance, predicted and
    observed time).  DSurfTomo's residualLast.dat belongs to the model of the previous iteration and its log mixes
    rejected rows in, so neither is used."""
    model = model or os.path.join(project, "DSurfTomo.inMeasure.dat")
    L, nx, ny, nz, periods = read_in(os.path.join(project, "DSurfTomo.in"))
    depz, vs = full_grid(project, model, nx, ny, nz)
    R = _forward_once(cfg, project, depz, vs, threads=threads)
    D = read_data(os.path.join(project, "surfdata.dat"), periods)
    if len(D) != len(R):
        raise RuntimeError("data file has %d travel times, the forward run %d" % (len(D), len(R)))
    if not np.allclose(R[:, 0] / D[:, 5], R[:, 2], rtol=2e-3):
        raise RuntimeError("row order check failed (distance / velocity differs from the program's observed time)")
    return D, R, vs


def precheck(cfg, project):
    """Before a project runs: all its data against its starting model, by DSurfTomo's own first iteration
    (checkerboards: against their synthetic data).  DSurfTomo's rejection rule (residual < threshold x q25 or >
    threshold x q75, on the signed quartiles; not changed) drops nearly every datum when the residuals all have
    one sign, e.g. a starting model that is too fast or too slow overall; such a project would 'finish' without
    using any data.  It is invalid when fewer than checks.precheck.min_kept of the data survive the rule, or when
    the median residual is further from 0 than checks.precheck.max_median_iqr x the interquartile range.
    Writes <project>/precheck.json; -> dict(valid, reasons, stats)."""
    import json
    pc = cfg["checks"]["precheck"]
    L, nx, ny, nz, periods = read_in(os.path.join(project, "DSurfTomo.in"))
    depz, vs = read_mod(os.path.join(project, "MOD"), nx, ny, nz)
    R = _forward_once(cfg, project, depz, vs, synthetic=os.path.exists(os.path.join(project, "MOD.true")))
    st = residual_stats(R)
    iqr = st["q75_ms"] - st["q25_ms"]
    reasons = []
    if st["kept_fraction"] < pc["min_kept"]:
        reasons.append("only %.1f %% of the %d travel times survive DSurfTomo's rejection rule (checks.precheck.min_kept "
                       "%.0f %%)" % (100 * st["kept_fraction"], st["n"], 100 * pc["min_kept"]))
    if abs(st["median_ms"]) > pc["max_median_iqr"] * iqr:
        reasons.append("median residual %.2f ms is beyond %.1f x the interquartile range %.2f ms: the starting model is "
                       "too fast or too slow overall" % (st["median_ms"], pc["max_median_iqr"], iqr))
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
    return dict(n=len(r), rms_ms=1e3 * float(np.sqrt(np.mean(r ** 2))), rel_rms_percent=100 * float(np.sqrt(np.mean((r / R[:, 2]) ** 2))),
                mean_ms=1e3 * float(r.mean()), median_ms=1e3 * float(q50), q25_ms=1e3 * float(q25), q75_ms=1e3 * float(q75),
                rejected_by_rule=int((~kept).sum()), kept_fraction=float(kept.mean()),
                kept_median_ms=1e3 * float(np.median(r[kept])) if kept.any() else None,
                zero_predicted=int((R[:, 1] <= 0).sum()))


def residuals(cfg, project, out_prefix, iterations=True):
    """Final-model residuals (csv + stats) and the residual of every iteration's model on the same data."""
    D, R, _ = evaluate(cfg, project)
    with open(out_prefix + ".csv", "w") as f:
        f.write("src_lat,src_lon,rcv_lat,rcv_lon,period_s,dist_km,t_obs_s,t_pred_s,kept_by_rule\n")
        for d, r in zip(D, R):
            f.write("%.6f,%.6f,%.6f,%.6f,%.4f,%.5f,%.5f,%.5f,%d\n" % (d[0], d[1], d[2], d[3], d[4], r[0], r[2], r[1], r[5]))
    final = residual_stats(R)
    iters = []
    prev = None
    models = ["start"] + sorted(glob.glob(os.path.join(project, "DSurfTomo.inMeasure.dat.iter???"))) if iterations else ["start"]
    for m in models:
        Dm, Rm, vs = evaluate(cfg, project, m)
        if m == "start":                          # the residuals of the starting model, for the residual figure
            with open(out_prefix.replace("residuals_final", "residuals_start") + ".csv", "w") as f:
                f.write("src_lat,src_lon,rcv_lat,rcv_lon,period_s,dist_km,t_obs_s,t_pred_s,kept_by_rule\n")
                for d, r in zip(Dm, Rm):
                    f.write("%.6f,%.6f,%.6f,%.6f,%.4f,%.5f,%.5f,%.5f,%d\n" % (d[0], d[1], d[2], d[3], d[4], r[0], r[2], r[1], r[5]))
        s = residual_stats(Rm)
        s["model"] = "start" if m == "start" else m[-7:]
        if prev is not None:
            s["dm_rms"] = float(np.sqrt(np.mean((vs - prev) ** 2)))
            s["dm_max"] = float(np.abs(vs - prev).max())
        prev = vs
        iters.append(s)
        print("    %-8s rms %.2f ms (%.3f %%), median %.2f ms, kept by the rejection rule %.1f %%"
              % (s["model"], s["rms_ms"], s["rel_rms_percent"], s["median_ms"], 100 * s["kept_fraction"]), flush=True)
    return dict(final=final, iterations=iters)


# ------------------------------------------------------------------------------------------------ models
def read_model(path):
    """DSurfTomo "lon lat depth vs" rows -> lon, lat, depth axes and v[lon, lat, depth]."""
    m = np.loadtxt(path)
    lon, lat, dep = np.unique(m[:, 0]), np.unique(m[:, 1]), np.unique(np.round(m[:, 2], 4))
    v = np.full((lon.size, lat.size, dep.size), np.nan)
    if len(m) != lon.size * lat.size * dep.size:
        raise RuntimeError("%s: %d rows for %d x %d x %d nodes (missing or repeated nodes)"
                           % (path, len(m), lon.size, lat.size, dep.size))
    v[np.searchsorted(lon, m[:, 0]), np.searchsorted(lat, m[:, 1]), np.searchsorted(dep, np.round(m[:, 2], 4))] = m[:, 3]
    if np.isnan(v).any():
        raise RuntimeError("%s: %d nodes have no row (repeated rows in their place)" % (path, int(np.isnan(v).sum())))
    return lon, lat, dep, v


def checker_recovery(project, stations, edge_km):
    """r, slope, sd of the recovered against the true anomaly per depth, at nodes within edge_km of a station."""
    syn = np.loadtxt(os.path.join(project, "DSurfTomo.inSyn.dat"))
    tru = np.loadtxt(os.path.join(project, "Vs_model.real"))
    M = open(os.path.join(project, "MOD")).read().split("\n")
    zmod = [round(float(v), 3) for v in M[0].split()]
    nrow = len([l for l in M[1:] if l.strip()]) // len(zmod)
    bg = {z: float(M[1 + k * nrow].split()[0]) for k, z in enumerate(zmod)}
    st = np.array([(s["lon"], s["lat"]) for s in stations.values()])
    out = []
    for z in sorted(set(np.round(syn[:, 2], 3))):
        ms, mt = np.round(syn[:, 2], 3) == z, np.round(tru[:, 2], 3) == z
        lon, lat = syn[ms, 0], syn[ms, 1]
        kx = 111.0 * np.cos(np.radians(lat.mean()))
        d = np.hypot((lon[:, None] - st[None, :, 0]) * kx, (lat[:, None] - st[None, :, 1]) * 111.0).min(axis=1)
        keep = d <= edge_km
        ar, at = syn[ms, 3][keep] - bg[z], tru[mt, 3][keep] - bg[z]
        if np.std(at) < 1e-6:
            continue
        out.append(dict(depth_km=float(z), r=float(np.corrcoef(ar, at)[0, 1]), slope=float(np.dot(ar, at) / np.dot(at, at)),
                        sd=float(np.std(ar) / np.std(at)), n=int(keep.sum())))
    return out
