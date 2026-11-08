"""Running DSurfTomo / pDSurfTomo projects, and the forward evaluation of a model on a fixed data set.

A project counts as finished only if its log ends with "Program finishes successfully" and none of its inputs is
newer than the log; the output of an interrupted run is never overwritten silently (status FAILED, move it away).
pDSurfTomo (Zhu et al. 2026) runs in its Mode 2 (parallel travel times and kernels, the original LSMR), which gives
the same model as DSurfTomo.
"""
import glob
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from ..tools import solver_binary

INPUTS = ("DSurfTomo.in", "MOD", "MOD.true", "surfdata.dat")


def status(d):
    log = os.path.join(d, "DSurfTomo.in.log")
    if not os.path.exists(os.path.join(d, "DSurfTomo.in")):
        return "missing"
    if not os.path.exists(log) and not os.path.exists(os.path.join(d, "run.log")):
        return "new"
    for f in INPUTS:
        p = os.path.join(d, f)
        if os.path.exists(p) and os.path.exists(log) and os.path.getmtime(p) > os.path.getmtime(log):
            return "stale"
    if os.path.exists(log) and "Program finishes successfully" in open(log, errors="ignore").read():
        return "done"
    return "unfinished"


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
    with open(os.path.join(d, "run.log"), "w") as fh:
        rc = subprocess.run([binary, "DSurfTomo.in"], cwd=d, env=env, stdout=fh, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        return d, "FAILED: exit status %d (see %s/run.log)" % (rc, d)
    if status(d) != "done":
        return d, "FAILED: no 'Program finishes successfully' in the log"
    return d, "ok"


def run_projects(cfg, dirs):
    inv = cfg["inversion"]
    binary = solver_binary(inv["solver"])
    with ThreadPoolExecutor(max_workers=int(inv["njobs"])) as ex:
        res = list(ex.map(lambda d: run_one(d, binary, int(inv["threads"])), dirs))
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


def full_grid(project, model, nx, ny, nz):
    """Starting MOD with the inner nodes (not boundary, not the fixed bottom node) taken from `model`."""
    depz, vs = read_mod(os.path.join(project, "MOD"), nx, ny, nz)
    if model != "start":
        m = np.loadtxt(model)
        if len(m) != (nx - 2) * (ny - 2) * (nz - 1):
            raise RuntimeError("%s: %d rows, expected %d" % (model, len(m), (nx - 2) * (ny - 2) * (nz - 1)))
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


def evaluate(cfg, project, model="", threads=None):
    """Predicted travel times of a model for every row of the project's data file, by the program's own forward
    operator (one iteration from the model in a scratch directory; residualFirst.dat holds distance, predicted and
    observed time).  DSurfTomo's residualLast.dat belongs to the model of the previous iteration and its log mixes
    rejected rows in, so neither is used."""
    binary = solver_binary(cfg["inversion"]["solver"])
    threads = threads or int(cfg["inversion"]["threads"])
    model = model or os.path.join(project, "DSurfTomo.inMeasure.dat")
    L, nx, ny, nz, periods = read_in(os.path.join(project, "DSurfTomo.in"))
    depz, vs = full_grid(project, model, nx, ny, nz)
    tmp = tempfile.mkdtemp(prefix="lant_fwd_")
    try:
        with open(os.path.join(tmp, "MOD"), "w") as f:
            f.write(" ".join("%9.4f" % d for d in depz) + "\n")
            for k in range(nz):
                for j in range(ny):
                    f.write(" ".join("%9.5f" % v for v in vs[k, j]) + "\n")
        shutil.copyfile(os.path.join(project, "surfdata.dat"), os.path.join(tmp, "surfdata.dat"))
        L2 = list(L)
        for key, val in (("maximum iteration", "1"), ("synthetic flag", "0")):
            k = next(i for i, l in enumerate(L2) if "c: " + key in l)
            L2[k] = "%-33s" % val + L2[k][L2[k].index("c: "):]
        open(os.path.join(tmp, "DSurfTomo.in"), "w").write("\n".join(L2))
        env = dict(os.environ, OMP_NUM_THREADS=str(threads))
        if os.path.basename(binary) == "pDSurfTomo":
            open(os.path.join(tmp, "ParallelConfig.in"), "w").write(
                "tt_senK_lsmr = parallel, parallel, default\nThreadNum = %d\nInvResultDir = InvResult\n" % threads)
        r = subprocess.run([binary, "DSurfTomo.in"], cwd=tmp, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        res = os.path.join(tmp, "residualFirst.dat")
        if not os.path.exists(res):
            raise RuntimeError("forward run failed:\n" + r.stdout.decode(errors="ignore")[-2000:])
        R = np.loadtxt(res, ndmin=2)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    D = read_data(os.path.join(project, "surfdata.dat"), periods)
    if len(D) != len(R):
        raise RuntimeError("data file has %d travel times, the forward run %d" % (len(D), len(R)))
    if not np.allclose(R[:, 0] / D[:, 5], R[:, 2], rtol=2e-3):
        raise RuntimeError("row order check failed (distance / velocity differs from the program's observed time)")
    return D, R, vs


def residual_stats(R):
    r = R[:, 2] - R[:, 1]
    return dict(n=len(r), rms_ms=1e3 * float(np.sqrt(np.mean(r ** 2))), rel_rms_percent=100 * float(np.sqrt(np.mean((r / R[:, 2]) ** 2))),
                mean_ms=1e3 * float(r.mean()), rejected_by_rule=int((R[:, 5] == 0).sum()), zero_predicted=int((R[:, 1] <= 0).sum()))


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
        _, Rm, vs = evaluate(cfg, project, m)
        s = residual_stats(Rm)
        s["model"] = "start" if m == "start" else m[-7:]
        if prev is not None:
            s["dm_rms"] = float(np.sqrt(np.mean((vs - prev) ** 2)))
            s["dm_max"] = float(np.abs(vs - prev).max())
        prev = vs
        iters.append(s)
        print("    %-8s rms %.2f ms (%.3f %%)" % (s["model"], s["rms_ms"], s["rel_rms_percent"]), flush=True)
    return dict(final=final, iterations=iters)


# ------------------------------------------------------------------------------------------------ models
def read_model(path):
    """DSurfTomo "lon lat depth vs" rows -> lon, lat, depth axes and v[lon, lat, depth]."""
    m = np.loadtxt(path)
    lon, lat, dep = np.unique(m[:, 0]), np.unique(m[:, 1]), np.unique(np.round(m[:, 2], 4))
    v = np.full((lon.size, lat.size, dep.size), np.nan)
    v[np.searchsorted(lon, m[:, 0]), np.searchsorted(lat, m[:, 1]), np.searchsorted(dep, np.round(m[:, 2], 4))] = m[:, 3]
    return lon, lat, dep, v


def checker_recovery(project, stations, edge_km=0.8):
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
