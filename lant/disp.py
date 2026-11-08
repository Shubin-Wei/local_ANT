"""Two-station group and phase velocity on a log-frequency axis (the measurement of steps 06 and 07).

For every pair of the notched stack:
  * symmetric component (causal + time-reversed acausal) / 2 (Bensen et al. 2007), Gaussian filter bank with
    constant relative width alpha on the frequencies of the reference curve -> noise-normalised envelope images
    (noise window from dist / vmin + 2 periods to the end of the lag window);
  * group ridge by dynamic programming guided by the group velocity implied by the phase reference,
    U = c / (1 - dln c / dln f) (Bensen et al. 2007 eq. 7): search window 0.65-1.35 U_ref, penalty beyond 25 %,
    smoothness 3 % between neighbouring frequencies; causal and acausal branches are picked as well (diagnostic);
  * phase travel time from the analytic phase at the group arrival, t_p = t_g - (psi - phi0) / (2 pi f),
    phi0 = pi/4 (far-field phase of the 2-D Green's function, Lin et al. 2008 eq. 7);
  * cycle count fixed at the lowest frequency of each contiguous band where the path is between lam_c and
    anchor_max reference wavelengths long, then followed upwards on a 48/octave grid with the measured group
    slowness (AFTAN phtovel.f); tracking stops at the first phase slip (> 0.3 cycle);
  * delete-a-group jackknife over the sub-stacks (njack groups) for both velocities;
  * per-frequency QC (thresholds: dispersion.qc); contiguous bands of >= min_run samples only.
The thresholds of the group velocity are applied again afterwards by lant.qc (group_qc).
"""
import os
from multiprocessing import Pool

import h5py
import numpy as np

from . import signal as S

VEL = np.round(np.arange(1.2, 4.001, 0.01), 3)
PHI0 = np.pi / 4
QC = dict(snr_g=8.0, lam_g=2.0, side_g=0.9, width_g=2.0, jk_g_rel=0.05, dev_g=0.25,
          snr_c=6.0, lam_c=2.0, anchor_max=3.5, side_c=1.0, jk_c_rel=0.02, dev_c=0.20, smooth_c=0.03,
          ucons_c=0.25, min_run=6, dev_pre=0.30)
G = {}


def off_band_mask(freqs, spec):
    m = np.zeros(len(freqs), bool)
    if spec == "all":
        m[:] = True
    elif spec:
        for b in spec.split(","):
            lo, hi = (float(x) for x in b.split("-"))
            m |= (freqs >= lo) & (freqs < hi)
    return m


def group_reference(freqs, c_ref):
    return c_ref / (1.0 - np.gradient(np.log(c_ref), np.log(freqs)))


def reference_work(freqs, c_ref, valid):
    """Reference curves where defined; the implied group velocity is differentiated inside each contiguous valid
    segment only.  Returns (valid, u_ref, c_fill, u_fill)."""
    valid = valid & np.isfinite(c_ref)
    u_ref = np.full(len(freqs), np.nan)
    for (k0, k1) in S.runs_of(valid):
        if k1 - k0 < 2:
            valid[k0:k1] = False
            continue
        u_ref[k0:k1] = group_reference(freqs[k0:k1], c_ref[k0:k1])
    if not valid.any():
        raise RuntimeError("reference curve has no valid frequency")
    idx = np.where(valid)[0]
    near = idx[np.argmin(np.abs(np.arange(len(freqs))[:, None] - idx[None, :]), axis=1)]
    return valid, u_ref, c_ref[near], u_ref[near]


def load_reference(path):
    ref = np.load(path)
    out = {k: ref[k] for k in ref.files}
    c = np.asarray(out["c_ref"], float)
    v = np.asarray(out["ref_valid"], bool) if "ref_valid" in out else np.ones(len(c), bool)
    out["ref_valid"], out["u_ref"], out["c_fill"], out["u_fill"] = reference_work(out["freqs"], c, v.copy())
    return out


def init_worker(opt):
    global VEL, PHI0
    VEL = np.round(np.arange(opt["vmin"], opt["vmax"] + 1e-9, opt["dv"]), 3)
    PHI0 = opt["phi0"]
    QC.update(opt["qc"])
    for k, v in opt["qc_off"].items():
        QC[k] = v if v > 0 else np.inf
    QC["lam_g"] = QC["lam_c"] = opt["lam_min"]
    QC["anchor_max"] = max(opt["anchor_max"], opt["lam_min"] + 1.5)
    G.update(opt)
    G["ref"] = load_reference(opt["reference"])
    G["c_phase"] = G["ref"]["c_fill"]
    with h5py.File(opt["stack"], "r") as h:
        G["fs"] = float(h.attrs["fs"])
        G["maxlag"] = float(h.attrs["maxlag"])
        G["dist"] = h["dist_km"][:]
        G["nseg"] = h["nseg"][:]


def smooth_f(x, n):
    k = np.ones(n) / n
    pad = np.pad(x, [(0, 0)] * (x.ndim - 1) + [(n // 2, n // 2)], mode="edge")
    return np.apply_along_axis(lambda v: np.convolve(v, k, mode="valid"), -1, pad)


def branch_images(x2, dist, fs, maxlag, freqs, alpha):
    """[(analytic, snr image, snr)] for the causal, acausal and symmetric component."""
    pos, neg = S.branches(x2, fs, maxlag)
    out = []
    t = np.arange(pos.shape[1]) / fs
    for s in (pos, neg, 0.5 * (pos + neg)):
        an = S.analytic_bank(s, fs, freqs, alpha)
        env = np.abs(an)
        times = np.broadcast_to(dist[:, None, None] / VEL[None, None, :], (len(dist), len(freqs), len(VEL))).copy()
        img = S.sample_at_times(env, fs, times)
        nrm = np.empty((len(dist), len(freqs)))
        for q in range(len(dist)):
            for k, fc in enumerate(freqs):
                m = t > min(dist[q] / VEL[0] + 2.0 / fc, maxlag - 3.0)
                nrm[q, k] = np.sqrt(np.mean(np.real(an[q, k, m]) ** 2)) + 1e-30
        simg = img / nrm[:, :, None]
        out.append((an, simg, simg.max(axis=2)))
    return out


def lobe_metrics(prof, j):
    """Height of the strongest distinct competing peak relative to the pick (valley below 70 % of the lower of
    the two), and the FWHM of the picked lobe."""
    p = prof
    lo = j
    while lo > 0 and p[lo - 1] <= p[lo]:
        lo -= 1
    hi = j
    while hi < len(p) - 1 and p[hi + 1] <= p[hi]:
        hi += 1
    pj = max(p[j], 1e-30)
    side = 0.0
    maxima = np.where((p[1:-1] >= p[:-2]) & (p[1:-1] >= p[2:]) & (p[1:-1] > 0))[0] + 1
    for m in maxima:
        if m == j:
            continue
        a, b = (m, j) if m < j else (j, m)
        if p[a:b + 1].min() <= 0.7 * min(p[m], pj):
            side = max(side, p[m] / pj)
    half = 0.5 * p[j]
    a = j
    while a > lo and p[a] > half:
        a -= 1
    b = j
    while b < hi and p[b] > half:
        b += 1
    width = VEL[b] - VEL[a]
    if p[a] > half or p[b] > half:
        width = VEL[hi] - VEL[lo]
    return side, width


def theoretical_width(v, dist, freqs, alpha):
    t0 = dist / v
    h = np.sqrt(alpha * np.log(2)) / (np.pi * freqs)
    vhi = np.where(t0 > h, dist / np.maximum(t0 - h, 1e-6), VEL[-1])
    vhi = np.minimum(vhi, VEL[-1])
    vlo = np.maximum(dist / (t0 + h), VEL[0])
    return vhi - vlo


def pick_group(img, u_ref, valid, segments=None):
    if segments is not None and not (len(segments) == 1 and segments[0] == (0, len(u_ref))):
        vg, side, width = (np.full(len(u_ref), np.nan) for _ in range(3))
        for (k0, k1) in segments:
            if valid[k0:k1].any():
                vg[k0:k1], side[k0:k1], width[k0:k1] = pick_group(img[k0:k1], u_ref[k0:k1], valid[k0:k1])
        return vg, side, width
    prof = img / (img.max(axis=1, keepdims=True) + 1e-30)
    inwin = (VEL[None, :] >= 0.65 * u_ref[:, None]) & (VEL[None, :] <= 1.35 * u_ref[:, None])
    prof_dp = np.where(inwin, prof, 1e-3)
    path = S.ridge_dp(prof_dp, VEL, prior=u_ref, prior_tol=0.25, smooth=0.03, valid=valid)
    nf = len(path)
    vg, side, width = (np.full(nf, np.nan) for _ in range(3))
    plaus = (VEL[None, :] >= 0.6 * u_ref[:, None]) & (VEL[None, :] <= 1.45 * u_ref[:, None])
    for k in range(nf):
        if not valid[k]:
            continue
        j = S.local_peak(prof_dp[k], path[k], 3)
        jr = S.refine_peak(prof[k], j)
        vg[k] = np.interp(jr, np.arange(len(VEL)), VEL)
        side[k], width[k] = lobe_metrics(np.where(plaus[k], prof[k], 0.0), j)
    return vg, side, width


def phase_travel_times(an_sel, fs, ffine, tg_fine, phi0):
    psi = np.angle(S.sample_at_times(an_sel, fs, tg_fine[:, None])[:, 0])
    return tg_fine - (psi - phi0) / (2 * np.pi * ffine)


def track_phase(tp0, dist, ffine, c_ref_f, start, stop, slip=0.3, tg_f=None):
    """Cycle count at `start` from c_ref, then the phase followed upward with the measured group slowness
    (w(k) tp(k) = w(k-1) tp(k-1) + (tg(k) + tg(k-1)) / 2 (w(k) - w(k-1)); AFTAN phtovel.f) until a slip."""
    T = 1.0 / ffine
    om = 2 * np.pi * ffine
    nf = len(ffine)
    tp = np.full(nf, np.nan)
    m = np.round((dist / c_ref_f[start] - tp0[start]) / T[start])
    tp[start] = tp0[start] + m * T[start]
    for k in range(start + 1, stop):
        if tg_f is None:
            pred = tp[k - 1] + (dist / c_ref_f[k] - dist / c_ref_f[k - 1])
        else:
            pred = (om[k - 1] * tp[k - 1] + 0.5 * (tg_f[k] + tg_f[k - 1]) * (om[k] - om[k - 1])) / om[k]
        m = np.round((pred - tp0[k]) / T[k])
        cand = tp0[k] + m * T[k]
        if abs(cand - pred) > slip * T[k]:
            break
        tp[k] = cand
    with np.errstate(divide="ignore", invalid="ignore"):
        return dist / tp


def process_chunk(idx):
    fs, maxlag, alpha = G["fs"], G["maxlag"], G["alpha"]
    ref = G["ref"]
    freqs = ref["freqs"]
    ref_valid = ref["ref_valid"]
    c_ref = ref["c_fill"]
    u_ref = ref["u_fill"]
    segs = S.runs_of(ref_valid)
    side_off = off_band_mask(freqs, G["side_off_bands"])
    phi0 = PHI0
    nf = len(freqs)
    with h5py.File(G["stack"], "r") as h:
        x2 = h["stack"][idx].astype(np.float64)
        daily = h["daily"][:, idx].astype(np.float64)
        dcount = h["daily_count"][:, idx].astype(np.float64)
    dist = G["dist"][idx]
    npair = len(idx)
    bim = branch_images(x2, dist, fs, maxlag, freqs, alpha)
    snr_b = np.stack([bim[0][2], bim[1][2]])
    w_b = smooth_f(snr_b, 7) ** 2
    img_sym = bim[2][1]
    snr_sym = bim[2][2]
    keys = ("vg", "vc", "snr", "side", "width_ratio", "jk_g", "jk_c", "u_from_c", "vg_pos", "vg_neg",
            "snr_pos", "snr_neg", "snr_c")
    res = {k: np.full((npair, nf), np.nan, np.float32) for k in keys}
    res["branch"] = np.full((npair, nf), -1, np.int8)
    res["w_pos"] = (w_b[0] / (w_b.sum(0) + 1e-30)).astype(np.float32)
    res["ok_g"] = np.zeros((npair, nf), bool)
    res["ok_c"] = np.zeros((npair, nf), bool)
    res["pre_c"] = np.zeros((npair, nf), bool)
    res["cycle"] = np.full((npair, nf), -32768, np.int16)
    res["anchor_f"] = np.full((npair, nf), np.nan, np.float32)
    img_keep = np.zeros((npair, nf, len(VEL)), np.float16)
    nday = daily.shape[0]
    groups = [np.arange(nday)[g::G["njack"]] for g in range(G["njack"])]
    tot_w = dcount.sum(axis=0)
    tot_sum = (daily * dcount[:, :, None]).sum(axis=0)
    po = int(round(1.0 / np.log2(freqs[1] / freqs[0])))          # points per octave of the measurement grid
    ffine = S.log_freqs(freqs[0], freqs[-1], 4 * po)              # phase tracked on a 4x finer grid
    lff = np.log(ffine)
    lf = np.log(freqs)
    c_phase_f = np.interp(lff, lf, G["c_phase"])
    c_ref_f = np.interp(lff, lf, c_ref)
    kc = np.clip(np.round(np.log2(ffine / freqs[0]) * po).astype(int), 0, nf - 1)
    for q in range(npair):
        if G["nseg"][idx[q]] == 0 or dist[q] < G["min_dist"]:
            continue
        lam_g = u_ref / freqs
        lam_c = c_ref / freqs
        img = img_sym[q]
        snr = snr_sym[q]
        valid = (dist[q] >= lam_c) & (snr >= 4.0) & ref_valid
        if valid.sum() < QC["min_run"]:
            continue
        vg, side, width = pick_group(img, u_ref, valid, segs)
        for b, key in ((0, "vg_pos"), (1, "vg_neg")):
            res[key][q] = pick_group(bim[b][1][q], u_ref, valid, segs)[0]
        res["snr_pos"][q], res["snr_neg"][q] = snr_b[0, q], snr_b[1, q]
        wratio = width / np.maximum(theoretical_width(vg, dist[q], freqs, alpha), 1e-3)
        img_keep[q] = (img / (img.max(axis=1, keepdims=True) + 1e-30)).astype(np.float16)
        subsets = []
        for g in groups:
            w = tot_w[q] - dcount[g, q].sum()
            if w <= 0 or dcount[g, q].sum() <= 0:
                continue
            subsets.append((tot_sum[q] - (daily[g, q] * dcount[g, q][:, None]).sum(axis=0)) / w)
        jk_v = []
        for sub in subsets:
            sb = branch_images(sub[None, :], dist[q:q + 1], fs, maxlag, freqs, alpha)
            im = sb[2][1][0]
            vk = np.full(nf, np.nan)
            for k in range(nf):
                if not np.isfinite(vg[k]):
                    continue
                j0 = int(np.argmin(np.abs(VEL - vg[k])))
                hw = max(2, int(round(0.10 * vg[k] / 0.01)))
                j = S.local_peak(im[k], j0, hw)
                vk[k] = np.interp(S.refine_peak(im[k], j), np.arange(len(VEL)), VEL)
            jk_v.append(vk)
        K = len(jk_v)
        if K >= G["min_jackknife_groups"]:
            jk_v = np.array(jk_v)
            jk_g = np.sqrt((K - 1) / K * np.sum((jk_v - jk_v.mean(0)) ** 2, axis=0))
        else:
            jk_g = np.full(nf, np.inf)
        dev_g = np.abs(vg - u_ref) / u_ref
        ok_g = (snr >= QC["snr_g"]) & (dist[q] >= QC["lam_g"] * lam_g) & ((side <= QC["side_g"]) | side_off) \
            & (wratio <= QC["width_g"]) & (jk_g <= QC["jk_g_rel"] * vg) & (dev_g <= QC["dev_g"]) & ref_valid
        ok_g = S.contiguous(ok_g, QC["min_run"])
        res["vg"][q], res["snr"][q], res["side"][q] = vg, snr, side
        res["width_ratio"][q], res["jk_g"][q], res["ok_g"][q] = wratio, jk_g, ok_g
        pre = (snr >= QC["snr_c"]) & (dist[q] >= QC["lam_c"] * lam_c) & (dev_g <= QC["dev_pre"]) \
            & ((side <= QC["side_c"]) | side_off) & ref_valid & np.isfinite(vg)
        pre = S.contiguous(pre, QC["min_run"])
        res["pre_c"][q] = pre
        if not pre.any():
            continue
        vg_fill = np.where(np.isfinite(vg), vg, u_ref)
        tg_f = dist[q] / np.interp(lff, lf, vg_fill)
        lam_ratio_f = dist[q] * ffine / c_ref_f
        sigsets = [x2[q:q + 1]] + [sub[None, :] for sub in subsets]
        an_sets = []
        for sig2 in sigsets:
            p, n = S.branches(sig2, fs, maxlag)
            an_sets.append(S.analytic_bank(0.5 * (p + n), fs, ffine, alpha)[0])
        c_fine_all = np.full((len(sigsets), len(ffine)), np.nan)
        m_fine = np.full(len(ffine), np.nan)
        anchor_run = np.full(nf, np.nan)
        branch_run = np.full(nf, -1)
        for (k0, k1) in S.runs_of(pre):
            branch_run[k0:k1] = 2
            fine = np.where((kc >= k0) & (kc < k1))[0]
            anchor = fine[(lam_ratio_f[fine] >= QC["lam_c"]) & (lam_ratio_f[fine] <= QC["anchor_max"])]
            if anchor.size == 0:
                continue
            st, sp = int(anchor[0]), int(fine[-1]) + 1
            anchor_run[k0:k1] = ffine[st]
            for s_i, an_s in enumerate(an_sets):
                tp0 = phase_travel_times(an_s, fs, ffine, tg_f, phi0)
                c_fine_all[s_i, st:sp] = track_phase(tp0, dist[q], ffine, c_phase_f, st, sp, tg_f=tg_f)[st:sp]
                if s_i == 0:
                    with np.errstate(invalid="ignore", divide="ignore"):
                        m_fine[st:sp] = np.round((dist[q] / c_fine_all[0, st:sp] - tp0[st:sp]) * ffine[st:sp])
        res["branch"][q] = branch_run
        res["anchor_f"][q] = anchor_run
        for k in range(nf):
            fk = np.where(kc == k)[0]
            if fk.size and np.all(np.isfinite(m_fine[fk])):
                res["cycle"][q, k] = int(m_fine[fk[np.argmin(np.abs(lff[fk] - lf[k]))]])
        snr_c_q = np.full(nf, np.nan)
        for k in range(nf):
            b = branch_run[k]
            if b >= 0 and np.isfinite(vg[k]):
                snr_c_q[k] = np.interp(vg[k], VEL, bim[b][1][q, k])
        res["snr_c"][q] = snr_c_q
        c_fine = c_fine_all[0]
        jc = c_fine_all[1:]
        if jc.shape[0] >= G["min_jackknife_groups"]:
            Kc = jc.shape[0]
            with np.errstate(invalid="ignore"):
                jk_c_f = np.sqrt((Kc - 1) / Kc * np.sum((jc - jc.mean(0)) ** 2, axis=0))
        else:
            jk_c_f = np.full(len(ffine), np.inf)
        vc = np.full(nf, np.nan)
        jk_c = np.full(nf, np.nan)
        for k in range(nf):
            fk = kc == k
            if fk.any() and np.all(np.isfinite(c_fine[fk])) and np.all(np.isfinite(jk_c_f[fk])):
                vc[k] = np.interp(lf[k], lff[fk], c_fine[fk]) if fk.sum() > 1 else c_fine[fk][0]
                jk_c[k] = np.max(jk_c_f[fk])
        with np.errstate(divide="ignore", invalid="ignore"):
            dlnc = np.gradient(np.log(vc), lf)
            u_from_c = vc / (1.0 - dlnc)
            step = np.abs(np.diff(np.log(vc)))
            smooth_ok = np.r_[step <= QC["smooth_c"], True] & np.r_[True, step <= QC["smooth_c"]]
            ucons = np.abs(u_from_c - vg) / vg <= QC["ucons_c"]
            ok_c = pre & np.isfinite(vc) & (jk_c <= QC["jk_c_rel"] * vc) \
                & (np.abs(vc - c_ref) / c_ref <= QC["dev_c"]) & smooth_ok & ucons \
                & np.isfinite(snr_c_q) & (snr_c_q >= QC["snr_c"])
        ok_c = S.contiguous(ok_c, QC["min_run"])
        res["vc"][q], res["jk_c"][q], res["ok_c"][q], res["u_from_c"][q] = vc, jk_c, ok_c, u_from_c
    return idx, res, img_keep


def auto_njack(stack_path):
    """Lines: the number of sub-stacks holding at least 25 % of the segments of the fullest one."""
    with h5py.File(stack_path, "r") as h:
        c = h["daily_count"][:].sum(1)
    return int((c >= 0.25 * c.max()).sum())


def measure(cfg, stack, reference, out, first_pass=False, save_images=0):
    """Run the measurement on every pair with data -> out (HDF5)."""
    from .config import phase_constant
    d = cfg["dispersion"]
    njack = auto_njack(stack) if d["njack"] == "auto" else int(d["njack"])
    if njack < 2:
        raise RuntimeError("only %d sub-stacks hold data: the jackknife needs at least 2 (stack.block_hours)" % njack)
    q = dict(d["qc"])
    lam_c_final = q.pop("lam_c")
    q.pop("min_run")
    qc_off = dict(dev_g=d["group_qc"]["dev_g"], dev_c=q.pop("dev_c"), dev_pre=q.pop("dev_pre"), ucons_c=q.pop("ucons_c"))
    if first_pass:
        qc_off.update(dev_g=0.0, dev_c=0.0, dev_pre=0.0)
    opt = dict(stack=stack, reference=reference, alpha=d["alpha"], njack=njack,
               min_jackknife_groups=min(d["min_jackknife_groups"], njack), lam_min=d["lam_min"],
               anchor_max=d["anchor_max_lambda"], side_off_bands=d["side_off_bands"], min_dist=d["min_dist_km"],
               vmin=d["vmin"], vmax=d["vmax"], dv=d["dv"], phi0=phase_constant(cfg),
               qc=dict(snr_c=q["snr_c"], jk_c_rel=q["jk_c_rel"], smooth_c=q["smooth_c"], side_c=q["side_c"],
                       min_run=d["qc"]["min_run"]), qc_off=qc_off)
    with h5py.File(stack, "r") as h:
        pairs = np.array([p.decode() for p in h["pairs"][:]])
        nseg = h["nseg"][:]
        dist = h["dist_km"][:]
    all_idx = np.where(nseg > 0)[0]
    init_worker(opt)
    ref = G["ref"]
    freqs = ref["freqs"]
    nf = len(freqs)
    chunks = [all_idx[i:i + d["chunk"]] for i in range(0, len(all_idx), d["chunk"])]
    rng = np.random.default_rng(1)
    keep_img = set(rng.choice(all_idx, min(save_images, len(all_idx)), replace=False).tolist()) if save_images else set()
    npair_all = len(pairs)
    out_d = {k: np.full((npair_all, nf), np.nan, np.float32) for k in
             ("vg", "vc", "snr", "side", "width_ratio", "jk_g", "jk_c", "u_from_c", "w_pos", "vg_pos", "vg_neg",
              "snr_pos", "snr_neg", "snr_c")}
    out_d["branch"] = np.full((npair_all, nf), -1, np.int8)
    out_d["ok_g"] = np.zeros((npair_all, nf), bool)
    out_d["ok_c"] = np.zeros((npair_all, nf), bool)
    out_d["pre_c"] = np.zeros((npair_all, nf), bool)
    out_d["cycle"] = np.full((npair_all, nf), -32768, np.int16)
    out_d["anchor_f"] = np.full((npair_all, nf), np.nan, np.float32)
    img_idx, img_data = [], []
    os.environ["OMP_NUM_THREADS"] = "1"
    with Pool(int(d["nproc"]), initializer=init_worker, initargs=(opt,)) as pool:
        done = 0
        for idx, res, imgs in pool.imap_unordered(process_chunk, chunks):
            for k in out_d:
                out_d[k][idx] = res[k]
            for j, i in enumerate(idx):
                if i in keep_img:
                    img_idx.append(i)
                    img_data.append(imgs[j])
            done += len(idx)
            if done % (20 * d["chunk"]) < d["chunk"] or done == len(all_idx):
                print("  pairs done %d / %d" % (done, len(all_idx)), flush=True)
    with h5py.File(out + ".tmp", "w") as h:
        for k, v in out_d.items():
            h.create_dataset(k, data=v)
        h.create_dataset("freqs", data=freqs)
        h.create_dataset("u_ref", data=ref["u_ref"])
        h.create_dataset("c_ref", data=np.where(ref["ref_valid"], ref["c_fill"], np.nan))
        h.create_dataset("ref_valid", data=ref["ref_valid"])
        h.create_dataset("side_off", data=off_band_mask(freqs, d["side_off_bands"]))
        h.create_dataset("vel_axis", data=VEL)
        h.create_dataset("pairs", data=pairs.astype("S"))
        h.create_dataset("dist_km", data=dist)
        order = np.argsort(img_idx)
        h.create_dataset("img_index", data=np.array(img_idx, int)[order])
        h.create_dataset("img", data=np.array(img_data)[order] if img_data else np.zeros((0, nf, len(VEL)), np.float16))
        h.attrs.update(dict(alpha=d["alpha"], njack=njack, min_jackknife_groups=opt["min_jackknife_groups"],
                            reference=os.path.abspath(reference), first_pass=bool(first_pass),
                            side_off_bands=d["side_off_bands"] or "none", phase_constant=opt["phi0"],
                            lam_c_final=lam_c_final,
                            **{"qc_" + k: (0.0 if np.isinf(v) else v) for k, v in QC.items()}))
    os.replace(out + ".tmp", out)
    return dict(pairs=len(all_idx), njack=njack, phase_points=int(out_d["ok_c"].sum()),
                group_points=int(out_d["ok_g"].sum()))
