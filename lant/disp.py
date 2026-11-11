"""Two-station group and phase velocity on a log-frequency axis (the measurement of steps 06 and 07).

For every pair of the stack:
  * symmetric component (causal + time-reversed acausal) / 2 (Bensen et al. 2007), Gaussian filter bank with
    constant relative width alpha on the frequencies of the reference curve -> noise-normalised envelope images
    (noise window from dist / vmin + 2 periods, or the last 15 % of the lag window, to its end);
  * group ridge by dynamic programming guided by the group velocity implied by the phase reference,
    U = c / (1 - dln c / dln f) (Bensen et al. 2007 eq. 7): search window 0.65-1.35 U_ref, penalty beyond 25 %,
    smoothness 3 % between neighbouring frequencies;
  * phase from the phase-matched filter (Herrin & Goforth 1977; AFTAN, Levshin & Ritzwoller 2001): the spectrum
    is multiplied by exp(i Psi), Psi(w) = integral of the measured group delay t_g dw, which removes the dispersion
    of the arrival, and the Gaussian-weighted sum around each frequency gives the spectral phase there,
    Phi(wc) = arg sum X(w) G(w) exp(i [Psi(w) - Psi(wc)]),  t_p = -(Phi - phi0) / wc  (matched_phase).
    Reading the phase of the narrow-band signal at the group arrival is the special case Psi = w t_g of this sum;
    it is biased when the group delay changes inside the filter band (the systematic error of Shapiro & Singh 1999
    for group velocities has the same origin).  phi0 = pi/4: far-field phase of the 2-D Green's function (Lin et
    al. 2008 eq. 7);
  * a frequency is not measured where the amplitude spectrum of the coherent signal (all pairs together; signal
    window minus the later lags, so that it does not depend on the noise level) at that frequency is below
    qc.band_edge (0.5 = -6 dB) of its mean over the filter band: the filter then passes mostly energy of other
    frequencies (the edge of the band the noise sources excite; band_amplitude);
  * cycle count fixed at the lowest frequency of each contiguous band where the path is between lam_c and
    anchor_max reference wavelengths long, then followed upwards on a 48/octave grid with the measured group
    slowness (AFTAN phtovel.f); tracking stops at the first phase slip (> 0.3 cycle);
  * delete-a-group jackknife over the sub-stacks (njack groups) for the phase velocity;
  * per-frequency QC (every threshold in dispersion.qc; none hidden); contiguous bands of >= min_run samples only;
    the reason of every rejection is kept as a bit mask ("reason", lant.qc.REASONS).
The group velocity is measured only to find the arrival the phase is read at (not a product).
"""

import os
from multiprocessing import Pool

from .progress import track

import h5py
import numpy as np

from . import signal as S

VEL = None  # velocity axis, set from dispersion.vmin/vmax/dv by init_worker
PHI0 = np.pi / 4
QC = {}  # dispersion.qc of the configuration, set by init_worker
G = {}
DATA_KEYS = ("vg", "vc", "snr", "side", "jk_c", "u_from_c")  # float arrays [npair, nf] of the measurement file


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
    for k0, k1 in S.runs_of(valid):
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
    sup = np.asarray(out["ref_support"], bool) if "ref_support" in out else out["ref_valid"].copy()
    out["ref_support"] = sup & out["ref_valid"]  # backed by measurements (F-J picks, medians, a borrowed curve)
    return out


def init_worker(opt):
    global VEL, PHI0
    VEL = np.round(np.arange(opt["vmin"], opt["vmax"] + opt["dv"] / 2, opt["dv"]), 6)
    PHI0 = opt["phi0"]
    QC.update(opt["qc"])
    for k, v in opt["qc_off"].items():
        QC[k] = v if v > 0 else np.inf
    QC["lam_c"] = opt["lam_min"]
    QC["anchor_max"] = opt["anchor_max"]
    G.update(opt)
    G["ref"] = load_reference(opt["reference"])
    G["refs"] = {}
    zp = os.path.join(os.path.dirname(opt["reference"]), "zones.npz")
    if os.path.exists(zp):  # one reference per zone (reference.zone_km)
        z = np.load(zp)
        for k in range(len(z["c_ref"])):
            ref = dict(G["ref"])
            ref["ref_valid"], ref["u_ref"], ref["c_fill"], ref["u_fill"] = reference_work(
                ref["freqs"], z["c_ref"][k].astype(float), np.ones(len(ref["freqs"]), bool)
            )
            ref["ref_support"] = z["ref_support"][k] & ref["ref_valid"]
            G["refs"][k] = ref
    with h5py.File(opt["stack"], "r") as h:
        G["fs"] = float(h.attrs["fs"])
        G["maxlag"] = float(h.attrs["maxlag"])
        G["dist"] = h["dist_km"][:]
        G["nseg"] = h["nseg"][:]


def symmetric_images(x2, dist, fs, maxlag, freqs, alpha):
    """(analytic signal, SNR image [npair, nf, nvel], largest SNR [npair, nf]) of the symmetric component."""
    pos, neg = S.branches(x2, fs, maxlag)
    t = np.arange(pos.shape[1]) / fs
    an = S.analytic_bank(0.5 * (pos + neg), fs, freqs, alpha)
    times = np.broadcast_to(dist[:, None, None] / VEL[None, None, :], (len(dist), len(freqs), len(VEL))).copy()
    img = S.sample_at_times(np.abs(an), fs, times)
    nrm = np.empty((len(dist), len(freqs)))
    for q in range(len(dist)):
        for k, fc in enumerate(freqs):
            m = t > min(dist[q] / VEL[0] + 2.0 / fc, 0.85 * maxlag)
            nrm[q, k] = np.sqrt(np.mean(np.real(an[q, k, m]) ** 2)) + 1e-30
    simg = img / nrm[:, :, None]
    return an, simg, simg.max(axis=2)


def at_pick(img, v):
    """Values of the SNR image [nf, nvel] at the picked velocities v [nf] (NaN where nothing was picked)."""
    out = np.full(len(v), np.nan)
    for k in range(len(v)):
        if np.isfinite(v[k]):
            out[k] = np.interp(v[k], VEL, img[k])
    return out


def lobe_metrics(prof, j):
    """Height of the strongest distinct competing peak relative to the pick (valley below 70 % of the lower of
    the two)."""
    p = prof
    pj = max(p[j], 1e-30)
    side = 0.0
    maxima = np.where((p[1:-1] >= p[:-2]) & (p[1:-1] >= p[2:]) & (p[1:-1] > 0))[0] + 1
    for m in maxima:
        if m == j:
            continue
        a, b = (m, j) if m < j else (j, m)
        if p[a : b + 1].min() <= 0.7 * min(p[m], pj):
            side = max(side, p[m] / pj)
    return side


def pick_group(img, u_ref, valid, segments=None, support=None):
    """Group-velocity ridge.  Where the reference is backed by measurements (support) the search is limited to
    0.65-1.35 U_ref and pulled towards it beyond 25 %; elsewhere (an extrapolated reference) the whole velocity axis
    is searched, as FTAN picks the envelope maximum, and only the continuity between neighbouring frequencies
    (3 %) guides the ridge."""
    support = np.ones(len(u_ref), bool) if support is None else support
    if segments is not None and not (len(segments) == 1 and segments[0] == (0, len(u_ref))):
        vg, side = np.full(len(u_ref), np.nan), np.full(len(u_ref), np.nan)
        for k0, k1 in segments:
            if valid[k0:k1].any():
                vg[k0:k1], side[k0:k1] = pick_group(img[k0:k1], u_ref[k0:k1], valid[k0:k1], support=support[k0:k1])
        return vg, side
    prof = img / (img.max(axis=1, keepdims=True) + 1e-30)
    s = support[:, None]
    inwin = ~s | ((VEL[None, :] >= 0.65 * u_ref[:, None]) & (VEL[None, :] <= 1.35 * u_ref[:, None]))
    prof_dp = np.where(inwin, prof, 1e-3)
    path = S.ridge_dp(prof_dp, VEL, prior=np.where(support, u_ref, np.nan), prior_tol=0.25, smooth=0.03, valid=valid)
    nf = len(path)
    vg, side = np.full(nf, np.nan), np.full(nf, np.nan)
    plaus = ~s | ((VEL[None, :] >= 0.6 * u_ref[:, None]) & (VEL[None, :] <= 1.45 * u_ref[:, None]))
    for k in range(nf):
        if not valid[k]:
            continue
        j = S.local_peak(prof_dp[k], path[k], 3)
        jr = S.refine_peak(prof[k], j)
        vg[k] = np.interp(jr, np.arange(len(VEL)), VEL)
        side[k] = lobe_metrics(np.where(plaus[k], prof[k], 0.0), j)
    return vg, side


def matched_phase(sig, fs, ffine, tg_fine, alpha, phi0, taper_s=1.0):
    """Phase travel times (any cycle) at the frequencies ffine from the one-sided signal sig (t = 0 first), by the
    phase-matched filter built from the group delays tg_fine.  -> t_p [nf]"""
    nt = len(sig)
    w = np.ones(nt)
    ntap = max(1, min(nt // 2, int(taper_s * fs)))
    w[-ntap:] = 0.5 + 0.5 * np.cos(np.pi * np.arange(ntap) / ntap)
    nfft = int(2 ** np.ceil(np.log2(max(4 * nt, 4096))))
    X = np.fft.rfft(sig * w, n=nfft)
    f = np.fft.rfftfreq(nfft, 1 / fs)
    tg = np.interp(f, ffine, tg_fine)
    psi = 2 * np.pi * np.concatenate([[0.0], np.cumsum(0.5 * (tg[1:] + tg[:-1]) * np.diff(f))])
    psi_c = np.interp(ffine, f, psi)
    tp = np.full(len(ffine), np.nan)
    half = np.sqrt(20.0 / alpha)  # the Gaussian is below 2e-9 beyond
    df = f[1]
    for k, fc in enumerate(ffine):
        i0, i1 = max(1, int(fc * (1 - half) / df)), min(len(f), int(fc * (1 + half) / df) + 2)
        g = np.exp(-alpha * ((f[i0:i1] - fc) / fc) ** 2)
        z = np.sum(X[i0:i1] * g * np.exp(1j * (psi[i0:i1] - psi_c[k])))
        if z != 0:
            tp[k] = -(np.angle(z) - phi0) / (2 * np.pi * fc)
    return tp


def excess_spectrum(sym, dist, fs, maxlag, vmin, fmin):
    """Sum over the pairs of the excess energy spectrum of the symmetric correlations sym [npair, nt] (lags 0 ..
    maxlag): the spectrum of the signal window (|t| up to the arrival at vmin plus two periods of fmin, the window
    of the SNR) minus that of the later lags (noise) scaled to the same length.  Incoherent noise, which is the
    same at every lag, cancels in the mean, so the result does not depend on how much data were stacked."""
    nt = sym.shape[1]
    t = np.arange(nt) / fs
    nfft = int(2 ** np.ceil(np.log2(max(4 * nt, 4096))))
    E = np.zeros(nfft // 2 + 1)
    for q in range(len(sym)):
        t1 = min(dist[q] / vmin + 2.0 / fmin, 0.85 * maxlag)
        ws = np.where(
            t <= 0.8 * t1, 1.0, np.where(t < t1, 0.5 + 0.5 * np.cos(np.pi * (t - 0.8 * t1) / (0.2 * t1)), 0.0)
        )
        late = t >= t1
        wn = np.zeros(nt)
        wn[late] = np.hanning(late.sum())
        Fs = np.fft.rfft(sym[q] * ws, n=nfft)
        Fn = np.fft.rfft(sym[q] * wn, n=nfft)
        # the correlation is even in lag: the two-sided signal window has the real spectrum 2 Re(Fs) - x(0)
        E += (2 * Fs.real - sym[q, 0]) ** 2 - 2 * np.abs(Fn) ** 2 * (2 * (ws**2).sum() - 1) / (2 * (wn**2).sum())
    return E


def band_amplitude(stack_path, idx, freqs, alpha, vmin, block=500):
    """Band-edge measure of every frequency, from all the pairs idx of the stack together: the amplitude of the
    coherent signal (square root of the excess energy spectrum, excess_spectrum) at the frequency (mean over a
    band a quarter as wide as the measurement filter) relative to its mean over the measurement filter.  Near 1
    inside the band the noise excites, falling towards 0 beyond its edge."""
    with h5py.File(stack_path, "r") as h:
        fs, maxlag = float(h.attrs["fs"]), float(h.attrs["maxlag"])
        dist = h["dist_km"][:]
        E = 0.0
        for i in range(0, len(idx), block):
            j = idx[i : i + block]
            pos, neg = S.branches(h["stack"][j].astype(np.float64), fs, maxlag)
            E = E + excess_spectrum(0.5 * (pos + neg), dist[j], fs, maxlag, vmin, freqs[0])
    A = np.sqrt(np.clip(E, 0.0, None))
    f = np.fft.rfftfreq(2 * (len(A) - 1), 1 / fs)
    out = np.full(len(freqs), np.nan)
    for k, fc in enumerate(freqs):
        x = ((f - fc) / fc) ** 2
        wide, narrow = np.exp(-alpha * x), np.exp(-16 * alpha * x)
        if narrow.sum() > 0 and (A * wide).sum() > 0:
            out[k] = (A * narrow).sum() / narrow.sum() / ((A * wide).sum() / wide.sum())
    return out


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


def process_chunk(task):
    zone, idx = task
    fs, maxlag, alpha = G["fs"], G["maxlag"], G["alpha"]
    ref = G["refs"].get(zone, G["ref"])
    freqs = ref["freqs"]
    ref_valid = ref["ref_valid"]
    support = ref["ref_support"]
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
    _, img_sym, snr_sym = symmetric_images(x2, dist, fs, maxlag, freqs, alpha)
    band_amp = G["band_amp"]
    res = {k: np.full((npair, nf), np.nan, np.float32) for k in DATA_KEYS}
    res["reason"] = np.zeros((npair, nf), np.int32)
    res["ok_c"] = np.zeros((npair, nf), bool)
    res["pre_c"] = np.zeros((npair, nf), bool)
    res["cycle"] = np.full((npair, nf), -32768, np.int16)
    res["anchor_f"] = np.full((npair, nf), np.nan, np.float32)
    tot_w = dcount.sum(axis=0)
    tot_sum = (daily * dcount[:, :, None]).sum(axis=0)
    po = int(round(1.0 / np.log2(freqs[1] / freqs[0])))  # points per octave of the measurement grid
    ffine = S.log_freqs(freqs[0], freqs[-1], 4 * po)  # phase tracked on a 4x finer grid
    lff = np.log(ffine)
    lf = np.log(freqs)
    c_ref_f = np.interp(lff, lf, c_ref)  # the reference on the fine tracking grid
    kc = np.clip(np.round(np.log2(ffine / freqs[0]) * po).astype(int), 0, nf - 1)
    from .qc import REASONS as RS

    for q in range(npair):
        if G["nseg"][idx[q]] == 0:
            res["reason"][q] = RS["no_data"]
            continue
        if dist[q] < G["min_dist"]:
            res["reason"][q] = RS["path_short_measure"]
            continue
        lam_c = c_ref / freqs
        img = img_sym[q]
        snr = snr_sym[q]
        reason = np.zeros(nf, np.int32)
        reason[~ref_valid] |= RS["no_reference"]
        reason[dist[q] < QC["lam_c"] * lam_c] |= RS["path_short_measure"]
        reason[snr < QC["snr_g"]] |= RS["snr_group"]
        valid = (dist[q] >= QC["lam_c"] * lam_c) & (snr >= QC["snr_g"]) & ref_valid
        if valid.sum() < QC["min_run"]:
            res["reason"][q] = reason | RS["short_run"]
            continue
        vg, side = pick_group(img, u_ref, valid, segs, support)
        snr = at_pick(img, vg)  # the SNR of the picked arrival, not of the strongest one
        # jackknife groups from this pair's own sub-stacks with data, dealt out in turn
        have = np.where(dcount[:, q] > 0)[0]
        groups = [have[g :: G["njack"]] for g in range(G["njack"])]
        subsets = []
        for g in groups:
            w = tot_w[q] - dcount[g, q].sum()
            if g.size == 0 or w <= 0 or dcount[g, q].sum() <= 0:
                continue
            subsets.append((tot_sum[q] - (daily[g, q] * dcount[g, q][:, None]).sum(axis=0)) / w)
        dev_g = np.where(support, np.abs(vg - u_ref) / u_ref, 0.0)  # no deviation screen off the measured reference
        res["vg"][q], res["snr"][q], res["side"][q] = vg, snr, side
        reason[~np.isfinite(vg)] |= RS["no_group_pick"]
        with np.errstate(invalid="ignore"):
            reason[snr < QC["snr_c"]] |= RS["snr_group"]
            reason[dev_g > QC["dev_pre"]] |= RS["dev_pre"]
            reason[(side > QC["side_c"]) & ~side_off] |= RS["side_lobe"]
            pre0 = (
                (snr >= QC["snr_c"])
                & (dist[q] >= QC["lam_c"] * lam_c)
                & (dev_g <= QC["dev_pre"])
                & ((side <= QC["side_c"]) | side_off)
                & ref_valid
                & np.isfinite(vg)
            )
        pre = S.contiguous(pre0, QC["min_run"])
        reason[pre0 & ~pre] |= RS["short_run"]
        res["pre_c"][q] = pre
        if not pre.any():
            res["reason"][q] = reason
            continue
        vg_fill = np.where(np.isfinite(vg), vg, u_ref)
        tg_f = dist[q] / np.interp(lff, lf, vg_fill)
        lam_ratio_f = dist[q] * ffine / c_ref_f
        sigsets = [x2[q : q + 1]] + [sub[None, :] for sub in subsets]
        tp_sets = []
        for sig2 in sigsets:
            p, n = S.branches(sig2, fs, maxlag)
            tp_sets.append(matched_phase(0.5 * (p[0] + n[0]), fs, ffine, tg_f, alpha, phi0))
        c_fine_all = np.full((len(sigsets), len(ffine)), np.nan)
        m_fine = np.full(len(ffine), np.nan)
        anchor_run = np.full(nf, np.nan)
        for k0, k1 in S.runs_of(pre):
            fine = np.where((kc >= k0) & (kc < k1))[0]
            anchor = fine[
                (lam_ratio_f[fine] >= QC["lam_c"]) & (lam_ratio_f[fine] <= QC["anchor_max"]) & support[kc[fine]]
            ]  # the cycle count only where the reference is measured
            if anchor.size == 0:
                reason[k0:k1] |= RS["no_anchor"]
                continue
            st, sp = int(anchor[0]), int(fine[-1]) + 1
            anchor_run[k0:k1] = ffine[st]
            for s_i, tp0 in enumerate(tp_sets):
                c_fine_all[s_i, st:sp] = track_phase(tp0, dist[q], ffine, c_ref_f, st, sp, tg_f=tg_f)[st:sp]
                if s_i == 0:
                    with np.errstate(invalid="ignore", divide="ignore"):
                        m_fine[st:sp] = np.round((dist[q] / c_fine_all[0, st:sp] - tp0[st:sp]) * ffine[st:sp])
        res["anchor_f"][q] = anchor_run
        for k in range(nf):
            fk = np.where(kc == k)[0]
            if fk.size and np.all(np.isfinite(m_fine[fk])):
                res["cycle"][q, k] = int(m_fine[fk[np.argmin(np.abs(lff[fk] - lf[k]))]])
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
            if fk.any() and np.all(np.isfinite(c_fine[fk])):
                vc[k] = np.interp(lf[k], lff[fk], c_fine[fk]) if fk.sum() > 1 else c_fine[fk][0]
                # inf: too few groups with data (no error estimate) -> rejected as "jackknife"; nan: a group lost
                # the phase track -> rejected as well
                jk_c[k] = np.max(np.nan_to_num(jk_c_f[fk], nan=np.inf))
        with np.errstate(divide="ignore", invalid="ignore"):
            dlnc = S.gradient_in_runs(np.log(vc), lf)  # within each run of measured samples, one-sided at its ends
            u_from_c = vc / (1.0 - dlnc)
            ucons = np.abs(u_from_c - vg) / vg <= QC["ucons_c"]
            jk_ok = jk_c <= QC["jk_c_rel"] * vc
            band_ok = ~(band_amp < G["band_edge"])
            ok0 = pre & np.isfinite(vc) & jk_ok & ucons & band_ok
        tracked = pre & np.isfinite(anchor_run)
        reason[tracked & ~np.isfinite(vc)] |= RS["phase_slip"]
        meas = pre & np.isfinite(vc)
        reason[meas & ~band_ok] |= RS["off_band"]
        reason[meas & ~jk_ok] |= RS["jackknife"]
        reason[meas & ~ucons] |= RS["group_phase"]
        ok_c = S.contiguous(ok0, QC["min_run"])
        reason[ok0 & ~ok_c] |= RS["short_run"]
        res["vc"][q], res["jk_c"][q], res["ok_c"][q], res["u_from_c"][q] = vc, jk_c, ok_c, u_from_c
        res["reason"][q] = np.where(ok_c, 0, reason)
    return idx, res


def auto_njack(stack_path):
    """Lines: the number of sub-stacks holding at least 25 % of the segments of the fullest one."""
    with h5py.File(stack_path, "r") as h:
        c = h["daily_count"][:].sum(1)
    return int((c >= 0.25 * c.max()).sum())


def measure(cfg, stack, reference, out):
    """Run the measurement on every pair with data -> out (HDF5)."""
    from .config import phase_constant

    d = cfg["dispersion"]
    njack = auto_njack(stack) if d["njack"] == "auto" else int(d["njack"])
    if njack < d["min_jackknife_groups"]:
        raise RuntimeError(
            "%d jackknife groups (sub-stacks holding data: stack.block_hours, dispersion.njack) but "
            "dispersion.min_jackknife_groups is %d: no error could be estimated; use shorter sub-stacks "
            "(stack.block_hours) or more groups" % (njack, d["min_jackknife_groups"])
        )
    q = dict(d["qc"])
    lam_c_final = q.pop("lam_c")
    q.pop("min_run")
    qc_off = dict(dev_pre=q.pop("dev_pre"), ucons_c=q.pop("ucons_c"))
    opt = dict(
        stack=stack,
        reference=reference,
        alpha=d["alpha"],
        njack=njack,
        min_jackknife_groups=d["min_jackknife_groups"],
        lam_min=d["lam_min"],
        anchor_max=d["anchor_max_lambda"],
        side_off_bands=d["side_off_bands"],
        min_dist=d["min_dist_km"],
        vmin=d["vmin"],
        vmax=d["vmax"],
        dv=d["dv"],
        phi0=phase_constant(cfg),
        qc=dict(
            snr_c=q["snr_c"], snr_g=q["snr_g"], jk_c_rel=q["jk_c_rel"], side_c=q["side_c"], min_run=d["qc"]["min_run"]
        ),
        qc_off=qc_off,
        band_edge=float(q["band_edge"]),
    )
    with h5py.File(stack, "r") as h:
        pairs = np.array([p.decode() for p in h["pairs"][:]])
        nseg = h["nseg"][:]
        dist = h["dist_km"][:]
    all_idx = np.where(nseg > 0)[0]
    init_worker(opt)
    ref = G["ref"]
    freqs = ref["freqs"]
    # one band-edge measure for the whole array (it must not depend on how the pairs are split among the workers)
    opt["band_amp"] = band_amplitude(stack, all_idx, freqs, d["alpha"], float(VEL[0]))
    nf = len(freqs)
    zone = np.full(len(pairs), -1, int)
    zp = os.path.join(os.path.dirname(reference), "zones.npz")
    if os.path.exists(zp):
        zone = np.load(zp)["zone"]
    chunks = []
    for zk in sorted(set(zone[all_idx].tolist())):  # the pairs of a zone share its reference
        zi = all_idx[zone[all_idx] == zk]
        chunks += [(zk, zi[i : i + d["chunk"]]) for i in range(0, len(zi), d["chunk"])]
    npair_all = len(pairs)
    out_d = {k: np.full((npair_all, nf), np.nan, np.float32) for k in DATA_KEYS}
    out_d["reason"] = np.zeros((npair_all, nf), np.int32)
    out_d["ok_c"] = np.zeros((npair_all, nf), bool)
    out_d["pre_c"] = np.zeros((npair_all, nf), bool)
    out_d["cycle"] = np.full((npair_all, nf), -32768, np.int16)
    out_d["anchor_f"] = np.full((npair_all, nf), np.nan, np.float32)
    os.environ["OMP_NUM_THREADS"] = "1"
    with Pool(int(d["nproc"]), initializer=init_worker, initargs=(opt,)) as pool:
        done = 0
        for idx, res in track(pool.imap_unordered(process_chunk, chunks), len(chunks), "pair batches"):
            for k in out_d:
                out_d[k][idx] = res[k]
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
        h.create_dataset("ref_support", data=ref["ref_support"])
        h.create_dataset("side_off", data=off_band_mask(freqs, d["side_off_bands"]))
        h.create_dataset("vel_axis", data=VEL)
        h.create_dataset("band_amp", data=opt["band_amp"])
        h.create_dataset("zone", data=zone)
        h.create_dataset("pairs", data=pairs.astype("S"))
        h.create_dataset("dist_km", data=dist)
        h.attrs.update(
            dict(
                alpha=d["alpha"],
                njack=njack,
                min_jackknife_groups=opt["min_jackknife_groups"],
                reference=os.path.abspath(reference),
                side_off_bands=d["side_off_bands"] or "none",
                phase_constant=opt["phi0"],
                lam_c_final=lam_c_final,
                **{"qc_" + k: (0.0 if np.isinf(v) else v) for k, v in QC.items()},
            )
        )
    os.replace(out + ".tmp", out)
    return dict(pairs=len(all_idx), njack=njack, phase_points=int(out_d["ok_c"].sum()))
