"""Two-station Rayleigh phase velocity of every pair (step 2, with lant.reference).

For each pair of the stack, at the measurement frequencies f (disp.per_octave per octave):
  1. symmetric component: the average of the positive and the time-reversed negative lags (Bensen et al. 2007);
  2. narrow-band analytic signals, Gaussian filters exp(-alpha ((f' - f) / f)^2) (FTAN; Levshin & Ritzwoller 2001);
  3. group arrival t_g: the envelope maximum between dist / ((1 + w) U_ref) and dist / ((1 - w) U_ref), w =
     disp.group_window, U_ref the reference group velocity; a maximum on the edge of this window is no arrival;
     signal-to-noise ratio: the envelope at t_g over the rms of the filtered signal after the window (two periods
     after its end to the end of the lag window, at least the last 15 %);
  4. phase by the phase-matched filter (Herrin & Goforth 1977; FTAN): the spectrum is multiplied by exp(i Psi),
     Psi(w) = integral of t_g dw, which takes the dispersion out of the arrival, and the Gaussian-weighted sum
     around f gives its phase Phi(f); the phase travel time t_p = -(Phi - pi/4) / (2 pi f) + N / f, where pi/4 is
     the far-field phase of the 2-D Green's function in the correlation (Lin et al. 2008) and N the cycle count;
  5. cycle count: at the lowest frequency of each run of usable frequencies, the N that brings dist / t_p closest to
     the reference; this is unambiguous when the reference is within half a branch spacing, i.e. when the path is
     at most 1 / (2 e) wavelengths long for a reference error e (disp.reference_error); then N is followed to the
     next frequency with the measured group delay, w_k t_k = w_(k-1) t_(k-1) + (t_g,k + t_g,k-1) / 2 (w_k - w_(k-1))
     (AFTAN), on a grid four times finer than the measurement grid; a jump of more than disp.slip_cycles of a cycle
     ends the run;
  6. error: delete-one-group jackknife over the sub-stacks of the cc step (each sub-stack left out in turn; Efron &
     Stein 1981), sigma^2 = (g - 1) / g sum (c_i - mean)^2; it measures the repeatability over time only;
  7. quality control (every threshold in disp): SNR >= snr_min; path >= min_wavelengths wavelengths of the measured
     c (Bensen et al. 2007; Luo et al. 2015); sigma <= max_jackknife_rel c; then of the accepted frequencies the
     longest consecutive run, at least min_points long, is the pair's curve (one segment per pair, as the manual
     of EGFAnalysisTimeFreq advises).
Every rejected value keeps its reasons (REASONS, a bit mask).
"""

import numpy as np

from .signal import analytic_bank, branches, log_freqs, runs_of

PHI0 = np.pi / 4
FINE = 4  # phase tracked on a grid FINE times denser than the measurement grid

REASONS = {
    "no_data": 1,
    "no_reference": 2,
    "path_short": 4,
    "no_arrival": 8,
    "snr": 16,
    "cycle_ambiguous": 32,
    "phase_slip": 64,
    "jackknife": 128,
    "short_run": 256,
    "not_longest_run": 512,
}


def matched_phase(sig, fs, ffine, tg_fine, alpha, phi0=PHI0, taper_s=1.0):
    """Phase travel times (cycle count 0) at ffine of the one-sided signals sig [n, nt] (t = 0 first), with the
    phase-matched filter of the group delays tg_fine -> [n, nf]."""
    sig = np.atleast_2d(sig)
    nt = sig.shape[-1]
    w = np.ones(nt)
    ntap = max(1, min(nt // 2, int(taper_s * fs)))
    w[-ntap:] = 0.5 + 0.5 * np.cos(np.pi * np.arange(ntap) / ntap)
    nfft = int(2 ** np.ceil(np.log2(max(4 * nt, 4096))))
    X = np.fft.rfft(sig * w, n=nfft, axis=-1)
    f = np.fft.rfftfreq(nfft, 1 / fs)
    tg = np.interp(f, ffine, tg_fine)
    psi = 2 * np.pi * np.concatenate([[0.0], np.cumsum(0.5 * (tg[1:] + tg[:-1]) * np.diff(f))])
    psi_c = np.interp(ffine, f, psi)
    tp = np.full((len(sig), len(ffine)), np.nan)
    half = np.sqrt(20.0 / alpha)  # the Gaussian is below 2e-9 beyond
    df = f[1]
    for k, fc in enumerate(ffine):
        i0, i1 = max(1, int(fc * (1 - half) / df)), min(len(f), int(fc * (1 + half) / df) + 2)
        g = np.exp(-alpha * ((f[i0:i1] - fc) / fc) ** 2)
        z = np.sum(X[:, i0:i1] * g * np.exp(1j * (psi[i0:i1] - psi_c[k])), axis=-1)
        ok = z != 0
        tp[ok, k] = -(np.angle(z[ok]) - phi0) / (2 * np.pi * fc)
    return tp


def group_arrival(env, fs, dist, u_ref, window):
    """Envelope maximum inside the group window -> (t_g, index) or (nan, -1) when the maximum is on its edge."""
    t0, t1 = dist / ((1 + window) * u_ref), dist / ((1 - window) * u_ref)
    i0, i1 = int(np.floor(t0 * fs)), min(env.size - 1, int(np.ceil(t1 * fs)))
    if i1 - i0 < 2:
        return np.nan, -1
    j = i0 + int(np.argmax(env[i0 : i1 + 1]))
    if j <= i0 or j >= i1:
        return np.nan, -1
    y0, y1, y2 = env[j - 1], env[j], env[j + 1]
    den = y0 - 2 * y1 + y2
    frac = 0.5 * (y0 - y2) / den if den < 0 else 0.0
    return (j + frac) / fs, j


def track(tp0, ffine, tg, dist, c_ref, slip):
    """Cycle count at ffine[0] from the reference, then followed upwards -> phase travel times (NaN after a slip)."""
    T = 1.0 / ffine
    om = 2 * np.pi * ffine
    tp = np.full(len(ffine), np.nan)
    tp[0] = tp0[0] + np.round((dist / c_ref - tp0[0]) / T[0]) * T[0]
    for k in range(1, len(ffine)):
        pred = (om[k - 1] * tp[k - 1] + 0.5 * (tg[k] + tg[k - 1]) * (om[k] - om[k - 1])) / om[k]
        cand = tp0[k] + np.round((pred - tp0[k]) / T[k]) * T[k]
        if abs(cand - pred) > slip * T[k]:
            break
        tp[k] = cand
    return tp


def measure_pair(sym, subsets, dist, fs, freqs, c_ref, u_ref, p):
    """One pair.  sym: one-sided symmetric component of the stack; subsets: those of the delete-one-group stacks
    [g, nt]; p: the disp settings.  -> dict of arrays over freqs: c, u, snr, sigma, reason."""
    nf = len(freqs)
    out = dict(c=np.full(nf, np.nan), u=np.full(nf, np.nan), snr=np.full(nf, np.nan), sigma=np.full(nf, np.nan),
               reason=np.zeros(nf, np.int32))  # fmt: skip
    R = REASONS
    have_ref = np.isfinite(c_ref) & np.isfinite(u_ref)
    out["reason"][~have_ref] |= R["no_reference"]
    short = dist < p["min_wavelengths"] * c_ref / freqs
    out["reason"][have_ref & short] |= R["path_short"]
    an = analytic_bank(sym[None, :], fs, freqs, p["alpha"])[0]
    env = np.abs(an)
    nt = sym.size
    t = np.arange(nt) / fs
    tg = np.full(nf, np.nan)
    for k in np.where(have_ref & ~short)[0]:
        tg[k], j = group_arrival(env[k], fs, dist, u_ref[k], p["group_window"])
        if j < 0:
            out["reason"][k] |= R["no_arrival"]
            continue
        t_end = dist / ((1 - p["group_window"]) * u_ref[k]) + 2.0 / freqs[k]
        m = t > min(t_end, 0.85 * t[-1])
        noise = np.sqrt(np.mean(np.real(an[k, m]) ** 2)) + 1e-30
        out["snr"][k] = env[k, j] / noise
        out["u"][k] = dist / tg[k]
        if out["snr"][k] < p["snr_min"]:
            out["reason"][k] |= R["snr"]
    usable = out["reason"] == 0
    po = int(round(1.0 / np.log2(freqs[1] / freqs[0])))
    for k0, k1 in runs_of(usable):
        if k1 - k0 < p["min_points"]:
            out["reason"][k0:k1] |= R["short_run"]
            continue
        if dist * freqs[k0] / c_ref[k0] > 1.0 / (2.0 * p["reference_error"]):
            out["reason"][k0:k1] |= R["cycle_ambiguous"]
            continue
        ffine = log_freqs(freqs[k0], freqs[k1 - 1], FINE * po)
        ffine = ffine[ffine <= freqs[k1 - 1] * (1 + 1e-9)]
        tgf = np.interp(np.log(ffine), np.log(freqs[k0:k1]), tg[k0:k1])
        tp0 = matched_phase(np.vstack([sym[None, :], subsets]) if len(subsets) else sym[None, :], fs, ffine, tgf,
                            p["alpha"])  # fmt: skip
        tp = track(tp0[0], ffine, tgf, dist, c_ref[k0], p["slip_cycles"])
        at = np.arange(k1 - k0) * FINE  # the measurement frequencies on the fine grid
        c = dist / tp[at]
        out["reason"][k0:k1][~np.isfinite(c)] |= R["phase_slip"]
        out["c"][k0:k1] = c
        if len(subsets) >= 3:
            T = 1.0 / ffine[at]
            tj = tp0[1:, at]
            tj = tj + np.round((tp[at][None, :] - tj) / T) * T  # each delete-one-group value on the stack's cycle
            cj = dist / tj
            g = cj.shape[0]
            out["sigma"][k0:k1] = np.sqrt((g - 1) / g * np.sum((cj - cj.mean(0)) ** 2, axis=0))
    ok = np.isfinite(out["c"]) & (out["reason"] == 0)
    with np.errstate(invalid="ignore"):
        lam_short = ok & (dist < p["min_wavelengths"] * out["c"] / freqs)
        out["reason"][lam_short] |= R["path_short"]
        if p["max_jackknife_rel"] is not None:
            bad = ok & ~(out["sigma"] <= p["max_jackknife_rel"] * out["c"])
            out["reason"][bad] |= R["jackknife"]
    ok = np.isfinite(out["c"]) & (out["reason"] == 0)
    runs = [(k0, k1) for k0, k1 in runs_of(ok) if k1 - k0 >= p["min_points"]]
    for k0, k1 in runs_of(ok):
        if k1 - k0 < p["min_points"]:
            out["reason"][k0:k1] |= R["short_run"]
    if len(runs) > 1:
        best = max(runs, key=lambda r: (r[1] - r[0], np.nanmean(out["snr"][r[0] : r[1]])))
        for k0, k1 in runs:
            if (k0, k1) != best:
                out["reason"][k0:k1] |= R["not_longest_run"]
    return out


def accepted(reason, c):
    return (reason == 0) & np.isfinite(c)


def reason_counts(reason):
    return {k: int(((reason & v) > 0).sum()) for k, v in REASONS.items()}


# ------------------------------------------------------------------------------------------------- all pairs
G = {}


def _init(shared):
    from . import resources

    G.update(shared)
    resources.limit_threads(1)


def _chunk(idx):
    """Measure the pairs idx (increasing indices of the stack)."""
    import h5py

    fs, maxlag = G["fs"], G["maxlag"]
    with h5py.File(G["stack"], "r") as h:
        x = h["ccf"][idx].astype(np.float64)
        dist = h["dist_km"][idx]
        cnt = h["count"][idx]
    subs, wts = _substacks(idx)
    total = x * cnt[:, None]
    res = []
    for q in range(len(idx)):
        if cnt[q] == 0:
            nf = len(G["freqs"])
            res.append(dict(c=np.full(nf, np.nan), u=np.full(nf, np.nan), snr=np.full(nf, np.nan),
                            sigma=np.full(nf, np.nan), reason=np.full(nf, REASONS["no_data"], np.int32)))  # fmt: skip
            continue
        pos, neg = branches(x[q : q + 1], fs, maxlag)
        sym = 0.5 * (pos + neg)[0]
        sets = []
        for s in range(subs.shape[0]):  # delete-one-group stacks
            w = cnt[q] - wts[s, q]
            if wts[s, q] > 0 and w > 0:
                y = (total[q] - wts[s, q] * subs[s, q]) / w
                p2, n2 = branches(y[None, :], fs, maxlag)
                sets.append(0.5 * (p2 + n2)[0])
        res.append(measure_pair(sym, np.array(sets).reshape(len(sets), -1), dist[q], fs, G["freqs"], G["c_ref"],
                                G["u_ref"], G["p"]))  # fmt: skip
    return idx, res


def _substacks(idx):
    """Sub-stacks of the pairs idx from the day files -> (ccf [s, n, nlag], count [s, n])."""
    import h5py

    out, wts = [], []
    for f in G["day_files"]:
        with h5py.File(f, "r") as h:
            gp = h["pair"][:]
            pos = np.searchsorted(gp, idx)
            pos = np.clip(pos, 0, max(0, gp.size - 1))
            hit = (gp.size > 0) & (gp[pos] == idx) if gp.size else np.zeros(len(idx), bool)
            for s in range(h["count"].shape[0]):
                a = np.zeros((len(idx), G["nlag"]))
                c = np.zeros(len(idx))
                if hit.any():
                    rows = np.unique(pos[hit])
                    blk = h["ccf"][s, rows, :].astype(np.float64)
                    cc = h["count"][s, rows].astype(np.float64)
                    where = {r: k for k, r in enumerate(rows)}
                    for q in np.where(hit)[0]:
                        a[q], c[q] = blk[where[pos[q]]], cc[where[pos[q]]]
                out.append(a)
                wts.append(c)
    return np.array(out).reshape(len(out), len(idx), G["nlag"]), np.array(wts).reshape(len(wts), len(idx))


def measure_all(stack_path, day_files, freqs, c_ref, u_ref, p, n_cores, log=print):
    """Every pair of the stack -> dict of arrays [npair, nf]."""
    from multiprocessing import Pool

    import h5py

    with h5py.File(stack_path, "r") as h:
        npair, nlag = h["ccf"].shape
        fs, maxlag = float(h.attrs["fs"]), float(h.attrs["maxlag"])
    shared = dict(stack=stack_path, day_files=day_files, freqs=freqs, c_ref=c_ref, u_ref=u_ref, p=p, fs=fs,
                  maxlag=maxlag, nlag=nlag)  # fmt: skip
    nsub = 0
    for f in day_files:
        with h5py.File(f, "r") as h:
            nsub += h["count"].shape[0]
    # pairs per task: enough tasks for the cores, and the sub-stacks of a task within about 200 MB
    size = max(1, min(50, npair // (4 * n_cores) or 1, int(2e8 // (8 * max(1, nsub) * nlag))))
    chunks = [np.arange(i, min(npair, i + size)) for i in range(0, npair, size)]
    nf = len(freqs)
    out = {k: np.full((npair, nf), np.nan) for k in ("c", "u", "snr", "sigma")}
    out["reason"] = np.zeros((npair, nf), np.int32)
    done = 0
    with Pool(max(1, min(n_cores, len(chunks))), initializer=_init, initargs=(shared,)) as pool:
        for idx, res in pool.imap_unordered(_chunk, chunks):
            for q, r in zip(idx, res):
                for k in out:
                    out[k][q] = r[k]
            done += len(idx)
            if done == npair or done % max(size, npair // 10) < len(idx):
                log("  pairs measured: %d / %d" % (done, npair))
    return out
