"""Step-specific narratives from recorded summaries; report assembly stays in lant.report."""

import json

from . import auto


def _round(v):
    if isinstance(v, float):
        return float("%.6g" % v)
    if isinstance(v, (list, tuple)):
        return [_round(x) for x in v]
    if isinstance(v, dict):
        return {k: _round(x) for k, x in v.items()}
    return v


def _short(v, n=70):
    s = json.dumps(_round(v), ensure_ascii=False) if not isinstance(v, str) else v
    return s if len(s) <= n else s[: n - 3] + "..."


def did(cfg, step, s):
    """What the step did, in the order it did it, with the counts (from the step summary)."""
    s = s or {}
    f = _DID.get(step)
    return f(cfg, lambda *k: _dig(s, *k)) if f else []


def _did_clean(cfg, g):
    L = []
    t = g("totals") or {}
    L.append(
        "read the headers of %s candidate files (%s not seismic data or unreadable); formats %s"
        % (g("files"), g("unreadable"), _short(g("formats")))
    )
    n = g("names") or {}
    L.append(
        "matched the header names to the station table: %s matched (%s under another spelling), %s ambiguous, "
        "%s table stations without data, %s records without a station"
        % (
            n.get("matched"),
            n.get("renamed"),
            len(n.get("ambiguous") or []),
            len(n.get("table_without_data") or []),
            len(n.get("data_without_table") or []),
        )
    )
    L.append(
        "chose the vertical channel and one sensor / rate per station: %s stations, raw rates %s Hz"
        % (g("stations"), _short(g("sampling_rates")))
    )
    L.append(
        "cut the records into processing days starting at %02d h UTC (labels: %s): %s days, %s station-days"
        % (int(g("day_start_utc_hour") or 0), g("day_label"), g("days"), g("station_days"))
    )
    if t:
        L.append(
            "put %s records on the UTC sample grid: %s moved by < 0.05 sample, %s interpolated (Lanczos); "
            "removed %s duplicates / identical overlaps; left out %.1f s of overlaps with differing samples "
            "and %s pieces too short; %s non-finite samples became gaps; %s records flipped (sensor pointing "
            "down); gaps are not filled (%s gaps)"
            % (
                t.get("traces_in"),
                t.get("shifted"),
                t.get("interpolated"),
                t.get("duplicates_or_identical_overlaps"),
                float(t.get("conflicting_overlap_s") or 0),
                t.get("dropped_short"),
                t.get("nonfinite_samples"),
                t.get("polarity_reversed"),
                t.get("gaps"),
            )
        )
        L.append(
            "wrote the day archive: %s station-days, %.0f station-hours; no filtering"
            % (g("archive_station_days"), float(g("station_hours") or 0))
        )
    L.append("derived the parameters marked 'rule' below from the geometry, the records and data.velocity_range_kms")
    return L


def _did_cc(cfg, g):
    L = []
    c = cfg["correlation"]
    L.append(
        "per station-day: zero-phase low-pass and decimation to %s Hz%s"
        % (
            c["fs"],
            "; 1-s periodic waveform subtracted hour by hour (%s)"
            % (
                g("pulse_removed")
                if isinstance(g("pulse_removed"), str)
                else "%d stations where the clean step detected it" % len(g("pulse_removed"))
            )
            if g("pulse_removed")
            else "; no 1-s periodic waveform subtracted (none detected)"
            if c["remove_1s_comb"] == "auto"
            else "",
        )
    )
    L.append(
        "segments of %s s every %s s; a segment is dropped when it holds a gap or its peak exceeds %s x the "
        "standard deviation of the station-day; linear detrend, 5 %% cosine taper"
        % (c["seg_len_s"], c["seg_step_s"], c["max_over_std"])
    )
    L.append(
        "cross-coherence conj(Fi) Fj / (|Fi| |Fj|) in %.4g-%.4g Hz (quarter-octave cosine edges), averaged per "
        "hour, lags +-%s s" % (c["fmin_hz"], c["fmax_hz"], c["maxlag_s"])
    )
    L.append(
        "%s processing days: %s computed, %s computed again (settings or input files changed), %s reused "
        "unchanged; %s days skipped (fewer than 2 stations); %s station-days failed to read"
        % (g("days"), g("computed"), g("recomputed"), g("reused"), len(g("skipped_days") or []), len(g("failed") or []))
    )
    return L


def _did_stack(cfg, g):
    L = []
    a = g("all_hours") or {}
    L.append(
        "stacked every hour, weighted by its segments: %s pairs (%s with data), %s segments; sub-stacks of %s h "
        "for the jackknife: %s"
        % (a.get("pairs"), a.get("pairs_with_data"), a.get("segments"), cfg["stack"]["block_hours"], a.get("blocks"))
    )
    q = g("station_qc")
    if q:
        lim = q.get("limits")
        L.append(
            "diagnosis: every sub-stack of a pair against the stack of its other sub-stacks (correlation and "
            "time shift), per station the median over its pairs with the nearest stations: %s station "
            "sub-stacks judged, %s not judged (too few pairs)" % (q.get("evaluated"), q.get("not_judged"))
        )
        if lim:
            L.append(
                "candidates = outliers of the two distributions (median and 1.4826 MAD, Chauvenet's "
                "criterion for N = %s): atanh(correlation) more than %.2f sigma (%.3f) below the array median of the "
                "sub-stack, or a time shift beyond %.2f sigma (%.4f s): %s candidates"
                % (
                    lim["n"],
                    lim["k_corr"],
                    lim["k_corr"] * lim["corr_sigma"],
                    lim["k_shift"],
                    lim["shift_max"],
                    q.get("candidates", 0),
                )
            )
        else:
            L.append("NOT CHECKED: too few station sub-stacks to set thresholds from their distribution")
    r = g("removed") or {}
    L.append(
        "removed %s station sub-stacks and %s station-hours (your lists%s); product: %s"
        % (
            r.get("station_blocks", 0),
            r.get("station_hours", 0),
            " and the confirmed candidates" if cfg["stack"]["exclude_candidates"] else "; candidates are kept",
            _short(g("product"), 80) if isinstance(g("product"), str) else "the stack without them",
        )
    )
    return L


def _did_reference(cfg, g):
    L = []
    if g("source") == "fj":
        npf = g("pairs_per_frequency") or [0]
        L.append(
            "frequency-Bessel spectrum of the stack, one frequency at a time with the pairs within %s "
            "wavelengths (%s-%s pairs per frequency), velocities %s-%s km/s in steps of %s"
            % (
                cfg["reference"]["fj_rmax_lambda"],
                min(npf),
                max(npf),
                cfg["reference"]["fj_cmin"],
                cfg["reference"]["fj_cmax"],
                cfg["reference"]["fj_dc"],
            )
        )
        L.append(
            "picked the ridge (%s): %s picks used, %.3g-%.3g Hz; fitted a cubic in log f, held constant outside"
            % (g("mode"), g("picks_used"), *(g("picked_hz") or [float("nan")] * 2))
        )
    else:
        L.append(
            "read the borrowed curve %s: defined at %s Hz; %s measurement frequencies have no reference"
            % (g("file"), _short(g("defined_hz")), g("frequencies_without_reference"))
        )
    L.append(
        "set the measurement band %s Hz and the velocity axis %s km/s from the curve"
        % (_short(g("band_hz")), _short(g("velocity_axis_kms")))
    )
    return L


def _did_disp(cfg, g):
    L = []
    d = cfg["dispersion"]
    m = g("measure") or {}
    L.append(
        "symmetric component of every pair, Gaussian filters (alpha %s) at %s frequencies per octave, "
        "%.4g-%.4g Hz: %s pairs" % (d["alpha"], d["per_octave"], d["fmin_hz"], d["fmax_hz"], m.get("pairs"))
    )
    L.append(
        "group arrival by ridge tracking guided by the reference; phase-matched filtering with that group delay (phase "
        "constant %s); cycle count fixed where the path is %s-%s reference wavelengths, then followed in "
        "frequency" % (d["phase_constant"], d["lam_min"], d["anchor_max_lambda"])
    )
    L.append("delete-a-group jackknife over %s groups of sub-stacks" % m.get("njack"))
    rc = g("rejections") or {}
    L.append(
        "quality control: %s samples accepted by the measurement; rejected samples by reason (a sample can "
        "have several): %s" % (m.get("phase_points"), ", ".join("%s %d" % kv for kv in rc.items() if kv[1]) or "-")
    )
    q = g("lambda_qc") or {}
    L.append(
        "wavelength rule (path >= %s of the sample's own wavelengths) and runs of >= %s samples: %s -> %s"
        % (d["qc"]["lam_c"], d["qc"]["min_run"], q.get("phase_points_before"), q.get("phase_points"))
    )
    e = g("export") or {}
    L.append(
        "exported %s curves with %s points (periods %.4g-%.4g s; one segment per pair: %s points of other "
        "segments left out)"
        % (
            e.get("phase_curves"),
            e.get("phase_points"),
            cfg["export"]["tmin_s"],
            cfg["export"]["tmax_s"],
            e.get("points_dropped_by_one_segment"),
        )
    )
    return L


def _did_inversion(cfg, g):
    L = []
    p, st, gr = g("project") or {}, g("start") or {}, g("grid") or {}
    L.append(
        "average curve of the array (periods with >= %s pairs), 1-D inversion for the starting model"
        % cfg["inversion"]["median_min_pairs"]
    )
    L.append(
        "main project: %s curves, %s travel times, %s periods (%.4g-%.4g s), %s depth nodes down to %s km"
        % (
            p.get("curves"),
            p.get("travel_times"),
            len(p.get("periods") or []),
            min(p.get("periods") or [0]),
            max(p.get("periods") or [0]),
            len(p.get("depths") or []),
            (p.get("depths") or [None])[-1],
        )
    )
    L.append(
        "grid %s x %s nodes, spacing %s x %s deg; weight %s, damp %s; Vs bounds %s-%s km/s"
        % (
            gr.get("nx"),
            gr.get("ny"),
            gr.get("dvxd"),
            gr.get("dvzd"),
            p.get("weight"),
            p.get("damp"),
            st.get("vmin"),
            st.get("vmax"),
        )
    )
    L.append("check projects built: %s" % (", ".join(n for n, _ in (g("checks", "built") or [])) or "none"))
    deferred = g("checks", "deferred") or []
    if deferred:
        L.append(
            "checkerboards built by the run step (their noise follows from the main fit): %s" % ", ".join(deferred)
        )
    dz = g("checks", "design") or {}
    if dz:
        L.append(
            "checkerboard cells (auto): %s; measured wavelengths %s-%s km over %s stations with data"
            % (dz.get("rule"), *(dz.get("wavelength_range_km") or [None, None]), dz.get("stations_with_data"))
        )
        if dz.get("reversal_rule"):
            L.append("depth-varying checkerboard: one sign reversal at %s" % dz["reversal_rule"])
        for k in ("note", "vertical_skipped"):
            if dz.get(k):
                L.append(dz[k])
    sp = g("short_paths") or {}
    if sp:
        L.append("shortest pair %s km%s" % (sp.get("shortest_pair_km"), "; " + sp["note"] if sp.get("note") else ""))
    return L


def _did_run(cfg, g):
    L = []
    res = g("results") or {}
    L.append(
        "precheck of every project (DSurfTomo's first-iteration rejection rule on the starting model), then "
        "the Python DSurfTomo: %s projects finished, %s invalid (not run), %s failed"
        % (
            sum(1 for v in res.values() if v == "ok"),
            len(g("invalid") or {}),
            len(g("failed") or []),
        )
    )
    nz = g("checkerboard_noise") or {}
    if nz:
        L.append(
            "Gaussian noise of the checkerboard travel times: %.2f %% of each time (%s)%s"
            % (
                100 * nz["value"],
                nz["reason"],
                "; the checkerboards ran after the main project" if auto.noise_is_auto(cfg) else "",
            )
        )
    return L


def _did_evaluate(cfg, g):
    import numpy as np

    L = []
    f = g("residuals", "final") or {}
    L.append(
        "travel-time residuals of the final model: %s times, rms %.2f ms (%.2f %%), %.0f %% kept by the "
        "rejection rule"
        % (
            f.get("n"),
            f.get("rms_ms", float("nan")),
            f.get("rel_rms_percent", float("nan")),
            100 * f.get("kept_fraction", float("nan")),
        )
    )
    s0 = next((r for r in (g("residuals", "iterations") or []) if r.get("model") == "start"), {})
    if s0.get("rms_ms") is not None:
        L.append(
            "the same for the starting model: rms %.2f ms (%.2f %%): the 3-D inversion lowered the residual "
            "by %.1f %%"
            % (
                s0["rms_ms"],
                s0.get("rel_rms_percent", float("nan")),
                100 * (1 - f.get("rms_ms", float("nan")) / s0["rms_ms"]),
            )
        )
    for name, rows in sorted((g("checkerboards") or {}).items()):
        ok = [r for r in rows if r.get("slope") is not None and np.isfinite(r["slope"])]
        if ok:
            best = max(ok, key=lambda r: r["slope"])
            L.append(
                "checkerboard %s: recovered amplitude / input amplitude (regression slope) %.2f at best (%.2f "
                "km depth), %.2f-%.2f over the depths; pattern correlation %.2f-%.2f (on the nodes with DWS > 0: "
                "%d of %d per depth at most).  An anomaly of this size in the model is damped by about that factor"
                % (
                    name,
                    best["slope"],
                    best["depth_km"],
                    min(r["slope"] for r in ok),
                    max(r["slope"] for r in ok),
                    min(r["r"] for r in ok),
                    max(r["r"] for r in ok),
                    max(r["n"] for r in ok),
                    max(r.get("n_all", r["n"]) for r in ok),
                )
            )
    for name, c in sorted((g("controls") or {}).items()):
        if name.startswith("start_"):
            rc, rm = c.get("residual_std_ms") or [], c.get("main_residual_std_ms") or []
            L.append(
                "control %s (another starting model): final models differ by %.3f km/s rms, %.3f km/s at "
                "most%s"
                % (
                    name,
                    c.get("rms_km_s", float("nan")),
                    c.get("max_abs_km_s", float("nan")),
                    ".  Weighted residual (standard deviation) of the control from the first "
                    "to the last iteration: %.1f to %.1f ms; of the main inversion: %.1f to %.1f ms "
                    "(a control that still fits worse has not finished converging: part of the "
                    "difference is then the starting model)" % (rc[0], rc[-1], rm[0], rm[-1])
                    if rc and rm
                    else "",
                )
            )
        else:
            L.append(
                "control %s (smoothing weight %s, damping %s): final models differ by %.3f km/s rms, %.3f km/s at most"
                % (name, c.get("weight"), c.get("damp"), c.get("rms_km_s", float("nan")),
                   c.get("max_abs_km_s", float("nan")))
            )  # fmt: skip
    pr = g("pruned_solver_files") or {}
    if pr.get("files") or pr.get("kept"):
        L.append(
            "removed %d solver files lant does not read from the check projects (iteration models, first/last "
            "residuals, 2-D maps; %.1f MB)%s"
            % (
                pr["files"],
                pr["bytes"] / 1e6,
                "; kept as changed or not written by lant: %s" % ", ".join(pr["kept"][:5]) if pr.get("kept") else "",
            )
        )
    return L


def _did_figures(cfg, g):
    L = []
    L.append("drew %s figures; not drawn: %s" % (len(g("drawn") or {}), _short(g("not_drawn") or "none", 100)))
    return L


_DID = {
    "clean": _did_clean,
    "cc": _did_cc,
    "stack": _did_stack,
    "reference": _did_reference,
    "disp": _did_disp,
    "inversion": _did_inversion,
    "run": _did_run,
    "evaluate": _did_evaluate,
    "figures": _did_figures,
}


def _dig(d, *keys):
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    return d
