"""Inversion check figures (matplotlib): travel-time residuals, 1-D starting model, checkerboard recovery numbers."""

import os

import numpy as np

from . import style
from .style import plt


def residuals(cfg, csv_path, evaluation, out):
    """All-observation residuals: milliseconds histogram, relative RMS by period, and separate mean/std/RMS.
    Show iteration convergence only when independently evaluated iteration predictions already exist."""
    from matplotlib.ticker import FuncFormatter

    style.apply(cfg["plots"]["font_scale"])
    sets = []
    start_csv = csv_path.replace("residuals_final", "residuals_start")
    for path, name, col in ((start_csv, "Start", "0.75"), (csv_path, "Final", style.RED)):
        if os.path.exists(path):
            a = np.genfromtxt(path, delimiter=",", names=True, ndmin=1)
            sets.append(
                (
                    name,
                    col,
                    a,
                    100 * (a["t_obs_s"] - a["t_pred_s"]) / a["t_obs_s"],
                    1e3 * (a["t_obs_s"] - a["t_pred_s"]),
                )
            )
    if not sets:
        return None
    it = evaluation.get("residuals", {}).get("iterations", [])
    many = len(it) > 2
    fig, axs = plt.subplots(1, 3 if many else 2, figsize=(12.5 if many else 9.5, 5.0))
    lim = max(1.0, max(np.max(np.abs(r)) for _, _, _, _, r in sets))
    bins = np.linspace(-lim, lim, 65)  # include every tail observation; no percentile clipping
    statistics = []
    for name, col, a, rel, r in sets:
        if name == "Final":
            axs[0].hist(r, bins=bins, histtype="step", color=col, lw=1.6, label=name)
        else:
            axs[0].hist(r, bins=bins, color=col, label=name)
        kept = int(np.count_nonzero(a["kept_by_rule"]))
        statistics.append(
            [name, "%.1f" % r.mean(), "%.1f" % r.std(), "%.1f" % np.sqrt(np.mean(r**2)), "%d / %d" % (kept, len(r))]
        )
        T = np.unique(a["period_s"])
        axs[1].semilogx(
            T, [np.sqrt(np.mean(rel[a["period_s"] == t] ** 2)) for t in T], ".-", color=col, ms=5, label=name
        )
    axs[0].set_title("(a) All observations", loc="left", fontsize=11)
    axs[1].set_title("(b) Fit by period", loc="left", fontsize=11)
    axs[0].set_xlabel("Observed - predicted travel time (ms)")
    axs[0].axvline(0, color="0.4", ls=":", lw=0.8)
    axs[0].set_ylabel("Travel times")
    axs[0].legend(loc="upper left", frameon=False, fontsize="small")
    axs[1].set_xlabel("Period (s)")
    axs[1].set_ylabel("Relative rms residual (%)")
    Tall = np.concatenate([np.unique(a["period_s"]) for _, _, a, _, _ in sets])
    ticks = [
        m * 10.0**k for k in range(-3, 3) for m in (1, 2, 5) if Tall.min() * 0.95 <= m * 10.0**k <= Tall.max() * 1.05
    ]
    if len(ticks) < 3:
        ticks = [
            m * 10.0**k
            for k in range(-3, 3)
            for m in (1, 1.5, 2, 3, 5, 7)
            if Tall.min() * 0.95 <= m * 10.0**k <= Tall.max() * 1.05
        ]
    axs[1].set_xticks(ticks)
    axs[1].xaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    axs[1].minorticks_off()
    lo0, hi0 = axs[0].get_ylim()
    axs[0].set_ylim(lo0, hi0 * 1.3)  # room for the legend above the histogram
    axs[1].legend(loc="best", frameon=False)
    if many:
        k = np.arange(len(it))
        axs[2].plot(k, [x["rms_ms"] for x in it], "ko-", ms=4)
        axs[2].set_xticks(k)
        axs[2].set_xticklabels(
            ["start" if x["model"] == "start" else x["model"].replace("iter", "").lstrip("0") or "0" for x in it]
        )
        axs[2].set_xlabel("Model (start, iterations)")
        axs[2].set_ylabel("RMS residual (ms)")
        axs[2].set_title("(c) Evaluated iterations", loc="left", fontsize=11)
    fig.tight_layout(rect=(0, 0.22, 1, 1))
    table_ax = fig.add_axes([0.10, 0.045, 0.84, 0.13])
    table_ax.axis("off")
    table = table_ax.table(
        cellText=statistics,
        colLabels=["Model", "Mean (ms)", "Std. dev. (ms)", "RMS (ms)", "Kept / all"],
        cellLoc="center",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9 * cfg["plots"]["font_scale"])
    for cell in table.get_celld().values():
        cell.set_edgecolor("0.85")
        cell.set_linewidth(0.4)
    fig.text(
        0.10,
        0.008,
        "All observations shown, including rejected rows. Residuals measure data fit, not model accuracy.",
        fontsize=9,
    )
    return style.save(fig, out, cfg)


def start_model(cfg, inv_dir, out):
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

    style.apply(cfg["plots"]["font_scale"])
    m = np.loadtxt(os.path.join(inv_dir, "vs1d.txt"))
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.6))
    fit = os.path.join(inv_dir, "vs1d_fit.txt")
    if os.path.exists(fit):
        f = np.loadtxt(fit)
        axs[0].semilogx(f[:, 0], f[:, 1], "ko", ms=3.5, label="Median observations")
        axs[0].semilogx(f[:, 0], f[:, 2], "-", c=style.RED, lw=1.8, label="Prediction of 1-D start")
        rms = 1000 * np.sqrt(np.mean((f[:, 1] - f[:, 2]) ** 2))
        axs[0].text(0.96, 0.05, "RMS difference: %.1f m/s" % rms, transform=axs[0].transAxes, ha="right", fontsize=10)
        axs[0].legend(frameon=False)
    axs[0].set_xlabel("Period (s)")
    axs[0].xaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
    axs[0].xaxis.set_major_formatter(FuncFormatter(lambda value, _: "%g" % value))
    axs[0].xaxis.set_minor_formatter(NullFormatter())
    axs[0].set_ylabel("Phase velocity (km/s)")
    axs[0].set_title("(a) Initial dispersion fit", loc="left")
    axs[1].set_title("(b) Starting Vs", loc="left")
    axs[1].plot(m[:, 1], m[:, 0], "k-", lw=1.1)
    axs[1].plot(m[:-1, 1], m[:-1, 0], "ko", ms=3.5, label="Can update in 3-D")
    axs[1].plot(m[-1, 1], m[-1, 0], "s", mfc="w", mec="k", label="Bottom value held fixed")
    axs[1].axhline(m[-1, 0], c="0.6", ls="--", lw=0.7)
    axs[1].set_ylim(m[-1, 0] * 1.12, -0.02 * m[-1, 0])
    axs[1].set_xlabel("Vs (km/s)")
    axs[1].set_ylabel("Depth (km)")
    axs[1].legend(frameon=False, loc="lower left", bbox_to_anchor=(0, 0.14), fontsize=10)
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    fig.text(
        0.06,
        0.02,
        "Initialization, not the final 3-D model. The fixed bottom is a model boundary, not a detected interface.",
        fontsize=9,
    )
    return style.save(fig, out, cfg)


def _checker_label(name):
    import re

    if name == "cb_vertical_thick":
        return "Thick vertical"
    if name == "cb_vertical":
        return "Depth-varying"
    m = re.match(r"^cb_(?:\d+_)?([\d.]+km)$", name)
    if m:
        return "%s cells" % m.group(1).replace("km", " km")
    return name.replace("cb_", "").replace("_", " ")


def recovery(cfg, evaluation, out):
    cb = evaluation.get("checkerboards") or {}
    if not cb:
        return None
    style.apply(cfg["plots"]["font_scale"])
    names = list(cb)
    labels = [_checker_label(n) for n in names]
    for name in names:
        z = np.array([r["depth_km"] for r in cb[name]], float)
        if not np.all(np.isfinite(z)) or np.any(np.diff(z) <= 0):
            raise ValueError("recovery depths of %s must be finite and strictly increasing" % name)
    # a depth can be missing in one test (an input pattern of zero there, e.g. at a sign reversal): blank cell
    z_km = np.array(sorted({round(float(r["depth_km"]), 6) for n in names for r in cb[n]}))
    if not len(z_km):
        raise ValueError("no recovery depths")
    z_m = z_km * 1000.0
    if len(z_m) == 1:
        edges = np.array([z_m[0] - 0.5, z_m[0] + 0.5])
    else:
        mid = (z_m[:-1] + z_m[1:]) / 2
        edges = np.r_[z_m[0] - (mid[0] - z_m[0]), mid, z_m[-1] + (z_m[-1] - mid[-1])]
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 7.2))
    for ax, field, lims, cblabel, title, cmap in (
        (axs[0, 0], "slope", (-20, 100), "Amplitude (%); 100% is ideal", "(a) Amplitude recovery", "Blues"),
        (axs[0, 1], "r", (-1, 1), "Correlation r; 1.0 is ideal", "(b) Pattern recovery", "RdBu"),
    ):
        vals = np.full((len(names), len(z_km)), np.nan)
        for i, n in enumerate(names):
            for r in cb[n]:
                vals[i, int(np.argmin(np.abs(z_km - r["depth_km"])))] = r[field]
        if field == "slope":
            vals = 100.0 * vals
        finite = vals[np.isfinite(vals)]
        if field == "slope" and finite.size:
            lims = (min(lims[0], float(finite.min())), max(lims[1], float(finite.max())))
        im = ax.pcolormesh(
            edges,
            np.arange(len(names) + 1) - 0.5,
            np.ma.masked_invalid(vals),
            shading="flat",
            vmin=lims[0],
            vmax=lims[1],
            cmap=cmap,
        )
        ax.set_yticks(np.arange(len(names)))
        ax.set_yticklabels(labels)
        ax.invert_yaxis()
        ax.set_xlabel("Depth (m)")
        ax.set_title(title, loc="left", fontsize=11)
        fig.colorbar(im, ax=ax, label=cblabel, shrink=0.86)
    markers = ("o", "s", "^", "D", "v", "P")
    for k, (name, lbl) in enumerate(zip(names, labels)):
        rows = cb[name]
        z = [r["depth_km"] for r in rows]
        mk = markers[k % len(markers)]
        axs[1, 0].plot([r["slope"] for r in rows], z, marker=mk, ls="-", ms=3.5, lw=1.3, label=lbl)
        axs[1, 1].plot([r["r"] for r in rows], z, marker=mk, ls="-", ms=3.5, lw=1.3, label=lbl)
    axs[1, 0].set_title("(c) Amplitude by depth", loc="left")
    axs[1, 1].set_title("(d) Correlation by depth", loc="left")
    axs[1, 0].set_xlabel(r"Amplitude factor $\beta$")
    axs[1, 1].set_xlabel("Pattern correlation r")
    axs[1, 0].set_ylabel("Depth (km)")
    axs[1, 1].set_ylabel("Depth (km)")
    for ax in (axs[1, 0], axs[1, 1]):
        ax.invert_yaxis()
        values = np.concatenate([np.asarray(line.get_xdata(), float) for line in ax.lines])
        finite = values[np.isfinite(values)]
        lo, hi = (
            min(-1.0 if ax is axs[1, 1] else 0.0, float(finite.min()) if finite.size else 0.0),
            max(1.0, float(finite.max()) if finite.size else 1.0),
        )
        margin = 0.05 * (hi - lo)
        ax.set_xlim(lo - margin, hi + margin)
        ax.axvline(0, c="0.65", lw=0.7, ls="--")
        ax.axvline(1, c="0.45", lw=0.8, ls=":")
    axs[1, 1].legend(frameon=False, fontsize=8.5, loc="best")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.text(
        0.04,
        0.012,
        "Nodes with nonzero data sensitivity (DWS > 0); cells are samples, not resolved depth intervals.\n"
        "High r with small beta: similar pattern, damped amplitude; r < 0: reversed sign. Ideal recovery: beta = r = 1.",
        fontsize=8.5,
    )
    return style.save(fig, out, cfg)


def kernels(cfg, inv_dir, out, periods=None, dz=None, step=0.05):
    """Illustrative 1-D sensitivity, not the inversion Jacobian: a separate fine model extended to 1.5 x the
    inversion bottom. Joint Vs +/-step perturbations carry Vp/rho through Brocher, divided by node spacing (per km).
    The inversion instead combines separate material partial derivatives on its own parameter nodes."""
    from ..inversion.start import forward

    style.apply(cfg["plots"]["font_scale"])
    m = np.loadtxt(os.path.join(inv_dir, "vs1d.txt"))
    avg = np.loadtxt(os.path.join(inv_dir, "average_curve.txt"), ndmin=2)
    T = np.asarray(periods or sorted({float(avg[int(round(i)), 0]) for i in np.linspace(0, len(avg) - 1, 6)}))
    zmax = 1.5 * m[-1, 0]
    dz = dz or zmax / 45  # DSurfTomo's forward allows about 50 nodes (lant.inversion.start.MAX_LAYERS)
    z = np.round(np.arange(0, zmax + 1e-9, dz), 4)
    vs = np.interp(z, m[:, 0], m[:, 1])
    K = np.zeros((len(z), len(T)))
    for j in range(len(z) - 1):
        vp, vm = vs.copy(), vs.copy()
        vp[j] *= 1 + step
        vm[j] *= 1 - step
        cp, cm = forward(z, vp, T, 0, 3), forward(z, vm, T, 0, 3)
        K[j] = (np.log(cp) - np.log(cm)) / (2 * step) / dz  # per km of node spacing (the surface node covers half)
    fig, ax = plt.subplots(figsize=(5.2, 6.0))
    for k, t in enumerate(T):
        ax.plot(K[:-1, k], z[:-1], "-", lw=1.6, label="%.2f s" % t)
    ax.axhline(m[-1, 0], color="0.5", ls="--", lw=0.8)
    ax.text(
        0.02,
        m[-1, 0],
        "Inversion bottom boundary",
        ha="left",
        va="bottom",
        fontsize=9,
        color="0.4",
        transform=ax.get_yaxis_transform(),
    )
    ax.invert_yaxis()
    ax.set_title("Illustrative 1-D depth sensitivity", loc="left", fontsize=11)
    ax.set_xlabel(r"$\partial\ln c / \partial\ln V_s$ (per km)")
    ax.set_ylabel("Depth (km)")
    ax.legend(frameon=False, title="Period", loc="center right", fontsize=9)
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    fig.text(
        0.05,
        0.02,
        "Separate fine model extended below the inversion.\nJoint Vs +/-%.0f%%; Vp/rho follow. Not the inversion Jacobian."
        % (100 * step),
        fontsize=8,
        wrap=True,
    )
    np.savetxt(
        os.path.join(inv_dir, "kernels.txt"),
        np.c_[z, K],
        fmt="%8.3f",
        header="depth_km  illustrative joint d ln c / d ln Vs per km; dz_km=%.9g; inversion_bottom_km=%.9g; fractional_step=%.9g; periods "
        % (dz, m[-1, 0], step)
        + " ".join("%.4f" % t for t in T),
    )
    return style.save(fig, out, cfg)
