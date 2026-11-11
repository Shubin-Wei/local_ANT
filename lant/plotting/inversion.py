"""Inversion check figures (matplotlib): travel-time residuals, 1-D starting model, checkerboard recovery numbers."""

import os

import numpy as np

from . import style
from .style import plt


def residuals(cfg, csv_path, evaluation, out):
    """Relative travel-time residuals of the starting (grey) and final (black) model: histogram and rms against
    period; the rms of every iteration's model only when checks.iteration_residuals has computed them."""
    from matplotlib.ticker import FuncFormatter

    style.apply(cfg["plots"]["font_scale"])
    sets = []
    start_csv = csv_path.replace("residuals_final", "residuals_start")
    for path, name, col in ((start_csv, "Starting model", "0.65"), (csv_path, "Final model", "k")):
        if os.path.exists(path):
            a = np.genfromtxt(path, delimiter=",", names=True)
            sets.append(
                (
                    name,
                    col,
                    a,
                    100 * (a["t_obs_s"] - a["t_pred_s"]) / a["t_obs_s"],
                    1e3 * (a["t_obs_s"] - a["t_pred_s"]),
                )
            )
    it = evaluation["residuals"]["iterations"]
    many = len(it) > 2
    fig, axs = plt.subplots(1, 3 if many else 2, figsize=(14 if many else 10.5, 4.2))
    lim = max(np.percentile(np.abs(rel), 99.5) for _, _, _, rel, _ in sets)
    bins = np.linspace(-lim, lim, 81)
    for name, col, a, rel, r in sets:
        txt = "%s: rms %.1f ms (%.1f %%)" % (name, np.sqrt(np.mean(r**2)), np.sqrt(np.mean(rel**2)))
        if col == "k":
            axs[0].hist(rel, bins=bins, histtype="step", color="k", lw=1.4, label=txt)
        else:
            axs[0].hist(rel, bins=bins, color=col, label=txt)
        T = np.unique(a["period_s"])
        axs[1].semilogx(
            T, [np.sqrt(np.mean(rel[a["period_s"] == t] ** 2)) for t in T], ".-", color=col, ms=5, label=name
        )
    axs[0].set_xlabel("Relative residual, observed - predicted (%)")
    axs[0].set_ylabel("Travel times")
    axs[0].legend(loc="upper left", frameon=False, fontsize="small")
    axs[0].text(0.97, 0.95, "%d travel times" % len(sets[-1][2]), transform=axs[0].transAxes, ha="right", va="top")
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
        axs[2].set_ylabel("rms residual (ms)")
    fig.tight_layout()
    return style.save(fig, out, cfg)


def start_model(cfg, inv_dir, out):
    style.apply(cfg["plots"]["font_scale"])
    m = np.loadtxt(os.path.join(inv_dir, "vs1d.txt"))
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.6))
    fit = os.path.join(inv_dir, "vs1d_fit.txt")
    if os.path.exists(fit):
        f = np.loadtxt(fit)
        axs[0].semilogx(f[:, 0], f[:, 1], "ko", ms=3.5, label="Average curve (median of the pairs)")
        axs[0].semilogx(f[:, 0], f[:, 2], "-", c=style.RED, lw=1.8, label="1-D model prediction")
        axs[0].legend(frameon=False)
    axs[0].set_xlabel("Period (s)")
    axs[0].set_ylabel("Phase velocity (km/s)")
    z = np.r_[m[:, 0], m[-1, 0] + 0.3]
    v = np.r_[m[:, 1], m[-1, 1]]
    axs[1].plot(v, z, "k-", drawstyle="steps-pre")
    axs[1].plot(m[:-1, 1], m[:-1, 0], "ko", ms=3.5, label="Updated nodes")
    axs[1].plot(m[-1, 1], m[-1, 0], "s", mfc="w", mec="k", label="Fixed bottom node")
    axs[1].invert_yaxis()
    axs[1].set_xlabel("Vs (km/s)")
    axs[1].set_ylabel("Depth (km)")
    axs[1].legend(frameon=False, loc="lower left")
    fig.tight_layout()
    return style.save(fig, out, cfg)


def recovery(cfg, evaluation, out):
    cb = evaluation.get("checkerboards") or {}
    if not cb:
        return None
    style.apply(cfg["plots"]["font_scale"])
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.8), sharey=True)
    for name, rows in cb.items():
        z = [r["depth_km"] for r in rows]
        axs[0].plot([r["slope"] for r in rows], z, "o-", ms=3, label=name)
        axs[1].plot([r["r"] for r in rows], z, "o-", ms=3, label=name)
    axs[0].set_xlabel("Recovery slope (recovered on true anomaly)")
    axs[1].set_xlabel("Correlation of recovered and true anomaly")
    axs[0].set_ylabel("Depth (km)")
    axs[0].invert_yaxis()
    for ax in axs:
        ax.set_xlim(-0.2, 1.2)
        ax.axvline(0, c="0.7", lw=0.6)
    axs[1].legend(frameon=False, fontsize=9)
    fig.tight_layout()
    return style.save(fig, out, cfg)


def kernels(cfg, inv_dir, out, periods=None, dz=None, step=0.05):
    """Depth sensitivity of the Rayleigh phase velocity to Vs for the starting model, d ln c / d ln Vs per km of depth:
    central differences (+-5 % at one node; 45 nodes down to 1.5 x the deepest node of the model, the linear layering of DSurfTomo between them) with fwd1d,
    i.e. DSurfTomo's own forward chain in which Vp and density follow Vs (Brocher 2005) -- the sensitivity the inversion
    uses.  Smaller perturbations or node spacings run into the single-precision root search of surfdisp96."""
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
        "fixed bottom node",
        ha="left",
        va="bottom",
        fontsize=9,
        color="0.4",
        transform=ax.get_yaxis_transform(),
    )
    ax.invert_yaxis()
    ax.set_xlabel("d ln c / d ln Vs (per km)")
    ax.set_ylabel("Depth (km)")
    ax.legend(frameon=False, title="Period", loc="lower right")
    np.savetxt(
        os.path.join(inv_dir, "kernels.txt"),
        np.c_[z, K],
        fmt="%8.3f",
        header="depth_km  d ln c / d ln Vs per km, periods " + " ".join("%.4f" % t for t in T),
    )
    return style.save(fig, out, cfg)
