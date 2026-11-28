"""Shared matplotlib style of the figures that are not maps (journal style: serif STIX type, ticks inside)."""

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

RED, BLUE, GREY = "#d62728", "#1f77b4", "0.55"
FAST, SLOW = "#e8201a", "#00a0d8"  # lines of the fastest and slowest velocity on correlations
ENERGY = LinearSegmentedColormap.from_list(
    "energy", ["#ffffff", "#cfe6fa", "#8cc8f0", "#ffe94a", "#ff9a1f", "#e8201a", "#8c0a0a"]
)
# GMT "seis" (red = slow, blue = fast), continuous
_SEIS = [
    (193.91, 0, 0),
    (241.72, 0, 0),
    (255, 34.531, 0),
    (255, 82.344, 0),
    (255, 130.16, 0),
    (255, 177.97, 0),
    (255, 225.78, 0),
    (255, 255, 0),
    (255, 255, 0),
    (198.28, 255, 10.312),
    (105.47, 255, 27.187),
    (47.813, 247.97, 67.5),
    (0, 235, 114.53),
    (0, 145, 196.09),
    (0, 67.5, 247.19),
    (0, 22.5, 219.06),
]
SEIS = LinearSegmentedColormap.from_list("seis", [tuple(v / 255 for v in c) for c in _SEIS])
PERTURBATION = LinearSegmentedColormap.from_list("vs_anomaly", ["#c50000", "#ffffff", "#0032cc"], N=257)


def apply(scale=1.0):
    """Journal style: serif (STIX) type, ticks inside on all four sides, no frame around legends."""
    plt.rcParams.update(
        {
            "font.family": "STIXGeneral",
            "mathtext.fontset": "stix",
            "font.size": 13 * scale,
            "axes.labelsize": 15 * scale,
            "xtick.labelsize": 13 * scale,
            "ytick.labelsize": 13 * scale,
            "legend.fontsize": 12 * scale,
            "legend.frameon": False,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 5,
            "ytick.major.size": 5,
            "xtick.minor.size": 2.5,
            "ytick.minor.size": 2.5,
            "xtick.major.pad": 6,
            "axes.linewidth": 0.9,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.05,
        }
    )


def period_axis(ax):
    """Period (s) on the top axis of a logarithmic frequency axis."""
    sec = ax.secondary_xaxis(
        "top", functions=(lambda x: 1 / np.clip(x, 1e-6, None), lambda x: 1 / np.clip(x, 1e-6, None))
    )
    lo, hi = sorted(1 / np.array(ax.get_xlim()))
    ticks = [v for v in nice_log_ticks(lo, hi)]
    sec.set_xticks(ticks)
    sec.set_xticklabels(["%g" % v for v in ticks])
    sec.minorticks_off()
    sec.tick_params(direction="in")
    sec.set_xlabel("Period (s)", labelpad=8)
    ax.tick_params(axis="x", which="both", top=False)
    return sec


def nice_log_ticks(lo, hi):
    """1, 2, 5 x 10^n within [lo, hi]; 1-9 x 10^n when that leaves fewer than three."""
    for subs in ((1, 2, 5), range(1, 10)):
        t = [m * 10.0**e for e in range(int(np.floor(np.log10(lo))) - 1, int(np.ceil(np.log10(hi))) + 1) for m in subs]
        t = [v for v in t if lo * 0.999 <= v <= hi * 1.001]
        if len(t) >= 3:
            return t
    return t


def freq_axis(ax, f):
    """Logarithmic frequency axis over the frequencies f, labelled 0.5, 1, 2, 5 ...; period on the top axis."""
    log_axis(ax.xaxis)
    ax.set_xlim(f[0] / 1.02, f[-1] * 1.02)
    ax.set_xlabel("Frequency (Hz)")
    period_axis(ax)


def log_axis(axis):
    """Logarithmic axis labelled with plain numbers (0.5, 1, 2, 5 ...) instead of powers of ten."""
    from matplotlib.ticker import FuncFormatter, NullFormatter

    axis.axes.set_xscale("log") if axis is axis.axes.xaxis else axis.axes.set_yscale("log")
    axis.set_major_locator(_PlainLogLocator())
    axis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    axis.set_minor_formatter(NullFormatter())


def plain_log(axis):
    """Plain numbers on an axis that is already logarithmic."""
    from matplotlib.ticker import FuncFormatter, NullFormatter

    axis.set_major_locator(_PlainLogLocator())
    axis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    axis.set_minor_formatter(NullFormatter())


def _PlainLogLocator():
    """Ticks at 1, 2, 5 x 10^n; over less than a decade (where that leaves one or two ticks) at 1-9 x 10^n."""
    import numpy as np
    from matplotlib.ticker import LogLocator

    class Locator(LogLocator):
        def tick_values(self, vmin, vmax):
            lo, hi = sorted((vmin, vmax))
            self._subs = np.array([1.0, 2.0, 5.0]) if hi / max(lo, 1e-30) > 4 else np.arange(1.0, 10.0)
            return super().tick_values(vmin, vmax)

    return Locator(subs=(1.0, 2.0, 5.0))


def save(fig, out, dpi=300):
    """out without extension -> out.png."""
    import os

    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out + ".png", dpi=dpi)
    plt.close(fig)
    return out + ".png"
