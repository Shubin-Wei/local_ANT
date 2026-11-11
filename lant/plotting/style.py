"""Shared matplotlib style of the non-map figures (journal style: serif STIX type, ticks inside, no titles)."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

RED, BLUE, GREY = "#d62728", "#1f77b4", "0.55"
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


def apply(scale=1.0):
    plt.rcParams.update(
        {
            "font.family": "STIXGeneral",
            "mathtext.fontset": "stix",
            "font.size": 12 * scale,
            "axes.labelsize": 13 * scale,
            "xtick.labelsize": 11 * scale,
            "ytick.labelsize": 11 * scale,
            "legend.fontsize": 11 * scale,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "axes.linewidth": 0.8,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.05,
        }
    )


def save(fig, out, cfg):
    import os

    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    files = []
    for f in cfg["plots"]["formats"]:
        fig.savefig(out + "." + f, dpi=cfg["plots"]["dpi"])
        files.append(out + "." + f)
    plt.close(fig)
    return files
