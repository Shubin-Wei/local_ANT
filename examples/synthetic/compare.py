"""Figures comparing lant's results with the true models (README of the synthetic example).

python compare.py dispersion OUT.png A_CONFIG [A_CONFIG ...]   part A: measured phase and group velocities
                                                                of every model against its theory
python compare.py model OUT B_DIR [z1,z2,..]                    part B: true and recovered Vs perturbation in depth
                                                                slices (OUT_slices.png) and sections
                                                                (OUT_sections_true.png, OUT_sections_lant.png)
"""

import os
import sys

import modes as M
import numpy as np

sys.path.insert(0, os.path.join(M.HERE, "..", ".."))


def dispersion(out, configs):
    from lant import config as C
    from lant.inversion.run import read_dispersion
    from lant.plotting import style
    from lant.plotting.style import plt

    style.apply()
    fig, axs = plt.subplots(2, len(configs), figsize=(6.2 * len(configs), 8.6), squeeze=False,
                            gridspec_kw=dict(height_ratios=[3, 1.3], hspace=0.08, wspace=0.22))  # fmt: skip
    for col, path in enumerate(configs):
        cfg = C.load(path)
        name = open(os.path.join(os.path.dirname(cfg["output_dir"]), "info.txt")).read().split()[1].rstrip(";")
        tab = read_dispersion(os.path.join(cfg["output_dir"], "3_disp", "dispersion.csv"))
        modes = M.read_table(os.path.join(M.HERE, "modes_%s.csv" % name))
        f0, c0, u0, _ = modes[0]
        ax, bx = axs[0, col], axs[1, col]
        freqs = np.unique(tab["frequency_hz"])
        for key, theo, colour, lab in (("c_kms", c0, "k", "phase"), ("u_kms", u0, style.BLUE, "group")):
            v = tab[key]
            med, lo, hi, dev, sd = [], [], [], [], []
            for f in freqs:
                x = v[(tab["frequency_hz"] == f) & np.isfinite(v)]
                t = np.interp(f, f0, theo)
                if x.size >= 10:
                    med.append(np.median(x))
                    lo.append(np.percentile(x, 16))
                    hi.append(np.percentile(x, 84))
                    dev.append(100 * (np.median(x) / t - 1))
                    sd.append(100 * 0.5 * (np.percentile(x, 84) - np.percentile(x, 16)) / t)
                else:
                    med.append(np.nan), lo.append(np.nan), hi.append(np.nan), dev.append(np.nan), sd.append(np.nan)
            ax.fill_between(freqs, lo, hi, color=colour, alpha=0.18, lw=0)
            ax.plot(freqs, med, "o", ms=3.5, color=colour, label="lant, %s (median, 16-84 %%)" % lab)
            m = (f0 >= freqs[0] / 1.2) & (f0 <= freqs[-1] * 1.2)
            ax.plot(f0[m], theo[m], "-", color=style.RED if lab == "phase" else "#ff7f0e", lw=1.6,
                    label="theory, %s" % lab)  # fmt: skip
            bx.fill_between(freqs, np.array(dev) - np.array(sd), np.array(dev) + np.array(sd), color=colour,
                            alpha=0.18, lw=0)  # fmt: skip
            bx.plot(freqs, dev, "o-", ms=3, lw=1, color=colour)
        for m_, (fm, cm, um, _) in modes.items():  # the first higher mode, for orientation
            if m_ == 1:
                ax.plot(fm, cm, ":", color="0.5", lw=1.2, label="theory, first higher mode (phase)")
        style.log_axis(ax.xaxis)
        style.log_axis(bx.xaxis)
        ax.set_xlim(freqs[0] / 1.05, freqs[-1] * 1.05)
        bx.set_xlim(ax.get_xlim())
        ax.tick_params(labelbottom=False)
        ax.set_ylabel("Velocity (km/s)")
        ax.text(0.03, 0.04, "model: %s" % name, transform=ax.transAxes, fontsize=14)
        ax.legend(loc="upper right", fontsize=10)
        bx.axhline(0, color="0.3", lw=0.8)
        bx.set_ylabel("lant / theory - 1 (%)")
        bx.set_xlabel("Frequency (Hz)")
        bx.set_ylim(-6, 6)
    style.save(fig, os.path.splitext(out)[0])


def model(out, bdir, depths=(0.05, 0.15, 0.3, 0.45)):
    """Vs perturbation from the mean of each depth, of the true model and of lant's model (inside the station hull,
    where lant has data): the usual way of showing a synthetic recovery test."""
    from matplotlib.path import Path

    from lant import config as C
    from lant import stations as ST
    from lant.plotting import maps
    from lant.plotting.inv_figures import grid_model

    cfg = C.load(os.path.join(bdir, "config.yaml"))
    sta = ST.read_csv(os.path.join(cfg["output_dir"], "1_clean", "stations_used.csv"))
    hull = Path(maps.hull_ring([s["lon"] for s in sta.values()], [s["lat"] for s in sta.values()]))

    def perturbation(rows):
        lon, lat, dep, V = grid_model(rows)
        LO, LA = np.meshgrid(lon, lat, indexing="ij")
        inside = hull.contains_points(np.c_[LO.ravel(), LA.ravel()]).reshape(LO.shape)
        mean = np.array([np.nanmean(V[:, :, k][inside]) for k in range(len(dep))])
        return lon, lat, dep, 100 * (V / mean - 1)

    t = np.loadtxt(os.path.join(bdir, "true_model.csv"), delimiter=",", skiprows=1)
    r = np.loadtxt(os.path.join(cfg["output_dir"], "4_inv", "vs_model.csv"), delimiter=",", skiprows=1)
    T, R = perturbation(t[:, :4]), perturbation(r[:, :4])
    zmax = float(R[2][-2])  # the deepest inverted node (the last one is the fixed half-space)
    layers, labels = [], []
    for name, (lon, lat, dep, V) in (("True", T), ("lant", R)):
        LO, LA = np.meshgrid(lon, lat, indexing="ij")
        for z in depths:
            kk = int(np.clip(np.searchsorted(dep, z), 1, len(dep) - 1))  # linear in depth between the nodes
            w = (z - dep[kk - 1]) / (dep[kk] - dep[kk - 1])
            layers.append((LO.ravel(), LA.ravel(), ((1 - w) * V[:, :, kk - 1] + w * V[:, :, kk]).ravel()))
            labels.append("%s, %g km" % (name, z))
    from scipy.interpolate import RegularGridInterpolator

    ft = RegularGridInterpolator(T[:3], T[3], bounds_error=False)
    lon, lat, dep, V = R
    LO, LA = np.meshgrid(lon, lat, indexing="ij")
    inside = hull.contains_points(np.c_[LO.ravel(), LA.ravel()]).reshape(LO.shape)
    with open(out + "_recovery.csv", "w") as fh:  # per inverted node depth: correlation and slope, lant on true
        fh.write("depth_km,correlation,slope\n")
        for k, z in enumerate(dep[:-1]):
            x = ft(np.c_[LO[inside], LA[inside], np.full(inside.sum(), z)])
            y = V[:, :, k][inside]
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() > 10 and np.std(x[ok]) > 0.5:
                fh.write(
                    "%.3f,%.3f,%.3f\n"
                    % (z, np.corrcoef(x[ok], y[ok])[0, 1], np.dot(x[ok], y[ok]) / np.dot(x[ok], x[ok]))
                )
    lines = maps.section_lines(sta)
    top = 5 * np.ceil(0.7 * np.nanmax(np.abs(T[3])) / 5)  # colour range: 70 % of the largest true perturbation
    maps.panels(cfg, layers, labels, sta, out + "_slices", "polar -I", "-%g/%g" % (top, top), "a%gf%g" % (top / 2,
                top / 10), "@~D@~Vs (%)",
                overlay=maps.line_overlay(lines))  # fmt: skip
    for name, (lon, lat, dep, V) in (("true", T), ("lant", R)):
        maps.sections(cfg, (lon, lat, dep, V), sta, lines, (-top, top), out + "_sections_" + name, zmax=zmax,
                      cpt="polar -I", label="@~D@~Vs (%)", contours=False)  # fmt: skip


if __name__ == "__main__":
    if sys.argv[1] == "dispersion":
        dispersion(sys.argv[2], sys.argv[3:])
    else:
        model(
            sys.argv[2], sys.argv[3], *([tuple(float(v) for v in sys.argv[4].split(","))] if len(sys.argv) > 4 else [])
        )
