#!/usr/bin/env python3
"""
Multi-panel cross-correlation moveout plots for journal publication.

Generates a 2×N figure (N = number of stack methods):
  - Top row: image (pcolormesh) moveout
  - Bottom row: wiggle with adaptive transparency, OR symmetrised image

Supports CLI arguments and YAML configuration.

Usage:
    python plot_stack_panels.py --stack-dir ./my_stack
    python plot_stack_panels.py --config config.yaml --methods linear pws robust
    python plot_stack_panels.py --stack-dir ./my_stack --bottom-style symmetric
"""

import argparse
import logging
import os
import tempfile
from fnmatch import fnmatch
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import h5py
import numpy as np
import yaml
from scipy.signal import butter, detrend, sosfiltfilt
from tqdm import tqdm

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

logging.getLogger("fontTools").setLevel(logging.WARNING)
logging.getLogger("fontTools.subset").setLevel(logging.WARNING)

LOG = logging.getLogger("plot_stack_panels")

JOURNAL_RCPARAMS = {
    "font.family": "serif",
    "font.serif": ["STIXGeneral", "Times New Roman", "Times"],
    "mathtext.fontset": "stix",
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.minor.size": 2,
    "ytick.minor.size": 2,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "axes.linewidth": 0.6,
    "lines.linewidth": 0.4,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}

METHOD_LABELS = {
    "linear": "Linear",
    "pws": "Phase-weighted",
    "robust": "Robust",
    "nroot": "Nth-root",
    "selective": "Selective",
    "auto_covariance": "Auto-cov",
}


# ──────────────── CLI / config ────────────────

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    parser.add_argument("--config", type=Path, default=None,
                        help="YAML config file. CLI flags override config values.")
    parser.add_argument("--stack-dir", type=Path, default=None,
                        help="Stack directory (NoisePy layout: stack_dir/SRC/REC/*.h5).")
    parser.add_argument("--component", default="ZZ", help="Component (default: ZZ).")
    parser.add_argument("--methods", nargs="+", default=["linear", "pws", "robust"],
                        help="Stack methods to compare as columns. Use 'all' to auto-discover (default: linear pws robust).")
    parser.add_argument("--station", default=None,
                        help="Comma-separated station filter. Supports glob (CF.*).")
    parser.add_argument("--stations", default=None,
                        help="Text file with one station code per line.")
    parser.add_argument("--freqmin", type=float, default=0.1, help="Bandpass low corner Hz.")
    parser.add_argument("--freqmax", type=float, default=0.2, help="Bandpass high corner Hz.")
    parser.add_argument("--min-dist", type=float, default=None, help="Min distance km.")
    parser.add_argument("--max-dist", type=float, default=None, help="Max distance km.")
    parser.add_argument("--max-records", type=int, default=None,
                        help="Max station pairs (shortest distances first).")
    parser.add_argument("--maxlag", type=float, default=None,
                        help="Half time window in seconds.")

    plot_grp = parser.add_argument_group("plot options")
    plot_grp.add_argument("--bottom-style", choices=["wiggle", "symmetric"],
                          default="wiggle",
                          help="Bottom row: wiggle traces or symmetrised image (default: wiggle).")
    plot_grp.add_argument("--wiggle-scale", type=float, default=50.0,
                          help="Wiggle amplitude multiplier relative to distance spacing (default: 50).")
    plot_grp.add_argument("--wiggle-every", type=int, default=1,
                          help="Draw one wiggle every N traces (default: 1).")
    plot_grp.add_argument("--cmap", default="RdBu_r", help="Colormap (default: RdBu_r).")
    plot_grp.add_argument("--clip-percentile", type=float, default=99.0,
                          help="Amplitude clip percentile (default: 99).")
    plot_grp.add_argument("--figure-width", type=float, default=7.2,
                          help="Figure width inches (default: 7.2).")
    plot_grp.add_argument("--figure-height", type=float, default=5.5,
                          help="Figure height inches (default: 5.5).")
    plot_grp.add_argument("--min-figure-height", type=float, default=5.4,
                          help="Min figure height in inches (default: 5.4).")
    plot_grp.add_argument("--max-figure-height", type=float, default=9.0,
                          help="Max figure height in inches (default: 9.0).")
    plot_grp.add_argument("--max-distance-ticks", type=int, default=8,
                          help="Max Y-axis ticks (default: 8).")
    plot_grp.add_argument("--max-time-ticks", type=int, default=5,
                          help="Max X-axis ticks (default: 5).")
    plot_grp.add_argument("--dpi", type=int, default=600, help="Raster DPI (default: 600).")
    plot_grp.add_argument("--no-title", action="store_true", help="Suppress suptitle.")

    parser.add_argument("--outdir", type=Path, default=None,
                        help="Output directory (default: figures/ next to stack-dir).")
    parser.add_argument("--output-format", nargs="+", default=["pdf", "png"],
                        choices=["png", "pdf", "svg"],
                        help="Output formats (default: pdf png).")
    parser.add_argument("--outname", default=None,
                        help="Output filename stem (auto-generated if omitted).")
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite existing output files.")
    parser.add_argument("--allow-missing-distance", action="store_true",
                        help="Plot records with missing distance using index order.")
    parser.add_argument("--validate-only", action="store_true",
                        help="Check data without generating figures.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO")
    return parser


VALID_CONFIG_KEYS = {
    "stack_dir", "component", "station", "stations",
    "freqmin", "freqmax", "min_dist", "max_dist", "max_records", "maxlag",
    "methods", "bottom_style",
    "wiggle_scale", "wiggle_every",
    "cmap", "clip_percentile",
    "figure_width", "figure_height",
    "min_figure_height", "max_figure_height",
    "max_distance_ticks", "max_time_ticks",
    "dpi", "no_title",
    "outdir", "output_format", "outname", "overwrite",
    "allow_missing_distance", "validate_only", "log_level",
}


def load_config(config_path: Path) -> Dict[str, Any]:
    resolved = config_path.expanduser().resolve()
    with resolved.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    if not isinstance(data, dict):
        raise SystemExit(f"Error: config must be a YAML mapping: {resolved}")

    norm = {k.strip().replace("-", "_"): v for k, v in data.items()}
    unknown = sorted(set(norm) - VALID_CONFIG_KEYS)
    if unknown:
        raise SystemExit(f"Error: unknown config keys in {resolved}: {', '.join(unknown)}")

    for pk in ("stack_dir", "outdir"):
        v = norm.get(pk)
        if v is None:
            continue
        p = Path(v).expanduser()
        if not p.is_absolute():
            p = resolved.parent / p
        norm[pk] = p.resolve()

    if "station" in norm and isinstance(norm["station"], list):
        norm["station"] = ",".join(str(s) for s in norm["station"])

    return norm


def parse_args() -> argparse.Namespace:
    parser = build_parser()
    pre, _ = parser.parse_known_args()
    if pre.config is not None:
        parser.set_defaults(**load_config(pre.config))
    args = parser.parse_args()

    if args.stack_dir is None:
        raise SystemExit("Error: --stack-dir is required (or set stack_dir in config).")
    if not Path(args.stack_dir).is_dir():
        raise SystemExit(f"Error: stack directory does not exist: {args.stack_dir}")
    if args.freqmin <= 0:
        raise SystemExit("Error: --freqmin must be > 0.")
    if args.freqmax <= args.freqmin:
        raise SystemExit("Error: --freqmax must be > --freqmin.")
    if args.outdir is None:
        args.outdir = Path(args.stack_dir).parent / "figures"
    return args


# ──────────────── Station matching ────────────────

def station_matches(name: str, patterns: List[str]) -> bool:
    for pat in patterns:
        if name == pat:
            return True
        if "." in name and name.split(".", 1)[1] == pat:
            return True
        if fnmatch(name, pat):
            return True
    return False


def parse_station_filter(station_str: Optional[str], stations_file: Optional[str]) -> Optional[List[str]]:
    result: List[str] = []
    if station_str is not None:
        for part in station_str.split(","):
            part = part.strip()
            if part:
                result.append(part)
    if stations_file is not None:
        path = Path(stations_file)
        if not path.is_file():
            raise SystemExit(f"Error: stations file not found: {path}")
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                result.append(line)
    return result if result else None


# ──────────────── I/O ────────────────

def _infer_stations(file_path: Path, h5: h5py.File, dataset_path: str) -> Tuple[str, str]:
    ds = h5[dataset_path]
    src_attr = ds.attrs.get("src") or ds.attrs.get("source")
    rec_attr = ds.attrs.get("rec") or ds.attrs.get("receiver")
    if src_attr and rec_attr:
        return str(src_attr), str(rec_attr)
    parts = file_path.parts
    if len(parts) >= 3:
        return parts[-3], parts[-2]
    return file_path.stem, file_path.stem


def load_all_methods(
    stack_dir: Path,
    methods: List[str],
    component: str,
    station_filter: Optional[List[str]] = None,
    min_dist: Optional[float] = None,
    max_dist: Optional[float] = None,
    max_records: Optional[int] = None,
    allow_missing_distance: bool = False,
) -> Dict[str, List[dict]]:
    """Load records for all requested methods in a single pass over h5 files."""
    stack_names = [f"Allstack_{m}" for m in methods]
    records: Dict[str, List[dict]] = {m: [] for m in methods}

    h5_files = sorted(stack_dir.rglob("*.h5"))
    if not h5_files:
        raise SystemExit(f"Error: no .h5 files in {stack_dir}")

    method_tag = ", ".join(methods)
    for fpath in tqdm(h5_files, desc=f"Loading [{method_tag}]", unit="file"):
        try:
            with h5py.File(fpath, "r") as h5:
                for method, sname in zip(methods, stack_names):
                    dp = f"AuxiliaryData/{sname}/{component}"
                    if dp not in h5:
                        continue
                    ds = h5[dp]
                    if "dt" not in ds.attrs or "maxlag" not in ds.attrs:
                        LOG.warning("Skipping %s: missing dt/maxlag", fpath)
                        continue

                    src, rec = _infer_stations(fpath, h5, dp)
                    if src == rec:
                        continue

                    records[method].append({
                        "src": src, "rec": rec,
                        "data": np.asarray(ds[()], dtype=np.float64),
                        "dt": float(ds.attrs["dt"]),
                        "maxlag": float(ds.attrs["maxlag"]),
                        "dist": float(ds.attrs.get("dist", np.nan)),
                    })
        except Exception as exc:
            LOG.warning("Failed to read %s: %s", fpath, exc)

    for method in methods:
        recs = records[method]
        if station_filter is not None:
            recs = [r for r in recs
                    if station_matches(r["src"], station_filter)
                    or station_matches(r["rec"], station_filter)]

        if allow_missing_distance:
            for idx, r in enumerate(recs):
                if not np.isfinite(r["dist"]):
                    r["dist"] = float(idx)
        else:
            recs = [r for r in recs if np.isfinite(r["dist"])]

        if min_dist is not None:
            recs = [r for r in recs if r["dist"] >= min_dist]
        if max_dist is not None:
            recs = [r for r in recs if r["dist"] <= max_dist]

        recs.sort(key=lambda r: (r["dist"], r["src"], r["rec"]))
        if max_records is not None:
            recs = recs[:max_records]
        records[method] = recs

    return records


def discover_stack_names(stack_dir: Path, component: str) -> List[str]:
    names: set = set()
    for fp in sorted(stack_dir.rglob("*.h5")):
        try:
            with h5py.File(fp, "r") as h5:
                aux = h5.get("AuxiliaryData")
                if aux is None:
                    continue
                for name in aux.keys():
                    if component in aux[name]:
                        names.add(str(name))
        except Exception:
            continue
    return sorted(names)


# ──────────────── Signal processing ────────────────

def bandpass(data: np.ndarray, dt: float, fmin: float, fmax: float) -> np.ndarray:
    nyquist = 0.5 / dt
    if fmax >= nyquist:
        raise ValueError(
            f"freqmax={fmax} Hz >= Nyquist={nyquist:.3f} Hz (dt={dt:.5f} s). "
            f"Lower --freqmax or check your data sampling rate."
        )
    sos = butter(4, [fmin, fmax], btype="bandpass", fs=1.0 / dt, output="sos")
    return sosfiltfilt(sos, detrend(data, type="constant"))


def prepare(
    records: List[dict], fmin: float, fmax: float, maxlag_plot: Optional[float],
    method_name: str = "",
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]:
    if not records:
        return None, None, None
    dt = records[0]["dt"]
    full_lag = records[0]["maxlag"]
    use_lag = full_lag if maxlag_plot is None else min(maxlag_plot, full_lag)
    npts = records[0]["data"].size
    tax = np.linspace(-full_lag, full_lag, npts)
    mask = np.abs(tax) <= use_lag + 0.5 * dt
    t = tax[mask]
    dists = np.array([r["dist"] for r in records])
    mat = np.zeros((len(records), t.size))
    desc = f"Filtering [{method_name}]" if method_name else "Filtering"
    for i, r in enumerate(tqdm(records, desc=desc, unit="trace")):
        tr = bandpass(r["data"], dt, fmin, fmax)[mask]
        pk = np.max(np.abs(tr))
        if pk > 0:
            tr /= pk
        mat[i] = tr
    return t, dists, mat


def symmetrise(
    t: np.ndarray, mat: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Fold causal + acausal: sym(t) = 0.5 * [f(t) + f(-t)]."""
    dt = t[1] - t[0]
    mask_pos = t >= -0.5 * dt
    mat_pos = mat[:, mask_pos]
    mask_neg = t <= 0.5 * dt
    mat_neg = mat[:, mask_neg][:, ::-1]
    n = min(mat_pos.shape[1], mat_neg.shape[1])
    half_sym = 0.5 * (mat_pos[:, :n] + mat_neg[:, :n])
    t_half = t[mask_pos][:n]
    t_full = np.concatenate([-t_half[:0:-1], t_half])
    mat_full = np.concatenate([half_sym[:, :0:-1], half_sym], axis=1)
    return t_full, mat_full


# ──────────────── Geometry helpers ────────────────

def time_edges(t: np.ndarray) -> np.ndarray:
    dt = t[1] - t[0]
    e = np.empty(len(t) + 1)
    e[1:-1] = 0.5 * (t[:-1] + t[1:])
    e[0] = t[0] - 0.5 * dt
    e[-1] = t[-1] + 0.5 * dt
    return e


def dist_edges(d: np.ndarray) -> np.ndarray:
    if len(d) == 1:
        h = max(abs(d[0]) * 0.1, 0.05)
        return np.array([d[0] - h, d[0] + h])
    m = 0.5 * (d[:-1] + d[1:])
    e = np.empty(len(d) + 1)
    e[1:-1] = m
    e[0] = d[0] - (d[1] - d[0]) * 0.5
    e[-1] = d[-1] + (d[-1] - d[-2]) * 0.5
    return e


def shared_ylim(prepared: Dict[str, tuple], methods: List[str]) -> List[float]:
    lo, hi = np.inf, -np.inf
    for m in methods:
        _, dists, _ = prepared[m]
        if dists is not None:
            lo = min(lo, dists.min())
            hi = max(hi, dists.max())
    pad = (hi - lo) * 0.03
    return [lo - pad, hi + pad]


def method_label(m: str) -> str:
    return METHOD_LABELS.get(m, m.replace("_", " ").title())


def freq_label(f: float) -> str:
    return f"{f:g}"


# ──────────────── Plotting ────────────────

def plot_panels(
    data_by_method: Dict[str, List[dict]],
    methods: List[str],
    args: argparse.Namespace,
) -> List[Path]:
    plt.rcParams.update(JOURNAL_RCPARAMS)

    ncols = len(methods)
    prepared = {}
    for m in methods:
        prepared[m] = prepare(data_by_method[m], args.freqmin, args.freqmax, args.maxlag, method_name=m)

    ylim = shared_ylim(prepared, methods)
    labels = [chr(ord("a") + i) for i in range(2 * ncols)]

    fig, axes = plt.subplots(
        2, ncols,
        figsize=(args.figure_width, args.figure_height),
        gridspec_kw={"height_ratios": [1, 1], "hspace": 0.28, "wspace": 0.08},
        squeeze=False,
    )

    if not args.no_title:
        fig.suptitle(
            f"{args.component}   {freq_label(args.freqmin)}–{freq_label(args.freqmax)} Hz"
            + (f"   lag ±{args.maxlag:.0f} s" if args.maxlag else ""),
            fontsize=11, fontweight="bold", y=0.995,
        )

    pcm_top = None
    pcm_bot = None

    for col, m in enumerate(methods):
        t, dists, mat = prepared[m]
        ax_top = axes[0, col]
        ax_bot = axes[1, col]

        if t is None:
            for ax in (ax_top, ax_bot):
                ax.text(0.5, 0.5, "No data", ha="center", va="center",
                        transform=ax.transAxes)
            continue

        n_pairs = len(dists)
        te = time_edges(t)
        de = dist_edges(dists)
        clip = np.percentile(np.abs(mat), args.clip_percentile) or 1.0

        # ── Top row: image ──
        pcm_top = ax_top.pcolormesh(
            te, de, mat, cmap=args.cmap,
            vmin=-clip, vmax=clip, shading="flat", rasterized=True,
        )
        ax_top.axvline(0, color="k", lw=0.4, alpha=0.5)
        ax_top.set_title(
            f"({labels[col]}) {method_label(m)}  ({n_pairs} pairs)",
            fontsize=9, pad=4,
        )
        ax_top.set_xlim(t[0], t[-1])
        ax_top.set_ylim(ylim)
        ax_top.xaxis.set_major_locator(ticker.MaxNLocator(args.max_time_ticks, integer=True))

        # ── Bottom row ──
        if args.bottom_style == "symmetric":
            t_sym, mat_sym = symmetrise(t, mat)
            te_sym = time_edges(t_sym)
            clip_sym = np.percentile(np.abs(mat_sym), args.clip_percentile) or 1.0
            pcm_bot = ax_bot.pcolormesh(
                te_sym, de, mat_sym, cmap=args.cmap,
                vmin=-clip_sym, vmax=clip_sym, shading="flat", rasterized=True,
            )
            ax_bot.set_title(
                f"({labels[col + ncols]}) {method_label(m)} (symmetric)",
                fontsize=9, pad=4,
            )
            ax_bot.set_xlim(t_sym[0], t_sym[-1])
        else:
            sel = slice(None, None, args.wiggle_every)
            dists_sel = dists[sel]
            mat_sel = mat[sel]
            uniq = np.unique(dists_sel)
            ref_sp = float(np.median(np.diff(uniq))) if len(uniq) > 1 else 0.1
            wscale = ref_sp * args.wiggle_scale if ref_sp > 0.001 else 0.1
            alpha = min(0.75, max(0.10, 55.0 / len(dists_sel)))

            for d, tr in zip(dists_sel, mat_sel):
                ax_bot.plot(t, d + tr * wscale, color="k", lw=0.4,
                            alpha=alpha, rasterized=True)
            ax_bot.axvline(0, color="r", lw=0.5, alpha=0.6)
            n_shown = len(dists_sel)
            ax_bot.set_title(
                f"({labels[col + ncols]}) {method_label(m)}  ({n_shown}/{n_pairs})",
                fontsize=9, pad=4,
            )
            ax_bot.set_xlim(t[0], t[-1])

        ax_bot.set_ylim(ylim)
        ax_bot.xaxis.set_major_locator(ticker.MaxNLocator(args.max_time_ticks, integer=True))

        if col > 0:
            ax_top.set_yticklabels([])
            ax_bot.set_yticklabels([])

    axes[0, 0].set_ylabel("Distance (km)")
    axes[1, 0].set_ylabel("Distance (km)")
    axes[1, ncols // 2].set_xlabel("Lag time (s)")

    if pcm_top is not None:
        cbar = fig.colorbar(pcm_top, ax=axes[0, :].tolist(), location="right",
                            fraction=0.025, pad=0.015, aspect=25)
        cbar.set_label("Normalized amplitude", fontsize=8)
        cbar.ax.tick_params(labelsize=7)

    if pcm_bot is not None and args.bottom_style == "symmetric":
        cbar2 = fig.colorbar(pcm_bot, ax=axes[1, :].tolist(), location="right",
                             fraction=0.025, pad=0.015, aspect=25)
        cbar2.set_label("Normalized amplitude", fontsize=8)
        cbar2.ax.tick_params(labelsize=7)

    # ── Save ──
    if args.outname is None:
        method_tag = "_".join(methods)
        stem = (f"{args.component}_{method_tag}"
                f"_{freq_label(args.freqmin)}-{freq_label(args.freqmax)}Hz"
                f"_{args.bottom_style}")
    else:
        stem = args.outname

    args.outdir.mkdir(parents=True, exist_ok=True)
    saved: List[Path] = []
    for fmt in args.output_format:
        out = args.outdir / f"{stem}.{fmt}"
        if out.exists() and not args.overwrite:
            raise SystemExit(
                f"Error: {out} already exists. Use --overwrite to replace."
            )
        fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
        saved.append(out)
    plt.close(fig)
    return saved


# ──────────────── Validate-only ────────────────

def run_validate_only(args: argparse.Namespace) -> None:
    stack_dir = Path(args.stack_dir)
    print(f"Validating: {stack_dir}")
    h5_files = sorted(stack_dir.rglob("*.h5"))
    print(f"  HDF5 files: {len(h5_files)}")
    if not h5_files:
        return

    all_names: set = set()
    all_comps: set = set()
    sample_dt = None
    for fp in tqdm(h5_files, desc="Scanning", unit="file"):
        try:
            with h5py.File(fp, "r") as h5:
                aux = h5.get("AuxiliaryData")
                if aux is None:
                    continue
                for sn in aux.keys():
                    all_names.add(str(sn))
                    for comp in aux[sn].keys():
                        all_comps.add(str(comp))
                        if sample_dt is None and "dt" in aux[sn][comp].attrs:
                            sample_dt = float(aux[sn][comp].attrs["dt"])
        except Exception as exc:
            print(f"  WARNING: {fp}: {exc}")

    print(f"  Stack methods: {sorted(all_names)}")
    print(f"  Components: {sorted(all_comps)}")
    if args.methods == ["all"]:
        print(f"  methods=all: will use all Allstack_* found above")
    else:
        for m in args.methods:
            full = f"Allstack_{m}"
            status = "FOUND" if full in all_names else "NOT FOUND"
            print(f"  {full}: {status}")
    status = "FOUND" if args.component in all_comps else "NOT FOUND"
    print(f"  Component '{args.component}': {status}")
    if sample_dt is not None:
        nyq = 0.5 / sample_dt
        ok = "OK" if args.freqmax < nyq else f"ERROR (Nyquist={nyq:.3f} Hz)"
        print(f"  dt={sample_dt}s, Nyquist={nyq:.3f} Hz, freq band: {ok}")
    print("Validation complete.")


# ──────────────── Main ────────────────

def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.validate_only:
        run_validate_only(args)
        return

    station_filter = parse_station_filter(args.station, args.stations)
    stack_dir = Path(args.stack_dir)

    if args.methods == ["all"]:
        discovered = discover_stack_names(stack_dir, args.component)
        args.methods = [
            name.replace("Allstack_", "")
            for name in discovered if name.startswith("Allstack_")
        ]
        if not args.methods:
            raise SystemExit(
                f"Error: no Allstack_* entries found in {stack_dir} for {args.component}. "
                f"Available: {discovered}"
            )
        print(f"Discovered methods: {args.methods}")

    print(f"Stack directory : {stack_dir}")
    print(f"Component       : {args.component}")
    print(f"Methods         : {', '.join(args.methods)}")
    print(f"Frequency band  : {args.freqmin}-{args.freqmax} Hz")
    if args.maxlag:
        print(f"Max lag         : ±{args.maxlag} s")
    if station_filter:
        print(f"Station filter  : {station_filter}")

    data = load_all_methods(
        stack_dir,
        args.methods,
        args.component,
        station_filter=station_filter,
        min_dist=args.min_dist,
        max_dist=args.max_dist,
        max_records=args.max_records,
        allow_missing_distance=args.allow_missing_distance,
    )

    any_data = False
    for m in args.methods:
        n = len(data[m])
        print(f"  {m}: {n} pairs")
        if n > 0:
            any_data = True
    if not any_data:
        available = discover_stack_names(stack_dir, args.component)
        raise SystemExit(
            f"Error: no records loaded for methods {args.methods}. "
            f"Available stack names: {available}"
        )

    print(f"Plotting ({args.bottom_style} bottom row) ...")
    saved = plot_panels(data, args.methods, args)
    print("Done.")
    for p in saved:
        print(f"  Saved: {p}")


if __name__ == "__main__":
    main()
