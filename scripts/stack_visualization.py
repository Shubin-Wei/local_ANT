"""Visualize NoisePy stack results.

Reads HDF5 stack files produced by NoisePy and generates moveout plots
(image and/or wiggle style). Supports CLI arguments and YAML configuration.

Expected HDF5 directory layout (NoisePy ASDFStackStore):
    stack_dir/SRC/REC/*.h5
    Each HDF5 contains: AuxiliaryData/Allstack_<method>/<component>
    Required dataset attributes: dt, maxlag
    Optional dataset attributes: dist, latS, lonS, latR, lonR

Dependencies: h5py, numpy, pyyaml, scipy, matplotlib.
noisepy is NOT required — it is only used for StackMethod enum if available.

Usage:
    python stack_visualization.py --stack-dir ./my_stack
    python stack_visualization.py --config config.yaml --freqmin 0.1 --freqmax 0.5
    python stack_visualization.py --stack-dir ./my_stack --validate-only
"""

import argparse
import logging
import os
import tempfile
from dataclasses import dataclass
from fnmatch import fnmatch
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import h5py
import numpy as np
import yaml
from scipy.signal import butter, detrend, sosfiltfilt

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

try:
    from noisepy.seis.io.datatypes import StackMethod
    STACK_METHOD_CHOICES = [method.value for method in StackMethod]
    STACK_METHOD_ALL = StackMethod.ALL.value
except ModuleNotFoundError:
    StackMethod = None
    STACK_METHOD_CHOICES = [
        "linear",
        "pws",
        "robust",
        "auto_covariance",
        "nroot",
        "selective",
        "all",
    ]
    STACK_METHOD_ALL = "all"

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

LOG = logging.getLogger("stack_visualization")


@dataclass
class StackRecord:
    src: str
    rec: str
    file_path: Path
    stack_name: str
    component: str
    data: np.ndarray
    dt: float
    maxlag: float
    dist: float


def is_autocorrelation(record: StackRecord) -> bool:
    return record.src == record.rec


def normalize_stack_method(name: str) -> str:
    """Accept both 'linear' and 'Allstack_linear' forms, return 'Allstack_linear'."""
    if name.lower() == "all":
        return "all"
    if name.startswith("Allstack_"):
        return name
    return f"Allstack_{name}"


def station_matches(src_or_rec: str, patterns: List[str]) -> bool:
    """Check if a station identifier matches any pattern.

    Matching rules:
    - Exact match: CF.17933
    - Station-only match (ignoring network prefix): 17933 matches CF.17933
    - Glob match: CF.* matches CF.17933
    """
    for pattern in patterns:
        if src_or_rec == pattern:
            return True
        if "." in src_or_rec and src_or_rec.split(".", 1)[1] == pattern:
            return True
        if fnmatch(src_or_rec, pattern):
            return True
    return False


def parse_station_arg(station_str: Optional[str], stations_file: Optional[str]) -> Optional[List[str]]:
    """Build station filter list from --station and/or --stations arguments."""
    result = []
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help=(
            "Optional YAML config file. Relative paths inside the YAML are resolved "
            "relative to that file. CLI flags override values from the config."
        ),
    )
    parser.add_argument(
        "--stack-dir",
        type=Path,
        default=None,
        help="Stack directory (NoisePy output layout: stack_dir/SRC/REC/*.h5).",
    )
    parser.add_argument("--component", default="ZZ", help="Component to plot (default: ZZ).")
    parser.add_argument(
        "--station",
        default=None,
        help=(
            "Station filter. Comma-separated list of station codes. "
            "Supports exact match (CF.17933), station-only (17933), or glob (CF.*)."
        ),
    )
    parser.add_argument(
        "--stations",
        default=None,
        help="Path to a text file with one station code per line.",
    )
    parser.add_argument("--freqmin", type=float, default=0.1, help="Bandpass low corner in Hz.")
    parser.add_argument("--freqmax", type=float, default=0.2, help="Bandpass high corner in Hz.")
    parser.add_argument(
        "--min-dist",
        type=float,
        default=None,
        help="Minimum inter-station distance in km.",
    )
    parser.add_argument(
        "--max-dist",
        type=float,
        default=None,
        help="Maximum inter-station distance in km.",
    )
    parser.add_argument(
        "--max-records",
        type=int,
        default=None,
        help="Maximum number of station pairs to plot (shortest distances first).",
    )
    parser.add_argument(
        "--maxlag",
        type=float,
        default=None,
        help="Half time window in seconds. Defaults to each dataset's full maxlag.",
    )
    parser.add_argument(
        "--stack-method",
        default="linear",
        help=(
            "Stack method name. Accepts both short (linear, pws, robust) and full "
            "(Allstack_linear) forms. Use 'all' to plot every Allstack_* result."
        ),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory for figures. Defaults to 'figures' next to --stack-dir.",
    )
    parser.add_argument(
        "--style",
        choices=["image", "wiggle", "both"],
        default="both",
        help="Plot style (default: both).",
    )
    parser.add_argument(
        "--output-format",
        choices=["png", "pdf", "svg"],
        default="png",
        help="Output figure format (default: png).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files without warning.",
    )

    plot_group = parser.add_argument_group("plot options")
    plot_group.add_argument(
        "--wiggle-scale",
        type=float,
        default=0.35,
        help="Waveform scale relative to median distance spacing (default: 0.35).",
    )
    plot_group.add_argument(
        "--wiggle-every",
        type=int,
        default=1,
        help="Draw one wiggle trace every N records (default: 1).",
    )
    plot_group.add_argument("--dpi", type=int, default=200, help="Figure DPI (default: 200).")
    plot_group.add_argument("--cmap", default="RdBu_r", help="Colormap for image plot (default: RdBu_r).")
    plot_group.add_argument(
        "--clip-percentile",
        type=float,
        default=99.0,
        help="Amplitude clip percentile for image plot (default: 99.0).",
    )
    plot_group.add_argument("--figure-width", type=float, default=10.5, help="Figure width in inches (default: 10.5).")
    plot_group.add_argument("--min-figure-height", type=float, default=5.4, help="Minimum figure height (default: 5.4).")
    plot_group.add_argument("--max-figure-height", type=float, default=9.0, help="Maximum figure height (default: 9.0).")
    plot_group.add_argument("--max-distance-ticks", type=int, default=8, help="Max Y-axis ticks (default: 8).")
    plot_group.add_argument("--max-time-ticks", type=int, default=5, help="Max X-axis ticks (default: 5).")
    plot_group.add_argument("--no-title", action="store_true", help="Suppress figure title.")

    parser.add_argument(
        "--allow-missing-distance",
        action="store_true",
        help="Plot records with missing distance using index order instead of filtering them out.",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Check data compatibility without generating figures.",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging verbosity (default: INFO).",
    )
    return parser


def normalize_config_key(key: str) -> str:
    return key.strip().replace("-", "_")


VALID_CONFIG_KEYS = {
    "stack_dir",
    "component",
    "station",
    "stations",
    "freqmin",
    "freqmax",
    "min_dist",
    "max_dist",
    "max_records",
    "maxlag",
    "stack_method",
    "outdir",
    "style",
    "output_format",
    "overwrite",
    "wiggle_scale",
    "wiggle_every",
    "dpi",
    "cmap",
    "clip_percentile",
    "figure_width",
    "min_figure_height",
    "max_figure_height",
    "max_distance_ticks",
    "max_time_ticks",
    "no_title",
    "allow_missing_distance",
    "validate_only",
    "log_level",
}


def load_config(config_path: Path) -> Dict[str, Any]:
    resolved_path = config_path.expanduser().resolve()
    with resolved_path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}

    if not isinstance(data, dict):
        raise ValueError(f"Config file must contain a YAML mapping at the top level: {resolved_path}")

    normalized = {normalize_config_key(key): value for key, value in data.items()}
    unknown_keys = sorted(set(normalized) - VALID_CONFIG_KEYS)
    if unknown_keys:
        raise ValueError(f"Unknown config keys in {resolved_path}: {', '.join(unknown_keys)}")

    for path_key in ("stack_dir", "outdir"):
        value = normalized.get(path_key)
        if value is None:
            continue
        path_value = Path(value).expanduser()
        if not path_value.is_absolute():
            path_value = resolved_path.parent / path_value
        normalized[path_key] = path_value.resolve()

    if "station" in normalized and isinstance(normalized["station"], list):
        normalized["station"] = ",".join(str(s) for s in normalized["station"])

    return normalized


def parse_args() -> argparse.Namespace:
    parser = build_parser()
    preliminary_args, _ = parser.parse_known_args()
    if preliminary_args.config is not None:
        parser.set_defaults(**load_config(preliminary_args.config))
    return validate_args(parser.parse_args())


def validate_args(args: argparse.Namespace) -> argparse.Namespace:
    if args.stack_dir is None:
        raise SystemExit(
            "Error: --stack-dir is required (or set stack_dir in YAML config).\n"
            "Example: python stack_visualization.py --stack-dir ./my_stack"
        )
    if not Path(args.stack_dir).is_dir():
        raise SystemExit(f"Error: stack directory does not exist: {args.stack_dir}")
    if args.freqmin <= 0:
        raise SystemExit("Error: --freqmin must be > 0.")
    if args.freqmax <= args.freqmin:
        raise SystemExit("Error: --freqmax must be larger than --freqmin.")
    if args.min_dist is not None and args.min_dist < 0:
        raise SystemExit("Error: --min-dist must be >= 0.")
    if args.max_dist is not None and args.max_dist < 0:
        raise SystemExit("Error: --max-dist must be >= 0.")
    if args.min_dist is not None and args.max_dist is not None and args.max_dist < args.min_dist:
        raise SystemExit("Error: --max-dist must be >= --min-dist.")
    if args.max_records is not None and args.max_records <= 0:
        raise SystemExit("Error: --max-records must be > 0 when provided.")
    if args.maxlag is not None and args.maxlag <= 0:
        raise SystemExit("Error: --maxlag must be > 0 when provided.")
    if args.wiggle_scale <= 0:
        raise SystemExit("Error: --wiggle-scale must be > 0.")
    if args.wiggle_every <= 0:
        raise SystemExit("Error: --wiggle-every must be > 0.")
    if args.dpi <= 0:
        raise SystemExit("Error: --dpi must be > 0.")
    if args.clip_percentile <= 0 or args.clip_percentile > 100:
        raise SystemExit("Error: --clip-percentile must be in (0, 100].")
    if args.outdir is None:
        args.outdir = Path(args.stack_dir).parent / "figures"
    return args


def _freq_label(freq: float) -> str:
    return f"{freq:g}"


def resolve_output_path(
    outdir: Path,
    stack_name: str,
    component: str,
    freqmin: float,
    freqmax: float,
    fmt: str = "png",
) -> Path:
    method = stack_method_label(stack_name)
    filename = f"{method}_{component}_{_freq_label(freqmin)}-{_freq_label(freqmax)}Hz.{fmt}"
    return outdir / filename


def _infer_station_names(file_path: Path, h5: h5py.File, dataset_path: str) -> Tuple[str, str]:
    """Infer source/receiver station names from HDF5 attrs or directory structure."""
    ds = h5[dataset_path]
    src_attr = ds.attrs.get("src") or ds.attrs.get("source")
    rec_attr = ds.attrs.get("rec") or ds.attrs.get("receiver")
    if src_attr and rec_attr:
        return str(src_attr), str(rec_attr)

    parts = file_path.parts
    if len(parts) >= 3:
        return parts[-3], parts[-2]

    return file_path.stem, file_path.stem


def load_records(
    stack_dir: Path,
    stack_name: str,
    component: str,
    station_filter: Optional[List[str]] = None,
    min_dist: Optional[float] = None,
    max_dist: Optional[float] = None,
    max_records: Optional[int] = None,
    allow_missing_distance: bool = False,
) -> Tuple[List[StackRecord], Dict[str, int]]:
    if not stack_dir.exists():
        raise FileNotFoundError(f"Stack directory not found: {stack_dir}")

    scan_stats: Dict[str, int] = {
        "h5_scanned": 0,
        "h5_read_ok": 0,
        "h5_read_error": 0,
        "h5_no_dataset": 0,
        "h5_missing_attrs": 0,
        "records_autocorrelation": 0,
        "records_station_filtered": 0,
        "records_distance_invalid": 0,
        "records_distance_filtered": 0,
        "records_final": 0,
    }

    records: List[StackRecord] = []
    missing_attrs_files: List[str] = []

    h5_files = sorted(stack_dir.rglob("*.h5"))
    scan_stats["h5_scanned"] = len(h5_files)

    for file_path in h5_files:
        try:
            with h5py.File(file_path, "r") as h5:
                dataset_path = f"AuxiliaryData/{stack_name}/{component}"
                if dataset_path not in h5:
                    scan_stats["h5_no_dataset"] += 1
                    LOG.debug("No dataset %s in %s", dataset_path, file_path)
                    continue

                dataset = h5[dataset_path]

                missing = []
                if "dt" not in dataset.attrs:
                    missing.append("dt")
                if "maxlag" not in dataset.attrs:
                    missing.append("maxlag")

                if missing:
                    scan_stats["h5_missing_attrs"] += 1
                    missing_attrs_files.append(f"{file_path} (missing: {', '.join(missing)})")
                    LOG.warning("Skipping %s: missing attributes %s", file_path, missing)
                    continue

                src, rec = _infer_station_names(file_path, h5, dataset_path)

                records.append(
                    StackRecord(
                        src=src,
                        rec=rec,
                        file_path=file_path,
                        stack_name=stack_name,
                        component=component,
                        data=np.asarray(dataset[()], dtype=np.float64),
                        dt=float(dataset.attrs["dt"]),
                        maxlag=float(dataset.attrs["maxlag"]),
                        dist=float(dataset.attrs.get("dist", np.nan)),
                    )
                )
                scan_stats["h5_read_ok"] += 1
        except Exception as e:
            scan_stats["h5_read_error"] += 1
            LOG.warning("Failed to read %s: %s", file_path, e)

    if missing_attrs_files:
        LOG.info("Files skipped due to missing attributes:\n  %s", "\n  ".join(missing_attrs_files))

    if not records:
        available = discover_stack_names(stack_dir, component)
        hint = f" Available stack names: {available}" if available else ""
        raise ValueError(
            f"No stacks found in {stack_dir} for {stack_name}/{component}.{hint}"
        )

    pre_auto = len(records)
    records = [r for r in records if not is_autocorrelation(r)]
    scan_stats["records_autocorrelation"] = pre_auto - len(records)
    if not records:
        raise ValueError(
            f"Only autocorrelation stacks found in {stack_dir} for {stack_name}/{component}."
        )

    if station_filter is not None:
        pre_station = len(records)
        records = [
            r for r in records
            if station_matches(r.src, station_filter) or station_matches(r.rec, station_filter)
        ]
        scan_stats["records_station_filtered"] = pre_station - len(records)
        if not records:
            raise ValueError(
                f"No stacks remain after station filtering with {station_filter}."
            )

    if allow_missing_distance:
        for idx, r in enumerate(records):
            if not np.isfinite(r.dist):
                r.dist = float(idx)
    else:
        pre_dist = len(records)
        records = [r for r in records if np.isfinite(r.dist)]
        scan_stats["records_distance_invalid"] = pre_dist - len(records)
        if not records:
            raise ValueError(
                f"All records have invalid distance. Use --allow-missing-distance to plot by index."
            )

    if min_dist is not None:
        pre = len(records)
        records = [r for r in records if r.dist >= min_dist]
        scan_stats["records_distance_filtered"] += pre - len(records)
    if max_dist is not None:
        pre = len(records)
        records = [r for r in records if r.dist <= max_dist]
        scan_stats["records_distance_filtered"] += pre - len(records)
    if not records:
        raise ValueError(
            f"No stacks remain after distance filtering [{min_dist}, {max_dist}] km."
        )

    records.sort(key=lambda r: (r.dist, r.src, r.rec))
    if max_records is not None:
        records = records[:max_records]

    scan_stats["records_final"] = len(records)
    return records, scan_stats


def discover_stack_names(stack_dir: Path, component: str) -> List[str]:
    if not stack_dir.exists():
        return []
    names = set()
    for file_path in sorted(stack_dir.rglob("*.h5")):
        try:
            with h5py.File(file_path, "r") as h5:
                aux = h5.get("AuxiliaryData")
                if aux is None:
                    continue
                for name in aux.keys():
                    if component in aux[name]:
                        names.add(str(name))
        except Exception:
            continue
    return sorted(names)


def design_filter(dt: float, freqmin: float, freqmax: float):
    sampling_rate = 1.0 / dt
    nyquist = 0.5 * sampling_rate
    if freqmax >= nyquist:
        raise ValueError(
            f"freqmax={freqmax} Hz exceeds Nyquist={nyquist:.3f} Hz (dt={dt:.5f} s). "
            f"Lower --freqmax or check your data sampling rate."
        )
    return butter(4, [freqmin, freqmax], btype="bandpass", fs=sampling_rate, output="sos")


def build_time_axis(npts: int, dt: float, maxlag: float) -> np.ndarray:
    return np.linspace(-maxlag, maxlag, npts, dtype=np.float64)


def prepare_matrix(
    records: List[StackRecord],
    freqmin: float,
    freqmax: float,
    requested_maxlag: Optional[float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    reference = records[0]
    npts = reference.data.size
    dt = reference.dt
    full_maxlag = reference.maxlag

    for record in records[1:]:
        if record.data.size != npts:
            raise ValueError("Inconsistent stack length across files.")
        if not np.isclose(record.dt, dt):
            raise ValueError("Inconsistent dt across files.")
        if not np.isclose(record.maxlag, full_maxlag):
            raise ValueError("Inconsistent maxlag across files.")

    time_axis = build_time_axis(npts, dt, full_maxlag)
    use_maxlag = full_maxlag if requested_maxlag is None else min(requested_maxlag, full_maxlag)
    time_mask = np.abs(time_axis) <= use_maxlag + 0.5 * dt
    filtered_time = time_axis[time_mask]

    sos = design_filter(dt, freqmin, freqmax)
    matrix = np.zeros((len(records), filtered_time.size), dtype=np.float64)
    distances = np.zeros(len(records), dtype=np.float64)

    for index, record in enumerate(records):
        trace = detrend(record.data, type="constant")
        trace = sosfiltfilt(sos, trace)
        trace = trace[time_mask]
        peak = np.max(np.abs(trace))
        if peak > 0:
            trace = trace / peak
        matrix[index] = trace
        distances[index] = record.dist

    return filtered_time, distances, matrix


def stack_method_label(stack_name: str) -> str:
    return stack_name.replace("Allstack_", "").replace("_", " ").strip() or stack_name


def compute_distance_edges(distances: np.ndarray) -> np.ndarray:
    distances = np.asarray(distances, dtype=np.float64)
    if distances.size == 1:
        half = 0.5 if distances[0] == 0 else max(abs(distances[0]) * 0.1, 0.1)
        return np.array([distances[0] - half, distances[0] + half], dtype=np.float64)

    mids = 0.5 * (distances[:-1] + distances[1:])
    first_half = max((distances[1] - distances[0]) * 0.5, 1e-6)
    last_half = max((distances[-1] - distances[-2]) * 0.5, 1e-6)
    edges = np.empty(distances.size + 1, dtype=np.float64)
    edges[1:-1] = mids
    edges[0] = distances[0] - first_half
    edges[-1] = distances[-1] + last_half
    return edges


def choose_distance_ticks(distances: np.ndarray, max_ticks: int = 8) -> Tuple[np.ndarray, List[str]]:
    tick_rows = np.linspace(0, len(distances) - 1, min(max_ticks, len(distances)), dtype=int)
    tick_values = distances[tick_rows]
    tick_labels = [f"{value:.2f}" for value in tick_values]
    return tick_values, tick_labels


def compute_wiggle_scale(distances: np.ndarray, wiggle_scale: float) -> float:
    unique_distances = np.unique(distances)
    if len(unique_distances) > 1:
        reference_spacing = float(np.median(np.diff(unique_distances)))
    else:
        reference_spacing = 0.1
    return max(reference_spacing * wiggle_scale, 0.02)


def build_summary_title(
    stack_name: str,
    component: str,
    freqmin: float,
    freqmax: float,
    record_count: int,
) -> str:
    method = stack_method_label(stack_name)
    return f"{method}  {component}  {_freq_label(freqmin)}–{_freq_label(freqmax)} Hz  ({record_count} pairs)"


def choose_time_ticks(time_axis: np.ndarray, max_ticks: int = 5) -> np.ndarray:
    if time_axis.size == 1:
        return np.asarray([time_axis[0]], dtype=np.float64)
    return np.linspace(time_axis[0], time_axis[-1], max_ticks)


def compute_time_edges(time_axis: np.ndarray) -> np.ndarray:
    if time_axis.size == 1:
        half = 0.5
        return np.array([time_axis[0] - half, time_axis[0] + half], dtype=np.float64)

    dt = time_axis[1] - time_axis[0]
    edges = np.empty(time_axis.size + 1, dtype=np.float64)
    edges[1:-1] = 0.5 * (time_axis[:-1] + time_axis[1:])
    edges[0] = time_axis[0] - 0.5 * dt
    edges[-1] = time_axis[-1] + 0.5 * dt
    return edges


def _check_overwrite(output_path: Path, overwrite: bool) -> None:
    if output_path.exists() and not overwrite:
        raise SystemExit(
            f"Error: output file already exists: {output_path}\n"
            f"Use --overwrite to replace, or delete it manually."
        )


def plot_stack_image(
    time_axis: np.ndarray,
    distances: np.ndarray,
    matrix: np.ndarray,
    stack_name: str,
    component: str,
    freqmin: float,
    freqmax: float,
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    clip = np.percentile(np.abs(matrix), args.clip_percentile)
    if clip == 0:
        clip = 1.0
    time_edges = compute_time_edges(time_axis)
    distance_edges = compute_distance_edges(distances)
    xticks = choose_time_ticks(time_axis, args.max_time_ticks)
    tick_values, tick_labels = choose_distance_ticks(distances, args.max_distance_ticks)

    fig_height = max(args.min_figure_height, min(args.max_figure_height, 0.11 * matrix.shape[0]))
    fig, ax = plt.subplots(figsize=(args.figure_width, fig_height))

    if not args.no_title:
        title = build_summary_title(stack_name, component, freqmin, freqmax, len(distances))
        fig.suptitle(title, fontsize=12, y=0.97)

    mesh = ax.pcolormesh(
        time_edges,
        distance_edges,
        matrix,
        cmap=args.cmap,
        vmin=-clip,
        vmax=clip,
        shading="flat",
    )
    ax.axvline(0.0, color="k", lw=0.8, alpha=0.4)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("distance [km]")
    ax.set_xlim(time_edges[0], time_edges[-1])
    ax.set_xticks(xticks)
    ax.set_yticks(tick_values)
    ax.set_yticklabels(tick_labels)
    ax.set_ylim(distance_edges[0], distance_edges[-1])
    ax.margins(x=0)
    ax.xaxis.set_ticks_position("bottom")
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(f"{component} amplitude")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0, 0, 1, 0.95) if not args.no_title else (0, 0, 1, 1))
    fig.savefig(output_path, dpi=args.dpi)
    plt.close(fig)


def plot_wiggle(
    time_axis: np.ndarray,
    distances: np.ndarray,
    matrix: np.ndarray,
    component: str,
    stack_name: str,
    freqmin: float,
    freqmax: float,
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    distance_edges = compute_distance_edges(distances)
    tick_values, tick_labels = choose_distance_ticks(distances, args.max_distance_ticks)
    offset_scale = compute_wiggle_scale(distances, args.wiggle_scale)
    wiggle_margin = max(offset_scale * 2.4, 0.05)
    xticks = choose_time_ticks(time_axis, args.max_time_ticks)

    fig_height = max(args.min_figure_height, min(args.max_figure_height, 0.13 * matrix.shape[0]))
    fig, ax = plt.subplots(figsize=(args.figure_width, fig_height))

    if not args.no_title:
        title = build_summary_title(stack_name, component, freqmin, freqmax, len(distances))
        fig.suptitle(title, fontsize=12, y=0.97)

    plot_distances = distances[::args.wiggle_every]
    plot_matrix = matrix[::args.wiggle_every]
    for dist in plot_distances:
        ax.hlines(dist, time_axis[0], time_axis[-1], color="0.88", lw=0.45, zorder=0)
    for dist, trace in zip(plot_distances, plot_matrix):
        ax.plot(time_axis, dist + trace * offset_scale, color="0.05", lw=0.8, zorder=2)

    ax.axvline(0.0, color="k", lw=0.8, alpha=0.4)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("distance [km]")
    ax.set_xlim(time_axis[0], time_axis[-1])
    ax.set_ylim(distance_edges[0] - wiggle_margin, distance_edges[-1] + wiggle_margin)
    ax.set_yticks(tick_values)
    ax.set_yticklabels(tick_labels)
    ax.set_xticks(xticks)
    ax.margins(x=0)
    ax.xaxis.set_ticks_position("bottom")
    ax.grid(True, axis="x", alpha=0.16)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0, 0, 1, 0.95) if not args.no_title else (0, 0, 1, 1))
    fig.savefig(output_path, dpi=args.dpi)
    plt.close(fig)


def add_suffix(path: Path, suffix: str) -> Path:
    return path.with_name(f"{path.stem}_{suffix}{path.suffix}")


def resolve_requested_stack_names(args: argparse.Namespace, stack_dir: Path) -> List[str]:
    method = normalize_stack_method(args.stack_method)
    if method == "all":
        stack_names = [
            name for name in discover_stack_names(stack_dir, args.component)
            if name.startswith("Allstack_")
        ]
        if not stack_names:
            available = discover_stack_names(stack_dir, args.component)
            raise ValueError(
                f"No Allstack_* entries found in {stack_dir} for component {args.component}. "
                f"Available names: {available}"
            )
        return stack_names
    return [method]


def resolve_run_output_path(
    args: argparse.Namespace,
    stack_name: str,
) -> Path:
    return resolve_output_path(
        args.outdir,
        stack_name,
        args.component,
        args.freqmin,
        args.freqmax,
        fmt=args.output_format,
    )


def run_validate_only(args: argparse.Namespace) -> None:
    stack_dir = Path(args.stack_dir)
    station_filter = parse_station_arg(args.station, args.stations)

    print(f"Validating stack directory: {stack_dir}")
    h5_files = sorted(stack_dir.rglob("*.h5"))
    print(f"  HDF5 files found: {len(h5_files)}")
    if not h5_files:
        print("  ERROR: No .h5 files found.")
        return

    all_stack_names = set()
    all_components = set()
    files_with_aux = 0
    files_without_aux = 0
    sample_dt = None
    sample_maxlag = None
    attrs_issues = []

    for fp in h5_files:
        try:
            with h5py.File(fp, "r") as h5:
                aux = h5.get("AuxiliaryData")
                if aux is None:
                    files_without_aux += 1
                    continue
                files_with_aux += 1
                for sname in aux.keys():
                    all_stack_names.add(str(sname))
                    for comp in aux[sname].keys():
                        all_components.add(str(comp))
                        ds = aux[sname][comp]
                        missing = [a for a in ("dt", "maxlag") if a not in ds.attrs]
                        if missing:
                            attrs_issues.append(f"{fp}: {sname}/{comp} missing {missing}")
                        else:
                            if sample_dt is None:
                                sample_dt = float(ds.attrs["dt"])
                                sample_maxlag = float(ds.attrs["maxlag"])
        except Exception as e:
            print(f"  WARNING: Cannot read {fp}: {e}")

    print(f"  Files with AuxiliaryData: {files_with_aux}")
    print(f"  Files without AuxiliaryData: {files_without_aux}")
    print(f"  Stack methods found: {sorted(all_stack_names)}")
    print(f"  Components found: {sorted(all_components)}")

    if attrs_issues:
        print(f"  Attribute issues ({len(attrs_issues)}):")
        for issue in attrs_issues[:10]:
            print(f"    {issue}")
        if len(attrs_issues) > 10:
            print(f"    ... and {len(attrs_issues) - 10} more")

    requested_method = normalize_stack_method(args.stack_method)
    if requested_method != "all":
        if requested_method in all_stack_names:
            print(f"  Requested stack method '{requested_method}': FOUND")
        else:
            print(f"  Requested stack method '{requested_method}': NOT FOUND")
            print(f"    Available: {sorted(all_stack_names)}")

    if args.component in all_components:
        print(f"  Requested component '{args.component}': FOUND")
    else:
        print(f"  Requested component '{args.component}': NOT FOUND")
        print(f"    Available: {sorted(all_components)}")

    if sample_dt is not None:
        nyquist = 0.5 / sample_dt
        print(f"  Sample dt={sample_dt}s, Nyquist={nyquist:.3f} Hz")
        if args.freqmax >= nyquist:
            print(f"  ERROR: freqmax={args.freqmax} Hz >= Nyquist={nyquist:.3f} Hz")
        else:
            print(f"  Frequency band {args.freqmin}-{args.freqmax} Hz: OK")

    if station_filter:
        print(f"  Station filter: {station_filter}")

    print("Validation complete.")


def render_outputs(
    args: argparse.Namespace,
    stack_dir: Path,
    stack_name: str,
    output_path: Path,
    station_filter: Optional[List[str]],
) -> Tuple[List[Path], int, np.ndarray, np.ndarray, Dict[str, int]]:
    records, scan_stats = load_records(
        stack_dir,
        stack_name,
        args.component,
        station_filter=station_filter,
        min_dist=args.min_dist,
        max_dist=args.max_dist,
        max_records=args.max_records,
        allow_missing_distance=args.allow_missing_distance,
    )
    time_axis, distances, matrix = prepare_matrix(records, args.freqmin, args.freqmax, args.maxlag)

    saved_paths: List[Path] = []
    if args.style in ("image", "both"):
        img_path = add_suffix(output_path, "image") if args.style == "both" else output_path
        _check_overwrite(img_path, args.overwrite)
        plot_stack_image(
            time_axis, distances, matrix,
            stack_name, args.component, args.freqmin, args.freqmax,
            img_path, args,
        )
        saved_paths.append(img_path)

    if args.style in ("wiggle", "both"):
        wig_path = add_suffix(output_path, "wiggle") if args.style == "both" else output_path
        _check_overwrite(wig_path, args.overwrite)
        plot_wiggle(
            time_axis, distances, matrix,
            args.component, stack_name, args.freqmin, args.freqmax,
            wig_path, args,
        )
        saved_paths.append(wig_path)

    return saved_paths, len(records), distances, time_axis, scan_stats


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

    station_filter = parse_station_arg(args.station, args.stations)
    stack_dir = args.stack_dir
    stack_names = resolve_requested_stack_names(args, stack_dir)

    for stack_name in stack_names:
        output_path = resolve_run_output_path(args, stack_name)
        saved_paths, record_count, distances, time_axis, scan_stats = render_outputs(
            args, stack_dir, stack_name, output_path, station_filter
        )
        LOG.debug("Scan: %s", scan_stats)
        for p in saved_paths:
            print(f"Saved: {p}  ({record_count} pairs, {_freq_label(args.freqmin)}-{_freq_label(args.freqmax)} Hz)")


if __name__ == "__main__":
    main()
