import argparse
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import h5py
import numpy as np
from tqdm import tqdm


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT_ROOT = SCRIPT_DIR / "egfanalysis_input"
DEFAULT_COMPONENT = "ZZ"


def _apply_yaml_config_defaults(parser, argv):
    """Seed parser defaults from a --config YAML file; explicit CLI flags win.

    YAML keys use the same spelling as the long options (with or without
    dashes). Unknown keys are reported and ignored.
    """
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config")
    known, _ = pre.parse_known_args(argv)
    if not known.config:
        return
    import yaml
    with open(known.config, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    actions = {a.dest: a for a in parser._actions}
    defaults = {}
    for key, value in data.items():
        dest = key.replace("-", "_")
        if dest not in actions:
            print(f"WARNING: unknown key in {known.config} ignored: {key}")
            continue
        if value is None:
            continue
        action = actions[dest]
        if callable(action.type) and isinstance(value, str):
            value = action.type(value)
        defaults[dest] = value
        action.required = False
    parser.set_defaults(**defaults)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Export NoisePy stack HDF5 files to the text format required by "
            "EGFAnalysisTimeFreq. The script writes one .dat waveform per "
            "cross-station stack."
        )
    )
    parser.add_argument(
        "--stack-dir",
        type=Path,
        required=True,
        help="Stack directory to export from.",
    )
    parser.add_argument(
        "--stack-method",
        default="all",
        help=(
            "Stack method to export. Examples: linear, pws, robust, all, or full "
            "dataset name like Allstack_linear."
        ),
    )
    parser.add_argument(
        "--component",
        default=DEFAULT_COMPONENT,
        help="Stack component to export. Default: ZZ.",
    )
    parser.add_argument(
        "--station",
        default=None,
        help="Only export pairs where either source or receiver station matches this station code.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help=(
            "Root directory for exported EGFAnalysisTimeFreq inputs. Default: "
            "Data_noisepy/egfanalysis_input."
        ),
    )
    parser.add_argument(
        "--config",
        help=(
            "YAML file providing defaults for the options above "
            "(e.g. configs/04_export_egf.yaml); explicit CLI flags win."
        ),
    )
    _apply_yaml_config_defaults(parser, sys.argv[1:])
    return validate_args(parser.parse_args())


def validate_args(args: argparse.Namespace) -> argparse.Namespace:
    if not args.component.strip():
        raise ValueError("--component must not be empty.")
    if not args.stack_method.strip():
        raise ValueError("--stack-method must not be empty.")
    if args.station is not None:
        args.station = args.station.strip()
        if not args.station:
            raise ValueError("--station must not be empty when provided.")
    return args


def normalize_stack_name(stack_method: str) -> str:
    stack_method = stack_method.strip()
    if stack_method == "all":
        return "all"
    return stack_method if stack_method.startswith("Allstack_") else f"Allstack_{stack_method}"


def iter_stack_files(stack_dir: Path) -> Iterable[Path]:
    return sorted(stack_dir.rglob("*.h5"))


def discover_stack_names(stack_dir: Path, component: str) -> List[str]:
    if not stack_dir.exists():
        raise FileNotFoundError(f"Stack directory not found: {stack_dir}")
    names: Set[str] = set()
    for h5_path in tqdm(iter_stack_files(stack_dir), desc="Scanning stacks", unit="file"):
        with h5py.File(h5_path, "r") as h5:
            aux = h5.get("AuxiliaryData")
            if aux is None:
                continue
            for stack_name in aux.keys():
                if component in aux[stack_name]:
                    names.add(stack_name)
    return sorted(names)


def split_station_name(full_name: str) -> Tuple[str, str]:
    if "." in full_name:
        network, code = full_name.split(".", 1)
    else:
        network, code = "", full_name
    return network, code


def station_alias(full_name: str) -> str:
    network, code = split_station_name(full_name)
    return f"{network}{code}" if network else code



def read_dataset_attrs(dataset) -> Dict[str, float]:
    return {key: dataset.attrs[key] for key in dataset.attrs.keys()}


def build_time_columns(data: np.ndarray, dt: float, maxlag: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    npts = data.size
    if npts % 2 == 0:
        raise ValueError(f"Expected odd number of points for symmetric stack, got {npts}")

    expected_npts = int(round((2 * maxlag) / dt)) + 1
    if expected_npts != npts:
        raise ValueError(
            f"Inconsistent npts/dt/maxlag: npts={npts}, dt={dt}, maxlag={maxlag}, "
            f"expected_npts={expected_npts}"
        )

    center = npts // 2
    times = np.arange(center + 1, dtype=np.float64) * dt
    gab = np.asarray(data[center:], dtype=np.float64)
    gba = np.asarray(data[: center + 1][::-1], dtype=np.float64)
    return times, gab, gba


def write_dat_file(
    output_path: Path,
    attrs: Dict[str, float],
    times: np.ndarray,
    gab: np.ndarray,
    gba: np.ndarray,
) -> None:
    with output_path.open("w", encoding="utf-8") as f:
        f.write(f"{float(attrs['lonS']):16.7e} {float(attrs['latS']):16.7e} {float(attrs.get('elevS', 0)):16.7e}\n")
        f.write(f"{float(attrs['lonR']):16.7e} {float(attrs['latR']):16.7e} {float(attrs.get('elevR', 0)):16.7e}\n")
        for t, ab, ba in zip(times, gab, gba):
            f.write(f"{t:16.7e} {ab:16.7e} {ba:16.7e}\n")


def export_file(
    h5_path: Path,
    stack_name: str,
    component: str,
    output_dir: Path,
) -> Optional[Path]:
    src = h5_path.parent.parent.name
    rec = h5_path.parent.name
    if src == rec:
        return None

    dataset_path = f"AuxiliaryData/{stack_name}/{component}"
    with h5py.File(h5_path, "r") as h5:
        if dataset_path not in h5:
            return None
        dataset = h5[dataset_path]
        data = np.asarray(dataset[()], dtype=np.float64)
        attrs = read_dataset_attrs(dataset)

    dt = float(attrs["dt"])
    maxlag = float(attrs["maxlag"])
    times, gab, gba = build_time_columns(data, dt=dt, maxlag=maxlag)

    src_alias = station_alias(src)
    rec_alias = station_alias(rec)
    out_name = f"{src_alias}_{rec_alias}.dat"
    out_path = output_dir / out_name
    write_dat_file(out_path, attrs=attrs, times=times, gab=gab, gba=gba)

    return out_path


def export_dataset_stack(
    stack_dir: Path,
    stack_name: str,
    component: str,
    station: Optional[str],
    output_root: Path,
) -> Dict[str, object]:
    output_dir = output_root / stack_name / component
    output_dir.mkdir(parents=True, exist_ok=True)

    stack_files = list(iter_stack_files(stack_dir))
    total_files = len(stack_files)
    if total_files == 0:
        raise ValueError(f"No .h5 files found in stack directory: {stack_dir}")

    exported: List[Path] = []

    desc = f"Exporting {stack_name}/{component}"
    for h5_path in tqdm(stack_files, desc=desc, unit="file"):
        src = h5_path.parent.parent.name
        rec = h5_path.parent.name
        if src == rec:
            continue
        if station is not None and src != station and rec != station:
            continue

        out_path = export_file(
            h5_path=h5_path,
            stack_name=stack_name,
            component=component,
            output_dir=output_dir,
        )
        if out_path is not None:
            exported.append(out_path)

    if not exported:
        raise ValueError(
            f"No cross-station stacks were exported from {stack_dir} for {stack_name}/{component}."
        )

    return {
        "stack_name": stack_name,
        "component": component,
        "output_dir": str(output_dir),
        "exported": len(exported),
    }


def main() -> None:
    args = parse_args()
    stack_name_arg = normalize_stack_name(args.stack_method)

    summaries: List[Dict[str, object]] = []
    stack_dir = args.stack_dir
    if not stack_dir.exists():
        raise FileNotFoundError(f"Stack directory not found: {stack_dir}")

    if stack_name_arg == "all":
        stack_names = [name for name in discover_stack_names(stack_dir, args.component) if name.startswith("Allstack_")]
    else:
        stack_names = [stack_name_arg]

    if not stack_names:
        raise ValueError(f"No Allstack_* entries found in {stack_dir} for component {args.component}.")

    for stack_name in stack_names:
        summary = export_dataset_stack(
            stack_dir=stack_dir,
            stack_name=stack_name,
            component=args.component,
            station=args.station,
            output_root=args.output_root,
        )
        summaries.append(summary)

    total_pairs = sum(s['exported'] for s in summaries)
    print(f"{total_pairs} pairs exported to {args.output_root}/")
    name_width = max(len(f"{s['stack_name']}/{s['component']}") for s in summaries)
    for s in summaries:
        label = f"{s['stack_name']}/{s['component']}"
        print(f"  {label:<{name_width}}  {s['exported']} pairs")


if __name__ == "__main__":
    main()
