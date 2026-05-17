"""
Prepare EGFAnalysisTimeFreq output for DisperPicker input.

Reads C.*.dat / G.*.dat dispersion images from an EGFAnalysisTimeFreq output
directory, renames them to the layout DisperPicker expects, and updates
DisperPicker's config.py to match the actual matrix dimensions.

Usage:
    python prepare_disperpicker_input.py \\
        --input-dir  ../TestData/EGFAnalysisOutputData/RobustDisperImage \\
        --picker-dir ../DisperPicker

What it does:
    1. Scan input-dir for C.*.dat (phase) and G.*.dat (group) files
    2. Copy them to picker-dir/data/TestData/phase_image/ and group_image/
       with renamed filenames: C.STA1_STA2.dat -> STA1.STA2.dat
    3. Auto-detect matrix dimensions (rows=velocity, cols=period)
    4. Patch picker-dir/config/config.py so range_T and range_V match
"""

import argparse
import re
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare EGFAnalysisTimeFreq dispersion images for DisperPicker.",
    )
    parser.add_argument(
        "--input-dir", type=Path, required=True,
        help="Directory containing C.*.dat and G.*.dat from EGFAnalysisTimeFreq.",
    )
    parser.add_argument(
        "--picker-dir", type=Path, required=True,
        help="DisperPicker root directory.",
    )
    parser.add_argument(
        "--start-period", type=float, default=0.2,
        help="Start period in seconds (default: 0.2).",
    )
    parser.add_argument(
        "--end-period", type=float, default=5.0,
        help="End period in seconds (default: 5.0).",
    )
    parser.add_argument(
        "--min-velocity", type=float, default=0.5,
        help="Minimum velocity in km/s (default: 0.5).",
    )
    parser.add_argument(
        "--max-velocity", type=float, default=4.0,
        help="Maximum velocity in km/s (default: 4.0).",
    )
    parser.add_argument(
        "--station-file", type=Path, default=None,
        help="stations.txt to copy into DisperPicker config/.",
    )
    parser.add_argument(
        "--clean", action="store_true",
        help="Remove existing files in TestData before copying.",
    )
    return parser.parse_args()


def rename_stem(matlab_name: str) -> str:
    """C.STA1_STA2.dat -> STA1.STA2"""
    stem = Path(matlab_name).stem
    stem = re.sub(r'^[CG]\.', '', stem)
    return stem.replace('_', '.')


def detect_shape(dat_path: Path) -> tuple:
    """Read first and last line to get (nrows, ncols)."""
    with dat_path.open() as f:
        first_line = f.readline()
        ncols = len(first_line.split())
        nrows = 1
        for _ in f:
            nrows += 1
    return nrows, ncols


def copy_images(input_dir: Path, dest_dir: Path, prefix: str) -> list:
    """Copy {prefix}.*.dat files, return list of (src, dst) and detected shape."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    pairs = []
    for src in sorted(input_dir.glob(f"{prefix}.*.dat")):
        new_name = rename_stem(src.name) + ".dat"
        dst = dest_dir / new_name
        shutil.copy2(src, dst)
        pairs.append((src.name, new_name))
    return pairs


def patch_config(config_path: Path, n_period: int, n_velocity: int,
                 start_T: float, end_T: float, min_V: float, max_V: float,
                 picker_dir: Path) -> None:
    """Update range_T, range_V, and root in config.py."""
    text = config_path.read_text()

    text = re.sub(
        r"self\.range_T\s*=\s*\[.*?\]",
        f"self.range_T = [{start_T}, {end_T}, {n_period}]",
        text,
    )
    text = re.sub(
        r"self\.range_V\s*=\s*\[.*?\]",
        f"self.range_V = [{min_V}, {max_V}, {n_velocity}]",
        text,
    )
    text = re.sub(
        r"self\.root\s*=\s*'.*?'",
        f"self.root = '{picker_dir.resolve()}'",
        text,
    )

    config_path.write_text(text)


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir
    picker_dir = args.picker_dir

    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    if not picker_dir.is_dir():
        raise FileNotFoundError(f"DisperPicker directory not found: {picker_dir}")

    test_data = picker_dir / "data" / "TestData"

    if args.clean:
        for sub in ("phase_image", "group_image"):
            d = test_data / sub
            if d.exists():
                shutil.rmtree(d)

    phase_pairs = copy_images(input_dir, test_data / "phase_image", "C")
    group_pairs = copy_images(input_dir, test_data / "group_image", "G")

    if not phase_pairs and not group_pairs:
        print("No C.*.dat or G.*.dat files found.")
        return

    sample = None
    for sub in ("phase_image", "group_image"):
        d = test_data / sub
        if d.exists():
            files = sorted(d.glob("*.dat"))
            if files:
                sample = files[0]
                break

    if sample is None:
        print("No .dat files copied.")
        return

    n_velocity, n_period = detect_shape(sample)

    config_path = picker_dir / "config" / "config.py"
    if config_path.exists():
        patch_config(
            config_path, n_period, n_velocity,
            args.start_period, args.end_period,
            args.min_velocity, args.max_velocity,
            picker_dir,
        )

    if args.station_file:
        if not args.station_file.is_file():
            raise FileNotFoundError(f"Station file not found: {args.station_file}")
        dst = picker_dir / "config" / "station.txt"
        shutil.copy2(args.station_file, dst)

    print(f"{len(phase_pairs)} phase + {len(group_pairs)} group images -> {test_data}/")
    print(f"Matrix: {n_velocity} x {n_period}  (velocity x period)")
    print(f"config.py: range_T=[{args.start_period}, {args.end_period}, {n_period}]"
          f"  range_V=[{args.min_velocity}, {args.max_velocity}, {n_velocity}]")
    if args.station_file:
        print(f"station.txt -> {picker_dir / 'config' / 'station.txt'}")


if __name__ == "__main__":
    main()
