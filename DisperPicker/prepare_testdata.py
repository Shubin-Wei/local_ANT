#!/usr/bin/env python3
"""Prepare EGFAnalysisTimeFreq outputs for DisperPicker.

This script converts files like:

    G.STA1_STA2.dat
    C.STA1_STA2.dat
    GDisp.STA1_STA2.dat
    CDisp.T.STA1_STA2.dat

into the directory layout required by DisperPicker:

    <output>/
      group_image/STA1.STA2.dat
      phase_image/STA1.STA2.dat
      group_velocity/STA1.STA2.txt
      phase_velocity/STA1.STA2.txt

The velocity text outputs are optional and are only generated when the
corresponding source files exist.
"""

import argparse
import shutil
from pathlib import Path
from typing import Dict, Iterable, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reorganize EGFAnalysisTimeFreq outputs into the TestData layout "
            "required by DisperPicker."
        )
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing G.*, C.*, GDisp.*, and CDisp.T.* files.",
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help=(
            "Destination TestData directory. Subdirectories group_image, "
            "phase_image, group_velocity, and phase_velocity will be created."
        ),
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove existing output subdirectories before writing new files.",
    )
    return parser.parse_args()


def ensure_dirs(output_dir: Path, clean: bool) -> Dict[str, Path]:
    subdirs = {
        "group_image": output_dir / "group_image",
        "phase_image": output_dir / "phase_image",
        "group_velocity": output_dir / "group_velocity",
        "phase_velocity": output_dir / "phase_velocity",
    }
    for path in subdirs.values():
        if clean and path.exists():
            shutil.rmtree(path)
        path.mkdir(parents=True, exist_ok=True)
    return subdirs


def pair_to_name(pair: str, suffix: str) -> str:
    return pair.replace("_", ".") + suffix


def read_pick_curve(src: Path) -> Iterable[Tuple[float, float]]:
    with src.open("r", encoding="utf-8") as f:
        for i, line in enumerate(f):
            if i < 2:
                continue
            cols = line.split()
            if len(cols) < 2:
                continue
            yield float(cols[0]), float(cols[1])


def write_curve(src: Path, dst: Path) -> None:
    with dst.open("w", encoding="utf-8") as f:
        for period, velocity in read_pick_curve(src):
            f.write(f"{period:.3f}  {velocity:.3f}\n")


def collect_by_pair(files: Iterable[Path], prefix: str) -> Dict[str, Path]:
    mapping: Dict[str, Path] = {}
    for path in files:
        name = path.name
        if not name.startswith(prefix):
            continue
        pair = path.stem.split(".", 1)[1]
        mapping[pair] = path
    return mapping


def collect_phase_pick(files: Iterable[Path]) -> Dict[str, Path]:
    mapping: Dict[str, Path] = {}
    for path in files:
        name = path.name
        if not name.startswith("CDisp.T."):
            continue
        pair = path.stem.split(".", 2)[2]
        mapping[pair] = path
    return mapping


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()

    if not input_dir.is_dir():
        raise SystemExit(f"Input directory not found: {input_dir}")

    subdirs = ensure_dirs(output_dir, clean=args.clean)
    files = list(input_dir.iterdir())

    group_images = collect_by_pair(files, "G.")
    phase_images = collect_by_pair(files, "C.")
    pairs = sorted(set(group_images) & set(phase_images))

    copied_group = 0
    copied_phase = 0
    for pair in pairs:
        out_name = pair_to_name(pair, ".dat")
        shutil.copy2(group_images[pair], subdirs["group_image"] / out_name)
        shutil.copy2(phase_images[pair], subdirs["phase_image"] / out_name)
        copied_group += 1
        copied_phase += 1

    copied_group_labels = 0
    copied_phase_labels = 0
    group_labels = collect_by_pair(files, "GDisp.")
    phase_labels = collect_phase_pick(files)

    for pair, src in group_labels.items():
        if pair not in pairs:
            continue
        out_name = pair_to_name(pair, ".txt")
        write_curve(src, subdirs["group_velocity"] / out_name)
        copied_group_labels += 1

    for pair, src in phase_labels.items():
        if pair not in pairs:
            continue
        out_name = pair_to_name(pair, ".txt")
        write_curve(src, subdirs["phase_velocity"] / out_name)
        copied_phase_labels += 1

    print(f"input_dir: {input_dir}")
    print(f"output_dir: {output_dir}")
    print(f"paired_images: {len(pairs)}")
    print(f"group_image files: {copied_group}")
    print(f"phase_image files: {copied_phase}")
    print(f"group_velocity labels: {copied_group_labels}")
    print(f"phase_velocity labels: {copied_phase_labels}")


if __name__ == "__main__":
    main()
