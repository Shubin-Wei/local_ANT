#!/usr/bin/env python
"""Clean raw mseed data and prepare NoisePy input.

Scans raw miniSEED files, finds the target channel (default: Z component),
merges its segments, checks completeness against a full day, and writes
qualifying data in NoisePy naming format.
"""

import argparse
import os
import re
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

try:
    from obspy import Stream, UTCDateTime, read
except ImportError:
    sys.exit("Error: ObsPy is required. Install with: pip install obspy")

FILENAME_PAT = re.compile(r"^\d{12}$")
EXPECTED_SECONDS = 86400


def parse_args():
    p = argparse.ArgumentParser(
        description="Clean mseed data and convert to NoisePy format."
    )
    p.add_argument(
        "--input-root",
        default="local_ANT/TestData/RawData",
        help="Root directory of raw mseed files.",
    )
    p.add_argument(
        "--output-root",
        default="local_ANT/TestData/CleanData",
        help="Output directory for cleaned files.",
    )
    p.add_argument(
        "--time-zone",
        choices=["beijing", "utc"],
        default="beijing",
        help="How to interpret filename timestamps (default: beijing).",
    )
    p.add_argument(
        "--tolerance-seconds",
        type=float,
        default=60.0,
        help="Tolerance in seconds for full-day check (default: 60).",
    )
    p.add_argument(
        "--target-channel",
        default=None,
        help="Channel or component for completeness check and output, e.g. Z, BNZ, BHZ. "
        "If not set, defaults to the Z-component channel.",
    )
    p.add_argument(
        "--merge-fill",
        choices=["interpolate", "zero", "none"],
        default="interpolate",
        help="Gap fill strategy when merging traces (default: interpolate).",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print actions without writing output.",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files.",
    )
    return p.parse_args()


def parse_utc_from_name(name: str, tz: str) -> datetime:
    dt = datetime.strptime(name, "%y%m%d%H%M%S")
    if tz == "beijing":
        dt = dt - timedelta(hours=8)
    return dt.replace(tzinfo=timezone.utc)


def find_target_traces(st: Stream, target_spec) -> list:
    """Return all traces matching the target channel spec.

    target_spec can be a component letter (Z/N/E) or a full channel code (BNZ/BHZ).
    If not specified, defaults to Z-component traces.
    Returns a list of traces (may be multiple segments of the same channel).
    """
    if target_spec:
        spec = target_spec.upper()
        return [
            tr for tr in st
            if tr.stats.channel.upper() == spec
            or tr.stats.channel.upper().endswith(spec)
        ]
    return [tr for tr in st if tr.stats.channel.upper().endswith("Z")]


def is_full_day(trace, tolerance: float) -> tuple:
    """Check if a merged trace covers a full day.

    Returns (is_complete, duration_npts, reason).
    """
    sr = trace.stats.sampling_rate
    if sr <= 0:
        return False, 0.0, f"invalid sampling rate {sr}"

    if hasattr(trace.data, "mask"):
        valid_samples = int(np.sum(~trace.data.mask))
    else:
        valid_samples = trace.stats.npts

    duration_npts = valid_samples / sr
    diff = abs(duration_npts - EXPECTED_SECONDS)

    if diff > tolerance:
        return (
            False,
            duration_npts,
            f"duration {duration_npts:.1f}s, diff {diff:.1f}s > tolerance {tolerance:.1f}s",
        )
    return True, duration_npts, ""


def build_output_name(trace, year: int, doy: int) -> str:
    net = (trace.stats.network or "").strip().upper()
    sta = (trace.stats.station or "").strip().upper()
    cha = (trace.stats.channel or "").strip().upper()
    loc = (trace.stats.location or "").strip().upper()
    if len(loc) != 2:
        loc = "__"
    return f"{net}{sta}_{cha}{loc}_{year:04d}{doy:03d}.ms"


def main():
    args = parse_args()
    input_root = Path(args.input_root)
    output_root = Path(args.output_root)

    if not input_root.is_dir():
        sys.exit(f"Input directory does not exist: {input_root}")

    if not args.dry_run:
        output_root.mkdir(parents=True, exist_ok=True)

    merge_fill_value = None
    if args.merge_fill == "interpolate":
        merge_fill_value = "interpolate"
    elif args.merge_fill == "zero":
        merge_fill_value = 0

    stats = {
        "total_files": 0,
        "read_failed_files": 0,
        "skipped_bad_name_files": 0,
        "skipped_empty_files": 0,
        "skipped_no_target_channel": 0,
        "skipped_incomplete_channels": 0,
        "skipped_merge_failed_channels": 0,
        "written_channels": 0,
        "output_root": str(output_root),
    }

    for root, _, files in os.walk(input_root):
        root_path = Path(root)
        for fname in sorted(files):
            fpath = root_path / fname
            if not fpath.is_file():
                continue
            stats["total_files"] += 1
            stem = fpath.stem

            if not FILENAME_PAT.match(stem):
                print(f"[SKIP] {fpath} — filename not YYMMDDHHMMSS")
                stats["skipped_bad_name_files"] += 1
                continue

            try:
                t_utc = parse_utc_from_name(stem, args.time_zone)
            except Exception:
                print(f"[SKIP] {fpath} — cannot parse timestamp from filename")
                stats["skipped_bad_name_files"] += 1
                continue

            try:
                st = read(str(fpath))
            except Exception as e:
                print(f"[SKIP] {fpath} — read failed: {e}")
                stats["read_failed_files"] += 1
                continue

            if len(st) == 0:
                print(f"[SKIP] {fpath} — no traces in file")
                stats["skipped_empty_files"] += 1
                continue

            # Warn if filename time and trace starttime differ significantly
            try:
                trace_start = min(tr.stats.starttime for tr in st)
                fname_utc = UTCDateTime(t_utc)
                time_diff = abs(trace_start - fname_utc)
                if time_diff > 3600:
                    print(
                        f"[WARN] {fpath} — filename time and trace starttime "
                        f"differ by {time_diff:.0f}s"
                    )
            except Exception:
                pass

            # Find target channel traces (may be multiple segments)
            target_traces = find_target_traces(st, args.target_channel)
            if not target_traces:
                print(
                    f"[SKIP] {fpath} — target channel "
                    f"'{args.target_channel or 'Z'}' not found"
                )
                stats["skipped_no_target_channel"] += 1
                continue

            # Group target traces by full ID (handle multiple location codes)
            groups = {}
            for tr in target_traces:
                key = (
                    tr.stats.network,
                    tr.stats.station,
                    tr.stats.location,
                    tr.stats.channel,
                )
                groups.setdefault(key, []).append(tr)

            for key in sorted(groups):
                traces = groups[key]
                channel_st = Stream(traces)

                # Check sampling rate consistency
                srs = set(tr.stats.sampling_rate for tr in channel_st)
                if len(srs) > 1:
                    print(
                        f"[SKIP] {fpath} — inconsistent sampling rates: {srs}"
                    )
                    stats["skipped_merge_failed_channels"] += 1
                    continue

                # Merge segments
                try:
                    if args.merge_fill == "none":
                        channel_st.merge(method=1, fill_value=None)
                    else:
                        channel_st.merge(method=1, fill_value=merge_fill_value)
                except Exception as e:
                    print(f"[SKIP] {fpath} — merge failed: {e}")
                    stats["skipped_merge_failed_channels"] += 1
                    continue

                if len(channel_st) != 1:
                    print(
                        f"[SKIP] {fpath} — "
                        f"merge produced {len(channel_st)} traces (expected 1)"
                    )
                    stats["skipped_merge_failed_channels"] += 1
                    continue

                merged_tr = channel_st[0]

                # Determine date from actual trace starttime
                tr_start = merged_tr.stats.starttime
                year = tr_start.year
                doy = tr_start.julday

                # Check completeness
                complete, duration, reason = is_full_day(
                    merged_tr, args.tolerance_seconds
                )
                if not complete:
                    print(f"[SKIP] {fpath} — incomplete: {reason}")
                    stats["skipped_incomplete_channels"] += 1
                    continue

                # Build output path
                out_name = build_output_name(merged_tr, year, doy)
                out_path = output_root / out_name

                if out_path.exists() and not args.overwrite:
                    print(f"[SKIP] {fpath} -> {out_path} — already exists")
                    continue

                if args.dry_run:
                    print(f"[DRY]  {fpath} -> {out_path}")
                else:
                    Stream([merged_tr]).write(str(out_path), format="MSEED")
                    print(f"[OK]   {fpath} -> {out_path}")
                stats["written_channels"] += 1

    # Summary
    print("\n" + "=" * 60)
    print("Processing Summary")
    print("=" * 60)
    for k, v in stats.items():
        print(f"  {k:<35} {v}")
    print("=" * 60)


if __name__ == "__main__":
    main()
