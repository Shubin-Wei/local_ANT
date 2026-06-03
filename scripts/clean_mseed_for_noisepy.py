#!/usr/bin/env python
"""Clean raw mseed data and prepare NoisePy input.

Scans station folders under input-root, groups files by Beijing date,
merges target channel traces, checks completeness, and writes
qualifying data in NoisePy naming format.

Usage:
    python clean_mseed_for_noisepy.py --input-root RawData --output-root CleanData
    python clean_mseed_for_noisepy.py --input-root RawData --output-root CleanData --date-range 260316 260429
    python clean_mseed_for_noisepy.py --input-root RawData --output-root CleanData --verbose --dry-run
"""

import argparse
import os
import re
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore", message="readMSEEDBuffer", category=UserWarning)
warnings.filterwarnings("ignore", module="obspy.io.mseed")

try:
    from obspy import Stream, UTCDateTime, read
except ImportError:
    sys.exit("Error: ObsPy is required. Install with: pip install obspy")

FILENAME_PAT = re.compile(r"^\d{12}$")
EXPECTED_SECONDS = 86400
MAX_TRACES_PER_FILE = 30
TIME_DIFF_THRESHOLD = 86400

VERBOSE = False


def log(msg):
    if VERBOSE:
        print(f"    {msg}")


def parse_args():
    p = argparse.ArgumentParser(
        description="Clean mseed data and convert to NoisePy format."
    )
    p.add_argument("--input-root", default="local_ANT/TestData/RawData",
                    help="Root directory of raw mseed files (contains station folders).")
    p.add_argument("--output-root", default="local_ANT/TestData/CleanData",
                    help="Output directory for cleaned files.")
    p.add_argument("--min-completeness", type=float, default=0.8,
                    help="Min completeness ratio, 0~1 (default: 0.8)")
    p.add_argument("--target-channel", default=None,
                    help="Channel or component, e.g. Z, BNZ (default: Z)")
    p.add_argument("--date-range", nargs=2, metavar=("START", "END"),
                    help="Only process Beijing dates strictly between START and END "
                    "(both excluded). Format: YYMMDD, e.g. --date-range 260316 260429")
    p.add_argument("--overwrite", action="store_true",
                    help="Overwrite existing output files.")
    p.add_argument("--dry-run", action="store_true",
                    help="Print actions without writing output.")
    p.add_argument("--verbose", "-v", action="store_true",
                    help="Print detailed per-file messages.")
    return p.parse_args()


def find_target_traces(st: Stream, target_spec) -> list:
    if target_spec:
        spec = target_spec.upper()
        return [
            tr for tr in st
            if tr.stats.channel.upper() == spec
            or tr.stats.channel.upper().endswith(spec)
        ]
    return [tr for tr in st if tr.stats.channel.upper().endswith("Z")]


def check_completeness(trace, min_completeness: float) -> tuple:
    sr = trace.stats.sampling_rate
    if sr <= 0:
        return False, 0.0, f"invalid sampling rate {sr}"
    if hasattr(trace.data, "mask"):
        valid_samples = int(np.sum(~trace.data.mask))
    else:
        valid_samples = trace.stats.npts
    duration = valid_samples / sr
    ratio = duration / EXPECTED_SECONDS
    if ratio < min_completeness:
        return False, duration, f"completeness {ratio:.1%} < {min_completeness:.0%} ({duration:.0f}s / {EXPECTED_SECONDS}s)"
    return True, duration, ""


def build_output_name(trace, year: int, doy: int) -> str:
    net = (trace.stats.network or "").strip().upper()
    sta = (trace.stats.station or "").strip().upper()
    cha = (trace.stats.channel or "").strip().upper()
    loc = (trace.stats.location or "").strip().upper()
    if len(loc) != 2:
        loc = "__"
    return f"{net}{sta}_{cha}{loc}_{year:04d}{doy:03d}.ms"


def collect_files_by_beijing_date(folder, date_start, date_end, stats):
    date_files = {}
    for fname in sorted(os.listdir(folder)):
        fpath = folder / fname
        if not fpath.is_file() or fname == "PAR.INI":
            continue
        stats["total_files"] += 1
        stem = fpath.stem

        if not FILENAME_PAT.match(stem):
            stats["skipped_bad_name"] += 1
            continue

        beijing_date_str = stem[:6]
        try:
            beijing_date = datetime.strptime(beijing_date_str, "%y%m%d").date()
        except Exception:
            stats["skipped_bad_name"] += 1
            continue

        if date_start is not None:
            if beijing_date <= date_start or beijing_date >= date_end:
                stats["skipped_date_range"] += 1
                continue

        date_files.setdefault(beijing_date, []).append(fpath)
    return date_files


def read_and_filter_file(fpath, stats):
    try:
        st = read(str(fpath))
    except Exception as e:
        log(f"read failed: {fpath.name} — {e}")
        stats["read_failed"] += 1
        return None

    if len(st) == 0:
        log(f"empty: {fpath.name}")
        stats["skipped_empty"] += 1
        return None

    if len(st) > MAX_TRACES_PER_FILE:
        log(f"corrupt ({len(st)} traces): {fpath.name}")
        stats["skipped_too_many_traces"] += 1
        return None

    stem = fpath.stem
    beijing_dt = datetime.strptime(stem, "%y%m%d%H%M%S")
    fname_utc = UTCDateTime(beijing_dt.year, beijing_dt.month, beijing_dt.day,
                            beijing_dt.hour, beijing_dt.minute, beijing_dt.second) - 8 * 3600

    filtered = Stream()
    n_garbage = 0
    for tr in st:
        time_diff = abs(tr.stats.starttime - fname_utc)
        if time_diff > TIME_DIFF_THRESHOLD:
            n_garbage += 1
            continue
        filtered.append(tr)

    if n_garbage > 0:
        log(f"dropped {n_garbage} garbage trace(s): {fpath.name}")

    if len(filtered) == 0:
        stats["skipped_empty"] += 1
        return None

    return filtered


def process_day(day_files, station_label, args, stats):
    all_traces = Stream()
    for fpath in day_files:
        st = read_and_filter_file(fpath, stats)
        if st is not None:
            all_traces += st

    if len(all_traces) == 0:
        return None

    target_traces = find_target_traces(all_traces, args.target_channel)
    if not target_traces:
        log(f"{station_label} — no target channel '{args.target_channel or 'Z'}'")
        stats["skipped_no_channel"] += 1
        return None

    channel_st = Stream(target_traces)

    srs = set(tr.stats.sampling_rate for tr in channel_st)
    if len(srs) > 1:
        log(f"{station_label} — inconsistent sampling rates: {srs}")
        stats["skipped_merge_failed"] += 1
        return None

    try:
        channel_st.merge(method=1, fill_value=None)
    except Exception as e:
        log(f"{station_label} — merge failed: {e}")
        stats["skipped_merge_failed"] += 1
        return None

    if len(channel_st) != 1:
        log(f"{station_label} — merge produced {len(channel_st)} traces")
        stats["skipped_merge_failed"] += 1
        return None

    merged_tr = channel_st[0]
    tr_start = merged_tr.stats.starttime
    year = tr_start.year
    doy = tr_start.julday

    complete, duration, reason = check_completeness(merged_tr, args.min_completeness)
    if not complete:
        log(f"{station_label} — {reason}")
        stats["skipped_incomplete"] += 1
        return None

    if hasattr(merged_tr.data, "mask"):
        merged_tr.data = np.ma.filled(merged_tr.data, 0)

    out_name = build_output_name(merged_tr, year, doy)
    out_path = args._output_root / out_name

    if out_path.exists() and not args.overwrite:
        log(f"{out_name} — already exists")
        return None

    if args.dry_run:
        print(f"  [DRY] {out_name}  ({duration:.0f}s)")
    else:
        Stream([merged_tr]).write(str(out_path), format="MSEED")
        print(f"  [OK]  {out_name}  ({duration:.0f}s)")
    stats["written"] += 1
    return merged_tr


def main():
    args = parse_args()

    global VERBOSE
    VERBOSE = args.verbose

    input_root = Path(args.input_root)
    output_root = Path(args.output_root)
    args._output_root = output_root

    if not input_root.is_dir():
        sys.exit(f"Input directory does not exist: {input_root}")

    date_start = date_end = None
    if args.date_range:
        date_start = datetime.strptime(args.date_range[0], "%y%m%d").date()
        date_end = datetime.strptime(args.date_range[1], "%y%m%d").date()

    print(f"Input:  {input_root}")
    print(f"Output: {output_root}")
    print(f"Min completeness: {args.min_completeness:.0%}")
    if date_start:
        print(f"Date range (Beijing): ({date_start}, {date_end}) — both excluded")
    print()

    if not args.dry_run:
        output_root.mkdir(parents=True, exist_ok=True)

    stats = {
        "total_files": 0,
        "read_failed": 0,
        "skipped_bad_name": 0,
        "skipped_date_range": 0,
        "skipped_empty": 0,
        "skipped_too_many_traces": 0,
        "skipped_no_channel": 0,
        "skipped_incomplete": 0,
        "skipped_merge_failed": 0,
        "written": 0,
    }

    station_folders = sorted(
        d for d in input_root.iterdir()
        if d.is_dir()
    )

    for folder in station_folders:
        date_files = collect_files_by_beijing_date(folder, date_start, date_end, stats)
        if not date_files:
            continue

        n_days = len(date_files)
        n_files = sum(len(v) for v in date_files.values())
        print(f"{folder.name}  ({n_days} days, {n_files} files)")

        for day in sorted(date_files):
            label = f"{folder.name}/{day}"
            process_day(date_files[day], label, args, stats)

    print("\n" + "=" * 50)
    print("Summary")
    print("=" * 50)
    for k, v in stats.items():
        print(f"  {k:<28} {v}")
    print(f"  {'output':<28} {output_root}")
    print("=" * 50)


if __name__ == "__main__":
    main()
