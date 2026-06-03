#!/usr/bin/env python
"""Local NoisePy cross-correlation and stacking workflow.

A generic CLI tool for running ambient noise cross-correlation on flat
miniSEED datasets using NoisePy. Supports YAML configuration, per-stage
execution, and flexible input formats.

Usage:
    python local_cc_workflow.py --config config.yaml
    python local_cc_workflow.py --data-dir ./data --station-csv ./station.csv
    python local_cc_workflow.py --config config.yaml --stage stack
    python local_cc_workflow.py --config config.yaml --validate-only
"""

import argparse
import json
import logging
import os
import re
import statistics
import tempfile
import threading
from dataclasses import asdict, dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

try:
    import numpy as np
    import obspy
    import pandas as pd
    import yaml
    from datetimerange import DateTimeRange
except ModuleNotFoundError as exc:
    raise SystemExit(
        f"Missing required Python package: {exc.name}. "
        "Install with: pip install obspy numpy pandas pyyaml datetimerange"
    ) from exc

from noisepy.seis import __version__, cross_correlate, noise_module
from noisepy.seis.io.asdfstore import ASDFCCStore, ASDFStackStore
from noisepy.seis.io.channelcatalog import CSVChannelCatalog
from noisepy.seis.io.channel_filter_store import channel_filter
from noisepy.seis.io.datatypes import (
    CCMethod,
    Channel,
    ChannelData,
    ChannelType,
    ConfigParameters,
    FreqNorm,
    RmResp,
    Stack,
    StackMethod,
    Station,
    TimeNorm,
)
from noisepy.seis.io.mseedstore import MiniSeedDataStore
from noisepy.seis.stack import calc_segments, validate_pairs

LOG = logging.getLogger("local_cc_workflow")

FILENAME_RE = re.compile(
    r"^(?P<net>\w{2})(?P<sta>\w+?)_(?P<cha>\w{3})(?P<loc>\w{2})_(?P<yj>\d{7})\.(?P<ext>ms|mseed|miniseed)$"
)

STATION_CSV_ALIASES = {
    "net": "network",
    "sta": "station",
    "chan": "channel",
    "lat": "latitude",
    "lon": "longitude",
    "lng": "longitude",
    "elev": "elevation",
    "loc": "location",
}

REQUIRED_STATION_COLS = ["network", "station", "latitude", "longitude", "elevation"]


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class DatasetPaths:
    data_dir: Path
    station_csv: Path
    ccf_dir: Path
    stack_dir: Path


@dataclass
class ProcessingOptions:
    inc_hours: int = 24
    cc_len: int = 3600
    step: float = 1800.0
    maxlag: int = 100
    sampling_rate: float = 20.0
    freqmin: float = 1.0
    freqmax: float = 5.0
    max_over_std: int = 10
    smoothspect_n: int = 10
    smooth_n: int = 10
    substack: bool = True
    substack_windows: int = 1
    keep_substack: bool = False
    station_pattern: str = "*"
    stack_method: str = StackMethod.ALL.value
    time_norm: str = TimeNorm.ONE_BIT.value
    freq_norm: str = FreqNorm.RMA.value
    rm_resp: str = RmResp.NO.value
    ncomp: int = 1
    rotation: bool = False
    xcorr_only: bool = True
    acorr_only: bool = False


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def yj_to_datetime(yj: str) -> datetime:
    year = int(yj[:4])
    doy = int(yj[4:])
    return datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(days=doy - 1)


def is_cross_station_pair(src_sta: Station, rec_sta: Station) -> bool:
    return src_sta != rec_sta


def infer_chunk_offset(ms_files: List[Path]) -> timedelta:
    offsets = []
    for path in ms_files[:20]:
        match = FILENAME_RE.match(path.name)
        if not match:
            continue
        try:
            trace = obspy.read(str(path), headonly=True)[0]
        except Exception:
            continue
        file_day = yj_to_datetime(match.group("yj"))
        trace_start = trace.stats.starttime.datetime.replace(tzinfo=timezone.utc)
        offsets.append((trace_start - file_day).total_seconds() % 86400)

    if not offsets:
        return timedelta(0)

    offset_seconds = int(round(statistics.median(offsets) / 60.0) * 60)
    return timedelta(seconds=offset_seconds)


# ---------------------------------------------------------------------------
# Station CSV handling
# ---------------------------------------------------------------------------


def normalize_station_csv(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = [c.strip().lower() for c in df.columns]
    df.rename(columns=STATION_CSV_ALIASES, inplace=True)

    missing = [col for col in REQUIRED_STATION_COLS if col not in df.columns]
    if missing:
        raise ValueError(
            f"station.csv is missing required columns: {missing}. "
            f"Available columns: {list(df.columns)}. "
            f"Supported aliases: {STATION_CSV_ALIASES}"
        )

    for col in ["network", "station", "channel"]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()

    if "location" in df.columns:
        df["location"] = df["location"].fillna("").astype(str).str.strip()

    return df


class LocalCSVChannelCatalog(CSVChannelCatalog):
    def __init__(self, file: str):
        self.df = normalize_station_csv(
            pd.read_csv(file, dtype={"network": str, "station": str, "channel": str})
        )

    @staticmethod
    def _normalize_location(value: Any) -> str:
        if value is None:
            return ""
        text = str(value).strip()
        return "" if text.lower() == "nan" else text

    def get_full_channel(self, timespan: DateTimeRange, ch: Channel) -> Channel:
        location = self._normalize_location(ch.station.location)
        rows = self.df[
            (self.df["network"] == ch.station.network)
            & (self.df["station"] == ch.station.name)
            & (self.df["channel"] == ch.type.name)
        ]
        if "location" in self.df.columns:
            rows = rows[rows["location"] == location]
        if rows.empty:
            raise ValueError(
                f"Station/channel not found in station.csv: "
                f"{ch.station.network}.{ch.station.name}.{ch.type.name}"
            )
        row = rows.iloc[0]
        return Channel(
            ch.type,
            Station(
                network=ch.station.network,
                name=ch.station.name,
                lat=row["latitude"],
                lon=row["longitude"],
                elevation=row["elevation"],
                location=ch.station.location,
            ),
        )


# ---------------------------------------------------------------------------
# MiniSEED data store
# ---------------------------------------------------------------------------


class FlatMiniSeedStore(MiniSeedDataStore):
    """Read flat miniSEED files named like NETSTA_CHALOC_YYYYDOY.{ms,mseed,miniseed}."""

    def __init__(
        self,
        path: Path,
        chan_catalog: LocalCSVChannelCatalog,
        chunk_offset: timedelta,
        chan_filter=None,
    ):
        self.path = Path(path)
        self.chan_catalog = chan_catalog
        self.chunk_offset = chunk_offset
        self.chan_filter = chan_filter
        self.paths = {}
        self.channels = {}
        self._skipped_files: List[Tuple[str, str]] = []

        for fpath in sorted(self.path.iterdir()):
            if fpath.is_dir():
                continue
            match = FILENAME_RE.match(fpath.name)
            if not match:
                if fpath.suffix.lower() in (".ms", ".mseed", ".miniseed"):
                    self._skipped_files.append(
                        (fpath.name, "filename does not match NETSTA_CHALOC_YYYYDOY pattern")
                    )
                continue

            info = match.groupdict()
            ts_start = yj_to_datetime(info["yj"]) + self.chunk_offset
            ts_end = ts_start + timedelta(days=1)
            timespan = DateTimeRange(ts_start, ts_end)

            station = Station(info["net"], info["sta"], location=info["loc"])
            channel = Channel(ChannelType(info["cha"], info["loc"]), station)

            if self.chan_filter is not None and not self.chan_filter(channel):
                continue

            key = str(timespan)
            self.channels.setdefault(key, []).append(channel)
            self.paths[(key, station.network, station.name, channel.type.name)] = str(fpath)

        if self._skipped_files:
            LOG.warning(
                "Skipped %d file(s) with unrecognized names:", len(self._skipped_files)
            )
            for name, reason in self._skipped_files[:10]:
                LOG.warning("  %s — %s", name, reason)

    def get_timespans(self) -> List[DateTimeRange]:
        return sorted(
            [DateTimeRange.from_range_text(key) for key in self.channels.keys()],
            key=lambda ts: ts.start_datetime,
        )

    def get_channels(self, timespan: DateTimeRange) -> List[Channel]:
        key = str(timespan)
        channels = self.channels.get(key, [])
        return [self.chan_catalog.get_full_channel(timespan, chan) for chan in channels]

    def read_data(self, timespan: DateTimeRange, chan: Channel) -> ChannelData:
        key = (str(timespan), chan.station.network, chan.station.name, chan.type.name)
        fpath = self.paths.get(key)
        if fpath is None:
            return ChannelData(obspy.Stream())
        return ChannelData(obspy.read(fpath))

    def get_inventory(self, timespan: DateTimeRange, station: Station):
        return self.chan_catalog.get_inventory(timespan, station)


# ---------------------------------------------------------------------------
# Thread-safe ASDF store
# ---------------------------------------------------------------------------


class LockedASDFCCStore(ASDFCCStore):
    def __init__(self, directory: str, mode: str = "a"):
        super().__init__(directory, mode=mode)
        self._io_lock = threading.RLock()

    def contains(self, src: Station, rec: Station, timespan: DateTimeRange) -> bool:
        with self._io_lock:
            return super().contains(src, rec, timespan)

    def append(self, timespan: DateTimeRange, src: Station, rec: Station, ccs) -> None:
        with self._io_lock:
            super().append(timespan, src, rec, ccs)

    def get_station_pairs(self):
        with self._io_lock:
            return super().get_station_pairs()

    def get_timespans(self, src: Station, rec: Station):
        with self._io_lock:
            return super().get_timespans(src, rec)

    def read(self, timespan: DateTimeRange, src_sta: Station, rec_sta: Station):
        with self._io_lock:
            return super().read(timespan, src_sta, rec_sta)


# ---------------------------------------------------------------------------
# Stacking (local fixed version)
# ---------------------------------------------------------------------------


def stack_pair_fixed(src_sta, rec_sta, timespans, cc_store, fft_params) -> List[Stack]:
    fauto = int(src_sta == rec_sta)
    nccomp = fft_params.ncomp * fft_params.ncomp
    if fft_params.rotation and fft_params.correction:
        if not fft_params.correction_csv:
            LOG.warning("Missing correction_csv but rotation=True and correction=True")
            locs = []
        else:
            locs = pd.read_csv(fft_params.correction_csv)
    else:
        locs = []

    if fft_params.ncomp == 1:
        enz_system = ["ZZ"]
    else:
        enz_system = ["EE", "EN", "EZ", "NE", "NN", "NZ", "ZE", "ZN", "ZZ"]

    rtz_components = ["ZR", "ZT", "ZZ_R", "RR", "RT", "RZ", "TR", "TT", "TZ"]
    num_chunk = len(timespans) * nccomp
    num_segmts, npts_segmt = calc_segments(fft_params, num_chunk)
    cc_array = np.zeros((num_chunk * num_segmts, npts_segmt), dtype=np.float32)
    cc_time = np.zeros(num_chunk * num_segmts, dtype=np.float32)
    cc_ngood = np.zeros(num_chunk * num_segmts, dtype=np.int16)
    cc_comp = np.chararray(num_chunk * num_segmts, itemsize=2, unicode=True)

    iseg = 0
    tparameters = {}
    for ts in timespans:
        if ts.end_datetime > fft_params.end_date or ts.start_datetime < fft_params.start_date:
            continue
        cross_correlations = cc_store.read(ts, src_sta, rec_sta)
        if not validate_pairs(
            fft_params.ncomp, str((src_sta, rec_sta)), fauto, ts, len(cross_correlations)
        ):
            continue
        for cc in cross_correlations:
            src_chan, rec_chan = cc.src, cc.rec
            tparameters, tdata = cc.parameters, cc.data
            tcmp1 = src_chan.get_orientation()
            tcmp2 = rec_chan.get_orientation()
            if fft_params.substack:
                for ii in range(tdata.shape[0]):
                    cc_array[iseg] = tdata[ii]
                    cc_time[iseg] = tparameters["time"][ii]
                    cc_ngood[iseg] = tparameters["ngood"][ii]
                    cc_comp[iseg] = tcmp1 + tcmp2
                    iseg += 1
            else:
                cc_array[iseg] = tdata
                cc_time[iseg] = tparameters["time"]
                cc_ngood[iseg] = tparameters["ngood"]
                cc_comp[iseg] = tcmp1 + tcmp2
                iseg += 1

    stack_results: List[Stack] = []
    if iseg <= 1 and fauto == 1:
        return stack_results

    if fft_params.rotation:
        bigstack = np.zeros(shape=(9, npts_segmt), dtype=np.float32)
        if fft_params.stack_method == StackMethod.ALL:
            bigstack1 = np.zeros(shape=(9, npts_segmt), dtype=np.float32)
            bigstack2 = np.zeros(shape=(9, npts_segmt), dtype=np.float32)

    def append_stacks(comp, params, stack_data):
        for method, data in stack_data:
            name = f"Allstack_{method.value}" if isinstance(method, StackMethod) else str(method)
            stack_results.append(Stack(comp, name, dict(params), data))

    iflag = 1
    for icomp in range(nccomp):
        comp = enz_system[icomp]
        indx = np.where(cc_comp.lower() == comp.lower())[0]
        if len(indx) < 2:
            iflag = 0
            continue

        (
            cc_final, ngood_final, stamps_final,
            allstacks1, allstacks2, allstacks3, nstacks,
        ) = noise_module.stacking(cc_array[indx], cc_time[indx], cc_ngood[indx], fft_params)
        if not len(allstacks1):
            continue

        params = dict(tparameters)
        params["time"] = stamps_final[0]
        params["ngood"] = int(nstacks)

        if fft_params.rotation:
            bigstack[icomp] = allstacks1
            if fft_params.stack_method == StackMethod.ALL:
                bigstack1[icomp] = allstacks2
                bigstack2[icomp] = allstacks3
        else:
            if fft_params.stack_method != StackMethod.ALL:
                to_write = [(fft_params.stack_method, allstacks1)]
            else:
                to_write = [
                    (StackMethod.LINEAR, allstacks1),
                    (StackMethod.PWS, allstacks2),
                    (StackMethod.ROBUST, allstacks3),
                ]
            append_stacks(comp, params, to_write)

        if fft_params.keep_substack:
            for ii in range(cc_final.shape[0]):
                sub_params = dict(tparameters)
                sub_params["time"] = stamps_final[ii]
                sub_params["ngood"] = int(ngood_final[ii])
                stack_name = "T" + str(int(stamps_final[ii]))
                append_stacks(comp, sub_params, [(stack_name, cc_final[ii])])

    if fft_params.rotation and iflag:
        if np.all(bigstack == 0):
            return stack_results
        params = dict(tparameters)
        params["station_source"] = src_sta.name
        params["station_receiver"] = rec_sta.name

        if fft_params.stack_method != StackMethod.ALL:
            bigstack_rotated = noise_module.rotation(bigstack, params, locs)
            for icomp in range(nccomp):
                rot_params = dict(params)
                rot_params["time"] = stamps_final[0]
                rot_params["ngood"] = int(nstacks)
                append_stacks(
                    rtz_components[icomp], rot_params,
                    [(fft_params.stack_method, bigstack_rotated[icomp])],
                )
        else:
            bigstack_rotated = noise_module.rotation(bigstack, params, locs)
            bigstack_rotated1 = noise_module.rotation(bigstack1, params, locs)
            bigstack_rotated2 = noise_module.rotation(bigstack2, params, locs)
            for icomp in range(nccomp):
                rot_params = dict(params)
                rot_params["time"] = stamps_final[0]
                rot_params["ngood"] = int(nstacks)
                stacks = [
                    (StackMethod.LINEAR, bigstack_rotated[icomp]),
                    (StackMethod.PWS, bigstack_rotated1[icomp]),
                    (StackMethod.ROBUST, bigstack_rotated2[icomp]),
                ]
                append_stacks(rtz_components[icomp], rot_params, stacks)

    return stack_results


# ---------------------------------------------------------------------------
# Pipeline stages
# ---------------------------------------------------------------------------


def discover_dataset(paths: DatasetPaths) -> Tuple:
    ms_files = sorted(
        f for f in paths.data_dir.iterdir()
        if f.is_file() and FILENAME_RE.match(f.name)
    )
    if not ms_files:
        raise FileNotFoundError(
            f"No miniSEED files matching NETSTA_CHALOC_YYYYDOY.{{ms,mseed,miniseed}} "
            f"found in {paths.data_dir}"
        )
    if not paths.station_csv.exists():
        raise FileNotFoundError(f"Missing station.csv: {paths.station_csv}")

    parsed = []
    for path in ms_files:
        match = FILENAME_RE.match(path.name)
        if match:
            parsed.append(match.groupdict())

    networks = sorted({item["net"] for item in parsed})
    stations = sorted({item["sta"] for item in parsed})
    channels = sorted({item["cha"] for item in parsed})
    locations = sorted({item["loc"] for item in parsed})
    days = sorted({item["yj"] for item in parsed})

    chunk_offset = infer_chunk_offset(ms_files)
    start = yj_to_datetime(days[0]) + chunk_offset
    end = yj_to_datetime(days[-1]) + timedelta(days=1) + chunk_offset
    timerange = DateTimeRange(start, end)

    LOG.info("NoisePy version : %s", __version__)
    LOG.info("DATA_DIR        : %s", paths.data_dir)
    LOG.info("STATION_CSV     : %s", paths.station_csv)
    LOG.info("CCF_DIR         : %s", paths.ccf_dir)
    LOG.info("STACK_DIR       : %s", paths.stack_dir)
    LOG.info("miniSEED files  : %d", len(ms_files))
    LOG.info("networks        : %s", networks)
    LOG.info("stations        : %d stations — %s", len(stations), stations)
    LOG.info("channels        : %s", channels)
    LOG.info("locations       : %s", locations)
    LOG.info("chunk offset    : %s", chunk_offset)
    LOG.info("time range      : %s", timerange)

    return ms_files, networks, channels, start, end, chunk_offset


def build_noisepy_config(
    networks, channels, start, end, options: ProcessingOptions
) -> ConfigParameters:
    config = ConfigParameters()
    config.start_date = start
    config.end_date = end
    config.networks = networks
    config.stations = [options.station_pattern]
    config.channels = channels

    config.acorr_only = options.acorr_only
    config.xcorr_only = options.xcorr_only

    config.inc_hours = options.inc_hours
    config.cc_len = options.cc_len
    config.step = options.step
    config.maxlag = options.maxlag

    config.ncomp = options.ncomp
    config.rotation = options.rotation
    config.stationxml = False
    config.rm_resp = RmResp(options.rm_resp)

    config.sampling_rate = options.sampling_rate
    config.freqmin = options.freqmin
    config.freqmax = options.freqmax
    config.max_over_std = options.max_over_std

    config.freq_norm = FreqNorm(options.freq_norm)
    config.smoothspect_N = options.smoothspect_n
    config.time_norm = TimeNorm(options.time_norm)
    config.smooth_N = options.smooth_n

    config.cc_method = CCMethod.XCORR
    config.substack = options.substack
    config.substack_windows = options.substack_windows
    config.keep_substack = options.keep_substack
    config.stack_method = StackMethod(options.stack_method)

    LOG.info("NoisePy config: %s", config)
    return config


def run_cc_stage(paths: DatasetPaths, options: ProcessingOptions) -> None:
    ms_files, networks, channels, start, end, chunk_offset = discover_dataset(paths)
    config = build_noisepy_config(networks, channels, start, end, options)

    paths.ccf_dir.mkdir(parents=True, exist_ok=True)

    catalog = LocalCSVChannelCatalog(str(paths.station_csv))
    raw_store = FlatMiniSeedStore(
        paths.data_dir,
        catalog,
        chunk_offset,
        channel_filter(config.networks, config.stations, config.channels),
    )

    timespans = raw_store.get_timespans()
    LOG.info("timespans       : %d", len(timespans))
    if timespans:
        first_channels = raw_store.get_channels(timespans[0])
        LOG.info("channels in first timespan: %d", len(first_channels))

    cc_store = LockedASDFCCStore(str(paths.ccf_dir))

    LOG.info("Starting cross-correlation...")
    cross_correlate(
        raw_store,
        config,
        cc_store,
        pair_filter=lambda src, rec: is_cross_station_pair(src.station, rec.station),
    )
    LOG.info("Cross-correlation finished.")


def run_stack_stage(paths: DatasetPaths, options: ProcessingOptions) -> None:
    if not paths.ccf_dir.exists():
        raise FileNotFoundError(f"CCF directory not found: {paths.ccf_dir}")

    ccf_files = sorted(paths.ccf_dir.glob("*.h5"))
    if not ccf_files:
        raise FileNotFoundError(f"No .h5 files found in {paths.ccf_dir}")

    cc_store = LockedASDFCCStore(str(paths.ccf_dir), mode="r")
    pairs = cc_store.get_station_pairs()
    if not pairs:
        raise FileNotFoundError(f"No station pairs found in CCF store: {paths.ccf_dir}")

    ccf_timespans = []
    for ccf_file in ccf_files:
        try:
            start_text, end_text = ccf_file.stem.split("T", maxsplit=1)
            span = DateTimeRange(
                datetime.strptime(start_text, "%Y_%m_%d_%H_%M_%S").replace(tzinfo=timezone.utc),
                datetime.strptime(end_text, "%Y_%m_%d_%H_%M_%S").replace(tzinfo=timezone.utc),
            )
            ccf_timespans.append(span)
        except ValueError:
            LOG.warning("Skipping CCF file with unexpected name: %s", ccf_file.name)

    if not ccf_timespans:
        raise FileNotFoundError(f"No valid timespans in CCF store: {paths.ccf_dir}")

    start = min(span.start_datetime for span in ccf_timespans)
    end = max(span.end_datetime for span in ccf_timespans)

    sample_pair = pairs[0]
    sample_span = cc_store.get_timespans(*sample_pair)[0]
    sample_cc = cc_store.read(sample_span, *sample_pair)[0]
    src_channel = str(sample_cc.src).split("_")[0]
    rec_channel = str(sample_cc.rec).split("_")[0]
    channels = sorted({src_channel, rec_channel})
    networks = sorted({sta.network for pair in pairs for sta in pair})

    config = build_noisepy_config(networks, channels, start, end, options)
    paths.stack_dir.mkdir(parents=True, exist_ok=True)

    stack_store = ASDFStackStore(str(paths.stack_dir))

    LOG.info("CCF files       : %d", len(ccf_files))
    LOG.info("station pairs   : %d", len(pairs))
    LOG.info("Starting stacking...")

    stack_cross_correlations(cc_store, stack_store, config)
    LOG.info("Stacking finished.")


def stack_cross_correlations(cc_store, stack_store, fft_params) -> None:
    stations = set(fft_params.stations)
    networks = set(fft_params.networks)
    stack_timespan = DateTimeRange(fft_params.start_date, fft_params.end_date)

    if fft_params.stack_method == StackMethod.ALL:
        requested_names = {"Allstack_linear", "Allstack_pws", "Allstack_robust"}
    else:
        requested_names = {f"Allstack_{fft_params.stack_method.value}"}

    def sta_filter(sta: Station) -> bool:
        return ("*" in stations or sta.name in stations) and (
            "*" in networks or sta.network in networks
        )

    def pair_has_requested_stacks(src_sta: Station, rec_sta: Station) -> bool:
        if not is_cross_station_pair(src_sta, rec_sta):
            return True
        if not stack_store.contains(src_sta, rec_sta, stack_timespan):
            return False
        existing = {stack.name for stack in stack_store.read(stack_timespan, src_sta, rec_sta)}
        return requested_names.issubset(existing)

    pairs = sorted(
        filter(
            lambda p: sta_filter(p[0]) and sta_filter(p[1]) and is_cross_station_pair(p[0], p[1]),
            cc_store.get_station_pairs(),
        ),
        key=lambda x: str(x),
    )
    if not pairs:
        raise IOError(
            f"No station pairs for networks={fft_params.networks}, stations={fft_params.stations}"
        )

    failed = []
    LOG.info("stack pairs     : %d", len(pairs))
    for idx, (src_sta, rec_sta) in enumerate(pairs, start=1):
        LOG.info("[%d/%d] stacking %s-%s", idx, len(pairs), src_sta, rec_sta)
        try:
            if pair_has_requested_stacks(src_sta, rec_sta):
                LOG.debug("  already stacked, skipping")
                continue
            timespans = cc_store.get_timespans(src_sta, rec_sta)
            stacks = stack_pair_fixed(src_sta, rec_sta, timespans, cc_store, fft_params)
            if not stacks:
                LOG.warning("No stacks for %s_%s", src_sta, rec_sta)
                failed.append((src_sta, rec_sta))
                continue
            stack_store.append(stack_timespan, src_sta, rec_sta, stacks)
        except Exception as exc:
            LOG.error("Error stacking pair %s-%s: %s", src_sta, rec_sta, exc)
            failed.append((src_sta, rec_sta))

    if failed:
        LOG.error(
            "Failed %d/%d pairs: %s",
            len(failed), len(pairs),
            ", ".join(f"{s}-{r}" for s, r in failed),
        )


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_options(options: ProcessingOptions) -> List[str]:
    errors = []
    if options.freqmin >= options.freqmax:
        errors.append(f"freqmin ({options.freqmin}) must be < freqmax ({options.freqmax})")
    if options.maxlag <= 0:
        errors.append(f"maxlag must be > 0, got {options.maxlag}")
    if options.cc_len <= 0:
        errors.append(f"cc_len must be > 0, got {options.cc_len}")
    if options.step <= 0:
        errors.append(f"step must be > 0, got {options.step}")
    if options.sampling_rate <= 0:
        errors.append(f"sampling_rate must be > 0, got {options.sampling_rate}")
    if options.inc_hours <= 0:
        errors.append(f"inc_hours must be > 0, got {options.inc_hours}")
    if options.ncomp not in (1, 3):
        errors.append(f"ncomp must be 1 or 3, got {options.ncomp}")
    if options.ncomp > 1 and not options.rotation:
        errors.append("ncomp > 1 requires rotation=True (multi-component not fully tested without rotation)")
    if options.acorr_only and options.xcorr_only:
        errors.append("acorr_only and xcorr_only cannot both be True")
    if options.inc_hours != 24:
        errors.append(
            f"inc_hours={options.inc_hours} but FlatMiniSeedStore uses daily files. "
            "Only inc_hours=24 is supported with this data format."
        )
    return errors


def validate_paths(paths: DatasetPaths, stage: str) -> List[str]:
    errors = []
    if stage in ("cc", "all"):
        if not paths.data_dir.is_dir():
            errors.append(f"data_dir does not exist: {paths.data_dir}")
        if not paths.station_csv.is_file():
            errors.append(f"station_csv does not exist: {paths.station_csv}")
    if stage == "stack":
        if not paths.ccf_dir.is_dir():
            errors.append(f"ccf_dir does not exist: {paths.ccf_dir}")
    return errors


def run_validate(paths: DatasetPaths, options: ProcessingOptions, stage: str) -> bool:
    LOG.info("=== Validation Mode ===")
    all_errors = []

    all_errors.extend(validate_options(options))
    all_errors.extend(validate_paths(paths, stage))

    if not all_errors and stage in ("cc", "all") and paths.station_csv.is_file():
        try:
            normalize_station_csv(
                pd.read_csv(
                    str(paths.station_csv),
                    dtype={"network": str, "station": str, "channel": str},
                )
            )
            LOG.info("station.csv     : OK")
        except Exception as e:
            all_errors.append(f"station.csv error: {e}")

    if not all_errors and stage in ("cc", "all") and paths.data_dir.is_dir():
        ms_count = sum(
            1 for f in paths.data_dir.iterdir()
            if f.is_file() and FILENAME_RE.match(f.name)
        )
        LOG.info("miniSEED files  : %d", ms_count)
        if ms_count == 0:
            all_errors.append(f"No valid miniSEED files in {paths.data_dir}")

    if all_errors:
        LOG.error("Validation FAILED:")
        for err in all_errors:
            LOG.error("  - %s", err)
        return False

    LOG.info("Validation PASSED — all checks OK")
    LOG.info("Options: %s", json.dumps(asdict(options), indent=2, default=str))
    return True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    parser.add_argument(
        "--config", type=Path, default=None,
        help="YAML config file. CLI flags override YAML values.",
    )
    parser.add_argument(
        "--stage", choices=["cc", "stack", "all"], default="all",
        help="Pipeline stage to run (default: all).",
    )
    parser.add_argument(
        "--validate-only", action="store_true",
        help="Check inputs and parameters without running the pipeline.",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Remove existing CCF/stack output directories before running.",
    )
    parser.add_argument(
        "--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO",
        help="Logging verbosity (default: INFO).",
    )

    paths_group = parser.add_argument_group("paths")
    paths_group.add_argument("--data-dir", type=Path, help="Input miniSEED directory.")
    paths_group.add_argument("--station-csv", type=Path, help="Station CSV file path.")
    paths_group.add_argument("--ccf-dir", type=Path, help="Cross-correlation output directory.")
    paths_group.add_argument("--stack-dir", type=Path, help="Stack output directory.")

    proc_group = parser.add_argument_group("processing")
    proc_group.add_argument("--freqmin", type=float, default=1.0, help="Bandpass low corner Hz (default: 1.0).")
    proc_group.add_argument("--freqmax", type=float, default=5.0, help="Bandpass high corner Hz (default: 5.0).")
    proc_group.add_argument("--sampling-rate", type=float, default=20.0, help="Target sampling rate Hz (default: 20.0).")
    proc_group.add_argument("--maxlag", type=int, default=100, help="CC lag seconds (default: 100).")
    proc_group.add_argument("--inc-hours", type=int, default=24, help="Processing chunk hours (default: 24).")
    proc_group.add_argument("--cc-len", type=int, default=3600, help="Correlation window seconds (default: 3600).")
    proc_group.add_argument("--step", type=float, default=1800.0, help="Sliding step seconds (default: 1800).")
    proc_group.add_argument("--max-over-std", type=int, default=10, help="Max over std threshold (default: 10).")
    proc_group.add_argument("--smoothspect-n", type=int, default=10, help="Spectral smoothing half-width (default: 10).")
    proc_group.add_argument("--smooth-n", type=int, default=10, help="Time-domain smoothing half-width (default: 10).")
    proc_group.add_argument("--substack", action="store_true", default=True, help="Enable substack (default: True).")
    proc_group.add_argument("--no-substack", action="store_true", help="Disable substack.")
    proc_group.add_argument("--substack-windows", type=int, default=1, help="Substack windows (default: 1).")
    proc_group.add_argument("--keep-substack", action="store_true", help="Keep per-window substacks.")
    proc_group.add_argument("--station-pattern", default="*", help="Station filter pattern (default: *).")
    proc_group.add_argument(
        "--stack-method", choices=[m.value for m in StackMethod], default=StackMethod.ALL.value,
        help="Stack method (default: all).",
    )
    proc_group.add_argument(
        "--time-norm", choices=[m.value for m in TimeNorm], default=TimeNorm.ONE_BIT.value,
        help="Time normalization (default: one_bit).",
    )
    proc_group.add_argument(
        "--freq-norm", choices=[m.value for m in FreqNorm], default=FreqNorm.RMA.value,
        help="Frequency normalization (default: rma).",
    )
    proc_group.add_argument(
        "--rm-resp", choices=[m.value for m in RmResp], default=RmResp.NO.value,
        help="Response removal mode (default: no).",
    )
    proc_group.add_argument("--ncomp", type=int, default=1, help="Number of components: 1 or 3 (default: 1).")
    proc_group.add_argument("--rotation", action="store_true", help="Enable rotation (requires ncomp=3).")
    proc_group.add_argument("--xcorr-only", action="store_true", default=True, help="Cross-correlation only (default: True).")
    proc_group.add_argument("--acorr-only", action="store_true", help="Auto-correlation only.")

    return parser


def normalize_config_key(key: str) -> str:
    return key.strip().replace("-", "_")


def load_yaml_config(config_path: Path) -> Dict[str, Any]:
    resolved = config_path.expanduser().resolve()
    with resolved.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}

    if not isinstance(data, dict):
        raise ValueError(f"Config must be a YAML mapping: {resolved}")

    normalized = {normalize_config_key(k): v for k, v in data.items()}

    for path_key in ("data_dir", "station_csv", "ccf_dir", "stack_dir"):
        value = normalized.get(path_key)
        if value is None:
            continue
        path_value = Path(value).expanduser()
        if not path_value.is_absolute():
            path_value = resolved.parent / path_value
        normalized[path_key] = path_value.resolve()

    return normalized


def parse_cli() -> argparse.Namespace:
    parser = build_parser()
    preliminary, _ = parser.parse_known_args()

    if preliminary.config is not None:
        yaml_cfg = load_yaml_config(preliminary.config)
        parser.set_defaults(**yaml_cfg)

    return parser.parse_args()


def build_paths(args: argparse.Namespace) -> DatasetPaths:
    data_dir = args.data_dir
    station_csv = args.station_csv
    ccf_dir = args.ccf_dir
    stack_dir = args.stack_dir

    if data_dir is None:
        raise SystemExit(
            "Error: --data-dir is required (or set data_dir in YAML config).\n"
            "Example: python local_cc_workflow.py --data-dir ./my_data --station-csv ./station.csv"
        )

    data_dir = Path(data_dir)
    if station_csv is None:
        station_csv = data_dir / "station.csv"
    else:
        station_csv = Path(station_csv)

    if ccf_dir is None:
        ccf_dir = data_dir.parent / f"{data_dir.name}_ccf"
    else:
        ccf_dir = Path(ccf_dir)

    if stack_dir is None:
        stack_dir = data_dir.parent / f"{data_dir.name}_stack"
    else:
        stack_dir = Path(stack_dir)

    return DatasetPaths(
        data_dir=data_dir,
        station_csv=station_csv,
        ccf_dir=ccf_dir,
        stack_dir=stack_dir,
    )


def build_options(args: argparse.Namespace) -> ProcessingOptions:
    substack = not args.no_substack if args.no_substack else args.substack
    return ProcessingOptions(
        inc_hours=args.inc_hours,
        cc_len=args.cc_len,
        step=args.step,
        maxlag=args.maxlag,
        sampling_rate=args.sampling_rate,
        freqmin=args.freqmin,
        freqmax=args.freqmax,
        max_over_std=args.max_over_std,
        smoothspect_n=args.smoothspect_n,
        smooth_n=args.smooth_n,
        substack=substack,
        substack_windows=args.substack_windows,
        keep_substack=args.keep_substack,
        station_pattern=args.station_pattern,
        stack_method=args.stack_method,
        time_norm=args.time_norm,
        freq_norm=args.freq_norm,
        rm_resp=args.rm_resp,
        ncomp=args.ncomp,
        rotation=args.rotation,
        xcorr_only=args.xcorr_only and not args.acorr_only,
        acorr_only=args.acorr_only,
    )


def main():
    args = parse_cli()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    paths = build_paths(args)
    options = build_options(args)

    if args.validate_only:
        ok = run_validate(paths, options, args.stage)
        raise SystemExit(0 if ok else 1)

    errors = validate_options(options)
    if errors:
        for err in errors:
            LOG.error("Config error: %s", err)
        raise SystemExit(1)

    if args.overwrite:
        import shutil
        if args.stage in ("cc", "all") and paths.ccf_dir.exists():
            LOG.info("Overwrite: removing %s", paths.ccf_dir)
            shutil.rmtree(paths.ccf_dir)
        if args.stage in ("stack", "all") and paths.stack_dir.exists():
            LOG.info("Overwrite: removing %s", paths.stack_dir)
            shutil.rmtree(paths.stack_dir)

    if args.stage in ("cc", "all"):
        run_cc_stage(paths, options)

    if args.stage in ("stack", "all"):
        run_stack_stage(paths, options)

    LOG.info("Done.")


if __name__ == "__main__":
    main()
