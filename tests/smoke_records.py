"""Raw records for the continuous-integration smoke test: does the chain run from raw files to the inversion project?

    python tests/smoke_records.py make --out DIR     16 stations, 6 hours at 25 Hz, about a minute
    python tests/smoke_records.py run --dir DIR      lant run all to the inversion project, gates confirmed

Not a test of accuracy or resolution (lant's synthetic tests are at the inversion level, docs/validation.md): plane
Rayleigh waves from random azimuths in a laterally uniform 1-D medium, phase velocities from lant's 1-D forward,
5 % incoherent noise and one burst.  The records are deliberately mixed (miniSEED one folder per station, SAC without
network code in day folders) and the station table has no header and longitude first, as lant must read such input.
The gates are confirmed only because these records are generated here; never do this with real data.
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

Z = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.9, 2.4, 3.0]
VS = [2.30, 2.45, 2.62, 2.78, 2.90, 3.00, 3.08, 3.16, 3.25, 3.33, 3.40]
SIZE_KM, MIN_SPACING_KM, BAND_HZ = 9.0, 0.4, (0.3, 12.0)
LON0, LAT0 = 100.0, 30.0
YAML = """# Smoke test of the chain (tests/smoke_records.py): 16 stations, 6 hours; not for real data.
project: {name: smoke, output_dir: output}
data: {raw_dir: data, station_file: stations.txt, velocity_range_kms: [1.5, 3.5]}
inversion: {median_min_pairs: 20, iterations: 2}
checks: {checkerboards: [], vertical: [], controls: []}
"""
REASON = "Generated smoke-test records (tests/smoke_records.py)"


def make(out, stations=16, hours=6, fs=25.0, waves=100):
    import obspy

    from lant.inversion.start import forward
    from lant.signal import km_per_degree

    rng = np.random.default_rng(0)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    xy = []
    while len(xy) < stations:
        q = rng.uniform(0, SIZE_KM, 2)
        if all(np.hypot(*(q - r)) > MIN_SPACING_KM for r in xy):
            xy.append(q)
    xy = np.array(xy)
    kx, ky = km_per_degree(LAT0)
    lon, lat = LON0 + xy[:, 0] / kx, LAT0 + xy[:, 1] / ky
    names = ["S%02d" % (i + 1) for i in range(stations)]
    with open(out / "stations.txt", "w") as f:
        for n, lo, la in zip(names, lon, lat):
            f.write("%s  %.6f  %.6f  500\n" % (n, lo, la))
    x, y = (lon - LON0) * kx, (lat - LAT0) * ky
    freq = np.geomspace(*BAND_HZ, 60)
    c = forward(np.array(Z), np.array(VS), 1 / freq)
    nh = int(3600 * fs)
    f = np.fft.rfftfreq(nh, 1 / fs)
    band = (f >= BAND_HZ[0]) & (f <= BAND_HZ[1])
    cf = np.ones_like(f)
    cf[band] = np.interp(np.log(f[band]), np.log(freq), c)
    record = np.zeros((stations, hours * nh), np.float32)
    for h in range(hours):
        az = rng.uniform(0, 2 * np.pi, waves)
        spec = (rng.standard_normal((waves, f.size)) + 1j * rng.standard_normal((waves, f.size))) * band
        for s in range(stations):
            delay = (x[s] * np.sin(az) + y[s] * np.cos(az))[:, None] / cf[None, :]
            record[s, h * nh : (h + 1) * nh] = np.fft.irfft((spec * np.exp(-2j * np.pi * f * delay)).sum(0), n=nh)
    record += 0.05 * record.std() * rng.standard_normal(record.shape).astype(np.float32)
    record[3, 30000 : 30000 + int(600 * fs)] *= 40.0  # a burst: segment rejection
    t0 = obspy.UTCDateTime(2024, 6, 1)
    for s, n in enumerate(names):
        tr = obspy.Trace(record[s] * 1e4, header=dict(network="SY", station=n, channel="HHZ", sampling_rate=fs,
                                                       starttime=t0))  # fmt: skip
        if s % 2 == 0:
            d = out / "data" / "mseed" / n
            d.mkdir(parents=True, exist_ok=True)
            tr.write(str(d / ("%s_day1.mseed" % n)), format="MSEED", encoding="FLOAT32")
        else:
            d = out / "data" / "sac_files" / "day01"
            d.mkdir(parents=True, exist_ok=True)
            tr.stats.network = ""
            tr.write(str(d / ("rec_%s.sac" % n)), format="SAC")
    (out / "synthetic.yaml").write_text(YAML)
    print("wrote %s (stations.txt, data/, synthetic.yaml)" % out)


def run(directory):
    """lant run all to the inversion project, confirming each gate (the records are generated here)."""
    config = str((Path(directory) / "synthetic.yaml").resolve())
    base = [sys.executable, "-m", "lant"]
    for _ in range(6):
        code = subprocess.run(base + ["run", "all", "-c", config, "--to", "inversion"]).returncode
        if code == 0:
            return
        if code != 3:
            raise SystemExit(code)
        state = json.loads(subprocess.check_output(base + ["status", "-c", config, "--json"]))
        gate = next((s for s, v in state["steps"].items() if v["state"] == "to confirm"), None)
        if gate is None:
            raise SystemExit("stopped without a gate")
        subprocess.run(base + ["confirm", gate, "-c", config, "--reason", REASON], check=True)
    raise SystemExit("too many stops")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("make").add_argument("--out", required=True)
    sub.add_parser("run").add_argument("--dir", required=True)
    a = p.parse_args(argv)
    if a.cmd == "make":
        if os.path.exists(a.out) and any(Path(a.out).iterdir()):
            raise SystemExit("%s is not empty" % a.out)
        make(a.out)
    else:
        run(a.dir)


if __name__ == "__main__":
    main()
