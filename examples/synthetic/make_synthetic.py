#!/usr/bin/env python
"""Synthetic test data with a known answer: ambient noise of a laterally uniform medium.

The noise field is a sum of plane Rayleigh waves arriving from random azimuths (an isotropic diffuse field), each
with a random spectrum, propagating with the phase velocity c(f) of a known 1-D Vs model.  For such a field the
cross-correlation of two stations is the 2-D Green's function, so the measured phase velocities should follow c(f)
and the inversion should return the 1-D model (within its resolution).

    python examples/synthetic/make_synthetic.py [--out examples/synthetic] [--stations 16] [--days 3]
    python examples/synthetic/make_synthetic.py --model sediment --out /some/dir     # a slow basin, 2 km array
    python examples/synthetic/make_synthetic.py --layout line --stations 24 --out /some/dir   # a line array
    python examples/synthetic/make_synthetic.py --quick --responses --out /some/dir   # two sensor types + responses
    python examples/synthetic/make_synthetic.py --quick --directional 0.7 --out /some/dir   # directional noise

writes
    data/                       the "raw" records, deliberately mixed: miniSEED and SAC (SAC without network code),
                                in nested folders with arbitrary names
    stations.txt                the station table (no header, longitude before latitude: recognised automatically)
    true_model.txt, true_dispersion.txt   the answer (depth Vs; frequency c)
    synthetic.yaml              a ready configuration (output in examples/synthetic/output)
    responses.xml               with --responses: the StationXML responses of the two sensor types (alternate
                                stations: geophones of 2 Hz and 10 Hz natural frequency, damping 0.7), named in
                                synthetic.yaml as data.response; synthetic_noresp.yaml is the same without it
Needs fwd1d (bash tools/install_external.sh dsurftomo) to compute c(f) of the model.  About 5-20 minutes, 300 MB.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

MODELS = {   # depth (km), Vs (km/s); array size (km), minimum spacing (km), noise band (Hz)
    "rock": dict(z=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.9, 2.4, 3.0],
                 vs=[2.30, 2.45, 2.62, 2.78, 2.90, 3.00, 3.08, 3.16, 3.25, 3.33, 3.40], size=9.0, spacing=0.4, band=(0.3, 12.0),
                 vrange=(1.5, 3.5)),
    "sediment": dict(z=[0.0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2.0], vs=[0.25, 0.30, 0.40, 0.55, 0.80, 1.00, 1.40, 1.80, 2.10, 2.50],
                     size=2.0, spacing=0.08, band=(1.0, 15.0), vrange=(0.2, 1.0)),
}
LON0, LAT0 = 100.0, 30.0          # an arbitrary place (no real array)


YAML = """# Configuration of the synthetic example (examples/synthetic/make_synthetic.py).
# Only the three paths are required; the velocity range is given to keep the example fast.  The other lines
# (a) lower thresholds meant for >= 50 stations, because this example has only 16 stations (120 pairs) and 3 days,
# and (b) make the inversion part fast enough for a first test (6 iterations, one checkerboard and one control).
project: {name: %(name)s, output_dir: output}
data: {raw_dir: data, station_file: stations.txt, velocity_range_kms: [%(vlo)s, %(vhi)s]}
dispersion: {njack: 3, nproc: 4}
%(export)sinversion: {median_min_pairs: 20, iterations: 6, threads: 4, njobs: 3}
checks: {checkerboards: [{name: cb, cell_km: [%(cell)s, %(cell)s]}], vertical: [], controls: [start_grad0.05]}
"""


YAML_LINE = """# Configuration of the synthetic line (examples/synthetic/make_synthetic.py --layout line).
# Only the three paths are required; the array type is found from the station table.  The other lines keep the
# example fast.
project: {name: synthetic_line, output_dir: output}
data: {raw_dir: data, station_file: stations.txt, velocity_range_kms: [%(vlo)s, %(vhi)s]}
dispersion: {nproc: 4}
inversion: {iterations: 6, threads: 4, njobs: 3}
"""

YAML_QUICK = """# Configuration of the quick first test (examples/synthetic/make_synthetic.py --quick): 10 stations, 6 hours.
# Only the three paths are required.  The other lines keep the test short: a velocity range that fits the model,
# a coarse grid, few iterations and no check inversions (DSurfTomo takes most of the time, and its time grows
# with the number of grid nodes and periods).  Do not copy them to a real project.
project: {name: quick, output_dir: output}
data: {raw_dir: data, station_file: stations.txt, velocity_range_kms: [%(vlo)s, %(vhi)s]}
correlation: {fs_max: 25}
dispersion: {njack: 3, nproc: 4}
export: {tmin_s: 0.2}
inversion: {median_min_pairs: 10, grid_refine: 1, iterations: 3, threads: 4}
checks: {checkerboards: [], vertical: [], controls: []}
"""


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default=HERE)
    p.add_argument("--stations", type=int, default=16)
    p.add_argument("--days", type=int, default=3)
    p.add_argument("--fs", type=float, default=50.0)
    p.add_argument("--waves", type=int, default=150, help="plane waves per hour")
    p.add_argument("--layout", choices=("area", "line"), default="area",
                   help="area: stations at random in a square; line: along a straight line with small offsets")
    p.add_argument("--model", choices=sorted(MODELS), default="rock",
                   help="rock: 9 km array on hard rock (2.3-3.4 km/s); sediment: 2 km array on a slow basin (0.25-2.5 km/s)")
    p.add_argument("--hours", type=int, default=24, help="hours of record per day (the rest of the day is a gap)")
    p.add_argument("--responses", action="store_true",
                   help="record alternate stations with two geophone types (2 Hz and 10 Hz, damping 0.7) and write "
                        "their responses (responses.xml, data.response); synthetic_noresp.yaml leaves them in")
    p.add_argument("--directional", type=float, default=0.0,
                   help="fraction of the plane waves arriving from one direction (azimuth 45 degrees, standard "
                        "deviation 15 degrees); the rest from random azimuths (0: an isotropic field)")
    p.add_argument("--quick", action="store_true",
                   help="a small first test: 10 stations, one day with 6 hours of record at 25 Hz (about a minute to "
                        "make, a few minutes to process; the numbers are rougher than those of the full example)")
    a = p.parse_args()
    if a.quick:
        a.stations, a.days, a.hours, a.fs, a.waves = 10, 1, 6, 25.0, 100
    M = MODELS[a.model]
    TRUE_Z, TRUE_VS = np.array(M["z"]), np.array(M["vs"])
    f_lo, f_hi = M["band"]
    import obspy
    from lant.inversion.start import forward
    rng = np.random.default_rng()                 # no fixed seed: every run draws a new noise field and layout
    os.makedirs(a.out, exist_ok=True)

    # stations: random in a square (9 km rock, 2 km sediment), at least 0.4 / 0.08 km apart
    xy = []
    while len(xy) < a.stations and a.layout == "area":
        q = rng.uniform(0, M["size"], 2)
        if all(np.hypot(*(q - r)) > M["spacing"] for r in xy):
            xy.append(q)
    if a.layout == "line":      # a line across the square (azimuth 60 degrees), stations up to 1 % of its length off it
        s_ = np.linspace(0, M["size"], a.stations) + rng.uniform(-0.2, 0.2, a.stations) * M["size"] / a.stations
        off = rng.uniform(-0.01, 0.01, a.stations) * M["size"]
        xy = np.c_[s_ * np.sin(np.radians(60)) + off * np.cos(np.radians(60)),
                   s_ * np.cos(np.radians(60)) - off * np.sin(np.radians(60))]
    xy = np.array(xy)
    from lant.signal import km_per_degree     # WGS84, as lant measures the distances
    kx, ky = km_per_degree(LAT0)
    lon = LON0 + xy[:, 0] / kx
    lat = LAT0 + xy[:, 1] / ky
    names = ["S%02d" % (i + 1) for i in range(a.stations)]
    # station table deliberately without a header and with longitude before latitude: lant recognises it
    with open(os.path.join(a.out, "stations.txt"), "w") as f:
        for n, lo, la in zip(names, lon, lat):
            f.write("%s  %.6f  %.6f  %.0f\n" % (n, lo, la, 500 + 100 * rng.standard_normal()))
    # positions recomputed from the written coordinates (exactly what the processing will use), in km
    x = (lon - LON0) * kx
    y = (lat - LAT0) * ky

    # true phase velocity of the model (DSurfTomo's own forward), interpolated in log f
    freq = np.geomspace(f_lo, f_hi, 60)               # increasing frequency (np.interp needs increasing abscissae)
    c = forward(TRUE_Z, TRUE_VS, 1 / freq)
    np.savetxt(os.path.join(a.out, "true_model.txt"), np.c_[TRUE_Z, TRUE_VS], fmt="%.3f %.3f", header="depth_km Vs_km/s")
    np.savetxt(os.path.join(a.out, "true_dispersion.txt"), np.c_[freq, c], fmt="%.4f %.4f",
               header="frequency_Hz Rayleigh_phase_velocity_km/s (fundamental mode of true_model.txt)")
    fl, cl = np.log(freq), c

    nh = int(3600 * a.fs)
    f = np.fft.rfftfreq(nh, 1 / a.fs)
    band = (f >= f_lo) & (f <= f_hi)
    cf = np.ones_like(f)
    cf[band] = np.interp(np.log(f[band]), fl, cl)
    taper = np.zeros_like(f)
    taper[band] = 1.0
    t0 = obspy.UTCDateTime(2024, 6, 1)
    resp = {}
    if a.responses:              # geophones: two zeros at 0, two poles of natural frequency f0 and damping h
        from obspy.core.inventory import Inventory, Network, Station, Channel
        from obspy.core.inventory.response import Response
        stas = []
        for s, n in enumerate(names):
            f0, h = (2.0, 0.7) if s % 2 == 0 else (10.0, 0.7)
            w0 = 2 * np.pi * f0
            resp[n] = Response.from_paz(zeros=[0j, 0j], poles=[complex(-h * w0, w0 * np.sqrt(1 - h * h)),
                                                               complex(-h * w0, -w0 * np.sqrt(1 - h * h))],
                                        stage_gain=1e4, stage_gain_frequency=20.0, normalization_frequency=20.0)
            stas.append(Station(n, lat[s], lon[s], 500.0, channels=[Channel("HHZ", "", lat[s], lon[s], 500.0, 0.0,
                                                                            response=resp[n], sample_rate=a.fs)]))
        Inventory(networks=[Network("SY", stations=stas)], source="make_synthetic").write(
            os.path.join(a.out, "responses.xml"), format="STATIONXML")
    for d in range(a.days):
        day = np.zeros((a.stations, a.hours * nh), np.float32)
        for h in range(a.hours):
            az = rng.uniform(0, 2 * np.pi, a.waves)
            one = rng.random(a.waves) < a.directional          # waves from one direction (directional noise)
            az[one] = np.radians(45.0 + 15.0 * rng.standard_normal(one.sum()))
            spec = (rng.standard_normal((a.waves, f.size)) + 1j * rng.standard_normal((a.waves, f.size))) * taper
            for s in range(a.stations):
                proj = x[s] * np.sin(az) + y[s] * np.cos(az)                  # km along each propagation direction
                ph = np.exp(-2j * np.pi * f[None, :] * proj[:, None] / cf[None, :])
                day[s, h * nh:(h + 1) * nh] = np.fft.irfft((spec * ph).sum(0), n=nh).astype(np.float32)
        day += 0.05 * day.std() * rng.standard_normal(day.shape).astype(np.float32)      # local incoherent noise
        if d == 1 and a.hours >= 2:
            day[3, 30000:30000 + int(600 * a.fs)] *= 40.0                                # a burst (segment rejection)
        if resp:                 # ground velocity -> counts, through the response of each station
            nn = day.shape[1]
            fd = np.fft.rfftfreq(nn, 1 / a.fs)
            for s, n in enumerate(names):
                H = np.r_[0, resp[n].get_evalresp_response_for_frequencies(fd[1:], output="VEL")]
                day[s] = np.fft.irfft(np.fft.rfft(day[s]) * H / 1e4, nn).astype(np.float32)
        for s, n in enumerate(names):
            tr = obspy.Trace(day[s] * 1e4, header=dict(network="SY", station=n, location="", channel="HHZ",
                                                       sampling_rate=a.fs, starttime=t0 + 86400 * d))
            if s % 2 == 0:              # miniSEED, one folder per station
                dd = os.path.join(a.out, "data", "mseed", n)
                os.makedirs(dd, exist_ok=True)
                tr.write(os.path.join(dd, "%s_day%d.mseed" % (n, d + 1)), format="MSEED", encoding="FLOAT32")
            else:                       # SAC without network code, all in one folder, other names
                dd = os.path.join(a.out, "data", "sac_files", "day%02d" % (d + 1))
                os.makedirs(dd, exist_ok=True)
                tr.stats.network = ""
                tr.write(os.path.join(dd, "rec_%s.sac" % n), format="SAC")
        print("day %d written" % (d + 1), flush=True)
    with open(os.path.join(a.out, "synthetic.yaml"), "w") as fh:
        if a.quick:
            fh.write(YAML_QUICK % dict(vlo=M["vrange"][0], vhi=M["vrange"][1]))
        else:
            fh.write(YAML_LINE % dict(vlo=M["vrange"][0], vhi=M["vrange"][1]) if a.layout == "line" else YAML % dict(name="synthetic" if a.model == "rock" else "synthetic_" + a.model,
                             cell=round(M["size"] / 3, 2), vlo=M["vrange"][0], vhi=M["vrange"][1],
                             export="export: {tmin_s: 0.2}\n" if a.model == "rock" else ""))
    if a.responses:              # the same configuration with the responses, and one without them (output_noresp)
        y = open(os.path.join(a.out, "synthetic.yaml")).read()
        open(os.path.join(a.out, "synthetic_noresp.yaml"), "w").write(y.replace("output_dir: output", "output_dir: output_noresp"))
        open(os.path.join(a.out, "synthetic.yaml"), "w").write(
            y.replace("station_file: stations.txt,", "station_file: stations.txt, response: responses.xml,"))
    print("wrote %s (stations.txt, data/, true_model.txt, true_dispersion.txt, synthetic.yaml%s)"
          % (a.out, ", responses.xml, synthetic_noresp.yaml" if a.responses else ""))


if __name__ == "__main__":
    main()
