#!/usr/bin/env python
"""Synthetic test data with a known answer: ambient noise of a laterally uniform medium.

The noise field is a sum of plane Rayleigh waves arriving from random azimuths (an isotropic diffuse field), each
with a random spectrum, propagating with the phase velocity c(f) of a known 1-D Vs model.  For such a field the
cross-correlation of two stations is the 2-D Green's function, so the measured phase velocities should follow c(f)
and the inversion should return the 1-D model (within its resolution).

    python examples/synthetic/make_synthetic.py [--out examples/synthetic] [--stations 16] [--days 3]

writes
    data/                       the "raw" records, deliberately mixed: miniSEED and SAC (SAC without network code),
                                in nested folders with arbitrary names
    stations.txt                the station table (no header, longitude before latitude: recognised automatically)
    true_model.txt, true_dispersion.txt   the answer (depth Vs; frequency c)
    synthetic.yaml              a ready configuration (output in examples/synthetic/output)
Needs fwd1d (bash tools/install_external.sh dsurftomo) to compute c(f) of the model.  About 5-20 minutes, 300 MB.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

TRUE_Z = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.9, 2.4, 3.0])
TRUE_VS = np.array([2.30, 2.45, 2.62, 2.78, 2.90, 3.00, 3.08, 3.16, 3.25, 3.33, 3.40])
LON0, LAT0 = 100.0, 30.0          # an arbitrary place (no real array)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default=HERE)
    p.add_argument("--stations", type=int, default=16)
    p.add_argument("--days", type=int, default=3)
    p.add_argument("--fs", type=float, default=50.0)
    p.add_argument("--waves", type=int, default=150, help="plane waves per hour")
    p.add_argument("--seed", type=int, default=1)
    a = p.parse_args()
    import obspy
    from lant.inversion.start import forward
    rng = np.random.default_rng(a.seed)
    os.makedirs(a.out, exist_ok=True)

    # stations: random in a 9 x 9 km square, at least 0.4 km apart
    xy = []
    while len(xy) < a.stations:
        q = rng.uniform(0, 9.0, 2)
        if all(np.hypot(*(q - r)) > 0.4 for r in xy):
            xy.append(q)
    xy = np.array(xy)
    k = np.radians(1.0) * 6371.0
    lon = LON0 + xy[:, 0] / (k * np.cos(np.radians(LAT0)))
    lat = LAT0 + xy[:, 1] / k
    names = ["S%02d" % (i + 1) for i in range(a.stations)]
    # station table deliberately without a header and with longitude before latitude: lant recognises it
    with open(os.path.join(a.out, "stations.txt"), "w") as f:
        for n, lo, la in zip(names, lon, lat):
            f.write("%s  %.6f  %.6f  %.0f\n" % (n, lo, la, 500 + 100 * rng.standard_normal()))
    # positions recomputed from the written coordinates (exactly what the processing will use), in km
    x = (lon - LON0) * k * np.cos(np.radians(LAT0))
    y = (lat - LAT0) * k

    # true phase velocity of the model (DSurfTomo's own forward), interpolated in log f
    freq = np.geomspace(0.3, 12.0, 60)                 # increasing frequency (np.interp needs increasing abscissae)
    c = forward(TRUE_Z, TRUE_VS, 1 / freq)
    np.savetxt(os.path.join(a.out, "true_model.txt"), np.c_[TRUE_Z, TRUE_VS], fmt="%.3f %.3f", header="depth_km Vs_km/s")
    np.savetxt(os.path.join(a.out, "true_dispersion.txt"), np.c_[freq, c], fmt="%.4f %.4f",
               header="frequency_Hz Rayleigh_phase_velocity_km/s (fundamental mode of true_model.txt)")
    fl, cl = np.log(freq), c

    nh = int(3600 * a.fs)
    f = np.fft.rfftfreq(nh, 1 / a.fs)
    band = (f >= 0.3) & (f <= 12.0)
    cf = np.ones_like(f)
    cf[band] = np.interp(np.log(f[band]), fl, cl)
    taper = np.zeros_like(f)
    taper[band] = 1.0
    t0 = obspy.UTCDateTime(2024, 6, 1)
    for d in range(a.days):
        day = np.zeros((a.stations, 24 * nh), np.float32)
        for h in range(24):
            az = rng.uniform(0, 2 * np.pi, a.waves)
            spec = (rng.standard_normal((a.waves, f.size)) + 1j * rng.standard_normal((a.waves, f.size))) * taper
            for s in range(a.stations):
                proj = x[s] * np.sin(az) + y[s] * np.cos(az)                  # km along each propagation direction
                ph = np.exp(-2j * np.pi * f[None, :] * proj[:, None] / cf[None, :])
                day[s, h * nh:(h + 1) * nh] = np.fft.irfft((spec * ph).sum(0), n=nh).astype(np.float32)
        day += 0.05 * day.std() * rng.standard_normal(day.shape).astype(np.float32)      # local incoherent noise
        if d == 1:
            day[3, 30000:30000 + int(600 * a.fs)] *= 40.0                                # a burst (segment rejection)
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
        fh.write("""# Configuration of the synthetic example (examples/synthetic/make_synthetic.py).
# Only the three paths are needed for real data.  The other lines (a) lower thresholds meant for >= 50 stations,
# because this example has only 16 stations (120 pairs) and 3 days, and (b) make the inversion part fast enough
# for a first test (periods >= 0.2 s, no extra depth nodes, 6 iterations, one checkerboard and one control).
project: {name: synthetic, output_dir: output}
data: {raw_dir: data, station_file: stations.txt}
stack: {station_day_qc: {min_partners: 5}}
notch: {raw_nsta: 8, raw_ndays: 3}
reference: {pass2_min_pairs: 20, fj_threads: 4}
dispersion: {njack: 3, nproc: 4}
export: {tmin_s: 0.2}
inversion: {median_min_pairs: 20, fine_step_km: 0, iterations: 6, threads: 4, njobs: 3}
checks: {checkerboards: [{name: cb_3km, cell_km: [3.0, 3.0]}], vertical: [], controls: [start_x1.05], iteration_residuals: false}
""")
    print("wrote %s (stations.txt, data/, true_model.txt, true_dispersion.txt, synthetic.yaml)" % a.out)


if __name__ == "__main__":
    main()
