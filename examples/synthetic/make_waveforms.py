"""Continuous vertical-component noise records of a dense array over a 1-D model (part A of the synthetic example).

    python make_waveforms.py MODEL OUT [--hours 24] [--cores 8]     MODEL: layered | lvl (models.txt)

Writes OUT/raw/<station>.mseed (one file per station, 50 Hz, float32), OUT/stations.txt and OUT/info.txt.

The wavefield (after the noise simulations of the literature, e.g. Herrmann's CPS tutorial on high-frequency noise,
Tsai 2009, Yao & van der Hilst 2009): point sources with white random source time functions, uniformly distributed
over a disk three times the array radius around its centre (inside the array too, but not within 25 m of a station),
60 new sources every 100 s.  Each source's
vertical force gives at a station at distance r the sum over the Rayleigh modes of the model (modes_<MODEL>.csv,
fundamental and higher modes with their excitation AR):  u(f) = sum_m AR_m H0^(2)(k_m r) exp(-pi f r / (Q U_m)),
k_m = 2 pi f / c_m, Q = 50.  On every record: independent local noise (half the rms of the wavefield), a few
amplitude glitches, two gaps; one station of the table has no records.  These are the things lant clean and lant cc
have to cope with.
"""

import argparse
import os
import time
from multiprocessing import Pool

import modes as M
import numpy as np
from scipy.special import hankel2

FS = 50.0
NFFT = 8192  # 163.84 s: a 100 s source window plus the longest delay that still carries energy
CHUNK_S = 100.0
SOURCES = 60  # per chunk
DMIN = 0.025  # km: no source closer to a station (as in Herrmann's CPS tutorial)
Q = 50.0
BAND = (0.2, 22.0)
LON0, LAT0 = 105.0, 35.0
KM = 111.195
START = "2024-01-01T00:00:00"


def stations(seed=1):
    """10 x 10 nodes 0.55 km apart, each moved at random by up to a third of the spacing -> (x, y) km, names."""
    rng = np.random.default_rng(seed)
    g = (np.arange(10) - 4.5) * 0.55
    X, Y = np.meshgrid(g, g)
    xy = np.c_[X.ravel(), Y.ravel()] + rng.uniform(-0.18, 0.18, (100, 2))
    return xy, ["S%03d" % (k + 1) for k in range(len(xy))]


def to_lonlat(xy):
    return LON0 + xy[:, 0] / (KM * np.cos(np.radians(LAT0))), LAT0 + xy[:, 1] / KM


def h02(x):
    """H0^(2)(x): exact below 10, the asymptotic series with two terms above (relative error < 1e-4)."""
    out = np.empty(x.shape, np.complex128)
    small = x < 10
    out[small] = hankel2(0, x[small])
    xs = x[~small]
    out[~small] = np.sqrt(2 / (np.pi * xs)) * np.exp(-1j * (xs - np.pi / 4)) * (1 + 1j / (8 * xs) - 9 / (128 * xs**2))
    return out


def mode_spectra(table, f):
    """Per mode, on the frequencies f where it exists and is at least 1 % of the fundamental: the frequency mask,
    k (rad/km), AR and the attenuation (1/km)."""
    out = []
    ar0 = np.interp(f, table[0][0], table[0][3])
    for m, (fm, c, u, ar) in table.items():
        a = np.interp(f, fm, ar)
        use = (f >= fm[0]) & (f <= fm[-1]) & ((a >= 0.01 * ar0) | (m == 0))
        if use.any():
            fu = f[use]
            out.append((use, 2 * np.pi * fu / np.interp(fu, fm, c), a[use], np.pi * fu / (Q * np.interp(fu, fm, u))))
    return out


G = {}


def _init(table, xy, seed):
    f = np.fft.rfftfreq(NFFT, 1 / FS)
    band = np.where((f >= BAND[0]) & (f <= BAND[1]))[0]
    G.update(f=f, band=band, modes=mode_spectra(table, f[band]), xy=xy, seed=seed)
    lo = np.clip((f[band] - BAND[0]) / 0.1, 0, 1)  # cosine tapers over 0.1 Hz and 2 Hz at the band edges
    hi = np.clip((BAND[1] - f[band]) / 2.0, 0, 1)
    G["taper"] = 0.5 * (1 - np.cos(np.pi * lo)) * 0.5 * (1 - np.cos(np.pi * hi))


def chunk(i):
    """The records of the sources of chunk i -> (i, array [stations, NFFT])."""
    rng = np.random.default_rng((G["seed"], i + 10))
    xy, f, band = G["xy"], G["f"], G["band"]
    r0 = np.hypot(*xy.T).max()
    src = np.empty((0, 2))
    while len(src) < SOURCES:  # uniform over a disk 3 times the array radius, not within DMIN of a station
        rad = np.sqrt(rng.uniform(0, (3.0 * r0) ** 2, SOURCES))
        az = rng.uniform(0, 2 * np.pi, SOURCES)
        p = np.c_[rad * np.cos(az), rad * np.sin(az)]
        near = np.hypot(p[:, None, 0] - xy[None, :, 0], p[:, None, 1] - xy[None, :, 1]).min(axis=1) < DMIN
        src = np.r_[src, p[~near]]
    src = src[:SOURCES]
    n = int(CHUNK_S * FS)
    stf = rng.standard_normal((SOURCES, n))
    S = np.fft.rfft(stf, NFFT, axis=1)[:, band] * G["taper"]
    U = np.zeros((len(xy), band.size), np.complex128)
    for s in range(SOURCES):
        r = np.hypot(*(xy - src[s]).T)[:, None]
        g = np.zeros((len(xy), band.size), np.complex128)
        for use, k, ar, att in G["modes"]:
            g[:, use] += ar * h02(k * r) * np.exp(-att * r)
        U += S[s] * g
    full = np.zeros((len(xy), f.size), np.complex128)
    full[:, band] = U
    return i, np.fft.irfft(full, NFFT, axis=1).astype(np.float32)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("model", choices=sorted(M.read_models()))
    ap.add_argument("out")
    ap.add_argument("--hours", type=float, default=24.0)
    ap.add_argument("--cores", type=int, default=os.cpu_count())
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    import obspy

    table = M.read_table(os.path.join(M.HERE, "modes_%s.csv" % a.model))
    xy, names = stations()
    npts = int(round(a.hours * 3600 * FS))
    nchunk = int(np.ceil(a.hours * 3600 / CHUNK_S))
    data = np.zeros((len(xy), npts + NFFT), np.float32)
    t0 = time.time()
    with Pool(a.cores, initializer=_init, initargs=(table, xy, a.seed)) as pool:
        for k, (i, rec) in enumerate(pool.imap_unordered(chunk, range(-2, nchunk))):  # two chunks before the start
            j = int(i * CHUNK_S * FS)
            lo, hi = max(0, j), min(data.shape[1], j + NFFT)
            if hi > lo:
                data[:, lo:hi] += rec[:, lo - j : hi - j]
            if (k + 1) % max(1, nchunk // 10) == 0:
                print("  %d of %d source windows (%.0f s)" % (k + 1, nchunk + 2, time.time() - t0), flush=True)
    data = data[:, :npts]
    rng = np.random.default_rng(a.seed + 1)
    rms = np.sqrt(np.mean(data[:, ::50].astype(np.float64) ** 2, axis=1))
    from scipy.signal import butter, sosfiltfilt

    sos = butter(4, BAND, btype="band", fs=FS, output="sos")
    for k in range(len(xy)):  # local noise, independent at every station
        data[k] += (0.5 * rms[k] * sosfiltfilt(sos, rng.standard_normal(npts))).astype(np.float32)
    for k in rng.choice(len(xy), 6, replace=False):  # glitches: 2 s bursts 50 times the rms
        j = int(rng.uniform(0, npts - 2 * FS))
        data[k, j : j + int(2 * FS)] += (50 * rms[k] * rng.standard_normal(int(2 * FS))).astype(np.float32)
    lon, lat = to_lonlat(xy)
    os.makedirs(os.path.join(a.out, "raw"), exist_ok=True)
    gaps = {names[7]: (3.0, 5.5), names[42]: (13.0, 13.4)}  # hours
    dead = names[90]
    for k, nm in enumerate(names):
        if nm == dead:
            continue
        st = obspy.Stream()
        segs = [(0.0, a.hours)] if nm not in gaps else [(0.0, gaps[nm][0]), (gaps[nm][1], a.hours)]
        for h0, h1 in segs:
            j0, j1 = int(h0 * 3600 * FS), min(npts, int(h1 * 3600 * FS))
            if j1 - j0 < 2:
                continue
            tr = obspy.Trace(data[k, j0:j1].copy())
            tr.stats.network, tr.stats.station, tr.stats.channel = "SY", nm, "HHZ"
            tr.stats.sampling_rate = FS
            tr.stats.starttime = obspy.UTCDateTime(START) + h0 * 3600
            st.append(tr)
        st.write(os.path.join(a.out, "raw", "%s.mseed" % nm), format="MSEED", encoding="FLOAT32")
    with open(os.path.join(a.out, "stations.txt"), "w") as fh:
        fh.write("network station longitude latitude elevation\n")
        for nm, x, y in zip(names, lon, lat):
            fh.write("SY %s %.6f %.6f 0\n" % (nm, x, y))
    with open(os.path.join(a.out, "info.txt"), "w") as fh:
        fh.write("model %s; %d stations (%s without records); %.1f h at %g Hz; %d source windows of %g s with %d "
                 "sources; Q %g; gaps %s\n" % (a.model, len(names), dead, a.hours, FS, nchunk, CHUNK_S, SOURCES, Q,
                                                gaps))  # fmt: skip
    print("written: %s (%.0f s)" % (a.out, time.time() - t0))


if __name__ == "__main__":
    main()
