#!/usr/bin/env python
"""Synthetic correlations with a known answer, written directly as a stack file (03_stack/stack.h5).

make_synthetic.py makes raw noise records of a laterally UNIFORM medium (plane waves).  This script makes the
correlations themselves, for cases plane waves cannot give:

  * a laterally varying medium: Vs(x, y, z) = Vs0(z) (1 + d(x, y)), d a chequer pattern of amplitude --hetero and
    cell size --cell; the phase velocity c(f, d) comes from the 1-D forward calculation (fwd1d) for the scaled
    profile, and the travel time of a pair is the integral of 1 / c along the straight line between the stations;
  * a second, faster arrival (--higher-mode AMPLITUDE: 1.25 x the phase velocity, above --mode-fmin Hz), to see
    what the measurement and its quality control do with a higher mode;
  * incoherent noise (--noise: its standard deviation in the stack relative to the largest amplitude of each
    correlation), drawn anew for every sub-stack.

The correlation spectrum of a pair is J0(2 pi f t(f)) (the real part of the 2-D Green's function, far field
cos(2 pi f t - pi/4)), in the band of the template.  Geometry, sampling, lag window, sub-stacks and bookkeeping
are copied from an existing stack (--like), so the file can replace the 03_stack/stack.h5 of a project and the
chain be run again from the reference step:

    python examples/synthetic/make_stack.py --like OUT/03_stack/stack.h5 --out OUT/03_stack/stack.h5 \
        --hetero 0.06 --cell 3 --truth OUT/truth.npz
    lant run all --from reference -c CONFIG

The truth file holds, per pair and frequency, the path-average phase velocity, and the pattern d(x, y).
"""
import argparse
import os
import sys

import h5py
import numpy as np
from scipy.special import j0

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))


def pattern(x, y, cell, amp):
    """Chequer pattern of relative Vs perturbation (x, y in km from the array centre)."""
    return amp * np.sin(np.pi * x / cell) * np.sin(np.pi * y / cell)


def velocity_table(z, vs, freqs, amp):
    """c(f) of the profile scaled by 1 - amp, 1, 1 + amp (rows) at the frequencies freqs."""
    from lant.inversion.start import forward
    return np.array([forward(z, np.asarray(vs) * (1 + a), 1 / freqs) for a in (-amp, 0.0, amp)]) if amp > 0 \
        else np.array([forward(z, np.asarray(vs), 1 / freqs)] * 3)


def build(like, out, z, vs, hetero=0.0, cell=3.0, noise=0.0, higher_mode=0.0, mode_fmin=4.0, truth=None):
    from lant.signal import local_xy_km
    rng = np.random.default_rng()
    with h5py.File(like, "r") as h:
        meta = {k: h[k][:] for k in ("pairs", "days", "dist_km", "az_deg", "lon1", "lat1", "lon2", "lat2",
                                     "daily_count", "nseg")}
        attrs = dict(h.attrs)
        nt = h["stack"].shape[1]
    fs, fmin, fmax = float(attrs["fs"]), float(attrs["fmin"]), float(attrs["fmax"])
    L = nt // 2
    nfft = 2 * L
    f = np.fft.rfftfreq(nfft, 1 / fs)
    band = (f >= fmin) & (f <= fmax)
    fb = f[band]
    grid = np.geomspace(fb[0], fb[-1], 40)                    # c(f) on a coarse grid, interpolated in log f
    tab = velocity_table(z, vs, grid, hetero)
    lon0, lat0 = np.mean(np.r_[meta["lon1"], meta["lon2"]]), np.mean(np.r_[meta["lat1"], meta["lat2"]])
    x1, y1 = local_xy_km(meta["lon1"], meta["lat1"], lat0, lon0)
    x2, y2 = local_xy_km(meta["lon2"], meta["lat2"], lat0, lon0)
    s = np.linspace(0, 1, 201)
    npair = len(meta["dist_km"])
    taper = np.ones(len(fb))
    ne = max(2, len(fb) // 20)
    taper[:ne] = 0.5 - 0.5 * np.cos(np.pi * np.arange(ne) / ne)
    taper[-ne:] = taper[:ne][::-1]
    stack = np.zeros((npair, nt), np.float32)
    c_path = np.zeros((npair, len(grid)))
    for q in range(npair):
        d = pattern(x1[q] + s * (x2[q] - x1[q]), y1[q] + s * (y2[q] - y1[q]), cell, hetero) if hetero > 0 else np.zeros(len(s))
        w = d / hetero if hetero > 0 else d                   # -1 .. 1: position between the three tabulated profiles
        slow = np.where(w[:, None] < 0, (1 + w[:, None]) / tab[1] - w[:, None] / tab[0],
                        (1 - w[:, None]) / tab[1] + w[:, None] / tab[2])
        c_path[q] = 1.0 / slow.mean(axis=0)
        t = meta["dist_km"][q] / np.interp(np.log(fb), np.log(grid), c_path[q])
        S = np.zeros(len(f))
        S[band] = j0(2 * np.pi * fb * t) * taper
        if higher_mode > 0:
            S[band] += higher_mode * j0(2 * np.pi * fb * t / 1.25) * (fb >= mode_fmin) * taper
        xx = np.fft.irfft(S, nfft)
        stack[q] = np.r_[xx[-L:], xx[:L + 1]]
    stack /= np.abs(stack).max(axis=1, keepdims=True) + 1e-30
    nblk = len(meta["days"])
    S_n = np.zeros(len(f))
    S_n[band] = taper
    with h5py.File(out + ".tmp", "w") as o:
        dd = o.create_dataset("daily", (nblk, npair, nt), np.float32, chunks=(1, npair, nt))
        acc = np.zeros((npair, nt))
        cnt = meta["daily_count"].astype(np.float64)
        for b in range(nblk):
            blk = stack.astype(np.float64)
            if noise > 0:                                     # band-limited noise, new for every sub-stack
                w = np.fft.irfft(np.fft.rfft(rng.standard_normal((npair, nfft)), axis=1) * S_n, nfft, axis=1)
                w = np.c_[w[:, -L:], w[:, :L + 1]]
                blk = blk + noise * np.sqrt(nblk) * w / (w.std() + 1e-30)
            blk[cnt[b] == 0] = 0
            dd[b] = blk
            acc += blk * cnt[b][:, None]
        o["stack"] = (acc / np.maximum(cnt.sum(0), 1)[:, None]).astype(np.float32)
        for k, v in meta.items():
            o[k] = v
        for k, v in attrs.items():
            o.attrs[k] = v
    os.replace(out + ".tmp", out)
    if truth:
        np.savez(truth, freqs=grid, c_path=c_path, c0=tab[1], pairs=meta["pairs"].astype("U"), hetero=hetero, cell=cell,
                 lon0=lon0, lat0=lat0, z=np.asarray(z), vs=np.asarray(vs))
    return grid, c_path


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--like", required=True, help="an existing stack.h5 (geometry, sampling, sub-stacks)")
    p.add_argument("--out", required=True)
    p.add_argument("--model", default=os.path.join(HERE, "true_model.txt"), help="depth_km Vs_km/s")
    p.add_argument("--hetero", type=float, default=0.0, help="amplitude of the lateral Vs pattern (0.06 = 6 %%)")
    p.add_argument("--cell", type=float, default=3.0, help="cell size of the pattern (km)")
    p.add_argument("--noise", type=float, default=0.0)
    p.add_argument("--higher-mode", type=float, default=0.0)
    p.add_argument("--mode-fmin", type=float, default=4.0)
    p.add_argument("--truth", default=None)
    a = p.parse_args()
    m = np.loadtxt(a.model)
    build(a.like, a.out, m[:, 0], m[:, 1], a.hetero, a.cell, a.noise, a.higher_mode, a.mode_fmin, a.truth)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
