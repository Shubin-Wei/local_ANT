"""lant export-sac: the stacked correlations as SAC files, for other programs (e.g. EGFAnalysisTimeFreq).

One file per pair, STA1_STA2.sac, two-sided (b = -maxlag), with the SAC headers of a correlation: the first
station as the "event" (evlo, evla, kevnm) and the second as the station (stlo, stla, kstnm), dist (km), az, and
user0 = number of segments stacked.  Positive lags: energy travelling from station 1 to station 2.
"""

import os

import h5py
import numpy as np


def run(cfg, out_dir, min_count=1):
    import obspy
    from obspy.io.sac import SACTrace

    stack = os.path.join(cfg["output_dir"], "1_cc", "stack.h5")
    if not os.path.exists(stack):
        raise FileNotFoundError("%s not found: run `lant cc` first" % stack)
    os.makedirs(out_dir, exist_ok=True)
    n = 0
    with h5py.File(stack, "r") as h:
        fs, maxlag = float(h.attrs["fs"]), float(h.attrs["maxlag"])
        cnt = h["count"][:]
        s1, s2 = [x.decode() for x in h["sta1"][:]], [x.decode() for x in h["sta2"][:]]
        for q in np.where(cnt >= min_count)[0]:
            tr = obspy.Trace(np.asarray(h["ccf"][q], np.float32))
            tr.stats.delta = 1.0 / fs
            sac = SACTrace.from_obspy_trace(tr)
            sac.b = -maxlag
            sac.evlo, sac.evla = float(h["lon1"][q]), float(h["lat1"][q])
            sac.stlo, sac.stla = float(h["lon2"][q]), float(h["lat2"][q])
            sac.kevnm, sac.kstnm = s1[q][-16:], s2[q].split(".")[-1][-8:]
            sac.dist, sac.az = float(h["dist_km"][q]), float(h["az_deg"][q])
            sac.user0 = float(cnt[q])
            sac.lcalda = False
            sac.write(os.path.join(out_dir, "%s_%s.sac" % (s1[q], s2[q])))
            n += 1
    print("%d SAC files written to %s" % (n, out_dir))
    return n
