"""Instrument responses (data.response): read, matched to the records, removed in the clean step.

The cross-coherence keeps only the phase of the two records, so a response common to all stations cancels in
every pair and need not be removed; responses that differ between stations do not cancel and are removed first
(Bensen et al. 2007, section 2.1).  The removal is ObsPy's Trace.remove_response: deconvolution in the frequency
domain to ground velocity, with a cosine pre-filter that leaves the correlation band untouched (flat between
correlation.fmin_hz and correlation.fmax_hz, zero below half the lower and above 1.25 times the upper edge) in place
of a water level.
"""

import glob
import os

_CACHE = {}


class ResponseError(RuntimeError):
    pass


def load(path):
    """Inventory of a response file, or of every readable file of a directory (StationXML, dataless SEED, RESP,
    SACPZ ...: what obspy.read_inventory reads)."""
    if path in _CACHE:
        return _CACHE[path]
    import obspy

    files = (
        [path]
        if os.path.isfile(path)
        else sorted(f for f in glob.glob(os.path.join(path, "**", "*"), recursive=True) if os.path.isfile(f))
    )
    if not files:
        raise ResponseError("data.response: no file at %s" % path)
    inv, bad = None, []
    for f in files:
        try:
            part = obspy.read_inventory(f)
        except Exception as exc:
            bad.append("%s (%s)" % (f, str(exc)[:60]))
            continue
        inv = part if inv is None else inv + part
    if inv is None:
        raise ResponseError(
            "data.response: none of the %d files under %s is a response file ObsPy reads: %s"
            % (len(files), path, "; ".join(bad[:3]))
        )
    _CACHE[path] = inv
    return inv


def lookup(inv, net, sta, loc, cha, t):
    """The SEED id of the response for a record with these header codes at time t (UTC timestamp), or None.  The
    exact id first; otherwise the station code alone when it leaves one channel (the same channel code, or the only
    vertical one) - headers of node data often lack the network or location code."""
    import obspy

    t = obspy.UTCDateTime(t)
    sid = "%s.%s.%s.%s" % (net, sta, loc, cha)
    try:
        inv.get_response(sid, t)
        return sid
    except Exception:
        pass
    found = []
    for n in inv.select(station=sta, time=t):
        for s in n:
            for c in s:
                if c.response is not None:
                    found.append((n.code, s.code, c.location_code, c.code))
    same = [f for f in found if f[3] == cha] or [f for f in found if f[3][-1:] in ("Z", cha[-1:])]
    if net:
        same = [f for f in same if f[0] == net] or same
    ids = sorted({"%s.%s.%s.%s" % f for f in same})
    return ids[0] if len(ids) == 1 else None


def pre_filter(cfg):
    f2, f3 = float(cfg["correlation"]["fmin_hz"]), float(cfg["correlation"]["fmax_hz"])
    return [0.5 * f2, f2, f3, 1.25 * f3]


def remove(tr, inv, seed_id, pre_filt):
    """Remove the response of `seed_id` from a gap-free trace, in place (-> ground velocity, float64)."""
    tr.stats.response = inv.get_response(seed_id, tr.stats.starttime)
    dur = tr.stats.npts * tr.stats.delta
    nyq = 0.5 * tr.stats.sampling_rate
    pf = list(pre_filt)
    if pf[3] > nyq:
        pf[3] = nyq
        pf[2] = min(pf[2], 0.8 * nyq)
    frac = min(0.05, (2.0 / pf[0]) / dur)  # end taper: two periods of the lowest frequency kept
    tr.remove_response(output="VEL", pre_filt=pf, water_level=None, zero_mean=True, taper=True, taper_fraction=frac)
    del tr.stats.response
    return tr
