"""Helpers of the tests: configuration files, temporary directories, a small set of raw records."""
import textwrap

import numpy as np

REQ = "velocity_range_kms: [0.5, 3.0]"


def write_cfg(tmp_path, text):
    import pathlib
    tmp_path = pathlib.Path(tmp_path)
    p = tmp_path / "cfg.yaml"
    p.write_text(textwrap.dedent(text))
    return str(p)


def tmp_dir():
    import tempfile
    return tempfile.mkdtemp()


def make_raw(tmp_path):
    import obspy
    rng = np.random.default_rng(0)
    t0 = obspy.UTCDateTime(2024, 1, 1)
    raw = tmp_path / "raw"
    for i, name in enumerate(["A01", "A02", "A03"]):
        tr = obspy.Trace(rng.standard_normal(int(3 * 3600 * 100)).astype(np.float32),
                         header=dict(network="", station=name, channel="HHZ", sampling_rate=100.0,
                                     starttime=t0 + (0.004 if i == 2 else 0.0003)))
        d = raw / ("deep" if i else "flat") / "x"
        d.mkdir(parents=True, exist_ok=True)
        if i == 0:
            tr.stats.network = "ZZ"
            tr.write(str(d / "rec1.mseed"), format="MSEED", encoding="FLOAT32")
            tr.copy().write(str(d / "dup.mseed"), format="MSEED", encoding="FLOAT32")        # exact duplicate
        else:
            tr.write(str(d / ("%s.sac" % name)), format="SAC")
    (raw / "notes.txt").write_text("not data")
    (tmp_path / "sta.csv").write_text("network,station,latitude,longitude,elevation,instrument\n"
                                      "ZZ,A01,30.0,100.0,10,N1\nZZ,A02,30.01,100.0,20,N1\nZZ,A03,30.0,100.01,,N1\n")
    return raw
