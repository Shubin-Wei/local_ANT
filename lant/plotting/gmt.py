"""GMT 6 (modern mode) for the map figures.

A figure is one shell script run in a temporary folder.  Relief (plots.relief): auto uses GMT's remote
@earth_relief_03s (downloaded and cached by GMT; needs internet), a file name uses that grid, null draws no relief.
When the relief cannot be had the map is drawn without it and the log says so.  Set LANT_KEEP_GMT=1 to keep the
scripts.
"""

import os
import shlex
import shutil
import subprocess
import tempfile

import numpy as np

SET = ("gmt set FORMAT_GEO_MAP ddd.xxF FONT_ANNOT_PRIMARY 9p,Helvetica FONT_LABEL 10p,Helvetica FONT_TITLE 11p,Helvetica MAP_FRAME_TYPE plain "
       "MAP_FRAME_PEN 0.5p MAP_TICK_LENGTH_PRIMARY 0.07c MAP_ANNOT_OFFSET_PRIMARY 0.06c PS_LINE_JOIN round")  # fmt: skip
CHECK_CPT = "polar -I"  # red (slow) - white - blue (fast)


class GMTError(RuntimeError):
    pass


def available():
    exe = shutil.which("gmt")
    if not exe:
        return None
    try:
        return subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=30).stdout.strip()
    except Exception:
        return None


def run(tmp, pre, cmd, out, dpi=300):
    """Run the commands `pre` (before the session) and `cmd` (inside it) -> out.png."""
    if not available():
        raise GMTError("GMT 6 is not installed (conda install -c conda-forge gmt): map not drawn")
    q = shlex.quote
    script = os.path.join(tmp, "draw.sh")
    with open(script, "w") as f:
        f.write("set -e\nexport GMT_TMPDIR=%s GMT_USERDIR=%s\ncd %s\n%s\ngmt begin figure png A+m0.25c,E%d\n%s\n%s\n"
                "gmt end\n" % (q(tmp), q(tmp), q(tmp), "\n".join(pre), dpi, SET, "\n".join(cmd)))  # fmt: skip
    r = subprocess.run(["bash", script], capture_output=True, text=True)
    if r.returncode != 0 or "ERROR" in r.stderr:
        raise GMTError("GMT failed for %s:\n%s\n(script: %s)" % (out, r.stderr[-2000:], script))
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    shutil.move(os.path.join(tmp, "figure.png"), out + ".png")
    if not os.environ.get("LANT_KEEP_GMT"):
        shutil.rmtree(tmp, ignore_errors=True)
    return out + ".png"


def relief_source(setting):
    if setting in (None, False, "none", "null"):
        return None
    return "@earth_relief_03s" if setting == "auto" else str(setting)


def relief_cmds(src, tmp, region, palette="gray"):
    """Commands that cut the relief to the region, its shading and a colour table (pale grey under model maps,
    coloured under the station map)."""
    W, E, S, N = region
    pad = 0.01
    pre = [
        "gmt grdcut %s -R%.5f/%.5f/%.5f/%.5f -G%s/topo.nc"
        % (shlex.quote(src), W - pad, E + pad, S - pad, N + pad, tmp),
        "gmt grdgradient %s/topo.nc -A315 -Nt%g -G%s/shade.nc" % (tmp, 0.22 if palette != "gray" else 0.5, tmp),
        "read ZMIN ZMAX <<< \"$(gmt grdinfo -C %s/topo.nc | awk '{print $6, $7}')\"" % tmp,
    ]
    if palette == "gray":
        pre.append("gmt makecpt -Cgray -T$(awk -v a=$ZMIN -v b=$ZMAX 'BEGIN{printf \"%%d/%%d\", a-2.5*(b-a)-1, "
                   "b+0.3*(b-a)+1}') -Z > %s/relief.cpt" % tmp)  # fmt: skip
    else:
        pre.append("gmt makecpt -Cgeo -T$(awk -v a=$ZMIN -v b=$ZMAX 'BEGIN{printf \"%%d/%%d\", a, (b>a)?b:a+1}') -Z "
                   "> %s/relief.cpt" % tmp)  # fmt: skip
    return pre


def draw(cfg, make, out, log=print):
    """make(tmp, relief_src) -> (pre, cmd).  Drawn with relief if possible, else without (and said so)."""
    src = relief_source(cfg["plots"]["relief"])
    for attempt in [src, None] if src else [None]:
        tmp = tempfile.mkdtemp(prefix="lant_gmt_")
        pre, cmd = make(tmp, attempt)
        try:
            return run(tmp, pre, cmd, out, int(cfg["plots"]["dpi"]))
        except GMTError as exc:
            if attempt is None:
                raise
            log("  relief %s not available (%s): %s drawn without relief" % (attempt, str(exc).splitlines()[1:2],
                                                                              os.path.basename(out)))  # fmt: skip


def region(lon, lat, pad=0.1, min_pad=0.005, min_aspect=0.4):
    """Map region around the points, the shorter side widened to min_aspect of the longer (a line is not a strip)."""
    W, E, S, N = min(lon), max(lon), min(lat), max(lat)
    px, py = max((E - W) * pad, min_pad), max((N - S) * pad, min_pad)
    W, E, S, N = W - px, E + px, S - py, N + py
    kx = np.cos(np.radians(0.5 * (S + N)))
    wx, wy = (E - W) * kx, N - S
    if wx < min_aspect * wy:
        d = 0.5 * (min_aspect * wy - wx) / kx
        W, E = W - d, E + d
    elif wy < min_aspect * wx:
        d = 0.5 * (min_aspect * wx - wy)
        S, N = S - d, N + d
    return W, E, S, N


def mercator_height(width, W, E, S, N):
    ys = lambda la: np.log(np.tan(np.pi / 4 + np.radians(la) / 2))
    return width * (ys(N) - ys(S)) / np.radians(E - W)


def nice_step(span, n=4):
    raw = span / n
    p = 10 ** np.floor(np.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if m * p >= raw:
            return m * p
    return 10 * p


def scale_km(W, E, S, N):
    width_km = (E - W) * 111.195 * np.cos(np.radians(0.5 * (S + N)))
    for v in (0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100):
        if v >= width_km / 5:
            return v
    return 100


def geo_format(*steps):
    nd = max([2] + [len(("%.6f" % s).rstrip("0").split(".")[1]) for s in steps])
    return "--FORMAT_GEO_MAP=DF --FORMAT_FLOAT_OUT=%%.%df --FORMAT_FLOAT_MAP=%%.%df" % (nd, nd)
