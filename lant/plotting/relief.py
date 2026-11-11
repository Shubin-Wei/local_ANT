"""Relief for the maps (plots.topography).

auto: SRTM 1-arc-second tiles (about 30 m; Farr et al. 2007) from the GMT data server (the server set in GMT,
GMT_DATA_SERVER, by default https://oceania.generic-mapping-tools.org), one tile per square degree, kept in
~/.cache/lant/srtm1 (tiles GMT 6.0 already holds in ~/.gmt/server/srtm1 are used as they are).  Only the names of the
tiles that cover the map are requested.  Maps wider than MAX_AUTO_DEG, beyond 60 degrees latitude (no SRTM) or with
no tile available are drawn without relief, and the reason is returned for the figure report.
A file name or a GMT remote grid (@earth_relief_...) is used as given; none (null) draws no relief.
"""
import math
import os
import shlex
import subprocess
import urllib.request

CACHE = os.path.expanduser(os.environ.get("LANT_RELIEF_CACHE", "~/.cache/lant/srtm1"))
GMT6_CACHE = os.path.expanduser("~/.gmt/server/srtm1")
MAX_AUTO_DEG = 3.0          # wider maps: no automatic relief (too many tiles); give a grid file in plots.topography
_NOTE = {}


def server():
    import tempfile
    try:
        with tempfile.TemporaryDirectory() as d:           # GMT may leave gmt.history in its working directory
            s = subprocess.run(["gmt", "get", "GMT_DATA_SERVER"], capture_output=True, text=True, timeout=30,
                               cwd=d).stdout.strip()
    except Exception:
        s = ""
    if not s.startswith("http"):
        s = "https://oceania.generic-mapping-tools.org"
    return s.rstrip("/")


def tile_names(W, E, S, N):
    out = []
    for la in range(int(math.floor(S)), int(math.floor(N - 1e-9)) + 1):
        for lo in range(int(math.floor(W)), int(math.floor(E - 1e-9)) + 1):
            out.append("%s%02d%s%03d" % ("N" if la >= 0 else "S", abs(la), "E" if lo >= 0 else "W", abs(lo)))
    return out


def tile(name):
    """Local path of one SRTM1 tile (downloaded when missing), or None when the server has none (sea) or is not
    reachable."""
    for p in (os.path.join(GMT6_CACHE, name + ".SRTMGL1.nc"), os.path.join(CACHE, name + ".SRTMGL1.jp2")):
        if os.path.exists(p) and os.path.getsize(p) > 1000:
            return p
    if os.path.exists(os.path.join(CACHE, name + ".none")):
        return None
    os.makedirs(CACHE, exist_ok=True)
    dst = os.path.join(CACHE, name + ".SRTMGL1.jp2")
    url = "%s/srtm1/%s.SRTMGL1.jp2" % (server(), name)
    try:
        with urllib.request.urlopen(url, timeout=120) as r, open(dst + ".part", "wb") as f:
            f.write(r.read())
        os.replace(dst + ".part", dst)
        return dst
    except urllib.error.HTTPError as exc:
        if exc.code == 404:                       # no tile: open sea
            open(os.path.join(CACHE, name + ".none"), "w").close()
        return None
    except Exception:
        return None


def grid(cfg, W, E, S, N, tmp):
    """-> (path of a relief grid covering W/E/S/N, note).  path None: no relief; note says why."""
    topo = cfg["plots"]["topography"]
    if not topo:
        return None, "plots.topography is none"
    if topo != "auto":
        if not topo.startswith("@") and not os.path.exists(topo):
            return None, "plots.topography %s not found" % topo
        return topo, ""
    if E - W > MAX_AUTO_DEG or N - S > MAX_AUTO_DEG:
        return None, "map wider than %g degrees: give a relief grid in plots.topography" % MAX_AUTO_DEG
    if S < -56 or N > 60:
        return None, "no SRTM relief beyond 60 N / 56 S: give a relief grid in plots.topography"
    key = "%.4f_%.4f_%.4f_%.4f" % (W, E, S, N)
    if key in _NOTE and _NOTE[key][0] and os.path.exists(_NOTE[key][0]):
        return _NOTE[key]
    paths = [p for p in (tile(n) for n in tile_names(W, E, S, N)) if p]
    if not paths:
        out = (None, "no SRTM tile could be loaded from %s (offline, or open sea): no relief" % server())
        _NOTE[key] = out
        return out
    dst = os.path.join(CACHE, "cut_%s.nc" % key)
    if not os.path.exists(dst) and not _mosaic_gdal(paths, W, E, S, N, dst, tmp) and not _mosaic_gmt(paths, W, E, S, N, dst, tmp):
        return None, "the SRTM tiles could not be read or joined (GMT needs GDAL with JPEG 2000): no relief"
    _NOTE[key] = (dst, "SRTM 1-arc-second relief (%d tile%s)" % (len(paths), "s" if len(paths) > 1 else ""))
    return _NOTE[key]


def _mosaic_gdal(paths, W, E, S, N, dst, tmp):
    """Tiles joined by a GDAL virtual mosaic and cut to the map (seamless; needs gdalbuildvrt / gdal_translate)."""
    import shutil
    if not (shutil.which("gdalbuildvrt") and shutil.which("gdal_translate")):
        return False
    vrt = os.path.join(tmp, "relief.vrt")
    q = shlex.quote
    tif = os.path.join(tmp, "relief.tif")
    cmd = ("gdalbuildvrt -q %s %s && gdal_translate -q -projwin %.6f %.6f %.6f %.6f -ot Float32 %s %s && "
           "gmt grdconvert %s=gd %s=nf" % (q(vrt), " ".join(map(q, paths)), W, N, E, S, q(vrt), q(tif), q(tif), q(dst + ".part")))
    ok = subprocess.run(cmd, shell=True, capture_output=True, cwd=tmp).returncode == 0 and os.path.exists(dst + ".part")
    if ok:
        os.replace(dst + ".part", dst)
    return ok


def _mosaic_gmt(paths, W, E, S, N, dst, tmp):
    """Without the GDAL programs: one tile cut by GMT (maps across tile edges get no relief)."""
    for p in paths:
        nm = os.path.basename(p)[:7]
        la0 = int(nm[1:3]) * (1 if nm[0] == "N" else -1)
        lo0 = int(nm[4:7]) * (1 if nm[3] == "E" else -1)
        if lo0 <= W and E <= lo0 + 1 and la0 <= S and N <= la0 + 1:
            src = p + ("=gd" if p.endswith(".jp2") else "")
            r = subprocess.run("gmt grdcut %s -R%.6f/%.6f/%.6f/%.6f -G%s" % (shlex.quote(src), W, E, S, N, shlex.quote(dst)),
                               shell=True, capture_output=True, cwd=tmp)
            return r.returncode == 0 and os.path.exists(dst)
    return False
