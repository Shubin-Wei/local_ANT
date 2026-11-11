"""Running GMT 6 (modern mode) for the map and section figures.

Every figure is one shell script run in a temporary directory: colour tables and grids are made before
`gmt begin` (GMT 6.0 writes no file inside a session), the session writes <out>.<format> for every format of
plots.formats.  On failure the script is kept and its path reported (set LANT_KEEP_GMT=1 to keep it always).
"""

import os
import shutil
import shlex
import subprocess
import tempfile

SET = (
    "gmt set FORMAT_GEO_MAP ddd.xxF FONT_ANNOT_PRIMARY {a}p,Helvetica FONT_LABEL {l}p,Helvetica MAP_FRAME_TYPE plain "
    "MAP_FRAME_PEN 0.5p MAP_TICK_LENGTH_PRIMARY 0.07c MAP_ANNOT_OFFSET_PRIMARY 0.06c MAP_LABEL_OFFSET 0.12c PS_LINE_JOIN round"
)
# continuous colour tables of the figures (GMT colour lists): perturbation and checkerboards red (slow) - white -
# blue (fast)
PERT_CPT = "#7a0000,#e31a1c,#fd8d3c,#fee391,white,#9bdcff,#2b8cff,#0b3fbf,#081d6b"
CHECK_CPT = "red,white,blue"
ELEV_CPT = "#a8cfa0,#cfe0b0,#ece6bf,#e3cfa4,#cdb08a,#f4efe6"


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


def workdir():
    return tempfile.mkdtemp(prefix="lant_gmt_")


def run(tmp, pre, cmd, out, formats=("png", "pdf"), dpi=300, font=(7, 8)):
    if not available():
        raise GMTError(
            "GMT 6 is not installed or not on PATH (needed for the map and section figures): "
            "`conda install -c conda-forge gmt`"
        )
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    script = os.path.join(tmp, "draw.sh")
    # the figure is drawn under a plain name inside the scratch directory and moved to `out` afterwards, so that
    # an output path with spaces or shell characters never reaches the GMT command line
    q = shlex.quote
    with open(script, "w") as f:
        f.write(
            "set -e\nexport GMT_TMPDIR=%s GMT_USERDIR=%s\ncd %s\n%s\ngmt begin lant_figure %s A+m0.25c,E%d\n%s\n%s\n"
            "gmt end\n"
            % (
                q(tmp),
                q(tmp),
                q(tmp),
                "\n".join(pre),
                ",".join(formats),
                dpi,
                SET.format(a=font[0], l=font[1]),
                "\n".join(cmd),
            )
        )
    r = subprocess.run(["bash", script], capture_output=True, text=True)
    if r.returncode != 0 or "ERROR" in r.stderr:
        raise GMTError("GMT failed for %s:\n%s\n(script kept: %s)" % (out, r.stderr[-3000:], script))
    for fmt in formats:
        shutil.move(os.path.join(tmp, "lant_figure." + fmt), os.path.abspath(out) + "." + fmt)
    if not os.environ.get("LANT_KEEP_GMT"):
        shutil.rmtree(tmp, ignore_errors=True)
    return [out + "." + f for f in formats]


RELIEF_NOTES = []  # why a figure was drawn without relief (for the figure report)


def relief_cmds(cfg, tmp, W, E, S, N, palette="gray"):
    """Commands that cut the relief grid (plots.topography, lant.plotting.relief) to the map, its hill shade and a
    colour table.  Returns (commands, have_relief).  palette gray: the light grey under the model maps; elevation:
    the coloured relief of the station map."""
    from . import relief

    pad = 0.01
    topo, note = relief.grid(cfg, W - pad, E + pad, S - pad, N + pad, tmp)
    if not topo:
        if note not in RELIEF_NOTES:
            RELIEF_NOTES.append(note)
        return [], False
    pre = [
        "gmt grdcut %s -R%.5f/%.5f/%.5f/%.5f -G%s/topo.nc"
        % (shlex.quote(topo), W - pad, E + pad, S - pad, N + pad, tmp),
        "gmt grdgradient %s/topo.nc -A315 -Nt0.5 -G%s/shade.nc" % (tmp, tmp),
    ]
    if palette == "gray":
        # elevations mapped onto the pale end of the grey scale: light relief under semi-transparent values
        pre.append(
            "read ZMIN ZMAX <<< \"$(gmt grdinfo -C %s/topo.nc | awk '{print $6, $7}')\"; "
            "gmt makecpt -Cgray -T$(awk -v a=$ZMIN -v b=$ZMAX 'BEGIN{printf \"%%d/%%d\", a-2.5*(b-a)-1, b+0.3*(b-a)+1}') -Z "
            "> %s/relief.cpt" % (tmp, tmp)
        )
    else:
        pre.append(
            "read ZMIN ZMAX <<< \"$(gmt grdinfo -C %s/topo.nc | awk '{print $6, $7}')\"; "
            "gmt makecpt -C%s -T$(awk -v a=$ZMIN -v b=$ZMAX 'BEGIN{printf \"%%d/%%d\", a, (b>a)?b:a+1}') -Z > %s/relief.cpt"
            % (tmp, ELEV_CPT, tmp)
        )
    return pre, True


def mercator_height(width, W, E, S, N):
    import numpy as np

    ys = lambda la: np.log(np.tan(np.pi / 4 + np.radians(la) / 2))
    return width * (ys(N) - ys(S)) / np.radians(E - W)


def nice_step(span, n=4):
    import numpy as np

    raw = span / n
    p = 10 ** np.floor(np.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if m * p >= raw:
            return m * p
    return 10 * p


def scale_km(W, E, S, N):
    """A round scale-bar length of about a quarter of the map width (km) and the latitude it is drawn for."""
    import numpy as np

    width_km = (E - W) * 111.195 * np.cos(np.radians(0.5 * (S + N)))
    for v in (0.5, 1, 2, 5, 10, 20, 50, 100):
        if v >= width_km / 5:
            return v, 0.5 * (S + N)
    return 100, 0.5 * (S + N)
