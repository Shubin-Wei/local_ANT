"""Isolated PyGMT C-API renderer; no shell or external GMT command execution.

Existing layout builders supply GMT module options. Keeping GMT's native options
preserves projections/gridding while replacing shell sessions and preprocessing.
A Python worker isolates GMT's process-global session/settings from the workflow.
"""

import json
import os
import shlex
import sys
from pathlib import Path
from string import Template


def _call(session, statement, variables):
    from pygmt.helpers import build_arg_list

    statement = Template(statement).substitute(variables)
    tokens = shlex.split(statement)
    if len(tokens) < 2 or tokens[0] != "gmt" or any(t in ("|", "&&", ";", "<<<") for t in tokens):
        raise ValueError("PyGMT renderer accepts native GMT module operations only")
    module = tokens[1]
    # The C API's argv interface preserves spaces without interpreting shell
    # quotation marks as literal label text (PyGMT >=0.12 prefers argv).
    args = tokens[2:]
    if ">" in tokens:
        index = tokens.index(">")
        if index != len(tokens) - 2:
            raise ValueError("invalid GMT output redirection")
        args = tokens[2:index] + build_arg_list({}, outfile=tokens[-1])
        if module == "makecpt":
            args.append("-H")
    session.call_module(module, args)


def _operation(session, operation, variables):
    import numpy as np

    grid, palette, output = operation["grid"], operation.get("palette"), operation.get("output")
    stats = "pygmt_grid_stats.txt"
    _call(session, "gmt grdinfo %s -Cn > %s" % (shlex.quote(grid), stats), variables)
    values = np.loadtxt(stats, ndmin=1)
    low, high = float(values[4]), float(values[5])
    if not np.isfinite(low + high):
        raise ValueError("relief has no finite elevation range")
    variables.update(ZMIN=str(low), ZMAX=str(high))
    if operation["action"] == "range":
        return
    if palette == "gray":
        minimum, maximum = int(low - 2.5 * (high - low) - 1), int(high + 0.3 * (high - low) + 1)
        colors = "gray"
    else:
        minimum, maximum = low, high if high > low else low + 1
        from .gmt import ELEV_CPT

        colors = ELEV_CPT
    _call(session, "gmt makecpt -C%s -T%g/%g -Z > %s" % (colors, minimum, maximum, shlex.quote(output)), variables)
    span = max((high - low) / 2, 1e-30)
    power = 10 ** np.floor(np.log10(span))
    m = span / power
    step = (1 if m < 1.5 else 2 if m < 3.5 else 5 if m < 7.5 else 10) * power
    if high <= low:
        step = 1.0
    variables.update(ESTEP=str(step), ESTEP_MINOR=str(step / 5))


def worker(request):
    import atexit

    from pygmt.clib import Session
    from pygmt.session_management import end

    variables = {}
    with Session() as session:
        if request.get("mode") == "track":
            _call(
                session,
                "gmt grdtrack %s -G%s > %s"
                % (shlex.quote(request["points"]), shlex.quote(request["grid"]), shlex.quote(request["output"])),
                variables,
            )
            return
        for operation in request["pre"]:
            _operation(session, operation, variables) if isinstance(operation, dict) else _call(
                session, operation, variables
            )
        _call(
            session, "gmt begin lant_figure %s A+m0.25c,E%d" % (",".join(request["formats"]), request["dpi"]), variables
        )
        from .gmt import SET

        _call(session, SET.format(a=request["font"][0], l=request["font"][1]), variables)
        for statement in request["commands"]:
            _call(session, statement, variables)
        # PyGMT starts a global modern session on import. Our worker closes the
        # shared GMT session explicitly; do not close it again at interpreter exit.
        atexit.unregister(end)
        _call(session, "gmt end", variables)


def render(tmp, pre, commands, formats, dpi, font):
    import subprocess

    from ..processes import run
    from ..progress import detail

    request = Path(tmp) / "pygmt_request.json"
    request.write_text(json.dumps(dict(pre=pre, commands=commands, formats=formats, dpi=dpi, font=font)))
    env = dict(os.environ, GMT_TMPDIR=tmp, GMT_USERDIR=tmp)
    with (Path(tmp) / "pygmt.log").open("w") as log:
        result = run(
            [sys.executable, "-m", "lant.plotting.pygmt_backend", str(request)],
            cwd=tmp,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            on_poll=lambda: detail("PyGMT rendering; process isolated"),
        )
    if result.returncode:
        raise RuntimeError(
            "PyGMT failed; request and log kept at %s: %s" % (tmp, (Path(tmp) / "pygmt.log").read_text()[-3000:])
        )
    missing = [f for f in formats if not (Path(tmp) / ("lant_figure." + f)).is_file()]
    if missing:
        raise RuntimeError("PyGMT omitted formats %s; evidence kept at %s" % (missing, tmp))


def sample_grid(tmp, grid, points):
    import subprocess

    import numpy as np

    from ..processes import run

    request = Path(tmp) / "pygmt_request.json"
    output = Path(tmp) / "tracked.txt"
    request.write_text(json.dumps(dict(mode="track", grid=grid, points=points, output=str(output))))
    with (Path(tmp) / "pygmt.log").open("w") as log:
        result = run(
            [sys.executable, "-m", "lant.plotting.pygmt_backend", str(request)],
            cwd=tmp,
            env=dict(os.environ, GMT_TMPDIR=tmp, GMT_USERDIR=tmp),
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    if result.returncode:
        raise RuntimeError(
            "PyGMT grid sampling failed; evidence kept at %s: %s" % (tmp, (Path(tmp) / "pygmt.log").read_text()[-3000:])
        )
    return np.loadtxt(output, ndmin=2)[:, 2]


if __name__ == "__main__":
    worker(json.loads(Path(sys.argv[1]).read_text()))
