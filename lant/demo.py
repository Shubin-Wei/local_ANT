"""`lant demo`: generate a small synthetic array with a known answer and run the whole chain on it.

The demo confirms its own judgement points, which is right only because it made the data itself (a laterally
uniform medium whose model is written next to the data).  It refuses to work in a directory it did not create, so
that it can never confirm anything for real data.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

MARKER = "lant-demo.json"
TRUTH = ("true_model.txt", "true_dispersion.txt", "synthetic.yaml")
REASON = "Generated synthetic input with known truth"


def run_known_truth(directory, full=True):
    """Run a generated example to its end (or to the inversion project), confirming its gates as synthetic."""
    directory = Path(directory)
    if not all((directory / name).is_file() for name in TRUTH):
        raise RuntimeError(
            "%s is not a generated synthetic example (true model and synthetic.yaml missing)" % directory
        )
    config = str((directory / "synthetic.yaml").resolve())
    base = [sys.executable, "-m", "lant"]
    target = [] if full else ["--to", "inversion"]
    for _ in range(6):
        result = subprocess.run(base + ["run", "all", "-c", config] + target)
        if result.returncode == 0:
            return config
        if result.returncode != 3:
            raise SystemExit(result.returncode)
        state = json.loads(subprocess.check_output(base + ["status", "-c", config, "--json"]))
        gate = next((step for step, value in state["steps"].items() if value["state"] == "to confirm"), None)
        if gate is None:
            raise SystemExit("stopped without a scientific gate")
        print(
            "[lant demo] confirming %s: the data are synthetic with a known answer" % gate, file=sys.stderr, flush=True
        )
        subprocess.run(base + ["confirm", gate, "-c", config, "--reason", REASON], check=True)
    raise SystemExit("too many synthetic confirmation stops")


def run(out):
    from . import synthetic

    out = Path(out).resolve()
    if out.exists() and any(out.iterdir()) and not (out / MARKER).exists():
        raise RuntimeError("%s is not empty and was not made by lant demo; choose a new directory" % out)
    if not (out / MARKER).exists():
        print("[lant demo] writing synthetic records to %s (about a minute)" % out, file=sys.stderr, flush=True)
        synthetic.main(["--quick", "--out", str(out)])
        sys.stdout.flush()  # before the chain's own output (Python 3.8 buffers redirected output)
        (out / MARKER).write_text(json.dumps(dict(created_by="lant demo", pid=os.getpid())) + "\n")
    config = run_known_truth(out)
    print(
        "\nThe demo is complete. Open %s, and compare the model with the true one in %s."
        % (out / "output" / "REPORT.md", out / "true_model.txt")
    )
    print("For your own data: lant run all --raw RAW_DIR --stations STATION_TABLE --out OUTPUT_DIR")
    return config
