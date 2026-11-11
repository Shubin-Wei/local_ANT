"""Layout of the output directory (one place, so that every step and the documentation agree)."""
import os


def P(cfg: dict, *parts: str) -> str:
    if parts and parts[0] in ("01_clean", CCF, STACK, REFERENCE, DISP, INV):
        parts = ("work",) + parts
    elif parts and parts[0] == FIG:
        parts = ("results",) + parts
    return os.path.join(cfg["project"]["output_dir"], *parts)


def mk(path):
    os.makedirs(path, exist_ok=True)
    return path


INDEX = "01_index"
CCF = "02_ccf"
STACK = "03_stack"
REFERENCE = "04_reference"
DISP = "05_dispersion"
INV = "06_inversion"
FIG = "figures"


def stack_product(cfg):
    """The cross-correlation product: the stack without the removed station-days / station-hours."""
    return P(cfg, STACK, "stack.h5")


def ascii_names(names):
    import numpy as np
    n = max(len(s) for s in names) if len(names) else 1
    return np.array(names, dtype="S%d" % n)
