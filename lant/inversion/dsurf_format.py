"""Text output of DSurfTomo's main.f90, written the way gfortran writes it.

Two edit forms occur: list-directed output (``write(u,*)``) of REAL values,
which gfortran (libgfortran write.c, checked against gfortran 9.4 on 140 000
values) writes with 9 significant digits in a 17-character item, fixed form for
0.1 <= |x| < 1e9 and exponent form otherwise; and the explicit ``Fw.d`` / ``Iw``
descriptors of the model, log and screen lines. Both round the exact binary
value to nearest (ties to even), as Python's ``%`` formatting does.
"""

import numpy as np


def list_real(value):
    """One list-directed REAL(4) item, including its leading separator blank."""
    x = float(np.float32(value))
    if x != x or x in (float("inf"), float("-inf")):
        text = "NaN" if x != x else ("Infinity" if x > 0 else "-Infinity")
        return text.rjust(17)
    mantissa, exponent = ("%.8e" % x).split("e")
    e = int(exponent)
    if x == 0 or -1 <= e < 9:
        sign = "-" if mantissa.startswith("-") else ""
        digits = mantissa.lstrip("-").replace(".", "")
        if x == 0:
            e = 0
        if e < 0:
            body = "0." + digits
        else:
            body = digits[: e + 1] + "." + digits[e + 1 :]
        return (sign + body).rjust(13) + "    "
    return ("%sE%s%02d" % (mantissa, "-" if e < 0 else "+", abs(e))).rjust(17)


def list_record(*items):
    """A list-directed record of character and REAL items (character first, or reals only)."""
    out = []
    for i, item in enumerate(items):
        if isinstance(item, str):
            out.append((" " if i == 0 else "") + item)
        else:
            out.append(list_real(item))
    return "".join(out)


def fixed(value, width, decimals):
    """Fortran Fw.d of a value (REAL or double), asterisks when it does not fit."""
    x = float(value)
    if x != x:
        return "NaN".rjust(width)
    text = "%.*f" % (decimals, x)
    if len(text) > width and decimals > 0:
        # Fortran drops the optional leading zero before the decimal point first.
        if text.startswith("0."):
            text = text[1:]
        elif text.startswith("-0."):
            text = "-" + text[2:]
    if len(text) > width:
        return "*" * width
    return text.rjust(width)


def integer(value, width):
    text = "%d" % value
    return "*" * width if len(text) > width else text.rjust(width)


def fixed_record(values, width, decimals):
    return "".join(fixed(v, width, decimals) for v in values)
