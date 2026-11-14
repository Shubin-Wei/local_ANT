"""gfortran-compatible text of DSurfTomo's main.f90 (strings checked against gfortran 9.4 output)."""

import numpy as np

from lant.inversion import dsurf_format as fmt


def test_list_directed_real_matches_gfortran():
    cases = [  # a list: 0.0 and -0.0 would be one dictionary key
        (0.0, "   0.00000000    "),
        (-0.0, "  -0.00000000    "),
        (0.15, "  0.150000006    "),
        (99.99995, "   99.9999466    "),
        (1e8, "   100000000.    "),
        (1e9, "   1.00000000E+09"),
        (0.09999999, "   9.99999866E-02"),
        (-3.5e-3, "  -3.50000011E-03"),
        (3e38, "   3.00000001E+38"),
    ]
    for value, text in cases:
        assert fmt.list_real(np.float32(value)) == text
    assert fmt.list_record("DWS:", 0.0, -0.0) == " DWS:   0.00000000      -0.00000000    "
    assert fmt.list_record(1.0, 2.0) == "   1.00000000       2.00000000    "


def test_fixed_rounds_to_nearest_even_and_overflows_with_asterisks():
    assert fmt.fixed(np.float32(0.125), 7, 2) == "   0.12"
    assert fmt.fixed(np.float32(-0.00004), 10, 5) == "  -0.00004"
    assert fmt.fixed(np.float32(123456.789), 10, 5) == "*" * 10
    assert fmt.fixed(-1e-7, 8, 1) == "    -0.0"
    assert fmt.fixed(0.5, 4, 3) == ".500"
    assert fmt.integer(149, 5) == "  149" and fmt.integer(100, 2) == "**"
    assert fmt.fixed_record((99.97554, 30.10809), 10, 5) == "  99.97554  30.10809"
