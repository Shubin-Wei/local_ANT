"""One DSurfTomo iteration's system: quartile rejection, DWS, smoothing rows, residual statistics."""

import numpy as np
import pytest

pytest.importorskip("numba")
from lant.inversion import dsurf_system as S  # noqa: E402


def test_percentiles_use_truncated_one_based_indices():
    r = np.arange(10, 0, -1, dtype=np.float32)  # sorted 1..10
    q25, q75 = S.percentiles(r)
    assert (q25, q75) == (2.0, 7.0)  # int(2.5) = 2 and int(7.5) = 7, 1-based
    with pytest.raises(ValueError, match="at least 4"):
        S.percentiles(np.ones(3, np.float32))


def test_rejection_uses_signed_quartiles_times_threshold():
    r = np.array([-1.0, 0.5, 1.0, 2.0, 9.0], np.float32)
    weights, kept = S.reject(r, np.float32(0.5), np.float32(2.0), np.float32(3.0))
    # The lower limit is 3 x q25 = 1.5 on the positive side: small residuals are rejected too.
    assert weights.tolist() == [0, 0, 0, 1, 0]
    assert kept.tolist() == [0, 0, 0, 2, 0]


def test_smoothing_rows_follow_main_order():
    nx, ny, nz = 5, 5, 4  # 3 x 3 x 3 parameters, one interior node
    rows, cols, values = S.smoothing(nx, ny, nz, 2.0, 10)
    counts = np.bincount(rows - 10)
    assert counts.tolist() == [1] * 13 + [7] + [1] * 13
    centre = np.flatnonzero(rows == 10 + 13)
    assert cols[centre].tolist() == [13, 12, 14, 10, 16, 4, 22]
    assert values[centre].tolist() == [12.0] + [-2.0] * 6
    assert values[0] == 4.0 and values.dtype == np.float32


def test_dws_and_statistics_in_single_precision():
    largest, average = S.dws([0, 1, 1], np.array([-1.0, 2.0, 0.5], np.float32), 4)
    assert (largest, average) == (2.5, np.float32(3.5) / 4)
    mean, std, rms = S.residual_statistics(np.array([0.1, -0.1, 0.3], np.float32))
    assert mean == np.float32(np.float32(0.3) / 3) or abs(mean - 0.1) < 1e-7
    assert abs(std - np.std([0.1, -0.1, 0.3])) < 1e-6 and abs(rms - np.sqrt(0.11 / 3)) < 1e-6
