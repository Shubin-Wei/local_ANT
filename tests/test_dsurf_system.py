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


def test_mad_rejection_is_centred_and_symmetric():
    r = np.r_[np.linspace(-1, 1, 41), 10.0, -10.0].astype(np.float32)
    weights, kept = S.reject_mad(r, 3.0)
    assert weights[-2:].tolist() == [0, 0] and weights[:-2].all()
    assert kept[-2:].tolist() == [0, 0]


def test_mad_rejection_per_period_does_not_reject_the_noisier_period():
    rng = np.random.default_rng(1)
    short, long_ = rng.normal(0, 0.01, 200), rng.normal(0, 0.04, 200)  # spread proportional to the period
    r = np.r_[short, long_].astype(np.float32)
    periods = np.r_[np.full(200, 1), np.full(200, 2)]
    pooled, _ = S.reject_mad(r, 3.0)
    split, _ = S.reject_mad(r, 3.0, periods)
    assert (pooled[200:] == 0).sum() > 20  # one MAD for both rejects the long period's ordinary values
    assert (split == 0).sum() <= 6  # about 0.27 % of normal values lie beyond 3 sigma


def test_small_periods_borrow_their_nearest_neighbours():
    r = np.r_[np.zeros(30), [0.0, 1.0]].astype(np.float32) + np.r_[np.linspace(-0.1, 0.1, 30), [0, 0]].astype(
        np.float32
    )
    periods = np.r_[np.full(30, 1), [2, 2]]
    centre, scale = S.robust_centre_scale(r, periods, min_rows=20)
    assert scale[-1] == scale[0] > 0  # period 2 (2 rows) is pooled with period 1
    weights, _ = S.reject_mad(r, 3.0, periods, min_rows=20)
    assert weights[-1] == 0 and weights[-2] == 1
