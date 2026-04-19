"""Tests for gpmap.stats — bugs from v1 (hardcoded axis=1) are fixed."""

from __future__ import annotations

import math

import numpy as np
from gpmap.stats import (
    c4_correction,
    corrected_std,
    unbiased_std,
    unbiased_sterror,
    unbiased_var,
)


def test_c4_correction_known_values() -> None:
    assert c4_correction(1) == 1.0
    # c4(2) = sqrt(2/pi)
    assert math.isclose(c4_correction(2), math.sqrt(2 / math.pi), rel_tol=1e-12)


def test_unbiased_std_honors_axis_none() -> None:
    x = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    s = unbiased_std(x, axis=None)
    assert np.isscalar(s) or s.shape == ()


def test_unbiased_std_honors_axis_0() -> None:
    x = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    s = unbiased_std(x, axis=0)
    assert s.shape == (3,)


def test_unbiased_std_honors_axis_1() -> None:
    x = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    s = unbiased_std(x, axis=1)
    assert s.shape == (2,)


def test_unbiased_var_is_square_of_std() -> None:
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    s = unbiased_std(x)
    v = unbiased_var(x)
    assert math.isclose(v, s * s, rel_tol=1e-12)


def test_unbiased_sterror_divides_by_sqrt_n() -> None:
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    s = unbiased_std(x)
    se = unbiased_sterror(x)
    assert math.isclose(se, s / math.sqrt(5), rel_tol=1e-12)


def test_corrected_std_roundtrip() -> None:
    # Given a biased std (ddof=0), recover the c4-corrected unbiased std.
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    biased = np.std(x, ddof=0)
    recovered = corrected_std(biased, len(x))
    direct = unbiased_std(x)
    assert math.isclose(recovered, direct, rel_tol=1e-12)
