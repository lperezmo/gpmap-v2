"""Tests for gpmap.enumerate — size guard, Rust/Python parity."""

from __future__ import annotations

import numpy as np
import pytest
from gpmap import enumerate_genotypes_int, enumerate_genotypes_str
from gpmap.exceptions import SpaceTooLargeError


def test_enumerate_int_shape_and_content() -> None:
    arr = enumerate_genotypes_int([2, 3])
    assert arr.shape == (6, 2)
    assert arr.dtype == np.uint8
    # Canonical row ordering: last site varies fastest.
    expected = np.array(
        [
            [0, 0],
            [0, 1],
            [0, 2],
            [1, 0],
            [1, 1],
            [1, 2],
        ],
        dtype=np.uint8,
    )
    np.testing.assert_array_equal(arr, expected)


def test_enumerate_str_covers_full_product() -> None:
    gs = enumerate_genotypes_str(wildtype="AA", mutations={0: ["A", "T"], 1: ["A", "C", "G"]})
    assert len(gs) == 6
    assert set(gs) == {"AA", "AC", "AG", "TA", "TC", "TG"}


def test_enumerate_str_respects_frozen_site() -> None:
    gs = enumerate_genotypes_str(wildtype="AX", mutations={0: ["A", "T"], 1: None})
    assert gs == ["AX", "TX"]


def test_size_guard_blocks_by_default() -> None:
    # 20^20 is astronomically large; the guard must catch it at the 2^28 cap.
    with pytest.raises(SpaceTooLargeError):
        enumerate_genotypes_int([20] * 20)


def test_size_guard_allow_huge_not_tested_here() -> None:
    # We do not actually run a huge enumeration; assert the error path only.
    with pytest.raises(SpaceTooLargeError):
        enumerate_genotypes_int([20] * 20, max_genotypes=1_000_000, allow_huge=False)


def test_size_guard_custom_max() -> None:
    with pytest.raises(SpaceTooLargeError):
        enumerate_genotypes_int([3, 3], max_genotypes=5)
