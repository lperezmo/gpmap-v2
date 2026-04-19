"""Tests for gpmap.errors — upper_transform / lower_transform fixed from v1."""

from __future__ import annotations

import numpy as np
from gpmap import lower_transform, upper_transform
from gpmap.errors import StandardDeviationMap, StandardErrorMap


def test_upper_transform_linear() -> None:
    y = np.array([1.0, 2.0, 3.0])
    u = np.array([1.5, 2.1, 5.0])
    out = upper_transform(y, u)
    np.testing.assert_allclose(out, [0.5, 0.1, 2.0])


def test_lower_transform_linear_differs_from_upper() -> None:
    y = np.array([1.0, 2.0, 3.0])
    u = np.array([5.0, 10.0, 30.0])
    lo = np.array([0.5, 1.9, 1.0])
    up_out = upper_transform(y, u)
    lo_out = lower_transform(y, lo)
    assert not np.allclose(up_out, lo_out)


def test_lower_transform_linear_values() -> None:
    y = np.array([1.0, 2.0, 3.0])
    lo = np.array([0.7, 1.8, 2.5])
    out = lower_transform(y, lo)
    np.testing.assert_allclose(out, [0.3, 0.2, 0.5])


def test_upper_transform_log() -> None:
    y = np.array([1.0, 10.0, 100.0])
    u = np.array([10.0, 100.0, 1000.0])
    out = upper_transform(y, u, log_transform=True, logbase=10.0)
    np.testing.assert_allclose(out, [1.0, 1.0, 1.0], atol=1e-9)


def test_upper_transform_log_clamps_at_zero() -> None:
    y = np.array([10.0])
    u = np.array([1.0])  # upper less than y, nonsensical, clipped to 0.
    out = upper_transform(y, u, log_transform=True)
    assert out[0] == 0.0


def test_stdeviation_map(tiny_binary_gpm) -> None:
    m = StandardDeviationMap(tiny_binary_gpm)
    np.testing.assert_allclose(m.upper, tiny_binary_gpm.stdeviations)
    np.testing.assert_allclose(m.lower, tiny_binary_gpm.stdeviations)


def test_standard_error_map_divides_by_sqrt_n() -> None:
    import gpmap

    gpm = gpmap.GenotypePhenotypeMap(
        wildtype="A",
        genotypes=["A"],
        phenotypes=[0.0],
        stdeviations=[1.0],
        n_replicates=[4],
    )
    m = StandardErrorMap(gpm)
    np.testing.assert_allclose(m.upper, [0.5])
    np.testing.assert_allclose(m.lower, [0.5])


def test_transforms_accept_scalars_and_lists() -> None:
    out = upper_transform([1.0, 2.0], [1.5, 2.5])
    np.testing.assert_allclose(out, [0.5, 0.5])
