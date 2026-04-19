"""Statistical utilities (v1 bugs fixed).

v1 had `unbiased_var(x, axis=None)` hardcoding `axis=1` and ignoring the kwarg.
v2 honors `axis`, uses c4-correction for bias, and drops the dead `coverage()` stub.
"""

from __future__ import annotations

import math

import numpy as np


def c4_correction(n_samples: int) -> float:
    """Bessel c4 correction factor for small-sample unbiased std estimation.

    For n <= 1 returns 1.0 (no correction possible).
    For n > 170 falls back to lgamma to avoid factorial overflow.
    """
    if n_samples <= 1:
        return 1.0
    # lgamma-based form is numerically stable for every n we care about.
    n = n_samples
    log_c4 = 0.5 * math.log(2.0 / (n - 1)) + math.lgamma(n / 2.0) - math.lgamma((n - 1) / 2.0)
    return math.exp(log_c4)


def unbiased_std(x: np.ndarray, axis: int | None = None) -> np.ndarray | float:
    """Unbiased standard deviation estimator using c4 correction.

    `std = np.std(x, ddof=1) / c4(n)`.
    """
    a = np.asarray(x, dtype=np.float64)
    n = a.size if axis is None else a.shape[axis]
    biased = np.std(a, axis=axis, ddof=1)
    c4 = c4_correction(n)
    return biased / c4


def unbiased_var(x: np.ndarray, axis: int | None = None) -> np.ndarray | float:
    """Unbiased variance estimator (c4-corrected-std, squared)."""
    s = unbiased_std(x, axis=axis)
    return np.square(s)


def unbiased_sterror(x: np.ndarray, axis: int | None = None) -> np.ndarray | float:
    """Unbiased standard error: unbiased_std / sqrt(n)."""
    a = np.asarray(x, dtype=np.float64)
    n = a.size if axis is None else a.shape[axis]
    return unbiased_std(a, axis=axis) / math.sqrt(n)


def corrected_std(biased_std: float | np.ndarray, n_samples: int) -> float | np.ndarray:
    """Given a biased std (ddof=0 style) and n, recover the unbiased estimator."""
    c4 = c4_correction(n_samples)
    # biased std (ddof=0) => multiply by sqrt(n/(n-1)) to get ddof=1, then /c4.
    if n_samples <= 1:
        return biased_std
    adj = math.sqrt(n_samples / (n_samples - 1))
    return (biased_std * adj) / c4


def corrected_sterror(biased_sterror: float | np.ndarray, n_samples: int) -> float | np.ndarray:
    return corrected_std(biased_sterror, n_samples)
