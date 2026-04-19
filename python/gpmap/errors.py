"""Error-bar transforms and stdeviation/stderror wrapper maps.

Signatures are pinned by SCHEMA.md section 5. v1 had a copy-paste bug where
upper_transform and lower_transform were identical; v2 fixes that.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .core import GenotypePhenotypeMap


def upper_transform(
    bar_y: np.ndarray,
    upper: np.ndarray,
    *,
    log_transform: bool = False,
    logbase: float = 10.0,
) -> np.ndarray:
    """Matplotlib-style positive offset from bar_y up to the upper bound.

    Linear: returns ``abs(upper - bar_y)``.
    Log: returns ``log_b(upper) - log_b(bar_y)``, clipped at 0.
    """
    y = np.asarray(bar_y, dtype=np.float64)
    u = np.asarray(upper, dtype=np.float64)
    if log_transform:
        out = np.log(u) / np.log(logbase) - np.log(y) / np.log(logbase)
        return np.clip(out, 0.0, None)
    return np.abs(u - y)


def lower_transform(
    bar_y: np.ndarray,
    lower: np.ndarray,
    *,
    log_transform: bool = False,
    logbase: float = 10.0,
) -> np.ndarray:
    """Matplotlib-style positive offset from bar_y down to the lower bound.

    Linear: returns ``abs(bar_y - lower)``.
    Log: returns ``log_b(bar_y) - log_b(lower)``, clipped at 0.
    """
    y = np.asarray(bar_y, dtype=np.float64)
    lo = np.asarray(lower, dtype=np.float64)
    if log_transform:
        out = np.log(y) / np.log(logbase) - np.log(lo) / np.log(logbase)
        return np.clip(out, 0.0, None)
    return np.abs(y - lo)


class StandardDeviationMap:
    """View over a GenotypePhenotypeMap returning raw stdeviations as error bars."""

    def __init__(self, gpm: GenotypePhenotypeMap) -> None:
        self._gpm = gpm

    @property
    def upper(self) -> np.ndarray:
        return np.asarray(self._gpm.stdeviations, dtype=np.float64)

    @property
    def lower(self) -> np.ndarray:
        return np.asarray(self._gpm.stdeviations, dtype=np.float64)


class StandardErrorMap:
    """View over a GenotypePhenotypeMap returning standard errors (std / sqrt(n))."""

    def __init__(self, gpm: GenotypePhenotypeMap) -> None:
        self._gpm = gpm

    @property
    def _sterr(self) -> np.ndarray:
        std = np.asarray(self._gpm.stdeviations, dtype=np.float64)
        n = np.asarray(self._gpm.n_replicates, dtype=np.float64)
        with np.errstate(invalid="ignore", divide="ignore"):
            return std / np.sqrt(n)

    @property
    def upper(self) -> np.ndarray:
        return self._sterr

    @property
    def lower(self) -> np.ndarray:
        return self._sterr
