"""Random phenotypes simulation: draw phenotypes uniformly in [low, high]."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from .base import BaseSimulation


class RandomPhenotypesSimulation(BaseSimulation):
    def __init__(
        self,
        wildtype: str,
        mutations: Mapping[int, list[str] | None],
        *,
        low: float = 0.0,
        high: float = 1.0,
        rng: np.random.Generator | None = None,
        site_labels: list[str] | None = None,
        include_binary: bool = True,
    ) -> None:
        self._low = float(low)
        self._high = float(high)
        self._rng = rng if rng is not None else np.random.default_rng()
        super().__init__(
            wildtype=wildtype,
            mutations=mutations,
            site_labels=site_labels,
            include_binary=include_binary,
        )

    def build(self) -> None:
        self._phenotypes[:] = self._rng.uniform(self._low, self._high, size=self.n)
