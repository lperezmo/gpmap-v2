"""House of Cards: NK with K = n_bits - 1 (every site sees the full genotype)."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from .nk import NKSimulation


class HouseOfCardsSimulation(NKSimulation):
    def __init__(
        self,
        wildtype: str,
        mutations: Mapping[int, list[str] | None],
        *,
        rng: np.random.Generator | None = None,
        site_labels: list[str] | None = None,
        include_binary: bool = True,
    ) -> None:
        # Compute n_bits from mutations up-front to set K. Mirrors v1.
        n_bits = 0
        for i in range(len(wildtype)):
            alpha = mutations[i]
            if alpha is None:
                continue
            n_bits += max(len(alpha) - 1, 0)
        K = max(n_bits - 1, 0)
        super().__init__(
            wildtype=wildtype,
            mutations=mutations,
            K=K,
            rng=rng,
            site_labels=site_labels,
            include_binary=include_binary,
        )
