"""NK fitness landscape simulator.

Uses vectorized numpy indexing instead of v1's triple Python loop.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from .base import BaseSimulation


class NKSimulation(BaseSimulation):
    """NK model: each site contributes a K-neighborhood-dependent fitness term."""

    def __init__(
        self,
        wildtype: str,
        mutations: Mapping[int, list[str] | None],
        *,
        K: int = 1,
        rng: np.random.Generator | None = None,
        site_labels: list[str] | None = None,
        include_binary: bool = True,
    ) -> None:
        self._K = int(K)
        self._rng = rng if rng is not None else np.random.default_rng()
        super().__init__(
            wildtype=wildtype,
            mutations=mutations,
            site_labels=site_labels,
            include_binary=include_binary,
        )

    def build(self) -> None:
        # Work on the packed binary matrix (n, n_bits).
        B = self.binary_packed
        n, n_bits = B.shape
        if n_bits == 0:
            self._phenotypes[:] = 0.0
            return
        K = min(self._K, n_bits - 1)
        # For each site i in 0..n_bits-1, the neighborhood is columns [i, i+1, ..., i+K]
        # (wrapping around). Build a random table of size 2^(K+1) per site,
        # sampled once, then index into it per row.
        table_size = 2 ** (K + 1)
        tables = self._rng.uniform(0.0, 1.0, size=(n_bits, table_size))

        # Compute neighborhood index for every (row, site) with vectorized rolls.
        # idx[:, i] = sum over j in 0..K of B[:, (i + j) % n_bits] << (K - j)
        idx = np.zeros((n, n_bits), dtype=np.int64)
        for j in range(K + 1):
            shift = K - j
            rolled = np.roll(B, -j, axis=1).astype(np.int64)
            idx += rolled << shift

        # Sum contributions across sites.
        contrib = tables[np.arange(n_bits), idx]  # shape (n, n_bits)
        self._phenotypes[:] = contrib.sum(axis=1) / n_bits
