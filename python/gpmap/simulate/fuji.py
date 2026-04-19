"""Mount Fuji: single-peak fitness landscape with additive term + noise."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from .base import BaseSimulation

try:
    from gpmap import _rust as _r
except ImportError:  # pragma: no cover
    _r = None


class MountFujiSimulation(BaseSimulation):
    def __init__(
        self,
        wildtype: str,
        mutations: Mapping[int, list[str] | None],
        *,
        field_strength: float = 1.0,
        roughness_width: float = 0.0,
        roughness_dist: str = "normal",
        rng: np.random.Generator | None = None,
        site_labels: list[str] | None = None,
        include_binary: bool = True,
    ) -> None:
        self._field_strength = float(field_strength)
        self._roughness_width = float(roughness_width)
        self._roughness_dist = str(roughness_dist)
        self._rng = rng if rng is not None else np.random.default_rng()
        super().__init__(
            wildtype=wildtype,
            mutations=mutations,
            site_labels=site_labels,
            include_binary=include_binary,
        )

    def build(self) -> None:
        # Map genotypes -> per-site int codes for vectorized hamming to WT.
        per_site_alpha: list[list[str]] = []
        for i in range(self.length):
            alpha = self._mutations[i]
            per_site_alpha.append([self._wildtype[i]] if alpha is None else list(alpha))

        geno_ints = np.zeros((self.n, self.length), dtype=np.uint8)
        for row, g in enumerate(self._genotypes):
            for site, letter in enumerate(g):
                geno_ints[row, site] = per_site_alpha[site].index(letter)
        wt_ints = np.zeros(self.length, dtype=np.uint8)
        for site, letter in enumerate(self._wildtype):
            wt_ints[site] = per_site_alpha[site].index(letter)

        if _r is not None:
            hamming = np.asarray(_r.hamming_to_reference(geno_ints, wt_ints), dtype=np.int64)
        else:
            hamming = (geno_ints != wt_ints[None, :]).sum(axis=1).astype(np.int64)

        if self._roughness_dist == "normal":
            noise = self._rng.normal(0.0, self._roughness_width, size=self.n)
        elif self._roughness_dist == "uniform":
            w = self._roughness_width
            noise = self._rng.uniform(-w, w, size=self.n)
        else:
            raise ValueError(
                f"roughness_dist must be 'normal' or 'uniform', got {self._roughness_dist!r}"
            )

        self._phenotypes[:] = -self._field_strength * hamming.astype(np.float64) + noise
