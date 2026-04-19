"""Multi-peak Mount Fuji: max over several single-peak landscapes.

v1 had an unguarded while loop in peak selection. v2 adds a retry cap.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from .base import BaseSimulation

try:
    from gpmap import _rust as _r
except ImportError:  # pragma: no cover
    _r = None


class MultiPeakMountFujiSimulation(BaseSimulation):
    def __init__(
        self,
        wildtype: str,
        mutations: Mapping[int, list[str] | None],
        *,
        peak_n: int = 2,
        min_peak_distance: int = 1,
        max_peak_distance: int | None = None,
        field_strength: float = 1.0,
        roughness_width: float = 0.0,
        roughness_dist: str = "normal",
        max_proposal_retries: int = 10_000,
        rng: np.random.Generator | None = None,
        site_labels: list[str] | None = None,
        include_binary: bool = True,
    ) -> None:
        self._peak_n = int(peak_n)
        self._min_peak_distance = int(min_peak_distance)
        self._max_peak_distance = max_peak_distance
        self._field_strength = float(field_strength)
        self._roughness_width = float(roughness_width)
        self._roughness_dist = str(roughness_dist)
        self._max_retries = int(max_proposal_retries)
        self._rng = rng if rng is not None else np.random.default_rng()
        super().__init__(
            wildtype=wildtype,
            mutations=mutations,
            site_labels=site_labels,
            include_binary=include_binary,
        )

    def _geno_int_matrix(self) -> tuple[np.ndarray, list[list[str]]]:
        per_site_alpha: list[list[str]] = []
        for i in range(self.length):
            alpha = self._mutations[i]
            per_site_alpha.append([self._wildtype[i]] if alpha is None else list(alpha))
        geno_ints = np.zeros((self.n, self.length), dtype=np.uint8)
        for row, g in enumerate(self._genotypes):
            for site, letter in enumerate(g):
                geno_ints[row, site] = per_site_alpha[site].index(letter)
        return geno_ints, per_site_alpha

    def build(self) -> None:
        geno_ints, _ = self._geno_int_matrix()
        max_d = self._max_peak_distance if self._max_peak_distance is not None else self.length

        # Pick peaks by sampling genotype indices and checking pairwise Hamming.
        peak_rows: list[int] = []
        retries = 0
        while len(peak_rows) < self._peak_n:
            if retries >= self._max_retries:
                raise RuntimeError(
                    f"failed to find {self._peak_n} peaks with "
                    f"min_dist={self._min_peak_distance}, max_dist={max_d} "
                    f"after {self._max_retries} proposals"
                )
            candidate = int(self._rng.integers(0, self.n))
            if candidate in peak_rows:
                retries += 1
                continue
            ok = True
            for existing in peak_rows:
                d = int((geno_ints[candidate] != geno_ints[existing]).sum())
                if d < self._min_peak_distance or d > max_d:
                    ok = False
                    break
            if ok:
                peak_rows.append(candidate)
            retries += 1

        # Hamming from every peak to every genotype.
        if _r is not None:
            hamming_rows = np.stack(
                [
                    np.asarray(
                        _r.hamming_to_reference(geno_ints, geno_ints[p]),
                        dtype=np.int64,
                    )
                    for p in peak_rows
                ]
            )
        else:
            hamming_rows = np.stack(
                [
                    (geno_ints != geno_ints[p][None, :]).sum(axis=1).astype(np.int64)
                    for p in peak_rows
                ]
            )

        scale = -self._field_strength * hamming_rows.astype(np.float64)

        if self._roughness_dist == "normal":
            noise = self._rng.normal(0.0, self._roughness_width, size=self.n)
        elif self._roughness_dist == "uniform":
            w = self._roughness_width
            noise = self._rng.uniform(-w, w, size=self.n)
        else:
            raise ValueError(
                f"roughness_dist must be 'normal' or 'uniform', got {self._roughness_dist!r}"
            )

        # Phenotype = max over peaks of each single-peak landscape, plus shared noise.
        self._phenotypes[:] = scale.max(axis=0) + noise
        self._peak_rows: list[int] = peak_rows

    @property
    def peak_genotypes(self) -> list[str]:
        return [self._genotypes[i] for i in getattr(self, "_peak_rows", [])]
