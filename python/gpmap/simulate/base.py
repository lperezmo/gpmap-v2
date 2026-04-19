"""BaseSimulation: a GenotypePhenotypeMap populated by a landscape generator.

v1 bug fixed: `random_mutation_set` used to shuffle `utils.AMINO_ACIDS` in place,
corrupting module state. v2 slice-copies first.
"""

from __future__ import annotations

import random
from collections.abc import Mapping

import numpy as np

from ..core import GenotypePhenotypeMap
from ..enumerate import enumerate_genotypes_str

# Canonical alphabets. Single source of truth; never mutated.
AMINO_ACIDS: tuple[str, ...] = (
    "A",
    "R",
    "N",
    "D",
    "C",
    "E",
    "Q",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
)
DNA: tuple[str, ...] = ("A", "C", "G", "T")
RNA: tuple[str, ...] = ("A", "C", "G", "U")
BINARY: tuple[str, ...] = ("0", "1")


def random_mutation_set(
    length: int,
    alphabet_size: int = 2,
    kind: str = "AA",
    *,
    rng: random.Random | None = None,
) -> dict[int, list[str]]:
    """Random per-site alphabets of a given size from the named source alphabet.

    Fixed from v1: makes a copy of the source alphabet before shuffling so the
    module-level tuple is never mutated.
    """
    if kind == "AA":
        source = list(AMINO_ACIDS)
    elif kind == "DNA":
        source = list(DNA)
    elif kind == "RNA":
        source = list(RNA)
    elif kind == "BINARY":
        source = list(BINARY)
    else:
        raise ValueError(f"unknown alphabet kind {kind!r}")
    if alphabet_size < 2 or alphabet_size > len(source):
        raise ValueError(f"alphabet_size must be in [2, {len(source)}] for kind={kind!r}")
    r = rng or random
    out: dict[int, list[str]] = {}
    for i in range(length):
        pool = source[:]
        r.shuffle(pool)
        out[i] = sorted(pool[:alphabet_size])
    return out


class BaseSimulation(GenotypePhenotypeMap):
    """Parent class for simulators. Subclasses override `build()` to fill phenotypes."""

    def __init__(
        self,
        wildtype: str,
        mutations: Mapping[int, list[str] | None],
        *,
        site_labels: list[str] | None = None,
        include_binary: bool = True,
    ) -> None:
        genotypes = enumerate_genotypes_str(wildtype, mutations)
        n = len(genotypes)
        phenotypes = np.zeros(n, dtype=np.float64)
        super().__init__(
            wildtype=wildtype,
            genotypes=genotypes,
            phenotypes=phenotypes,
            mutations=dict(mutations),
            site_labels=site_labels,
            include_binary=include_binary,
        )
        # Hook for subclasses to fill phenotypes via self._phenotypes[:] = ...
        self.build()
        # Invalidate any cached data view so it reflects freshly-built phenotypes.
        self._data = None

    def build(self) -> None:  # pragma: no cover - abstract
        raise NotImplementedError("subclasses must implement build()")

    @classmethod
    def from_length(
        cls,
        length: int,
        alphabet_size: int = 2,
        kind: str = "BINARY",
        *,
        rng: random.Random | None = None,
        **kwargs: object,
    ) -> BaseSimulation:
        mutations = random_mutation_set(length, alphabet_size, kind, rng=rng)
        # Build a wildtype from the first letter of each site alphabet.
        wildtype = "".join(mutations[i][0] for i in range(length))
        return cls(wildtype=wildtype, mutations=mutations, **kwargs)

    def set_stdeviations(self, sigma: float) -> None:
        """Set uniform stdeviations across all genotypes."""
        if sigma < 0:
            raise ValueError("sigma must be non-negative")
        self._stdeviations = np.full(self.n, float(sigma), dtype=np.float64)
        self._data = None
