"""Subsample a GenotypePhenotypeMap.

v1 returned a (float, GPM) tuple. v2 returns a NamedTuple so callers index by name.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

from ..core import GenotypePhenotypeMap


class MaskedGPM(NamedTuple):
    fraction: float
    gpm: GenotypePhenotypeMap


def mask(
    gpm: GenotypePhenotypeMap,
    fraction: float,
    *,
    rng: np.random.Generator | None = None,
) -> MaskedGPM:
    """Keep `fraction` of genotypes uniformly at random."""
    if not 0.0 < fraction <= 1.0:
        raise ValueError("fraction must be in (0, 1]")
    r = rng if rng is not None else np.random.default_rng()
    n = gpm.n
    keep_count = max(1, round(fraction * n))
    chosen = np.sort(r.choice(n, size=keep_count, replace=False))
    new_gpm = GenotypePhenotypeMap(
        wildtype=gpm.wildtype,
        genotypes=gpm.genotypes[chosen].tolist(),
        phenotypes=gpm.phenotypes[chosen],
        stdeviations=gpm.stdeviations[chosen],
        n_replicates=gpm.n_replicates[chosen],
        mutations=gpm.mutations,
        site_labels=gpm.site_labels,
    )
    true_fraction = float(keep_count) / n
    return MaskedGPM(fraction=true_fraction, gpm=new_gpm)
