"""Genotype-space enumeration with size guards.

Wraps the Rust `enumerate_genotypes` primitive and lifts the uint8 2D result
back to string genotypes.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from .exceptions import SpaceTooLargeError

try:
    from gpmap import _rust as _r
except ImportError:  # pragma: no cover
    _r = None

DEFAULT_MAX_GENOTYPES: int = 1 << 28  # ~268M


def enumerate_genotypes_int(
    alphabet_sizes: list[int],
    *,
    max_genotypes: int = DEFAULT_MAX_GENOTYPES,
    allow_huge: bool = False,
) -> np.ndarray:
    """Emit the Cartesian product of per-site alphabet indices as a uint8 2D array."""
    if _r is not None:
        try:
            return _r.enumerate_genotypes(list(alphabet_sizes), max_genotypes, allow_huge)
        except ValueError as exc:
            if "exceeds max_genotypes" in str(exc):
                raise SpaceTooLargeError(str(exc)) from exc
            raise

    # Python fallback.
    l = len(alphabet_sizes)
    total = 1
    for s in alphabet_sizes:
        total *= s
    if not allow_huge and total > max_genotypes:
        raise SpaceTooLargeError(
            f"enumerate_genotypes would produce {total} rows, "
            f"exceeds max_genotypes = {max_genotypes}; pass allow_huge=True to override"
        )
    strides = [1] * l
    for i in range(l - 2, -1, -1):
        strides[i] = strides[i + 1] * alphabet_sizes[i + 1]
    out = np.zeros((total, l), dtype=np.uint8)
    for row in range(total):
        rem = row
        for site in range(l):
            letter = rem // strides[site]
            rem = rem % strides[site]
            out[row, site] = letter
    return out


def enumerate_genotypes_str(
    wildtype: str,
    mutations: Mapping[int, list[str] | None],
    *,
    max_genotypes: int = DEFAULT_MAX_GENOTYPES,
    allow_huge: bool = False,
) -> list[str]:
    """Enumerate all genotypes in the mutations dict as strings, WT-prefixed alphabetically."""
    L = len(wildtype)
    per_site: list[list[str]] = []
    for i in range(L):
        alpha = mutations[i]
        if alpha is None:
            per_site.append([wildtype[i]])
        else:
            per_site.append(list(alpha))
    sizes = [len(a) for a in per_site]
    ints = enumerate_genotypes_int(sizes, max_genotypes=max_genotypes, allow_huge=allow_huge)
    out: list[str] = ["".join(per_site[site][int(row[site])] for site in range(L)) for row in ints]
    return out
