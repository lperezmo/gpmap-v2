"""Type stubs for the Rust-backed accelerators exposed as gpmap._rust."""

from __future__ import annotations

import numpy as np

__version__: str

def genotypes_to_binary_packed(
    genotypes: list[str],
    n_sites: int,
    n_bits: int,
    letter_to_bits_flat: np.ndarray,
    letter_to_bits_offsets: np.ndarray,
    letter_index: list[dict[str, int]],
) -> np.ndarray: ...

def hamming_to_reference(
    genotypes_int: np.ndarray, reference: np.ndarray
) -> np.ndarray: ...

def enumerate_genotypes(
    alphabet_sizes: list[int],
    max_genotypes: int = ...,
    allow_huge: bool = ...,
) -> np.ndarray: ...
