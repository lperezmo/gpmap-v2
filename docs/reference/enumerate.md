---
title: "Enumeration"
description: "enumerate_genotypes_int and enumerate_genotypes_str: size-guarded Cartesian product enumeration over per-site alphabets."
---

# `gpmap.enumerate`

Two functions enumerate the full Cartesian product of per-site alphabets, both with a `SpaceTooLargeError` safety cap.

## `enumerate_genotypes_int`

```python
def enumerate_genotypes_int(
    alphabet_sizes: list[int],
    *,
    max_genotypes: int = 2**28,
    allow_huge: bool = False,
) -> np.ndarray
```

Emit the Cartesian product of per-site alphabet indices as a `uint8` 2D array.

### Parameters

`alphabet_sizes` (`list[int]`, required)
:   Per-site alphabet sizes. Each entry must be >= 1.

`max_genotypes` (`int`, default `2**28`)
:   Maximum number of rows to materialize. Raises `SpaceTooLargeError` if the product exceeds this.

`allow_huge` (`bool`, default `False`)
:   If `True`, bypass the cap. Use sparingly.

### Returns

`np.ndarray[uint8]` of shape `(prod(alphabet_sizes), len(alphabet_sizes))`. Row k has `out[k, site] = letter_index` for the k-th genotype in lexicographic order.

### Raises

`SpaceTooLargeError`
:   If `prod(alphabet_sizes) > max_genotypes` and `allow_huge=False`.

### Performance

Rust-backed when `gpmap._rust` is available. Pure-Python fallback walks the indices row by row.

## `enumerate_genotypes_str`

```python
def enumerate_genotypes_str(
    wildtype: str,
    mutations: Mapping[int, list[str] | None],
    *,
    max_genotypes: int = 2**28,
    allow_huge: bool = False,
) -> list[str]
```

Enumerate the Cartesian product as genotype strings.

### Parameters

`wildtype` (`str`, required)
:   Reference sequence, used to fix the per-site letter for frozen sites (`mutations[i] = None`).

`mutations` (`Mapping[int, list[str] | None]`, required)
:   Per-site alphabets. Frozen sites (`None`) contribute only the wildtype letter.

`max_genotypes`, `allow_huge`
:   Forwarded to `enumerate_genotypes_int`.

### Returns

A `list[str]` of length `prod(alphabet_sizes)`. Order is lexicographic by per-site alphabet index, which puts the wildtype-derived ordering at the front of each site's run.

### Example

```python
from gpmap import enumerate_genotypes_str

enumerate_genotypes_str(
    wildtype="AAT",
    mutations={0: ["A", "T"], 1: None, 2: ["T", "C", "G"]},
)
# ['AAT', 'AAC', 'AAG', 'TAT', 'TAC', 'TAG']
```

### Use cases

- Building a list of unobserved genotypes (`gpm.get_missing_genotypes()` calls this under the hood).
- Generating an exhaustive panel for an in-silico screen.
- Sanity-checking that the per-site alphabet sizes match the binary representation produced by `genotypes_to_binary_packed`.

!!! warning "Memory cost"

    A full L=20 amino-acid space has `20^20 ~ 10^26` entries. The default cap (`2**28 ~ 2.68e8`) is sized to fit comfortably in a few GB of `uint8` storage. Raising it is rarely the right move; consider whether you really need every genotype or just a sample.
