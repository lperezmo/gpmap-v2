---
title: "Encoding"
description: "get_encoding_table, genotypes_to_binary, genotypes_to_binary_packed, and the ENCODING_COLUMNS schema tuple."
---

# `gpmap.encoding`

Encoding-table construction and string-to-binary conversion. The encoding table itself is documented in detail under [Concepts: Encoding table](../concepts/encoding-table.md).

## `get_encoding_table`

```python
def get_encoding_table(
    wildtype: str,
    mutations: Mapping[int, list[str] | None],
    site_labels: list[str] | None = None,
) -> pd.DataFrame
```

Build the encoding lookup table. The result has the columns and dtypes locked by `SCHEMA.md` section 3.

### Parameters

`wildtype` (`str`, required)
:   Reference sequence.

`mutations` (`Mapping[int, list[str] | None]`, required)
:   Per-site alphabets. Keys must be exactly `0..len(wildtype) - 1`. Each non-`None` alphabet must contain the wildtype letter for that site.

`site_labels` (`list[str] | None`, default `None`)
:   Optional per-site labels. Defaults to `["0", "1", ..., "L-1"]`.

### Returns

A `pd.DataFrame` (specifically a `_LegacyAliasFrame` subclass that exposes a deprecated `genotype_index` alias for `site_index`).

### Raises

`ValueError`
:   If `wildtype` is empty, if `mutations` keys do not match `0..L-1`, or if any per-site alphabet is malformed (empty, multi-char letters, or missing the wildtype letter).

## `genotypes_to_binary_packed`

```python
def genotypes_to_binary_packed(
    genotypes: Iterable[str],
    encoding_table: pd.DataFrame,
) -> np.ndarray
```

Convert genotype strings to the packed binary representation.

### Returns

`np.ndarray[uint8]` of shape `(n_genotypes, n_bits)`. Cells are 0 or 1. Order matches the input.

### Performance

Rust-backed when the `gpmap._rust` extension is available; rayon-parallel over rows. Pure-Python fallback runs roughly 50x to 200x slower at L=16, n=65k, but is still correct.

### Raises

`UnknownLetterError`
:   If any genotype contains a letter not present in the corresponding per-site alphabet.

`SchemaError`
:   If `encoding_table` does not match the locked schema.

## `genotypes_to_binary`

```python
def genotypes_to_binary(
    genotypes: Iterable[str],
    encoding_table: pd.DataFrame,
) -> np.ndarray
```

Returns a 1-D `np.ndarray[object]` of binary strings (`"0010..."`). Calls `genotypes_to_binary_packed` under the hood and then joins each row.

!!! note "Prefer the packed form"

    `genotypes_to_binary` is the v1 surface; new code should call `genotypes_to_binary_packed` to skip the string-construction step. The string form is retained only for back-compat with downstream consumers that hash on the binary string.

## `ENCODING_COLUMNS`

```python
ENCODING_COLUMNS: tuple[str, ...] = (
    "site_index",
    "site_label",
    "wildtype_letter",
    "mutation_letter",
    "mutation_index",
    "binary_repr",
    "binary_index_start",
    "binary_index_stop",
)
```

Frozen tuple of required encoding-table column names, in canonical order. Use this if you are building an encoding table by hand and want to align column order with the schema.
