---
title: "GenotypePhenotypeMap"
description: "Reference for the GenotypePhenotypeMap container, validate_encoding_table, and the from_dataframe class constructor."
---

# `gpmap.core`

The core container module exposes `GenotypePhenotypeMap` and the schema validator `validate_encoding_table`.

## `GenotypePhenotypeMap`

```python
class GenotypePhenotypeMap:
    def __init__(
        self,
        wildtype: str,
        genotypes: list[str] | np.ndarray,
        phenotypes: list[float] | np.ndarray,
        *,
        stdeviations: list[float] | np.ndarray | None = None,
        n_replicates: list[int] | np.ndarray | None = None,
        mutations: Mapping[int, list[str] | None] | None = None,
        site_labels: list[str] | None = None,
        include_binary: bool = True,
    ) -> None: ...
```

### Parameters

`wildtype` (`str`, required)
:   Reference sequence. Must be a non-empty string. All `genotypes` must have this length.

`genotypes` (`list[str] | np.ndarray`, required)
:   Observed sequences, each of length `len(wildtype)`. Will be coerced to a numpy object array of strings.

`phenotypes` (`list[float] | np.ndarray`, required)
:   Measured phenotypes, length matching `genotypes`. Coerced to `float64`. If the coercion changes any values (e.g., a NumPy `Int64` overflow), a `UserWarning` is emitted once.

`stdeviations` (`list[float] | np.ndarray | None`, default `None`)
:   Per-genotype standard deviations. Defaults to `np.nan` for every genotype when omitted.

`n_replicates` (`list[int] | np.ndarray | None`, default `None`)
:   Per-genotype replicate counts. Defaults to 1 for every genotype.

`mutations` (`dict[int, list[str] | None] | None`, default `None`)
:   Per-site alphabets. Keys must be `0..L-1`. A `None` value marks a frozen site (only the wildtype letter is allowed). If omitted, the alphabet is inferred from observed letters.

`site_labels` (`list[str] | None`, default `None`)
:   Optional per-site labels. Defaults to `["0", "1", ..., "L-1"]`.

`include_binary` (`bool`, default `True`)
:   If `True`, eagerly build `binary_packed` at construction. Set to `False` to defer the cost until the first access.

### Attributes

| Attribute | Type | Notes |
|---|---|---|
| `wildtype` | `str` | Reference sequence |
| `genotypes` | `np.ndarray[object]` | Length n, dtype object of `str` |
| `phenotypes` | `np.ndarray[float64]` | Length n |
| `stdeviations` | `np.ndarray[float64]` | Length n, NaN by default |
| `n_replicates` | `np.ndarray[int64]` | Length n, 1 by default |
| `mutations` | `dict[int, list[str] \| None]` | Per-site alphabets |
| `site_labels` | `list[str]` | Length L |
| `length` | `int` | L, the wildtype length |
| `n` | `int` | n, the number of genotypes |
| `encoding_table` | `pd.DataFrame` | Lazy, schema-locked, cached |
| `binary_packed` | `np.ndarray[uint8]` | Shape (n, n_bits), lazy, cached |
| `binary` | `np.ndarray[object]` | Length n, dtype object of '0'/'1' strings, lazy, cached |
| `n_mutations` | `np.ndarray[int64]` | Per-genotype Hamming weight |
| `data` | `pd.DataFrame` | Lazy, cached, see [the data column](#data-columns) |
| `stdeviation_map` | `StandardDeviationMap` | View returning raw stdeviations |
| `standard_error_map` | `StandardErrorMap` | View returning std / sqrt(n_replicates) |

### `data` columns

The `data` property returns a DataFrame with these columns, in order:

| Column | dtype |
|---|---|
| `genotypes` | `string` |
| `phenotypes` | `float64` |
| `stdeviations` | `float64` |
| `n_replicates` | `Int64` |
| `binary` | `string` |
| `n_mutations` | `Int64` |

### `from_dataframe`

```python
@classmethod
def from_dataframe(
    cls,
    df: pd.DataFrame,
    *,
    wildtype: str | None = None,
    mutations: Mapping[int, list[str] | None] | None = None,
    site_labels: list[str] | None = None,
    include_binary: bool = True,
) -> GenotypePhenotypeMap
```

Hydrate a `GenotypePhenotypeMap` from a pandas DataFrame.

Required columns: `genotypes`, `phenotypes`. Optional: `stdeviations`, `n_replicates`. Other columns are ignored.

`wildtype` may be omitted if the DataFrame `attrs` carry a `"wildtype"` key. Otherwise it must be provided.

### `get_missing_genotypes`

```python
def get_missing_genotypes(self) -> np.ndarray
```

Returns a numpy object array of genotype strings that are in the full Cartesian product of `mutations` but not in `self.genotypes`. Respects the `SpaceTooLargeError` guard; for huge spaces, use `enumerate_genotypes_str(..., allow_huge=True)` directly.

## `validate_encoding_table`

```python
def validate_encoding_table(table: pd.DataFrame) -> None
```

Raises `SchemaError` if the table is missing required columns, has NaN in `site_index` or `binary_index_stop`, or if `binary_index_stop` disagrees with the implied total bit count. Returns `None` on success.

Use this when you build an encoding table externally (or read one from a third-party source) and want to confirm it matches the schema before feeding it into `genotypes_to_binary_packed`.
