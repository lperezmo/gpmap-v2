---
title: "The encoding table"
description: "Schema-locked DataFrame describing how each (site, mutation letter) maps to a binary slot. The contract that epistasis-v2 reads to build design matrices."
---

# The encoding table

The encoding table is a pandas DataFrame, one row per `(site, mutation_letter)` combination, that says how a string genotype maps to its packed binary form. It is the load-bearing contract between `gpmap-v2` and any package that builds models on top of it.

## Schema

| Column | dtype | Nullable | Meaning |
|---|---|---|---|
| `site_index` | `Int64` | no | 0..L-1, position in the wildtype |
| `site_label` | `string` | no | user-visible label, defaults to `str(site_index)` |
| `wildtype_letter` | `string` | no | single character, the WT at `site_index` |
| `mutation_letter` | `string` | yes (NaN on frozen sites) | the letter this row describes |
| `mutation_index` | `Int64` | yes (NaN on WT rows and frozen sites) | global mutation index starting at 1 |
| `binary_repr` | `string` | no (empty on frozen sites) | unary-minus-one encoding, length = `alphabet_size - 1` |
| `binary_index_start` | `Int64` | no | left edge of this site's slot in the concatenated binary |
| `binary_index_stop` | `Int64` | no | right edge (exclusive) |

## Worked example

For `wildtype="AX"` and `mutations={0: ["A","C","G"], 1: None}` (site 1 frozen):

| site_index | site_label | wildtype_letter | mutation_letter | mutation_index | binary_repr | binary_index_start | binary_index_stop |
|---|---|---|---|---|---|---|---|
| 0 | "0" | "A" | "A" | `<NA>` | "00" | 0 | 2 |
| 0 | "0" | "A" | "C" | 1 | "10" | 0 | 2 |
| 0 | "0" | "A" | "G" | 2 | "01" | 0 | 2 |
| 1 | "1" | "X" | `<NA>` | `<NA>` | "" | 2 | 2 |

Site 0 has alphabet size 3, so `n_bits = 2`. Site 1 is frozen and contributes 0 bits. The total `n_bits` for the map is 2.

## Unary-minus-one encoding

Each active site contributes `alphabet_size - 1` bits. The wildtype letter is the all-zero string. Each non-WT letter sets exactly one bit:

| Letter at a site with alphabet `["A", "C", "G"]` (WT=A) | binary_repr |
|---|---|
| A | "00" |
| C | "10" |
| G | "01" |

The full binary representation of a genotype is the concatenation of the per-site `binary_repr` strings, in `site_index` order. `binary_index_start` and `binary_index_stop` give you the slice indices into that concatenation.

## Reading from the encoding table

```python
gpm.encoding_table[
    ["site_index", "site_label", "mutation_letter", "mutation_index"]
].dropna()
```

This is the canonical query pattern for downstream consumers like `epistasis-v2`: drop frozen-site rows and WT rows, then group by `site_index` to build the per-site mutation index used in the design matrix.

## Legacy alias

The v1 schema had a column named `genotype_index` that was actually a site index. The v2 schema renames it to `site_index`. The old name is still available as a read-only alias for one minor version:

```python
gpm.encoding_table["genotype_index"]  # DeprecationWarning, returns site_index
```

Writes only go to `site_index`. The alias issues a `DeprecationWarning` and will be removed in a future minor release.

!!! warning "Update your consumers"

    If you maintain a package that reads `gpmap` encoding tables, switch to `site_index` now. The alias is a transition aid, not a permanent feature.

## Building one manually

For most uses you do not need to build the encoding table yourself; it is computed lazily by `GenotypePhenotypeMap` when you first access `gpm.encoding_table`. If you need it as a standalone artifact:

```python
from gpmap import get_encoding_table

table = get_encoding_table(
    wildtype="AXG",
    mutations={0: ["A", "T"], 1: None, 2: ["G", "C"]},
)
```

## Validating an externally-built table

```python
from gpmap import validate_encoding_table

validate_encoding_table(table)  # raises SchemaError if anything is off
```

`validate_encoding_table` checks that all required columns are present, that `site_index` has no NaN, and that `binary_index_stop` agrees with the implied `n_bits`.

## Going from strings to packed binary

```python
from gpmap import genotypes_to_binary_packed

packed = genotypes_to_binary_packed(["ATG", "ATC"], table)  # shape (2, n_bits), dtype uint8
```

This is the fastest path from a list of genotype strings to the packed binary representation. Unknown letters raise `UnknownLetterError`. Under the hood this dispatches to the Rust `genotypes_to_binary_packed` kernel, parallelized over rows with rayon.

The string-form sibling `genotypes_to_binary(...)` returns a NumPy object array of `'0'`/`'1'` strings. Prefer the packed form for any hot path; the string form is kept for back-compat with v1 consumers.
