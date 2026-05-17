---
title: "Genotype-phenotype maps"
description: "The conceptual model behind GenotypePhenotypeMap: wildtype, per-site alphabets, measured phenotypes, and the derived binary encoding."
---

# Genotype-phenotype maps

A genotype-phenotype map (GPM) is a table of measured organisms. Each row holds a sequence (the genotype) and one or more measured properties (the phenotype, plus optional standard deviation and replicate count). `gpmap-v2` wraps that table in a typed container that adds several derived views: per-site mutation counts, a binary encoding for high-order models, and a pandas DataFrame view for analysis.

## The five fields

A `GenotypePhenotypeMap` is defined by five inputs:

| Field | Type | Meaning |
|---|---|---|
| `wildtype` | `str` | Reference sequence of length L |
| `genotypes` | `list[str] \| np.ndarray` | Observed sequences, each of length L |
| `phenotypes` | `list[float] \| np.ndarray` | Measured phenotypes, one per genotype |
| `stdeviations` | `list[float] \| np.ndarray \| None` | Optional. Defaults to `nan` |
| `n_replicates` | `list[int] \| np.ndarray \| None` | Optional. Defaults to 1 |

Two derived fields are stored alongside:

| Field | Type | Meaning |
|---|---|---|
| `mutations` | `dict[int, list[str] \| None]` | Per-site alphabets. Defaults to inferring from observed letters |
| `site_labels` | `list[str]` | Per-site labels for plotting. Defaults to `["0", "1", ..., "L-1"]` |

## Sites and alphabets

A site is a position in the wildtype. A frozen site has alphabet `None` and only ever holds the wildtype letter. An active site has an explicit alphabet that always includes the wildtype letter:

```python
gpm = GenotypePhenotypeMap(
    wildtype="AXG",
    genotypes=["AXG", "TXG", "AXC"],
    phenotypes=[0.0, 0.5, 0.3],
    mutations={
        0: ["A", "T"],
        1: None,            # frozen
        2: ["G", "C"],
    },
)
```

When you do not pass `mutations`, gpmap-v2 infers an alphabet for every site from the union of letters seen in `genotypes`. Frozen sites are detected by having only the wildtype letter in observed data.

## Two views of the same genotype

Every genotype has three equally valid representations inside the container:

1. **String form** (`gpm.genotypes`): the human-readable sequence, e.g., `"ATG"`.
2. **Binary string form** (`gpm.binary`): a `'0'`/`'1'` string under the unary-minus-one encoding from the encoding table.
3. **Packed binary form** (`gpm.binary_packed`): the same encoding as a `uint8` 2D NumPy array.

Both binary views are lazy: they are computed the first time you ask, then cached. Downstream packages should prefer `binary_packed` because it is the closest thing to bytes-on-the-wire and feeds directly into Rust kernels.

```python
gpm.genotypes      # array(['AAA', 'AAT', ...], dtype=object)
gpm.binary         # array(['00', '01', '10', '11'], dtype=object)
gpm.binary_packed  # array([[0, 0], [0, 1], [1, 0], [1, 1]], dtype=uint8)
```

## The pandas view

`gpm.data` builds a DataFrame combining all of the above, suitable for analysis and notebook display:

| Column | dtype |
|---|---|
| `genotypes` | `string` |
| `phenotypes` | `float64` |
| `stdeviations` | `float64` |
| `n_replicates` | `Int64` |
| `binary` | `string` |
| `n_mutations` | `Int64` |

It is derived state. Mutating the DataFrame does not mutate the container. To round-trip a modified DataFrame back into a `GenotypePhenotypeMap`, use `GenotypePhenotypeMap.from_dataframe(df, wildtype=...)`.

## Invariants

The container guarantees, on every successful construction:

1. `len(wildtype) == L` and every genotype has length L.
2. Every letter in every genotype is in the per-site alphabet, or equals wildtype on frozen sites.
3. `phenotypes`, `stdeviations`, and `n_replicates` are all length n.
4. `encoding_table` has exactly the right number of rows for the alphabet sizes.
5. `binary_packed[i, :].sum()` equals the Hamming distance (in the unary-minus-one sense) from genotype i to wildtype.
6. Round-trip: `from_dataframe(gpm.data).data.equals(gpm.data)` is `True` for every map.

These six invariants are enforced by contract tests in `tests/schema/` on every release.

## What is not in the container

`gpmap-v2` deliberately does not own:

- **Models.** Linear, nonlinear, classifier, and Bayesian fits live in [`epistasis-v2`](https://github.com/lperezmo/epistasis-v2).
- **Graphs.** The NetworkX wrapping for evolutionary-path analysis lives in [`gpgraph-v2`](https://github.com/lperezmo/gpgraph-v2).
- **Plotting.** `gpmap-v2` ships no plotting code; the streamlit demo uses plotly directly.

This is intentional: the container is the load-bearing dependency, kept small.

## Next

- [Encoding table](encoding-table.md) for the contract that downstream consumers rely on.
- [Schema contract](schema-contract.md) for the load-bearing surface as of `SCHEMA.md`.
