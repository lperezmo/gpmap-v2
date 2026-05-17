---
title: "I/O"
description: "CSV, JSON, Excel, and pickle round-trip functions for GenotypePhenotypeMap."
---

# `gpmap.io`

CSV, JSON, Excel, and pickle round-trip functions. See [Loading and saving](../guides/loading-saving.md) for a guide-style walkthrough.

## `read_json` and `to_json`

```python
def read_json(path: str | Path) -> GenotypePhenotypeMap

def to_json(gpm: GenotypePhenotypeMap, path: str | Path) -> None
```

JSON writes carry a top-level `"schema_version": "1"`. Reads of files without that field emit a `UserWarning`.

```python
from gpmap import read_json, to_json

to_json(gpm, "map.json")
gpm2 = read_json("map.json")
```

## `read_csv` and `to_csv`

```python
def read_csv(path: str | Path, **kwargs: Any) -> GenotypePhenotypeMap

def to_csv(gpm: GenotypePhenotypeMap, path: str | Path, **kwargs: Any) -> None
```

`to_csv` writes two files: the data CSV (genotypes, phenotypes, stdeviations, n_replicates) and a sidecar `<basename>.csv.meta.json` (wildtype, mutations, site_labels, schema_version).

`read_csv` looks for the sidecar at either `<basename>.csv.meta.json` or the older `<basename>.meta.json` and raises `FileNotFoundError` if neither exists.

`**kwargs` are forwarded to `pd.read_csv` / `pd.to_csv`. `index=False` is set as a default for `to_csv`.

## `read_excel` and `to_excel`

```python
def read_excel(path: str | Path, **kwargs: Any) -> GenotypePhenotypeMap

def to_excel(gpm: GenotypePhenotypeMap, path: str | Path, **kwargs: Any) -> None
```

Writes a two-sheet workbook (`data` and `meta`). The `meta` sheet stores schema fields as `key`/`value` pairs.

`**kwargs` are forwarded to `pd.read_excel(... sheet_name="data")` / `pd.ExcelWriter.

!!! note "Engine"

    Excel I/O uses the engine pandas picks by default. Install `openpyxl` if it is not already present in your environment.

## `read_pickle` and `to_pickle`

```python
def read_pickle(path: str | Path) -> GenotypePhenotypeMap

def to_pickle(gpm: GenotypePhenotypeMap, path: str | Path) -> None
```

Standard `pickle.dump` and `pickle.load` under the hood. `read_pickle` checks the loaded object's type and raises `TypeError` if it is not a `GenotypePhenotypeMap`.

Pickle preserves cached binary arrays, so a re-loaded map is immediately ready to use without re-running `binary_packed` construction.

!!! warning "Cross-version pickle"

    Pickle files are only guaranteed to load on the same major version. Class paths inside `gpmap` are stable within a major version; an internal refactor that moves a class between modules constitutes a breaking change.

## File-format choice

| Format | Pros | Cons |
|---|---|---|
| JSON | Self-contained, human-readable, schema-versioned | Slowest large-map write (~150 ms at n=65k) |
| CSV + sidecar | Pandas-friendly, fast | Two files to track |
| Excel | Single file, collaborator-friendly | Slowest format overall |
| Pickle | Fastest round-trip | Version-fragile, opaque |

For sharing maps with anyone outside the `gpmap-v2` ecosystem, prefer CSV or JSON.
