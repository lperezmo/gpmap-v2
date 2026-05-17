---
title: "Loading and saving"
description: "JSON, CSV, Excel, and pickle round-trips for GenotypePhenotypeMap. When to use each format."
---

# Loading and saving

`gpmap-v2` supports four file formats out of the box: JSON, CSV (with sidecar metadata), Excel, and pickle. Each round-trips losslessly: load a saved map and it is byte-for-byte equivalent to the one you stored.

## Picking a format

=== "JSON"

    Self-contained, human-readable, schema-versioned. The default choice for sharing maps.

    ```python
    from gpmap import to_json, read_json

    to_json(gpm, "map.json")
    gpm2 = read_json("map.json")
    ```

=== "CSV"

    Pandas-friendly. Stores the data columns in CSV and the wildtype, alphabets, and site labels in a sidecar `<basename>.meta.json`. Keep both files together.

    ```python
    from gpmap import to_csv, read_csv

    to_csv(gpm, "map.csv")        # writes map.csv + map.csv.meta.json
    gpm2 = read_csv("map.csv")
    ```

=== "Excel"

    Two-sheet workbook: `data` for the table, `meta` for the schema. Good for handing a map to a non-programmer collaborator.

    ```python
    from gpmap import to_excel, read_excel

    to_excel(gpm, "map.xlsx")
    gpm2 = read_excel("map.xlsx")
    ```

=== "Pickle"

    Python-native, fastest, version-fragile. Use it only for caches you control.

    ```python
    from gpmap import to_pickle, read_pickle

    to_pickle(gpm, "map.pkl")
    gpm2 = read_pickle("map.pkl")
    ```

## What each format stores

| Format | Genotypes | Phenotypes | Stdevs | Replicates | Wildtype | Alphabets | Site labels |
|---|---|---|---|---|---|---|---|
| JSON | yes | yes | yes | yes | yes | yes | yes |
| CSV + sidecar | yes (CSV) | yes (CSV) | yes (CSV) | yes (CSV) | yes (sidecar) | yes (sidecar) | yes (sidecar) |
| Excel | yes (data) | yes (data) | yes (data) | yes (data) | yes (meta) | yes (meta) | yes (meta) |
| Pickle | full container, including caches |

## JSON schema version

JSON files written by `to_json` carry a top-level `"schema_version": "1"`. Legacy v1 files without this field still load, but issue a `UserWarning`:

```python
>>> read_json("v1_legacy_map.json")
UserWarning: JSON file has no schema_version; treating as v1 legacy format
```

If you want to migrate a v1 file to the v2 format, just re-save it after loading:

```python
from gpmap import read_json, to_json

gpm = read_json("v1_legacy_map.json")
to_json(gpm, "v2_map.json")
```

## CSV sidecar layout

`to_csv(gpm, "map.csv")` writes two files:

```
map.csv               # genotypes, phenotypes, stdeviations, n_replicates
map.csv.meta.json     # {schema_version, wildtype, mutations, site_labels}
```

`read_csv` will look for the sidecar at either `map.csv.meta.json` (current convention) or `map.meta.json` (older convention) and raises `FileNotFoundError` if neither exists.

!!! warning "Keep CSV and sidecar together"

    The CSV file alone is not enough to reconstruct the map: wildtype and per-site alphabets live in the sidecar. If you only have the CSV, you can still load it via `pd.read_csv` and rebuild manually with `GenotypePhenotypeMap.from_dataframe(df, wildtype="...", mutations=...)`.

## Excel layout

`to_excel` writes a workbook with two sheets:

| Sheet | Columns |
|---|---|
| `data` | genotypes, phenotypes, stdeviations, n_replicates |
| `meta` | key, value (rows for schema_version, wildtype, mutations, site_labels) |

`read_excel` expects this layout. JSON-encoded `mutations` and `site_labels` are stored as strings in the `meta` sheet so Excel does not mangle them.

## Pickle compatibility

Pickle preserves the full Python object, including cached binary arrays. This is the fastest format and also the most fragile:

- A pickle written with one version of `gpmap-v2` is only guaranteed to load on the same major version.
- Class paths inside the `gpmap` package are stable within a major version. Internal refactors that move classes between modules will break pickles across minor versions; we treat that as a breaking change.

For long-term storage, prefer JSON or CSV.

## DataFrame round-trip

If you already have a `pd.DataFrame`, hydrate directly:

```python
gpm = GenotypePhenotypeMap.from_dataframe(df, wildtype="AAA")
```

This is also how to round-trip through any format pandas supports natively (parquet, feather, HDF, ...). You provide the wildtype on the way in; everything else (alphabets, site labels) is inferred or defaulted.

When you save a DataFrame externally, store the wildtype and per-site alphabets in your own sidecar; `gpmap-v2` only owns the JSON/CSV/Excel sidecar conventions.

## Performance notes

| Format | n=65k (L=16) round-trip | Notes |
|---|---|---|
| Pickle | ~50 ms | full object dump |
| JSON | ~150 ms | one large indented dictionary |
| CSV | ~80 ms | data only, sidecar is trivial |
| Excel | ~1 s | openpyxl overhead dominates |

For programmatic pipelines, CSV is the sweet spot: pandas-friendly, version-stable, and fast. JSON is the best human-readable archival format.
