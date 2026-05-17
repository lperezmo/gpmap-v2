---
title: "Quickstart"
description: "Construct a GenotypePhenotypeMap, inspect its packed binary representation, and round-trip it through JSON or CSV."
---

# Quickstart

Build a `GenotypePhenotypeMap` from a wildtype string, a list of genotypes, and their measured phenotypes. Everything else (binary encoding, encoding table, pandas view) is derived on demand.

## Construct a map

```python
from gpmap import GenotypePhenotypeMap

gpm = GenotypePhenotypeMap(
    wildtype="AAA",
    genotypes=["AAA", "AAT", "ATA", "TAA", "ATT", "TAT", "TTA", "TTT"],
    phenotypes=[0.1, 0.2, 0.2, 0.6, 0.4, 0.6, 1.0, 1.1],
    stdeviations=[0.05] * 8,
)
```

The constructor figures the per-site alphabet out from observed letters. Pass `mutations={i: [...]}` explicitly when you want to lock in an alphabet that does not appear in every column of your measured genotypes.

## Inspect what you got back

```python
gpm.genotypes        # np.ndarray of strings, shape (8,)
gpm.phenotypes       # np.ndarray[float64], shape (8,)
gpm.binary_packed    # np.ndarray[uint8], shape (8, 3), the fast path
gpm.binary           # np.ndarray of '0'/'1' strings, back-compat accessor
gpm.n_mutations      # per-genotype Hamming weight
gpm.encoding_table   # pandas DataFrame per SCHEMA.md
gpm.data             # pandas DataFrame view for Jupyter
```

`binary_packed` is the recommended input for downstream consumers that need the encoded representation. It is computed once at construction (when `include_binary=True`, the default), then cached.

!!! tip "Inspecting in a notebook"

    `gpm.data` is the right thing to print in a notebook. It returns a `pd.DataFrame` view with `genotypes`, `phenotypes`, `stdeviations`, `n_replicates`, `binary`, and `n_mutations` columns, suitable for displaying inline.

## Round-trip through JSON

```python
from gpmap import to_json, read_json

to_json(gpm, "map.json")
gpm2 = read_json("map.json")

assert (gpm.genotypes == gpm2.genotypes).all()
assert (gpm.phenotypes == gpm2.phenotypes).all()
```

The JSON file carries `schema_version` so legacy v1 files still load with a `UserWarning`.

## Round-trip through CSV

```python
from gpmap import to_csv, read_csv

to_csv(gpm, "map.csv")          # writes map.csv + map.csv.meta.json
gpm3 = read_csv("map.csv")
```

CSV stores the data columns; the sidecar `map.csv.meta.json` carries the wildtype, per-site alphabets, and site labels. Keep both files together.

## Hydrate from a DataFrame

```python
import pandas as pd
from gpmap import GenotypePhenotypeMap

df = pd.DataFrame({
    "genotypes": ["AAA", "AAT", "ATA", "TAA"],
    "phenotypes": [0.1, 0.2, 0.2, 0.6],
})

gpm = GenotypePhenotypeMap.from_dataframe(df, wildtype="AAA")
```

`from_dataframe` is also how you load the output of `gpm.data` back into a fresh map after manipulating it externally.

## Simulate a landscape

For toy data, the `gpmap.simulate` subpackage provides a small zoo of generators:

```python
import numpy as np
from gpmap.simulate import NKSimulation

sim = NKSimulation(
    wildtype="AAAA",
    mutations={i: ["A", "T"] for i in range(4)},
    K=2,
    rng=np.random.default_rng(0),
)
sim.phenotypes.shape  # (16,)
```

Each simulator is a `GenotypePhenotypeMap` subclass, so anything you can do with `gpm` you can do with `sim`. See the [simulators guide](guides/simulators.md) for the full menagerie (Mount Fuji, multi-peak Fuji, House of Cards, random, mask).

## Where to next

- The container surface and invariants: [Genotype-phenotype maps](concepts/genotype-phenotype-maps.md).
- The encoding-table contract for downstream consumers: [Encoding table](concepts/encoding-table.md).
- File formats and round-tripping: [Loading and saving](guides/loading-saving.md).
- Toy landscapes: [Simulators](guides/simulators.md).
