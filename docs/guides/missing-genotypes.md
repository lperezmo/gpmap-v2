---
title: "Missing genotypes"
description: "Find all unobserved genotypes in the full Cartesian product. Respects the SpaceTooLargeError safety guard."
---

# Missing genotypes

For a partial measurement, you often want the list of genotypes in the full Cartesian product of per-site alphabets that are not in your observed set. `gpm.get_missing_genotypes()` returns that complement.

![A genotype graph with observed nodes filled and missing genotypes drawn hollow](../assets/missing-genotypes-light.png#only-light)
![A genotype graph with observed nodes filled and missing genotypes drawn hollow](../assets/missing-genotypes-dark.png#only-dark)

The hollow nodes are the genotypes `get_missing_genotypes()` returns: present in the full combinatorial space but absent from the measured set.

```python
from gpmap import GenotypePhenotypeMap

gpm = GenotypePhenotypeMap(
    wildtype="AA",
    genotypes=["AA", "AT", "TA"],
    phenotypes=[0.1, 0.2, 0.3],
    mutations={0: ["A", "T"], 1: ["A", "T"]},
)

gpm.get_missing_genotypes()  # array(['TT'], dtype=object)
```

The method enumerates the full space defined by `mutations` and removes the observed `genotypes`. Order matches `enumerate_genotypes_str`: WT-prefixed alphabetical.

## Size-guarded

The full enumeration goes through `enumerate_genotypes_str`, which respects the `2**28`-row safety cap. For a 20-residue amino-acid space (20^20 = 10^26), this raises `SpaceTooLargeError`:

```python
from gpmap import SpaceTooLargeError

try:
    gpm.get_missing_genotypes()
except SpaceTooLargeError as e:
    print(e)
```

If you actually do want a huge enumeration, build the full-space genotypes manually with `allow_huge=True`:

```python
from gpmap import enumerate_genotypes_str

full = enumerate_genotypes_str(
    wildtype=gpm.wildtype,
    mutations=gpm.mutations,
    allow_huge=True,
)
missing = set(full) - set(gpm.genotypes.tolist())
```

## When to use it

The typical use is to write the complement to a CSV so a wet-lab collaborator can target the unobserved genotypes:

```python
import pandas as pd

missing = gpm.get_missing_genotypes()
pd.DataFrame({"genotype": missing}).to_csv("to_measure.csv", index=False)
```

Or as a sanity check on coverage:

```python
n_observed = gpm.n
n_missing = len(gpm.get_missing_genotypes())
print(f"{n_observed} measured, {n_missing} missing ({n_observed / (n_observed + n_missing):.1%} coverage)")
```
