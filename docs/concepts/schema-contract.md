---
title: "Schema contract"
description: "Locked surface that downstream consumers depend on. Versioning policy and migration guide from v1."
---

# Schema contract

`gpmap-v2` is a load-bearing dependency for `epistasis-v2` and `gpgraph-v2`. The contract between producer and consumers lives in [`SCHEMA.md`](https://github.com/lperezmo/gpmap-v2/blob/main/SCHEMA.md) at the repo root and is versioned with the package: breaking changes bump the major version, additive changes bump the minor.

## What is locked

The schema document pins:

1. The public surface of `GenotypePhenotypeMap`: attribute names, return types, lifecycle.
2. The dtype contract for every attribute (`float64` for phenotypes, `uint8` 2D for `binary_packed`, etc.).
3. The `encoding_table` column names, dtypes, and nullability.
4. The signatures of `genotypes_to_binary` and `genotypes_to_binary_packed`.
5. Error-bar transforms (`upper_transform`, `lower_transform`).
6. The pandas `data` view shape.
7. Serialization formats (JSON `schema_version`, CSV sidecar layout, pickle compatibility window).
8. Size-guard behavior (`SpaceTooLargeError` on `enumerate_genotypes` and friends).
9. Container invariants enforced on every construction.

Anything not in `SCHEMA.md` is internal and may change without notice. If you find yourself reaching past the schema surface, file an issue: that is a sign the schema needs an additive change.

## What downstream consumers rely on

The minimal set of imports from `gpmap` that `epistasis-v2` and `gpgraph-v2` depend on:

- `GenotypePhenotypeMap`, `GenotypePhenotypeMap.from_dataframe`
- `get_encoding_table`, `genotypes_to_binary`, `genotypes_to_binary_packed`
- `upper_transform`, `lower_transform`
- `StandardDeviationMap`, `StandardErrorMap`
- `SpaceTooLargeError`, `SchemaError`, `UnknownLetterError`
- `read_csv`, `read_json`, `read_pickle`, `read_excel` (and the `to_*` counterparts)
- All simulators under `gpmap.simulate`

These are the load-bearing exports.

## Migration from v1

`gpmap-v2` is not wire-compatible with `harmslab/gpmap`. The deltas that matter for code that consumed v1:

### Distribution

- The PyPI distribution is now `gpmap-v2`. The import path stays `gpmap`.
- Python 3.10+ required.

### Encoding table

- The column `genotype_index` is renamed to `site_index` (the v1 name was a misnomer; it was always a site index). Alias is live with a `DeprecationWarning` for one minor version.
- The new sibling `binary_packed` is exposed alongside the string-form `binary`. Prefer the packed form for any hot-path consumer.

### Constructors

- `GenotypePhenotypeMap.read_dataframe` is renamed to `from_dataframe`.

### I/O

- JSON files must carry `"schema_version": "1"`. Legacy files are readable with a `UserWarning`.

### Error transforms

- `upper_transform` and `lower_transform` now do different things. v1 had a copy-paste bug where they were identical. `lower_transform` is now genuinely the lower-bound distance.

### Stats

- `stats.unbiased_var` now honors the `axis` kwarg. v1 ignored it and hardcoded `axis=1`.

### Simulators

- `simulate.random_mutation_set` no longer mutates the module-level amino-acid list (v1 shuffled it in place).
- `simulate.MultiPeakMountFujiSimulation` peak search has a retry cap; it raises `RuntimeError` instead of spinning forever on infeasible constraints.

### Size guards

- Cartesian-product enumeration (`enumerate_genotypes_int`, `enumerate_genotypes_str`, and anything that materializes the full space like `get_missing_genotypes`) refuses to allocate beyond `2**28` rows by default. Pass `allow_huge=True` to override. v1 would silently attempt the allocation.

See [`CHANGELOG.md`](https://github.com/lperezmo/gpmap-v2/blob/main/CHANGELOG.md) for the line-by-line list.

## Versioning policy

The schema document moves in lockstep with the package version:

- A `feat: ...` commit that breaks the schema bumps the major version.
- A `feat: ...` commit that only adds bumps the minor version.
- A `fix: ...` commit that does not touch the schema bumps the patch.

`python-semantic-release` handles the version bump, changelog entry, GitHub release, wheel build, and PyPI upload on every merge to `main`.

## Size guards

Size limits are enforced by `gpmap.SpaceTooLargeError` (subclass of `ValueError`). The default cap on full-space enumeration is `2**28` rows (about 268 million), which fits comfortably in a few GB of `uint8` storage. Callers who really do want a larger space pass `allow_huge=True`:

```python
from gpmap import enumerate_genotypes_int

ints = enumerate_genotypes_int([4, 4, 4, 4, 4], max_genotypes=2**30, allow_huge=True)
```

The guard catches the most common v1 footgun: an unintentional 10^26 allocation from a 30-residue amino-acid space.
