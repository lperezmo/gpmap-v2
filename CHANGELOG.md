# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Subsequent entries are generated automatically by [`python-semantic-release`](https://python-semantic-release.readthedocs.io/) from Conventional Commit messages on merges to `main`.

## [0.1.0] - 2026-04-19

Initial release of `gpmap-v2`. Clean-break rewrite of `harmslab/gpmap`; no back-compat with v1 wire formats or packaging.

### Added

- `GenotypePhenotypeMap` container with typed numpy-backed internals (`genotypes`, `phenotypes`, `stdeviations`, `n_replicates`, `encoding_table`, `binary`, `binary_packed`, `n_mutations`, `data`).
- `GenotypePhenotypeMap.from_dataframe` classmethod (replaces v1 `read_dataframe`).
- `binary_packed`: `np.uint8` 2D representation exposed alongside the string-form `binary`. Hot-path consumers should prefer the packed form.
- `gpmap-core` Rust crate exposed as `gpmap._rust`. PyO3 bindings for `genotypes_to_binary_packed` (rayon-parallel), `hamming_to_reference`, `enumerate_genotypes`.
- `SpaceTooLargeError` size guard on all full-space enumerations (default cap `2**28`, override with `allow_huge=True`).
- Schema contract in `SCHEMA.md` at repo root. Tests in `tests/schema/` enforce the contract.
- JSON I/O carries `"schema_version": "1"`. Legacy v1 JSON is read with a `UserWarning`.
- Simulators: `BaseSimulation`, `RandomPhenotypesSimulation`, `NKSimulation`, `HouseOfCardsSimulation`, `MountFujiSimulation`, `MultiPeakMountFujiSimulation`. `mask()` returns a `MaskedGPM` NamedTuple.

### Fixed

- `stats.unbiased_var` now honors the `axis` keyword argument (v1 ignored it and hardcoded `axis=1`).
- `errors.upper_transform` and `errors.lower_transform` are now genuinely distinct functions with correct semantics. Accept `log_transform=False` kwarg.
- `simulate.random_mutation_set` no longer mutates the module-level `AMINO_ACIDS` list in place.
- `MultiPeakMountFujiSimulation` peak search now has a configurable retry cap (`max_proposal_retries`) and raises `RuntimeError` instead of spinning forever on infeasible constraints.

### Changed

- Python 3.10+ required.
- Build system: `uv` + `maturin` + `pyproject.toml`. No `setup.py`, no `Pipfile`, no `requirements.txt`, no Travis.
- `encoding_table` column `genotype_index` renamed to `site_index` (its actual meaning). The old name is still readable via a deprecated alias for one minor release.
- `.values` attribute access replaced by `.to_numpy()` internally; public properties return numpy arrays directly.
- Exceptions use specific types (`ValueError`, `NotImplementedError`, `SchemaError`, `SpaceTooLargeError`, `UnknownLetterError`) instead of bare `Exception`.

### Removed

- Dead `stats.coverage()` stub.
- `setup.py`, `MANIFEST.in`, `Pipfile`, `Pipfile.lock`, `requirements.txt`, `.travis.yml`.
- `gpmap.simulate.mask.mask` no longer returns a bare `(float, GPM)` tuple.

[0.1.0]: https://github.com/lperezmo/gpmap-v2/releases/tag/v0.1.0
