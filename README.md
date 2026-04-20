# gpmap-v2

[![CI](https://github.com/lperezmo/gpmap-v2/actions/workflows/ci.yml/badge.svg)](https://github.com/lperezmo/gpmap-v2/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/gpmap-v2.svg)](https://pypi.org/project/gpmap-v2/)
[![Python](https://img.shields.io/pypi/pyversions/gpmap-v2.svg)](https://pypi.org/project/gpmap-v2/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://gpmap-v2.streamlit.app)

Live interactive tour: [gpmap-v2.streamlit.app](https://gpmap-v2.streamlit.app) (source under [`examples/streamlit/`](examples/streamlit/)).

A typed, Rust-accelerated container and simulator toolkit for genotype-phenotype maps.

gpmap-v2 is a clean-break rewrite of `harmslab/gpmap`. It exposes the same conceptual model (a `GenotypePhenotypeMap` backed by a pandas DataFrame) with a locked schema, vectorized hot paths, and a PyO3 Rust core for the operations that used to be inner Python loops.

## Why

- **Fast.** String-encoded genotypes are replaced with packed `uint8` matrices. The encoding step (`genotypes_to_binary`) runs rayon-parallel in Rust, delivering two orders of magnitude over the pure-Python v1 at L >= 16.
- **Typed.** Full type hints, mypy-checked, strict mode.
- **Safe.** Cartesian-product enumeration is size-guarded out of the box (`SpaceTooLargeError`). No more silent 10^26 allocations.
- **Stable surface.** The container and `encoding_table` schema are locked in [`SCHEMA.md`](SCHEMA.md) for downstream consumers.
- **Modern tooling.** `uv` + `maturin` + `pyproject.toml`. Automated releases via `python-semantic-release`. OIDC-based PyPI publishing.

## Install

```bash
pip install gpmap-v2
```

Or with uv:

```bash
uv add gpmap-v2
```

Python 3.10+. Prebuilt wheels ship for Linux (x86_64, aarch64), macOS (x86_64, aarch64), and Windows (x64).

## Quick start

```python
from gpmap import GenotypePhenotypeMap

gpm = GenotypePhenotypeMap(
    wildtype="AAA",
    genotypes=["AAA", "AAT", "ATA", "TAA", "ATT", "TAT", "TTA", "TTT"],
    phenotypes=[0.1, 0.2, 0.2, 0.6, 0.4, 0.6, 1.0, 1.1],
    stdeviations=[0.05] * 8,
)

gpm.genotypes        # np.ndarray of strings, shape (8,)
gpm.phenotypes       # np.ndarray[float64], shape (8,)
gpm.binary_packed    # np.ndarray[uint8], shape (8, 3) - the fast path
gpm.binary           # np.ndarray of '0'/'1' strings - back-compat accessor
gpm.n_mutations      # per-genotype Hamming weight
gpm.encoding_table   # pandas DataFrame per SCHEMA.md
gpm.data             # pandas DataFrame view for Jupyter
```

### Simulating landscapes

```python
from gpmap.simulate import NKSimulation, MountFujiSimulation
import numpy as np

sim = NKSimulation(
    wildtype="AAAA",
    mutations={i: ["A", "T"] for i in range(4)},
    K=2,
    rng=np.random.default_rng(0),
)
sim.phenotypes.shape  # (16,)
```

### I/O round-trips

```python
from gpmap import to_json, read_json, to_csv, read_csv

to_json(gpm, "map.json")
gpm2 = read_json("map.json")

to_csv(gpm, "map.csv")   # writes map.csv + map.csv.meta.json
gpm3 = read_csv("map.csv")
```

## Public API surface

The full list of stable exports lives in `gpmap.__all__`. The ones downstream consumers depend on (notably `epistasis-v2`) are:

- `GenotypePhenotypeMap`, `GenotypePhenotypeMap.from_dataframe`
- `get_encoding_table`, `genotypes_to_binary`, `genotypes_to_binary_packed`
- `upper_transform`, `lower_transform`
- `StandardDeviationMap`, `StandardErrorMap`
- `SpaceTooLargeError`, `SchemaError`, `UnknownLetterError`
- `read_csv`, `read_json`, `read_pickle`, `read_excel` (and `to_*` counterparts)
- All simulators under `gpmap.simulate`

The load-bearing schema contract (column names, dtypes, invariants) is in [`SCHEMA.md`](SCHEMA.md). Breaking changes to that document bump the major version.

## Development

```bash
git clone https://github.com/lperezmo/gpmap-v2
cd gpmap-v2
uv sync
uv run maturin develop --release
uv run pytest
uv run ruff check python/gpmap tests
```

A typical dev inner loop after editing Rust:

```bash
uv run maturin develop --release && uv run pytest
```

## Consuming from another local project

`gpmap-v2` is designed to be consumed as an editable dependency during co-development with sister packages (e.g. `epistasis-v2`). In the consumer's `pyproject.toml`:

```toml
[tool.uv.sources]
gpmap-v2 = { path = "/absolute/path/to/gpmap-v2", editable = true }

[project]
dependencies = ["gpmap-v2"]
```

Then `uv sync` or `uv add gpmap-v2` in the consumer. Imports remain `import gpmap`.

## Migration from v1 (`harmslab/gpmap`)

gpmap-v2 is not wire-compatible with v1. Key differences:

- Distribution name is `gpmap-v2` on PyPI; import path is still `gpmap`.
- The `encoding_table` column `genotype_index` has been renamed to `site_index` to match its actual meaning. The old name is still readable via a deprecated alias.
- `read_dataframe` is now `from_dataframe`.
- `binary_packed` (uint8 2D) is exposed alongside the string-form `binary`. Prefer the packed form for any hot-path consumer.
- JSON files must carry `"schema_version": "1"`. Legacy files are readable with a warning.
- `upper_transform` and `lower_transform` are now genuinely distinct (v1 had a copy-paste bug where they were identical).
- `stats.unbiased_var` honors the `axis` kwarg (v1 ignored it and hardcoded `axis=1`).
- `simulate.random_mutation_set` no longer mutates the module-level amino-acid list.
- `simulate.MultiPeakMountFujiSimulation` peak search has a retry cap; it raises instead of spinning forever on infeasible constraints.

See [`CHANGELOG.md`](CHANGELOG.md) for the full list.

## Releases

Releases are driven by [`python-semantic-release`](https://python-semantic-release.readthedocs.io/) on merge to `main`. Commit messages follow [Conventional Commits](https://www.conventionalcommits.org/):

- `fix: ...` -> patch
- `feat: ...` -> minor
- `feat!: ...` or `BREAKING CHANGE:` footer -> major

`CHANGELOG.md`, version bumps, Git tags, GitHub Releases, wheel builds, and PyPI uploads all happen automatically.

## License

MIT. See `LICENSE`.
