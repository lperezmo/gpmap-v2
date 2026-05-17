---
title: "Install gpmap-v2"
description: "Install gpmap-v2 from PyPI with uv or pip. Build from source with maturin if you need an unreleased commit."
---

# Install gpmap-v2

gpmap-v2 ships prebuilt wheels for Linux (x86_64, aarch64), macOS (x86_64, aarch64), and Windows (x64). Python 3.10 or newer is required.

## Stable release

=== "uv"

    ```bash
    uv add gpmap-v2
    ```

=== "pip"

    ```bash
    pip install gpmap-v2
    ```

The package installs the import path `gpmap`, regardless of the distribution name being `gpmap-v2`.

```python
import gpmap
print(gpmap.__version__)
```

!!! info "Why two names"

    The PyPI distribution is published as `gpmap-v2` so it does not collide with the legacy `gpmap` package. The import path stays `gpmap` so downstream code that already wrote `from gpmap import GenotypePhenotypeMap` keeps working.

## Optional extras

| Extra | What it adds |
|---|---|
| `gpmap-v2[parquet]` | `pyarrow` for the optional parquet I/O codec |
| `gpmap-v2[dev]` | `pytest`, `pytest-benchmark`, `mypy`, `ruff`, `maturin` |

## Build from source

You need a Rust toolchain (stable, `>= 1.70`) and `maturin`. The recommended workflow uses `uv`:

```bash
git clone https://github.com/lperezmo/gpmap-v2
cd gpmap-v2
uv sync
uv run maturin develop --release
uv run pytest
```

After editing Rust under `crates/`:

```bash
uv run maturin develop --release && uv run pytest
```

!!! warning "Editable installs without the Rust extension"

    If `maturin develop` is skipped, the Python package still imports, but the encoding and enumeration hot paths fall back to a pure-Python implementation. You will see correct results, just measured in seconds instead of milliseconds for large maps.

## Consuming from another local project

`gpmap-v2` is designed to be consumed as an editable dependency during co-development with sister packages like [`epistasis-v2`](https://github.com/lperezmo/epistasis-v2) and [`gpgraph-v2`](https://github.com/lperezmo/gpgraph-v2). In the consumer's `pyproject.toml`:

```toml
[tool.uv.sources]
gpmap-v2 = { path = "/absolute/path/to/gpmap-v2", editable = true }

[project]
dependencies = ["gpmap-v2"]
```

Then `uv sync` or `uv add gpmap-v2` in the consumer. Imports remain `import gpmap`.

## Verify the install

```python
from gpmap import GenotypePhenotypeMap

gpm = GenotypePhenotypeMap(
    wildtype="AA",
    genotypes=["AA", "AT", "TA", "TT"],
    phenotypes=[0.0, 0.5, 0.5, 1.0],
)
print(gpm.binary_packed)
```

If the import succeeds and `binary_packed` returns a 4x2 `uint8` array, the Rust extension is wired in correctly.
