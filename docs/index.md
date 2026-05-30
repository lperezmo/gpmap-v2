---
title: "gpmap-v2"
description: "Typed, Rust-accelerated container and simulator toolkit for genotype-phenotype maps."
---

# gpmap-v2

A typed, Rust-accelerated container and simulator toolkit for genotype-phenotype maps. Clean-break rewrite of [`harmslab/gpmap`](https://github.com/harmslab/gpmap), with a locked schema, vectorized hot paths, and a PyO3 core for the operations that used to be inner Python loops.

![An L=4 biallelic genotype-phenotype map drawn as a hypercube laid out by Hamming distance, nodes colored by phenotype](assets/hypercube-light.png#only-light)
![An L=4 biallelic genotype-phenotype map drawn as a hypercube laid out by Hamming distance, nodes colored by phenotype](assets/hypercube-dark.png#only-dark)

<div class="grid cards" markdown>

-   :material-rocket-launch: **Quickstart**

    Build a `GenotypePhenotypeMap`, inspect its packed binary representation, and round-trip through JSON or CSV.

    [:octicons-arrow-right-24: Read the quickstart](quickstart.md)

-   :material-package-down: **Installation**

    Install from PyPI with `uv` or `pip`. Source builds need a Rust toolchain.

    [:octicons-arrow-right-24: Install gpmap-v2](installation.md)

-   :material-lightbulb-on: **Concepts**

    The container model, the encoding table, and the schema contract that downstream consumers depend on.

    [:octicons-arrow-right-24: Learn the model](concepts/genotype-phenotype-maps.md)

-   :material-book-open-page-variant: **Reference**

    Per-module API reference for the container, simulators, I/O, statistics, and exceptions.

    [:octicons-arrow-right-24: Browse the reference](reference/core.md)

</div>

## What you get

- **Fast.** String-encoded genotypes are replaced with packed `uint8` matrices. The encoding step (`genotypes_to_binary`) runs rayon-parallel in Rust. Construction is 5x faster than v1 at L=16 (65k genotypes), and `binary_packed` is a free cached lookup afterward, with no repeated re-encoding on each access.
- **Typed.** Full type hints, mypy-checked, strict mode.
- **Safe.** Cartesian-product enumeration is size-guarded out of the box via `SpaceTooLargeError`. No more silent 10^26 allocations.
- **Stable surface.** The container and `encoding_table` schema are locked in [`SCHEMA.md`](https://github.com/lperezmo/gpmap-v2/blob/main/SCHEMA.md) for downstream consumers.
- **Modern tooling.** `uv` plus `maturin` plus `pyproject.toml`. Automated releases via `python-semantic-release`. OIDC-based PyPI publishing.

## Live demo

A multi-page Streamlit tour is published at [gpmap-v2.streamlit.app](https://gpmap-v2.streamlit.app); the source lives under [`examples/streamlit/`](https://github.com/lperezmo/gpmap-v2/tree/main/examples/streamlit) in the repo.

## Status

Phase 1 port complete. Phase 2 Rust kernel (`gpmap._rust`) covers the dominant hot paths: encoding, enumeration, hamming reference distances. The schema in `SCHEMA.md` is the load-bearing contract that `epistasis-v2` and `gpgraph-v2` consume.

## Next steps

Pick the page that matches what you want to do:

- Just want to wire it up and play? [Quickstart](quickstart.md).
- Already familiar with v1 and want the diffs? [Migration notes](concepts/schema-contract.md#migration-from-v1).
- Building a downstream consumer that needs the binary representation? [Encoding table reference](concepts/encoding-table.md).
- Need toy landscapes? [Simulators guide](guides/simulators.md).
