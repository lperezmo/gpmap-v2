# CHANGELOG


## v0.0.1 (2026-04-19)

### Bug Fixes

- **ci**: Drop spurious build_command from semantic_release config
  ([`59a99d8`](https://github.com/lperezmo/gpmap-v2/commit/59a99d80d4d3c28923aa99a01a0b3fa7c9ffcdff))

The release job only needs to cut the version and tag; wheel building is handled by dedicated matrix
  jobs in release.yml. The inline build_command was running on a Docker image without Rust
  installed, which caused maturin to fail.

### Chores

- Initial scaffold of gpmap-v2
  ([`967c0a7`](https://github.com/lperezmo/gpmap-v2/commit/967c0a7af9dac74ac631ef9bfc67bedb8dd0fea3))

Clean-break rewrite of harmslab/gpmap. Typed, Rust-accelerated GenotypePhenotypeMap container and
  toy-landscape simulators.

Stack: uv + maturin + PyO3 + rayon. Python 3.10+. Pandas stays. Schema contract locked in SCHEMA.md.
  Release automation via python-semantic-release on Conventional Commits; PyPI via OIDC trusted
  publisher.

See CHANGELOG.md for the full list of behaviors ported or fixed.
