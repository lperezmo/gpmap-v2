# CHANGELOG


## v1.0.1 (2026-04-20)

### Bug Fixes

- **streamlit**: Restore missing utils.ui imports on nine pages
  ([`04b03b4`](https://github.com/lperezmo/gpmap-v2/commit/04b03b4dacb726dd1b8247b729951f1f162eb62a))

The previous refactor added stats_row() call sites across the app but the formatter that runs after
  each Edit stripped the import on its first pass (at that moment the import had no usage yet, so
  the unused-import pruner removed it). The second Edit introduced the call, but by then the import
  was already gone. Every page that renders a stats_row() row now re-imports it explicitly.

### Chores

- Add streamlit demo app under examples/streamlit
  ([`269d668`](https://github.com/lperezmo/gpmap-v2/commit/269d6685370ad270c370cb2698b449ad5475711e))

Multi-page Streamlit tour of the public API: container, encoding table, enumerate + size guard, all
  five simulators, masking, error bar transforms and unbiased stats, and JSON/CSV/pickle
  round-trips. Structured after references/streamlit-aggrid-v2. Hosted target is
  gpmap.streamlit.app; README gets the open-in-streamlit badge.

Not added as a project dependency. examples/streamlit ships its own requirements.txt for streamlit
  cloud (streamlit, plotly, gpmap-v2). uv.lock also syncs the editable project version with the
  pyproject bump to 1.0.0.

- **streamlit**: Guard letters-per-site slider against min==max
  ([`5ef2662`](https://github.com/lperezmo/gpmap-v2/commit/5ef26629bbd811c3ebc96375682899e2ca700789))

When the selected alphabet is BINARY, len(source) == 2 so the slider was rendered with min_value ==
  max_value == 2, which Streamlit rejects with a StreamlitAPIException. Branch on len(source) > 2:
  show the slider only when the alphabet actually has choices; otherwise render a static caption and
  fix alpha_size to the full alphabet size. Same fix applied in enumerate.py, which had the same
  pattern.

- **streamlit**: Replace st.subheader and st.metric with tighter UI
  ([`02228d8`](https://github.com/lperezmo/gpmap-v2/commit/02228d8cc0c093a36d5ee0d4ca7e61e113519c45))

Streamlit's default st.subheader / st.header / st.title render with heavy default sizing, and
  st.metric draws a bordered box. Both feel noisy for this kind of docs-meets-demo app. Swap for a
  small utils/ui.py module exposing a stat() helper (tiny uppercase label above a bold value) and
  stats_row() for horizontal sets, and use st.markdown("#### ...") for section headings instead of
  st.subheader. Applied across every page and the showcase entry point.

- **streamlit**: Replace use_container_width with width='stretch'
  ([`ebde7bb`](https://github.com/lperezmo/gpmap-v2/commit/ebde7bb89883fa6940074cd601f29698610ee77c))

Streamlit deprecated use_container_width; it is removed after 2025-12-31. Switch every
  st.plotly_chart call to the new width='stretch' equivalent to silence the runtime warning and stay
  ahead of the removal.

### Documentation

- Point streamlit badge at gpmap-v2.streamlit.app
  ([`5e1e939`](https://github.com/lperezmo/gpmap-v2/commit/5e1e9393f5f0dfbcb6d1f8740a2fbc2dce6a84b7))


## v1.0.0 (2026-04-19)

### Code Style

- Apply ruff format to match CI check
  ([`f04987e`](https://github.com/lperezmo/gpmap-v2/commit/f04987e7de83c73c71f5703bad07eaa8f17ff252))

### Features

- Mark 1.0.0 as the clean-break release baseline
  ([`cbd854f`](https://github.com/lperezmo/gpmap-v2/commit/cbd854f2bcc738a4b09b9d26370816513d3787c7))

BREAKING CHANGE: gpmap-v2 is a full rewrite of harmslab/gpmap with no wire-compat and no API-compat.
  The preceding 0.0.1 release was an artifact of python-semantic-release's default initial version;
  1.0.0 is the intended baseline for the stable v2 surface documented in SCHEMA.md. The
  encoding_table schema, GenotypePhenotypeMap public attributes, I/O formats, and simulator APIs are
  locked from this release forward; breaking changes bump the major version.

### Breaking Changes

- Gpmap-v2 is a full rewrite of harmslab/gpmap with no wire-compat and no API-compat. The preceding
  0.0.1 release was an artifact of python-semantic-release's default initial version; 1.0.0 is the
  intended baseline for the stable v2 surface documented in SCHEMA.md. The encoding_table schema,
  GenotypePhenotypeMap public attributes, I/O formats, and simulator APIs are locked from this
  release forward; breaking changes bump the major version.


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
