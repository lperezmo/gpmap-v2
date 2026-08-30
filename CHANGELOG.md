# CHANGELOG


## v1.0.2 (2026-08-30)

### Bug Fixes

- **security**: Declare least-privilege workflow permissions
  ([`f96b824`](https://github.com/lperezmo/gpmap-v2/commit/f96b82464eefc7f277f340e88f0a6ce0145df6f7))

### Chores

- Add docs badge to README ([#2](https://github.com/lperezmo/gpmap-v2/pull/2),
  [`c13848c`](https://github.com/lperezmo/gpmap-v2/commit/c13848c1fd1ac58ad1fd815c7a079ad27767681c))

Link to the new GitHub Pages docs site at lperezmo.github.io/gpmap-v2 via the docs workflow status
  badge, between the CI and PyPI badges.

- Add Zensical docs site ([#1](https://github.com/lperezmo/gpmap-v2/pull/1),
  [`da4f401`](https://github.com/lperezmo/gpmap-v2/commit/da4f401c87014a4cb6f462ea5b37b241df0e7caf))

Stand up GitHub Pages documentation for gpmap-v2 using Zensical (modern theme, orange accents
  matching zensical.org).

Site contents - Quickstart and installation pages - Concepts: genotype-phenotype maps, encoding
  table, schema contract - Guides: loading and saving, simulators, missing genotypes - Per-module
  API reference: core, encoding, enumerate, io, simulate, stats, errors, exceptions, changelog

Configuration - zensical.toml with modern variant, light and dark palette toggle, navigation tabs
  and sections, content edit and view actions - docs/stylesheets/extra.css for the orange accent
  palette - docs/** ignored from semantic-release scope (chore: prefix)

Deployment - .github/workflows/docs.yml builds with pip-installed zensical and publishes to GitHub
  Pages via actions/deploy-pages on every push to main that touches docs/**, zensical.toml, or the
  workflow file - site/ added to .gitignore so local zensical build does not pollute the working
  tree

Pages will need to be enabled with build_type=workflow after this merges; the first run after
  enablement may need a workflow_dispatch retrigger because of the configure-pages race window.

- Bump pyo3 and numpy to 0.29 to resolve Dependabot alerts
  ([`b088265`](https://github.com/lperezmo/gpmap-v2/commit/b088265529439fe2f09a9d08b5df875fc68abfcf))

Bumps pyo3 0.22.6 -> 0.29.0 and numpy 0.22.1 -> 0.29.0 (numpy tracks the pyo3 version) to clear four
  open Dependabot alerts:

- GHSA-36hh-v3qg-5jq4 (high): out-of-bounds read in nth/nth_back for PyList and PyTuple iterators (<
  0.29.0) - GHSA-chgr-c6px-7xpp (medium): missing Sync bound on PyCFunction::new_closure closures (<
  0.29.0)

0.29.0 is the smallest version that patches all four alerts.

Required API migration for the major bump: - Python::allow_threads -> Python::detach (renamed in
  pyo3 0.27) - ndarray into_pyarray_bound -> into_pyarray (numpy dropped the _bound transition
  methods)

cargo build passes; maturin develop + pytest green (83 passed, 88% coverage). The crate does not
  call the vulnerable pyo3 code paths (no nth/nth_back on Py list/tuple iterators, no new_closure),
  so this is a dependency hardening bump, not a shipped-code fix, hence chore.

- Fix inaccurate v1 speed claim and add v1 vs v2 benchmark table
  ([`5afbf95`](https://github.com/lperezmo/gpmap-v2/commit/5afbf9598bf3af54dc6ad7973187f770d7059b32))

- Replace broken static.streamlit.io badge with shields.io
  ([`7b75aed`](https://github.com/lperezmo/gpmap-v2/commit/7b75aed844c94819035dc568aff9f7d67e976df0))

- **streamlit**: Add phenotype landscape heatmap and 3D surface to all sim pages
  ([`c94f417`](https://github.com/lperezmo/gpmap-v2/commit/c94f4178fbcdd726be0b8c26c687939eea110774))

### Documentation

- Add light/dark gallery images to docs and README
  ([#3](https://github.com/lperezmo/gpmap-v2/pull/3),
  [`91776b3`](https://github.com/lperezmo/gpmap-v2/commit/91776b37e6f6e4c7099284d202ee45a711a177e2))

Add transparent-background figures that adapt to light and dark themes: docs pages pair them with
  #only-light / #only-dark, the README uses <picture> with prefers-color-scheme (absolute raw URLs
  so they also resolve on PyPI, where <picture> degrades to the light <img>).

Images: hypercube hero, the encoding_table rendered, the genotype binary matrix, NK K=0 vs K=3,
  Mount Fuji vs House of Cards (line plots and graph views), smooth vs rugged landscape graphs, and
  a missing-genotype graph with hollow nodes. Transparent backgrounds blend into any page background
  without a visible seam.

Docs-only change; no package code touched.

- Enable MathJax rendering for inline and display math
  ([`2f9a3a5`](https://github.com/lperezmo/gpmap-v2/commit/2f9a3a51ff99b564117520aca33f06dfebecaacc))

Wire up MathJax via extra_javascript so any LaTeX in the docs renders as proper math rather than as
  code-styled text.

- Make hypercube graph node labels legible in light and dark mode
  ([#4](https://github.com/lperezmo/gpmap-v2/pull/4),
  [`153b1d0`](https://github.com/lperezmo/gpmap-v2/commit/153b1d01e1782a6ce73224b22875c17286529122))

The genotype-graph figures (Hamming hypercube hero, ruggedness graph, landscape graphs,
  missing-genotypes map) colored each node by phenotype but drew the node label and outline in a
  single per-variant ink (black in the light PNG, white in the dark PNG). That ink matched some
  fills exactly and made those labels vanish: dark text on the dark low-phenotype nodes in light
  mode, light text on the bright peaks in dark mode.

Label and outline color are now chosen per node from each node's own fill luminance, so labels stay
  readable on any page background in both PNG variants with no manual overrides. Missing (hollow)
  nodes are filled with the page background, so their labels land on the correct page ink.


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
