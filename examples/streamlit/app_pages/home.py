"""Home page: what a genotype-phenotype map is and why gpmap-v2 exists."""

from __future__ import annotations

import gpmap
import streamlit as st

st.markdown(
    """
A **genotype-phenotype map** (GP-map) is a function from every reachable
sequence in some local space to a measured phenotype. `gpmap-v2` gives you
a typed container for that data, plus a handful of Rust-accelerated
hot paths for encoding, enumeration, and Hamming-distance work.
"""
)

c1, c2, c3 = st.columns(3)
c1.metric("Version", gpmap.__version__)
c2.metric("Schema", "locked (SCHEMA.md)")
c3.metric("Python", "3.10+")

st.subheader("What this app covers")
st.markdown(
    """
- **Core** - the `GenotypePhenotypeMap` container, the encoding table, and
  the safe enumerator with its `SpaceTooLargeError` guard.
- **Simulators** - Mount Fuji, NK, House of Cards, Multi-peak Fuji, and
  Random. Each page is interactive so you can feel the parameters.
- **Utilities** - subsampling via `mask()`, error-bar transforms, and
  JSON / CSV / pickle round-trips.

Every page ends with a `Code` expander showing the exact snippet it ran,
so you can copy it into a notebook and keep going.
"""
)

st.subheader("Quick start")
st.code(
    """from gpmap import GenotypePhenotypeMap

gpm = GenotypePhenotypeMap(
    wildtype="AAA",
    genotypes=["AAA", "AAT", "ATA", "TAA", "ATT", "TAT", "TTA", "TTT"],
    phenotypes=[0.1, 0.2, 0.2, 0.6, 0.4, 0.6, 1.0, 1.1],
    stdeviations=[0.05] * 8,
)

gpm.binary_packed   # np.ndarray[uint8], shape (8, 3)
gpm.n_mutations     # per-genotype Hamming weight
gpm.encoding_table  # pandas DataFrame per SCHEMA.md
""",
    language="python",
)

st.subheader("Links")
st.markdown(
    """
- [GitHub](https://github.com/lperezmo/gpmap-v2)
- [PyPI](https://pypi.org/project/gpmap-v2/)
- [Schema contract](https://github.com/lperezmo/gpmap-v2/blob/main/SCHEMA.md)
"""
)
