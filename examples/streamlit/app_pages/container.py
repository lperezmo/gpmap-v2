"""Container page: build a GenotypePhenotypeMap from raw inputs."""

from __future__ import annotations

import numpy as np
import streamlit as st
from gpmap import GenotypePhenotypeMap
from gpmap.simulate import RandomPhenotypesSimulation
from utils.controls import pick_wildtype_and_alphabet, rng_seed_input

st.markdown(
    "The central object is `GenotypePhenotypeMap`. Hand it a wildtype, a list "
    "of genotypes, and matching phenotypes, and it becomes a typed view over "
    "the whole space (encoding table, packed binary, Hamming weight, and a "
    "pandas `data` frame, all lazy and cached)."
)

with st.container(border=True):
    st.caption("Choose the space")
    wildtype, mutations = pick_wildtype_and_alphabet("container", default_length=4)
    rng = rng_seed_input("container")

# Build a small space of real phenotypes (RandomPhenotypesSimulation enumerates).
sim = RandomPhenotypesSimulation(wildtype=wildtype, mutations=mutations, rng=rng)

# Now rebuild as a plain GenotypePhenotypeMap to show the container API directly.
gpm = GenotypePhenotypeMap(
    wildtype=sim.wildtype,
    genotypes=sim.genotypes.tolist(),
    phenotypes=sim.phenotypes,
    stdeviations=np.full(sim.n, 0.05),
    mutations=sim.mutations,
)

stats_row(
    [
        ("wildtype", gpm.wildtype),
        ("length (L)", gpm.length),
        ("genotypes (n)", gpm.n),
        ("binary bits", gpm.binary_packed.shape[1]),
    ]
)

st.markdown("#### .data  (pandas view)")
st.dataframe(gpm.data, width="stretch", hide_index=True)

left, right = st.columns(2)
with left:
    st.markdown("#### .phenotypes")
    st.code(np.array2string(gpm.phenotypes, precision=3, max_line_width=60))
with right:
    st.markdown("#### .n_mutations")
    st.code(np.array2string(gpm.n_mutations, max_line_width=60))

with st.expander("Code", icon=":material/code:"):
    st.code(
        """from gpmap import GenotypePhenotypeMap

gpm = GenotypePhenotypeMap(
    wildtype="00",
    genotypes=["00", "01", "10", "11"],
    phenotypes=[0.1, 0.4, 0.5, 0.9],
    stdeviations=[0.05] * 4,
)

gpm.data            # pandas DataFrame of the whole thing
gpm.phenotypes      # np.ndarray[float64]
gpm.n_mutations     # per-row Hamming weight
gpm.binary_packed   # np.ndarray[uint8], shape (n, n_bits)
""",
        language="python",
    )
