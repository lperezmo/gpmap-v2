"""RandomPhenotypesSimulation: phenotypes drawn uniformly in [low, high]."""

from __future__ import annotations

import streamlit as st
from gpmap.simulate import RandomPhenotypesSimulation
from utils.charts import phenotype_histogram, phenotype_vs_hamming
from utils.controls import pick_wildtype_and_alphabet, rng_seed_input

st.markdown(
    "**Random** phenotypes, no landscape structure at all. Useful as a "
    "baseline or as a quick way to get a valid `GenotypePhenotypeMap` of a "
    "given shape."
)

with st.container(border=True):
    st.caption("Space")
    wildtype, mutations = pick_wildtype_and_alphabet("rand", default_length=5)
    c1, c2, c3 = st.columns(3)
    with c1:
        low = st.number_input("low", value=0.0, step=0.1)
    with c2:
        high = st.number_input("high", value=1.0, step=0.1)
    with c3:
        rng = rng_seed_input("rand")

if high <= low:
    st.warning("`high` must be greater than `low`.")
    st.stop()

sim = RandomPhenotypesSimulation(
    wildtype=wildtype, mutations=mutations, low=low, high=high, rng=rng
)

c1, c2 = st.columns(2)
c1.metric("n genotypes", sim.n)
c2.metric("range", f"[{low}, {high}]")

st.subheader("Phenotype vs Hamming distance")
st.plotly_chart(phenotype_vs_hamming(sim), use_container_width=True)

st.subheader("Phenotype distribution")
st.plotly_chart(phenotype_histogram(sim), use_container_width=True)

with st.expander("Code", icon=":material/code:"):
    st.code(
        """import numpy as np
from gpmap.simulate import RandomPhenotypesSimulation

sim = RandomPhenotypesSimulation(
    wildtype="00000",
    mutations={i: ["0", "1"] for i in range(5)},
    low=0.0,
    high=1.0,
    rng=np.random.default_rng(0),
)
""",
        language="python",
    )
