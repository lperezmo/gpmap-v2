"""HouseOfCardsSimulation: maximally rugged landscape."""

from __future__ import annotations

import streamlit as st
from gpmap.simulate import HouseOfCardsSimulation
from utils.charts import phenotype_histogram, phenotype_vs_hamming
from utils.controls import pick_wildtype_and_alphabet, rng_seed_input

st.markdown(
    "**House of Cards** is NK with `K = n_bits - 1`: every site sees the full "
    "genotype, so there is no useful local structure. Expect phenotypes to "
    "look like noise, with Hamming distance providing no signal."
)

with st.container(border=True):
    st.caption("Space")
    wildtype, mutations = pick_wildtype_and_alphabet("hoc", default_length=5)
    rng = rng_seed_input("hoc")

sim = HouseOfCardsSimulation(wildtype=wildtype, mutations=mutations, rng=rng)

c1, c2 = st.columns(2)
c1.metric("n genotypes", sim.n)
c2.metric("max phenotype", f"{sim.phenotypes.max():.3f}")

st.subheader("Phenotype vs Hamming distance")
st.plotly_chart(phenotype_vs_hamming(sim), width='stretch')

st.subheader("Phenotype distribution")
st.plotly_chart(phenotype_histogram(sim), width='stretch')

with st.expander("Code", icon=":material/code:"):
    st.code(
        """import numpy as np
from gpmap.simulate import HouseOfCardsSimulation

sim = HouseOfCardsSimulation(
    wildtype="00000",
    mutations={i: ["0", "1"] for i in range(5)},
    rng=np.random.default_rng(0),
)
""",
        language="python",
    )
