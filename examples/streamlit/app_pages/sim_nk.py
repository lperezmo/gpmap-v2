"""NKSimulation: tunable ruggedness via K."""

from __future__ import annotations

import streamlit as st
from gpmap.simulate import NKSimulation
from utils.charts import phenotype_histogram, phenotype_vs_hamming
from utils.controls import pick_wildtype_and_alphabet, rng_seed_input

st.markdown(
    "**Kauffman's NK** model places a random fitness contribution at every "
    "site, each of which depends on the focal site plus its `K` neighbors. "
    "`K=0` is additive (smooth). `K=L-1` is House of Cards (fully random). "
    "Pull the slider to feel the phase transition."
)

with st.container(border=True):
    st.caption("Space")
    wildtype, mutations = pick_wildtype_and_alphabet("nk", default_length=6)
    max_K = max(len(wildtype) - 1, 1)
    c1, c2 = st.columns(2)
    with c1:
        K = st.slider("K (neighborhood size)", 0, max_K, value=min(2, max_K))
    with c2:
        rng = rng_seed_input("nk")

sim = NKSimulation(wildtype=wildtype, mutations=mutations, K=K, rng=rng)

c1, c2, c3 = st.columns(3)
c1.metric("n genotypes", sim.n)
c2.metric("K", K)
c3.metric("max phenotype", f"{sim.phenotypes.max():.3f}")

st.subheader("Phenotype vs Hamming distance")
st.plotly_chart(phenotype_vs_hamming(sim), use_container_width=True)

st.subheader("Phenotype distribution")
st.plotly_chart(phenotype_histogram(sim), use_container_width=True)

with st.expander("Code", icon=":material/code:"):
    st.code(
        '''import numpy as np
from gpmap.simulate import NKSimulation

sim = NKSimulation(
    wildtype="000000",
    mutations={i: ["0", "1"] for i in range(6)},
    K=2,
    rng=np.random.default_rng(0),
)
sim.phenotypes.shape   # (64,)
''',
        language="python",
    )
