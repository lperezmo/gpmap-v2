"""NKSimulation: tunable ruggedness via K."""

from __future__ import annotations

import streamlit as st
from gpmap.simulate import NKSimulation
from utils.charts import phenotype_histogram, phenotype_vs_hamming
from utils.controls import pick_wildtype_and_alphabet, rng_seed_input
from utils.ui import render_landscape, stats_row

st.markdown(
    "**Kauffman's NK** model places a random fitness contribution at every "
    "site, each of which depends on the focal site plus its `K` neighbors. "
    "`K=0` is additive (smooth). `K=L-1` is House of Cards (fully random). "
    "Pull the slider to feel the phase transition."
)

with st.container(border=True):
    st.caption("Space")
    wildtype, mutations = pick_wildtype_and_alphabet("nk", default_length=6)
    max_k = max(len(wildtype) - 1, 1)
    c1, c2 = st.columns(2)
    with c1:
        K = st.slider("K (neighborhood size)", 0, max_k, value=min(2, max_k))
    with c2:
        rng = rng_seed_input("nk")

sim = NKSimulation(wildtype=wildtype, mutations=mutations, K=K, rng=rng)

stats_row(
    [
        ("n genotypes", sim.n),
        ("K", K),
        ("max phenotype", f"{sim.phenotypes.max():.3f}"),
    ]
)

st.markdown("#### Phenotype vs Hamming distance")
st.plotly_chart(phenotype_vs_hamming(sim), width="stretch")

st.markdown("#### Phenotype distribution")
st.plotly_chart(phenotype_histogram(sim), width="stretch")

render_landscape(sim)

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
