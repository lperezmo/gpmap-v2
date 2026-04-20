"""MountFujiSimulation: single-peak additive landscape."""

from __future__ import annotations

import streamlit as st
from gpmap.simulate import MountFujiSimulation
from utils.charts import phenotype_histogram, phenotype_vs_hamming
from utils.controls import pick_wildtype_and_alphabet, rng_seed_input
from utils.ui import stats_row

st.markdown(
    "**Mount Fuji** is a smooth, single-peak landscape: phenotype is a linear "
    "function of Hamming distance to wildtype, optionally perturbed by noise. "
    "The wildtype sits at the top; every additional mutation costs "
    "`field_strength`."
)

with st.container(border=True):
    st.caption("Space")
    wildtype, mutations = pick_wildtype_and_alphabet("fuji", default_length=5)
    c1, c2, c3 = st.columns(3)
    with c1:
        field_strength = st.slider("field_strength", 0.0, 5.0, value=1.0, step=0.1)
    with c2:
        roughness_width = st.slider("roughness_width", 0.0, 2.0, value=0.2, step=0.05)
    with c3:
        roughness_dist = st.radio("roughness_dist", options=["normal", "uniform"], horizontal=True)
    rng = rng_seed_input("fuji")

sim = MountFujiSimulation(
    wildtype=wildtype,
    mutations=mutations,
    field_strength=field_strength,
    roughness_width=roughness_width,
    roughness_dist=roughness_dist,
    rng=rng,
)

stats_row(
    [
        ("n genotypes", sim.n),
        ("peak phenotype", f"{sim.phenotypes.max():.3f}"),
    ]
)

st.markdown("#### Phenotype vs Hamming distance")
st.plotly_chart(phenotype_vs_hamming(sim), width="stretch")

st.markdown("#### Phenotype distribution")
st.plotly_chart(phenotype_histogram(sim), width="stretch")

with st.expander("Code", icon=":material/code:"):
    st.code(
        """import numpy as np
from gpmap.simulate import MountFujiSimulation

sim = MountFujiSimulation(
    wildtype="00000",
    mutations={i: ["0", "1"] for i in range(5)},
    field_strength=1.0,
    roughness_width=0.2,
    roughness_dist="normal",
    rng=np.random.default_rng(0),
)
sim.phenotypes.shape   # (32,)
sim.data               # pandas view
""",
        language="python",
    )
