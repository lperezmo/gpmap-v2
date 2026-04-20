"""mask(): uniformly subsample a GenotypePhenotypeMap."""

from __future__ import annotations

import streamlit as st
from gpmap.simulate import MountFujiSimulation, mask
from utils.charts import phenotype_vs_hamming
from utils.controls import pick_wildtype_and_alphabet, rng_seed_input
from utils.ui import stats_row

st.markdown(
    "`mask(gpm, fraction)` returns a `MaskedGPM(fraction, gpm)` NamedTuple "
    "containing a new `GenotypePhenotypeMap` keeping `fraction` of the rows "
    "uniformly at random. Useful for simulating incomplete measurements or "
    "building train/test splits."
)

with st.container(border=True):
    st.caption("Source landscape")
    wildtype, mutations = pick_wildtype_and_alphabet("mask", default_length=5)
    c1, c2 = st.columns(2)
    with c1:
        fraction = st.slider("fraction", 0.05, 1.0, value=0.3, step=0.05)
    with c2:
        rng = rng_seed_input("mask")

full = MountFujiSimulation(
    wildtype=wildtype,
    mutations=mutations,
    field_strength=1.0,
    roughness_width=0.2,
    rng=rng,
)
masked = mask(full, fraction=fraction, rng=rng)

stats_row(
    [
        ("full n", full.n),
        ("masked n", masked.gpm.n),
        ("actual fraction", f"{masked.fraction:.3f}"),
    ]
)

left, right = st.columns(2)
with left:
    st.caption("Full landscape")
    st.plotly_chart(phenotype_vs_hamming(full), width="stretch")
with right:
    st.caption("After mask")
    st.plotly_chart(phenotype_vs_hamming(masked.gpm), width="stretch")

with st.expander("Code", icon=":material/code:"):
    st.code(
        """import numpy as np
from gpmap.simulate import MountFujiSimulation, mask

full = MountFujiSimulation(
    wildtype="00000",
    mutations={i: ["0", "1"] for i in range(5)},
    rng=np.random.default_rng(0),
)
masked = mask(full, fraction=0.3, rng=np.random.default_rng(1))
masked.fraction   # actual fraction kept, e.g. 0.3125
masked.gpm        # new GenotypePhenotypeMap
""",
        language="python",
    )
