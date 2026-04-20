"""MultiPeakMountFujiSimulation: max over several single-peak landscapes."""

from __future__ import annotations

import streamlit as st
from gpmap.simulate import MultiPeakMountFujiSimulation
from utils.charts import peaks_scatter
from utils.controls import pick_wildtype_and_alphabet, rng_seed_input
from utils.ui import stats_row

st.markdown(
    "**Multi-peak Fuji** takes the max over `peak_n` single-peak Fujis, "
    "spacing the peaks by at least `min_peak_distance` in Hamming distance. "
    "v2 caps the retry loop so infeasible constraints raise instead of "
    "spinning forever."
)

with st.container(border=True):
    st.caption("Space")
    wildtype, mutations = pick_wildtype_and_alphabet("mp", default_length=6)
    c1, c2, c3 = st.columns(3)
    with c1:
        peak_n = st.slider("peak_n", 2, 5, value=2)
    with c2:
        min_peak_distance = st.slider(
            "min_peak_distance", 1, max(1, len(wildtype)), value=min(2, len(wildtype))
        )
    with c3:
        field_strength = st.slider("field_strength", 0.1, 3.0, value=1.0, step=0.1)
    c4, c5 = st.columns(2)
    with c4:
        roughness_width = st.slider("roughness_width", 0.0, 1.0, value=0.1, step=0.05)
    with c5:
        rng = rng_seed_input("mp")

try:
    sim = MultiPeakMountFujiSimulation(
        wildtype=wildtype,
        mutations=mutations,
        peak_n=peak_n,
        min_peak_distance=min_peak_distance,
        field_strength=field_strength,
        roughness_width=roughness_width,
        rng=rng,
    )
except RuntimeError as exc:
    st.error(f"Peak search failed: {exc}")
    st.stop()

stats_row(
    [
        ("n genotypes", sim.n),
        ("peaks", peak_n),
    ]
)

peak_rows = getattr(sim, "_peak_rows", [])
st.markdown("#### Peak genotypes")
st.code("\n".join(sim.peak_genotypes))

st.markdown("#### Phenotype vs Hamming distance (peaks in red)")
st.plotly_chart(peaks_scatter(sim, peak_rows), width="stretch")

with st.expander("Code", icon=":material/code:"):
    st.code(
        """import numpy as np
from gpmap.simulate import MultiPeakMountFujiSimulation

sim = MultiPeakMountFujiSimulation(
    wildtype="000000",
    mutations={i: ["0", "1"] for i in range(6)},
    peak_n=2,
    min_peak_distance=3,
    field_strength=1.0,
    roughness_width=0.1,
    rng=np.random.default_rng(0),
)
sim.peak_genotypes     # list of peak strings
""",
        language="python",
    )
