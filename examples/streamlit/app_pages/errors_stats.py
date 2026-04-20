"""Error-bar transforms, stdev/sterror views, and unbiased stats."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
import streamlit as st
from gpmap import StandardDeviationMap, StandardErrorMap, lower_transform, upper_transform
from gpmap.simulate import MountFujiSimulation
from gpmap.stats import c4_correction, unbiased_std, unbiased_var
from utils.controls import rng_seed_input
from utils.ui import stats_row

st.markdown(
    "`gpmap-v2` ships matplotlib/plotly-friendly error-bar transforms, two "
    "view objects (`StandardDeviationMap`, `StandardErrorMap`), and a small "
    "stats module with c4-corrected unbiased estimators."
)

with st.container(border=True):
    st.caption("Simulated measurements")
    c1, c2, c3 = st.columns(3)
    with c1:
        length = st.slider("length", 3, 5, value=4)
    with c2:
        sigma = st.slider("noise sigma", 0.0, 1.0, value=0.2, step=0.05)
    with c3:
        rng = rng_seed_input("errs")

mutations = {i: ["0", "1"] for i in range(length)}
sim = MountFujiSimulation(wildtype="0" * length, mutations=mutations, field_strength=1.0, rng=rng)
sim.set_stdeviations(sigma)

st.markdown("#### Error bars with upper_transform / lower_transform")
st.markdown(
    "Matplotlib and plotly want a *positive offset* from the bar height to "
    "the upper/lower bound. `upper_transform` and `lower_transform` do that "
    "(with optional log support). In v1 they were accidentally identical; v2 "
    "makes them genuinely distinct."
)

bar_y = sim.phenotypes
upper = bar_y + sigma
lower = bar_y - sigma
upper_off = upper_transform(bar_y, upper)
lower_off = lower_transform(bar_y, lower)

fig = go.Figure(
    go.Bar(
        x=list(sim.genotypes),
        y=bar_y,
        error_y=dict(type="data", symmetric=False, array=upper_off, arrayminus=lower_off),
        marker_color="steelblue",
    )
)
fig.update_layout(
    xaxis_title="genotype",
    yaxis_title="phenotype",
    height=360,
    margin=dict(l=40, r=20, t=20, b=40),
)
st.plotly_chart(fig, width="stretch")

st.markdown("#### StandardDeviationMap vs StandardErrorMap")
stdev_map = StandardDeviationMap(sim)
err_map = StandardErrorMap(sim)
c1, c2 = st.columns(2)
with c1:
    st.caption("StandardDeviationMap.upper  (raw stdev)")
    st.code(np.array2string(stdev_map.upper[:8], precision=3))
with c2:
    st.caption("StandardErrorMap.upper  (std / sqrt(n_replicates))")
    st.code(np.array2string(err_map.upper[:8], precision=3))

st.markdown("#### Unbiased statistics")
samples = rng.normal(loc=0.0, scale=sigma, size=(sim.n, 8))
stats_row(
    [
        ("unbiased_std", f"{float(unbiased_std(samples.ravel())):.4f}"),
        ("unbiased_var", f"{float(unbiased_var(samples.ravel())):.4f}"),
        ("c4(n=8)", f"{c4_correction(8):.4f}"),
    ]
)

with st.expander("Code", icon=":material/code:"):
    st.code(
        """from gpmap import StandardDeviationMap, upper_transform, lower_transform
from gpmap.stats import unbiased_std, unbiased_var, c4_correction

upper_off = upper_transform(bar_y, upper_bounds)
lower_off = lower_transform(bar_y, lower_bounds)

stdev_map = StandardDeviationMap(gpm)
stdev_map.upper   # np.ndarray of positive offsets

unbiased_std(samples)         # c4-corrected
unbiased_var(samples)
c4_correction(n_samples=8)
""",
        language="python",
    )
