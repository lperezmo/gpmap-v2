"""Enumerate the Cartesian product safely."""

from __future__ import annotations

import streamlit as st
from gpmap import enumerate_genotypes_str
from gpmap.exceptions import SpaceTooLargeError
from utils.controls import ALPHABETS

st.markdown(
    "`enumerate_genotypes_str` walks the full Cartesian product of per-site "
    "alphabets. It is guarded against accidental 10^20 allocations via "
    "`SpaceTooLargeError`. Toggle `allow_huge=True` to override, at your own risk."
)

with st.container(border=True):
    st.caption("Small space")
    c1, c2, c3 = st.columns(3)
    with c1:
        kind = st.selectbox("Alphabet", options=list(ALPHABETS), index=list(ALPHABETS).index("DNA"))
    source = ALPHABETS[kind]
    with c2:
        length = st.slider("Length", 2, 6, value=3)
    with c3:
        if len(source) > 2:
            alpha_size = st.slider("Letters per site", 2, len(source), value=min(3, len(source)))
        else:
            st.markdown("**Letters per site**")
            st.caption(f"{len(source)} (fixed for {kind})")
            alpha_size = len(source)

letters = list(source[:alpha_size])
wildtype = letters[0] * length
mutations = {i: letters[:] for i in range(length)}

total = alpha_size**length
st.info(f"Theoretical size = {alpha_size}^{length} = **{total:,}** genotypes.")

genotypes = enumerate_genotypes_str(wildtype, mutations)
st.code("\n".join(genotypes[:32]) + ("\n..." if len(genotypes) > 32 else ""))

st.divider()
st.subheader("The size guard")
st.markdown(
    "The default `max_genotypes` is **2**^28 (~268M rows). Ask for more without "
    "`allow_huge=True` and you get a hard error before numpy ever tries to "
    "allocate. Hit the button to watch it trip on purpose."
)

if st.button("Request a 4^40 space (impossible)", type="primary"):
    try:
        enumerate_genotypes_str("A" * 40, {i: ["A", "C", "G", "T"] for i in range(40)})
    except SpaceTooLargeError as exc:
        st.error(f"SpaceTooLargeError raised as expected:\n\n`{exc}`")

with st.expander("Code", icon=":material/code:"):
    st.code(
        """from gpmap import enumerate_genotypes_str
from gpmap.exceptions import SpaceTooLargeError

mutations = {i: ["A", "C", "G", "T"] for i in range(3)}
genos = enumerate_genotypes_str("AAA", mutations)
# len(genos) == 64

try:
    enumerate_genotypes_str("A" * 40, {i: ["A", "C", "G", "T"] for i in range(40)})
except SpaceTooLargeError as e:
    print("caught:", e)
""",
        language="python",
    )
