"""Encoding table and binary conversion."""

from __future__ import annotations

import streamlit as st
from gpmap import genotypes_to_binary_packed, get_encoding_table
from utils.controls import pick_wildtype_and_alphabet

st.markdown(
    "`gpmap-v2` encodes genotypes with a **unary-minus-one** scheme: a site "
    "with alphabet size *k* gets *k-1* binary bits. The wildtype row is all "
    "zeros; each non-WT letter owns a single bit. The full scheme is locked "
    "in [SCHEMA.md](https://github.com/lperezmo/gpmap-v2/blob/main/SCHEMA.md) "
    "section 3."
)

with st.container(border=True):
    st.caption("Build an encoding table")
    wildtype, mutations = pick_wildtype_and_alphabet(
        "encoding", default_length=3, default_kind="DNA"
    )

table = get_encoding_table(wildtype, mutations)
st.markdown("#### encoding_table")
st.dataframe(table, width="stretch", hide_index=True)

st.markdown("#### genotypes_to_binary_packed")
st.markdown(
    "Pass any list of genotype strings through the same table to get a "
    "`(n_genotypes, n_bits)` uint8 matrix. This is the fast path every "
    "downstream consumer should prefer."
)

sample_genotypes = [wildtype]
for i in range(len(wildtype)):
    letters = mutations[i]
    if len(letters) > 1:
        alt = letters[1]
        g = list(wildtype)
        g[i] = alt
        sample_genotypes.append("".join(g))

packed = genotypes_to_binary_packed(sample_genotypes, table)

left, right = st.columns([1, 1])
with left:
    st.caption("genotypes")
    st.code("\n".join(sample_genotypes))
with right:
    st.caption(f"binary_packed  shape={packed.shape}")
    st.code(str(packed))

with st.expander("Code", icon=":material/code:"):
    st.code(
        """from gpmap import get_encoding_table, genotypes_to_binary_packed

mutations = {0: ["A", "C", "G", "T"], 1: ["A", "T"], 2: ["A", "C"]}
table = get_encoding_table(wildtype="AAA", mutations=mutations)

packed = genotypes_to_binary_packed(["AAA", "CAA", "ATA", "AAC"], table)
# packed.shape -> (4, n_bits), packed.dtype -> uint8
""",
        language="python",
    )
