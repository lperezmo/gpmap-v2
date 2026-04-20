"""Shared widgets used across demo pages."""

from __future__ import annotations

import numpy as np
import streamlit as st

ALPHABETS: dict[str, tuple[str, ...]] = {
    "BINARY": ("0", "1"),
    "DNA": ("A", "C", "G", "T"),
    "RNA": ("A", "C", "G", "U"),
    "AA": (
        "A",
        "R",
        "N",
        "D",
        "C",
        "E",
        "Q",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
    ),
}


def pick_wildtype_and_alphabet(
    key_prefix: str,
    *,
    default_length: int = 4,
    max_length: int = 8,
    default_kind: str = "BINARY",
) -> tuple[str, dict[int, list[str]]]:
    """Render length + alphabet controls and return (wildtype, mutations)."""
    c1, c2, c3 = st.columns(3)
    with c1:
        kind = st.selectbox(
            "Alphabet",
            options=list(ALPHABETS),
            index=list(ALPHABETS).index(default_kind),
            key=f"{key_prefix}_kind",
        )
    source = ALPHABETS[kind]
    with c2:
        length = st.slider(
            "Wildtype length",
            min_value=2,
            max_value=max_length,
            value=default_length,
            key=f"{key_prefix}_length",
        )
    with c3:
        if len(source) > 2:
            alpha_size = st.slider(
                "Letters per site",
                min_value=2,
                max_value=len(source),
                value=2,
                key=f"{key_prefix}_alpha",
            )
        else:
            st.markdown("**Letters per site**")
            st.caption(f"{len(source)} (fixed for {kind})")
            alpha_size = len(source)

    letters = list(source[:alpha_size])
    wildtype = letters[0] * length
    mutations = {i: letters[:] for i in range(length)}
    return wildtype, mutations


def rng_seed_input(key_prefix: str, *, default: int = 0) -> np.random.Generator:
    seed = st.number_input(
        "Random seed",
        min_value=0,
        max_value=1_000_000,
        value=default,
        step=1,
        key=f"{key_prefix}_seed",
    )
    return np.random.default_rng(int(seed))
