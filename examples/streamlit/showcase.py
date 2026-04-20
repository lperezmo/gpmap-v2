"""
gpmap-v2 Showcase
=================
Interactive demo of the gpmap-v2 container and simulators.
Structure follows streamlit-aggrid-v2's multi-page pattern.
"""

from __future__ import annotations

import sys
from pathlib import Path

import streamlit as st

# Make `utils` importable when pages run (Streamlit executes pages with their
# own module context).
_HERE = Path(__file__).parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import gpmap  # noqa: E402

st.set_page_config(
    page_title="gpmap-v2 | genotype-phenotype map toolkit",
    page_icon=":material/hub:",
    layout="wide",
    initial_sidebar_state="collapsed",
)

st.markdown(
    """<style>.block-container { padding-top: 1rem; }</style>""",
    unsafe_allow_html=True,
)

page = st.navigation(
    {
        "": [
            st.Page("app_pages/home.py", title="Home", icon=":material/home:"),
        ],
        "Core": [
            st.Page(
                "app_pages/container.py",
                title="Container",
                icon=":material/inventory_2:",
            ),
            st.Page(
                "app_pages/encoding.py",
                title="Encoding table",
                icon=":material/table_view:",
            ),
            st.Page(
                "app_pages/enumerate.py",
                title="Enumerate",
                icon=":material/list:",
            ),
        ],
        "Simulators": [
            st.Page(
                "app_pages/sim_fuji.py",
                title="Mount Fuji",
                icon=":material/landscape:",
            ),
            st.Page("app_pages/sim_nk.py", title="NK", icon=":material/terrain:"),
            st.Page(
                "app_pages/sim_hoc.py",
                title="House of Cards",
                icon=":material/shuffle:",
            ),
            st.Page(
                "app_pages/sim_multipeak.py",
                title="Multi-peak Fuji",
                icon=":material/filter_hdr:",
            ),
            st.Page(
                "app_pages/sim_random.py",
                title="Random",
                icon=":material/casino:",
            ),
        ],
        "Utilities": [
            st.Page(
                "app_pages/masking.py",
                title="Masking",
                icon=":material/filter_alt:",
            ),
            st.Page(
                "app_pages/errors_stats.py",
                title="Errors and stats",
                icon=":material/functions:",
            ),
            st.Page("app_pages/io.py", title="I/O round-trips", icon=":material/save:"),
        ],
    },
    position="top",
)

_IS_DARK = st.context.theme.type == "dark"
_HEADER_GRADIENT = (
    "linear-gradient(135deg, #7dd3fc, #c4b5fd)"
    if _IS_DARK
    else "linear-gradient(135deg, #0369a1, #6d28d9)"
)
st.html(f"""
<div style="text-align:center; padding:1.25rem 0 0.25rem;">
    <h2 style="
        margin:0; font-size:1.5rem; font-weight:700;
        background:{_HEADER_GRADIENT};
        -webkit-background-clip:text; -webkit-text-fill-color:transparent;
        background-clip:text;
    ">gpmap-v2</h2>
    <p style="margin:0.3rem 0 0; font-size:0.85rem; opacity:0.7;">
        Typed, Rust-accelerated genotype-phenotype map toolkit
    </p>
</div>
""")

with st.sidebar:
    st.divider()
    st.markdown(
        "**gpmap-v2** is a clean-break rewrite of `harmslab/gpmap`. Same "
        "conceptual model (`GenotypePhenotypeMap` + simulators), locked "
        "schema, vectorized hot paths, PyO3 Rust core."
    )
    st.caption(f"package: gpmap-v2 {gpmap.__version__}")

st.markdown(f"#### {page.title}")
page.run()

st.divider()
st.caption(
    "Built with [gpmap-v2](https://github.com/lperezmo/gpmap-v2) - "
    "schema contract in [SCHEMA.md](https://github.com/lperezmo/gpmap-v2/blob/main/SCHEMA.md)"
)
