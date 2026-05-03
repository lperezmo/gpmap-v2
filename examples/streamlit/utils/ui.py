"""Small UI primitives. Avoid Streamlit's heavy default headers and metric boxes."""

from __future__ import annotations

import streamlit as st


def stat(label: str, value: object) -> None:
    """Compact stat: tiny uppercase label, bold value. No box."""
    st.html(
        f"""
<div style="line-height:1.2; padding:0.1rem 0;">
  <div style="font-size:0.7rem; opacity:0.55; text-transform:uppercase;
              letter-spacing:0.06em; font-weight:500;">{label}</div>
  <div style="font-size:1.05rem; font-weight:600; margin-top:0.15rem;">{value}</div>
</div>
"""
    )


def stats_row(items: list[tuple[str, object]]) -> None:
    """Render a horizontal row of `stat` tiles."""
    cols = st.columns(len(items))
    for col, (label, value) in zip(cols, items, strict=False):
        with col:
            stat(label, value)


def section(title: str, caption: str | None = None) -> None:
    """Section heading without Streamlit's oversized default."""
    st.markdown(f"#### {title}")
    if caption:
        st.caption(caption)


def render_landscape(sim) -> None:
    """Render phenotype landscape charts for a simulation.

    Shows a 3D surface as the primary visual when the space is small enough
    (N <= 64, biallelic). A 2D heatmap is shown below the surface, or on its
    own when the landscape is too large for a surface (up to N=512).
    """
    from utils.charts import phenotype_landscape_heatmap, phenotype_landscape_surface

    surface = phenotype_landscape_surface(sim)
    heatmap = phenotype_landscape_heatmap(sim)

    if surface is None and heatmap is None:
        return

    st.markdown("#### Phenotype landscape")
    if surface is not None:
        st.plotly_chart(surface, width="stretch")
    if heatmap is not None:
        if surface is not None:
            st.markdown("**2D view**")
        else:
            st.caption(
                "Genotype space split into two halves; colour encodes phenotype."
            )
        st.plotly_chart(heatmap, width="stretch")
