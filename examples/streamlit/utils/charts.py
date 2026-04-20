"""Plotly chart helpers shared across demo pages."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import plotly.express as px
import plotly.graph_objects as go


def phenotype_vs_hamming(gpm) -> go.Figure:
    fig = px.strip(
        x=gpm.n_mutations,
        y=gpm.phenotypes,
        labels={"x": "Hamming distance to wildtype", "y": "Phenotype"},
        stripmode="overlay",
    )
    fig.update_traces(jitter=0.35, marker=dict(size=7, opacity=0.75))
    fig.update_layout(height=360, margin=dict(l=40, r=20, t=20, b=40))
    return fig


def phenotype_histogram(gpm, *, nbins: int = 30) -> go.Figure:
    fig = px.histogram(
        x=gpm.phenotypes,
        nbins=nbins,
        labels={"x": "Phenotype"},
    )
    fig.update_layout(height=320, margin=dict(l=40, r=20, t=20, b=40), bargap=0.05)
    return fig


def peaks_scatter(gpm, peak_indices: Iterable[int]) -> go.Figure:
    peaks = set(int(i) for i in peak_indices)
    colors = np.where(
        np.fromiter((i in peaks for i in range(gpm.n)), dtype=bool),
        "crimson",
        "steelblue",
    )
    sizes = np.where(
        np.fromiter((i in peaks for i in range(gpm.n)), dtype=bool),
        14,
        7,
    )
    fig = go.Figure(
        go.Scatter(
            x=gpm.n_mutations,
            y=gpm.phenotypes,
            mode="markers",
            marker=dict(color=colors, size=sizes, opacity=0.8),
            hovertext=gpm.genotypes,
        )
    )
    fig.update_layout(
        xaxis_title="Hamming distance to wildtype",
        yaxis_title="Phenotype",
        height=380,
        margin=dict(l=40, r=20, t=20, b=40),
    )
    return fig
