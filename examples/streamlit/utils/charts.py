"""Plotly chart helpers shared across demo pages."""

from __future__ import annotations

import itertools
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


def phenotype_landscape_heatmap(gpm) -> go.Figure | None:
    """2D phenotype heatmap splitting sites into two halves.

    Returns None when the landscape is too large to display (N > 512).
    """
    if gpm.n > 512:
        return None
    L = gpm.length
    mutations = gpm.mutations
    l1 = (L + 1) // 2
    l2 = L - l1
    if l2 == 0:
        return None
    alleles1 = [mutations[s] for s in range(l1)]
    alleles2 = [mutations[s] for s in range(l1, L)]
    rows = ["".join(c) for c in itertools.product(*alleles1)]
    cols = ["".join(c) for c in itertools.product(*alleles2)]
    geno_pheno = {str(g): float(p) for g, p in zip(gpm.genotypes, gpm.phenotypes, strict=False)}
    z = np.full((len(rows), len(cols)), float("nan"))
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            val = geno_pheno.get(r + c)
            if val is not None:
                z[i, j] = val
    row_title = f"sites 1-{l1}"
    col_title = f"sites {l1 + 1}-{L}"
    height = max(280, min(540, len(rows) * 28 + 120))
    fig = go.Figure(
        go.Heatmap(
            z=z,
            x=cols,
            y=rows,
            colorscale="Viridis",
            colorbar=dict(title="phenotype", thickness=14),
            hovertemplate="<b>%{y}%{x}</b><br>phenotype=%{z:.4f}<extra></extra>",
        )
    )
    fig.update_layout(
        xaxis=dict(title=col_title, tickangle=-45, side="bottom"),
        yaxis=dict(title=row_title, autorange="reversed"),
        height=height,
        margin=dict(l=70, r=20, t=20, b=60),
    )
    return fig


def phenotype_landscape_surface(gpm) -> go.Figure | None:
    """3D surface of the phenotype landscape for small landscapes (N <= 64).

    Returns None when the landscape is too large or has non-biallelic sites.
    """
    if gpm.n > 64:
        return None
    mutations = gpm.mutations
    if not all(len(v) == 2 for v in mutations.values()):
        return None
    L = gpm.length
    l1 = (L + 1) // 2
    l2 = L - l1
    if l2 == 0:
        return None
    alleles1 = [mutations[s] for s in range(l1)]
    alleles2 = [mutations[s] for s in range(l1, L)]
    rows = ["".join(c) for c in itertools.product(*alleles1)]
    cols = ["".join(c) for c in itertools.product(*alleles2)]
    geno_pheno = {str(g): float(p) for g, p in zip(gpm.genotypes, gpm.phenotypes, strict=False)}
    z = np.full((len(rows), len(cols)), float("nan"))
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            val = geno_pheno.get(r + c)
            if val is not None:
                z[i, j] = val
    fig = go.Figure(
        go.Surface(
            z=z,
            x=list(range(len(cols))),
            y=list(range(len(rows))),
            colorscale="Viridis",
            showscale=True,
            colorbar=dict(title="phenotype", thickness=14),
        )
    )
    fig.update_layout(
        scene=dict(
            xaxis=dict(
                title=f"sites {l1 + 1}-{L}",
                tickvals=list(range(len(cols))),
                ticktext=cols,
            ),
            yaxis=dict(
                title=f"sites 1-{l1}",
                tickvals=list(range(len(rows))),
                ticktext=rows,
            ),
            zaxis=dict(title="phenotype"),
        ),
        height=480,
        margin=dict(l=0, r=0, t=10, b=0),
    )
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
