"""
components/jaccard_heatmap.py

Jaccard similarity heatmap for selected variants.

Shows pairwise mutation signature overlap between panel variants.
High overlap (coral) = variants share many mutations = harder for LolliPop
to discriminate between them in deconvolution.
Low overlap (teal) = distinct signatures = easier to discriminate.

Pure pango signature math — no cooc or deconv results needed.
Computable immediately when variants are selected.
"""
from __future__ import annotations

import plotly.graph_objects as go
import streamlit as st


def compute_jaccard_matrix(
    variants: list[str],
    pango_loader,
) -> list[list[float]]:
    """
    Compute pairwise Jaccard similarity matrix.

    Jaccard(A, B) = |sig(A) ∩ sig(B)| / |sig(A) ∪ sig(B)|

    Returns a 2D list (n x n) of overlap values 0.0–1.0.
    Diagonal is always 1.0 (variant compared to itself).
    """
    n = len(variants)
    matrix = [[0.0] * n for _ in range(n)]

    for i, v1 in enumerate(variants):
        for j, v2 in enumerate(variants):
            if i == j:
                matrix[i][j] = 1.0
                continue
            s1 = pango_loader.get_signature(v1)
            s2 = pango_loader.get_signature(v2)
            if not s1 or not s2:
                matrix[i][j] = 0.0
                continue
            intersection = len(s1 & s2)
            union = len(s1 | s2)
            matrix[i][j] = intersection / union if union > 0 else 0.0

    return matrix


def render_jaccard_heatmap(
    variants: list[str],
    pango_loader,
) -> None:
    """
    Render a Jaccard similarity heatmap for the selected variants.

    Args:
        variants: List of selected pango lineage names
        pango_loader: PangoLoader instance
    """
    if len(variants) < 2:
        st.caption("Select at least 2 variants to see similarity matrix.")
        return

    matrix = compute_jaccard_matrix(variants, pango_loader)

    fig = go.Figure(data=go.Heatmap(
        z=matrix,
        x=variants,
        y=variants,
        text=[[f"{matrix[i][j]:.2f}" for j in range(len(variants))]
              for i in range(len(variants))],
        texttemplate="%{text}",
        textfont={"size": 11},
        hovertemplate="%{y} ↔ %{x}<br>Jaccard: %{z:.2f}<extra></extra>",
        colorscale=[
            [0.0,  "#1D9E75"],
            [0.3,  "#E1F5EE"],
            [0.5,  "#F1EFE8"],
            [0.7,  "#F5C4B3"],
            [1.0,  "#712B13"],
        ],
        zmin=0,
        zmax=1,
        showscale=True,
        colorbar=dict(
            title=dict(text="overlap", side="right"),
            thickness=12,
            len=0.8,
            tickvals=[0, 0.25, 0.5, 0.75, 1.0],
            ticktext=["0.0<br>easy", "0.25", "0.5", "0.75", "1.0<br>hard"],
        ),
    ))

    fig.update_layout(
        title=dict(
            text="Signature similarity (Jaccard index)",
            font=dict(size=13),
            x=0,
        ),
        margin=dict(l=0, r=0, t=40, b=0),
        height=max(200, len(variants) * 44 + 80),
        xaxis=dict(
            tickfont=dict(size=11),
            side="bottom",
        ),
        yaxis=dict(
            tickfont=dict(size=11),
            autorange="reversed",
        ),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    st.plotly_chart(fig, use_container_width=True)
    st.caption(
        "Teal = low overlap, easy to discriminate. "
        "Coral = high overlap, harder for LolliPop to distinguish. "
        "Diagonal = variant compared to itself (always 1.0)."
    )
    