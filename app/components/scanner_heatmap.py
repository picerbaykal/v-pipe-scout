"""Scanner heatmap component.

For each missing variant found by the scanner, shows a heatmap of
mutation frequencies over time using the queriesOverTime endpoint.

x-axis: date ranges (weeks)
y-axis: discriminating mutations observed in unexplained patterns
cells: frequency = count / coverage per mutation per week

Helps the user understand the strength and timing of the signal —
when did the missing variant start rising, and which mutations drive it.
"""

import asyncio
from datetime import datetime, timedelta
from typing import List, Dict

import pandas as pd
import plotly.graph_objects as go
import streamlit as st


def _build_queries(mutations: List[str]) -> List[Dict[str, str]]:
    """
    Convert a list of '{pos}{alt}' mutation strings into queriesOverTime
    query objects with count and coverage queries.

    e.g. '23018T' → {
        countQuery: 'main:23018T',
        coverageQuery: '!main:23018N',
        displayLabel: '23018T'
    }
    """
    queries = []
    for mut in mutations:
        pos = mut[:-1]
        queries.append({
            "displayLabel": mut,
            "countQuery": f"main:{pos}{mut[-1]}",
            "coverageQuery": f"!main:{pos}N",
        })
    return queries


def _fetch_frequencies(
    client,
    location: str,
    date_range: tuple,
    mutations: List[str],
) -> pd.DataFrame:
    """
    Fetch mutation frequencies over time for the given mutations.
    Returns a DataFrame with columns: mutation, dateFrom, dateTo, frequency, count, coverage.
    """
    queries = _build_queries(mutations)

    async def _run():
        return await client.get_queries_over_time(
            locationName=location,
            date_range=date_range,
            queries=queries,
            date_granularity="week",
        )

    result = asyncio.run(_run())

    if not result:
        return pd.DataFrame()

    query_labels = result.get("queries", [])
    date_ranges = result.get("dateRanges", [])
    data = result.get("data", [])

    rows = []
    for q_idx, label in enumerate(query_labels):
        for d_idx, dr in enumerate(date_ranges):
            cell = data[q_idx][d_idx] if q_idx < len(data) and d_idx < len(data[q_idx]) else {}
            count = cell.get("count", 0)
            coverage = cell.get("coverage", 0)
            freq = count / coverage if coverage > 0 else 0.0
            rows.append({
                "mutation": label,
                "dateFrom": dr["dateFrom"],
                "dateTo": dr["dateTo"],
                "frequency": freq,
                "count": count,
                "coverage": coverage,
            })

    return pd.DataFrame(rows)


def render_scanner_heatmap(
    variant: str,
    mutations: List[str],
    client,
    location: str,
    date_range: tuple,
    max_mutations: int = 20,
) -> None:
    """
    Render a mutation frequency heatmap for a missing variant.

    Args:
        variant: Variant name e.g. "XFG"
        mutations: List of "{pos}{alt}" mutation strings observed in unexplained patterns
        client: WiseLoculusLapis instance
        location: Location name e.g. "Lugano (TI)"
        date_range: (start_datetime, end_datetime) tuple
        max_mutations: Cap the number of mutations shown (top by position)
    """
    if not mutations:
        st.caption(f"No discriminating mutations found for {variant}.")
        return

    # cap to avoid huge heatmaps
    shown = sorted(mutations)[:max_mutations]
    if len(mutations) > max_mutations:
        st.caption(
            f"Showing {max_mutations} of {len(mutations)} observed mutations for {variant}."
        )

    with st.spinner(f"Fetching {variant} signal over time…"):
        df = _fetch_frequencies(client, location, date_range, shown)

    if df.empty:
        st.caption("No frequency data returned.")
        return

    # pivot: rows = mutations, cols = date ranges
    pivot = df.pivot(index="mutation", columns="dateFrom", values="frequency")
    hover = df.pivot(index="mutation", columns="dateFrom", values="count")
    cov = df.pivot(index="mutation", columns="dateFrom", values="coverage")

    # x-axis labels: "May 01" style
    col_labels = [
        datetime.strptime(c, "%Y-%m-%d").strftime("%b %d")
        for c in pivot.columns
    ]

    # custom hover text
    hover_text = []
    for mut in pivot.index:
        row_hover = []
        for col in pivot.columns:
            f = pivot.loc[mut, col]
            c = hover.loc[mut, col] if mut in hover.index else 0
            v = cov.loc[mut, col] if mut in cov.index else 0
            row_hover.append(
                f"{mut}<br>freq: {f:.1%}<br>{int(c):,} / {int(v):,} reads"
            )
        hover_text.append(row_hover)

    fig = go.Figure(go.Heatmap(
        z=pivot.values.tolist(),
        x=col_labels,
        y=list(pivot.index),
        text=hover_text,
        hovertemplate="%{text}<extra></extra>",
        colorscale="YlOrRd",
        zmin=0,
        zmax=max(pivot.values.max(), 0.01),
        colorbar=dict(title="frequency", tickformat=".0%"),
    ))

    fig.update_layout(
        title=f"{variant} — mutation signal over time",
        height=max(200, 30 * len(shown) + 80),
        margin=dict(t=40, b=40, l=80, r=20),
        template="plotly_white",
        xaxis=dict(side="bottom"),
        yaxis=dict(autorange="reversed"),
    )

    st.plotly_chart(fig, use_container_width=True)
    st.caption(
        f"Frequency = reads carrying the mutation / reads covering that position. "
        f"Only mutations observed in unexplained co-occurrence patterns are shown."
    )