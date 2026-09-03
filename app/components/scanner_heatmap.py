"""Scanner heatmap component.
For each missing variant or emerging sublineage found by the scanner,
shows a heatmap of mutation frequencies over time using queriesOverTime.

x-axis: date ranges (weeks)
y-axis: discriminating mutations observed in unexplained patterns
cells:  frequency = count / coverage per mutation per week

Mutations are color-coded by uniqueness:
  Blue  (unique)        — not in any panel or cowwid variant → real signal
  Orange (shared)       — present in a known variant (e.g. XFG not in panel)
  Gray  (shared panel)  — present in a panel variant

Helps the user understand: is this signal truly new, or just a known
variant the user didn't include in their panel?
"""
import asyncio
from datetime import datetime
from typing import List, Dict, Set, Optional
import pandas as pd
import plotly.graph_objects as go
import streamlit as st


def _build_queries(mutations: List[str]) -> List[Dict[str, str]]:
    """
    Convert '{pos}{alt}' mutation strings into queriesOverTime query objects.
    e.g. '23018T' → {countQuery: 'main:23018T', coverageQuery: '!main:23018N', ...}
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
    """Fetch mutation frequencies over time. Returns DataFrame with columns:
    mutation, dateFrom, dateTo, frequency, count, coverage."""
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


def _classify_mutations(
    mutations: List[str],
    lineage_sig: Optional[Set[str]],
    panel_variants: List[str],
    all_lineage_signatures: Dict[str, Set[str]],
    cowwid_signatures: Dict[str, Set[str]],
) -> Dict[str, str]:
    """
    Classify each mutation as:
    - 'unique':        not in any panel or cowwid variant signature
    - 'shared_cowwid': in a cowwid variant (e.g. XFG not in panel)
    - 'shared_panel':  in a panel variant

    If lineage_sig is None (e.g. for missing cowwid variants where we only
    have observed_mutations), classify only from cowwid/panel membership.
    """
    panel_union = set().union(
        *(all_lineage_signatures.get(p, set()) for p in panel_variants)
    ) if all_lineage_signatures else set()

    cowwid_union = set().union(
        *cowwid_signatures.values()
    ) if cowwid_signatures else set()

    result = {}
    for mut in mutations:
        if mut in panel_union:
            result[mut] = 'shared_panel'
        elif mut in cowwid_union:
            result[mut] = 'shared_cowwid'
        else:
            result[mut] = 'unique'
    return result


def _make_figure(
    pivot: pd.DataFrame,
    hover: pd.DataFrame,
    cov: pd.DataFrame,
    mut_classes: Dict[str, str],
    title: str,
) -> go.Figure:
    """
    Build a Plotly heatmap figure with color-coded mutation rows.
    Unique mutations use Blues colorscale, shared use Oranges, panel use Greys.
    Renders as three overlapping traces with masks so each row uses its own scale.
    """
    col_labels = [
        datetime.strptime(c, "%Y-%m-%d").strftime("%b %d")
        for c in pivot.columns
    ]

    # sort rows: unique first, then shared_cowwid, then shared_panel
    priority = {'unique': 0, 'shared_cowwid': 1, 'shared_panel': 2}
    sorted_muts = sorted(pivot.index, key=lambda m: (priority.get(mut_classes.get(m, 'shared_panel'), 2), m))

    # build per-group traces
    groups = [
        ('unique',       'Blues',   'unique to lineage'),
        ('shared_cowwid','Oranges', 'shared with known variant'),
        ('shared_panel', 'Greys',   'shared with panel variant'),
    ]

    fig = go.Figure()

    for group_key, colorscale, label in groups:
        group_muts = [m for m in sorted_muts if mut_classes.get(m) == group_key]
        if not group_muts:
            continue

        z = pivot.loc[group_muts].values.tolist()
        h = hover.loc[group_muts].values if group_muts[0] in hover.index else None
        v = cov.loc[group_muts].values if group_muts[0] in cov.index else None

        hover_text = []
        for i, mut in enumerate(group_muts):
            row_hover = []
            for j, col in enumerate(pivot.columns):
                f = pivot.loc[mut, col]
                c = hover.loc[mut, col] if mut in hover.index else 0
                cv = cov.loc[mut, col] if mut in cov.index else 0
                cls = mut_classes.get(mut, '')
                cls_label = {'unique': '🔵 unique', 'shared_cowwid': '🟠 shared (known variant)', 'shared_panel': '⚫ shared (panel)'}.get(cls, '')
                row_hover.append(
                    f"<b>{mut}</b> {cls_label}<br>"
                    f"freq: {f:.1%}<br>"
                    f"{int(c):,} / {int(cv):,} reads"
                )
            hover_text.append(row_hover)

        zmax = max((max(max(r) for r in z), 0.01))

        fig.add_trace(go.Heatmap(
            z=z,
            x=col_labels,
            y=group_muts,
            text=hover_text,
            hovertemplate="%{text}<extra></extra>",
            colorscale=colorscale,
            zmin=0,
            zmax=zmax,
            showscale=(group_key == 'unique'),  # only show scale for unique
            colorbar=dict(title="freq", tickformat=".0%", x=1.02) if group_key == 'unique' else None,
            name=label,
        ))

    fig.update_layout(
        title=dict(text=title, font=dict(size=14)),
        height=max(200, 32 * len(sorted_muts) + 90),
        margin=dict(t=50, b=40, l=100, r=60),
        template="plotly_white",
        xaxis=dict(side="bottom", title=None),
        yaxis=dict(autorange="reversed", title=None, tickfont=dict(family="monospace", size=11)),
        legend=dict(orientation="h", y=-0.12, x=0, font=dict(size=11)),
    )
    return fig


def render_scanner_heatmap(
    variant: str,
    mutations: List[str],
    client,
    location: str,
    date_range: tuple,
    max_mutations: int = 30,
    panel_variants: Optional[List[str]] = None,
    all_lineage_signatures: Optional[Dict[str, Set[str]]] = None,
    cowwid_signatures: Optional[Dict[str, Set[str]]] = None,
    lineage_sig: Optional[Set[str]] = None,
) -> None:
    """
    Render a mutation frequency heatmap for a scanner candidate.

    Works for both:
    - missing_from_panel variants (cowwid variants not in panel)
    - emerging_sublineage candidates

    Args:
        variant:               Variant/lineage name e.g. "RF.5.2"
        mutations:             List of "{pos}{alt}" observed mutations
        client:                WiseLoculusLapis instance
        location:              Location name e.g. "Zürich (ZH)"
        date_range:            (start_datetime, end_datetime) tuple
        max_mutations:         Cap shown mutations (top by position)
        panel_variants:        Currently selected panel variant names
        all_lineage_signatures: All pango lineage signatures (for classification)
        cowwid_signatures:     All cowwid variant signatures
        lineage_sig:           Full signature of this lineage (for classification)
    """
    if not mutations:
        st.caption(f"No discriminating mutations found for {variant}.")
        return

    shown = sorted(mutations)[:max_mutations]
    if len(mutations) > max_mutations:
        st.caption(f"Showing {max_mutations} of {len(mutations)} observed mutations.")

    # classify mutations if signatures available
    mut_classes = {}
    if all_lineage_signatures is not None and cowwid_signatures is not None and panel_variants:
        mut_classes = _classify_mutations(
            shown, lineage_sig, panel_variants, all_lineage_signatures, cowwid_signatures
        )
    else:
        # fallback: all unknown
        mut_classes = {m: 'unique' for m in shown}

    n_unique = sum(1 for c in mut_classes.values() if c == 'unique')
    n_shared = sum(1 for c in mut_classes.values() if c != 'unique')

    # legend pills
    legend_parts = []
    if n_unique:
        legend_parts.append(
            f'<span style="background:#dbeafe;color:#1e40af;padding:2px 8px;'
            f'border-radius:10px;font-size:11px;margin-right:6px">'
            f'🔵 {n_unique} unique</span>'
        )
    if n_shared:
        legend_parts.append(
            f'<span style="background:#ffedd5;color:#9a3412;padding:2px 8px;'
            f'border-radius:10px;font-size:11px;margin-right:6px">'
            f'🟠 {n_shared} shared with known variant</span>'
        )
    if legend_parts:
        st.markdown(" ".join(legend_parts), unsafe_allow_html=True)

    with st.spinner(f"Fetching {variant} signal over time…"):
        df = _fetch_frequencies(client, location, date_range, shown)

    if df.empty:
        st.caption("No frequency data returned.")
        return

    pivot = df.pivot(index="mutation", columns="dateFrom", values="frequency")
    hover_df = df.pivot(index="mutation", columns="dateFrom", values="count")
    cov_df = df.pivot(index="mutation", columns="dateFrom", values="coverage")

    # auto-verdict based on unique mutations
    unique_muts_in_pivot = [m for m in shown if mut_classes.get(m) == 'unique' and m in pivot.index]
    if unique_muts_in_pivot:
        unique_pivot = pivot.loc[unique_muts_in_pivot]
        max_unique_freq = unique_pivot.values.max()
        # count unique mutations rising in last 2 weeks
        if len(pivot.columns) >= 2:
            last2 = unique_pivot.iloc[:, -2:].values
            n_rising = sum(1 for row in last2 if row[-1] > row[0])
        else:
            n_rising = 0
            max_unique_freq = 0

        if max_unique_freq > 0.05 and n_rising >= 2:
            st.success(f"**Real signal** — {n_rising} unique mutations rising together (max {max_unique_freq:.1%})")
        elif max_unique_freq > 0.005:
            st.warning(f"**Weak signal** — unique mutations present but low ({max_unique_freq:.1%} max). Monitor over time.")
        else:
            st.info("**Likely noise** — unique mutations not rising. Signal may be from a related known variant.")

    title = f"{variant} — mutation frequencies by week in {location}"
    fig = _make_figure(pivot, hover_df, cov_df, mut_classes, title)
    st.plotly_chart(fig, use_container_width=True)

    st.caption(
        "Frequency = reads with mutation / reads covering that position. "
        "Blue = mutations unique to this lineage (true signal). "
        "Orange = shared with a known variant not in your panel. "
        "Only mutations observed in unexplained co-occurrence patterns are shown."
    )