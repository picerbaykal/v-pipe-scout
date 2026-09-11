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
    """Fetch mutation frequencies over time, cached in session_state so that
    Streamlit re-runs (expanding, slider, city switch) don't re-query LAPIS
    for the same location/dates/mutations."""
    try:
        import streamlit as _st
        _key = (location, str(date_range[0]), str(date_range[1]),
                tuple(sorted(mutations)))
        _cache = _st.session_state.setdefault("_scanner_freq_cache", {})
        if _key in _cache:
            return _cache[_key]
        _df = _fetch_frequencies_uncached(client, location, date_range, mutations)
        _cache[_key] = _df
        return _df
    except Exception:
        # if streamlit/session unavailable (e.g. worker), fetch directly
        return _fetch_frequencies_uncached(client, location, date_range, mutations)


def _fetch_frequencies_uncached(
    client,
    location: str,
    date_range: tuple,
    mutations: List[str],
) -> pd.DataFrame:
    """Fetch mutation frequencies over time via /sample/aggregated — the SAME
    endpoint the co-occurrence pipeline uses, so the heatmap and the scanner
    see identical data. For each real sampling date we query all positions at
    once; frequency = reads with the mutation / reads covering that position
    (non-N). Returns columns: mutation, dateFrom, dateTo, frequency, count,
    coverage."""
    import re as _re
    positions = sorted({int(_re.match(r'^(\d+)', m).group(1))
                        for m in mutations if _re.match(r'^(\d+)', m)})
    alt_of = {int(_re.match(r'^(\d+)', m).group(1)): m[-1]
              for m in mutations if _re.match(r'^(\d+)', m)}

    async def _run():
        import aiohttp
        start, end = date_range
        dates = await client._get_sampling_dates(location, (start, end))
        dates = sorted(dates)
        out = []
        async with aiohttp.ClientSession() as session:
            for d in dates:
                rows = await client._fetch_cooccurrence_for_date(
                    session, location, d, positions)
                # per position: total covered (non-N) and mutation count
                cov = {p: 0 for p in positions}
                mut = {p: 0 for p in positions}
                for row in (rows or []):
                    c = row.get('count', 0)
                    for p in positions:
                        b = row.get(f'[{p}]')
                        if b is not None and b != 'N':
                            cov[p] += c
                            if b == alt_of[p]:
                                mut[p] += c
                for m in mutations:
                    mm = _re.match(r'^(\d+)', m)
                    if not mm:
                        continue
                    p = int(mm.group(1))
                    freq = mut[p] / cov[p] if cov[p] > 0 else 0.0
                    out.append({
                        'mutation': m, 'dateFrom': d, 'dateTo': d,
                        'frequency': freq, 'count': mut[p], 'coverage': cov[p],
                    })
        return out

    rows = asyncio.run(_run())
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
    member_blocks: Optional[List[dict]] = None,
    shared_mutations: Optional[List[str]] = None,
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

    # member-block guide: which mutations distinguish which clade member
    if member_blocks:
        _guide = []
        if shared_mutations:
            _sm = ", ".join(shared_mutations[:4])
            _guide.append(f"<b>shared (clade):</b> {_sm}")
        for _b in member_blocks[:6]:
            _dm = ", ".join(_b.get("discriminating", [])[:3])
            if _dm:
                _guide.append(f"<b>{_b['member']}:</b> {_dm}")
        if _guide:
            st.markdown(
                "<div style='font-size:11px;color:#6b7280;margin:4px 0;'>"
                "Rows below, grouped by member — a member is present when its "
                "own row(s) light up alongside the shared rows.<br>"
                + " &nbsp;·&nbsp; ".join(_guide) + "</div>",
                unsafe_allow_html=True,
            )

    title = f"{variant} — mutation frequencies by week in {location}"
    fig = _make_figure(pivot, hover_df, cov_df, mut_classes, title)
    st.plotly_chart(fig, use_container_width=True)

    st.caption(
        "Frequency = reads with mutation / reads covering that position. "
        "Blue = mutations unique to this lineage (true signal). "
        "Orange = shared with a known variant not in your panel. "
        "Only mutations observed in unexplained co-occurrence patterns are shown."
    )

def render_clade_heatmap(
    clade_node: str,
    shared_mutations: list,
    member_blocks: list,
    client,
    location: str,
    date_range: tuple,
    max_members: int = 12,
    max_muts_per_block: int = 25,
    reads_threshold: int = 0,
) -> None:
    """Signal-over-time for a clade: one small heatmap per family block.

    The shared block confirms the clade is present. Each member/family block
    below shows that family's discriminating co-occurrence group — when its
    block lights up alongside shared, that family is driving the signal.
    Each block is a separate small heatmap with its own title, so labels are
    always clear and different-amplicon groups aren't misread as one.
    """
    import numpy as np
    import plotly.graph_objects as go
    from datetime import datetime

    # blocks to render: split by the reads threshold. Above-threshold blocks
    # show their full heatmap; below-threshold ones collapse to a dimmed line
    # (still listed, nothing removed). shared block always shown.
    render_blocks = []
    if shared_mutations:
        render_blocks.append(("shared (clade)",
                              "present in all members — confirms the clade",
                              list(shared_mutations)[:max_muts_per_block]))
    below = []
    for b in member_blocks[:max_members]:
        name = b.get("member", "?")
        mc = b.get("member_count", 1)
        reads = b.get("reads", 0)
        subtitle = (f"{mc} lineages · {reads:,} reads" if mc > 1
                    else f"{reads:,} reads")
        if reads >= reads_threshold:
            render_blocks.append((name, subtitle,
                                  list(b.get("discriminating", []))[:max_muts_per_block]))
        else:
            below.append((name, reads))

    if not render_blocks and not below:
        st.caption("No co-occurrence groups to display for this clade.")
        return

    # fetch all muts once
    all_muts = list(dict.fromkeys(
        m for _, _, muts in render_blocks for m in muts))
    df = None
    if all_muts:
        with st.spinner(f"Fetching {clade_node} signal over time…"):
            df = _fetch_frequencies(client, location, date_range, all_muts)
    if df is None or df.empty:
        if below:
            _lines = " · ".join(f"{n} ({r:,})" for n, r in below)
            st.caption(f"All regions below the threshold: {_lines}")
        else:
            st.caption("No frequency data returned.")
        return
    df = df.drop_duplicates(subset=["mutation", "dateFrom"], keep="first")
    freq = df.pivot(index="mutation", columns="dateFrom", values="frequency")
    cov = df.pivot(index="mutation", columns="dateFrom", values="coverage")
    cols = list(freq.columns)
    col_labels = [datetime.strptime(c, "%Y-%m-%d").strftime("%b %d") for c in cols]

    st.caption(
        "Each block is a co-occurrence group for this clade in one amplicon region. "
        "Row labels show [number of lineages carrying that mutation]: a low number "
        "with ★ = discriminating (few variants have it); a high number = backbone "
        "(broadly shared, not variant-specific). The ★ rows are the ones that "
        "confirm the variant — if they are dark the variant is present; if only the "
        "high-number backbone rows are dark, it is just shared mutations. "
        "Hatched = no coverage that week (not absence)."
    )

    # carrier-count lookup across all blocks: how many lineages carry each mut.
    # Used to annotate rows so the user sees which are discriminating (few
    # carriers) vs backbone (many). Discriminating threshold ~30 lineages.
    _mut_car = {}
    for _b in member_blocks:
        _mut_car.update(_b.get("mut_carriers", {}))

    def _row_label(m):
        n = _mut_car.get(m)
        if n is None:
            return m
        if n <= 30:
            return f"{m} [{n}]★"          # discriminating (star marks it)
        return f"{m} [{n}]"                # backbone (high count)

    def _one_block(title, subtitle, muts, key):
        rows = [m for m in muts if m in freq.index]
        if not rows:
            return
        row_labels = [_row_label(m) for m in rows]
        z, hatch = [], []
        for m in rows:
            frow, hrow = [], []
            for c in cols:
                fv = freq.loc[m, c]
                cvv = cov.loc[m, c] if m in cov.index else 0
                if cvv is None or cvv == 0 or (isinstance(cvv, float) and np.isnan(cvv)):
                    frow.append(np.nan); hrow.append(True)
                else:
                    frow.append(float(fv) if fv == fv else 0.0); hrow.append(False)
            z.append(frow); hatch.append(hrow)

        fig = go.Figure()
        fig.add_trace(go.Heatmap(
            z=z, x=col_labels, y=row_labels, colorscale="Blues", zmin=0, zmax=1,
            hoverongaps=False, showscale=False, xgap=1, ygap=1,
        ))
        hz = [[1 if hatch[i][j] else np.nan for j in range(len(cols))]
              for i in range(len(rows))]
        fig.add_trace(go.Heatmap(
            z=hz, x=col_labels, y=row_labels,
            colorscale=[[0, "#e5e7eb"], [1, "#e5e7eb"]],
            showscale=False, hoverongaps=False,
            hovertemplate="no coverage<extra></extra>", xgap=1, ygap=1,
        ))
        fig.update_layout(
            height=max(90, 24 * len(rows) + 40),
            margin=dict(l=90, r=10, t=6, b=24),
            template="plotly_white",
            yaxis=dict(autorange="reversed", tickfont=dict(size=10)),
            xaxis=dict(tickfont=dict(size=9)),
        )
        st.markdown(
            f"<div style='font-size:13px;font-weight:500;margin:8px 0 0;'>{title}"
            f"<span style='font-size:11px;color:#6b7280;font-weight:400;'> · {subtitle}</span></div>",
            unsafe_allow_html=True,
        )
        st.plotly_chart(fig, use_container_width=True, key=key)

    for i, (title, subtitle, muts) in enumerate(render_blocks):
        _one_block(title, subtitle, muts,
                   key=f"clade_hm_{clade_node}_{location}_{i}")

    # dimmed list of below-threshold regions — kept, not removed
    if below:
        _lines = "<br>".join(
            f"<span style='color:#9ca3af;'>{n} · {r:,} reads · below threshold</span>"
            for n, r in below)
        st.markdown(
            f"<div style='font-size:12px;margin-top:6px;'>{_lines}</div>",
            unsafe_allow_html=True,
        )