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

    Within each block the discriminating (★) rows are shown FIRST, in the
    Blues scale; backbone (high-carrier) rows follow in a muted Greys scale
    so they read as context, not signal. A block with no ★ row is flagged
    'backbone only'. Cell tooltips are labelled with date, mutation, carrier
    count, discriminating-vs-backbone and frequency as a %% (count/coverage).
    """
    import numpy as np
    import re as _re
    import plotly.graph_objects as go
    from datetime import datetime

    _STAR_MAX = 30   # carriers <= this → discriminating (★)

    # carrier-count lookup across all blocks
    _mut_car = {}
    for _b in member_blocks:
        _mut_car.update(_b.get("mut_carriers", {}))

    def _is_star(m):
        n = _mut_car.get(m)
        return n is not None and n <= _STAR_MAX

    def _pos(m):
        mm = _re.match(r"^(\d+)", m)
        return int(mm.group(1)) if mm else 0

    # ── block selection: DEFAULT shows only regions that contain a
    # discriminating (★) mutation — those are what confirm the variant.
    # Backbone-only regions are hidden behind a per-finding toggle (nothing is
    # dropped permanently). reads_threshold is kept in the signature for
    # backward compatibility but is no longer used to gate blocks.
    def _blk_has_star(b):
        return any(_is_star(m) for m in b.get("discriminating", []))

    def _blk_min_carrier(b):
        # lowest carrier count among the block's discriminating muts
        cs = [_mut_car.get(m, 10**9) for m in b.get("discriminating", [])
              if _mut_car.get(m) is not None]
        return min(cs) if cs else 10**9

    _star_blocks = [b for b in member_blocks if _blk_has_star(b)]
    _back_blocks = [b for b in member_blocks if not _blk_has_star(b)]

    # Group blocks by family (family_root): all blocks of one family are shown
    # together, families ordered by their MOST discriminating block (lowest
    # carrier count), and within a family by genome position — so the same
    # family's regions read left-to-right instead of being scattered.
    def _blk_fam(b):
        return b.get("family_root") or b.get("member", "").split(" @ ")[0]

    def _blk_startpos(b):
        ps = [_pos(m) for m in b.get("discriminating", [])]
        return min(ps) if ps else 0

    _fam_rank = {}
    for _b in _star_blocks:
        _f = _blk_fam(_b)
        _c = _blk_min_carrier(_b)
        if _f not in _fam_rank or _c < _fam_rank[_f]:
            _fam_rank[_f] = _c
    _star_blocks.sort(key=lambda b: (_fam_rank.get(_blk_fam(b), 10**9),
                                     _blk_fam(b), _blk_startpos(b)))
    # backbone-only regions: same family grouping, then position
    _bfam_rank = {}
    for _b in _back_blocks:
        _f = _blk_fam(_b)
        _r = -_b.get("reads", 0)
        if _f not in _bfam_rank or _r < _bfam_rank[_f]:
            _bfam_rank[_f] = _r
    _back_blocks.sort(key=lambda b: (_bfam_rank.get(_blk_fam(b), 0),
                                     _blk_fam(b), _blk_startpos(b)))

    _toggle_key = f"clade_showall_{clade_node}_{location}"
    _show_all = st.session_state.get(_toggle_key, False)
    _SAFETY_CAP = 40   # only to guard a pathological finding, not a real limit

    # DEFAULT: ALL regions that have a discriminating (★) mutation (these also
    # contain backbone rows — that's fine, the ★ is what matters). Backbone-ONLY
    # regions (no ★ at all) are appended only when the user asks. Nothing is
    # capped away in practice; _SAFETY_CAP just bounds a runaway finding.
    render_blocks = []
    if shared_mutations:
        render_blocks.append(("shared (clade)",
                              "present in all members — confirms the clade",
                              list(shared_mutations)[:max_muts_per_block]))
    _selected = list(_star_blocks)
    if _show_all:
        _selected += _back_blocks
    for b in _selected[:_SAFETY_CAP]:
        name = b.get("member", "?")
        mc = b.get("member_count", 1)
        reads = b.get("reads", 0)
        subtitle = (f"{mc} lineages · {reads:,} reads" if mc > 1
                    else f"{reads:,} reads")
        render_blocks.append((name, subtitle,
                              list(b.get("discriminating", []))[:max_muts_per_block]))
    below = []  # kept for the "no data" branch below

    # header line: what's shown, what's hidden, and the toggle — at the TOP so it
    # isn't lost below a long stack of grids.
    _n_star = len(_star_blocks)
    _n_back = len(_back_blocks)
    _hdr = (f"Showing {_n_star} region(s) with a discriminating (\u2605) mutation."
            if _n_star else "No regions have a discriminating (\u2605) mutation.")
    if _n_back:
        _hdr += (f"  {_n_back} backbone-only region(s) "
                 + ("shown below." if _show_all else "hidden."))
    st.markdown(f"<div style='font-size:12px;color:#6b7280;margin:2px 0 4px;'>{_hdr}"
                f"</div>", unsafe_allow_html=True)
    if _n_back:
        _lbl = (f"Hide backbone-only regions ({_n_back})" if _show_all
                else f"Show backbone-only regions (+{_n_back})")
        if st.button(_lbl, key=f"btn_top_{_toggle_key}"):
            st.session_state[_toggle_key] = not _show_all
            st.rerun()

    # all-backbone finding (no ★ region anywhere): with backbone hidden there is
    # nothing to draw — say so and stop (the toggle above reveals them).
    if not _star_blocks and not shared_mutations and not _show_all:
        st.caption("Every co-occurrence group here is backbone (shared) only, so "
                   "co-occurrence can't confirm a specific variant. Use the button "
                   "above to inspect the backbone regions.")
        return

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
    cnt = df.pivot(index="mutation", columns="dateFrom", values="count")
    cols = list(freq.columns)
    col_labels = [datetime.strptime(c, "%Y-%m-%d").strftime("%b %d") for c in cols]

    st.caption(
        "Each block is a co-occurrence group for this clade in one amplicon region. "
        "Rows are labelled [number of lineages carrying that mutation]: a low number "
        "with ★ = discriminating (few variants have it), shown first; a high number = "
        "backbone (broadly shared, context only). The ★ rows are the ones that confirm "
        "the variant — if they are dark on recent dates the variant is present; if only "
        "the high-number backbone rows are dark, it is just shared mutations. "
        "Hatched = no coverage that week (not absence)."
    )

    def _row_label(m):
        n = _mut_car.get(m)
        if n is None:
            return m
        if n <= _STAR_MAX:
            # discriminating — bold blue so ★ rows stand out even on one scale
            return f"<b><span style='color:#185FA5'>{m} [{n}] \u2605</span></b>"
        # backbone — muted grey label
        return f"<span style='color:#9ca3af'>{m} [{n}]</span>"

    def _one_block(title, subtitle, muts, key):
        rows = [m for m in muts if m in freq.index]
        if not rows:
            return
        # discriminating (★) rows first, backbone after; each sorted by position
        star_rows = sorted([m for m in rows if _is_star(m)], key=_pos)
        back_rows = sorted([m for m in rows if not _is_star(m)], key=_pos)
        ordered = star_rows + back_rows
        n_star = len(star_rows)
        row_labels = [_row_label(m) for m in ordered]

        z, hatch, htxt = [], [], []
        for m in ordered:
            frow, hrow, trow = [], [], []
            star = _is_star(m)
            ncar = _mut_car.get(m, "?")
            tag = "★ discriminating" if star else "backbone"
            for c in cols:
                fv = freq.loc[m, c]
                cvv = cov.loc[m, c] if m in cov.index else 0
                ct = cnt.loc[m, c] if m in cnt.index else 0
                dl = datetime.strptime(c, "%Y-%m-%d").strftime("%b %d")
                if cvv is None or cvv == 0 or (isinstance(cvv, float) and np.isnan(cvv)):
                    frow.append(np.nan); hrow.append(True)
                    trow.append(f"{dl} · {m} ({tag})<br>no coverage")
                else:
                    fval = float(fv) if fv == fv else 0.0
                    frow.append(fval); hrow.append(False)
                    trow.append(
                        f"{dl} · {m} [{ncar}] {tag}<br>"
                        f"{fval*100:.0f}% ({int(ct):,} / {int(cvv):,} reads)")
            z.append(frow); hatch.append(hrow); htxt.append(trow)

        fig = go.Figure()
        # single heatmap trace for ALL rows (★ first, then backbone) so rows
        # pack tightly with no seam. Two separate traces left a visible gap /
        # a lone detached backbone row. The ★-vs-backbone distinction is carried
        # by the row labels (★ marker + [carrier count]) and by ordering, not by
        # a second colour scale. One muted blue ramp for everything.
        fig.add_trace(go.Heatmap(
            z=z, x=col_labels, y=row_labels,
            text=htxt, hovertemplate="%{text}<extra></extra>",
            colorscale="Blues", zmin=0, zmax=1, hoverongaps=False,
            showscale=False, xgap=1, ygap=1))
        # hatch overlay for no-coverage cells
        hz = [[1 if hatch[i][j] else np.nan for j in range(len(cols))]
              for i in range(len(ordered))]
        fig.add_trace(go.Heatmap(
            z=hz, x=col_labels, y=row_labels,
            colorscale=[[0, "#e5e7eb"], [1, "#e5e7eb"]],
            showscale=False, hoverongaps=False,
            hovertemplate="no coverage<extra></extra>", xgap=1, ygap=1,
        ))
        fig.update_layout(
            height=max(90, 24 * len(ordered) + 40),
            margin=dict(l=110, r=10, t=6, b=24),
            template="plotly_white",
            yaxis=dict(autorange="reversed", tickfont=dict(size=10)),
            xaxis=dict(tickfont=dict(size=9)),
        )
        # divider between the ★ group (top) and the backbone group (below) so the
        # two groups read as intentional even though they share one trace/scale.
        # y is reversed, so the boundary sits below the last ★ row at n_star-0.5.
        if n_star and len(ordered) > n_star:
            fig.add_shape(
                type="line", xref="paper", yref="y",
                x0=0, x1=1, y0=n_star - 0.5, y1=n_star - 0.5,
                line=dict(color="#9ca3af", width=1, dash="dot"),
            )
        _flag = ("" if n_star else
                 " · ⚠ backbone only — shared mutations, not variant-specific")
        st.markdown(
            f"<div style='font-size:13px;font-weight:500;margin:8px 0 0;'>{title}"
            f"<span style='font-size:11px;color:#6b7280;font-weight:400;'> · "
            f"{subtitle}{_flag}</span></div>",
            unsafe_allow_html=True,
        )
        st.plotly_chart(fig, use_container_width=True, key=key)

    for i, (title, subtitle, muts) in enumerate(render_blocks):
        _one_block(title, subtitle, muts,
                   key=f"clade_hm_{clade_node}_{location}_{i}")