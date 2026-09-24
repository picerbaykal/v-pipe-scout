"""Scanner heatmap component.

Two views of a scanner finding's signal over time:

1. DEFAULT — read-level HAPLOTYPES.  For each co-occurrence block that carries a
   discriminating (★) mutation, we show the actual base COMBINATIONS that
   co-occur on reads at that block's positions, ordered by read count.  A read
   that carries all of the block's mutations is the VARIANT haplotype; one that
   carries a ★ but not the rest is a partial★ (sublineage / dropout); one that
   carries none of the ★ is a co-circulating relative (shares the backbone, not
   the variant).  Striking rows are shown, the long tail is folded.

2. DETAIL (behind a button) — the older per-MUTATION heatmap: one row per
   mutation showing its individual frequency trajectory and coverage.  Useful
   for reading coverage/dropout, but noisy (25+ marginal rows), so it is opt-in.

Both read the SAME /sample/aggregated endpoint the co-occurrence pipeline uses.
"""
import asyncio
from datetime import datetime
from collections import defaultdict
from typing import List, Dict, Set, Optional, Tuple
import re as _re

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

_STAR_MAX = 30          # carriers <= this → discriminating (★)
_MAX_POS = 5            # cap positions per haplotype block (keeps coverage high)
_TOP_N = 5             # dominant combos always shown by default
_MIN_READS = 50         # combos below this fold into "other"
_UNCOVERED = {"N", "-"}


def _pos(m: str) -> int:
    mm = _re.match(r"^(\d+)", m)
    return int(mm.group(1)) if mm else 0


# ─────────────────────────────────────────────────────────────────────────────
# per-MUTATION frequency fetch  (used by the DETAIL view + render_scanner_heatmap)
# ─────────────────────────────────────────────────────────────────────────────
def _build_queries(mutations: List[str]) -> List[Dict[str, str]]:
    queries = []
    for mut in mutations:
        pos = mut[:-1]
        queries.append({
            "displayLabel": mut,
            "countQuery": f"main:{pos}{mut[-1]}",
            "coverageQuery": f"!main:{pos}N",
        })
    return queries


def _fetch_frequencies(client, location, date_range, mutations) -> pd.DataFrame:
    """Cached wrapper around the per-mutation fetch."""
    try:
        _key = (location, str(date_range[0]), str(date_range[1]),
                tuple(sorted(mutations)))
        _cache = st.session_state.setdefault("_scanner_freq_cache", {})
        if _key in _cache:
            return _cache[_key]
        _df = _fetch_frequencies_uncached(client, location, date_range, mutations)
        _cache[_key] = _df
        return _df
    except Exception:
        return _fetch_frequencies_uncached(client, location, date_range, mutations)


def _fetch_frequencies_uncached(client, location, date_range, mutations) -> pd.DataFrame:
    positions = sorted({int(_re.match(r'^(\d+)', m).group(1))
                        for m in mutations if _re.match(r'^(\d+)', m)})
    alt_of = {int(_re.match(r'^(\d+)', m).group(1)): m[-1]
              for m in mutations if _re.match(r'^(\d+)', m)}

    async def _run():
        import aiohttp
        start, end = date_range
        dates = sorted(await client._get_sampling_dates(location, (start, end)))
        out = []
        async with aiohttp.ClientSession() as session:
            for d in dates:
                rows = await client._fetch_cooccurrence_for_date(
                    session, location, d, positions)
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
                    out.append({'mutation': m, 'dateFrom': d, 'dateTo': d,
                                'frequency': freq, 'count': mut[p], 'coverage': cov[p]})
        return out

    return pd.DataFrame(asyncio.run(_run()))


# ─────────────────────────────────────────────────────────────────────────────
# per-READ haplotype fetch  (used by the DEFAULT view)
# ─────────────────────────────────────────────────────────────────────────────
def _fetch_haplotypes(client, location, date_range, positions: Tuple[int, ...]):
    """Return (dates, per_date, covered) for a block's positions.

    per_date : {date: {base_tuple: reads}}   base_tuple aligned to `positions`
    covered  : {date: reads fully covering ALL positions}
    Only reads that are non-N at every position count (a partial-coverage read
    can't be assigned to a haplotype).  Cached per (location, dates, positions).
    """
    positions = tuple(positions)
    try:
        _key = ("hap", location, str(date_range[0]), str(date_range[1]), positions)
        _cache = st.session_state.setdefault("_scanner_hap_cache", {})
        if _key in _cache:
            return _cache[_key]
    except Exception:
        _cache, _key = None, None

    async def _run():
        import aiohttp
        start, end = date_range
        dates = sorted(await client._get_sampling_dates(location, (start, end)))
        per_date, covered = {}, {}
        async with aiohttp.ClientSession() as session:
            for d in dates:
                rows = await client._fetch_cooccurrence_for_date(
                    session, location, d, list(positions))
                cc, tot = defaultdict(int), 0
                for row in (rows or []):
                    cnt = int(row.get("count", 0) or 0)
                    if cnt <= 0:
                        continue
                    bases = tuple(row.get(f"[{p}]", "N") for p in positions)
                    if any(b in _UNCOVERED for b in bases):
                        continue
                    cc[bases] += cnt
                    tot += cnt
                per_date[d] = dict(cc)
                covered[d] = tot
        return dates, per_date, covered

    result = asyncio.run(_run())
    if _cache is not None and _key is not None:
        _cache[_key] = result
    return result


def _hap_positions(blk_muts: List[str], mut_car: Dict[str, int]):
    """Choose <= _MAX_POS positions for a block: all ★ positions, then fill with
    the lowest-carrier (most informative) backbone positions.  Returns
    (positions, mut_base{pos->alt}, star_pos set)."""
    pos_alt, star_pos = {}, set()
    for m in blk_muts:
        p = _pos(m)
        pos_alt[p] = m[-1]
        n = mut_car.get(m)
        if n is not None and n <= _STAR_MAX:
            star_pos.add(p)
    all_pos = sorted(pos_alt)
    if len(all_pos) <= _MAX_POS:
        keep = all_pos
    else:
        stars = sorted(star_pos)[:_MAX_POS]
        back = sorted((p for p in all_pos if p not in star_pos),
                      key=lambda p: mut_car.get(f"{p}{pos_alt[p]}", 10**9))
        keep = sorted(stars + back[:max(0, _MAX_POS - len(stars))])
    return keep, {p: pos_alt[p] for p in keep}, {p for p in keep if p in star_pos}


def _one_haplotype_block(positions, mut_base, star_pos, dates, per_date, covered,
                         title, subtitle, key):
    """Render one block's haplotype table as a compact heatmap.

    Rows = base combinations (striking set), ordered by reads.  The row label
    carries the combo, tag (VARIANT / partial★ / relative), read count and %.
    Cells = that combo's share of covered reads per date (trajectory).  Minor
    haplotypes are folded into a caption; a ▸ expander reveals them.
    """
    npos = len(positions)
    totals = defaultdict(int)
    for d in dates:
        for combo, n in per_date.get(d, {}).items():
            totals[combo] += n
    grand = sum(totals.values())
    if grand == 0:
        st.markdown(f"<div style='font-size:13px;font-weight:500;margin:8px 0 0;'>{title}"
                    f"<span style='font-size:11px;color:#6b7280;font-weight:400;'> · "
                    f"{subtitle}</span></div>", unsafe_allow_html=True)
        st.caption("No reads fully cover this block's positions in the window.")
        return

    def _classify(combo):
        n_mut = sum(1 for i, p in enumerate(positions) if combo[i] == mut_base[p])
        has_star = any(combo[i] == mut_base[p]
                       for i, p in enumerate(positions) if p in star_pos)
        is_var = n_mut == npos
        return n_mut, has_star, is_var

    ranked = sorted(totals.items(), key=lambda kv: -kv[1])
    keep, minor = [], []
    for i, (combo, n) in enumerate(ranked):
        n_mut, has_star, is_var = _classify(combo)
        striking = is_var or (n >= _MIN_READS and (i < _TOP_N or has_star))
        (keep if striking else minor).append((combo, n, n_mut, has_star, is_var))

    # header
    _pos_str = " · ".join(f"{p}{'★' if p in star_pos else ''}" for p in positions)
    st.markdown(
        f"<div style='font-size:13px;font-weight:600;margin:10px 0 0;'>{title}"
        f"<span style='font-size:11px;color:#6b7280;font-weight:400;'> · "
        f"{subtitle} · positions {_pos_str}</span></div>",
        unsafe_allow_html=True,
    )

    col_labels = [datetime.strptime(d, "%Y-%m-%d").strftime("%b %d") for d in dates]

    def _row_label(combo, n, n_mut, has_star, is_var):
        pct = n / grand * 100
        combo_s = "".join(combo)
        star = " ★" if has_star else ""
        if is_var:
            return (f"<b><span style='color:#185FA5'>{combo_s}{star} · "
                    f"{n:,} ({pct:.1f}%) VARIANT</span></b>")
        if has_star:
            return (f"<span style='color:#b45309'>{combo_s}{star} · "
                    f"{n:,} ({pct:.1f}%) partial</span>")
        return (f"<span style='color:#9ca3af'>{combo_s} · "
                f"{n:,} ({pct:.1f}%) relative</span>")

    row_labels, z, hz, htxt = [], [], [], []
    for combo, n, n_mut, has_star, is_var in keep:
        row_labels.append(_row_label(combo, n, n_mut, has_star, is_var))
        frow, hrow, trow = [], [], []
        tag = "VARIANT" if is_var else ("partial★" if has_star else "relative")
        for d, dl in zip(dates, col_labels):
            tot = covered.get(d, 0)
            c = per_date.get(d, {}).get(combo, 0)
            if not tot:
                frow.append(np.nan); hrow.append(1)
                trow.append(f"{dl} · {''.join(combo)}<br>no coverage")
            else:
                share = c / tot
                frow.append(share); hrow.append(np.nan)
                trow.append(f"{dl} · {''.join(combo)} ({tag})<br>"
                            f"{share*100:.0f}% ({c:,} / {tot:,} reads)")
        z.append(frow); hz.append(hrow); htxt.append(trow)

    fig = go.Figure()
    fig.add_trace(go.Heatmap(
        z=z, x=col_labels, y=row_labels, text=htxt,
        hovertemplate="%{text}<extra></extra>", colorscale="Blues",
        zmin=0, zmax=1, hoverongaps=False, showscale=False, xgap=1, ygap=1))
    fig.add_trace(go.Heatmap(
        z=hz, x=col_labels, y=row_labels,
        colorscale=[[0, "#e5e7eb"], [1, "#e5e7eb"]], showscale=False,
        hoverongaps=False, hovertemplate="no coverage<extra></extra>",
        xgap=1, ygap=1))
    fig.update_layout(
        height=max(90, 26 * len(keep) + 40),
        margin=dict(l=190, r=10, t=6, b=24),
        template="plotly_white",
        yaxis=dict(autorange="reversed", tickfont=dict(size=10,
                   family="monospace")),
        xaxis=dict(tickfont=dict(size=9)),
    )
    st.plotly_chart(fig, use_container_width=True, key=key)

    if minor:
        msum = sum(r[1] for r in minor)
        with st.expander(f"▸ {len(minor)} minor haplotypes  "
                         f"({msum:,} reads, {msum/grand*100:.1f}%)"):
            _lines = []
            for combo, n, n_mut, has_star, is_var in minor[:60]:
                _lines.append(f"`{''.join(combo)}` — {n:,} ({n/grand*100:.2f}%)"
                              + (" ★" if has_star else ""))
            st.markdown("  \n".join(_lines))


def _render_haplotype_blocks(clade_node, star_blocks, mut_car, client,
                             location, date_range):
    """DEFAULT view: one haplotype table per ★-block, strongest (most reads)
    first and open, the rest collapsed.  Blocks are ordered by read count."""
    blocks = sorted(star_blocks, key=lambda b: -b.get("reads", 0))
    st.caption(
        "Read-level haplotypes — the actual base combinations that co-occur on "
        "reads. **VARIANT** (blue) carries every block mutation; **partial ★** "
        "(orange) carries the discriminating marker on an incomplete background "
        "(sublineage / dropout); **relative** (grey) shares the backbone but not "
        "the ★, so it is a co-circulating lineage, not this variant. Hatched = "
        "no coverage that week."
    )
    for i, b in enumerate(blocks):
        muts = list(b.get("discriminating", []))
        positions, mut_base, star_pos = _hap_positions(muts, mut_car)
        if len(positions) < 2:
            continue
        name = b.get("member", "?")
        mc = b.get("member_count", 1)
        reads = b.get("reads", 0)
        subtitle = (f"{mc} lineages · {reads:,} reads" if mc > 1
                    else f"{reads:,} reads")
        with st.spinner(f"Fetching {name} haplotypes…"):
            dates, per_date, covered = _fetch_haplotypes(
                client, location, date_range, positions)
        key = f"hap_{clade_node}_{location}_{i}"
        if i == 0:
            _one_haplotype_block(positions, mut_base, star_pos, dates, per_date,
                                 covered, name, subtitle, key)
        else:
            with st.expander(f"{name} · {subtitle}", expanded=False):
                _one_haplotype_block(positions, mut_base, star_pos, dates,
                                     per_date, covered, name, subtitle, key)


# ─────────────────────────────────────────────────────────────────────────────
# DETAIL view: the older per-mutation heatmap (opt-in, behind a button)
# ─────────────────────────────────────────────────────────────────────────────
def _render_per_mutation_blocks(clade_node, shared_mutations, member_blocks,
                                client, location, date_range,
                                max_muts_per_block=25):
    _mut_car = {}
    for _b in member_blocks:
        _mut_car.update(_b.get("mut_carriers", {}))

    def _is_star(m):
        n = _mut_car.get(m)
        return n is not None and n <= _STAR_MAX

    def _blk_has_star(b):
        return any(_is_star(m) for m in b.get("discriminating", []))

    def _blk_min_carrier(b):
        cs = [_mut_car.get(m, 10**9) for m in b.get("discriminating", [])
              if _mut_car.get(m) is not None]
        return min(cs) if cs else 10**9

    _star_blocks = [b for b in member_blocks if _blk_has_star(b)]
    _back_blocks = [b for b in member_blocks if not _blk_has_star(b)]
    _star_blocks.sort(key=lambda b: (_blk_min_carrier(b), -b.get("reads", 0)))
    _back_blocks.sort(key=lambda b: -b.get("reads", 0))

    _toggle_key = f"clade_showall_{clade_node}_{location}"
    _show_all = st.session_state.get(_toggle_key, False)
    _SAFETY_CAP = 40

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

    _n_back = len(_back_blocks)
    if _n_back:
        _lbl = (f"Hide backbone-only regions ({_n_back})" if _show_all
                else f"Show backbone-only regions (+{_n_back})")
        if st.button(_lbl, key=f"btn_pm_{_toggle_key}"):
            st.session_state[_toggle_key] = not _show_all
            st.rerun()

    if not render_blocks:
        st.caption("No co-occurrence groups to display for this clade.")
        return

    all_muts = list(dict.fromkeys(m for _, _, muts in render_blocks for m in muts))
    df = None
    if all_muts:
        with st.spinner(f"Fetching {clade_node} per-mutation detail…"):
            df = _fetch_frequencies(client, location, date_range, all_muts)
    if df is None or df.empty:
        st.caption("No frequency data returned.")
        return
    df = df.drop_duplicates(subset=["mutation", "dateFrom"], keep="first")
    freq = df.pivot(index="mutation", columns="dateFrom", values="frequency")
    cov = df.pivot(index="mutation", columns="dateFrom", values="coverage")
    cnt = df.pivot(index="mutation", columns="dateFrom", values="count")
    cols = list(freq.columns)
    col_labels = [datetime.strptime(c, "%Y-%m-%d").strftime("%b %d") for c in cols]

    st.caption(
        "Per-mutation frequency (count / coverage) by week. Rows labelled "
        "[carriers]: low + ★ = discriminating (shown first); high = backbone "
        "(context). Hatched = no coverage that week."
    )

    def _row_label(m):
        n = _mut_car.get(m)
        if n is None:
            return m
        if n <= _STAR_MAX:
            return f"<b><span style='color:#185FA5'>{m} [{n}] ★</span></b>"
        return f"<span style='color:#9ca3af'>{m} [{n}]</span>"

    def _one_block(title, subtitle, muts, key):
        rows = [m for m in muts if m in freq.index]
        if not rows:
            return
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
                    trow.append(f"{dl} · {m} [{ncar}] {tag}<br>"
                                f"{fval*100:.0f}% ({int(ct):,} / {int(cvv):,} reads)")
            z.append(frow); hatch.append(hrow); htxt.append(trow)

        fig = go.Figure()
        fig.add_trace(go.Heatmap(
            z=z, x=col_labels, y=row_labels, text=htxt,
            hovertemplate="%{text}<extra></extra>", colorscale="Blues",
            zmin=0, zmax=1, hoverongaps=False, showscale=False, xgap=1, ygap=1))
        hz = [[1 if hatch[i][j] else np.nan for j in range(len(cols))]
              for i in range(len(ordered))]
        fig.add_trace(go.Heatmap(
            z=hz, x=col_labels, y=row_labels,
            colorscale=[[0, "#e5e7eb"], [1, "#e5e7eb"]], showscale=False,
            hoverongaps=False, hovertemplate="no coverage<extra></extra>",
            xgap=1, ygap=1))
        fig.update_layout(
            height=max(90, 24 * len(ordered) + 40),
            margin=dict(l=110, r=10, t=6, b=24), template="plotly_white",
            yaxis=dict(autorange="reversed", tickfont=dict(size=10)),
            xaxis=dict(tickfont=dict(size=9)))
        if n_star and len(ordered) > n_star:
            fig.add_shape(type="line", xref="paper", yref="y", x0=0, x1=1,
                          y0=n_star - 0.5, y1=n_star - 0.5,
                          line=dict(color="#9ca3af", width=1, dash="dot"))
        _flag = ("" if n_star else
                 " · ⚠ backbone only — shared mutations, not variant-specific")
        st.markdown(
            f"<div style='font-size:13px;font-weight:500;margin:8px 0 0;'>{title}"
            f"<span style='font-size:11px;color:#6b7280;font-weight:400;'> · "
            f"{subtitle}{_flag}</span></div>", unsafe_allow_html=True)
        st.plotly_chart(fig, use_container_width=True, key=key)

    for i, (title, subtitle, muts) in enumerate(render_blocks):
        _one_block(title, subtitle, muts,
                   key=f"clade_pm_{clade_node}_{location}_{i}")


# ─────────────────────────────────────────────────────────────────────────────
# public entry point for a clade finding
# ─────────────────────────────────────────────────────────────────────────────
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
    """Signal-over-time for a clade.

    DEFAULT = read-level haplotype tables (one per ★-block, strongest open).
    A button reveals the older per-mutation heatmap for coverage detail.
    """
    _mut_car = {}
    for _b in member_blocks:
        _mut_car.update(_b.get("mut_carriers", {}))

    def _is_star(m):
        n = _mut_car.get(m)
        return n is not None and n <= _STAR_MAX

    def _blk_has_star(b):
        return any(_is_star(m) for m in b.get("discriminating", []))

    _star_blocks = [b for b in member_blocks if _blk_has_star(b)]

    if not _star_blocks:
        st.caption(
            "Every co-occurrence group here is backbone (shared) only, so no "
            "discriminating haplotype can be shown. See per-mutation detail below."
        )
    else:
        _render_haplotype_blocks(clade_node, _star_blocks, _mut_car, client,
                                 location, date_range)

    # ── DETAIL: old per-mutation heatmap, opt-in ────────────────────────────
    _detail_key = f"pm_detail_{clade_node}_{location}"
    if st.checkbox("Show per-mutation heatmap (coverage detail)",
                   key=_detail_key, value=False):
        _render_per_mutation_blocks(
            clade_node, shared_mutations, member_blocks, client, location,
            date_range, max_muts_per_block=max_muts_per_block)


# ─────────────────────────────────────────────────────────────────────────────
# unchanged: single-variant heatmap (used elsewhere)
# ─────────────────────────────────────────────────────────────────────────────
def _classify_mutations(mutations, lineage_sig, panel_variants,
                        all_lineage_signatures, cowwid_signatures) -> Dict[str, str]:
    panel_union = set().union(
        *(all_lineage_signatures.get(p, set()) for p in panel_variants)
    ) if all_lineage_signatures else set()
    cowwid_union = set().union(*cowwid_signatures.values()) if cowwid_signatures else set()
    result = {}
    for mut in mutations:
        if mut in panel_union:
            result[mut] = 'shared_panel'
        elif mut in cowwid_union:
            result[mut] = 'shared_cowwid'
        else:
            result[mut] = 'unique'
    return result


def _make_figure(pivot, hover, cov, mut_classes, title) -> go.Figure:
    col_labels = [datetime.strptime(c, "%Y-%m-%d").strftime("%b %d")
                  for c in pivot.columns]
    priority = {'unique': 0, 'shared_cowwid': 1, 'shared_panel': 2}
    sorted_muts = sorted(pivot.index, key=lambda m: (
        priority.get(mut_classes.get(m, 'shared_panel'), 2), m))
    groups = [('unique', 'Blues', 'unique to lineage'),
              ('shared_cowwid', 'Oranges', 'shared with known variant'),
              ('shared_panel', 'Greys', 'shared with panel variant')]
    fig = go.Figure()
    for group_key, colorscale, label in groups:
        group_muts = [m for m in sorted_muts if mut_classes.get(m) == group_key]
        if not group_muts:
            continue
        z = pivot.loc[group_muts].values.tolist()
        hover_text = []
        for mut in group_muts:
            row_hover = []
            for col in pivot.columns:
                f = pivot.loc[mut, col]
                c = hover.loc[mut, col] if mut in hover.index else 0
                cv = cov.loc[mut, col] if mut in cov.index else 0
                cls = mut_classes.get(mut, '')
                cls_label = {'unique': '🔵 unique',
                             'shared_cowwid': '🟠 shared (known variant)',
                             'shared_panel': '⚫ shared (panel)'}.get(cls, '')
                row_hover.append(f"<b>{mut}</b> {cls_label}<br>freq: {f:.1%}<br>"
                                 f"{int(c):,} / {int(cv):,} reads")
            hover_text.append(row_hover)
        zmax = max((max(max(r) for r in z), 0.01))
        fig.add_trace(go.Heatmap(
            z=z, x=col_labels, y=group_muts, text=hover_text,
            hovertemplate="%{text}<extra></extra>", colorscale=colorscale,
            zmin=0, zmax=zmax, showscale=(group_key == 'unique'),
            colorbar=dict(title="freq", tickformat=".0%", x=1.02)
            if group_key == 'unique' else None, name=label))
    fig.update_layout(
        title=dict(text=title, font=dict(size=14)),
        height=max(200, 32 * len(sorted_muts) + 90),
        margin=dict(t=50, b=40, l=100, r=60), template="plotly_white",
        xaxis=dict(side="bottom", title=None),
        yaxis=dict(autorange="reversed", title=None,
                   tickfont=dict(family="monospace", size=11)),
        legend=dict(orientation="h", y=-0.12, x=0, font=dict(size=11)))
    return fig


def render_scanner_heatmap(
    variant, mutations, client, location, date_range, max_mutations=30,
    panel_variants=None, all_lineage_signatures=None, cowwid_signatures=None,
    lineage_sig=None, member_blocks=None, shared_mutations=None,
) -> None:
    if not mutations:
        st.caption(f"No discriminating mutations found for {variant}.")
        return
    shown = sorted(mutations)[:max_mutations]
    if len(mutations) > max_mutations:
        st.caption(f"Showing {max_mutations} of {len(mutations)} observed mutations.")
    if all_lineage_signatures is not None and cowwid_signatures is not None and panel_variants:
        mut_classes = _classify_mutations(shown, lineage_sig, panel_variants,
                                          all_lineage_signatures, cowwid_signatures)
    else:
        mut_classes = {m: 'unique' for m in shown}
    n_unique = sum(1 for c in mut_classes.values() if c == 'unique')
    n_shared = sum(1 for c in mut_classes.values() if c != 'unique')
    legend_parts = []
    if n_unique:
        legend_parts.append(
            f'<span style="background:#dbeafe;color:#1e40af;padding:2px 8px;'
            f'border-radius:10px;font-size:11px;margin-right:6px">'
            f'🔵 {n_unique} unique</span>')
    if n_shared:
        legend_parts.append(
            f'<span style="background:#ffedd5;color:#9a3412;padding:2px 8px;'
            f'border-radius:10px;font-size:11px;margin-right:6px">'
            f'🟠 {n_shared} shared with known variant</span>')
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
    unique_muts_in_pivot = [m for m in shown
                            if mut_classes.get(m) == 'unique' and m in pivot.index]
    if unique_muts_in_pivot:
        unique_pivot = pivot.loc[unique_muts_in_pivot]
        max_unique_freq = unique_pivot.values.max()
        if len(pivot.columns) >= 2:
            last2 = unique_pivot.iloc[:, -2:].values
            n_rising = sum(1 for row in last2 if row[-1] > row[0])
        else:
            n_rising, max_unique_freq = 0, 0
        if max_unique_freq > 0.05 and n_rising >= 2:
            st.success(f"**Real signal** — {n_rising} unique mutations rising "
                       f"together (max {max_unique_freq:.1%})")
        elif max_unique_freq > 0.005:
            st.warning(f"**Weak signal** — unique mutations present but low "
                       f"({max_unique_freq:.1%} max). Monitor over time.")
        else:
            st.info("**Likely noise** — unique mutations not rising. Signal may "
                    "be from a related known variant.")
    if member_blocks:
        _guide = []
        if shared_mutations:
            _guide.append(f"<b>shared (clade):</b> {', '.join(shared_mutations[:4])}")
        for _b in member_blocks[:6]:
            _dm = ", ".join(_b.get("discriminating", [])[:3])
            if _dm:
                _guide.append(f"<b>{_b['member']}:</b> {_dm}")
        if _guide:
            st.markdown("<div style='font-size:11px;color:#6b7280;margin:4px 0;'>"
                        "Rows below, grouped by member — a member is present when "
                        "its own row(s) light up alongside the shared rows.<br>"
                        + " &nbsp;·&nbsp; ".join(_guide) + "</div>",
                        unsafe_allow_html=True)
    title = f"{variant} — mutation frequencies by week in {location}"
    fig = _make_figure(pivot, hover_df, cov_df, mut_classes, title)
    st.plotly_chart(fig, use_container_width=True)
    st.caption(
        "Frequency = reads with mutation / reads covering that position. "
        "Blue = mutations unique to this lineage (true signal). "
        "Orange = shared with a known variant not in your panel.")