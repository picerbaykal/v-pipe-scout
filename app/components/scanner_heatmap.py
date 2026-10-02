"""Scanner heatmap component.

Two views of a scanner finding's signal over time:

1. CO-OCCURRENCE (default). For each group of nearby positions that carries a
   ★ marker (block["mut_star"], decided by the scanner), the reads covering
   all of them, per sample, in three fixed rows: all mutations (the finding) ·
   ★ but not all (a sublineage or a dropout) · no ★ (relatives / reference).

2. MUTATIONS (checkbox): one row per mutation, ★ first — its frequency per
   sample.

Both use one grid (2026-10-02): same margins, date columns and row height, so
the heatmaps of a finding line up; grey = fewer reads than the check's min_cov
(a share from a handful of reads jumps between 0 and 100 %); the right column
is the whole window.

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

_STAR_MAX = 30          # legacy fallback only: results without block["mut_star"]
_MAX_POS = 5            # cap positions per haplotype block (keeps coverage high)
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


def _star_info(member_blocks):
    """(is_star(m), outside(m)) from the scanner's blocks.

    ★ = rare outside the finding's own family (block["mut_star"], set by the
    scanner with the same rule it confirms findings with). outside(m) = number of
    lineages outside the family carrying m. Older results without these fields
    fall back to the legacy global carrier count <= _STAR_MAX."""
    star, outside, car = {}, {}, {}
    for b in member_blocks or []:
        star.update(b.get("mut_star", {}) or {})
        outside.update(b.get("mut_outside", {}) or {})
        car.update(b.get("mut_carriers", {}) or {})

    def is_star(m):
        if m in star:
            return bool(star[m])
        n = car.get(m)
        return n is not None and n <= _STAR_MAX

    def n_outside(m):
        return outside.get(m, car.get(m))

    return is_star, n_outside


def _hap_positions(blk_muts: List[str], is_star, n_outside):
    """Choose <= _MAX_POS positions for a block: all ★ positions, then fill with
    the shared mutations carried by the fewest lineages outside the family.
    Returns (positions, mut_base{pos->alt}, star_pos set)."""
    pos_alt, star_pos = {}, set()
    for m in blk_muts:
        p = _pos(m)
        pos_alt[p] = m[-1]
        if is_star(m):
            star_pos.add(p)
    all_pos = sorted(pos_alt)
    if len(all_pos) <= _MAX_POS:
        keep = all_pos
    else:
        stars = sorted(star_pos)[:_MAX_POS]
        back = sorted((p for p in all_pos if p not in star_pos),
                      key=lambda p: n_outside(f"{p}{pos_alt[p]}") or 10**9)
        keep = sorted(stars + back[:max(0, _MAX_POS - len(stars))])
    return keep, {p: pos_alt[p] for p in keep}, {p for p in keep if p in star_pos}


# ─────────────────────────────────────────────────────────────────────────────
# one grid for both views (2026-10-02): same left margin, same date columns,
# same row height, same grey rule — so every heatmap of a finding lines up
# ─────────────────────────────────────────────────────────────────────────────
_LEFT, _RIGHT, _ROW_H = 170, 64, 24


def _min_cov() -> int:
    try:
        from process.cooc import _check_cfg
        return int(_check_cfg()["min_cov"])
    except Exception:
        return 100


def _date_label(d):
    x = datetime.strptime(str(d)[:10], "%Y-%m-%d")
    return f"{x.day} {x.strftime('%b')}"


def _grid(row_labels, right_notes, dates, cells, key, sep_after=None):
    """rows × sampling dates. cells[i][j] = (count, coverage). A cell is the
    share count / coverage; grey when coverage < min_cov (a share from a handful
    of reads jumps between 0 and 100 %)."""
    mc = _min_cov()
    x = [_date_label(d) for d in dates]
    z, grey, txt = [], [], []
    for lab, row in zip(row_labels, cells):
        zr, gr, tr = [], [], []
        for (n, cov), xl in zip(row, x):
            if cov < mc:
                zr.append(np.nan); gr.append(1)
                tr.append(f"{xl} · {lab}<br>too few reads ({cov:,} < {mc})")
            else:
                zr.append(n / cov); gr.append(np.nan)
                tr.append(f"{xl} · {lab}<br>{n / cov * 100:.1f} % ({n:,} / {cov:,} reads)")
        z.append(zr); grey.append(gr); txt.append(tr)
    fig = go.Figure()
    fig.add_trace(go.Heatmap(z=z, x=x, y=row_labels, text=txt, hovertemplate="%{text}<extra></extra>",
                             colorscale="Blues", zmin=0, zmax=1, showscale=False,
                             hoverongaps=False, xgap=2, ygap=2))
    fig.add_trace(go.Heatmap(z=grey, x=x, y=row_labels, colorscale=[[0, "#eceef1"], [1, "#eceef1"]],
                             showscale=False, hoverinfo="skip", xgap=2, ygap=2))
    for i, note in enumerate(right_notes):
        fig.add_annotation(x=1.0, xref="paper", xanchor="left", xshift=6, y=row_labels[i],
                           text=note, showarrow=False, font=dict(size=11, color="#6b7280"))
    if sep_after:
        fig.add_shape(type="line", xref="paper", x0=0, x1=1, y0=sep_after - 0.5,
                      y1=sep_after - 0.5, line=dict(color="#9ca3af", width=1, dash="dot"))
    fig.update_layout(
        height=_ROW_H * len(row_labels) + 46, template="plotly_white",
        margin=dict(l=_LEFT, r=_RIGHT, t=4, b=34),
        yaxis=dict(autorange="reversed", tickfont=dict(size=12, family="monospace"),
                   automargin=False),
        xaxis=dict(tickfont=dict(size=10), tickangle=0, automargin=False, side="bottom"))
    st.plotly_chart(fig, use_container_width=True, key=key, config={"displayModeBar": False})


def _title(text, sub=""):
    st.markdown(f"<div style='font-size:13px;font-weight:600;margin:12px 0 0;'>{text}"
                + (f" <span style='font-weight:400;color:#6b7280'>· {sub}</span>" if sub else "")
                + "</div>", unsafe_allow_html=True)


def _one_haplotype_block(positions, mut_base, star_pos, dates, per_date, covered,
                         title, subtitle, key):
    """One ★ block as three fixed rows, so every block lines up:
      all mutations   the block's full combination (the finding's haplotype)
      ★, not all      its ★ marker without the rest (a sublineage / a dropout)
      no ★            relatives that carry only shared mutations, or reference
    Cell = share of the reads covering every position of the block."""
    combo_all = "·".join(f"{p}{mut_base[p]}" for p in positions)
    rows = [("all mutations", lambda c: all(c[i] == mut_base[p] for i, p in enumerate(positions))),
            ("★, not all", None),
            ("no ★", lambda c: not any(c[i] == mut_base[p]
                                       for i, p in enumerate(positions) if p in star_pos))]
    counts = {r: {} for r, _ in rows}
    for d in dates:
        for combo, n in per_date.get(d, {}).items():
            if rows[0][1](combo):
                r = rows[0][0]
            elif rows[2][1](combo):
                r = rows[2][0]
            else:
                r = rows[1][0]
            counts[r][d] = counts[r].get(d, 0) + n
    cov_tot = sum(covered.get(d, 0) for d in dates)
    labels = [r for r, _ in rows]
    notes = [(f"{sum(counts[r].values()) / cov_tot * 100:.0f} %" if cov_tot else "–")
             for r in labels]
    cells = [[(counts[r].get(d, 0), covered.get(d, 0)) for d in dates] for r in labels]
    pos_s = " + ".join(f"{p}{mut_base[p]}{'★' if p in star_pos else ''}" for p in positions)
    _title(pos_s, subtitle)
    if not cov_tot:
        st.caption("No read covers all of these positions in the window.")
        return
    _grid(labels, notes, dates, cells, key)


def _render_haplotype_blocks(clade_node, star_blocks, is_star, n_outside, client,
                             location, date_range):
    """Co-occurrence: one grid per ★ block, the block with most reads open."""
    blocks = sorted(star_blocks, key=lambda b: -b.get("reads", 0))
    st.caption(f"Co-occurrence — share of the reads covering all of a group's positions, "
               f"per sample · grey = under {_min_cov()} reads · right: whole window")
    for i, b in enumerate(blocks):
        muts = list(b.get("discriminating", []))
        positions, mut_base, star_pos = _hap_positions(muts, is_star, n_outside)
        if len(positions) < 2:
            continue
        reads = b.get("reads", 0)
        subtitle = f"{reads:,} reads in the scan" if reads else ""
        with st.spinner("Reading co-occurrence…"):
            dates, per_date, covered = _fetch_haplotypes(client, location, date_range, positions)
        key = f"hap_{clade_node}_{location}_{i}"
        if i == 0:
            _one_haplotype_block(positions, mut_base, star_pos, dates, per_date, covered,
                                 "", subtitle, key)
        else:
            with st.expander(" + ".join(str(p) for p in positions)
                             + (f" · {subtitle}" if subtitle else ""), expanded=False):
                _one_haplotype_block(positions, mut_base, star_pos, dates, per_date,
                                     covered, "", subtitle, key)


# ─────────────────────────────────────────────────────────────────────────────
# per-mutation view (opt-in): one row per mutation, ★ first
# ─────────────────────────────────────────────────────────────────────────────
def _render_per_mutation_blocks(clade_node, shared_mutations, member_blocks,
                                client, location, date_range,
                                max_muts_per_block=25):
    """One grid per genome region (the same groups as the co-occurrence view),
    ★ markers first in each; groups without a ★ folded below."""
    _is_star, _n_out = _star_info(member_blocks)

    def _order(muts):
        stars = sorted((m for m in muts if _is_star(m)), key=_pos)
        rest = sorted((m for m in muts if not _is_star(m)), key=_pos)
        return (stars + rest)[:max_muts_per_block], len(stars)

    groups = []
    if shared_mutations:
        groups.append(("shared by the whole clade", list(shared_mutations)))
    for b in member_blocks or []:
        muts = list(dict.fromkeys(b.get("discriminating", [])))
        if muts:
            groups.append((None, muts))
    if not groups:
        st.caption("No mutations to show.")
        return
    all_muts = list(dict.fromkeys(m for _t, g in groups for m in g))
    with st.spinner("Reading mutation frequencies…"):
        df = _fetch_frequencies(client, location, date_range, all_muts)
    if df is None or df.empty:
        st.caption("No frequency data returned.")
        return
    df = df.drop_duplicates(subset=["mutation", "dateFrom"], keep="first")
    dates = sorted(df["dateFrom"].unique())
    by = {(r.mutation, r.dateFrom): (int(r["count"]), int(r.coverage)) for _, r in df.iterrows()}

    def _one(title, muts, key):
        ordered, n_star = _order(muts)
        labels, notes, cells = [], [], []
        for m in ordered:
            row = [by.get((m, d), (0, 0)) for d in dates]
            n, c = sum(x[0] for x in row), sum(x[1] for x in row)
            labels.append(f"{m} ★" if _is_star(m) else m)
            notes.append(f"{n / c * 100:.0f} %" if c else "–")
            cells.append(row)
        ps = [_pos(m) for m in muts]
        _title(title or f"region {min(ps):,}–{max(ps):,}",
               f"{n_star} ★" if n_star else "no ★ — shared mutations only")
        _grid(labels, notes, dates, cells, key=key,
              sep_after=n_star if 0 < n_star < len(ordered) else None)

    st.caption(f"Mutations — reads with it / reads covering its position, per sample · "
               f"one heatmap per genome region · ★ = marker of this finding · "
               f"grey = under {_min_cov()} reads")
    star_g = [(t, g) for t, g in groups if t or any(_is_star(m) for m in g)]
    back_g = [(t, g) for t, g in groups if not t and not any(_is_star(m) for m in g)]
    star_g.sort(key=lambda tg: (tg[0] is None, min(_pos(m) for m in tg[1])))
    for i, (t, g) in enumerate(star_g):
        _one(t, g, f"clade_pm_{clade_node}_{location}_{i}")
    if back_g:
        with st.expander(f"Regions without a ★ marker ({len(back_g)})", expanded=False):
            for i, (t, g) in enumerate(sorted(back_g, key=lambda tg: min(_pos(m) for m in tg[1]))):
                _one(t, g, f"clade_pmb_{clade_node}_{location}_{i}")


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
    _is_star, _n_out = _star_info(member_blocks)

    def _blk_has_star(b):
        return any(_is_star(m) for m in b.get("discriminating", []))

    _star_blocks = [b for b in member_blocks if _blk_has_star(b)]

    if not _star_blocks:
        st.caption("No ★ marker: no combination specific to this finding — see the "
                   "mutations below.")
    else:
        _render_haplotype_blocks(clade_node, _star_blocks, _is_star, _n_out, client,
                                 location, date_range)

    # ── DETAIL: old per-mutation heatmap, opt-in ────────────────────────────
    _detail_key = f"pm_detail_{clade_node}_{location}"
    if st.checkbox("Show mutations one by one", key=_detail_key, value=not _star_blocks):
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