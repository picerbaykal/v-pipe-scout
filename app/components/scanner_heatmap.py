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
_LEFT_W = 300          # px of the left part (bases + total, or the mutation name):
                       #   the same in every table, so the week columns line up


def _min_cov() -> int:
    """Reads at one position for a share to count (the check's min_cov)."""
    try:
        from process.cooc import _check_cfg
        return int(_check_cfg()["min_cov"])
    except Exception:
        return 100


def _min_link() -> int:
    """Reads covering several positions at once for a share to count (the
    check's link test, min_link)."""
    try:
        from process.cooc import _check_cfg
        return int(_check_cfg()["min_link"])
    except Exception:
        return 20


def _week(d) -> str:
    """ISO week start (Monday) of a sampling date, 'YYYY-MM-DD'."""
    x = datetime.strptime(str(d)[:10], "%Y-%m-%d")
    return (x - pd.Timedelta(days=x.weekday())).strftime("%Y-%m-%d")


def _short(d) -> str:
    x = datetime.strptime(str(d)[:10], "%Y-%m-%d")
    return f"{x.day} {x.strftime('%b')}"


def _weeks(dates):
    """{week start: [sampling dates]} in order."""
    out: Dict[str, List[str]] = {}
    for d in sorted(dates):
        out.setdefault(_week(d), []).append(str(d)[:10])
    return out


# colour steps (2026-10-02): tied to the check's shares, so a minority lineage
# at a few % is visible — on a 0–100 % scale 8 % was almost white
_STEPS = [(0.01, "#ffffff", "#9ca3af"), (0.05, "#dbeafe", "#1f2937"), (0.20, "#93c5fd", "#1f2937"),
          (0.50, "#3b82f6", "#ffffff"), (1.01, "#1e3a8a", "#ffffff")]
_GREY = "#eceef1"


def _colour(share):
    for top, bg, fg in _STEPS:
        if share < top:
            return bg, fg
    return _STEPS[-1][1], _STEPS[-1][2]


def _pct(share):
    if share <= 0:
        return "0%"
    return "<1%" if share < 0.01 else f"{share * 100:.0f}%"


_CSS = """<style>
.hm{border-collapse:separate;border-spacing:2px;font-size:11.5px;table-layout:fixed;margin:4px 0 6px}
.hm th{font-weight:600;color:#6b7280;font-size:10.5px;padding:0 2px;white-space:nowrap}
.hm td{height:22px;padding:0;white-space:nowrap}
.hm td.b{font-family:ui-monospace,Menlo,monospace;font-size:13px;text-align:center}
.hm td.n{text-align:right;padding-right:8px;color:#374151}
.hm td.l{font-family:ui-monospace,Menlo,monospace;font-size:12.5px;padding-left:4px}
.hm td.c{text-align:center;border-radius:2px;cursor:default}
.hm .v{color:#1d4ed8;font-weight:700}.hm .r{color:#9ca3af}
.hm .dot{font-size:8px;margin-left:2px}
.hmnote{font-size:11.5px;color:#6b7280;margin:0 0 4px}
</style>"""


def _legend(extra=""):
    sw = "".join(f"<span style='display:inline-block;width:16px;height:10px;background:{bg};"
                 f"border:1px solid #e5e7eb;margin:0 3px -1px 8px'></span>{lab}"
                 for (bg, lab) in (("#ffffff", "&lt; 1 %"), ("#dbeafe", "1–5 %"),
                                   ("#93c5fd", "5–20 %"), ("#3b82f6", "20–50 %"),
                                   ("#1e3a8a", "≥ 50 %"), (_GREY, "too few reads")))
    st.markdown(f"<div class='hmnote'>{sw}{extra}</div>", unsafe_allow_html=True)


def _cell(n, cov, days, min_reads, what, dot=False):
    """One week cell: share n / cov, with the numbers on hover."""
    tip = (f"week of {_short(days[0])} · {len(days)} sample{'s' if len(days) != 1 else ''} "
           f"({', '.join(_short(d) for d in days)})")
    if cov < min_reads:
        tip += f" · {cov:,} covering reads — too few (< {min_reads})"
        return f"<td class='c' style='background:{_GREY}' title=\"{_esc(tip)}\"></td>"
    sh = n / cov
    bg, fg = _colour(sh)
    tip += f" · {what}: {n:,} of {cov:,} covering reads = {sh * 100:.2f} %"
    if dot:
        tip += " · the scanner counted a day of evidence this week"
    return (f"<td class='c' style='background:{bg};color:{fg}' title=\"{_esc(tip)}\">{_pct(sh)}"
            + ("<span class='dot'>●</span>" if dot else "") + "</td>")


def _esc(t):
    import html as _h
    return _h.escape(str(t), quote=True)


def _region_title(muts):
    ps = [_pos(m) for m in muts]
    return f"region {min(ps):,}–{max(ps):,}"


# ─────────────────────────────────────────────────────────────────────────────
# co-occurrence: the actual base combinations on reads, per week
# ─────────────────────────────────────────────────────────────────────────────
_TOP = 10


def _cooc_table(positions, mut_base, star_pos, dates, per_date, covered, finding,
                evidence_weeks):
    weeks = _weeks(dates)
    full = tuple(mut_base[p] for p in positions)
    totals = defaultdict(int)
    for d in dates:
        for combo, n in per_date.get(d, {}).items():
            totals[combo] += n
    ranked = [c for c, _n in sorted(totals.items(), key=lambda kv: -kv[1])]
    rows = [full] + [c for c in ranked if c != full][:_TOP - 1]
    rest = [c for c in ranked if c not in rows]
    wk_cov = {w: sum(covered.get(d, 0) for d in ds) for w, ds in weeks.items()}
    mr = _min_link()
    # left width = _LEFT_W in every table (border-spacing included) -> weeks line up
    pw = max(30, (_LEFT_W - 72 - 2 * len(positions)) // max(1, len(positions)))
    cols = "".join(f"<col style='width:{pw}px'>" for _ in positions) + "<col style='width:72px'>" \
        + "".join("<col style='width:50px'>" for _ in weeks)
    head = ("<tr>" + "".join(f"<th title=\"position {p}: {mut_base[p]} = the variant's base\">{p}"
                             f"{' ★' if p in star_pos else ''}</th>" for p in positions)
            + "<th style='text-align:right;padding-right:8px'>reads</th>"
            + "".join(f"<th>{_short(w)}</th>" for w in weeks) + "</tr>")
    body = []
    for c in rows:
        me = c == full
        bases = "".join(f"<td class='b'><span class='{'v' if b == mut_base[p] else 'r'}'>{b}</span></td>"
                        for p, b in zip(positions, c))
        what = (f"{finding}'s combination {''.join(c)}" if me else f"combination {''.join(c)}")
        cells = "".join(_cell(sum(per_date.get(d, {}).get(c, 0) for d in ds), wk_cov[w], ds, mr,
                              what, dot=me and w in evidence_weeks)
                        for w, ds in weeks.items())
        body.append(f"<tr class='{'me' if me else ''}'>{bases}<td class='n'>{totals.get(c, 0):,}"
                    f"</td>{cells}</tr>")
    note = (f"+ {len(rest)} more combinations ({sum(totals[c] for c in rest):,} reads) · "
            if rest else "")
    return (f"<table class='hm'>{cols}{head}{''.join(body)}</table>"
            f"<div class='hmnote'>{note}first row = {finding}'s combination (all its "
            f"mutations) · blue base = the variant's, grey = reference · grey cell = fewer than "
            f"{mr} reads cover all positions that week</div>")


def _render_haplotype_blocks(clade_node, star_blocks, is_star, n_outside, client,
                             location, date_range, evidence_days=None):
    """Co-occurrence: one table per genome region with a ★ marker, folded."""
    ev_weeks = {_week(d) for d in (evidence_days or [])}
    blocks = []
    for b in star_blocks:
        positions, mut_base, star_pos = _hap_positions(list(b.get("discriminating", [])),
                                                       is_star, n_outside)
        if len(positions) >= 2:
            blocks.append((positions, mut_base, star_pos))
    if not blocks:
        st.caption("No group of ★ positions that one read can cover.")
        return
    blocks.sort(key=lambda x: x[0][0])
    st.markdown(_CSS + "<div style='font-size:13.5px;font-weight:600;margin-top:8px'>"
                "Co-occurrence</div><div class='hmnote'>the base combinations on reads covering "
                "all positions of a region · cell = share of those reads, per week · ● = the "
                "scanner counted a day of evidence that week</div>", unsafe_allow_html=True)
    _legend()
    for positions, mut_base, star_pos in blocks:
        muts = [f"{p}{mut_base[p]}" for p in positions]
        label = (_region_title(muts) + " · "
                 + " + ".join(f"{m}{'★' if _pos(m) in star_pos else ''}" for m in muts))
        with st.expander(label, expanded=False):
            with st.spinner("Reading co-occurrence…"):
                dates, per_date, covered = _fetch_haplotypes(client, location, date_range,
                                                             tuple(positions))
            if not sum(covered.values()):
                st.caption("No read covers all of these positions in the window.")
                continue
            st.markdown(_CSS + _cooc_table(positions, mut_base, star_pos, dates, per_date,
                                           covered, clade_node, ev_weeks),
                        unsafe_allow_html=True)


# ─────────────────────────────────────────────────────────────────────────────
# mutations: one row per mutation, per week, one table per region
# ─────────────────────────────────────────────────────────────────────────────
def _render_per_mutation_blocks(clade_node, shared_mutations, member_blocks,
                                client, location, date_range,
                                max_muts_per_block=25):
    _is_star, _n_out = _star_info(member_blocks)
    groups = []
    for b in member_blocks or []:
        muts = list(dict.fromkeys(b.get("discriminating", [])))
        if muts:
            groups.append(muts)
    if shared_mutations:
        groups.append(list(shared_mutations))
    if not groups:
        st.caption("No mutations to show.")
        return
    groups.sort(key=lambda g: (not any(_is_star(m) for m in g), min(_pos(m) for m in g)))
    st.markdown(_CSS + "<div style='font-size:13.5px;font-weight:600;margin-top:14px'>"
                "Mutations</div><div class='hmnote'>reads with the mutation / reads covering "
                "its position, per week · ★ = marker of this finding"
                "</div>", unsafe_allow_html=True)
    _legend()
    mc = _min_cov()
    for g in groups:
        stars = sorted((m for m in g if _is_star(m)), key=_pos)
        rest = sorted((m for m in g if not _is_star(m)), key=_pos)
        ordered = (stars + rest)[:max_muts_per_block]
        label = (_region_title(g) + (f" · {len(stars)} ★" if stars
                                     else " · no ★ — shared mutations only"))
        with st.expander(label, expanded=False):
            with st.spinner("Reading mutation frequencies…"):
                df = _fetch_frequencies(client, location, date_range, ordered)
            if df is None or df.empty:
                st.caption("No reads for these positions.")
                continue
            df = df.drop_duplicates(subset=["mutation", "dateFrom"], keep="first")
            weeks = _weeks(df["dateFrom"].unique())
            by = {(r.mutation, str(r.dateFrom)[:10]): (int(r["count"]), int(r.coverage))
                  for _, r in df.iterrows()}
            cols = (f"<col style='width:{_LEFT_W}px'>"
                    + "".join("<col style='width:50px'>" for _ in weeks))
            head = ("<tr><th></th>" + "".join(f"<th>{_short(w)}</th>" for w in weeks)
                    + "</tr>")
            body = []
            for m in ordered:
                cells = []
                for w, ds in weeks.items():
                    n = sum(by.get((m, d), (0, 0))[0] for d in ds)
                    c = sum(by.get((m, d), (0, 0))[1] for d in ds)
                    cells.append(_cell(n, c, ds, mc, f"{m}"))
                star = " <b style='color:#1d4ed8'>★</b>" if _is_star(m) else ""
                n_out = _n_out(m)
                tip = (f"{m}: carried by {n_out} lineages outside this family" if n_out is not None
                       else m)
                body.append(f"<tr><td class='l' title=\"{_esc(tip)}\">{m}{star}</td>"
                            f"{''.join(cells)}</tr>")
            st.markdown(_CSS + f"<table class='hm'>{cols}{head}{''.join(body)}</table>"
                        f"<div class='hmnote'>grey cell = fewer than {mc} reads cover the "
                        f"position that week</div>", unsafe_allow_html=True)


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
    evidence_days=None,
) -> None:
    """Signal over time for a lineage or pattern in one city (2026-10-02):
    co-occurrence tables and mutation tables, one per genome region, per
    week, all folded. evidence_days: the scanner's counted days (● marks)."""
    _is_star, _n_out = _star_info(member_blocks)
    _star_blocks = [b for b in member_blocks
                    if any(_is_star(m) for m in b.get("discriminating", []))]
    if not _star_blocks:
        st.caption("No ★ marker: no combination specific to this finding — see the mutations.")
    else:
        _render_haplotype_blocks(clade_node, _star_blocks, _is_star, _n_out, client,
                                 location, date_range, evidence_days)
    _render_per_mutation_blocks(clade_node, shared_mutations, member_blocks, client,
                                location, date_range, max_muts_per_block=max_muts_per_block)


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