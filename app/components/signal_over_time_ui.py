"""components/signal_over_time_ui.py

'Signal over time' (2026-10-02): the co-occurrence and mutation heatmaps of one
lineage or novel pattern in one of the run's cities, per week
(components.scanner_heatmap). Its own section, so the co-occurrence results
and 'Investigate a variant' stay short.

Picker order: lineages the deep scan found (most reads first), the panel,
the run's recurring novel patterns, then any pango lineage. A lineage the
scanner reported in the chosen city uses the scanner's groups of positions
(and its evidence days as ● marks); any other lineage, groups of its ★
markers (process.variant_explorer.marker_blocks).
"""
from __future__ import annotations

from datetime import date

import streamlit as st

from process.variant_explorer import marker_blocks

NOVEL_PREFIX = "novel · "


def _scan_index(scanner_results):
    """From the deep scan, per city: {node: blocks}, {node: counted days};
    the found nodes by reads; the recurring novel patterns."""
    blocks, days, reads, novel = {}, {}, {}, {}
    for city, r in (scanner_results or {}).items():
        for k in ("resolved_clade", "one_day"):
            for f in r.get(k) or []:
                n = f.get("node")
                if not n:
                    continue
                blocks.setdefault(city, {})[n] = f.get("member_blocks") or []
                days.setdefault(city, {})[n] = list(f.get("counted_days") or [])
                reads[n] = reads.get(n, 0) + int(f.get("total_reads", 0) or 0)
        for g in (r.get("novel") or {}).get("groups") or []:
            if len(g.get("days") or []) >= 2 and not g.get("likely_error"):
                novel[" ".join(g["mutations"])] = list(g["mutations"])
    found = [n for n, _r in sorted(reads.items(), key=lambda kv: -kv[1])]
    return blocks, days, found, novel


def render_signal_over_time(pango_loader, client, cities, start_date, end_date,
                            panel=None, scanner_results=None, default_city=None,
                            key_prefix="acooc_sot"):
    st.markdown("#### 📈 Signal over time")
    st.caption("Per week, in one city: which base combinations sit on the same reads "
               "(co-occurrence) and how frequent each mutation is. Regions are folded — "
               "open the ones you want.")
    if not (cities and start_date and end_date and client):
        st.info("Available once the run has cities and dates.")
        return

    blocks_by, days_by, found, novel = _scan_index(scanner_results)
    panel = list(panel or [])
    every = sorted(pango_loader.get_raw_data().keys())
    first = found + [p for p in panel if p not in found]
    options = (first + [NOVEL_PREFIX + k for k in novel]
               + [x for x in every if x not in set(first)])

    c1, c2 = st.columns([3, 2])
    with c1:
        pick = st.selectbox("Lineage or pattern", options, index=0 if options else None,
                            placeholder="Search a lineage…", key=f"{key_prefix}_pick",
                            label_visibility="collapsed")
    with c2:
        idx = cities.index(default_city) if default_city in cities else 0
        city = st.selectbox("City", cities, index=idx, key=f"{key_prefix}_city",
                            label_visibility="collapsed")
    if not pick:
        return
    if found:
        st.caption(f"Found by the deep scan: {', '.join(found[:8])}"
                   + (" …" if len(found) > 8 else "") + " — listed first.")

    from components.scanner_heatmap import render_clade_heatmap
    dr = (date.fromisoformat(str(start_date)[:10]), date.fromisoformat(str(end_date)[:10]))
    if pick.startswith(NOVEL_PREFIX):
        lab = pick[len(NOVEL_PREFIX):]
        muts = novel[lab]
        blocks = [{"member": "novel pattern", "discriminating": muts, "reads": 0,
                   "mut_star": {m: True for m in muts}}]
        name, ev = lab, None
    else:
        name = pick
        blocks = (blocks_by.get(city) or {}).get(pick) or marker_blocks(pick, pango_loader, panel)
        ev = (days_by.get(city) or {}).get(pick)
    if not blocks:
        st.caption(f"{name} has no ★ marker — nothing specific to show.")
        return
    render_clade_heatmap(clade_node=name, shared_mutations=[], member_blocks=blocks,
                         client=client, location=city, date_range=dr, evidence_days=ev)