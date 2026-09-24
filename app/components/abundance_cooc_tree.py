"""
components/abundance_cooc_tree.py

Phylogenetic variant tree for the Abundance & Co-occurrence tab.

Shows the global variant panel rooted at B. Selected (panel) variants are
coloured by their co-occurrence verdict:

    green  (#16a34a)  confirmed      — co-occurring in the wastewater
    orange (#b45309)  not found in WW — distinctive haplotype looked for, not seen
    black  (#1f2430)  selected       — chosen, verdict still pending / inconclusive

Officially-tracked (cowwid) variants get an "OT" badge; tracked-but-not-selected
variants are hollow blue; spine/structural nodes are small grey circles.

Verdicts come from the Co-occurrence check (session_state["acooc_verdicts"],
{variant: status}); until a scan has run every selected node shows black.
"""
from __future__ import annotations
from collections import defaultdict

import streamlit as st
import streamlit.components.v1 as components

C = {
    "panel_ot": "#185FA5",
    "panel":    "#185FA5",
    "yaml":     "#185FA5",
    "spine":    "#C8C6BE",
}

# verdict → colour for selected panel nodes. Green is deliberately brighter than
# the check-table teal so it is clearly distinct from the black "pending" nodes.
# Orange matches the scanner list's "Not found in wastewater" band (#b45309).
STATUS_C = {
    "confirmed": "#16a34a",   # green  — present in WW (confirmed)
    "not_found": "#b45309",   # orange — not found in WW (looked for, absent)
    "pending":   "#1f2430",   # black  — chosen, verdict pending
}


def _status_color(v: str, variant_status: dict | None) -> str:
    s = (variant_status or {}).get(v)
    if s == "confirmed":
        return STATUS_C["confirmed"]
    if s == "not_found":
        return STATUS_C["not_found"]
    return STATUS_C["pending"]


def _ancestors(v: str, parent_map: dict) -> list[str]:
    path = []
    while v in parent_map:
        v = parent_map[v]
        path.append(v)
    return path


def _build_spine(selected_set: set, yaml_set: set, parent_map: dict, recomb_set: set = None):
    recomb_set = recomb_set or set()
    # recombinants are shown in a separate section, not the main bifurcating tree
    selected_set = {v for v in selected_set if v not in recomb_set}
    if not selected_set:
        return None

    needed: set = set()
    for v in selected_set:
        needed.add(v)
        for a in _ancestors(v, parent_map):
            needed.add(a)
    for v in yaml_set:
        if v in recomb_set:
            continue
        if any(a in selected_set for a in _ancestors(v, parent_map) + [v]):
            needed.add(v)
    needed.add("B")

    root = "B"

    def chain_root(v):
        """the top of v's parent chain (where parent is empty/absent)."""
        x = v
        seen = set()
        while x in parent_map and x not in seen:
            seen.add(x)
            x = parent_map[x]
        return x

    def keeps(v):
        r = chain_root(v)
        return r == root or (r.startswith("X"))

    needed = {v for v in needed if keeps(v)}
    needed.add(root)

    _recomb_roots = {chain_root(v) for v in needed
                     if chain_root(v) != root and chain_root(v).startswith("X")}
    needed |= _recomb_roots

    children = defaultdict(list)
    for v in needed:
        p = parent_map.get(v)
        if p and p in needed:
            children[p].append(v)
        elif v in _recomb_roots:
            children[root].append(v)  # recombinant root hangs under B

    def kind_of(v):
        if v in selected_set:
            return "panel_ot" if v in yaml_set else "panel"
        if v in yaml_set:
            return "yaml"
        return "spine"

    def collapse(v):
        parts = [v]
        cur = v
        while True:
            ch = children.get(cur, [])
            if len(ch) != 1:
                break
            child = ch[0]
            if kind_of(child) != "spine":
                break
            parts.append(child)
            cur = child
        return " › ".join(parts), cur

    return children, root, needed, kind_of, collapse


def _build_svg(children, root, kind_of, collapse, width=340,
               recombinant_set=None, variant_status=None):
    recombinant_set = recombinant_set or set()
    ROW_H, INDENT, X0 = 26, 20, 16
    rows = []
    visited = set()

    def assign_rows(v, depth):
        if v in visited:
            return
        visited.add(v)
        kind = kind_of(v)
        if kind == "spine" and depth > 0:
            label, real_v = collapse(v)
        else:
            label, real_v = v.replace("_", " "), v
        rows.append((v, label, real_v, depth, kind))
        ch = sorted(children.get(real_v, []), key=lambda x: (
            0 if kind_of(x) in ("panel", "panel_ot") else
            1 if kind_of(x) == "yaml" else 2, x))
        for c in ch:
            assign_rows(c, depth + 1)

    assign_rows(root, 0)
    row_y = {rows[i][0]: i * ROW_H + ROW_H // 2 for i in range(len(rows))}
    total_h = len(rows) * ROW_H + 74   # extra room for the two-row legend

    def node_color(v, kind):
        if kind in ("panel", "panel_ot"):
            return _status_color(v, variant_status)
        return C[kind]

    lines_svg, nodes_svg = [], []

    for v, _, real_v, depth, kind in rows:
        ch = sorted(children.get(real_v, []), key=lambda x: (
            0 if kind_of(x) in ("panel", "panel_ot") else
            1 if kind_of(x) == "yaml" else 2, x))
        if not ch:
            continue
        x = X0 + depth * INDENT
        r_parent = 5 if kind not in ("spine",) else 3
        y_start = row_y[v] + r_parent + 1
        y_last = row_y[ch[-1]]
        lines_svg.append(
            f'<line x1="{x}" y1="{y_start}" x2="{x}" y2="{y_last}" '
            f'stroke="#D3D1C7" stroke-width="1.5"/>'
        )

    for v, label, real_v, depth, kind in rows:
        x = X0 + depth * INDENT
        y = row_y[v]
        color = node_color(v, kind)
        is_spine = kind == "spine"
        filled = kind in ("panel", "panel_ot")
        fw = "600" if filled else "400"
        fsize = "11" if is_spine else "13"
        fcolor = color if not is_spine else "#B4B2A9"
        r = 5 if not is_spine else 3

        if depth > 0:
            px = X0 + (depth - 1) * INDENT
            branch_color = color if not is_spine else "#D3D1C7"
            lines_svg.append(
                f'<line x1="{px+1}" y1="{y}" x2="{x-r-2}" y2="{y}" '
                f'stroke="{branch_color}" stroke-width="1.8"/>'
            )

        if filled:
            nodes_svg.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{color}"/>')
        elif is_spine:
            nodes_svg.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="#D3D1C7"/>')
        else:
            nodes_svg.append(
                f'<circle cx="{x}" cy="{y}" r="{r}" fill="white" '
                f'stroke="{color}" stroke-width="1.8"/>'
            )

        tx = x + r + 6
        nodes_svg.append(
            f'<text x="{tx}" y="{y}" dy="0.35em" font-size="{fsize}" '
            f'font-weight="{fw}" fill="{fcolor}" '
            f'font-family="-apple-system,BlinkMacSystemFont,\'Segoe UI\',sans-serif">'
            f'{label}</text>'
        )

        # OT badge
        _badge_x_end = tx + len(label) * (7 if not is_spine else 6) + 16
        if kind in ("panel_ot", "yaml"):
            lw = len(label) * (7 if not is_spine else 6) + 16
            bx = tx + lw
            bw = 6.2 * 2 + 14
            bh = 15
            nodes_svg.append(
                f'<rect x="{bx}" y="{y - bh//2}" width="{bw}" height="{bh}" '
                f'rx="7" fill="#F1EFE8"/>'
                f'<text x="{bx + bw/2:.1f}" y="{y}" dy="0.35em" '
                f'text-anchor="middle" font-size="10" fill="#5F5E5A" '
                f'font-family="-apple-system,BlinkMacSystemFont,\'Segoe UI\',sans-serif">'
                f'OT</text>'
            )
            _badge_x_end = bx + bw + 4

        # recombinant badge (mixed ancestry — no single parent in the tree)
        if real_v in recombinant_set:
            bw2 = 6.2 * 6 + 14
            bh = 15
            bx2 = _badge_x_end
            nodes_svg.append(
                f'<rect x="{bx2}" y="{y - bh//2}" width="{bw2}" height="{bh}" '
                f'rx="7" fill="#F3E8F1"/>'
                f'<text x="{bx2 + bw2/2:.1f}" y="{y}" dy="0.35em" '
                f'text-anchor="middle" font-size="10" fill="#8A4E82" '
                f'font-family="-apple-system,BlinkMacSystemFont,\'Segoe UI\',sans-serif">'
                f'recomb</text>'
            )

    # ── legend (two rows: status swatches, then tracked/OT) ──────────────────
    leg_svg = [
        f'<line x1="0" y1="{total_h - 60}" x2="{width}" y2="{total_h - 60}" '
        f'stroke="#E8E6E0" stroke-width="1"/>'
    ]

    def _dot(cx, cy, fill, hollow=False):
        if hollow:
            return (f'<circle cx="{cx}" cy="{cy}" r="4" fill="white" '
                    f'stroke="{fill}" stroke-width="1.5"/>')
        return f'<circle cx="{cx}" cy="{cy}" r="4" fill="{fill}"/>'

    def _txt(x, y, s, fill="#888"):
        return (f'<text x="{x}" y="{y}" dy="0.35em" font-size="10" fill="{fill}" '
                f'font-family="-apple-system,BlinkMacSystemFont,\'Segoe UI\',sans-serif">'
                f'{s}</text>')

    # row 1: verdict colours — wording matches the scanner list / check labels
    row1_y = total_h - 44
    row1 = [
        (STATUS_C["confirmed"], False, "confirmed"),
        (STATUS_C["not_found"], False, "not found in WW"),
        (STATUS_C["pending"], False, "selected"),
    ]
    lx = 0
    for lc, hollow, ltxt in row1:
        leg_svg.append(_dot(lx + 5, row1_y, lc, hollow))
        leg_svg.append(_txt(lx + 13, row1_y, ltxt))
        lx += len(ltxt) * 6 + 26

    # row 2: tracked (hollow) + OT chip
    row2_y = total_h - 26
    lx = 0
    leg_svg.append(_dot(lx + 5, row2_y, "#185FA5", hollow=True))
    leg_svg.append(_txt(lx + 13, row2_y, "tracked, not selected"))
    lx += len("tracked, not selected") * 6 + 26
    leg_svg.append(
        f'<rect x="{lx}" y="{row2_y-7}" width="26" height="14" rx="7" fill="#F1EFE8"/>'
        f'<text x="{lx+13}" y="{row2_y}" dy="0.35em" text-anchor="middle" '
        f'font-size="9" fill="#5F5E5A" '
        f'font-family="-apple-system,BlinkMacSystemFont,\'Segoe UI\',sans-serif">OT</text>'
    )
    leg_svg.append(_txt(lx + 32, row2_y, "officially tracked"))

    leg_svg.append(
        f'<text x="0" y="{total_h - 8}" font-size="9" fill="#B4B2A9" '
        f'font-family="-apple-system,BlinkMacSystemFont,\'Segoe UI\',sans-serif">'
        f'colour = co-occurrence verdict · black until a scan has run</text>'
    )

    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{total_h}" style="display:block">'
        + "".join(lines_svg) + "".join(nodes_svg) + "".join(leg_svg) + "</svg>"
    )
    html = (
        f'<!DOCTYPE html><html><head>'
        f'<style>body{{margin:0;padding:8px 4px;background:transparent}}</style>'
        f'</head><body>{svg}</body></html>'
    )
    return html, total_h


def render_panel_tree(
    selected_variants: list[str],
    yaml_variants: list[str],
    pango_loader,
    scanner_results: dict | None = None,
    cooc_only: set | None = None,
    scanner_added_for: dict | None = None,
    variant_status: dict | None = None,
):
    """
    Render the global variant tree, colouring selected nodes by their
    co-occurrence verdict.

    variant_status: {variant: status} from the Co-occurrence check
        ("confirmed" -> green, "not_found" -> orange, anything else / missing
        -> black "pending"). Falls back to session_state["acooc_verdicts"].
    (scanner_results, cooc_only, scanner_added_for kept for signature
    compatibility but unused.)
    """
    if variant_status is None:
        variant_status = st.session_state.get("acooc_verdicts", {}) or {}

    # Reset the colouring whenever the panel differs from the last completed run
    # (or nothing has run): every selected node shows black (pending) until a new
    # run recomputes verdicts. Stops stale colours lingering after the panel,
    # dates, or locations change.
    _ran_panel = st.session_state.get("acooc_ran_panel")
    _stale = (_ran_panel is None
              or sorted(set(selected_variants)) != sorted(set(_ran_panel)))
    if _stale:
        variant_status = {}

    selected_set = set(selected_variants)
    yaml_set = set(yaml_variants)

    raw = pango_loader.get_raw_data()
    parent_map = {v: e.get("parent", "") for v, e in raw.items() if e.get("parent")}

    all_known = set(raw.keys())
    for v in list(selected_set) + list(yaml_set):
        if v in parent_map or "." not in v:
            continue
        name_par = v.rsplit(".", 1)[0]
        if name_par in all_known:
            parent_map[v] = name_par
        elif "." in name_par:
            parent_map[v] = name_par
            cur = name_par
            while cur not in parent_map and "." in cur:
                par = cur.rsplit(".", 1)[0]
                if par in all_known or par in parent_map:
                    parent_map[cur] = par
                    break
                parent_map[cur] = par
                cur = par

    if not selected_set and not yaml_set:
        st.caption("Select at least one variant to build the tree.")
        return

    def _recomb_root(v):
        head = v.split(".")[0]
        if head.startswith("X") and head in raw and not raw.get(head, {}).get("parent"):
            return head
        return None

    _all_panel = set(selected_variants) | set(yaml_variants)
    recomb_members = {v for v in _all_panel if _recomb_root(v) is not None}

    result = _build_spine(selected_set, yaml_set, parent_map, recomb_set=recomb_members)

    # ── main bifurcating tree (non-recombinants) ──
    if result:
        children, root, needed, kind_of, collapse = result
        _html, _h = _build_svg(children, root, kind_of, collapse,
                               variant_status=variant_status)
        components.html(_html, height=_h + 20, scrolling=True)
    elif not recomb_members:
        st.caption("Select at least one variant to build the tree.")
        return

    # ── separate Recombinants section (grouped by family), coloured by verdict ──
    selected_recomb = {v for v in selected_variants if v in recomb_members}
    if selected_recomb:
        from collections import defaultdict as _dd
        fam = _dd(list)
        for v in selected_recomb:
            fam[_recomb_root(v)].append(v)
        st.caption("Recombinants (mixed ancestry — shown separately from the tree):")
        for froot in sorted(fam):
            members = sorted(fam[froot])
            lines = []
            for m in members:
                is_ot = m in set(yaml_variants)
                badge = " · OT" if is_ot else ""
                indent = "&nbsp;&nbsp;&nbsp;" if m != froot else ""
                _mc = _status_color(m, variant_status)
                lines.append(
                    f"{indent}<span style='color:{_mc};'>●</span> "
                    f"<b style='color:{_mc};'>{m}</b>{badge}")
            st.markdown(
                f"<div style='border:0.5px solid #e5e7eb;border-radius:8px;"
                f"padding:6px 10px;margin:4px 0;font-size:13px;'>"
                f"<span style='color:#8A4E82;font-weight:600;'>⧉ {froot} recombinant</span><br>"
                + "<br>".join(lines) + "</div>",
                unsafe_allow_html=True,
            )