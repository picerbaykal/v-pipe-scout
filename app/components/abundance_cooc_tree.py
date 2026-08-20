"""
components/abundance_cooc_tree.py

Phylogenetic variant tree for the Abundance & Co-occurrence tab.

Shows the global variant panel rooted at B:
- Selected variants — filled blue circles
- Officially tracked (cowwid OT) variants get an "OT" badge
- Spine/structural nodes — small grey circles connecting the hierarchy

No per-city labels — the panel is global. The only distinction is
OT (officially tracked) vs not.
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


def _ancestors(v: str, parent_map: dict) -> list[str]:
    path = []
    while v in parent_map:
        v = parent_map[v]
        path.append(v)
    return path


def _build_spine(selected_set: set, yaml_set: set, parent_map: dict):
    if not selected_set:
        return None

    needed: set = set()
    for v in selected_set:
        needed.add(v)
        for a in _ancestors(v, parent_map):
            needed.add(a)
    for v in yaml_set:
        if any(a in selected_set for a in _ancestors(v, parent_map) + [v]):
            needed.add(v)
    needed.add("B")

    root = "B"

    def is_desc(v, r):
        x = v
        while x in parent_map:
            if x == r:
                return True
            x = parent_map[x]
        return x == r or v == r

    needed = {v for v in needed if is_desc(v, root)}
    needed.add(root)

    children = defaultdict(list)
    for v in needed:
        p = parent_map.get(v)
        if p and p in needed:
            children[p].append(v)

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


def _build_svg(children, root, kind_of, collapse, width=340):
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
    total_h = len(rows) * ROW_H + 60

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
        color = C[kind]
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

        # OT badge only
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

    # legend
    leg_y = total_h - 34
    leg_svg = [
        f'<line x1="0" y1="{total_h - 48}" x2="{width}" y2="{total_h - 48}" '
        f'stroke="#E8E6E0" stroke-width="1"/>'
    ]
    items = [
        ("#185FA5", True, "selected"),
        ("#185FA5", False, "not selected"),
    ]
    lx = 0
    for lc, lf, ltxt in items:
        if lf:
            leg_svg.append(f'<circle cx="{lx+5}" cy="{leg_y}" r="4" fill="{lc}"/>')
        else:
            leg_svg.append(
                f'<circle cx="{lx+5}" cy="{leg_y}" r="4" fill="white" '
                f'stroke="{lc}" stroke-width="1.5"/>'
            )
        leg_svg.append(
            f'<text x="{lx+13}" y="{leg_y}" dy="0.35em" font-size="10" fill="#888" '
            f'font-family="-apple-system,BlinkMacSystemFont,\'Segoe UI\',sans-serif">'
            f'{ltxt}</text>'
        )
        lx += len(ltxt) * 6 + 26
    # OT chip in legend
    leg_svg.append(
        f'<rect x="{lx}" y="{leg_y-7}" width="26" height="14" rx="7" fill="#F1EFE8"/>'
        f'<text x="{lx+13}" y="{leg_y}" dy="0.35em" text-anchor="middle" '
        f'font-size="9" fill="#5F5E5A" '
        f'font-family="-apple-system,BlinkMacSystemFont,\'Segoe UI\',sans-serif">OT</text>'
    )
    leg_svg.append(
        f'<text x="{lx+32}" y="{leg_y}" dy="0.35em" font-size="10" fill="#888" '
        f'font-family="-apple-system,BlinkMacSystemFont,\'Segoe UI\',sans-serif">'
        f'officially tracked</text>'
    )
    leg_svg.append(
        f'<text x="0" y="{total_h - 12}" font-size="10" fill="#B4B2A9" '
        f'font-family="-apple-system,BlinkMacSystemFont,\'Segoe UI\',sans-serif">'
        f'OT = Officially Tracked (cowwid surveillance panel)</text>'
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
):
    """
    Render the global variant tree. Only distinction: OT vs not-OT.
    (scanner_results, cooc_only, scanner_added_for kept for signature
    compatibility but unused.)
    """
    selected_set = set(selected_variants)
    yaml_set = set(yaml_variants)

    raw = pango_loader.get_raw_data()
    parent_map = {v: e.get("parent", "") for v, e in raw.items() if e.get("parent")}

    RECOMBINANT_PARENT = {
        "XDV": "JN.1", "XFG": "JN.1", "XEC": "JN.1",
        "XBB": "BA.2", "XBB.1": "XBB", "XBB.2": "XBB",
    }
    for alias, par in RECOMBINANT_PARENT.items():
        if alias in raw and alias not in parent_map:
            parent_map[alias] = par

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

    result = _build_spine(selected_set, yaml_set, parent_map)
    if not result:
        st.caption("Select at least one variant to build the tree.")
        return

    children, root, needed, kind_of, collapse = result
    _html, _h = _build_svg(children, root, kind_of, collapse)
    components.html(_html, height=_h + 20, scrolling=True)