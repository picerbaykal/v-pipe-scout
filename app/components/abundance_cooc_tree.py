"""
components/abundance_cooc_tree.py

Interactive phylogenetic variant tree for the Abundance & Co-occurrence tab.

Renders a pruned pango tree rooted at B, showing:
- Selected variants (in deconv panel) — filled blue circles
- Cowwid OT variants (officially tracked, not selected) — open blue circles
  with "OT" badge
- Spine/structural nodes — small grey circles connecting the hierarchy
- Scanner signal placeholders (v1 — requires co-occurrence/SaneQL data)

Adapted from lollijac's core/panel_tree_component.py.
Scanner/cooc logic stripped for v0 — parameters kept as None defaults
so the function signature is forward-compatible with v1.
"""
from __future__ import annotations
from collections import defaultdict

import streamlit as st
import streamlit.components.v1 as components

# ── Visual config ─────────────────────────────────────────────────────────────
C = {
    "panel_ot": "#185FA5",
    "panel":    "#185FA5",
    "cooc":     "#1D9E75",
    "signal":   "#D85A30",
    "moderate": "#BA7517",
    "clade":    "#BA7517",
    "parallel": "#7F77DD",
    "yaml":     "#185FA5",
    "none":     "#C8C6BE",
    "spine":    "#C8C6BE",
}

BADGE = {
    "panel_ot": ("OT", "#F1EFE8", "#5F5E5A"),
    "panel":    ("", "", ""),
    "cooc":     ("cooc only",       "#E1F5EE", "#085041"),
    "signal":   ("strong signal",   "#FAECE7", "#712B13"),
    "moderate": ("moderate signal", "#FAEEDA", "#633806"),
    "clade":    ("unresolved",      "#FAEEDA", "#633806"),
    "parallel": ("parallel signal", "#EEEDFE", "#3C3489"),
    "yaml": ("OT", "#F1EFE8", "#5F5E5A"),
    "none":     ("", "", ""),
    "spine":    ("", "", ""),
}

# ── Ancestry helpers ──────────────────────────────────────────────────────────

def _ancestors(v: str, parent_map: dict) -> list[str]:
    path = []
    while v in parent_map:
        v = parent_map[v]
        path.append(v)
    return path


def _lca(variants: list[str], parent_map: dict) -> str:
    """Always root at B so the full spine is always visible."""
    return "B"


def _build_spine(
    selected_set: set,
    cooc_set: set,
    yaml_set: set,
    scanner_results: dict | None,
    parent_map: dict,
):
    if not selected_set:
        return None

    # Scanner nodes — v1, not yet available
    # Structure kept for forward compatibility when SaneQL lands
    scanner_nodes: dict[str, str] = {}
    if scanner_results:
        for v, res in scanner_results.items():
            if v in selected_set or v in yaml_set:
                continue
            sl = getattr(res, 'signal_level', None)
            sl = sl.value if sl else ""
            has_cooc = bool(
                getattr(res, 'cooc_signal', None) and
                any(f > 0.01 for f in res.cooc_signal)
            )
            if sl == "strong" and has_cooc:
                scanner_nodes[v] = "signal"
            elif sl == "moderate" and has_cooc:
                scanner_nodes[v] = "moderate"

    needed: set = set()
    for v in selected_set:
        needed.add(v)
        for a in _ancestors(v, parent_map):
            needed.add(a)
    for v in scanner_nodes:
        needed.add(v)
        for a in _ancestors(v, parent_map):
            needed.add(a)
    for v in yaml_set:
        if any(a in selected_set for a in _ancestors(v, parent_map) + [v]):
            needed.add(v)

    needed.add("B")

    children: dict = defaultdict(list)
    for v in needed:
        p = parent_map.get(v)
        if p and p in needed:
            children[p].append(v)

    root = _lca(list(selected_set), parent_map)

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
            if v in yaml_set:
                return "panel_ot"
            return "panel"
        if v in cooc_set:
            return "cooc"
        if v in yaml_set:
            return "yaml"
        if v in scanner_nodes:
            return scanner_nodes[v]
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

# ── SVG renderer ──────────────────────────────────────────────────────────────

def _build_svg(children, root, kind_of, collapse, width=320):
    ROW_H  = 26
    INDENT = 20
    X0     = 16

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
            label = v.replace("_", " ")
            real_v = v
        rows.append((v, label, real_v, depth, kind))
        ch = sorted(children.get(real_v, []), key=lambda x: (
            0 if kind_of(x) == "panel" else
            1 if kind_of(x) in ("signal", "clade", "parallel") else
            2 if kind_of(x) == "yaml" else 3, x))
        for c in ch:
            assign_rows(c, depth + 1)

    assign_rows(root, 0)

    row_y   = {rows[i][0]: i * ROW_H + ROW_H // 2 for i in range(len(rows))}
    total_h = len(rows) * ROW_H + 80

    lines_svg = []
    nodes_svg = []

    for v, _, real_v, depth, kind in rows:
        ch = sorted(children.get(real_v, []), key=lambda x: (
            0 if kind_of(x) == "panel" else
            1 if kind_of(x) in ("signal", "clade", "parallel") else
            2 if kind_of(x) == "yaml" else 3, x))
        if not ch:
            continue
        x = X0 + depth * INDENT
        r_parent = 5 if kind not in ("spine", "none") else 3
        y_start  = row_y[v] + r_parent + 1
        y_last   = row_y[ch[-1]]
        lines_svg.append(
            f'<line x1="{x}" y1="{y_start}" x2="{x}" y2="{y_last}" '
            f'stroke="#D3D1C7" stroke-width="1.5"/>'
        )

    for v, label, real_v, depth, kind in rows:
        x = X0 + depth * INDENT
        y = row_y[v]
        color = C[kind]
        is_spine = kind in ("spine", "none")
        filled = kind in ("panel", "panel_ot", "cooc", "signal")
        fw = "600" if kind in ("panel", "panel_ot", "cooc") else "400"
        fsize = "11" if is_spine else "13"
        fcolor = color if not is_spine else "#B4B2A9"
        r = 5 if not is_spine else 3

        if depth > 0:
            px = X0 + (depth - 1) * INDENT
            branch_color = color if not is_spine else "#D3D1C7"
            dash = ' stroke-dasharray="4,3"' if kind in ("parallel", "none") else ""
            lines_svg.append(
                f'<line x1="{px+1}" y1="{y}" x2="{x-r-2}" y2="{y}" '
                f'stroke="{branch_color}" stroke-width="1.8"{dash}/>'
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
            f'<text x="{tx}" y="{y}" dy="0.35em" '
            f'font-size="{fsize}" font-weight="{fw}" fill="{fcolor}" '
            f'font-family="-apple-system,BlinkMacSystemFont,\'Segoe UI\',sans-serif">'
            f'{label}</text>'
        )

        bdg, bbg, bfg = BADGE[kind]
        if bdg:
            lw  = len(label) * (7 if not is_spine else 6) + 16
            bx  = tx + lw
            bw  = len(bdg) * 6.2 + 14
            bh  = 15
            nodes_svg.append(
                f'<rect x="{bx}" y="{y - bh//2}" width="{bw}" height="{bh}" '
                f'rx="7" fill="{bbg}"/>'
                f'<text x="{bx + bw/2:.1f}" y="{y}" dy="0.35em" '
                f'text-anchor="middle" font-size="10" fill="{bfg}" '
                f'font-family="-apple-system,BlinkMacSystemFont,\'Segoe UI\',sans-serif">'
                f'{bdg}</text>'
            )

    legend_items_row1 = [
        ("#185FA5", True, "selected"),
        ("#185FA5", False, "not selected"),
        ("#D85A30", True, "signal"),
    ]
    legend_items_row2 = [
        ("#1D9E75", True, "cooc only"),
        ("#BA7517", True, "unresolved"),
    ]

    leg_y1 = total_h - 42
    leg_y2 = total_h - 28
    leg_svg = [
        f'<line x1="0" y1="{total_h - 56}" x2="{width}" y2="{total_h - 56}" '
        f'stroke="#E8E6E0" stroke-width="1"/>'
    ]

    leg_x = 0
    for lc, lf, ltxt in legend_items_row1:
        if lf:
            leg_svg.append(
                f'<circle cx="{leg_x + 5}" cy="{leg_y1}" r="4" '
                f'fill="{lc}" stroke="{lc}" stroke-width="1.5"/>'
            )
        else:
            leg_svg.append(
                f'<circle cx="{leg_x + 5}" cy="{leg_y1}" r="4" '
                f'fill="white" stroke="{lc}" stroke-width="1.5"/>'
            )
        leg_svg.append(
            f'<text x="{leg_x + 13}" y="{leg_y1}" dy="0.35em" '
            f'font-size="10" fill="#888" '
            f'font-family="-apple-system,BlinkMacSystemFont,\'Segoe UI\',sans-serif">'
            f'{ltxt}</text>'
        )
        leg_x += len(ltxt) * 6 + 26

    leg_x = 0
    for lc, lf, ltxt in legend_items_row2:
        if lf:
            leg_svg.append(
                f'<circle cx="{leg_x + 5}" cy="{leg_y2}" r="4" '
                f'fill="{lc}" stroke="{lc}" stroke-width="1.5"/>'
            )
        else:
            leg_svg.append(
                f'<circle cx="{leg_x + 5}" cy="{leg_y2}" r="4" '
                f'fill="white" stroke="{lc}" stroke-width="1.5"/>'
            )
        leg_svg.append(
            f'<text x="{leg_x + 13}" y="{leg_y2}" dy="0.35em" '
            f'font-size="10" fill="#888" '
            f'font-family="-apple-system,BlinkMacSystemFont,\'Segoe UI\',sans-serif">'
            f'{ltxt}</text>'
        )
        leg_x += len(ltxt) * 6 + 26

    leg_svg.append(
        f'<text x="0" y="{total_h - 12}" '
        f'font-size="10" fill="#B4B2A9" '
        f'font-family="-apple-system,BlinkMacSystemFont,\'Segoe UI\',sans-serif">'
        f'OT = Officially Tracked (cowwid surveillance panel)</text>'
    )

    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{width}" height="{total_h}" style="display:block">'
        + "".join(lines_svg)
        + "".join(nodes_svg)
        + "".join(leg_svg)
        + "</svg>"
    )
    html = (
        f'<!DOCTYPE html><html><head>'
        f'<style>body{{margin:0;padding:8px 4px;background:transparent}}</style>'
        f'</head><body>{svg}</body></html>'
    )
    return html, total_h

# ── Public entry point ────────────────────────────────────────────────────────

def render_panel_tree(
    selected_variants: list[str],
    yaml_variants: list[str],
    pango_loader,
    scanner_results: dict | None = None,
    cooc_only: set | None = None,
):
    """
    Render the phylogenetic variant tree for the Abundance & Co-occurrence tab.

    Args:
        selected_variants: Variants currently in the deconv panel
        yaml_variants: Officially tracked variants (cowwid OT list)
        pango_loader: PangoLoader instance
        scanner_results: Scanner candidate results (v1, pass None for now)
        cooc_only: Cooc-only variant set (v1, pass None for now)
    """
    selected_set = set(selected_variants)
    cooc_set     = cooc_only or set()
    yaml_set     = set(yaml_variants)

    raw = pango_loader.get_raw_data()
    parent_map = {
        v: e.get("parent", "")
        for v, e in raw.items()
        if e.get("parent")
    }

    # Patch recombinant aliases that have no parent in pango_summary
    RECOMBINANT_PARENT = {
        "XDV": "JN.1",
        "XFG": "JN.1",
        "XEC": "JN.1",
        "XBB": "BA.2",
        "XBB.1": "XBB",
        "XBB.2": "XBB",
    }
    for alias, par in RECOMBINANT_PARENT.items():
        if alias in raw and alias not in parent_map:
            parent_map[alias] = par

    # Patch missing intermediate nodes via name-prefix inference
    all_known = set(raw.keys())
    for v in list(selected_set) + list(yaml_set):
        if v in parent_map:
            continue
        if "." not in v:
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

    result = _build_spine(
        selected_set, cooc_set, yaml_set,
        scanner_results, parent_map,
    )
    if not result:
        st.caption("Select at least one variant to build the tree.")
        return

    children, root, needed, kind_of, collapse = result
    _html, _h = _build_svg(children, root, kind_of, collapse)
    components.html(_html, height=_h + 20, scrolling=True)
