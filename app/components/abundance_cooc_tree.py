"""
components/abundance_cooc_tree.py

Phylogenetic variant tree for the Abundance & Co-occurrence tab.

Shows the global variant panel rooted at B. Two uses:

  structural (left column): every selected variant is black.
  per city (Co-occurrence results): each selected variant is coloured by its
  read evidence in that city, and scanner findings are added as red nodes:

    green        confirmed        specific evidence on >= 2 days
    light green  seen on 1 day    possible artefact, watch
    orange       not found        markers covered, absent
    purple       inconsistent     markers disagree (goes with the old check)
    grey         no ★ marker      nothing specific to test (e.g. KP.2, KP.3)
    dashed grey  no data          too few reads on its markers
    red          not in panel     scanner finding

Officially-tracked (cowwid) variants get an "OT" badge; tracked-but-not-selected
variants are hollow blue; spine/structural nodes are small grey circles. The
legend lists only the states present in the tree; hovering a node names its state.
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

# state → (colour, legend / hover label). Orange matches the scanner list's
# "Not found in wastewater" band, red its "Not in your panel" band.
STATUS_C = {
    "confirmed":    ("#16a34a", "confirmed"),
    "one_day":      ("#86c98a", "seen on 1 day"),
    "not_found":    ("#b45309", "not found in WW"),
    "inconsistent": ("#7c3aed", "inconsistent"),
    "no_marker":    ("#9ca3af", "no ★ marker"),
    "no_data":      ("#9ca3af", "no data"),
    "finding":      ("#dc2626", "found, not in panel"),
    "pending":      ("#1f2430", "selected"),
}
_DASHED = {"no_data"}          # drawn as an empty dashed circle


def _status_of(v: str, variant_status: dict | None, default: str) -> str:
    s = (variant_status or {}).get(v, default)
    return s if s in STATUS_C else default


def _status_color(v: str, variant_status: dict | None, default: str = "pending") -> str:
    return STATUS_C[_status_of(v, variant_status, default)][0]


def _ancestors(v: str, parent_map: dict) -> list[str]:
    path = []
    while v in parent_map:
        v = parent_map[v]
        path.append(v)
    return path


def _build_spine(selected_set: set, yaml_set: set, parent_map: dict, recomb_set: set = None,
                 finding_set: set = None):
    recomb_set = recomb_set or set()
    finding_set = finding_set or set()
    selected_set = set(selected_set) | set(finding_set)
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
        if v in finding_set:
            return "finding"
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
               recombinant_set=None, variant_status=None, default_status="pending",
               per_city=False, legend_states=None, status_detail=None):
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
            0 if kind_of(x) in ("panel", "panel_ot", "finding") else
            1 if kind_of(x) == "yaml" else 2, x))
        for c in ch:
            assign_rows(c, depth + 1)

    assign_rows(root, 0)
    row_y = {rows[i][0]: i * ROW_H + ROW_H // 2 for i in range(len(rows))}
    # legend: the states used in this tree (plus any from the recombinant
    # section, passed in), wrapped onto as many rows as the width needs
    _order = list(STATUS_C)
    _leg_items, _leg_rows, _lx = [], 0, 0
    for st_ in sorted(set(legend_states or ()), key=_order.index):
        w = len(STATUS_C[st_][1]) * 6 + 26
        if _lx and _lx + w > width:
            _leg_rows, _lx = _leg_rows + 1, 0
        _leg_items.append((st_, _lx, _leg_rows))
        _lx += w
    n_leg = (_leg_rows + 1) if _leg_items else 0
    leg_top = len(rows) * ROW_H + 14
    total_h = leg_top + 18 * n_leg + 18 + 22

    def node_state(v, kind):
        if kind == "finding":
            return "finding"
        if kind in ("panel", "panel_ot"):
            return _status_of(v, variant_status, default_status)
        return None

    def node_color(v, kind):
        st_ = node_state(v, kind)
        return STATUS_C[st_][0] if st_ else C[kind]


    lines_svg, nodes_svg = [], []

    for v, _, real_v, depth, kind in rows:
        ch = sorted(children.get(real_v, []), key=lambda x: (
            0 if kind_of(x) in ("panel", "panel_ot", "finding") else
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
        state = node_state(v, kind)
        is_spine = kind == "spine"
        filled = kind in ("panel", "panel_ot", "finding")
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

        # hover names the state (per-city tree) — colour is never the only cue
        _det = (status_detail or {}).get(real_v, "")
        _det = f" — {_det}" if _det else ""
        nodes_svg.append(f'<g><title>{real_v}: {STATUS_C[state][1]}{_det}</title>'
                         if state else '<g>')
        if filled and state in _DASHED:
            nodes_svg.append(
                f'<circle cx="{x}" cy="{y}" r="{r}" fill="white" stroke="{color}" '
                f'stroke-width="1.5" stroke-dasharray="2,2"/>')
        elif filled:
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
        nodes_svg.append('</g>')

    # ── legend (two rows: status swatches, then tracked/OT) ──────────────────
    leg_svg = [
        f'<line x1="0" y1="{leg_top - 8}" x2="{width}" y2="{leg_top - 8}" '
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

    # states: only the ones that appear, wrapped
    for st_, lx, ri in _leg_items:
        lc, ltxt = STATUS_C[st_]
        ly = leg_top + 18 * ri + 4
        if st_ in _DASHED:
            leg_svg.append(f'<circle cx="{lx + 5}" cy="{ly}" r="4" fill="white" '
                           f'stroke="{lc}" stroke-width="1.5" stroke-dasharray="2,2"/>')
        else:
            leg_svg.append(_dot(lx + 5, ly, lc))
        leg_svg.append(_txt(lx + 13, ly, ltxt))

    # row 2: tracked (hollow) + OT chip
    row2_y = leg_top + 18 * n_leg + 4
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

    _foot = ("colour = evidence from reads in this city · hover a node"
             if per_city else "your panel (structure only)")
    leg_svg.append(
        f'<text x="0" y="{total_h - 8}" font-size="9" fill="#B4B2A9" '
        f'font-family="-apple-system,BlinkMacSystemFont,\'Segoe UI\',sans-serif">'
        f'{_foot}</text>'
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
    findings: list[str] | None = None,
    per_city: bool = False,
    status_detail: dict | None = None,
):
    """
    Render the global variant tree.

    Structural (per_city=False, left column): selected variants are black.
    Per city (per_city=True): variant_status {variant: state} colours each
        selected variant (states: see STATUS_C; a variant missing from it is
        "no data"), and `findings` (scanner nodes not in the panel) are added
        as red nodes under their place in the tree. Pass the panel that was RUN
        as selected_variants, so the colours always match the results.
        status_detail {variant: text} is added to the hover (the evidence).
    (scanner_results, cooc_only, scanner_added_for kept for signature
    compatibility but unused.)
    """
    default_status = "no_data" if per_city else "pending"
    if not per_city:
        variant_status = {}
    variant_status = variant_status or {}
    finding_set = {f for f in (findings or []) if f not in set(selected_variants)}

    selected_set = set(selected_variants)
    yaml_set = set(yaml_variants)

    raw = pango_loader.get_raw_data()
    parent_map = {v: e.get("parent", "") for v, e in raw.items() if e.get("parent")}

    all_known = set(raw.keys())
    for v in list(selected_set) + list(yaml_set) + list(finding_set):
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

    _all_panel = set(selected_variants) | set(yaml_variants) | finding_set
    recomb_members = {v for v in _all_panel if _recomb_root(v) is not None}

    result = _build_spine(selected_set, yaml_set, parent_map, recomb_set=recomb_members,
                          finding_set=finding_set)

    # ── main bifurcating tree (non-recombinants) ──
    if result:
        children, root, needed, kind_of, collapse = result
        _states = {("finding" if v in finding_set
                     else _status_of(v, variant_status, default_status))
                    for v in selected_set | finding_set}
        _html, _h = _build_svg(children, root, kind_of, collapse,
                               variant_status=variant_status,
                               default_status=default_status, per_city=per_city,
                               legend_states=_states, status_detail=status_detail)
        components.html(_html, height=_h + 20, scrolling=True)
    elif not recomb_members:
        st.caption("Select at least one variant to build the tree.")
        return

    # ── separate Recombinants section (grouped by family), coloured by verdict ──
    selected_recomb = {v for v in list(selected_variants) + sorted(finding_set)
                       if v in recomb_members}
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
                _st = ("finding" if m in finding_set
                       else _status_of(m, variant_status, default_status))
                _mc, _ml = STATUS_C[_st]
                _dot_ch = "◌" if _st in _DASHED else "●"
                _lbl = f" · {_ml}" if per_city else ""
                lines.append(
                    f"{indent}<span title='{m}: {_ml}"
                    + (f" — {(status_detail or {}).get(m, '')}".replace("'", "&#39;")
                       if (status_detail or {}).get(m) else "")
                    + f"'><span style='color:{_mc};'>{_dot_ch}</span> "
                    f"<b style='color:{_mc};'>{m}</b>{badge}"
                    f"<span style='color:#9ca3af;font-size:11px;'>{_lbl}</span></span>")
            st.markdown(
                f"<div style='border:0.5px solid #e5e7eb;border-radius:8px;"
                f"padding:6px 10px;margin:4px 0;font-size:13px;'>"
                f"<span style='color:#8A4E82;font-weight:600;'>⧉ {froot} recombinant</span><br>"
                + "<br>".join(lines) + "</div>",
                unsafe_allow_html=True,
            )