"""
components/abundance_cooc_tree.py

The variant tree of the Abundance & Co-occurrence tab, as rows:

  tree_rows()          the panel (and scanner findings) with their pango
                       ancestors, rooted at A; recombinant families (XEC, XFG,
                       XDV → NB.1.8.1 …) under their own "Recombinants" root,
                       each with its parents from pango's alias key.
  render_panel_tree()  the structural tree of the left column (panel only).

The Co-occurrence results page draws the same rows with per-city evidence
(components.variants_table.build), so both trees look alike.
"""
from __future__ import annotations
from collections import defaultdict

import streamlit as st
import streamlit.components.v1 as components


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
    # the root: the top of B's chain — A in pango naming (_parent_map puts B
    # under A), else B itself
    root = "B"
    _seen_r = set()
    while parent_map.get(root) and root not in _seen_r:
        _seen_r.add(root)
        root = parent_map[root]
    needed.add(root)

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
            children[root].append(v)  # recombinant roots: split off by tree_rows

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


def _parent_map(raw: dict, names: set) -> dict:
    """pango parent of every lineage; names missing from the summary are
    hung under their dotted-name parent."""
    parent_map = {v: e.get("parent", "") for v, e in raw.items() if e.get("parent")}
    # The Nextclade tree is rooted on the reference (Wuhan-Hu-1, lineage B)
    # and hangs A under B; in pango naming A is the root and B descends
    # from it — show it that way.
    if "A" in raw and parent_map.get("A") == "B":
        parent_map.pop("A")
        parent_map["B"] = "A"
    all_known = set(raw.keys())
    for v in names:
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
    return parent_map


def recombinant_parents(path=None) -> dict:
    """{recombinant: [parent, parent]} from pango-designation's alias_key.json
    (app/data/alias_key.json); {} if the file is missing. "*" = any
    descendant, as pango writes it."""
    import json
    from pathlib import Path
    path = Path(path) if path else Path(__file__).resolve().parents[1] / "data" / "alias_key.json"
    try:
        raw = json.loads(path.read_text())
    except Exception:
        return {}
    out = {}
    for k, v in raw.items():
        if k.startswith("X") and isinstance(v, list):
            out[k] = list(dict.fromkeys(x for x in v if x))   # de-duplicated, in order
    return out


def _parents_tip(x: str, parents: dict, depth: int = 0) -> str:
    """'XDV = XDE × JN.1 · XDE = GW.5.1 × FL.13.4' — recursing into parents
    that are recombinants themselves."""
    ps = parents.get(x)
    if not ps or depth > 3:
        return ""
    txt = f"{x} = {' × '.join(ps)}"
    for p in ps:
        sub = _parents_tip(p.rstrip("*"), parents, depth + 1)
        if sub:
            txt += " · " + sub
    return txt


def tree_rows(selected_variants, yaml_variants, pango_loader, findings=(),
              recomb_parents=None) -> list:
    """The panel tree as table rows (for the Variants table): one dict per
    row, depth-first — {node, label, depth, kind, last, guides, has_children,
    tip}. kind: panel / panel_ot / finding / yaml / spine / group. guides[k] = a
    vertical line continues at depth column k. Recombinant families get their
    own root row ("Recombinants", kind group) below the main tree (rooted at
    A); with
    recomb_parents (recombinant_parents()) the family's first row gets
    "parents" (e.g. "LF.7 × LP.8.1.2") and "parents_tip" (the full chain)."""
    raw = pango_loader.get_raw_data()
    selected_set, yaml_set = set(selected_variants), set(yaml_variants)
    finding_set = set(findings or ()) - selected_set
    parent_map = _parent_map(raw, selected_set | yaml_set | finding_set)
    res = _build_spine(selected_set, yaml_set, parent_map, recomb_set=set(),
                       finding_set=finding_set)
    if not res:
        return []
    children, root, _needed, kind_of, collapse = res

    def order(x):
        k = kind_of(x)
        return (0 if k in ("panel", "panel_ot", "finding") else 1 if k == "yaml" else 2, x)

    rows, seen = [], set()

    def walk(v, depth, guides, last):
        if v in seen:
            return
        seen.add(v)
        kind = kind_of(v)
        label, real = collapse(v) if (kind == "spine" and depth > 0) else (v, v)
        ch = sorted(children.get(real, []), key=order)
        rows.append({"node": real, "label": label, "depth": depth, "kind": kind,
                     "last": last, "guides": list(guides), "has_children": bool(ch)})
        for i, c in enumerate(ch):
            walk(c, depth + 1, guides + [not last] if depth >= 1 else [], i == len(ch) - 1)

    # recombinant families (XEC, XFG, XDV → NB.1.8.1 …) have two parents, so
    # they don't hang under the main tree: they get their own root below it
    recs = sorted([c for c in children.get(root, []) if c.startswith("X")
                   and not parent_map.get(c)], key=order)
    children[root] = [c for c in children.get(root, []) if c not in recs]
    if children[root]:
        walk(root, 0, [], True)
    if recs:
        rows.append({"node": "__recombinants__", "label": "Recombinants",
                     "tip": "Recombinant families — mixed ancestry from two parent lineages, "
                            "so they are not placed in the main tree",
                     "depth": 0, "kind": "group", "last": True, "guides": [],
                     "has_children": True})
        for i, c in enumerate(recs):
            n0 = len(rows)
            walk(c, 1, [], i == len(recs) - 1)
            ps = (recomb_parents or {}).get(c)
            if ps and len(rows) > n0:          # the family's first row
                rows[n0]["parents"] = " × ".join(ps)
                rows[n0]["parents_tip"] = ("Recombinant: " + _parents_tip(c, recomb_parents))
    return rows


def render_panel_tree(selected_variants: list[str], yaml_variants: list[str],
                      pango_loader, **_unused) -> None:
    """The structural tree of the left column: the panel with its ancestors,
    same look as the Variants view (recombinants under their own root, with
    their parents). Updates live with the selection; no results needed."""
    from components.variants_table import build_tree
    if not selected_variants:
        st.caption("Select at least one variant to build the tree.")
        return
    if "acooc_recomb_parents" not in st.session_state:
        st.session_state["acooc_recomb_parents"] = recombinant_parents()
    rows = tree_rows(selected_variants, yaml_variants, pango_loader,
                     recomb_parents=st.session_state["acooc_recomb_parents"])
    if not rows:
        st.caption("Select at least one variant to build the tree.")
        return
    html_, height = build_tree(rows, ot=yaml_variants)
    components.html(html_, height=height, scrolling=False)