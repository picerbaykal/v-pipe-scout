"""Variant explorer — on-demand look-up of ANY pango lineage.

Two parts (2026-10-02):
  1. investigate_variant — from the tree only: can co-occurrence tell this
     lineage apart (★ markers, the SAME rule as the panel check:
     process.cooc.specific_markers), and how it connects to others.
  2. check_in_data — the worker's read counts at its ★ markers in one city
     (tasks.run_cooc_variant_check_lapis) -> present / absent / mixed /
     not covered, plus one mark per sampling day for a tiny calendar.

Pure logic (no Streamlit, no IO). Not anchored to the panel or the tracked
list: works for any lineage the user picks, also an old one — its markers are
read and the answer is "absent" when the data covers them and they're not there.
"""

from typing import Dict, List, Optional, Set

_SUB = "ACGT"


def _pos(m: str) -> int:
    """'29774T' -> 29774."""
    d = "".join(ch for ch in m if ch.isdigit())
    return int(d) if d else 0


class _Tree:
    """Cheap parent/child index over the pango raw data."""
    def __init__(self, raw: dict):
        self.parent = {l: e.get("parent", "") for l, e in raw.items()}
        self.child: Dict[str, List[str]] = {}
        for l, p in self.parent.items():
            if p:
                self.child.setdefault(p, []).append(l)

    def children(self, v: str) -> List[str]:
        return sorted(self.child.get(v, []))

    def siblings(self, v: str) -> List[str]:
        p = self.parent.get(v, "")
        if not p:
            return []
        return sorted(c for c in self.child.get(p, []) if c != v)


# one index per pango loader (the page keeps one loader for its lifetime)
_CACHE: Dict[int, tuple] = {}


def _indices(pango_loader):
    k = id(pango_loader)
    if k not in _CACHE:
        raw = pango_loader.get_raw_data()
        sigs = {}
        for l in raw:
            s = {m for m in (pango_loader.get_signature(l) or ()) if m and m[-1] in _SUB}
            if len(s) >= 2:
                sigs[l] = s
        parent = {l: e.get("parent", "") for l, e in raw.items()}
        _CACHE.clear()
        _CACHE[k] = (raw, sigs, parent, _Tree(raw))
    return _CACHE[k]


def variant_markers(variant: str, pango_loader, panel: Optional[List[str]] = None) -> List[str]:
    """★ markers of `variant` as the panel check would compute them if it were
    added to `panel` (most specific first)."""
    from process.cooc import specific_markers
    _, sigs, parent, _ = _indices(pango_loader)
    vs = {v: sigs.get(v, set()) for v in list(dict.fromkeys(list(panel or []) + [variant]))}
    return specific_markers(vs, sigs, parent).get(variant, [])


def marker_blocks(variant: str, pango_loader, panel: Optional[List[str]] = None,
                  span: int = 300) -> List[Dict]:
    """Groups of the lineage's own mutations that one read can cover (within
    `span` bases), keeping the groups with a ★ marker — the input of the
    co-occurrence heatmap (components.scanner_heatmap) for ANY lineage, also
    one the scanner didn't report. Same shape as the scanner's member blocks:
    {member, discriminating, mut_star, mut_carriers, reads}."""
    _raw, sigs, _parent, _tree = _indices(pango_loader)
    sig = sorted(sigs.get(variant, set()), key=lambda m: int(_pos(m)))
    stars = set(variant_markers(variant, pango_loader, panel))
    groups, cur = [], []
    for m in sig:
        if cur and int(_pos(m)) - int(_pos(cur[0])) > span:
            groups.append(cur)
            cur = []
        cur.append(m)
    if cur:
        groups.append(cur)
    out = []
    for g in groups:
        if not stars & set(g):
            continue
        car = {m: sum(1 for s in sigs.values() if m in s) for m in g}
        out.append({"member": variant, "discriminating": g, "reads": 0,
                    "mut_star": {m: m in stars for m in g}, "mut_carriers": car})
    return out


def investigate_variant(variant: str, pango_loader,
                        panel: Optional[List[str]] = None) -> Dict:
    """Fact sheet for `variant`:
        {found, name, in_panel, is_recombinant, detectable, markers, reason,
         parent, siblings, children, oscillates_with}
    detectable = it has >= 1 ★ marker (same rule as the panel check)."""
    raw, sigs, _parent, tree = _indices(pango_loader)
    if variant not in raw:
        return {"found": False, "name": variant,
                "reason": "Not a known pango lineage in the current data."}
    sig = sigs.get(variant, set())
    panel = list(panel or [])
    try:
        markers = variant_markers(variant, pango_loader, panel)
    except Exception:
        markers = []
    parent = raw[variant].get("parent", "")
    siblings = tree.siblings(variant)

    if markers:
        reason = (f"{len(markers)} ★ marker{'s' if len(markers) > 1 else ''} "
                  f"(e.g. {', '.join(markers[:3])}) — mutations it and its "
                  f"sublineages carry, and at most a few other lineages.")
    elif sig:
        reason = ("no ★ marker — each of its mutations is also carried by another "
                  "panel variant or by many lineages outside its family, so reads "
                  "can't point to it. Quantify with deconvolution.")
    else:
        reason = "No signature mutations available for this lineage."

    # near-identical siblings (Jaccard >= 0.98): they can't be told apart
    osc = []
    for s in siblings:
        ss = sigs.get(s, set())
        if sig and ss and len(sig & ss) / len(sig | ss) >= 0.98:
            osc.append(s)

    return {
        "found": True, "name": variant,
        "in_panel": variant in set(panel),
        "is_recombinant": variant.startswith("X") and not parent,
        "detectable": bool(markers),
        "markers": markers,
        "reason": reason,
        "parent": parent,
        "siblings": siblings[:8],
        "children": tree.children(variant)[:8],
        "oscillates_with": osc,
    }


def check_in_data(markers: List[str], per_date: Dict[str, Dict[str, list]],
                  dates: List[str], recent: Optional[int] = None) -> Dict:
    """One city's answer from the worker's counts.

    per_date: {date: {marker: [cov, hit, link_n, link_ok]}} (accumulate_check_stats)
    dates:    every sampling date in the window (also those with no reads on it)

    recent:   when set, the vote (state, n_present, n_measured) uses only the
              last `recent` covered sampling days ("is it there now?"); the
              timeline still covers the whole window. recent_dates lists them.

    Returns {state, n_present, n_measured, n_markers, timeline, recent_dates}
      state: present / absent / mixed / not_covered / no_marker
             (the check's own vote, process.cooc.check_verdicts)
      timeline: [[date, mark]], mark per day pooled over the markers:
             present  >= present_freq of >= min_cov reads
             absent   <  absent_freq  of >= min_cov reads
             weak     in between
             uncovered < min_cov reads at its markers
    """
    from process.cooc import _check_cfg, check_verdicts
    c = _check_cfg()
    if not markers:
        return {"state": "no_marker", "n_present": 0, "n_measured": 0,
                "n_markers": 0, "timeline": [], "recent_dates": []}
    tl = []
    for d in sorted(set(dates or []) | set(per_date or {})):
        cells = (per_date or {}).get(d, {})
        cov = sum(cells.get(m, [0, 0])[0] for m in markers)
        hit = sum(cells.get(m, [0, 0])[1] for m in markers)
        # a day is covered when at least one marker has min_cov reads
        if not any(cells.get(m, [0])[0] >= c["min_cov"] for m in markers):
            mark = "uncovered"
        elif hit / cov >= c["present_freq"]:
            mark = "present"
        elif hit / cov < c["absent_freq"]:
            mark = "absent"
        else:
            mark = "weak"
        tl.append([d, mark])
    recent_dates = ([d for d, m in tl if m != "uncovered"][-recent:] if recent else [])
    pd_ = ({d: (per_date or {}).get(d, {}) for d in recent_dates} if recent
           else (per_date or {}))
    res = check_verdicts(markers, pd_)
    state = {"confirmed": "present", "not_found": "absent",
             "inconsistent": "mixed"}.get(res["verdict"], "not_covered")
    return {"state": state, "n_present": res["n_present"], "n_measured": res["n_measured"],
            "n_markers": res["n_markers"], "timeline": tl, "recent_dates": recent_dates}