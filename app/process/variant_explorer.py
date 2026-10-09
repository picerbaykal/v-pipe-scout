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


def marker_funnel(variant: str, pango_loader, panel: Optional[List[str]] = None) -> Dict:
    """The ★ marker rule step by step for one lineage (2026-10-05) — the same
    steps as process.cooc.specific_markers, with what each step removed:

      signature   every substitution of the lineage (deletions left out)
      by_panel    {mutation: [panel variants carrying it]} — step 2 drops these
                  (the lineage's own descendants in the panel don't count)
      rows        for the rest: [{mutation, outside, forgiven, examples}] —
                  outside = carriers outside the family after option C;
                  forgiven = recombinants made from it that inherited it
      markers     outside <= out_other_max (most specific first)
      near        rejected by 1–5 outsiders over the limit: where the threshold
                  might be cutting off a real marker
      settings    {out_other_max, panel}
    """
    from process.cooc import _carrier_index, _check_cfg, _children_map, _lineage_family
    from process.recombinants import get_recombinants
    _raw, sigs, parent, _tree = _indices(pango_loader)
    lim = int(_check_cfg()["out_other_max"])
    sig = sorted(sigs.get(variant, set()), key=_pos)
    fam = _lineage_family(variant, parent, _children_map(parent))
    panel = [w for w in dict.fromkeys(panel or []) if w != variant]
    others = [w for w in panel if w not in fam]
    by_panel = {}
    for m in sig:
        who = [w for w in others if m in sigs.get(w, set())]
        if who:
            by_panel[m] = who
    rest = [m for m in sig if m not in by_panel]
    idx = _carrier_index(sigs)
    rc = get_recombinants(parent, sigs)
    own = {rc.rec_root(variant)} if variant in parent else set()
    rows = []
    for m in rest:
        outside = idx.get(m, set()) - fam
        forgiven = sorted(l for l in outside if rc.tolerated(l, variant, m, own))
        bad = sorted(outside - set(forgiven))
        rows.append({"mutation": m, "outside": len(bad), "forgiven": forgiven,
                     "examples": bad[:4]})
    rows.sort(key=lambda r: (r["outside"], _pos(r["mutation"])))
    markers = [r["mutation"] for r in rows if r["outside"] <= lim]
    near = [r for r in rows if lim < r["outside"] <= lim + 5]
    return {"variant": variant, "signature": sig, "by_panel": by_panel, "rows": rows,
            "markers": markers, "near": near,
            "settings": {"out_other_max": lim, "panel": panel}}


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
                  dates: List[str]) -> Dict:
    """One city's answer from the worker's counts.

    per_date: {date: {marker: [cov, hit, link_n, link_ok]}} (accumulate_check_stats)
    dates:    every sampling date in the window (also those with no reads on it)

    The vote (state, n_present, n_measured) adds up the reads of every date
    (2026-10-09: the chosen window; from 2 Oct to 9 Oct the last 5 samples).

    Returns {state, n_present, n_measured, n_markers, timeline, day_detail}
      state: present / absent / mixed / not_covered / no_marker
             (the check's own vote, process.cooc.check_verdicts)
      timeline: [[date, mark]], one mark per day: the cells' vote on that
             day's reads alone (2026-10-09; until then the reads were pooled
             over the markers, so one marker carried by an unrelated lineage
             could make a day "present": XFZ and LP.8's 1954A):
             present    the vote says present (>= min_present markers)
             absent     the vote says absent
             mixed      >= 2 markers decided and they disagree
             thin       reads, but fewer than 2 markers (1 for a lineage with
                        one marker) decided that day: one day alone can't
                        tell; the cell pools the last samples
             uncovered  no marker with >= min_cov reads that day
             The link test ("are the reads carrying the marker this lineage's
             reads?") is a property of the marker, not of the day, and one day
             rarely has the 20 reads it needs; so a day uses each marker's link
             pooled over the whole window, and its own reads for the share.
      day_detail: {date: "14808C 8 % of 140 reads ✓ · …"} for the squares'
             hover: why the day got its mark
    """
    from process.cooc import _check_cfg, check_verdicts
    c = _check_cfg()
    if not markers:
        return {"state": "no_marker", "n_present": 0, "n_measured": 0,
                "n_markers": 0, "timeline": [], "day_detail": {}}
    tl, detail = [], {}
    link = {m: [0, 0] for m in markers}          # [link_n, link_ok] over the window
    for cells in (per_date or {}).values():
        for m in markers:
            x = cells.get(m)
            if x:
                link[m][0] += x[2]
                link[m][1] += x[3]
    need = min(2, len(markers))
    for d in sorted(set(dates or []) | set(per_date or {})):
        cells = (per_date or {}).get(d, {})
        if not any((cells.get(m) or [0])[0] >= c["min_cov"] for m in markers):
            tl.append([d, "uncovered"])
            continue
        # this day's reads, the window's link
        day = {m: [x[0], x[1], link[m][0], link[m][1]]
               for m, x in cells.items() if m in link and x}
        r = check_verdicts(markers, {d: day}, c)
        if r["n_measured"] < need:
            mark = "present" if r["verdict"] == "present" else "thin"
        else:
            mark = r["verdict"] if r["verdict"] in ("present", "absent") else "mixed"
        tl.append([d, mark])
        detail[d] = _day_detail(markers, day, r, c)
    pd_ = per_date or {}
    res = check_verdicts(markers, pd_)
    state = res["verdict"]          # present / absent / mixed / not_covered
    return {"state": state, "n_present": res["n_present"], "n_measured": res["n_measured"],
            "n_markers": res["n_markers"], "timeline": tl,
            "day_detail": detail, "detail": _marker_detail(markers, pd_, res, c)}


def _day_detail(markers, cells, r, c, n_max=8) -> str:
    """One day, every marker: '14808C 8 % of 140 reads ✓ · 21249A 0 % of 900
    reads ✗ · 24604G 3 % of 300 reads (1–5 %) · 4927T 12 reads (under 100)'."""
    out = []
    for m in markers:
        x = (r.get("markers") or {}).get(m) or {}
        cov, f, st = x.get("cov", 0), x.get("freq"), x.get("status")
        if not cov:
            out.append(f"{m} no reads")
            continue
        num = f"{m} {f * 100:.0f} % of {cov:,} reads"
        if st == "present":
            out.append(num + " ✓")
        elif st == "absent":
            out.append(num + " ✗")
        elif cov < c["min_cov"]:
            out.append(f"{m} {cov} reads (under {c['min_cov']})")
        elif f < c["present_freq"]:
            out.append(num + f" ({c['absent_freq'] * 100:g}–{c['present_freq'] * 100:g} %)")
        else:
            lk_n = (cells.get(m) or [0, 0, 0, 0])[2]
            out.append(num + (f" (link test: {lk_n} reads in the window, needs "
                              f"{c['min_link']})" if lk_n < c["min_link"]
                              else " (its reads don't look like this lineage)"))
    more = f" · +{len(out) - n_max} more" if len(out) > n_max else ""
    return " · ".join(out[:n_max]) + more


def _marker_detail(markers, per_date, res, c) -> List[Dict]:
    """Each marker's pooled numbers and which rule decided it, flagging the
    ones that just missed a rule ("near"): within 20 % of a read count or of
    the present share, or within 5 points of the link share."""
    out = []
    for m in markers:
        link_n = sum((per_date.get(d, {}).get(m) or [0, 0, 0, 0])[2] for d in per_date)
        r = (res.get("markers") or {}).get(m) or {}
        cov, f, link, st = r.get("cov", 0), r.get("freq"), r.get("link"), r.get("status")
        why, near = "", ""
        if cov < c["min_cov"]:
            why = f"under {c['min_cov']} reads"
            if cov >= 0.8 * c["min_cov"]:
                near = f"{cov} reads, needs {c['min_cov']}"
        elif f is not None and f < c["absent_freq"]:
            why = f"absent (< {c['absent_freq'] * 100:g} %)"
        elif f is not None and f < c["present_freq"]:
            why = f"between {c['absent_freq'] * 100:g} and {c['present_freq'] * 100:g} %"
            if f >= 0.8 * c["present_freq"]:
                near = f"{f * 100:.1f} %, present needs {c['present_freq'] * 100:g} %"
        elif link_n < c["min_link"]:
            why = f"only {link_n} reads cover a neighbour (needs {c['min_link']})"
            if link_n >= 0.8 * c["min_link"]:
                near = why
        elif link is not None and link < c["link_min"]:
            why = (f"{link * 100:.0f} % of its reads look like the lineage nearby "
                   f"(needs {c['link_min'] * 100:g} %)")
            if link >= c["link_min"] - 0.05:
                near = why
        else:
            why = "present"
        out.append({"marker": m, "cov": cov, "freq": f, "link": link, "link_n": link_n,
                    "status": st, "why": why, "near": near})
    return out