"""Completeness composition — per-date, normalized to 0-100%, split into
categories, with empty/low-read dates dropped.

  explained : reads the panel accounts for exactly (matched)      → dark green
  near-panel: a panel variant + <2 stray mutations (fingerprint<2)→ light green
              (counts toward completeness — no beyond-panel signal; usually
              homoplasy / sequencing error; watch a rising trend as an early
              sublineage hint)
  addable   : gap reads the scanner resolved to a clade           → red (add these)
  novel     : gap reads matching a novel pattern                  → blue (investigate)
  noise     : the rest of the gap (fingerprint>=2, unresolved)    → grey

Completeness for a date is explained + near-panel (both are the panel's variants).
The remaining bands show what the genuinely-beyond-panel gap is made of — the
same fingerprint>=2 signal the scanner acts on, so grey here matches the scanner.

`panel_union` is the set of the panel variants' signature mutations ("{pos}{alt}",
substitutions). When it's None the near-panel split is skipped (legacy behaviour:
those reads stay in noise). Built from the EXISTING cooc + scanner result —
`classify_pattern` and the scanner are untouched.
"""

from typing import Dict, List, Set, Optional

_NOVEL_MATCH_FRACTION = 0.8


def _resolved_groups(scanner_result: dict) -> List[Set[str]]:
    """Discriminating mutation groups of each resolved finding (from its member
    blocks) — the same matching the scanner uses."""
    groups: List[Set[str]] = []
    for c in scanner_result.get("resolved_clade", []):
        for b in c.get("member_blocks", []):
            g = set(b.get("discriminating", []))
            if len(g) >= 2:
                groups.append(g)
    return groups


def _novel_patterns(scanner_result: dict) -> List[Set[str]]:
    return [set(p.get("mutations", []))
            for p in scanner_result.get("novel", {}).get("top_patterns", [])]


def compute_completeness_composition(cooc_result: dict,
                                     scanner_result: dict = None,
                                     min_reads: int = 1000,
                                     panel_union: Optional[Set[str]] = None) -> List[Dict]:
    """Per-date normalized composition rows.

    Returns list of dicts (empty/low-read dates dropped):
      {date, explained, near, addable, novel, noise, total,
       explained_pct, near_pct, addable_pct, novel_pct, noise_pct}
    Percentages are of that date's total and sum to 1.0.
    """
    dates = cooc_result.get("dates", [])
    matched = cooc_result.get("matched_counts", [])
    unexpl = cooc_result.get("unexplained_counts", [])
    ups = cooc_result.get("unexplained_patterns", [])
    scanner_result = scanner_result or {}

    resolved_groups = _resolved_groups(scanner_result)
    novel_pats = _novel_patterns(scanner_result)

    def gap_category(pat: Set[str]) -> str:
        for g in resolved_groups:
            if len(pat & g) / len(g) >= _NOVEL_MATCH_FRACTION:
                return "addable"
        for np_ in novel_pats:
            if np_ and len(pat & np_) / len(np_) >= _NOVEL_MATCH_FRACTION:
                return "novel"
        return "noise"

    # markers of the scanner's findings — a lone beyond-panel mutation that is one
    # of these is a partial-coverage read of a REAL (close) variant, not a stray.
    addable_muts: Set[str] = set().union(*resolved_groups) if resolved_groups else set()
    novel_muts: Set[str] = set().union(*novel_pats) if novel_pats else set()

    add_d = {d: 0 for d in dates}
    nov_d = {d: 0 for d in dates}
    noi_d = {d: 0 for d in dates}
    near_d = {d: 0 for d in dates}
    for p in ups:
        d = p.get("date")
        if d not in add_d:
            continue
        pat = set(p.get("confirmed_present", []))
        if len(pat) < 2:
            continue
        cnt = int(p.get("count", 0))
        if panel_union is not None:
            fp = pat - panel_union
            if len(fp) < 2:
                # A single beyond-panel mutation. If it is a MARKER of a scanner
                # finding, this is a partial-coverage read of that real variant
                # (e.g. XFG, close to the panel, where most reads catch only one
                # of its discriminating mutations) -> attribute to addable/novel,
                # NOT near-panel. Only a mutation that matches no finding is a
                # benign stray (near-panel homoplasy) that counts to completeness.
                if fp & addable_muts:
                    add_d[d] += cnt
                elif fp & novel_muts:
                    nov_d[d] += cnt
                else:
                    near_d[d] += cnt
                continue
        cat = gap_category(pat)
        if cat == "addable":
            add_d[d] += cnt
        elif cat == "novel":
            nov_d[d] += cnt
        else:
            noi_d[d] += cnt

    rows: List[Dict] = []
    for i, d in enumerate(dates):
        m = int(matched[i]) if i < len(matched) else 0
        u = int(unexpl[i]) if i < len(unexpl) else 0
        total = m + u
        if total < min_reads:
            continue  # drop empty / low-read dates
        # near-panel (a panel variant + 1 NON-marker stray mutation) counts as
        # explained — it is one of the panel's variants, the stray is homoplasy /
        # batch noise. Marker-carrying single mutations already went to addable/
        # novel above, so nothing real is hidden here. Folded into green to keep
        # the graph to four categories that match the scanner list.
        near = near_d[d]
        expl = m + near
        add = add_d[d]
        nov = nov_d[d]
        noi = max(0, u - near - add - nov)   # unresolved gap (fingerprint>=2)
        rows.append({
            "date": d,
            "explained": expl, "near": near, "addable": add, "novel": nov,
            "noise": noi, "total": total,
            "explained_pct": expl / total,
            "near_pct": near / total,   # kept for diagnostics; not a band
            "addable_pct": add / total,
            "novel_pct": nov / total,
            "noise_pct": noi / total,
        })
    return rows