"""Completeness composition — per-date, normalized to 0-100%, split into four
categories, with empty/low-read dates dropped.

  explained : reads the panel accounts for (matched)          → green
  addable   : gap reads the scanner resolved to a clade       → amber (add these)
  novel     : gap reads matching a novel pattern              → violet (investigate)
  noise     : the rest of the gap (recurrent / homoplastic)   → grey

The green height IS the completeness for that date; the other bands show what the
unexplained gap is made of. Built from the EXISTING cooc result + scanner result
— nothing is re-classified, `classify_pattern` and the scanner are untouched.
The scanner is the arbiter of what in the gap is coherent (addable/novel) vs
noise; the noise band is simply the remainder (no threshold defined here).
"""

from typing import Dict, List, Set

_NOVEL_MATCH_FRACTION = 0.8


def _resolved_groups(scanner_result: dict) -> List[Set[str]]:
    """Discriminating mutation groups of each resolved finding (from its member
    blocks). A gap pattern is 'addable' if it contains such a group — the same
    matching the scanner uses, which works (subset-of-observed_mutations does
    not, because a single read is rarely a subset of the whole clade union)."""
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
                                     min_reads: int = 1000) -> List[Dict]:
    """Per-date normalized composition rows.

    Returns list of dicts (empty/low-read dates dropped):
      {date, explained, addable, novel, noise, total,
       explained_pct, addable_pct, novel_pct, noise_pct}
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
        # addable: contains a resolved finding's discriminating group (>=80%)
        for g in resolved_groups:
            if len(pat & g) / len(g) >= _NOVEL_MATCH_FRACTION:
                return "addable"
        for np_ in novel_pats:
            if np_ and len(pat & np_) / len(np_) >= _NOVEL_MATCH_FRACTION:
                return "novel"
        return "noise"

    add_d = {d: 0 for d in dates}
    nov_d = {d: 0 for d in dates}
    noi_d = {d: 0 for d in dates}
    for p in ups:
        d = p.get("date")
        if d not in add_d:
            continue
        pat = set(p.get("confirmed_present", []))
        if len(pat) < 2:
            continue
        cat = gap_category(pat)
        cnt = int(p.get("count", 0))
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
        add = add_d[d]
        nov = nov_d[d]
        # noise = the gap remainder (guard against double-counted overflow)
        noi = max(0, u - add - nov)
        rows.append({
            "date": d,
            "explained": m, "addable": add, "novel": nov, "noise": noi,
            "total": total,
            "explained_pct": m / total,
            "addable_pct": add / total,
            "novel_pct": nov / total,
            "noise_pct": noi / total,
        })
    return rows