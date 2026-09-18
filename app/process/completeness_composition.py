"""Completeness composition — per-date, normalized to 0-100%, split into four
categories, with empty/low-read dates dropped.

  explained : reads the panel accounts for (matched)          → green
  addable   : gap reads the scanner resolved to a clade       → red (add these)
  novel     : gap reads the scanner called novel              → violet (investigate)
  noise     : the rest of the gap (recurrent / homoplastic)   → grey

The green height IS the completeness for that date; the other bands show what the
unexplained gap is made of.

Addable is taken DIRECTLY from the scanner's own per-finding read accounting
(`reads_by_date` on each resolved clade + its sub-findings) and novel from the
scanner's novel `reads_by_date`. Completeness therefore reflects 100% of what the
scanner surfaced — a read the scanner counted into a confirmed finding is red, not
grey — and it does not re-classify anything: `classify_pattern` and the scanner's
detection are untouched. Reads the scanner did NOT attribute to a finding
(fingerprint < 2, unresolved, or resolved-without-a-discriminating-block) are not
in any `reads_by_date`, so they remain in noise — the graph never invents addable
signal the scanner didn't find.

`noise` is the gap remainder (unexplained − addable − novel), never a threshold,
so nothing leaks into or out of it.

Legacy fallback: if the scanner result predates `reads_by_date` (e.g. a stale
Redis result from an old worker), we fall back to the previous 80%-overlap
matching so nothing crashes. Each row is tagged `_addable_source` so a diagnostic
can tell which path ran.
"""

from typing import Dict, List, Set

_NOVEL_MATCH_FRACTION = 0.8


# ── scanner attribution (preferred path) ────────────────────────────────

def _has_reads_by_date(scanner_result: dict) -> bool:
    """True if the scanner result carries per-date read attribution."""
    for c in scanner_result.get("resolved_clade", []):
        if "reads_by_date" in c:
            return True
        for sf in c.get("sub_findings", []):
            if "reads_by_date" in sf:
                return True
    if "reads_by_date" in scanner_result.get("novel", {}):
        return True
    return False


def _addable_by_date(scanner_result: dict) -> Dict[str, int]:
    """Reads the scanner attributed to confirmed findings, per date.

    Sums each resolved clade's own `reads_by_date` plus those of its
    sub-findings. Parent and child slots are disjoint (each pattern is
    assigned to exactly one node in the scanner), so nothing is double
    counted; summing the whole confirmed branch gives its full read total.
    """
    add: Dict[str, int] = {}
    for c in scanner_result.get("resolved_clade", []):
        for d, n in (c.get("reads_by_date") or {}).items():
            add[d] = add.get(d, 0) + int(n)
        for sf in c.get("sub_findings", []):
            for d, n in (sf.get("reads_by_date") or {}).items():
                add[d] = add.get(d, 0) + int(n)
    return add


def _novel_by_date(scanner_result: dict) -> Dict[str, int]:
    return {d: int(n)
            for d, n in (scanner_result.get("novel", {}).get("reads_by_date") or {}).items()}


# ── legacy 80%-overlap matching (fallback only) ─────────────────────────

def _resolved_groups(scanner_result: dict) -> List[Set[str]]:
    groups: List[Set[str]] = []
    for c in scanner_result.get("resolved_clade", []):
        for b in c.get("member_blocks", []):
            g = set(b.get("discriminating", []))
            if len(g) >= 2:
                groups.append(g)
        for sf in c.get("sub_findings", []):
            for b in sf.get("member_blocks", []):
                g = set(b.get("discriminating", []))
                if len(g) >= 2:
                    groups.append(g)
    return groups


def _novel_patterns(scanner_result: dict) -> List[Set[str]]:
    return [set(p.get("mutations", []))
            for p in scanner_result.get("novel", {}).get("top_patterns", [])]


def _legacy_gap_counts(cooc_result: dict, scanner_result: dict,
                       dates: List[str]) -> (Dict[str, int], Dict[str, int]):
    """Old behaviour: re-classify each gap pattern by 80% overlap with a
    finding's discriminating group (addable) or a novel pattern (novel)."""
    resolved_groups = _resolved_groups(scanner_result)
    novel_pats = _novel_patterns(scanner_result)
    ups = cooc_result.get("unexplained_patterns", [])

    def gap_category(pat: Set[str]) -> str:
        for g in resolved_groups:
            if len(pat & g) / len(g) >= _NOVEL_MATCH_FRACTION:
                return "addable"
        for np_ in novel_pats:
            if np_ and len(pat & np_) / len(np_) >= _NOVEL_MATCH_FRACTION:
                return "novel"
        return "noise"

    add_d = {d: 0 for d in dates}
    nov_d = {d: 0 for d in dates}
    for p in ups:
        d = p.get("date")
        if d not in add_d:
            continue
        pat = set(p.get("confirmed_present", []))
        if len(pat) < 2:
            continue
        cnt = int(p.get("count", 0))
        cat = gap_category(pat)
        if cat == "addable":
            add_d[d] += cnt
        elif cat == "novel":
            nov_d[d] += cnt
    return add_d, nov_d


# ── main entry point ────────────────────────────────────────────────────

def compute_completeness_composition(cooc_result: dict,
                                     scanner_result: dict = None,
                                     min_reads: int = 1000) -> List[Dict]:
    """Per-date normalized composition rows.

    Returns list of dicts (empty/low-read dates dropped):
      {date, explained, addable, novel, noise, total,
       explained_pct, addable_pct, novel_pct, noise_pct, _addable_source}
    Percentages are of that date's total and sum to 1.0.
    """
    dates = cooc_result.get("dates", [])
    matched = cooc_result.get("matched_counts", [])
    unexpl = cooc_result.get("unexplained_counts", [])
    scanner_result = scanner_result or {}

    if _has_reads_by_date(scanner_result):
        source = "scanner_attribution"
        add_by_date = _addable_by_date(scanner_result)
        nov_by_date = _novel_by_date(scanner_result)
    else:
        source = "legacy_overlap"
        add_by_date, nov_by_date = _legacy_gap_counts(
            cooc_result, scanner_result, dates)

    rows: List[Dict] = []
    for i, d in enumerate(dates):
        m = int(matched[i]) if i < len(matched) else 0
        u = int(unexpl[i]) if i < len(unexpl) else 0
        total = m + u
        if total < min_reads:
            continue  # drop empty / low-read dates
        add = int(add_by_date.get(d, 0))
        nov = int(nov_by_date.get(d, 0))
        # addable + novel can't exceed the gap; clamp defensively so noise
        # (the remainder) never goes negative.
        add = min(add, u)
        nov = min(nov, max(0, u - add))
        noi = max(0, u - add - nov)
        rows.append({
            "date": d,
            "explained": m, "addable": add, "novel": nov, "noise": noi,
            "total": total,
            "explained_pct": m / total,
            "addable_pct": add / total,
            "novel_pct": nov / total,
            "noise_pct": noi / total,
            "_addable_source": source,
        })
    return rows