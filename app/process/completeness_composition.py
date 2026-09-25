"""Completeness composition — per-date, normalized to 0-100%, split into
categories, with empty/low-read dates dropped.

  explained  : reads the panel accounts for exactly (matched)       → green
  near-panel : "panel variant + 1 change" — a single panel variant
               explains the read except at exactly one position
               (process.cooc.near_panel_label, decided per read in
               the worker)                                          → teal
               Shown as its OWN band, not added to explained: a thin flat
               band is noise; a growing one is a sublineage of a panel
               variant spreading (e.g. "XFG + 22896C").
  addable    : gap reads the scanner resolved to a clade            → red (add these)
  novel      : gap reads matching a novel pattern                   → blue (investigate)
  noise      : the rest of the gap                                  → grey

`panel_union` is {variant: signature} for the panel ("{pos}{alt}" substitutions);
a read's beyond-panel mutations are those its BEST-MATCHING panel variant lacks
(process.scanner.beyond_panel — the same definition the scanner uses, so an
extinct control in the panel can't hide a real variant). A plain set is still
accepted and treated as the old pooled union. A pattern with fewer than 2
beyond-panel mutations is not seen by the scanner; if its one outside mutation is a scanner-finding marker it is a
partial read of that finding (addable / novel), otherwise its near-panel reads
(near_count) are teal and the rest is noise.

Legacy results without `near_count` (scans from before 2026-09-25) keep the old
rule: every such read counts as near-panel.
"""

from collections import defaultdict
from typing import Dict, List, Set, Optional, Union

try:
    from process.scanner import beyond_panel
except Exception:                                   # pragma: no cover
    def beyond_panel(present, panel_sigs):
        present = set(present)
        if isinstance(panel_sigs, dict):
            if not panel_sigs:
                return present
            return present - max(panel_sigs.values(), key=lambda s: len(present & s))
        return present - set(panel_sigs or ())

_NOVEL_MATCH_FRACTION = 0.8
_NEAR_TOP = 3


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


def _star_markers(scanner_result: dict) -> Set[str]:
    """★ markers of the confirmed findings (block["mut_star"] from the scanner;
    older results without it fall back to every block mutation)."""
    out: Set[str] = set()
    for c in scanner_result.get("resolved_clade", []):
        for b in c.get("member_blocks", []):
            flags = b.get("mut_star")
            if flags is None:
                out |= set(b.get("discriminating", []))
            else:
                out |= {m for m, f in flags.items() if f}
    return out


def _novel_patterns(scanner_result: dict) -> List[Set[str]]:
    return [set(p.get("mutations", []))
            for p in scanner_result.get("novel", {}).get("top_patterns", [])]


def compute_completeness_composition(cooc_result: dict,
                                     scanner_result: dict = None,
                                     min_reads: int = 1000,
                                     panel_union: Optional[Union[Set[str], Dict[str, Set[str]]]] = None
                                     ) -> List[Dict]:
    """Per-date normalized composition rows.

    Returns list of dicts (empty/low-read dates dropped):
      {date, explained, near, addable, novel, noise, total,
       explained_pct, near_pct, addable_pct, novel_pct, noise_pct,
       near_top: [(label, reads), ...]  — largest "variant + 1 change" groups}
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

    # A lone beyond-panel mutation counts for a finding only if it is one of that
    # finding's ★ markers (a partial read of a real variant). A SHARED mutation
    # that merely appears in a finding's group, or in a small novel pattern (e.g.
    # XFG's 21653C inside a 2-mutation novel pattern), would otherwise repaint
    # every read of the variant that carries it. Novel is about a combination,
    # so a lone mutation is never attributed to it.
    addable_muts: Set[str] = _star_markers(scanner_result)
    novel_muts: Set[str] = set()

    add_d = {d: 0 for d in dates}
    nov_d = {d: 0 for d in dates}
    near_d = {d: 0 for d in dates}
    near_lab = {d: defaultdict(int) for d in dates}
    for p in ups:
        d = p.get("date")
        if d not in add_d:
            continue
        pat = set(p.get("confirmed_present", []))
        if len(pat) < 2:
            continue
        cnt = int(p.get("count", 0))
        if panel_union is not None:
            fp = beyond_panel(pat, panel_union)
            if len(fp) < 2:
                if fp & addable_muts:
                    add_d[d] += cnt
                elif fp & novel_muts:
                    nov_d[d] += cnt
                elif "near_count" in p:
                    ncnt = int(p.get("near_count", 0) or 0)
                    lab = p.get("near_label") or ""
                    # the ONE change vs the best-matching panel variant; if it is
                    # a scanner-finding marker the read is a partial read of that
                    # finding -> red/blue, independent of who else is in the panel
                    chg = lab.split(" + ", 1)[1] if " + " in lab else ""
                    if ncnt and chg and chg in addable_muts:
                        add_d[d] += ncnt
                    elif ncnt and chg and chg in novel_muts:
                        nov_d[d] += ncnt
                    else:
                        near_d[d] += ncnt
                        if ncnt and lab:
                            near_lab[d][lab] += ncnt
                    # the pattern's other reads fall through to noise below
                else:
                    near_d[d] += cnt          # legacy result: old rule
                continue
        cat = gap_category(pat)
        if cat == "addable":
            add_d[d] += cnt
        elif cat == "novel":
            nov_d[d] += cnt

    rows: List[Dict] = []
    for i, d in enumerate(dates):
        m = int(matched[i]) if i < len(matched) else 0
        u = int(unexpl[i]) if i < len(unexpl) else 0
        total = m + u
        if total < min_reads:
            continue  # drop empty / low-read dates
        near = near_d[d]
        add = add_d[d]
        nov = nov_d[d]
        noi = max(0, u - near - add - nov)
        top = sorted(near_lab[d].items(), key=lambda kv: -kv[1])[:_NEAR_TOP]
        rows.append({
            "date": d,
            "explained": m, "near": near, "addable": add, "novel": nov,
            "noise": noi, "total": total,
            "explained_pct": m / total,
            "near_pct": near / total,
            "addable_pct": add / total,
            "novel_pct": nov / total,
            "noise_pct": noi / total,
            "near_top": top,
        })
    return rows