"""
process/near_changes.py

"Panel variant ± 1 change" per city, for the Variants view (2026-10-01).

Completeness labels every unexplained read that differs from ONE panel
variant at exactly one position: "XFG + 22896C" (the read carries a mutation
XFG lacks) or "XFG − 23021G" (the read shows the reference where XFG has a
mutation). These reads are the amber band of the completeness graph. Here
they are counted per change, city and day, with the same day rule as scanner
findings (>= EVIDENCE_MIN_READS reads and >= EVIDENCE_MIN_SHARE of the day's
informative reads), and each change is placed in the pango tree: the
sublineages of that variant where the mutation was gained (+) or lost (−).

A single mutation is the weakest kind of evidence (one convergent mutation or
one PCR jackpot can make it), so the view shows these as hints to watch,
never as findings.
"""
from __future__ import annotations

import re
from typing import Dict, List, Optional, Set

_LABEL = re.compile(r"^(.+?) ([+−]) (\d+[ACGT])$")


def parse_label(label: str):
    """'XFG + 22896C' -> ('XFG', '+', '22896C'); anything else -> None."""
    m = _LABEL.match(label or "")
    return (m.group(1), m.group(2), m.group(3)) if m else None


def _day_total(res: dict) -> Dict[str, int]:
    return {str(d)[:10]: int(m) + int(u) for d, m, u in zip(
        res.get("dates", []) or [], res.get("matched_counts", []) or [],
        res.get("unexplained_counts", []) or [])}


def near_changes(cooc_by_city: Dict[str, dict], panel: List[str],
                 min_reads: int, min_share: float) -> Dict[str, List[dict]]:
    """{variant: [{variant, sign, mut, label, per_city: {city: days},
    reads: {city: reads}}]} for the panel variants, from the completeness
    results of each city. A day counts when the change has >= min_reads reads
    and >= min_share of the day's informative reads. Changes are listed when
    they count on >= 1 day in some city; "+" first, then by days."""
    panel_set = set(panel)
    acc: Dict[tuple, dict] = {}
    for city, res in (cooc_by_city or {}).items():
        if not res:
            continue
        tot = _day_total(res)
        per_day: Dict[tuple, Dict[str, int]] = {}
        for row in res.get("unexplained_patterns", []) or []:
            n = int(row.get("near_count", 0) or 0)
            p = parse_label(row.get("near_label", ""))
            if not n or not p or p[0] not in panel_set:
                continue
            d = str(row.get("date", ""))[:10]
            per_day.setdefault(p, {})
            per_day[p][d] = per_day[p].get(d, 0) + n
        for p, days in per_day.items():
            ok = [d for d, n in days.items()
                  if n >= min_reads and n >= min_share * tot.get(d, 0)]
            slot = acc.setdefault(p, {"variant": p[0], "sign": p[1], "mut": p[2],
                                      "label": f"{p[0]} {p[1]} {p[2]}",
                                      "per_city": {}, "reads": {}})
            slot["reads"][city] = sum(days.values())
            if ok:
                slot["per_city"][city] = len(ok)
    out: Dict[str, List[dict]] = {}
    for slot in acc.values():
        if not slot["per_city"]:
            continue
        out.setdefault(slot["variant"], []).append(slot)
    for v in out:
        out[v].sort(key=lambda s: (s["sign"] != "+", -max(s["per_city"].values()),
                                   -sum(s["reads"].values()), s["mut"]))
    return out


def recurring_gains(changes: List[dict], min_days: int) -> int:
    """How many '+' changes count on >= min_days days in some city (the tag)."""
    return sum(1 for c in changes
               if c["sign"] == "+" and max(c["per_city"].values(), default=0) >= min_days)


def where_in_tree(variant: str, sign: str, mut: str, sigs: Dict[str, Set[str]],
                  children: Dict[str, List[str]], limit: int = 3) -> Optional[dict]:
    """The sublineages of `variant` where `mut` was gained ('+': carries it,
    its parent doesn't) or lost ('−': lacks it, its parent has it).
    Returns {"roots": [...up to limit, shallowest first], "n_roots": int,
    "n_lineages": int (roots and all their descendants)} or None when no
    designated sublineage of the variant has that change."""
    def has(l):
        return mut in (sigs.get(l) or ())

    roots, n_all = [], 0
    stack = [(c, variant, 1) for c in children.get(variant, [])]
    while stack:
        node, par, depth = stack.pop()
        changed = (has(node) and not has(par)) if sign == "+" else (not has(node) and has(par))
        if changed:
            roots.append((depth, node))
            sub, todo = 0, [node]
            while todo:
                x = todo.pop()
                sub += 1
                todo.extend(children.get(x, []))
            n_all += sub
            continue                       # its descendants follow it
        stack.extend((c, node, depth + 1) for c in children.get(node, []))
    if not roots:
        return None
    roots.sort()
    return {"roots": [n for _, n in roots[:limit]], "n_roots": len(roots),
            "n_lineages": n_all}