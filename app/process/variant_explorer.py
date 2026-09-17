"""Variant explorer — on-demand investigation of ANY variant's co-occurrence
detectability and its connections to other variants.

Purpose (not detection — the scanner does that): let a user look up any pango
lineage, including blind spots the scanner cannot surface (e.g. KP.2), to
understand WHY it's a blind spot, HOW it connects to other variants, and its
co-occurrence signal if any. Pure logic (no Streamlit, no IO): everything is
derived from the pango tree + signatures, plus an optional co-occurrence result
for "found in data". Not anchored to the officially-tracked list — works for any
variant the user types.
"""

from typing import Dict, List, Optional, Set

# a mutation is "distinctive" (variant-discriminating) if carried by at most this
# many lineages — the same notion the scanner uses for clade-distinctiveness
_DISC_CARRIERS = 30


def _pos(m: str) -> str:
    return m[:-1] if m and m[-1].isalpha() else m


class _Tree:
    """Cheap parent/child index over the pango raw data."""
    def __init__(self, raw: dict):
        self.parent = {l: e.get("parent", "") for l, e in raw.items()}
        self.child: Dict[str, List[str]] = {}
        for l, p in self.parent.items():
            if p:
                self.child.setdefault(p, []).append(l)

    def ancestors(self, v: str) -> List[str]:
        out, cur, seen = [], v, set()
        while cur in self.parent and self.parent[cur] and cur not in seen:
            seen.add(cur)
            cur = self.parent[cur]
            out.append(cur)
        return out

    def children(self, v: str) -> List[str]:
        return sorted(self.child.get(v, []))

    def siblings(self, v: str) -> List[str]:
        p = self.parent.get(v, "")
        if not p:
            return []
        return sorted(c for c in self.child.get(p, []) if c != v)


def _mutation_carrier_counts(all_sigs: Dict[str, Set[str]]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for s in all_sigs.values():
        for m in s:
            counts[m] = counts.get(m, 0) + 1
    return counts


def investigate_variant(
    variant: str,
    pango_loader,
    panel: Optional[List[str]] = None,
    cooc_result: Optional[dict] = None,
    all_sigs: Optional[Dict[str, Set[str]]] = None,
    mut_carriers: Optional[Dict[str, int]] = None,
) -> Dict:
    """Return a fact sheet + connections for `variant`.

    Args:
        variant: the pango lineage to investigate (any lineage, not just OT/panel).
        pango_loader: provides get_raw_data() and get_signature().
        panel: current panel variant names (for the in-panel flag).
        cooc_result: optional cooc-pipeline result (for found-in-data).
        all_sigs / mut_carriers: optional precomputed indices (built if omitted).

    Returns dict:
        {found, name, in_data, in_panel, is_recombinant,
         detectable, n_discriminating, reason,
         found_in_data, found_reads,
         parent, siblings, children, oscillates_with, nearest_detectable}
    """
    raw = pango_loader.get_raw_data()
    if variant not in raw:
        return {"found": False, "name": variant,
                "reason": "Not a known pango lineage in the current data."}

    if all_sigs is None:
        all_sigs = {l: pango_loader.get_signature(l) for l in raw}
        all_sigs = {l: s for l, s in all_sigs.items() if s}
    if mut_carriers is None:
        mut_carriers = _mutation_carrier_counts(all_sigs)

    tree = _Tree(raw)
    sig = pango_loader.get_signature(variant) or set()
    panel = panel or []

    # detectability: distinctive (few-carrier) mutations
    disc = sorted(m for m in sig if mut_carriers.get(m, 0) <= _DISC_CARRIERS)
    detectable = len(disc) >= 1
    is_recomb = variant.startswith("X") and not raw.get(variant, {}).get("parent")

    parent = raw.get(variant, {}).get("parent", "")
    siblings = tree.siblings(variant)
    children = tree.children(variant)

    # reason string
    if detectable:
        reason = (f"{len(disc)} distinctive mutation(s) carried by ≤{_DISC_CARRIERS} "
                  f"lineages (e.g. {', '.join(disc[:3])}) — co-occurrence can "
                  f"resolve this clade.")
    else:
        # find the single lowest-carrier mutation to explain the sharing
        if sig:
            m_least = min(sig, key=lambda m: mut_carriers.get(m, 0))
            n_least = mut_carriers.get(m_least, 0)
            reason = (f"0 distinctive mutations — its rarest mutation ({m_least}) is "
                      f"shared with {n_least} lineages, so no combination is unique "
                      f"to it. Co-occurrence can't confirm it; quantify with "
                      f"deconvolution.")
        else:
            reason = "No signature mutations available for this lineage."

    # oscillating partners: siblings that are near-identical (Jaccard ≥ 0.98).
    # Whole families (e.g. KP.*) can be mutually near-identical, so this may be
    # several — that's the honest picture (they can't be told apart).
    oscillates_with = []
    for sib in siblings:
        ss = all_sigs.get(sib, set())
        if not ss or not sig:
            continue
        inter = len(sig & ss)
        union = len(sig | ss)
        if union and inter / union >= 0.98:
            oscillates_with.append(sib)

    # For a blind spot: anchor to the PARENT clade, and look within that parent's
    # branch (its descendants) for any GENUINELY detectable variant the user could
    # track instead. Honest and simple — no wide "nearest relative" hunt.
    def _is_detectable(cand: str) -> bool:
        cs = all_sigs.get(cand, set())
        return sum(1 for m in cs if mut_carriers.get(m, 0) <= _DISC_CARRIERS) >= 1

    def _branch(v: str) -> List[str]:
        out, stack = [], [v]
        seen = {v}
        while stack:
            x = stack.pop()
            for c in tree.children(x):
                if c not in seen:
                    seen.add(c); out.append(c); stack.append(c)
        return out

    _osc = set(oscillates_with)
    detectable_in_branch = None
    if not detectable and parent:
        # scan the parent's whole branch for a detectable member (excluding this
        # variant's own near-identical oscillating siblings)
        for cand in _branch(parent):
            if cand == variant or cand in _osc:
                continue
            if _is_detectable(cand):
                detectable_in_branch = cand
                break
    # keep the field name used downstream
    nearest_detectable = detectable_in_branch

    # found in data (from cooc result's unexplained + matched patterns if present)
    found_in_data, found_reads = False, 0
    if cooc_result and sig:
        # a pattern "belongs" to this variant if it's a subset of its signature
        for p in cooc_result.get("unexplained_patterns", []):
            pat = set(p.get("confirmed_present", []))
            if len(pat) >= 2 and pat.issubset(sig):
                found_in_data = True
                found_reads += int(p.get("count", 0))

    return {
        "found": True,
        "name": variant,
        "in_data": True,
        "in_panel": variant in set(panel),
        "is_recombinant": is_recomb,
        "detectable": detectable,
        "n_discriminating": len(disc),
        "discriminating": disc[:8],
        "reason": reason,
        "found_in_data": found_in_data,
        "found_reads": found_reads,
        "parent": parent,
        "siblings": siblings[:8],
        "children": children[:8],
        "oscillates_with": oscillates_with,
        "nearest_detectable": nearest_detectable,
    }