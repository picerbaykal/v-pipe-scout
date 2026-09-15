"""Co-occurrence annotation logic for the integrated deconvolution view.

Given a deconvolution panel + scanner findings, this module computes, for each
variant, whether co-occurrence corroborates it — turning the standalone scanner
into per-variant annotations on the deconvolution results.

All functions are pure (no Streamlit, no IO) so they're unit-testable and can be
called from the worker or the UI. Signatures come from PangoLoader.

Key concepts (validated against real Lugano data + COJAC):
  - detectability: a variant is co-occurrence-detectable if its tightest
    amplicon-local mutation combination is carried by few lineages OUTSIDE its
    own clade (`outside` count). Low outside -> distinctive haplotype exists.
  - Only distinctive-haplotype variants can be confirmed by co-occurrence.
    Variants defined by few/broadly-shared mutations (KP.2, KP.3) are blind
    spots — deconvolution's domain.
"""

import re
from functools import lru_cache
from typing import Dict, List, Optional, Set, Tuple

# outside-clade carrier threshold: a combo carried by <= this many lineages
# outside the variant's own clade is "distinctive" (co-occurrence can resolve
# the clade). Derived empirically — most variants sit far from this boundary.
DETECTABLE_OUTSIDE_MAX = 15

# jaccard >= this between two panel variants => they oscillate in deconvolution
# (near-collinear signatures; deconvolution can't split them reliably).
OSCILLATION_JACCARD = 0.98

_POS_RE = re.compile(r"^(\d+)")


def _pos(m: str) -> int:
    mm = _POS_RE.match(m)
    return int(mm.group(1)) if mm else -1


class VariantIndex:
    """Precomputed indices over all lineage signatures for fast annotation.

    Build once per request; reuse across all variants. Holds the mutation ->
    lineages index (for combo carrier lookups) and the parent/child maps (for
    clade membership).
    """

    def __init__(self, all_sigs: Dict[str, Set[str]],
                 parent_map: Dict[str, str]):
        self.all_sigs = all_sigs
        self.parent = parent_map
        self.child: Dict[str, List[str]] = {}
        for lin, par in parent_map.items():
            self.child.setdefault(par, []).append(lin)
        self.mut_index: Dict[str, Set[str]] = {}
        for lin, s in all_sigs.items():
            for m in s:
                self.mut_index.setdefault(m, set()).add(lin)
        self._clade_cache: Dict[str, Set[str]] = {}
        self._outside_cache: Dict[str, Optional[int]] = {}

    def carriers(self, combo) -> Set[str]:
        it = iter(combo)
        try:
            acc = set(self.mut_index.get(next(it), ()))
        except StopIteration:
            return set()
        for m in it:
            acc &= self.mut_index.get(m, set())
            if not acc:
                break
        return acc

    def clade_of(self, v: str) -> Set[str]:
        cached = self._clade_cache.get(v)
        if cached is not None:
            return cached
        out = {v}
        stack = [v]
        while stack:
            x = stack.pop()
            for c in self.child.get(x, []):
                if c not in out:
                    out.add(c)
                    stack.append(c)
        self._clade_cache[v] = out
        return out

    def outside(self, v: str) -> Optional[int]:
        """Fewest lineages OUTSIDE v's clade carrying any amplicon-local combo
        of v's signature. None if v has no signature. Lower = more detectable."""
        if v in self._outside_cache:
            return self._outside_cache[v]
        vsig = self.all_sigs.get(v)
        if not vsig:
            self._outside_cache[v] = None
            return None
        positions = sorted(_pos(m) for m in vsig if _pos(m) >= 0)
        clusters, cur = [], []
        for p in positions:
            if cur and p - cur[-1] > 350:
                clusters.append(cur)
                cur = []
            cur.append(p)
        if cur:
            clusters.append(cur)
        clade = self.clade_of(v)
        best = None
        for cl in clusters:
            pos_set = set(cl)
            g = frozenset(m for m in vsig if _pos(m) in pos_set)
            if len(g) < 2:
                continue
            n_out = sum(1 for l in self.carriers(g) if l not in clade)
            if best is None or n_out < best:
                best = n_out
        self._outside_cache[v] = best
        return best

    def detectable(self, v: str) -> bool:
        o = self.outside(v)
        return o is not None and o <= DETECTABLE_OUTSIDE_MAX

    def branch_activity(self, v: str, found: Set[str], max_steps: int = 4) -> int:
        """How many confirmed-present variants are CLOSE relatives of v — within
        max_steps up the tree to a common ancestor. Ancient shared ancestors
        (JN.1, B.1) don't count: sharing B.1 doesn't make XFG a relevant relative
        of Alpha. A dead branch (0) means v is unlikely to circulate; an active
        branch means v has potential even if undetectable."""
        v_anc = self._ancestors_limited(v, max_steps)
        v_clade = self.clade_of(v)
        n = 0
        for f in found:
            if f == v or f in v_clade or v in self.clade_of(f):
                n += 1
                continue
            # close common ancestor within max_steps on both sides
            f_anc = self._ancestors_limited(f, max_steps)
            if v_anc & f_anc:
                n += 1
        return n

    def _ancestors_limited(self, v: str, max_steps: int) -> Set[str]:
        out, cur, steps = set(), v, 0
        while cur and steps < max_steps:
            cur = self.parent.get(cur, "")
            if cur:
                out.add(cur)
            steps += 1
        return out

    def _ancestors(self, v: str) -> Set[str]:
        out, cur, seen = set(), v, set()
        while cur and cur not in seen:
            seen.add(cur)
            cur = self.parent.get(cur, "")
            if cur:
                out.add(cur)
        return out


def jaccard(a: Set[str], b: Set[str]) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


def oscillating_pairs(panel: List[str],
                      all_sigs: Dict[str, Set[str]]) -> Dict[str, List[str]]:
    """Panel variants that are near-identical (jaccard >= threshold) and will
    oscillate in deconvolution. Returns variant -> [its near-identical siblings]."""
    out: Dict[str, List[str]] = {}
    for i, a in enumerate(panel):
        for b in panel[i + 1:]:
            sa, sb = all_sigs.get(a, set()), all_sigs.get(b, set())
            if jaccard(sa, sb) >= OSCILLATION_JACCARD:
                out.setdefault(a, []).append(b)
                out.setdefault(b, []).append(a)
    return out


def annotate_variant(variant: str,
                     idx: VariantIndex,
                     scanner_found: Set[str],
                     osc_pairs: Dict[str, List[str]]) -> Tuple[str, str]:
    """Co-occurrence verdict for one deconvolution variant.

    Returns (status, reason) where status is one of:
      "confirmed"    — scanner confirmed it (or its clade) present
      "oscillating"  — near-identical panel sibling; deconvolution unstable
      "cant_confirm" — no distinctive haplotype, or detectable but not seen
    """
    # A variant is confirmed if a scanner finding COVERS it: the finding's node
    # equals the variant, or the variant descends from the found clade (the
    # finding's node is an ancestor of the variant). We do NOT confirm a variant
    # just because one of ITS descendants was found — finding XFG does not
    # confirm the ancestral JN.1.
    found = any(f == variant or variant in idx.clade_of(f)
                for f in scanner_found)
    if found:
        return "confirmed", "co-occurrence confirms"
    sib = osc_pairs.get(variant, [])
    if sib:
        return "oscillating", f"oscillates with {', '.join(sib)} — trust the sum"
    o = idx.outside(variant)
    if o is not None and o > DETECTABLE_OUTSIDE_MAX:
        return "cant_confirm", "no distinctive haplotype (blind spot)"
    # detectable, but not surfaced as a scanner finding. For a panel variant
    # this is expected (the scanner reports non-panel findings), so we don't
    # claim it's absent — only that co-occurrence CAN resolve this variant and
    # its level is deconvolution's estimate.
    return "detectable", "has a distinctive haplotype (co-occurrence can resolve it)"


def variant_facts(variant: str,
                  idx: VariantIndex,
                  panel: List[str],
                  scanner_found: Set[str]) -> dict:
    """Systematic fact sheet for the 'learn about a variant' explainer.
    Pure quantitative fields — no prose interpretation."""
    o = idx.outside(variant)
    found = any(f == variant or variant in idx.clade_of(f)
                for f in scanner_found)
    return {
        "variant": variant,
        "in_panel": variant in panel,
        "detectable": (o is not None and o <= DETECTABLE_OUTSIDE_MAX),
        "outside": o,
        "found": found,
        "branch_activity": idx.branch_activity(variant, scanner_found),
    }