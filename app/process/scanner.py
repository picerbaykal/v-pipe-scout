"""Panel scanner: classify unexplained co-occurrence patterns (Option C).

Given the unexplained patterns from run_cooc_panel_completeness and the full
pango tree, assigns each co-occurrence pattern to the tightest pango node its
fingerprint supports, then categorizes by relationship to the user's panel.

Categories:
  resolved_lineage  — fingerprint matches exactly one pango lineage.
  resolved_clade    — fingerprint matches several lineages forming a tight
                      clade; labelled by their common ancestor.
  unresolved        — fingerprint matches many lineages across unrelated
                      clades (common ancestor too ancient to be meaningful).
  novel             — no pango lineage explains the fingerprint.

Within resolved_* each finding is tagged by panel relationship:
  in_panel     — the assigned node is a panel variant (explained; dropped).
  sublineage   — assigned node descends from a panel variant.
  new_lineage  — assigned node unrelated to the panel.

Co-occurrence requires >= 2 fingerprint mutations by definition (a single
mutation is allele frequency, not a haplotype), so patterns with fewer than
2 mutations beyond the panel are excluded.

Called by the worker Celery task (run_cooc_scanner_lapis) — pure computation.
"""

import logging
import re
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

# A clade label is only meaningful if the common ancestor of the candidate
# set is reasonably recent. If the common ancestor is at or above this depth
# threshold from the root it's too broad (e.g. BA.2) -> unresolved.
# Depth is measured as number of parent hops from the node to the tree root.
MIN_CLADE_DEPTH = 6

# Minimum fingerprint size for co-occurrence (2 = haplotype, 1 = allele freq).
MIN_FINGERPRINT = 2


def _sig_explains(present: Set[str], sig: Set[str]) -> bool:
    return bool(present) and present.issubset(sig)


class _Tree:
    """Lightweight pango tree helper built from a parent map."""

    def __init__(self, parent_map: Dict[str, str]):
        self.parent = parent_map
        self.children: Dict[str, List[str]] = {}
        for node, par in parent_map.items():
            if par:
                self.children.setdefault(par, []).append(node)
        self._depth: Dict[str, int] = {}

    def depth(self, node: str) -> int:
        if node in self._depth:
            return self._depth[node]
        d, cur = 0, node
        seen = set()
        while True:
            par = self.parent.get(cur, "")
            if not par or par in seen:
                break
            seen.add(par)
            d += 1
            cur = par
        self._depth[node] = d
        return d

    def ancestors(self, node: str) -> Set[str]:
        a, cur, seen = set(), node, set()
        while cur and cur not in seen:
            a.add(cur)
            seen.add(cur)
            cur = self.parent.get(cur, "")
        return a

    def is_descendant(self, node: str, ancestor: str) -> bool:
        cur, seen = node, set()
        while cur and cur not in seen:
            seen.add(cur)
            cur = self.parent.get(cur, "")
            if cur == ancestor:
                return True
        return False

    def lca(self, nodes: List[str]) -> Optional[str]:
        """Deepest common ancestor of all nodes (may be one of the nodes)."""
        if not nodes:
            return None
        common: Optional[Set[str]] = None
        for n in nodes:
            a = self.ancestors(n)
            common = a if common is None else (common & a)
        if not common:
            return None
        return max(common, key=self.depth)

    def dominant_clade(
        self, nodes: List[str], min_fraction: float = 0.6
    ) -> Optional[str]:
        """Deepest node that is an ancestor of >= min_fraction of `nodes`.

        Unlike strict LCA, this tolerates outliers — recombinants (which have
        no parent chain) and convergent lineages that share the fingerprint
        but sit outside the main clade don't drag the label up to the root.
        Returns the tightest (deepest) clade covering the bulk of candidates.
        """
        if not nodes:
            return None
        # count how many candidates each ancestor covers
        cover: Dict[str, int] = {}
        for n in nodes:
            for anc in self.ancestors(n):
                cover[anc] = cover.get(anc, 0) + 1
        threshold = max(2, int(len(nodes) * min_fraction))
        qualifying = [a for a, c in cover.items() if c >= threshold]
        if not qualifying:
            return None
        return max(qualifying, key=self.depth)


def _assign(
    fingerprint: Set[str],
    all_sigs: Dict[str, Set[str]],
    tree: _Tree,
) -> Tuple[Optional[str], str, List[str]]:
    """Assign a fingerprint to the tightest clade it supports.

    Returns (label_node, kind, candidates).
      kind: 'clade' | 'unresolved' | 'novel'
      label_node: the clade root (clade), else None.
      candidates: all lineages whose signature contains the fingerprint.

    We deliberately do NOT claim a single "exact" lineage. A fingerprint
    matching exactly one lineage is usually a coincidence of which lineages
    happen to carry those muts (e.g. {1722T,1895A} intersecting at PY.1.1),
    not evidence that lineage specifically is present. Co-occurrence resolves
    to clades; which member drives the signal is shown in the drill-down
    discriminating-mutation heatmap, not claimed as a label.
    """
    candidates = [l for l, s in all_sigs.items() if fingerprint.issubset(s)]
    if not candidates:
        return None, "novel", []
    if len(candidates) == 1:
        # single candidate: label it as a clade-of-one at that node, but only
        # if it's deep enough to be meaningful; else unresolved.
        node = candidates[0]
        if tree.depth(node) >= MIN_CLADE_DEPTH:
            return node, "clade", candidates
        return None, "unresolved", candidates

    clade = tree.dominant_clade(candidates)
    if clade is None or tree.depth(clade) < MIN_CLADE_DEPTH:
        return None, "unresolved", candidates
    return clade, "clade", candidates


def _panel_relationship(
    node: str,
    panel_set: Set[str],
    tree: _Tree,
) -> Tuple[str, Optional[str]]:
    """Return (relationship, panel_ancestor).

    relationship: 'in_panel' | 'sublineage' | 'new_lineage'.
    """
    if node in panel_set:
        return "in_panel", node
    for pv in panel_set:
        if tree.is_descendant(node, pv):
            return "sublineage", pv
    return "new_lineage", None


def scan_unexplained_patterns(
    unexplained_patterns: pd.DataFrame,
    panel_variants: List[str],
    all_lineage_signatures: Dict[str, Set[str]],
    panel_parent_map: Dict[str, str],
    min_read_count: int = 2,
    # kept for backwards-compat with the task signature; unused in Option C.
    cowwid_signatures: Optional[Dict[str, Set[str]]] = None,
    truly_private_muts: Optional[Dict[str, Set[str]]] = None,
) -> dict:
    """Classify unexplained co-occurrence patterns (Option C).

    Args:
        unexplained_patterns: DataFrame with columns date, count,
            confirmed_present (list of "{pos}{alt}" per row).
        panel_variants: currently selected panel variant names.
        all_lineage_signatures: {lineage: set of "{pos}{alt}"} for all pango.
        panel_parent_map: {lineage: parent} for the full pango tree.
        min_read_count: minimum reads for a pattern to count.

    Returns dict with keys:
        resolved_lineage:  [{node, relationship, panel_ancestor, total_reads,
                             pattern_count, observed_mutations, designation}]
        resolved_clade:    [{node, relationship, panel_ancestor, member_count,
                             members, total_reads, pattern_count,
                             observed_mutations, designation}]
        unresolved:        [{fingerprint, candidate_count, common_ancestor,
                             total_reads, pattern_count}]
        novel:             {total_reads, pattern_count, top_patterns}
        total_unexplained_reads: int
        summary: str
    """
    panel_set = set(panel_variants)
    tree = _Tree(panel_parent_map)

    panel_union: Set[str] = set()
    for pv in panel_set:
        panel_union |= all_lineage_signatures.get(pv, set())

    patterns = unexplained_patterns
    if patterns is None or patterns.empty:
        return _empty_result()

    # aggregate per assigned clade
    clade_hits: Dict[str, dict] = {}
    unresolved_hits: Dict[frozenset, dict] = {}
    novel_reads = 0
    novel_patterns: List[dict] = []
    total_unexplained = 0

    # collect all observed co-occurrence patterns (mut-sets) for member
    # filtering — a member is only plotted if its specific amplicon combo
    # actually appears co-occurring in the data. Stored WITH read counts so we
    # can report per-region co-occurrence strength for the UI threshold slider.
    observed_patterns: List = []   # list of (frozenset(muts), count)

    for _, row in patterns.iterrows():
        present = set(row["confirmed_present"])
        count = int(row["count"])
        if count < min_read_count:
            continue
        total_unexplained += count
        if len(present) >= 2:
            observed_patterns.append((frozenset(present), count))

        fingerprint = present - panel_union
        if len(fingerprint) < MIN_FINGERPRINT:
            continue  # not co-occurrence beyond panel

        node, kind, candidates = _assign(
            fingerprint, all_lineage_signatures, tree
        )

        if kind == "novel":
            novel_reads += count
            if len(novel_patterns) < 10:
                novel_patterns.append(
                    {"count": count, "date": row.get("date", ""),
                     "mutations": sorted(fingerprint)[:8]}
                )
            continue

        if kind == "unresolved":
            lca = tree.dominant_clade(candidates, min_fraction=0.9) or ""
            key = frozenset(fingerprint)
            slot = unresolved_hits.get(key)
            if slot is None:
                slot = unresolved_hits[key] = {
                    "fingerprint": sorted(fingerprint),
                    "candidate_count": len(candidates),
                    "common_ancestor": lca or "",
                    "total_reads": 0, "pattern_count": 0,
                }
            slot["total_reads"] += count
            slot["pattern_count"] += 1
            continue

        # clade (may still be in_panel -> explained, drop those)
        rel, panel_anc = _panel_relationship(node, panel_set, tree)
        if rel == "in_panel":
            continue

        slot = clade_hits.get(node)
        if slot is None:
            slot = clade_hits[node] = {
                "node": node,
                "relationship": rel,
                "panel_ancestor": panel_anc,
                "total_reads": 0, "pattern_count": 0,
                "observed_mutations": set(),
                "designation": "",
                "candidates": set(),
            }
        slot["total_reads"] += count
        slot["pattern_count"] += 1
        slot["observed_mutations"].update(fingerprint)
        slot["candidates"].update(candidates)

    # ── build output lists ────────────────────────────────────────────────
    _all_clades = [
        _finalize_clade(s, tree, all_lineage_signatures, observed_patterns)
        for s in clade_hits.values()
    ]

    # nest clades that are inside another reported clade. e.g. PY.1.1 sits
    # inside LF.7 — show it as a sub-finding of LF.7, not a separate top-level
    # entry, so the same branch isn't reported multiple times.
    _by_node = {c["node"]: c for c in _all_clades}
    _nodes = set(_by_node)
    for c in _all_clades:
        c["sub_findings"] = []
        c["parent_clade"] = None
        cur = tree.parent.get(c["node"], "")
        seen = set()
        while cur and cur not in seen:
            seen.add(cur)
            if cur in _nodes:
                c["parent_clade"] = cur
                break
            cur = tree.parent.get(cur, "")

    top_level = []
    for c in _all_clades:
        if c["parent_clade"]:
            _by_node[c["parent_clade"]]["sub_findings"].append(c)
        else:
            top_level.append(c)

    # flag whether each clade resolved to a specific member (has blocks) or
    # is only the clade backbone (panel-worthy at clade level, member unknown)
    for c in _all_clades:
        c["member_resolved"] = len(c.get("member_blocks", [])) > 0

    # A clade is a confirmed finding only if it has at least one discriminating
    # co-occurrence block (a specific + observed amplicon group). Clades that
    # were assigned only by a fingerprint match but have no discriminating
    # block (e.g. KW.1.2, PA.1 — flagged by a couple of broad mutations) are
    # NOT co-occurrence-confirmed. Drop them from the confident findings and
    # record them as unresolved so their reads aren't silently lost.
    confirmed = []
    for c in top_level:
        # a clade counts if it, OR any of its sub-findings, has a block
        has_block = bool(c.get("member_blocks")) or any(
            sf.get("member_blocks") for sf in c.get("sub_findings", []))
        if has_block:
            confirmed.append(c)
        else:
            unresolved_hits[frozenset(c["observed_mutations"])] = {
                "fingerprint": c["observed_mutations"][:6],
                "candidate_count": c.get("member_count", 1),
                "common_ancestor": c["node"],
                "total_reads": c["total_reads"],
                "pattern_count": c.get("pattern_count", 0),
            }

    resolved_clade = sorted(confirmed, key=lambda x: -x["total_reads"])
    unresolved = sorted(
        unresolved_hits.values(), key=lambda x: -x["total_reads"]
    )

    novel = {
        "total_reads": novel_reads,
        "top_patterns": sorted(
            novel_patterns, key=lambda x: -x["count"]
        )[:10],
    }
    novel["pattern_count"] = _count_novel(
        patterns, panel_union, all_lineage_signatures, tree, min_read_count
    )

    summary = _summary(resolved_clade, unresolved, novel)

    result = {
        "resolved_clade": resolved_clade,
        "unresolved": unresolved,
        "novel": novel,
        "total_unexplained_reads": total_unexplained,
        "summary": summary,
    }
    # ── backward-compat shim ──────────────────────────────────────────────
    # Map the new clade-only shape onto the legacy keys the current UI reads,
    # so nothing crashes during the UI transition. Legacy consumers see:
    #   missing_from_panel  <- new-lineage clades (not descended from panel)
    #   emerging_sublineage <- sublineage clades (descend from panel)
    #   possibly_new        <- novel
    result["missing_from_panel"] = [
        {"variant": c["node"], "total_reads": c["total_reads"],
         "pattern_count": c["pattern_count"],
         "observed_mutations": c["observed_mutations"],
         "cluster_key": c["node"]}
        for c in resolved_clade if c["relationship"] == "new_lineage"
    ]
    result["emerging_sublineage"] = [
        {"lineage": c["node"], "parent": c["panel_ancestor"] or "",
         "total_reads": c["total_reads"], "pattern_count": c["pattern_count"],
         "observed_mutations": c["observed_mutations"]}
        for c in resolved_clade if c["relationship"] == "sublineage"
    ]
    result["possibly_new"] = novel
    return result


def _clade_root_of(members: List[str], tree: "_Tree") -> str:
    """Tightest single label for a set of members: the deepest node that is an
    ancestor of (or equal to) all of them. Falls back to the shortest name
    when they don't share a common ancestor (e.g. recombinants)."""
    if len(members) == 1:
        return members[0]
    common = None
    for m in members:
        anc = tree.ancestors(m) | {m}
        common = anc if common is None else (common & anc)
    if not common:
        return sorted(members, key=lambda x: (len(x), x))[0]
    return max(common, key=tree.depth)


def _finalize_clade(s: dict, tree: "_Tree", all_sigs: Dict[str, Set[str]],
                    observed_patterns: List = None) -> dict:
    """Build the clade finding with per-member co-occurrence blocks.

    Selection rule for which members to plot:
      1. specific  — the member has an amplicon-local mutation combination
         carried by few other lineages (discriminating).
      2. observed  — that specific combination actually appears co-occurring
         in the data (is a subset of some observed co-occurrence pattern).
      3. collapse  — members sharing the same observed combination can't be
         told apart, so they're grouped into one "family" block.

    This reduces a clade of 100+ candidates to the handful of distinguishable
    groups the co-occurrence data can actually confirm.
    """
    observed_patterns = observed_patterns or []
    node = s["node"]
    candidates = sorted(s["candidates"])
    phylo = [c for c in candidates
             if c == node or tree.is_descendant(c, node)]
    associated = [c for c in candidates
                  if c not in phylo and not tree.parent.get(c, "")]
    members = phylo + associated
    if not members:
        members = candidates

    member_sigs = {m: all_sigs.get(m, set()) for m in members}
    shared = set.intersection(*member_sigs.values()) if member_sigs else set()

    def _pos(m):
        mm = re.match(r"^(\d+)", m)
        return int(mm.group(1)) if mm else -1

    def _amplicon_groups(member_sig):
        """Amplicon-local mutation groups (positions within ~350bp)."""
        positions = sorted(_pos(m) for m in member_sig if _pos(m) >= 0)
        clusters, cur = [], []
        for p in positions:
            if cur and p - cur[-1] > 350:
                clusters.append(cur)
                cur = []
            cur.append(p)
        if cur:
            clusters.append(cur)
        out = []
        for cl in clusters:
            pos_set = set(cl)
            g = frozenset(m for m in member_sig if _pos(m) in pos_set)
            if len(g) >= 2:
                out.append(g)
        return out

    def _is_observed(group):
        # fractional: the group is "observed" if a single read/pattern shows
        # most of it (>= 80%), not necessarily all. Dense discriminating groups
        # (e.g. NB.1.8.1's 24-mut spike) rarely appear in full on one pattern
        # due to coverage/variation, but a strong partial match is real signal.
        g = set(group)
        best = 0.0
        for op, _cnt in observed_patterns:
            inter = len(g & op)
            if inter >= 2:
                best = max(best, inter / len(g))
        return best >= 0.8

    def _group_reads(group):
        # total co-occurrence reads where this group appears (fractional >=80%)
        # — the region's signal strength, used for the UI threshold slider.
        g = set(group)
        total = 0
        for op, cnt in observed_patterns:
            inter = len(g & op)
            if inter >= 2 and inter / len(g) >= 0.8:
                total += cnt
        return total

    # for each member, find its best specific+observed amplicon group(s)
    # keyed by the observed combination so identical combos collapse.
    # Specificity counts carriers OUTSIDE the clade family — sibling
    # descendants (PQ.* for NB.1.8.1) sharing the group are the same family,
    # not "other" variants, so they don't count against specificity.
    family = set(members)
    combo_to_members = {}   # frozenset(group) -> {members, n_other, n_total}
    for m in members:
        for g in _amplicon_groups(member_sigs[m]):
            carriers = [l for l, sg in all_sigs.items() if g.issubset(sg)]
            n_outside = sum(1 for l in carriers if l not in family)
            if n_outside > 15:
                continue          # not specific (many non-family carriers)
            if not _is_observed(g):
                continue          # not seen co-occurring
            slot = combo_to_members.setdefault(
                g, {"members": [], "n_other": n_outside,
                    "n_total": len(carriers)})
            slot["members"].append(m)

    # build blocks: one per distinguishable observed amplicon group. Each is a
    # piece of discriminating co-occurrence evidence for this clade — no region
    # is privileged, all observed discriminating groups are shown. Groups are
    # labelled by the clade family + a region index (the amplicon they're in).
    # For each block we also compute ABSENT markers: mutations from competing
    # variants in the SAME amplicon window that the family lacks (present +
    # absent gives sharper discrimination than present-only).
    _raw_blocks = []
    for g, info in sorted(combo_to_members.items(),
                          key=lambda kv: (min(_pos(m) for m in kv[0]))):
        fam = sorted(info["members"], key=lambda x: (len(x), x))
        gpos = sorted(_pos(m) for m in g)
        lo, hi = gpos[0] - 20, gpos[-1] + 20
        fam_sig = set().union(*(all_sigs.get(m, set()) for m in fam))
        absent = set()
        for l, sl in all_sigs.items():
            if l in fam:
                continue
            for mm in sl:
                p = _pos(mm)
                if lo <= p <= hi and mm not in fam_sig:
                    n_car = sum(1 for x in all_sigs.values() if mm in x)
                    if n_car <= 60:
                        absent.add(mm)
        _raw_blocks.append({
            "members": fam,
            "n_outside": info["n_other"],
            "n_total": info["n_total"],
            "discriminating": sorted(g, key=_pos)[:25],
            "absent_markers": sorted(absent, key=_pos)[:6],
            "region_start": gpos[0],
            "reads": _group_reads(g),
        })

    # label: name each region by the tightest family its discriminating group
    # points to (the lineages carrying that exact group), not a generic index.
    # This surfaces the actual driver — XFG, XFP, PY.1.1 etc. — per region.
    blocks = []
    for b in _raw_blocks:
        fam = b["members"]
        ntot = b["n_total"]   # total lineages carrying this group (specificity)
        # label each region by the tightest family it points to, with the
        # lineage count so the user can judge specificity directly:
        #   1 carrier          -> "XFG @ 4184"              (specific variant)
        #   small (<=15)       -> "XFG family @ 4184 [9 lineages]"
        #   broad (>15)        -> "NB.1.8.1 clade @ 8299 [68 lineages]"
        root = fam[0] if len(fam) == 1 else _clade_root_of(fam, tree)
        if ntot == 1:
            label = f"{root} @ {b['region_start']}"
        elif ntot <= 15:
            label = f"{root} family @ {b['region_start']} [{ntot} lineages]"
        else:
            label = f"{root} clade @ {b['region_start']} [{ntot} lineages]"
        blocks.append({
            "member": label,
            "family_root": root,
            "members": fam[:10],
            "member_count": len(fam),
            "discriminating": b["discriminating"],
            "absent_markers": b["absent_markers"],
            "reads": b["reads"],
            # per-mutation carrier count so the UI can show which rows are
            # discriminating (few carriers) vs backbone (many carriers)
            "mut_carriers": {m: sum(1 for s in all_sigs.values() if m in s)
                             for m in b["discriminating"]},
        })
    # strongest regions first (by co-occurrence reads) for the UI slider
    blocks.sort(key=lambda x: -x["reads"])

    # plottable members = union of all collapsed families
    plottable = sorted({m for info in combo_to_members.values()
                        for m in info["members"]})

    # no separate "shared/backbone" block — the discriminating groups ARE the
    # evidence. Keep shared_mutations empty so the UI doesn't show backbone.
    clade_block_muts = []

    return {
        "node": node,
        "relationship": s["relationship"],
        "panel_ancestor": s["panel_ancestor"],
        "member_count": len(members),
        "members": (associated + phylo)[:30],
        "associated_members": associated[:10],
        "plottable_members": plottable,
        "total_reads": s["total_reads"],
        "pattern_count": s["pattern_count"],
        "observed_mutations": sorted(s["observed_mutations"]),
        "shared_mutations": clade_block_muts,
        "member_blocks": blocks,
        "designation": s["designation"],
    }


def _count_novel(patterns, panel_union, all_sigs, tree, min_read_count) -> int:
    n = 0
    all_sig_list = list(all_sigs.values())
    for _, row in patterns.iterrows():
        count = int(row["count"])
        if count < min_read_count:
            continue
        fp = set(row["confirmed_present"]) - panel_union
        if len(fp) < MIN_FINGERPRINT:
            continue
        if not any(fp.issubset(s) for s in all_sig_list):
            n += 1
    return n


def _summary(clade, unresolved, novel) -> str:
    parts = []
    if clade:
        top = clade[0]
        label = f"{top['node']} clade" if top["member_count"] > 1 else top["node"]
        extra = len(clade) - 1
        parts.append(label + (f" + {extra} more" if extra else ""))
    if novel["total_reads"] > 0:
        parts.append(f"{novel['pattern_count']} novel pattern(s)")
    if unresolved:
        parts.append(f"{len(unresolved)} unresolved")
    if not parts:
        return "No co-occurrence signal beyond the panel."
    return "; ".join(parts) + " not explained by panel."


def _empty_result() -> dict:
    return {
        "resolved_clade": [], "unresolved": [],
        "novel": {"total_reads": 0, "pattern_count": 0, "top_patterns": []},
        "total_unexplained_reads": 0,
        "summary": "No unexplained patterns.",
    }