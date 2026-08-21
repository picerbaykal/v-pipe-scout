"""Panel scanner: classify unexplained co-occurrence patterns.

Given the unexplained patterns from run_cooc_panel_completeness and
a set of known lineage signatures, classifies each pattern into three buckets:

  missing_from_panel  — a cowwid surveillance variant not in the panel
                        explains this pattern. Actionable: add it.

  emerging_sublineage — a pango lineage (descendant of a panel variant,
                        not in cowwid) explains this pattern. Worth watching.

  possibly_new        — no known lineage explains this pattern. Could be
                        a novel variant, recombinant, or sequencing artifact.

Called by the worker Celery task (run_cooc_scanner_lapis) — no Streamlit,
no I/O. Pure computation on DataFrames and signature sets.
"""

import logging
import re
from typing import Dict, List, Set, Tuple

import pandas as pd

logger = logging.getLogger(__name__)


def _sig_explains(present: Set[str], sig: Set[str]) -> bool:
    """True if the signature explains the pattern — present is a subset of sig."""
    return bool(present) and present.issubset(sig)


def _is_descendant_of_panel(
    lineage: str,
    panel_set: Set[str],
    parent_map: Dict[str, str],
    max_depth: int = 10,
) -> Tuple[bool, str]:
    """
    Walk the parent chain up to max_depth steps.
    Returns (True, panel_ancestor) if a panel variant is found, else (False, "").
    """
    current = lineage
    for _ in range(max_depth):
        parent = parent_map.get(current, "")
        if not parent:
            break
        if parent in panel_set:
            return True, parent
        current = parent
    return False, ""


def _cluster_missing_ot(
    missing_from_panel: List[dict],
    candidate_names,
    parent_map: Dict[str, str],
) -> List[dict]:
    """Tag each missing OT variant with a stable cluster_key = its clade.

    Siblings only: two missing variants group iff they share the SAME immediate
    parent on the pango tree, and that shared parent is not itself a candidate
    (so a broad ancestor never absorbs its own children). This keeps the clade
    one step up from real lineages (e.g. XBB.1, never a basal catch-all like
    BA.2), which is the only level where "add the clade" is a lineage worth
    tracking in deconvolution. Cousins and ancestor/descendant pairs do NOT
    group — their similarity, if any, is the Jaccard heatmap's job, not this.

    Grouping is computed over the FULL candidate set (all cowwid-not-in-panel
    variants), not just the ones that got hits here, so cluster_key is identical
    across every city in a run (all cities share the panel -> same candidates).
    The UI aggregates findings by cluster_key; a group of >= 2 renders as one
    clade row headlined by the parent, a group of 1 renders as that lineage.

    Mutates and returns missing_from_panel, adding 'cluster_key' to each item.
    A variant with no qualifying sibling keys to its own name.
    """
    cand_set = set(candidate_names)

    by_parent: Dict[str, List[str]] = {}
    for v in cand_set:
        p = parent_map.get(v)
        if p:
            by_parent.setdefault(p, []).append(v)

    key_of: Dict[str, str] = {}
    for p, siblings in by_parent.items():
        # a real sibling group: >= 2 candidates sharing a parent that is not
        # itself a candidate finding
        if len(siblings) >= 2 and p not in cand_set:
            for v in siblings:
                key_of[v] = p

    for item in missing_from_panel:
        item["cluster_key"] = key_of.get(item["variant"], item["variant"])
    return missing_from_panel


def scan_unexplained_patterns(
    unexplained_patterns: pd.DataFrame,
    panel_variants: List[str],
    cowwid_signatures: Dict[str, Set[str]],
    all_lineage_signatures: Dict[str, Set[str]],
    panel_parent_map: Dict[str, str],
    min_read_count: int = 2,
) -> dict:
    """
    Classify unexplained co-occurrence patterns into three buckets.

    Args:
        unexplained_patterns: DataFrame with columns date, count, confirmed_present.
            confirmed_present is a list of "{pos}{alt}" strings per row.
            Produced by run_cooc_panel_completeness result["unexplained_patterns"].
        panel_variants: Currently selected panel variant names.
        cowwid_signatures: All cowwid surveillance variant signatures.
            Format: {variant_name: set of "{pos}{alt}" strings}
        all_lineage_signatures: All known pango lineage signatures.
            Format: {lineage_name: set of "{pos}{alt}" strings}
        panel_parent_map: lineage -> parent lineage, for descendant check.
        min_read_count: Minimum read count for a pattern to be considered.

    Returns:
        Dict with keys:
            missing_from_panel:  list of {variant, total_reads, pattern_count,
                                 observed_mutations, cluster_key}
                                 cluster_key = the clade (shared immediate
                                 parent) for sibling groups, else the variant's
                                 own name. Identical across every location in a
                                 run so the UI can aggregate then collapse each
                                 2+-member clade to one row.
            emerging_sublineage: list of {lineage, parent, total_reads, pattern_count}
            possibly_new:        {total_reads, pattern_count, top_patterns}
            total_unexplained_reads: int
            summary: human-readable one-line summary
    """
    if unexplained_patterns.empty:
        return _empty_result("No unexplained patterns to classify.")

    panel_set = set(panel_variants)
    cowwid_set = set(cowwid_signatures.keys())
    cowwid_not_in_panel = {
        v: sig for v, sig in cowwid_signatures.items()
        if v not in panel_set and sig
    }

    total_unexplained = int(unexplained_patterns["count"].sum())
    patterns = unexplained_patterns[
        unexplained_patterns["count"] >= min_read_count
    ].copy()

    if patterns.empty:
        return {
            **_empty_result(
                f"All {total_unexplained:,} unexplained reads are in singleton "
                f"patterns (< {min_read_count} reads each)."
            ),
            "total_unexplained_reads": total_unexplained,
        }

    # ── Bucket 1: missing cowwid tracked variants ─────────────────────────────
    missing_hits: Dict[str, dict] = {}

    for idx, row in patterns.iterrows():
        present = set(row["confirmed_present"])
        count = int(row["count"])
        for variant, sig in cowwid_not_in_panel.items():
            if _sig_explains(present, sig):
                if variant not in missing_hits:
                    missing_hits[variant] = {
                        "total_reads": 0,
                        "pattern_count": 0,
                        "observed_mutations": set(),
                    }
                missing_hits[variant]["total_reads"] += count
                missing_hits[variant]["pattern_count"] += 1
                missing_hits[variant]["observed_mutations"].update(present & sig)

    missing_from_panel = sorted(
        [{
            "variant": v,
            "total_reads": s["total_reads"],
            "pattern_count": s["pattern_count"],
            "observed_mutations": sorted(s["observed_mutations"]),
        } for v, s in missing_hits.items()],
        key=lambda x: -x["total_reads"],
    )

    # group missing variants by clade (siblings sharing an immediate parent):
    # tag each with a cluster_key stable across all cities. The UI aggregates
    # then collapses each 2+-member clade to a single row headlined by the
    # parent (adding the clade, not an unresolvable leaf).
    missing_from_panel = _cluster_missing_ot(
        missing_from_panel, cowwid_not_in_panel.keys(), panel_parent_map
    )

    # ── Bucket 2: emerging sublineages ───────────────────────────────────────
    emerging_hits: Dict[str, dict] = {}

    for idx, row in patterns.iterrows():
        present = set(row["confirmed_present"])
        count = int(row["count"])
        for lineage, sig in all_lineage_signatures.items():
            if lineage in cowwid_set or lineage in panel_set:
                continue
            if not _sig_explains(present, sig):
                continue
            is_desc, panel_parent = _is_descendant_of_panel(
                lineage, panel_set, panel_parent_map
            )
            if not is_desc:
                continue
            parent_sig = all_lineage_signatures.get(panel_parent, set())
            private_sig = sig - parent_sig
            private_observed = present & private_sig
            # require >=2 co-occurring private mutations — matches the
            # process/cooc.py uninformative cutoff and eliminates single-
            # position sublineages that cannot be distinguished from noise
            # (e.g. a sublineage whose only private mutation is one position
            # appearing at high frequency — indistinguishable from drift)
            if len(private_observed) < 2:
                continue
            if lineage not in emerging_hits:
                emerging_hits[lineage] = {
                    "parent": panel_parent,
                    "total_reads": 0,
                    "pattern_count": 0,
                    "observed_mutations": set(),
                }
            emerging_hits[lineage]["total_reads"] += count
            emerging_hits[lineage]["pattern_count"] += 1
            emerging_hits[lineage]["observed_mutations"].update(private_observed)

    # deduplicate by observed-mutation fingerprint:
    # sublineages matching the exact same private mutations are explaining
    # the same reads (signature overlap). Keep the most specific lineage —
    # deepest pango designation (most dots = most specific).
    fingerprint_to_lineages: dict = {}
    for l, s in emerging_hits.items():
        fp = frozenset(s["observed_mutations"])
        fingerprint_to_lineages.setdefault(fp, []).append(l)

    deduped: dict = {}
    for fp, lineages in fingerprint_to_lineages.items():
        winner = max(lineages, key=lambda x: (len(x.split(".")), x))
        deduped[winner] = emerging_hits[winner]

    emerging_sublineage = sorted(
        [{
            "lineage": l,
            "parent": s["parent"],
            "total_reads": s["total_reads"],
            "pattern_count": s["pattern_count"],
            "observed_mutations": sorted(s["observed_mutations"]),
        } for l, s in deduped.items()],
        key=lambda x: -x["total_reads"],
    )

    # ── Bucket 3: possibly new ───────────────────────────────────────────────
    possibly_new_reads = 0
    possibly_new_patterns = []
    all_lineage_sigs = list(all_lineage_signatures.values())

    for idx, row in patterns.iterrows():
        present = set(row["confirmed_present"])
        count = int(row["count"])
        if any(_sig_explains(present, sig) for sig in all_lineage_sigs):
            continue
        possibly_new_reads += count
        possibly_new_patterns.append({
            "mutations": sorted(present),
            "count": count,
            "date": str(row["date"]),
        })

    possibly_new_patterns.sort(key=lambda x: -x["count"])

    # ── Summary ───────────────────────────────────────────────────────────────
    parts = []
    if missing_from_panel:
        n = len(missing_from_panel)
        top = missing_from_panel[0]["variant"]
        parts.append(
            f"{top} + {n-1} other tracked variant(s) missing from panel"
            if n > 1 else f"{top} missing from panel"
        )
    if emerging_sublineage:
        n = len(emerging_sublineage)
        top = emerging_sublineage[0]["lineage"]
        parts.append(
            f"{top} + {n-1} other sublineage(s) rising"
            if n > 1 else f"{top} sublineage rising"
        )
    if possibly_new_reads > 0:
        parts.append(f"{possibly_new_reads:,} reads match no known lineage")
    summary = "; ".join(parts) if parts else "All unexplained patterns are noise."

    return {
        "missing_from_panel": missing_from_panel,
        "emerging_sublineage": emerging_sublineage,
        "possibly_new": {
            "total_reads": possibly_new_reads,
            "pattern_count": len(possibly_new_patterns),
            "top_patterns": possibly_new_patterns[:10],
        },
        "total_unexplained_reads": total_unexplained,
        "summary": summary,
    }


def _empty_result(summary: str) -> dict:
    return {
        "missing_from_panel": [],
        "emerging_sublineage": [],
        "possibly_new": {"total_reads": 0, "pattern_count": 0, "top_patterns": []},
        "total_unexplained_reads": 0,
        "summary": summary,
    }