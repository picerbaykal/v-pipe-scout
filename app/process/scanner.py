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
            missing_from_panel:  list of {variant, total_reads, pattern_count}
            emerging_sublineage: list of {lineage, parent, total_reads, pattern_count}
            possibly_new:        {total_reads, pattern_count, top_patterns}
            total_unexplained_reads: int
            summary: human-readable one-line summary
    """
    # NOTE on "possibly new" coverage:
    # Currently checks patterns only at ~421 positions from known variant signatures.
    # When SaneQL deploys, expand to all covered positions (non-N reads > threshold)
    # for comprehensive novel variant detection. Zero-coverage positions should be
    # excluded from the query set — no coverage means no signal to detect.
    # See: worker/cooc.py TODO for the SaneQL migration plan.
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
    bucket1_explained = set()  # pattern indices explained by bucket 1

    for idx, row in patterns.iterrows():
        present = set(row["confirmed_present"])
        count = int(row["count"])
        for variant, sig in cowwid_not_in_panel.items():
            if _sig_explains(present, sig):
                if variant not in missing_hits:
                    missing_hits[variant] = {
                        "total_reads": 0,
                        "pattern_count": 0,
                        "observed_mutations": set(),  # mutations seen in unexplained patterns
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
            # only count if pattern contains at least one private mutation
            # (not just inherited from parent panel variant)
            parent_sig = all_lineage_signatures.get(panel_parent, set())
            private_sig = sig - parent_sig
            if not (present & private_sig):
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
            emerging_hits[lineage]["observed_mutations"].update(present & private_sig)

    emerging_sublineage = sorted(
        [{
            "lineage": l,
            "parent": s["parent"],
            "total_reads": s["total_reads"],
            "pattern_count": s["pattern_count"],
            "observed_mutations": sorted(s["observed_mutations"]),
        } for l, s in emerging_hits.items()],
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