"""Amplicon coordinates from BED files.

Used for amplicon-aware co-occurrence queries. Positions on the same
amplicon can be queried together because they come from the same physical
read; positions on different amplicons cannot.
"""

import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)



def load_amplicons(bed_path: Path) -> List[Tuple[str, int, int]]:
    """
    Load amplicon coordinates from a BED file.

    Args:
        bed_path: Path to a BED file with amplicon inserts.
            Expected columns: chrom, start, end, name, ...
            Coordinates are 0-based half-open (standard BED).

    Returns:
        List of (name, start, end) tuples, coordinates converted to
        1-based inclusive to match VCF/LAPIS conventions.
    """
    amplicons = []
    with open(bed_path) as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                logger.warning(f"Skipping malformed BED line {i}: {line}")
                continue
            _chrom, start_bed, end_bed, name = parts[:4]
            # Convert BED 0-based half-open to 1-based inclusive
            start = int(start_bed) + 1
            end = int(end_bed)
            amplicons.append((name, start, end))
    logger.info(f"Loaded {len(amplicons)} amplicons from {bed_path}")
    return amplicons


def group_positions_by_amplicon(
    positions: set,
    amplicons: List[Tuple[str, int, int]],
) -> dict:
    """
    Group target positions by which amplicon they fall on.

    Positions in overlapping regions get assigned to ALL amplicons
    containing them (they can be queried on any of those amplicons).

    Args:
        positions: Set of 1-based positions of interest.
        amplicons: List of (name, start, end) tuples from load_amplicons.

    Returns:
        Dict mapping amplicon_name -> list of positions on that amplicon.
        Positions not covered by any amplicon are excluded (with a warning).
    """
    grouped: dict = {name: [] for name, _, _ in amplicons}
    uncovered = set(positions)
    assigned = set()  # track positions already assigned to an amplicon
    for name, start, end in amplicons:
        for pos in positions:
            if start <= pos <= end and pos not in assigned:
                grouped[name].append(pos)
                uncovered.discard(pos)
                assigned.add(pos)
    if uncovered:
        logger.warning(
            f"{len(uncovered)} positions not covered by any amplicon: {sorted(uncovered)[:10]}..."
        )
    # Drop amplicons with no positions of interest
    grouped = {name: sorted(pos) for name, pos in grouped.items() if pos}
    return grouped

def split_positions_by_distance(
    amplicon_groups: Dict[str, List[int]],
    max_distance: int = None,
) -> List[Tuple[str, List[int]]]:
    """
    Split each amplicon's positions into batches that fit within read length.

    Two positions can only appear co-called on the same read if they are within
    the read length of each other. This function greedily groups positions on
    the same amplicon such that (max - min) within each group ≤ max_distance.

    Args:
        amplicon_groups: Output from `group_positions_by_amplicon` — dict
            mapping amplicon_name -> sorted list of positions on that amplicon.
        max_distance: Maximum span within a batch (from group's first to last
            position, inclusive). Defaults to the value in cooc_config.yaml
            (`query.max_position_distance_bp`, currently 200bp).

    Returns:
        List of (amplicon_name, [positions]) tuples. Each tuple can be sent
        as a single LAPIS /aggregated query.
    """
    if max_distance is None:
        from utils.config import get_cooc_setting
        max_distance = get_cooc_setting("query.max_position_distance_bp", default=200)

    batches: List[Tuple[str, List[int]]] = []
    for amplicon_name, positions in amplicon_groups.items():
        if not positions:
            continue
        sorted_positions = sorted(positions)
        current_batch = [sorted_positions[0]]
        for pos in sorted_positions[1:]:
            span = pos - current_batch[0] + 1
            if span <= max_distance:
                current_batch.append(pos)
            else:
                batches.append((amplicon_name, current_batch))
                current_batch = [pos]
        batches.append((amplicon_name, current_batch))
    logger.info(
        f"Split into {len(batches)} query batches "
        f"(max_distance={max_distance}bp)"
    )
    return batches


def build_amp_dict_from_variants(
        variants: List[str],
        pango_loader,
        cowwid_variants=None,
) -> Dict[int, List[str]]:
    """
    Build an amp_dict mapping position -> list of tracked alt bases,
    aggregated across the signatures of selected panel variants.

    Args:
        variants: List of variant names (e.g. ['KP.2', 'NB.1.8.1']).
        pango_loader: PangoLoader instance for signature lookup.
        cowwid_variants: Optional fallback for reconstructed centroid nodes.

    Returns:
        Dict mapping position (int) -> list of alt bases (str), deduplicated.
        Example: {241: ['T', 'A'], 297: ['G'], ...}
    """
    amp_dict: Dict[int, set] = defaultdict(set)

    for variant in variants:
        # Always use pango_summary.json as the single source of truth.
        # pango_loader handles reconstructed nodes (e.g. BA.3.2) via
        # _fill_empty_node_signatures — no cowwid yaml fallback needed.
        signature = pango_loader.get_signature(variant)

        for mut in signature:
            # Parse "{pos}{alt}" — e.g. "241T"
            m = re.match(r"^(\d+)([ACGT-])$", mut)
            if not m:
                continue
            pos = int(m.group(1))
            alt = m.group(2)
            # Skip deletions if you want (matching remove_deletions=True)
            if alt == "-":
                continue
            amp_dict[pos].add(alt)

    return {pos: sorted(alts) for pos, alts in amp_dict.items()}
# ─────────────────────────────────────────────────────────────────────────
# Pair and stretch geometry (2026-11-10)
#
# process.cooc approximates "can one read carry both positions" with
# check.pair_span (300) and chains stretches with check.stretch_span (400).
# Both are distances standing in for a physical fact: a read comes off ONE
# amplicon. ArticV542inserts.bed says amplicon 74's insert is 22775-23119,
# 345 bases -- so pair_span 300 discards observable pairs -- while TA.1's
# stretch 23,148-23,604 is 456 bases spanning amplicons 75 and 76, which no
# read can produce and stretch_span 400 merged into one.
#
# group_positions_by_amplicon assigns each position to a single amplicon
# independently, so two positions of one pair can land in different groups
# when either sits in a tiling overlap. These helpers answer the pair
# question directly instead.
#
# Everything returns None / [] when the BED is unreadable, so callers keep
# their distance behaviour.
# ─────────────────────────────────────────────────────────────────────────

_INSERTS_DEFAULT = "/app_shared/data/ArticV542inserts.bed"
_inserts_cache: Dict[str, List[Tuple[str, int, int]]] = {}


def get_inserts(bed_path=None) -> List[Tuple[str, int, int]]:
    """load_amplicons on the inserts BED, cached. [] when unreadable.

    Resolution: the path given, then $AMPLICON_BED, then the container default.
    """
    import os
    # app/data next to app/process: the same file in the worker (/app_shared)
    # and in streamlit (/app), whatever the mount point (2026-11-10)
    _here = Path(__file__).resolve().parent.parent / "data" / "ArticV542inserts.bed"
    for p in (bed_path, os.environ.get("AMPLICON_BED"), _INSERTS_DEFAULT, _here):
        if not p:
            continue
        key = str(p)
        if key in _inserts_cache:
            if _inserts_cache[key]:
                return _inserts_cache[key]
            continue
        try:
            amps = load_amplicons(Path(key))
        except OSError:
            _inserts_cache[key] = []
            continue
        _inserts_cache[key] = amps
        if amps:
            return amps
    return []


def amplicon_of_pair(p1: int, p2: int,
                     amplicons: List[Tuple[str, int, int]]):
    """Name of an amplicon whose insert holds BOTH positions, else None.

    Lowest start wins, so a pair always groups into the same stretch however
    it is reached when amplicons overlap.
    """
    lo, hi = (p1, p2) if p1 <= p2 else (p2, p1)
    best = None
    for name, start, end in amplicons:
        if start <= lo and hi <= end:
            if best is None or start < best[1]:
                best = (name, start)
    return best[0] if best else None


def pair_reach(pos: int, amplicons: List[Tuple[str, int, int]]):
    """Highest position that could share a read with pos, else None.

    Used only to keep the early exit in specific_pairs' inner loop.
    """
    ends = [end for _n, start, end in amplicons if start <= pos <= end]
    return max(ends) if ends else None
