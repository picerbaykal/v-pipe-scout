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
    for name, start, end in amplicons:
        for pos in positions:
            if start <= pos <= end:
                grouped[name].append(pos)
                uncovered.discard(pos)
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
        # Get signature (same fallback logic as elsewhere)
        if variant in pango_loader._reconstructed_signatures and cowwid_variants and variant in cowwid_variants:
            signature = cowwid_variants[variant]
        else:
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