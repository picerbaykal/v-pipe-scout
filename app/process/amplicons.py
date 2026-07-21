"""Amplicon coordinates from BED files.

Used for amplicon-aware co-occurrence queries. Positions on the same
amplicon can be queried together because they come from the same physical
read; positions on different amplicons cannot.
"""

from pathlib import Path
from typing import List, Tuple
import logging

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