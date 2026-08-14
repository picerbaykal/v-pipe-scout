"""Cooc panel-completeness pipeline (worker side).

Given a location, date range, and panel variants, computes per-date panel
completeness by:
1. Building amp_dict + variant_signatures from panel variants.
2. Loading amplicons from BED, grouping panel positions by amplicon,
   splitting by read-length distance → query batches.
3. For each batch, calling LAPIS /aggregated for read-level co-occurrence.
4. Annotating results (confirmed_present/absent/classification).
5. Aggregating matched vs unexplained counts per date across batches.

Returns a dict shaped for JSON serialization back through Celery.
"""

import asyncio
import logging
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Callable, Dict, List, Optional

import pandas as pd

sys.path.insert(0, "/app_shared")

from api.pango_loader import PangoLoader, get_pango_summary_path
from api.wiseloculus import WiseLoculusLapis
from process.amplicons import (
    build_amp_dict_from_variants,
    group_positions_by_amplicon,
    load_amplicons,
    split_positions_by_distance,
)
from process.cooc import annotate_cooc_dataframe, panel_completeness_by_date
from utils.config import get_wiseloculus_url

logger = logging.getLogger(__name__)


def _load_cowwid_variants() -> dict:
    """
    Load cowwid signatures at first call; return empty dict if unavailable.

    Used as fallback for reconstructed centroid nodes (e.g. BA.3.2) where
    pango_summary has no direct sequences.
    """
    try:
        from api.signatures import get_variant_list
        variant_list = get_variant_list()
        result = {
            v.name: {m[1:] for m in v.signature_mutations if len(m) > 1}
            for v in variant_list.variants
        }
        logger.info(f"Loaded cowwid signatures for {len(result)} variants")
        return result
    except Exception as e:
        logger.warning(f"Could not load cowwid signatures: {e}")
        return {}


# existing
_COWWID_VARIANTS = _load_cowwid_variants()

# add these two after
def _load_all_lineage_signatures() -> dict:
    """Load all pango lineage signatures at startup. Cached for scanner use."""
    import re
    try:
        pl = PangoLoader(get_pango_summary_path())
        raw = pl.get_raw_data()
        result = {}
        for lineage in raw:
            try:
                sig = {m for m in pl.get_signature(lineage)
                       if re.match(r'^\d+[ACGT]$', m)}
                if len(sig) >= 2:
                    result[lineage] = sig
            except Exception:
                continue
        logger.info(f"Loaded signatures for {len(result)} pango lineages")
        return result
    except Exception as e:
        logger.warning(f"Could not load pango signatures: {e}")
        return {}


def _load_panel_parent_map() -> dict:
    """Load pango parent map at startup. Cached for scanner use."""
    try:
        pl = PangoLoader(get_pango_summary_path())
        return {lin: d.get("parent", "") for lin, d in pl.get_raw_data().items()}
    except Exception as e:
        logger.warning(f"Could not load parent map: {e}")
        return {}


_ALL_LINEAGE_SIGNATURES = None
_PANEL_PARENT_MAP = None


def get_all_lineage_signatures() -> dict:
    global _ALL_LINEAGE_SIGNATURES
    if _ALL_LINEAGE_SIGNATURES is None:
        _ALL_LINEAGE_SIGNATURES = _load_all_lineage_signatures()
    return _ALL_LINEAGE_SIGNATURES


def get_panel_parent_map() -> dict:
    global _PANEL_PARENT_MAP
    if _PANEL_PARENT_MAP is None:
        _PANEL_PARENT_MAP = _load_panel_parent_map()
    return _PANEL_PARENT_MAP

def _build_variant_signatures(
    variants: List[str],
    pango_loader: PangoLoader,
    cowwid_variants: Dict[str, set],
) -> Dict[str, set]:
    """
    Build per-variant signature sets in "{pos}{alt}" format matching amp_dict.

    Deletions are excluded (matching remove_deletions=True in deconvolution).
    """
    sigs: Dict[str, set] = {}
    for variant in variants:
        if variant in pango_loader._reconstructed_signatures and variant in cowwid_variants:
            sig = cowwid_variants[variant]
        else:
            sig = pango_loader.get_signature(variant)
        # Keep only substitution entries (skip deletions ending in "-")
        sigs[variant] = {m for m in sig if re.match(r"^\d+[ACGT]$", m)}
    return sigs


def run_cooc_panel_completeness(
    location: str,
    start_date: datetime,
    end_date: datetime,
    variants: List[str],
    bed_path: Optional[str] = None,
    progress_callback: Optional[Callable[[int, str], None]] = None,
) -> dict:
    """
    Compute per-date panel completeness for one location.

    Args:
        location: Location name e.g. "Lugano (TI)"
        start_date, end_date: Datetime bounds (inclusive).
        variants: Panel variant names.
        bed_path: Path to amplicon BED. Defaults to /app_shared/data/SARS-CoV-2.insert.bed.
        progress_callback: Optional callback(step, message) for progress reporting.

    Returns:
        Dict with keys: location, dates, matched_counts, unexplained_counts, completeness.
        All list values are aligned by index (one entry per date).
    """
    if bed_path is None:
        bed_path = "/app_shared/data/SARS-CoV-2.insert.bed"

    from utils.config import get_cooc_setting

    scope_weeks = get_cooc_setting("scope.weeks", default=None)
    if scope_weeks:
        from datetime import timedelta
        clamped = end_date - timedelta(weeks=int(scope_weeks))
        if clamped > start_date:
            logger.info(
                f"[cooc][{location}] scope.weeks={scope_weeks} — "
                f"clamping start {start_date.date()} → {clamped.date()}"
            )
            start_date = clamped

    def _progress(step: int, msg: str):
        if progress_callback:
            progress_callback(step, msg)
        logger.info(f"[cooc][{location}] step {step}: {msg}")

    _progress(1, f"Building amp_dict + signatures for {len(variants)} variants")
    pango_loader = PangoLoader(get_pango_summary_path())
    cowwid_variants = _COWWID_VARIANTS
    reference_variants = get_cooc_setting("scope.reference_variants", default=None)
    if reference_variants is None and get_cooc_setting("scope.use_tracked_variants", default=False):
        reference_variants = sorted(cowwid_variants.keys())
    amp_dict = build_amp_dict_from_variants(
        reference_variants or variants, pango_loader, cowwid_variants
    )
    variant_signatures = _build_variant_signatures(
        variants, pango_loader, cowwid_variants
    )
    logger.info(
        f"[cooc][{location}] amp_dict from "
        f"{'reference list' if reference_variants else 'panel'}: "
        f"{len(amp_dict)} positions"
    )

    if get_cooc_setting("scope.discriminating_positions_only", default=False):
        sig_list = [s for s in variant_signatures.values() if s]
        if len(sig_list) >= 2:
            shared = set.intersection(*sig_list)
            before = len(amp_dict)
            amp_dict = {
                pos: alts for pos, alts in amp_dict.items()
                if not all(f"{pos}{a}" in shared for a in alts)
            }
            logger.info(
                f"[cooc][{location}] discriminating filter: "
                f"{before} → {len(amp_dict)} positions"
            )

    positions = set(amp_dict.keys())
    logger.info(
        f"[cooc][{location}] Panel: {len(positions)} positions across "
        f"{len(variants)} variants"
    )

    _progress(2, "Grouping positions by amplicon + read-length distance")
    # TODO: when SaneQL REST endpoint deploys, replace per-amplicon batching with:
    # 1. Coverage scan: query positions with non-N reads > MIN_COVERAGE threshold
    #    (skip positions with no real coverage — no point querying them)
    # 2. Single SaneQL query for all covered positions in one request:
    #    default
    #    .filter(locationName == location & samplingDate in date_range)
    #    .map({"pos1" := main.at(pos1), "pos2" := main.at(pos2), ...covered_positions...})
    #    .groupBy({count := count()}, {pos1, pos2, ..., samplingDate})
    #    This gives full genome co-occurrence without BED file or amplicon batching.
    #    For "possibly new" detection: expand position set beyond known variant signatures
    #    to all covered positions (~500-1000 on a typical date) for comprehensive novel
    #    variant detection. Zero-coverage positions are excluded by the coverage scan,
    #    keeping the query fast (~2-5s estimated based on 130-position benchmark).
    amplicons = load_amplicons(Path(bed_path))
    amp_groups = group_positions_by_amplicon(positions, amplicons)
    batches = list(amp_groups.items())
    amplicons = load_amplicons(Path(bed_path))
    # new — one batch per amplicon, all positions together
    amp_groups = group_positions_by_amplicon(positions, amplicons)
    batches = list(amp_groups.items())  # (amplicon_name, positions) directly
    logger.info(
        f"[cooc][{location}] {len(batches)} query batches "
        f"(from {len(amp_groups)} amplicons)"
    )

    _progress(3, f"Querying LAPIS for {len(batches)} batches")
    client = WiseLoculusLapis(get_wiseloculus_url())

    async def _query_all_batches():
        dates = await client._get_sampling_dates(location, (start_date, end_date))
        logger.info(f"[cooc][{location}] {len(dates)} sampling dates")
        if not dates:
            return []

        sem = asyncio.Semaphore(8)  # cap concurrent batches

        async def _one_batch(i, amp_name, batch_positions):
            async with sem:
                df = await client.get_cooccurrence(
                    locationName=location,
                    date_range=(start_date, end_date),
                    positions=batch_positions,
                    dates=dates,
                )
            logger.info(f"[cooc][{location}] batch {i}/{len(batches)} — {len(df)} rows")
            if df.empty:
                return None
            annotated = annotate_cooc_dataframe(
                df, batch_positions, amp_dict, variant_signatures
            )
            if "classification" in annotated.columns:
                logger.info(
                    f"[cooc][{location}] batch {i} rows={annotated['classification'].value_counts().to_dict()} "
                    f"reads={annotated.groupby('classification')['count'].sum().to_dict()}"
                )
            else:
                logger.info(f"[cooc][{location}] batch {i} not annotated (early exit)")

            per_date = panel_completeness_by_date(annotated)
            unexplained = annotated[annotated["classification"] == "unexplained"][
                ["date", "count", "confirmed_present"]
            ].copy() if "classification" in annotated.columns else pd.DataFrame()
            if per_date.empty:
                return None
            return per_date, unexplained

        tasks = [
            _one_batch(i, amp, pos)
            for i, (amp, pos) in enumerate(batches, 1)
        ]
        done = await asyncio.gather(*tasks, return_exceptions=True)

        completeness_results = []
        pattern_results = []
        for r in done:
            if isinstance(r, Exception):
                logger.error(f"[cooc][{location}] batch failed: {r}")
                continue
            if r is not None:
                per_date, unexplained = r
                completeness_results.append(per_date)
                if not unexplained.empty:
                    pattern_results.append(unexplained)
        return completeness_results, pattern_results

    per_batch_results, pattern_results = asyncio.run(_query_all_batches())

    _progress(4, "Aggregating across batches")
    if not per_batch_results:
        logger.warning(f"[cooc][{location}] No batches returned data")
        return {
            "location": location,
            "dates": [],
            "matched_counts": [],
            "unexplained_counts": [],
            "completeness": [],
            "unexplained_patterns": [],
        }

    combined = pd.concat(per_batch_results, ignore_index=True)

    if pattern_results:
        all_patterns = pd.concat(pattern_results, ignore_index=True)
        all_patterns["pattern_key"] = all_patterns["confirmed_present"].apply(tuple)
        unexplained_agg = (
            all_patterns.groupby(["date", "pattern_key"])["count"]
            .sum().reset_index()
            .rename(columns={"pattern_key": "confirmed_present"})
        )
        unexplained_agg["confirmed_present"] = unexplained_agg["confirmed_present"].apply(list)
    else:
        unexplained_agg = pd.DataFrame(columns=["date", "count", "confirmed_present"])

    per_date = combined.groupby("date", as_index=False)[
        ["matched_count", "unexplained_count"]
    ].sum()
    per_date["completeness"] = per_date["matched_count"] / (
        per_date["matched_count"] + per_date["unexplained_count"]
    ).replace(0, pd.NA)
    per_date = per_date.sort_values("date").reset_index(drop=True)

    return {
        "location": location,
        "dates": per_date["date"].astype(str).tolist(),
        "matched_counts": per_date["matched_count"].astype(int).tolist(),
        "unexplained_counts": per_date["unexplained_count"].astype(int).tolist(),
        "completeness": per_date["completeness"].astype(float).tolist(),
        "unexplained_patterns": unexplained_agg.to_dict("records"),
    }