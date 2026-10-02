"""Cooc panel-completeness pipeline (worker side).

Given a location, date range, and panel variants, computes per-date panel
completeness by:
1. Building amp_dict + variant_signatures from panel variants.
2. Loading amplicons from BED, grouping panel positions by amplicon → query
   batches (one batch per amplicon; positions within an amplicon co-occur).
3. For each batch, calling LAPIS /aggregated for read-level co-occurrence.
4. Annotating results (confirmed_present/absent/classification).
5. Aggregating matched vs unexplained counts per date across batches.

Returns a dict shaped for JSON serialization back through Celery.

Note (2026-08): the LAPIS co-occurrence endpoint (PR #1768, `[position]`
bracket fields in /aggregated) is deployed on WASAP, and query cost is
confirmed independent of position count. The `scope.weeks` clamp has been
removed from the default path on that basis.

CAVEAT (2026-08-26, unverified): the "~5s full sweep" figure above has NOT
been reproduced end-to-end. Measured cold benchmarks for a full BED-based
sweep (this file's code path, one city, ~421 positions) were 11.7-17.9s,
not ~5s — see cooc_investigation_summary.md. Separately, this file's
_fetch_cooccurrence_for_date calls now use the `date` field instead of
`samplingDate` (temporary LAPIS-side hack, see api/wiseloculus.py), which
should help, but the combined fetch+classify wall time with this change has
not yet been measured for a full sweep. Treat "~5s" as aspirational until
re-benchmarked; do not assume the scope.weeks clamp removal is safe for the
worst-case scenario (all variants, all cities, 6 months) without testing it.

2026-09: adds the new per-variant co-occurrence check ("panel_check" in the
result): specific markers chosen with own-family-excluded carrier counts
(process.cooc.specific_markers) and raw per-date marker counts
(process.cooc.accumulate_check_stats). Verdicts are derived in the UI via
process.cooc.check_verdicts so thresholds can be retuned without a re-scan.
(The legacy "panel_presence" block was removed on 2026-09-29.)
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
)
from process.cooc import (annotate_cooc_dataframe, panel_completeness_by_date,
                          specific_markers, accumulate_check_stats,
                          aggregate_unexplained)
from utils.config import get_wiseloculus_url

logger = logging.getLogger(__name__)

# Concurrency cap for parallel batch queries. Benchmarks (2026-08) show
# near-linear scaling to 16 workers with no server-side rate-limiting;
# 16 roughly halves wall time vs 8 for a full sweep.
BATCH_CONCURRENCY = 8


def _load_cowwid_variants() -> dict:
    """
    Load cowwid signatures at first call; return empty dict if unavailable.

    Used as fallback for reconstructed centroid nodes (e.g. BA.3.2) where
    pango_summary has no direct sequences.
    """
    try:
        from api.signatures import get_variant_list
        variant_list = get_variant_list()
        result = {v.name: set() for v in variant_list.variants}
        logger.info(f"Loaded cowwid variant names: {len(result)}")
        return result
    except Exception as e:
        logger.warning(f"Could not load cowwid signatures: {e}")
        return {}


_COWWID_VARIANTS = _load_cowwid_variants()


def _load_all_lineage_signatures() -> dict:
    """Load all pango lineage signatures at startup. Cached for scanner use."""
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
        # Always use pango_summary.json — single source of truth.
        sig = pango_loader.get_signature(variant)
        # Keep only substitution entries (skip deletions ending in "-")
        sigs[variant] = {m for m in sig if re.match(r"^\d+[ACGT]$", m)}
    return sigs


def _sig_positions(sig: set) -> set:
    out = set()
    for m in sig or ():
        mm = re.match(r"^(\d+)[ACGT]$", m)
        if mm:
            out.add(int(mm.group(1)))
    return out


def add_panel_positions(amp_dict: Dict[int, list], variant_signatures: Dict[str, set]) -> int:
    """Add every panel variant's substitutions ("{pos}{alt}") to amp_dict
    in place; returns how many mutations were new."""
    added = 0
    for sig in variant_signatures.values():
        for m in sig or ():
            mm = re.match(r"^(\d+)([ACGT])$", m)
            if not mm:
                continue
            pos, alt = int(mm.group(1)), mm.group(2)
            alts = amp_dict.setdefault(pos, [])
            if alt in alts:
                continue
            if isinstance(alts, set):
                alts.add(alt)
            else:
                alts.append(alt)
            added += 1
    return added


def data_positions(location: str, start_date: datetime, end_date: datetime,
                   min_share: float, min_cov: int) -> Dict[int, set]:
    """Positions where the DATA shows a substitution: every mutation with >=
    min_share of >= min_cov covering reads in this city and window (LAPIS
    nucleotideMutations, pooled over the window). {position: {new bases}}.

    Used by the deep scan (2026-10-02): today's positions come from lists
    (cowwid names + panel), so a lineage's newest mutations — its best markers —
    and undesignated mutations were never fetched. The thresholds are the check's
    own (absent_freq, min_cov): a position is added when its mutation is "not
    absent" on enough reads. Deletions are left out, as everywhere."""
    try:
        from api.wiseloculus import MutationType
    except ImportError:
        from api.signatures import MutationType
    client = WiseLoculusLapis(get_wiseloculus_url())
    df = asyncio.run(client.sample_mutations(MutationType.NUCLEOTIDE, (start_date, end_date),
                                             locationName=location, min_proportion=min_share))
    if df is None or df.empty:
        return {}
    if len(df) >= 10000:
        logger.warning(f"[cooc][{location}] nucleotideMutations hit its 10,000-row limit — "
                       "some positions may be missing")
    out: Dict[int, set] = {}
    for _, r in df.iterrows():
        alt = str(r.get("mutationTo", ""))
        if alt in ("A", "C", "G", "T") and int(r.get("coverage", 0) or 0) >= min_cov:
            out.setdefault(int(r["position"]), set()).add(alt)
    return out


def coverage_from_rows(rows: List[dict], positions) -> Dict[str, int]:
    """{position: reads covering it} for one LAPIS co-occurrence answer (one
    date, one batch). A read covers a position when its base there is not N
    and not a deletion; rows are read patterns with a "count"."""
    if not rows:
        return {}
    df = pd.DataFrame(rows)
    cols = [f"[{p}]" for p in positions if f"[{p}]" in df.columns]
    if not cols or "count" not in df.columns:
        return {}
    cnt = pd.to_numeric(df["count"], errors="coerce").fillna(0)
    covered = df[cols].notna() & ~df[cols].isin(["N", "-"])
    tot = covered.mul(cnt, axis=0).sum()
    return {c[1:-1]: int(n) for c, n in tot.items() if n > 0}


def run_cooc_panel_completeness(
    location: str,
    start_date: datetime,
    end_date: datetime,
    variants: List[str],
    bed_path: Optional[str] = None,
    progress_callback: Optional[Callable[[int, str], None]] = None,
    extra_positions: Optional[Dict[int, set]] = None,
) -> dict:
    """
    Compute per-date panel completeness for one location.

    Args:
        location: Location name e.g. "Lugano (TI)"
        start_date, end_date: Datetime bounds (inclusive).
        variants: Panel variant names.
        bed_path: Path to amplicon BED. Defaults to /app_shared/data/ArticV542inserts.bed.
        progress_callback: Optional callback(step, message) for progress reporting.
        extra_positions: {position: {bases}} added to the completeness positions
            (the deep scan's positions from the data, see data_positions).

    Returns:
        Dict with keys: location, dates, matched_counts, unexplained_counts, completeness,
        unexplained_patterns, position_coverage ({date: {position: reads}})
        and panel_check (co-occurrence check:
        per variant {"markers": [...], "per_date": {date: {marker: [cov, hit,
        link_n, link_ok]}}}). List values are aligned by index (one per date).
    """
    if bed_path is None:
        bed_path = "/app_shared/data/ArticV542inserts.bed"

    from utils.config import get_cooc_setting


    # NOTE: the former `scope.weeks` start-date clamp has been removed.
    # Co-occurrence queries are now fast enough (~5s full sweep) to run over
    # the full requested date range without truncation. If a genuine upper
    # bound on range is ever needed again, enforce it in the UI/date-picker
    # rather than silently clamping here.

    def _progress(step: int, msg: str):
        if progress_callback:
            progress_callback(step, msg)
        logger.info(f"[cooc][{location}] step {step}: {msg}")

    _progress(1, f"Building amp_dict + signatures for {len(variants)} variants")
    pango_loader = PangoLoader(get_pango_summary_path())
    cowwid_variants = _COWWID_VARIANTS
    reference_variants = get_cooc_setting("scope.reference_variants", default=None)
    # amp_dict_cutoff: build amp_dict from ALL lineages designated on/after
    # this date (covers circulating variants + sublineages, not just cowwid).
    # Positions are the cheap axis (benchmarked ~2x for 5.7x more positions).
    amp_dict_cutoff = get_cooc_setting("scope.amp_dict_cutoff", default=None)
    if reference_variants is None and amp_dict_cutoff:
        _raw = pango_loader.get_raw_data()
        reference_variants = sorted(
            l for l, d in _raw.items()
            if (d.get("designationDate", "") or "") >= amp_dict_cutoff
            and pango_loader.get_signature(l)
        )
        logger.info(
            f"[cooc][{location}] amp_dict_cutoff={amp_dict_cutoff}: "
            f"{len(reference_variants)} lineages"
        )
    if reference_variants is None and get_cooc_setting("scope.use_tracked_variants", default=False):
        reference_variants = sorted(cowwid_variants.keys())
    amp_dict = build_amp_dict_from_variants(
        reference_variants or variants, pango_loader, cowwid_variants
    )
    variant_signatures = _build_variant_signatures(
        variants, pango_loader, cowwid_variants
    )

    # New check (Part A): specific markers per panel variant, carriers counted
    # outside the variant's own family. Degrades to "no markers" on error so
    # the completeness scan never breaks because of it.
    try:
        check_markers = specific_markers(
            variant_signatures, get_all_lineage_signatures(), get_panel_parent_map())
    except Exception as e:
        logger.warning(f"[cooc][{location}] specific_markers failed: {e}")
        check_markers = {v: [] for v in variants}
    check_stats: Dict[str, dict] = {}
    logger.info(
        f"[cooc][{location}] check markers: "
        + ", ".join(f"{v}={len(m)}" for v, m in check_markers.items())
    )
    logger.info(
        f"[cooc][{location}] amp_dict from "
        f"{'reference list' if reference_variants else 'panel'}: "
        f"{len(amp_dict)} positions"
    )

    # Panel variants' own mutations are completeness positions too
    # (2026-09-30): before, a panel variant outside the reference list (e.g.
    # PJ.2) had its own positions fetched for the check only — completeness
    # and the scanner never saw them, so its parent's reads counted as
    # "explained" by it and nothing beyond it could be named.
    _n_before = len(amp_dict)
    _added_muts = add_panel_positions(amp_dict, variant_signatures)
    logger.info(
        f"[cooc][{location}] panel signatures: +{_added_muts} mutations, "
        f"+{len(amp_dict) - _n_before} positions in completeness"
    )

    if extra_positions:
        _n_before = len(amp_dict)
        _n_mut = 0
        for _p, _alts in extra_positions.items():
            _cur = amp_dict.setdefault(int(_p), [])
            for _a in _alts:
                if _a not in _cur:
                    _cur.add(_a) if isinstance(_cur, set) else _cur.append(_a)
                    _n_mut += 1
        logger.info(f"[cooc][{location}] data positions: +{_n_mut} mutations, "
                    f"+{len(amp_dict) - _n_before} positions")

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

    # Query positions = completeness positions (amp_dict) PLUS every panel
    # variant's signature positions, so each check marker and its neighbours are
    # read even when amp_dict is built from the tracked list. The extra
    # positions are invisible to completeness: annotate only sees amp_dict ones.
    completeness_positions = set(amp_dict.keys())
    check_positions = set()
    for v, mks in check_markers.items():
        if mks:
            check_positions |= _sig_positions(variant_signatures.get(v, set()))
    positions = completeness_positions | check_positions
    logger.info(
        f"[cooc][{location}] Panel: {len(completeness_positions)} completeness "
        f"positions + {len(positions - completeness_positions)} check-only positions "
        f"across {len(variants)} variants"
    )

    _progress(2, "BED-free: preparing all positions for query")
    # BED-free: one batch containing all positions, one query per date.
    # Benchmarked 2026-09: 8-65x faster than BED-based (90 batches x dates),
    # identical completeness (0.9935 vs 0.9935), more correct (full
    # confirmed_absent set catches contradictions BED-based misses).
    # Cross-amplicon positions return N automatically — physics enforces
    # amplicon scoping, not the BED file. BED still used for UI labeling only.
    # Chunk positions under LAPIS's per-query field limit (~700 measured).
    # Reads span ~400bp so co-occurring positions are always genomically
    # adjacent; contiguous chunks preserve every real co-occurrence.
    _CHUNK = 500
    _sorted_pos = sorted(positions)
    if len(_sorted_pos) <= _CHUNK:
        batches = [("all", _sorted_pos)]
    else:
        batches = [
            (f"chunk{i//_CHUNK}", _sorted_pos[i:i + _CHUNK])
            for i in range(0, len(_sorted_pos), _CHUNK)
        ]
    logger.info(
        f"[cooc][{location}] BED-free: {len(batches)} batch(es), "
        f"{len(positions)} positions"
    )

    _progress(3, f"Querying LAPIS for {len(batches)} batches")
    client = WiseLoculusLapis(get_wiseloculus_url())

    async def _query_all_batches():
        dates = await client._get_sampling_dates(location, (start_date, end_date))
        logger.info(f"[cooc][{location}] {len(dates)} sampling dates")
        if not dates:
            return [], [], {}

        # Flat concurrent pool: pre-build ALL (batch, date) pairs and run them
        # in one shared session with a single semaphore. This avoids the nested
        # fan-out (16 batch-sessions each opening 35 connections) that fought
        # the per-host connection cap. Instead, all n_batches × n_dates queries
        # share one connector and are serialized only by BATCH_CONCURRENCY.
        # For 1 batch × N dates = N queries total (e.g. 54 dates = 54 queries).
        # BED-free: 8-65x faster than the previous 90-batch structure.
        import aiohttp
        from api.wiseloculus import MAX_CONCURRENT_CONNECTIONS, MAX_CONNECTIONS_PER_HOST

        connector = aiohttp.TCPConnector(
            limit=MAX_CONCURRENT_CONNECTIONS,
            limit_per_host=MAX_CONNECTIONS_PER_HOST,
        )
        timeout = aiohttp.ClientTimeout(total=120)
        sem = asyncio.Semaphore(BATCH_CONCURRENCY)

        # Accumulate results. Multiple chunks per date are merged before
        # classification (each chunk covers a genomic region; a read's
        # confirmed_present spans only its own chunk since reads are short).
        per_date_results = []
        per_date_unexplained = []
        # reads covering each fetched position per date (base not N / deletion),
        # ALL reads, explained or not: tells the scanner whether a later sample
        # could have shown a one-day finding again (2026-10-02)
        position_coverage: Dict[str, Dict[str, int]] = {}

        async def _one_query(session, batch_idx, batch_positions, date_str):
            async with sem:
                rows = await client._fetch_cooccurrence_for_date(
                    session, location, date_str, batch_positions
                )
            if not rows:
                return
            _slot = position_coverage.setdefault(str(date_str)[:10], {})
            for _k, _n in coverage_from_rows(rows, batch_positions).items():
                _slot[_k] = _slot.get(_k, 0) + _n
            # New check (Part B): raw per-date marker counts. Pure and
            # synchronous, so no interleaving between concurrent tasks.
            accumulate_check_stats(rows, batch_positions, variant_signatures,
                                   check_markers, check_stats)
            comp_positions = [p for p in batch_positions if p in completeness_positions]
            if not comp_positions:
                return
            import pandas as _pd
            df = _pd.DataFrame(rows)
            annotated = annotate_cooc_dataframe(
                df, comp_positions, amp_dict, variant_signatures
            )
            per_date = panel_completeness_by_date(annotated)
            if not per_date.empty:
                per_date_results.append(per_date)
            if "classification" in annotated.columns:
                _cols = ["date", "count", "confirmed_present", "confirmed_absent"]
                if "near" in annotated.columns:
                    _cols.append("near")
                unexp = annotated[annotated["classification"] == "unexplained"][
                    _cols
                ].copy()
                if not unexp.empty:
                    per_date_unexplained.append(unexp)
            del df, rows, annotated

        async with aiohttp.ClientSession(
            timeout=timeout, connector=connector
        ) as session:
            tasks = [
                _one_query(session, bi, bpos, d)
                for bi, (_, bpos) in enumerate(batches)
                for d in dates
            ]
            done = await asyncio.gather(*tasks, return_exceptions=True)

        n_errors = sum(1 for r in done if isinstance(r, Exception))
        if n_errors:
            logger.error(f"[cooc][{location}] {n_errors}/{len(tasks)} queries failed")
        logger.info(
            f"[cooc][{location}] processed {len(per_date_results)} dates, "
            f"{len(per_date_unexplained)} with unexplained patterns"
        )
        return per_date_results, per_date_unexplained, position_coverage

    try:
        loop = asyncio.get_event_loop()
        if loop.is_running():
            import concurrent.futures
            with concurrent.futures.ThreadPoolExecutor() as pool:
                future = pool.submit(asyncio.run, _query_all_batches())
                per_batch_results, pattern_results, position_coverage = future.result()
        else:
            per_batch_results, pattern_results, position_coverage = loop.run_until_complete(
                _query_all_batches())
    except RuntimeError:
        per_batch_results, pattern_results, position_coverage = asyncio.run(_query_all_batches())

    panel_check = {
        v: {"markers": list(check_markers.get(v, [])),
            "per_date": check_stats.get(v, {})}
        for v in variants
    }

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
            "panel_check": panel_check,
            "position_coverage": position_coverage,
        }

    combined = pd.concat(per_batch_results, ignore_index=True)

    # one row per (date, present, absent) — the scanner needs the absent side
    # to tell which lineages a read can come from (process.cooc)
    unexplained_agg = aggregate_unexplained(pattern_results)

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
        "panel_check": panel_check,
        # {date: {position: reads covering it}} — see _query_all_batches
        "position_coverage": position_coverage,
    }