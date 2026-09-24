"""Panel-completeness computation from LAPIS co-occurrence data.

Given a DataFrame of read-level base combinations at target positions
(from WiseLoculusLapis.get_cooccurrence), compute:

- Per-row translation to confirmed_present / confirmed_absent (matching
  the manual scan's semantics).
- Per-row classification as matched (some panel variant explains all
  observed mutations) or unexplained (no variant explains all).
- Per-date aggregation of panel completeness.

Two configurable knobs live in app/config/cooc_config.yaml:
- filter.require_all_positions_covered: skip rows where any queried position
  is uncovered (N).
- filter.include_deletion_states: whether to treat "-" as a valid observed
  base (kept in signatures) or as uncovered.
"""

import logging
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

from utils.config import get_cooc_setting

logger = logging.getLogger(__name__)


# ── Per-row translation ─────────────────────────────────────────────────

def row_to_confirmed_sets(
    row: dict,
    positions: List[int],
    amp_dict: Dict[int, List[str]],
    include_deletion_states: bool = False,
) -> Tuple[Set[str], Set[str]]:
    """
    Translate a single LAPIS cooc row into (confirmed_present, confirmed_absent).

    Matches the manual scan's semantics: for each position tracked in
    amp_dict, if the read shows a tracked alt base → all tracked alts at
    that position go into confirmed_present. If the read shows reference
    (non-N, non-tracked) → all tracked alts go into confirmed_absent.
    If the read shows N (or "-" when include_deletion_states is False),
    the position is uncovered and contributes to neither set.

    Args:
        row: One row from the LAPIS cooc DataFrame, e.g.
            {"[241]": "T", "[297]": "G", "count": 15000, "date": "2025-11-09"}
        positions: The positions queried in this batch.
        amp_dict: Panel amp_dict — position -> list of tracked alt bases.
        include_deletion_states: If False, "-" is treated as uncovered.

    Returns:
        (confirmed_present, confirmed_absent) as sets of "{pos}{alt}" strings.
    """
    confirmed_present: Set[str] = set()
    confirmed_absent: Set[str] = set()
    uncovered = {"N"} if include_deletion_states else {"N", "-"}

    for pos in positions:
        base = row.get(f"[{pos}]", "N")
        if base in uncovered:
            continue
        tracked_alts = amp_dict.get(pos, [])
        if not tracked_alts:
            continue
        if base in tracked_alts:
            # Record only the base actually observed. A read is one molecule
            # and carries one base per position — emitting every tracked alt
            # would describe a read that cannot exist, and such patterns match
            # no signature and inflate the unexplained fraction.
            confirmed_present.add(f"{pos}{base}")
        else:
            # Reference observed — mark all tracked alts as absent
            for alt in tracked_alts:
                confirmed_absent.add(f"{pos}{alt}")

    return confirmed_present, confirmed_absent


# ── Classification ──────────────────────────────────────────────────────

def classify_pattern(
    confirmed_present: Set[str],
    confirmed_absent: Set[str],
    variant_signatures: Dict[str, Set[str]],
) -> str:
    """
    Classify a pattern by whether SOME panel variant explains it.

    A variant explains a pattern when both hold:
      - every observed mutation is in its signature
      - none of its signature mutations were observed to be absent

    The second condition carries most of the discriminating power: a read
    covering a position where the variant requires a mutation, and showing
    reference there, is positive evidence against that variant.

    Args:
        confirmed_present: "{pos}{alt}" strings observed on this read.
        confirmed_absent: "{pos}{alt}" strings looked for and not found.
        variant_signatures: variant_name -> signature mutation set.

    Returns:
        "uninformative" (fewer than 2 observed mutations),
        "matched" (some variant explains it), or
        "unexplained" (no variant does — the panel-gap signal).
    """
    if len(confirmed_present) < 2:
        return "uninformative"
    for sig in variant_signatures.values():
        if confirmed_present.issubset(sig) and not (sig & confirmed_absent):
            return "matched"
    return "unexplained"


# ── Per-batch DataFrame processing ──────────────────────────────────────

def annotate_cooc_dataframe(
    df: pd.DataFrame,
    positions: List[int],
    amp_dict: Dict[int, List[str]],
    variant_signatures: Dict[str, Set[str]],
) -> pd.DataFrame:
    """
    Annotate each row of a cooc DataFrame with confirmed_present /
    confirmed_absent / classification.

    Args:
        df: DataFrame from WiseLoculusLapis.get_cooccurrence
            (columns: date, count, [pos1], [pos2], ...).
        positions: The positions queried in this batch.
        amp_dict: Panel amp_dict.
        variant_signatures: Dict of variant_name -> signature mutation set,
            in the same "{pos}{alt}" format as amp_dict produces.

    Returns:
        A new DataFrame with additional columns:
            confirmed_present, confirmed_absent, classification.
        Rows where require_all_positions_covered is True and any queried
        position is N/- are dropped.
    """
    if df.empty:
        return df.copy()

    require_covered = get_cooc_setting("filter.require_all_positions_covered", False)
    include_dels = get_cooc_setting("filter.include_deletion_states", False)
    uncovered_bases = {"N"} if include_dels else {"N", "-"}

    annotated_rows = []
    dropped = 0
    for row in df.to_dict("records"):
        if require_covered:
            if any(row.get(f"[{p}]", "N") in uncovered_bases for p in positions):
                dropped += 1
                continue
        cp, ca = row_to_confirmed_sets(row, positions, amp_dict, include_dels)
        classification = classify_pattern(cp, ca, variant_signatures)
        annotated_rows.append({
            "date": row["date"],
            "count": row["count"],
            "confirmed_present": sorted(cp),
            "confirmed_absent": sorted(ca),
            "classification": classification,
        })

    if dropped:
        logger.info(
            f"coverage filter dropped {dropped}/{len(df)} rows "
            f"({dropped / len(df):.1%}) for positions {positions}"
        )

    return pd.DataFrame(annotated_rows)


# ── Per-date panel completeness ─────────────────────────────────────────

def panel_completeness_by_date(
    annotated_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Aggregate an annotated cooc DataFrame into per-date panel completeness.

    Ignores "uninformative" patterns (0-1 mutations). Only the balance
    between "matched" and "unexplained" matters.

    Args:
        annotated_df: Output of annotate_cooc_dataframe.

    Returns:
        DataFrame with columns:
            date, matched_count, unexplained_count, completeness
        where completeness = matched / (matched + unexplained).
        Returns empty DataFrame if annotated_df is empty.
    """
    if annotated_df.empty:
        return pd.DataFrame(columns=[
            "date", "matched_count", "unexplained_count", "completeness"
        ])

    informative = annotated_df[annotated_df["classification"] != "uninformative"]
    if informative.empty:
        return pd.DataFrame(columns=[
            "date", "matched_count", "unexplained_count", "completeness"
        ])

    grouped = informative.groupby(["date", "classification"])["count"].sum().unstack(fill_value=0)
    if "matched" not in grouped.columns:
        grouped["matched"] = 0
    if "unexplained" not in grouped.columns:
        grouped["unexplained"] = 0

    result = grouped.reset_index().rename(columns={
        "matched": "matched_count",
        "unexplained": "unexplained_count",
    })
    result["completeness"] = result["matched_count"] / (
        result["matched_count"] + result["unexplained_count"]
    ).replace(0, pd.NA)
    return result[["date", "matched_count", "unexplained_count", "completeness"]]
# ── Per-panel-variant presence / co-coverage ("not found in WW") ───────────
# Empirical presence/absence for each PANEL variant, from the same confirmed_
# present / confirmed_absent computed above. Amplicon scoping is automatic:
# cross-amplicon positions come back N together, so they never both appear in
# confirmed_present / confirmed_absent on one read.
import re as _re_presence


# A distinctive mutation must be globally RARE to count — carried by at most this
# many pango lineages. Matches the scanner's ★ discriminating bar. Positions more
# common than this (E484K, the N-gene triplet, ...) are shared with currently-
# circulating lineages, so they cannot confirm a specific (possibly extinct)
# panel variant. Tune here if the scanner's STAR_CARRIER_MAX changes.
_DISTINCT_STAR_MAX = 30

_GLOBAL_CARRIER = {}


def _global_carrier_counts(pango_loader):
    """substitution -> number of pango lineages carrying it. Computed once per
    process from the pango summary (the same ~5k-lineage reference every scan
    sees), then cached. Degrades to {} on any error, which restores the old
    panel-relative behaviour rather than breaking the scan."""
    global _GLOBAL_CARRIER
    if _GLOBAL_CARRIER:
        return _GLOBAL_CARRIER
    carrier = {}
    try:
        for _lin in pango_loader.get_raw_data():
            for _m in (pango_loader.get_signature(_lin) or []):
                if _re_presence.match(r"^\d+[ACGT]$", _m):
                    carrier[_m] = carrier.get(_m, 0) + 1
    except Exception:
        return {}
    _GLOBAL_CARRIER = carrier
    return carrier


def distinctive_within_panel(variant_signatures, carrier_counts=None,
                             star_max=_DISTINCT_STAR_MAX):
    """variant -> mutations that both (a) NO OTHER panel member carries and, when
    `carrier_counts` is given, (b) are GLOBALLY RARE (carried by <= star_max
    lineages). These are the positions co-occurrence uses to test a variant's
    presence.

    Panel-relative uniqueness alone is unsafe in a small panel: a variant's
    "distinctive vs the other members" set can be dominated by mutations that
    modern circulating lineages carry, which then co-occur in the reads and
    falsely confirm an extinct variant (B.1.1.7, B.1.351). Intersecting with the
    global-rarity bar leaves only genuinely-defining positions, so an extinct
    variant is judged on mutations that are actually absent from today's WW.

    A variant with no globally-rare distinctive mutation gets an empty set — it
    cannot be told apart from the circulating background, so it is left
    unconfirmable (the verdict layer reads that as a blind spot / not found)
    rather than confirmed off shared signal. Without carrier_counts the old
    panel-relative behaviour is preserved."""
    panel = list(variant_signatures)
    out = {}
    for v in panel:
        others = set()
        for w in panel:
            if w != v:
                others |= (variant_signatures.get(w) or set())
        d = (variant_signatures.get(v) or set()) - others
        if carrier_counts:
            d = {m for m in d if carrier_counts.get(m, 0) <= star_max}
        out[v] = d
    return out


def _presence_pos(m: str) -> int:
    mm = _re_presence.match(r"^(\d+)", m)
    return int(mm.group(1)) if mm else -1


def accumulate_panel_presence(annotated_df: pd.DataFrame,
                              distinctive: Dict[str, Set[str]],
                              acc: Dict[str, Dict[str, int]]) -> None:
    """Update `acc` (variant -> {"present":int, "co_covered":int}) from one
    annotated batch (output of annotate_cooc_dataframe). Per read row, weighted
    by its `count`:
      present    += count if >=2 of the variant's distinctive mutations appear
                    TOGETHER in confirmed_present (the defining co-occurrence,
                    as the mutant base) — the actual evidence the variant is here.
      co_covered += count if >=2 of the variant's distinctive POSITIONS were
                    sequenced together (present in confirmed_present OR
                    confirmed_absent — non-N, any base) — "did we get to look".
    Pure; no IO. Accumulate across every batch/date, then decide the verdict
    from the totals (freq = present / co_covered)."""
    if annotated_df is None or getattr(annotated_df, "empty", True):
        return
    dpos = {v: {_presence_pos(m) for m in muts if _presence_pos(m) >= 0}
            for v, muts in distinctive.items()}
    pos2mut = {v: {_presence_pos(m): m for m in muts if _presence_pos(m) >= 0}
               for v, muts in distinctive.items()}
    for row in annotated_df.to_dict("records"):
        cnt = int(row.get("count", 0) or 0)
        if cnt <= 0:
            continue
        cp = set(row.get("confirmed_present", []) or [])
        ca = set(row.get("confirmed_absent", []) or [])
        covered_pos = ({_presence_pos(m) for m in cp}
                       | {_presence_pos(m) for m in ca})
        for v, muts in distinctive.items():
            a = acc.get(v)
            if a is None:
                a = acc[v] = {"present": 0, "co_covered": 0,
                              "mut_cov": {}, "mut_pres": {}}
            dcov = dpos[v] & covered_pos
            if len(dcov) >= 2:
                a["co_covered"] += cnt
            if len(muts & cp) >= 2:
                a["present"] += cnt
            # per-mutation constellation: cover / present for each distinctive
            # position sequenced on this read (works even when they never
            # co-occur, i.e. co_covered stays 0)
            _mc = a["mut_cov"]; _mp = a["mut_pres"]
            for _p in dcov:
                _m = pos2mut[v][_p]
                _mc[_m] = _mc.get(_m, 0) + cnt
                if _m in cp:
                    _mp[_m] = _mp.get(_m, 0) + cnt


def constellation_counts(mut_cov, mut_pres, cov_min: int = 3000,
                         freq_min: float = 0.01):
    """(#present, #testable) over a variant's distinctive mutations.

    testable = distinctive mutations with >= cov_min covering reads (readable);
    present  = of those, mutations whose per-mutation frequency
               (mut_pres / mut_cov) is >= freq_min. A variant that is genuinely
               present shows most of its constellation; an extinct one shows
               only stray homoplastic sites."""
    mut_cov = mut_cov or {}
    mut_pres = mut_pres or {}
    testable = [m for m, c in mut_cov.items() if c >= cov_min]
    present = [m for m in testable
               if (mut_pres.get(m, 0) / mut_cov[m]) >= freq_min]
    return len(present), len(testable)