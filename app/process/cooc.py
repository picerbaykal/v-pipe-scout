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


def near_panel_label(
    confirmed_present: Set[str],
    confirmed_absent: Set[str],
    variant_signatures: Dict[str, Set[str]],
) -> str:
    """For an UNEXPLAINED read: "variant + 1 change" label, else "".

    Compared with each panel variant SEPARATELY (not against the pooled panel
    mutations), the read differs from the best-matching variant at exactly one
    position:
      extra   a mutation on the read the variant doesn't have (a new mutation,
              or another base where the variant has its own) -> "XFG + 22896C"
      missing a mutation of the variant, read shows reference there
                                                        -> "XFG − 23021G"
    Needs >= 2 of the variant's mutations on the read so the match means
    something. Reads differing at 2+ positions from every variant return "".
    Independent of which other variants are in the panel."""
    best = None
    for name, sig in variant_signatures.items():
        if len(confirmed_present & sig) < 2:
            continue
        extra = confirmed_present - sig
        missing = sig & confirmed_absent
        n = len(extra) + len(missing)
        if n == 1 and best is None:
            best = (name, extra, missing)
        elif n == 1 and best is not None:
            # tie between variants: keep the one sharing more mutations
            if len(confirmed_present & sig) > len(confirmed_present &
                                                   variant_signatures[best[0]]):
                best = (name, extra, missing)
    if best is None:
        return ""
    name, extra, missing = best
    if extra:
        return f"{name} + {next(iter(extra))}"
    return f"{name} − {next(iter(missing))}"


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
        near = (near_panel_label(cp, ca, variant_signatures)
                if classification == "unexplained" else "")
        annotated_rows.append({
            "date": row["date"],
            "count": row["count"],
            "confirmed_present": sorted(cp),
            "confirmed_absent": sorted(ca),
            "classification": classification,
            "near": near,
        })

    if dropped:
        logger.info(
            f"coverage filter dropped {dropped}/{len(df)} rows "
            f"({dropped / len(df):.1%}) for positions {positions}"
        )

    return pd.DataFrame(annotated_rows)


# ── Unexplained reads handed to the scanner ─────────────────────────────

def aggregate_unexplained(frames: List[pd.DataFrame]) -> pd.DataFrame:
    """Merge the unexplained rows of all batches into one row per
    (date, present mutations, absent mutations), with the read count.

    `confirmed_absent` (positions where the read shows the reference base) is
    kept: reads with the same present mutations but different absent ones are
    different evidence — a read showing the reference where a lineage has its
    mutation cannot come from that lineage. (Until 2026-09-30 rows were keyed
    by present mutations only and the absent side was dropped.)

    Also kept: the "panel variant + 1 change" summary per row — near_count
    (reads with a label) and near_label (the most common label).

    Returns columns: date, count, confirmed_present, confirmed_absent,
    near_count, near_label (lists for the mutation columns)."""
    cols = ["date", "count", "confirmed_present", "confirmed_absent",
            "near_count", "near_label"]
    frames = [f for f in (frames or []) if f is not None and not f.empty]
    if not frames:
        return pd.DataFrame(columns=cols)
    df = pd.concat(frames, ignore_index=True)
    if "confirmed_absent" not in df.columns:
        df["confirmed_absent"] = [[] for _ in range(len(df))]
    if "near" not in df.columns:
        df["near"] = ""
    df["near"] = df["near"].fillna("")
    df["_p"] = df["confirmed_present"].apply(lambda x: tuple(sorted(x)))
    df["_a"] = df["confirmed_absent"].apply(lambda x: tuple(sorted(x or ())))
    df["near_count"] = df["count"].where(df["near"] != "", 0)
    keys = ["date", "_p", "_a"]
    agg = (df.groupby(keys).agg(count=("count", "sum"), near_count=("near_count", "sum"))
             .reset_index())
    lab = (df[df["near"] != ""].groupby(keys + ["near"])["count"].sum().reset_index()
             .sort_values("count", ascending=False).drop_duplicates(keys)
             .rename(columns={"near": "near_label"})[keys + ["near_label"]])
    agg = agg.merge(lab, on=keys, how="left")
    agg["near_label"] = agg["near_label"].fillna("")
    agg["near_count"] = agg["near_count"].astype(int)
    agg["count"] = agg["count"].astype(int)
    agg["confirmed_present"] = agg["_p"].apply(list)
    agg["confirmed_absent"] = agg["_a"].apply(list)
    return agg[cols]


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
import re as _re_presence


# ══ New co-occurrence check (2026-09): specific markers + neighbour linkage ══
#
# Part A — which markers: panel-private substitutions whose carriers OUTSIDE the
#   variant's own family are few. Own sublineages never count against a variant
#   (that is what broke the global <=30 rule on the granular Nextclade tree).
#   Outside carriers count, except those under a recombinant that INHERITED the
#   marker from this variant (XFV = LP.8.1 × XFG.3.3.1 for XFG's 8350C) — see
#   process.recombinants (2026-10-02, "option C"; before, every carrier under any
#   other recombinant root was tolerated up to 60, including unrelated circulating
#   ones such as SV.4 for B.1.1.7's 3267T).
# Part B — reads: per date, each marker's coverage/frequency, plus "link": of
#   reads carrying the marker that also cover other positions of the variant's
#   signature, the share where every covered neighbour carries the variant base
#   (marker sits on a variant-like haplotype, not an isolated error). Markers
#   are NOT required to co-occur with each other (often ~kb apart).
# Verdicts are computed from raw counts by check_verdicts() with thresholds from
#   cooc_config.yaml `check:`, so retuning needs no worker re-scan.

_REC_ROOT_RE = _re_presence.compile(r"^X[A-Z]+$")
_SUB_RE = _re_presence.compile(r"^(\d+)([ACGT])$")
_CHECK_UNCOVERED = {"N", "-"}

CHECK_DEFAULTS = {
    # Part A — marker selection (worker; changing needs a re-scan)
    "out_other_max": 5,      # carriers outside the family that didn't inherit it
    # Part B — the vote (UI; a streamlit restart is enough)
    "min_cov": 100,          # reads covering a marker (whole window) to measure it
    "present_freq": 0.05,    # present: >= this share of covering reads carry it ...
    "link_min": 0.8,         # ... and >= this share of those reads match the
    "min_link": 20,          #     variant at neighbouring positions (>= 20 reads)
    "absent_freq": 0.01,     # absent: < this share carry it
    "confirm_share": 0.75,   # present / measurable >= this -> confirmed
    "notfound_share": 0.25,  # present / measurable <= this -> not found
}


# "What makes a mutation a ★ marker" is shared with the scanner and lives under
# `markers:` in cooc_config.yaml; `check.*` is still read as a fallback so older
# config files keep working.
_MARKER_KEYS = ("out_other_max",)


def _check_cfg(cfg: Optional[dict] = None) -> dict:
    out = {}
    for k, v in CHECK_DEFAULTS.items():
        if cfg and k in cfg:
            out[k] = cfg[k]
        elif k in _MARKER_KEYS:
            out[k] = get_cooc_setting(f"markers.{k}",
                                      get_cooc_setting(f"check.{k}", v))
        else:
            out[k] = get_cooc_setting(f"check.{k}", v)
    return out


def _children_map(parent_map: Dict[str, str]) -> Dict[str, List[str]]:
    kids: Dict[str, List[str]] = {}
    for c, p in parent_map.items():
        if p:
            kids.setdefault(p, []).append(c)
    return kids


def _lineage_family(v: str, parent_map: Dict[str, str],
                    kids: Dict[str, List[str]]) -> Set[str]:
    """v plus all descendants. A panel variant absent from the file (built by the
    children-intersection fallback) is seeded with its name-prefix children."""
    roots = {v}
    if v not in parent_map:
        roots |= {l for l in parent_map if l.startswith(v + ".")}
    fam, stack = set(roots), list(roots)
    while stack:
        for c in kids.get(stack.pop(), []):
            if c not in fam:
                fam.add(c)
                stack.append(c)
    return fam


def _recombinant_root(lin: str, parent_map: Dict[str, str]) -> Optional[str]:
    """Topmost X.. ancestor-or-self (recombinant roots have an empty parent)."""
    x, found, seen = lin, None, set()
    while x and x not in seen:
        seen.add(x)
        if _REC_ROOT_RE.match(x):
            found = x
        x = parent_map.get(x, "")
    return found


_CARRIER_INDEX: Dict[int, Dict[str, Set[str]]] = {}


def _carrier_index(lineage_signatures: Dict[str, Set[str]]) -> Dict[str, Set[str]]:
    key = id(lineage_signatures)
    if key not in _CARRIER_INDEX:
        idx: Dict[str, Set[str]] = {}
        for lin, sig in lineage_signatures.items():
            for m in sig:
                idx.setdefault(m, set()).add(lin)
        _CARRIER_INDEX.clear()
        _CARRIER_INDEX[key] = idx
    return _CARRIER_INDEX[key]


def specific_markers(variant_signatures: Dict[str, Set[str]],
                     lineage_signatures: Dict[str, Set[str]],
                     parent_map: Dict[str, str],
                     cfg: Optional[dict] = None) -> Dict[str, List[str]]:
    """Part A. variant -> specific markers ("{pos}{alt}"), most specific first.

    variant_signatures: panel variant -> substitutions ("241T" form).
    lineage_signatures: every pango lineage -> substitutions (the carrier pool).
    parent_map:         lineage -> parent lineage ("" for roots/recombinants).

    A variant with no specific marker gets [] and is reported "no ★ marker"
    (KP.2/KP.3: their only own mutations are shared with unrelated variants).
    The variant's own sublineages in the panel don't count as "other" panel
    variants."""
    from process.recombinants import get_recombinants
    c = _check_cfg(cfg)
    idx = _carrier_index(lineage_signatures)
    kids = _children_map(parent_map)
    rc = get_recombinants(parent_map, lineage_signatures)

    out: Dict[str, List[str]] = {}
    for v, sig in variant_signatures.items():
        fam = _lineage_family(v, parent_map, kids)
        # panel-private: not carried by another panel variant — except v's own
        # sublineages in the panel, which carry all of v's mutations by descent
        # (PQ.16.1.1 in the panel left NB.1.8.1 with no marker at all). With
        # both in the panel, v's markers read "v or its sublineages".
        others: Set[str] = set()
        for w, s in variant_signatures.items():
            if w != v and w not in fam:
                others |= (s or set())
        private = {m for m in (sig or set()) - others if _SUB_RE.match(m)}
        own = {rc.rec_root(v)} if v in parent_map else set()
        ranked = []
        for m in private:
            outside = idx.get(m, set()) - fam
            n_inh = sum(1 for l in outside if rc.tolerated(l, v, m, own))
            n_out = len(outside) - n_inh
            if n_out <= c["out_other_max"]:
                ranked.append((n_out, n_inh, int(_SUB_RE.match(m).group(1)), m))
        out[v] = [m for *_, m in sorted(ranked)]
    return out


def _sig_map(sig: Set[str]) -> Dict[int, str]:
    out = {}
    for m in sig or ():
        mm = _SUB_RE.match(m)
        if mm:
            out[int(mm.group(1))] = mm.group(2)
    return out


def accumulate_check_stats(rows: List[dict],
                           positions: List[int],
                           variant_signatures: Dict[str, Set[str]],
                           markers: Dict[str, List[str]],
                           stats: Dict[str, Dict[str, Dict[str, List[int]]]]) -> None:
    """Part B tally from RAW LAPIS rows (one batch of one date).

    stats[variant][date][marker] = [cov, hit, link_n, link_ok]
      cov     reads covering the marker position (non-N, non-deletion)
      hit     of those, reads carrying the marker base
      link_n  of hit, reads that also cover >=1 other signature position
      link_ok of link_n, reads where every covered neighbour has the variant base
    Pure; no IO. Positions outside this batch are ignored."""
    if not rows:
        return
    pos_set = set(positions)
    prepared = []
    for v, mks in markers.items():
        if not mks:
            continue
        sm = {p: b for p, b in _sig_map(variant_signatures.get(v, set())).items()
              if p in pos_set}
        mk = [(int(_SUB_RE.match(m).group(1)), _SUB_RE.match(m).group(2), m)
              for m in mks if _SUB_RE.match(m) and int(_SUB_RE.match(m).group(1)) in pos_set]
        if mk:
            prepared.append((v, sm, mk))
    if not prepared:
        return
    for row in rows:
        cnt = int(row.get("count", 0) or 0)
        if cnt <= 0:
            continue
        date = str(row.get("date", ""))[:10]
        for v, sm, mk in prepared:
            cov_pos = None
            for p, b, m in mk:
                base = row.get(f"[{p}]", "N")
                if base in _CHECK_UNCOVERED:
                    continue
                cell = stats.setdefault(v, {}).setdefault(date, {}).setdefault(m, [0, 0, 0, 0])
                cell[0] += cnt
                if base != b:
                    continue
                cell[1] += cnt
                if cov_pos is None:
                    cov_pos = {q: row.get(f"[{q}]", "N") for q in sm}
                nb = [q for q, bq in cov_pos.items()
                      if q != p and bq not in _CHECK_UNCOVERED]
                if nb:
                    cell[2] += cnt
                    if all(cov_pos[q] == sm[q] for q in nb):
                        cell[3] += cnt


def check_verdicts(markers: List[str],
                   per_date: Dict[str, Dict[str, List[int]]],
                   cfg: Optional[dict] = None) -> dict:
    """Verdict for one variant in one city: a vote among its specific markers.

    Counts are pooled over all dates in the window (the heatmaps show dates).
    Each marker is
      present     >= present_freq of covering reads carry it, and those reads
                  match the variant at neighbouring positions (link >= link_min)
      absent      <  absent_freq of covering reads carry it
      unmeasured  fewer than min_cov covering reads, or anything in between
    Verdict from present / (present + absent):
      confirmed >= confirm_share · not_found <= notfound_share · else inconsistent
      cant_confirm when the variant has no markers or none is measurable.
    Uses reads only — never the deconvolution abundance."""
    c = _check_cfg(cfg)
    empty = {"verdict": "cant_confirm", "n_markers": len(markers or []),
             "n_present": 0, "n_measured": 0, "markers": {}}
    if not markers:
        return empty
    pooled = {m: [0, 0, 0, 0] for m in markers}
    for cells in (per_date or {}).values():
        for m in markers:
            cell = cells.get(m)
            if cell:
                for i in range(4):
                    pooled[m][i] += cell[i]
    out, n_p, n_a = {}, 0, 0
    for m in markers:
        cov, hit, lk_n, lk_ok = pooled[m]
        f = hit / cov if cov else None
        link = lk_ok / lk_n if lk_n else None
        if cov < c["min_cov"]:
            st = "unmeasured"
        elif f < c["absent_freq"]:
            st = "absent"
        elif (f >= c["present_freq"] and lk_n >= c["min_link"]
              and link >= c["link_min"]):
            st = "present"
        else:
            st = "unmeasured"
        n_p += st == "present"
        n_a += st == "absent"
        out[m] = {"cov": cov, "freq": f, "link": link, "status": st}
    n_meas = n_p + n_a
    if n_meas == 0:
        verdict = "cant_confirm"
    elif n_p / n_meas >= c["confirm_share"]:
        verdict = "confirmed"
    elif n_p / n_meas <= c["notfound_share"]:
        verdict = "not_found"
    else:
        verdict = "inconsistent"
    return {"verdict": verdict, "n_markers": len(markers), "n_present": n_p,
            "n_measured": n_meas, "markers": out}