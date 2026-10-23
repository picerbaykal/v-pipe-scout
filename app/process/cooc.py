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
    in_tree=None,
    support: Optional[Dict[str, float]] = None,
    stats: Optional[Dict[str, int]] = None,
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
    Independent of which other variants are in the panel.

    Several variants one change away: the one sharing more mutations with the
    read. Still a tie (2026-10-09; was: first in panel order, so the label
    changed with the panel's order):
      1. pango decides: in_tree(variant, sign, mut) -> True when a designated
         sublineage of that variant has exactly this change (where_in_tree);
         if exactly one tied variant has it, it is that one;
      2. else the dominant variant: support = {variant: reads it explains
         that day}, the most supported wins;
      3. else alphabetical (stable, but arbitrary: counted as unresolved).
    stats, if given, counts {"ties", "pango", "dominant", "unresolved"}."""
    cands = []
    for name, sig in variant_signatures.items():
        shared = len(confirmed_present & sig)
        if shared < 2:
            continue
        extra = confirmed_present - sig
        missing = sig & confirmed_absent
        if len(extra) + len(missing) == 1:
            if extra:
                cands.append((name, "+", next(iter(extra)), shared))
            else:
                cands.append((name, "−", next(iter(missing)), shared))
    if not cands:
        return ""
    top = max(c[3] for c in cands)
    tied = sorted((c for c in cands if c[3] == top), key=lambda c: c[0])
    best = tied[0]
    if len(tied) > 1:
        how = "unresolved"
        if in_tree is not None:
            known = [c for c in tied if in_tree(c[0], c[1], c[2])]
            if len(known) == 1:
                best, how = known[0], "pango"
            elif known:
                tied = known
        if how == "unresolved" and support:
            def _sup(c):
                return float(support.get(c[0], 0) or 0)
            sup = sorted(tied, key=lambda c: (-_sup(c), c[0]))
            if _sup(sup[0]) > _sup(sup[1]):
                best, how = sup[0], "dominant"
        if how == "unresolved":
            best = tied[0]
        if stats is not None:
            stats["ties"] = stats.get("ties", 0) + 1
            stats[how] = stats.get(how, 0) + 1
    return f"{best[0]} {best[1]} {best[2]}"


# ── Per-batch DataFrame processing ──────────────────────────────────────

def annotate_cooc_dataframe(
    df: pd.DataFrame,
    positions: List[int],
    amp_dict: Dict[int, List[str]],
    variant_signatures: Dict[str, Set[str]],
    in_tree=None,
    stats: Optional[Dict[str, int]] = None,
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
    parsed = []
    # reads each panel variant explains, per date: the "dominant variant"
    # tie-breaker of near_panel_label (2026-10-09)
    support: Dict[str, Dict[str, float]] = {}
    for row in df.to_dict("records"):
        if require_covered:
            if any(row.get(f"[{p}]", "N") in uncovered_bases for p in positions):
                dropped += 1
                continue
        cp, ca = row_to_confirmed_sets(row, positions, amp_dict, include_dels)
        classification = classify_pattern(cp, ca, variant_signatures)
        if classification == "matched":
            day = support.setdefault(str(row["date"]), {})
            for name, sig in variant_signatures.items():
                if cp.issubset(sig) and not (sig & ca):
                    day[name] = day.get(name, 0) + float(row["count"] or 0)
        parsed.append((row, cp, ca, classification))
    for row, cp, ca, classification in parsed:
        near = (near_panel_label(cp, ca, variant_signatures, in_tree=in_tree,
                                 support=support.get(str(row["date"])), stats=stats)
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
    # the vote is a simple majority of the measured markers (2026-10-20; until
    # then present >= 75 % and absent <= 25 %, `confirm_share` / `notfound_share`):
    # more present than absent -> present, more absent -> absent, a tie -> mixed.
    # Same false alarms, one more circulating lineage found, two numbers fewer
    # (bench_params 16 Oct; bench_scorecard 20 Oct, 6 cities, 2 windows)
    "min_present": 2,        # present also needs at least this many markers
                             #     present (one marker can come from an unrelated
                             #     lineage that carries it: XFZ carries LP.8's
                             #     1954A; RV.1.3's 15982A at 52 % on one Zürich
                             #     sample while its 4444A was at 1 %)
    "calendar_samples": 12,  # Evidence squares shown per row (UI only)
    "cap_reads": 1000,       # a sample counts at most this many reads per marker
                             #   (2026-10-13): in low-virus samples a few amplicons
                             #   reach 100,000s of reads from few molecules and
                             #   would outvote every other sample
    "cov_factor": 10,        # ... and at most this × its median depth (2026-10-16):
                             #   every sample votes, weighted by its coverage. A
                             #   sample with median depth 2 and 400,000 reads at one
                             #   marker (a jackpot: few virus copies, copied many
                             #   times) counts as 20 reads; a broadly covered one as
                             #   cap_reads. ×5 and ×20 gave the same verdicts (6
                             #   cities, 3 windows: bench_samples, 16 Oct 2026)
    "own_marker": True,      # present also needs one of the mutations the lineage
                             #   gained itself (not its parent's) to be present,
                             #   when one of them is measured not seen (2026-10-23).
                             #   A child carries its parent's mutations: SB.1.1
                             #   2:1 in Zürich, but both seen markers are SB.1's
                             #   and its own 1501A is at 0.8 % of 2,587 reads
    # A found lineage, graded on all cities pooled (process.national.grade,
    # 2026-10-23; bench_scorecard 6 cities, summer and spring):
    "find_pooled_freq": 0.025,   # found also when all cities pooled say seen at this
                                 #   share, though no city alone does (summer: +2
                                 #   circulating lineages, no false alarm)
    "confirm_pooled_freq": 0.05,  # confirmed: all cities pooled say seen at this share,
    "confirm_min_cities": 2,     #   ... and it's seen in at least this many cities
    "confirm_min_weeks": 2,      #   ... or in this many weeks of the all-cities line,
    "confirm_max_missing": 0.10,  #   ... and at most this share of its measured
                                 #   mutations is clearly missing (< absent_freq) pooled
                                 #   (PJ.2.1, Sep–Oct: 4 of 18 at 0.0–0.2 %). Confirmed
                                 #   held no gone or not-in-Europe lineage in either
                                 #   window; "early" (found, not confirmed): summer 4
                                 #   of 5 real
    "thin_well_min": 2,          # a city with fewer well-covered samples: a seen cell
                                 #   is striped (display only; counting thin cities
                                 #   or not gave the same grades, 23 Oct)
}


def sample_cap(c: dict, depth: Optional[Dict[str, float]], date: str) -> float:
    """Reads one sample counts at most per marker: check.cap_reads, and at
    most check.cov_factor × the sample's median depth when it is known
    (never under 1). 0 = no cap."""
    cap = float(c.get("cap_reads") or 0)
    if depth is None:
        return cap
    md = depth.get(str(date)[:10])
    if md is None:
        return cap
    lim = max(1.0, float(c.get("cov_factor") or 0) * float(md))
    return min(cap, lim) if cap else lim


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


def own_markers(variant: str, markers: Optional[List[str]],
                lineage_signatures: Dict[str, Set[str]],
                parent_map: Dict[str, str]) -> Optional[dict]:
    """The mutations `variant` gained itself — not carried by its parent in
    the tree (2026-10-23): {"parent": name, "markers": [...]}, kept to
    `markers` when given (check_verdicts keeps them to the lineage's markers
    anyway), or None when it has no parent in the tree (a recombinant root:
    every marker is its own). SB.1.1: parent SB.1; of its markers 1501A,
    21414C (SB.1's) and 18585T (XFG.6.2.2's) only 1501A is its own."""
    p = (parent_map or {}).get(variant) or ""
    if not p or p not in (lineage_signatures or {}):
        return None
    gained = set(lineage_signatures.get(variant) or ()) - set(lineage_signatures[p] or ())
    keep = sorted(gained) if markers is None else [m for m in markers if m in gained]
    return {"parent": p, "markers": keep}


def check_verdicts(markers: List[str],
                   per_date: Dict[str, Dict[str, List[int]]],
                   cfg: Optional[dict] = None,
                   no_link: Optional[set] = None,
                   depth: Optional[Dict[str, float]] = None,
                   own: Optional[dict] = None) -> dict:
    """Verdict for one variant in one city: a vote among its specific markers.

    Counts are pooled over all dates in the window (the heatmaps show dates),
    each date counting at most cap_reads reads per marker (2026-10-13), and
    at most cov_factor × its median depth when depth = {date: median depth}
    is given (sample_cap, 2026-10-16). Shares are kept, counts scaled.
    Each marker is
      present     >= present_freq of covering reads carry it, and those reads
                  match the variant at neighbouring positions (link >= link_min)
      absent      <  absent_freq of covering reads carry it
      unmeasured  fewer than min_cov covering reads, or anything in between
    A marker with no neighbour in reach (no_link: at least min_link reads carry
    it, and not one of them covers another position of the variant — its
    nearest mutation is further than a read, e.g. NB.1.8.1's 2095A, 695 bases)
    can't take the link test; it is judged on its share alone (2026-10-13;
    before, it could never be present). Computed from the reads unless given
    (check_in_data passes the window's set to each day).
    Verdict: a majority vote of the measured markers (2026-10-20)
      present  more markers present than absent, AND at least min_present
               present (one marker alone can come from an unrelated
               lineage, 2026-10-09)
      absent   more markers absent than present
      mixed    a tie, or more present than absent but fewer than min_present
      not_covered when the variant has no markers or none is measurable.
    (Until 2026-10-20: present >= 75 %, absent <= 25 %, mixed in between.)
    (Until 2026-10-08 the verdicts were named confirmed / not_found /
    inconsistent / cant_confirm.)
    own (own_markers(), 2026-10-23, setting check.own_marker): a present vote
    becomes mixed ("own_missing") when none of the lineage's own markers is
    present and at least one is absent — every seen marker is one it shares
    with its parent, and its own say not seen. Own markers with too few reads
    leave the vote as it is.
    Uses reads only — never the deconvolution abundance."""
    c = _check_cfg(cfg)
    empty = {"verdict": "not_covered", "n_markers": len(markers or []),
             "n_present": 0, "n_measured": 0, "markers": {}}
    if not markers:
        return empty
    pooled = {m: [0.0, 0.0, 0.0, 0.0] for m in markers}
    raw = {m: 0 for m in markers}
    for dt, cells in (per_date or {}).items():
        cap = sample_cap(c, depth, dt)
        for m in markers:
            cell = cells.get(m)
            if cell:
                raw[m] += cell[0]
                # each sample counts at most `cap` reads, its share unchanged
                k = cap / cell[0] if cap and cell[0] > cap else 1.0
                for i in range(4):
                    pooled[m][i] += cell[i] * k
    if no_link is None:
        no_link = {m for m in markers
                   if pooled[m][2] == 0 and pooled[m][1] >= c["min_link"]}
    out, n_p, n_a = {}, 0, 0
    for m in markers:
        cov, hit, lk_n, lk_ok = pooled[m]
        f = hit / cov if cov else None
        link = lk_ok / lk_n if lk_n else None
        if cov < c["min_cov"]:
            st = "unmeasured"
        elif f < c["absent_freq"]:
            st = "absent"
        elif f >= c["present_freq"] and (m in no_link or (lk_n >= c["min_link"]
                                                          and link >= c["link_min"])):
            st = "present"
        else:
            st = "unmeasured"
        n_p += st == "present"
        n_a += st == "absent"
        out[m] = {"cov": int(round(cov)), "cov_raw": raw[m], "freq": f, "link": link,
                  "status": st, "no_link": m in no_link}
    n_meas = n_p + n_a
    if n_meas == 0:
        verdict = "not_covered"
    elif n_p > n_a and n_p >= int(c["min_present"]):
        verdict = "present"
    elif n_a > n_p:
        verdict = "absent"
    else:
        verdict = "mixed"
    own_l = [m for m in ((own or {}).get("markers") or []) if m in out]
    own_missing = bool(
        verdict == "present" and c.get("own_marker", True) and own_l
        and not any(out[m]["status"] == "present" for m in own_l)
        and any(out[m]["status"] == "absent" for m in own_l))
    if own_missing:
        verdict = "mixed"
    return {"verdict": verdict, "n_markers": len(markers), "n_present": n_p,
            "n_measured": n_meas, "markers": out, "no_link": sorted(no_link),
            "own": own_l if own else None,
            "parent": (own or {}).get("parent"), "own_missing": own_missing}

# ══ Lineages found below a found lineage: how much of it is itself ═════════
# (2026-10-14) An ancestor's ★ markers are carried by all its descendants, so
# its ★ check measures the whole family (BA.2.86 "present" in Zürich because
# JN.1 and everything after it carry its markers). Per city, for each found
# lineage X, compare it with the nearest lineages below it that are present
# and will be in the panel (in the panel already, or found and offered + Add):
#
#     Jaccard(X, D) = reads with both / reads with either
#                   = share of D's family ÷ share of X's markers D carries
#     X itself      = 1 − Σ Jaccard over those lineages
#
# Every read of D's family also has X's markers (D descends from X), so
# "both" = D's family and "either" = X's. Each share is read at its own
# positions; the two need not be on one read. (A first version counted both
# on the same reads; in Zürich a marker of X lay within one read of a
# mutation of D in 6 of 37 cases only.)
#
# A family's share is read only at markers that every lineage below it with
# reads in the city still carries: descendants lose markers again (JN.1's
# 3565C is on 58 % of reads in every city, 27810C on 100 %), and the median of
# all markers would leave a false "itself" (BA.2.86 22 % instead of 0 %).
#
# "Below" comes from the Pango tree: X is an ancestor of D. A recombinant (and
# its descendants) is also below a lineage outside its own line whose ★
# markers it carries: XFG below LF.7 and LP.8.1.


def _ancestors(lin: str, parent_map: Dict[str, str]) -> List[str]:
    out, x, seen = [], parent_map.get(lin, ""), {lin}
    while x and x not in seen:
        out.append(x)
        seen.add(x)
        x = parent_map.get(x, "")
    return out


def lineages_below_all(nodes: List[str], markers: Dict[str, List[str]],
                       sigs: Dict[str, Set[str]],
                       parent_map: Dict[str, str]) -> Dict[str, List[str]]:
    """{X: [every lineage of `nodes` below X]}. Pure.

    D is below X when X is a Pango ancestor of D, or when D is (or descends
    from) a recombinant, X is not in its line, and D carries at least one of
    X's ★ markers. Two lineages each "below" the other (two recombinants
    sharing markers) are left unrelated."""
    nodes = list(dict.fromkeys(nodes))
    anc = {n: set(_ancestors(n, parent_map)) for n in nodes}
    rec = {n: _recombinant_root(n, parent_map) is not None for n in nodes}

    def is_below(d, x):
        if d == x:
            return False
        if x in anc[d]:
            return True
        sd = sigs.get(d) or set()
        return (rec[d] and d not in anc[x]
                and any(m in sd for m in markers.get(x) or []))

    below = {x: [d for d in nodes if is_below(d, x)] for x in nodes}
    return {x: [d for d in ds if x not in below.get(d, [])] for x, ds in below.items()}


def nearest_below(x: str, below_all: Dict[str, List[str]], take) -> List[str]:
    """The lineages below X for which take(d) is true and with no other such
    lineage between X and them. Pure."""
    cand = [d for d in below_all.get(x) or [] if take(d)]
    return [d for d in cand
            if not any(o != d and d in (below_all.get(o) or []) for o in cand)]


def kept_markers(mks: List[str], d: str, below_all: Dict[str, List[str]],
                 sigs: Dict[str, Set[str]], has_reads) -> List[str]:
    """The markers of `mks` that D and every lineage below D with
    has_reads(e) still carry (by their signatures). Pure."""
    keep = [m for m in mks or [] if m in (sigs.get(d) or set())]
    for e in below_all.get(d) or []:
        if has_reads(e):
            se = sigs.get(e) or set()
            keep = [m for m in keep if m in se]
    return keep


def lineage_share(markers: List[str], per_date: Dict[str, Dict[str, List[int]]],
                  dates=None, cfg: Optional[dict] = None,
                  depth: Optional[Dict[str, float]] = None) -> Optional[float]:
    """Share of reads carrying a lineage's markers in one city: the median,
    over the markers with ≥ check.min_cov reads, of the share of reads at the
    marker that carry it. Reads pooled over `dates` (all when None), each date
    counting at most sample_cap reads per marker, as the ★ check does.
    None when no marker has enough reads."""
    c = _check_cfg(cfg)
    vals = []
    for m in markers or []:
        cov = hit = 0.0
        for dt, cells in (per_date or {}).items():
            if dates is not None and dt not in dates:
                continue
            cell = (cells or {}).get(m)
            if not cell:
                continue
            cap = sample_cap(c, depth, dt)
            k = cap / cell[0] if cap and cell[0] > cap else 1.0
            cov += cell[0] * k
            hit += cell[1] * k
        if cov >= c["min_cov"]:
            vals.append(hit / cov)
    if not vals:
        return None
    vals.sort()
    n = len(vals)
    return vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2


def family_jaccard(share_x: Optional[float], share_d: Optional[float]) -> Optional[float]:
    """Jaccard(X, D) = share of D ÷ share of X (at the markers D carries),
    at most 1. None when either share can't be measured or X's is 0."""
    if share_x is None or share_d is None or share_x <= 0:
        return None
    return min(1.0, share_d / share_x)
