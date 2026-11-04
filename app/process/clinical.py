"""Clinical sequences per week, for the details' week columns (2026-10-27).

One request per region for the chosen dates: every sequence's collection date
and Nextclade Pango lineage, counted per week (Monday). A lineage's count
includes its descendants in the Pango tree the app uses, so aliased children
(TK.1 under RV.2.1.1) count with their parent's family.

Open CoV-Spectrum LAPIS (no key): Europe was 213 sequences in Jun-Sep 2026,
171 of them from the UK and 42 from Germany, so these are thin and late (the
last 2-3 weeks fill in as sequences are uploaded). They are shown as counts,
never used for a decision.
"""
import json
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import date, timedelta
from typing import Dict, Iterable, Optional

LAPIS = "https://lapis.cov-spectrum.org/open/v2/sample/aggregated"


def _url():
    try:
        from utils.config import get_cooc_setting
        return get_cooc_setting("clinical.lapis_url", default=LAPIS) or LAPIS
    except Exception:
        return LAPIS


def _monday(d) -> str:
    x = date.fromisoformat(str(d)[:10])
    return (x - timedelta(days=x.weekday())).isoformat()


def fetch_weekly(d0: str, d1: str, region: Optional[str] = None,
                 timeout: int = 60) -> Optional[Dict[str, Dict[str, int]]]:
    """{monday: {lineage: n, "": total}} for the dates, or None if LAPIS
    can't be reached. region None = the world."""
    q = {"dateFrom": str(d0)[:10], "dateTo": str(d1)[:10],
         "fields": "date,nextcladePangoLineage"}
    if region:
        q["region"] = region
    try:
        with urllib.request.urlopen(_url() + "?" + urllib.parse.urlencode(q),
                                    timeout=timeout) as r:
            rows = json.loads(r.read()).get("data") or []
    except Exception:
        return None
    out: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for x in rows:
        if not x.get("date"):
            continue
        wk, n = _monday(x["date"]), int(x.get("count") or 0)
        out[wk][""] += n
        if x.get("nextcladePangoLineage"):
            out[wk][x["nextcladePangoLineage"]] += n
    return {k: dict(v) for k, v in out.items()}


def family(lineage: str, parent: Dict[str, str]) -> set:
    """lineage and all its descendants in the Pango tree"""
    kids = defaultdict(list)
    for c, p in (parent or {}).items():
        if p:
            kids[p].append(c)
    out, stack = {lineage}, [lineage]
    while stack:
        for c in kids.get(stack.pop(), []):
            if c not in out:
                out.add(c)
                stack.append(c)
    return out


def lineage_weeks(weekly: Optional[Dict[str, Dict[str, int]]], members: Iterable[str],
                  weeks: Iterable[str]) -> Optional[Dict[str, tuple]]:
    """{monday: (sequences of the lineage family, all sequences)} for the weeks"""
    if weekly is None:
        return None
    m = set(members)
    return {w: (sum(n for k, n in (weekly.get(w) or {}).items() if k in m),
                (weekly.get(w) or {}).get("", 0)) for w in weeks}


# ── who carries a pattern (2026-11-04) ─────────────────────────────────────
# A recurring unassigned pattern is a combination no lineage's signature has.
# Clinical sequences can still say whose it is: 25409A 25455T 25528T was ST.1
# in 90 of 90 sequences (Aug–Oct 2026), a branch of ST.1 with no name yet
# (69 % of ST.1). Shown as words and an extra tree row, never used for a
# decision.
BRANCH_DEFAULTS = {"branch_min_sequences": 20, "branch_min_share": 0.9}


def branch_cfg(cfg: Optional[dict] = None) -> dict:
    out = dict(BRANCH_DEFAULTS)
    try:
        from utils.config import get_cooc_setting
        for k, v in BRANCH_DEFAULTS.items():
            out[k] = type(v)(get_cooc_setting(f"clinical.{k}", v))
    except Exception:
        pass
    out.update(cfg or {})
    return out


def fetch_carriers(mutations: Iterable[str], d0: str, d1: str, region: Optional[str] = None,
                   timeout: int = 15) -> Optional[Dict[str, int]]:
    """{lineage: n} of the clinical sequences with all the mutations in the
    dates ({} = none), or None if LAPIS can't be reached"""
    q = {"nucleotideMutations": ",".join(mutations), "dateFrom": str(d0)[:10],
         "dateTo": str(d1)[:10], "fields": "nextcladePangoLineage"}
    if region:
        q["region"] = region
    try:
        with urllib.request.urlopen(_url() + "?" + urllib.parse.urlencode(q),
                                    timeout=timeout) as r:
            rows = json.loads(r.read()).get("data") or []
    except Exception:
        return None
    out: Dict[str, int] = defaultdict(int)
    for x in rows:
        out[x.get("nextcladePangoLineage") or ""] += int(x.get("count") or 0)
    return dict(out)


def carrier(counts: Dict[str, int], parent: Dict[str, str], nodes: Iterable[str],
            cfg: Optional[dict] = None):
    """(node, n in its family, total): the deepest of `nodes` (the tree's rows)
    whose family holds >= branch_min_share of >= branch_min_sequences
    sequences with the pattern; None otherwise"""
    c = branch_cfg(cfg)
    total = sum((counts or {}).values())
    if total < int(c["branch_min_sequences"]):
        return None
    nodes = set(nodes)
    fam: Dict[str, int] = defaultdict(int)
    for lin, n in counts.items():
        x, seen = lin, set()
        while x and x not in seen:
            seen.add(x)
            if x in nodes:
                fam[x] += n
            x = (parent or {}).get(x, "")

    def depth(x):
        d, seen = 0, set()
        while x and x not in seen:
            seen.add(x)
            x = (parent or {}).get(x, "")
            d += 1
        return d
    ok = [x for x, n in fam.items() if n >= float(c["branch_min_share"]) * total]
    if not ok:
        return None
    best = max(ok, key=lambda x: (depth(x), -fam[x]))
    return best, fam[best], total


def carrier_words(counts: Optional[Dict[str, int]], parent: Dict[str, str], nodes: Iterable[str],
                  cfg: Optional[dict] = None):
    """(words, node or None) for a pattern's clinical carriers: 'clinical: ST.1
    (55 of 55 sequences)', 'clinical: mixed · most XFG (8 of 30)', 'clinical:
    no sequences with all of them', 'clinical: 7 sequences, too few (ST.1 7)';
    ("", None) when LAPIS couldn't be reached"""
    if counts is None:
        return "", None
    c = branch_cfg(cfg)
    total = sum(counts.values())
    if not total:
        return "clinical: no sequences with all of them", None
    top, k = max(counts.items(), key=lambda kv: kv[1])
    top = top or "no lineage"
    hit = carrier(counts, parent, nodes, c)
    if hit:
        x, n, t = hit
        return f"clinical: {x} ({n:,} of {t:,} {'sequence' if t == 1 else 'sequences'})", x
    if total < int(c["branch_min_sequences"]):
        return (f"clinical: {total:,} {'sequence' if total == 1 else 'sequences'}, too few "
                f"(most {top}, {k:,})"), None
    return f"clinical: mixed · most {top} ({k:,} of {total:,})", None
