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
