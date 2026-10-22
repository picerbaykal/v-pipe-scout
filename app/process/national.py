"""One lineage, all cities pooled (2026-10-22): the national answer and one
square per week, shown on top of the per-city sample squares in a row's
details.

Per city most summer samples have low coverage (Jun–Sep 2026: 20 of 138
BA.3.2 samples well covered; Chur none), so a city's squares rest on thin
reads. Pooled over the cities, every week but one gave an answer and the
trend shows: BA.3.2 seen every week to 20 Apr, not seen from 1 Jun; TA.1 not
seen until the week of 3 Aug, seen from the week of 10 Aug (7:0, 4:1, 5:2,
3:0). In winter, with well-covered samples, the city squares carry the
detail (TA.1 reached Genève before Zürich) and this line is the summary.

Input per city: process.variant_explorer.check_in_data's result, which
carries "samples" = {date: {marker: [cov, hit, link_n, link_ok]}} — the
COUNTED reads (each sample at most min(cap_reads, cov_factor × its median
depth) per marker, shares unchanged) — and "quality".

Every answer here is the cells' own rule (process.cooc.check_verdicts) on
the counted reads it pools: no new number, no new setting. Pure, no IO."""
from datetime import date, timedelta
from typing import Dict, List, Optional


def _monday(d: str) -> date:
    x = date.fromisoformat(str(d)[:10])
    return x - timedelta(days=x.weekday())


def _vote(markers: List[str], samples: List[dict]) -> dict:
    """check_verdicts on already-counted reads: keys are unique per sample
    (date#city) and no second cap is applied."""
    from process.cooc import check_verdicts
    return check_verdicts(markers, {f"{s['date']}#{s['city']}": s["cells"] for s in samples},
                          dict(cap_reads=0))


def national(per_city: Dict[str, dict]) -> Optional[dict]:
    """The pooled answer for the chosen dates and one per week, or None
    without samples.

    Returns {answer: check_verdicts(...) for all samples, n_cities,
             weeks: [{start (Monday, ISO), answer, n_samples, n_well,
                      cities: {city: verdict of that city's samples that week}}]}"""
    try:
        from process.sample_qc import well_covered
    except Exception:                                   # pragma: no cover
        well_covered = lambda _q: False                 # noqa: E731
    markers, samples = [], []
    for city, d in (per_city or {}).items():
        if not d:
            continue
        for m in d.get("markers") or []:
            if m not in markers:
                markers.append(m)
        q = d.get("quality") or {}
        for dt, cells in (d.get("samples") or {}).items():
            day = str(dt)[:10]
            samples.append({"date": day, "city": city, "cells": cells,
                            "well": bool(well_covered(q.get(day)))})
    if not markers or not samples:
        return None
    by_week: Dict[date, List[dict]] = {}
    for s in samples:
        by_week.setdefault(_monday(s["date"]), []).append(s)
    weeks = []
    for wk in sorted(by_week):
        ss = by_week[wk]
        cities = {}
        for c in sorted({s["city"] for s in ss}):
            cities[c] = _vote(markers, [s for s in ss if s["city"] == c])["verdict"]
        weeks.append({"start": wk.isoformat(), "answer": _vote(markers, ss),
                      "n_samples": len(ss), "n_well": sum(s["well"] for s in ss),
                      "cities": cities})
    return {"answer": _vote(markers, samples),
            "n_cities": len({s["city"] for s in samples}), "weeks": weeks}
