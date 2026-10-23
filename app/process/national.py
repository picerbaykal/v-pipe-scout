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


def _vote(markers: List[str], samples: List[dict], own: Optional[dict] = None,
          present_freq: Optional[float] = None) -> dict:
    """check_verdicts on already-counted reads: keys are unique per sample
    (date#city) and no second cap is applied. own: the lineage's own markers,
    as the cells use them (2026-10-23); present_freq: another "seen" share
    for grade() (else check.present_freq)."""
    from process.cooc import check_verdicts
    cfg = dict(cap_reads=0)
    if present_freq is not None:
        cfg["present_freq"] = present_freq
    return check_verdicts(markers, {f"{s['date']}#{s['city']}": s["cells"] for s in samples},
                          cfg, own=own)


def national(per_city: Dict[str, dict], present_freq: Optional[float] = None) -> Optional[dict]:
    """The pooled answer for the chosen dates and one per week, or None
    without samples. present_freq: the "seen" share (default check.present_freq).

    Returns {answer: check_verdicts(...) for all samples, n_cities,
             weeks: [{start (Monday, ISO), answer, n_samples, n_well,
                      cities: {city: verdict of that city's samples that week}}]}"""
    try:
        from process.sample_qc import well_covered
    except Exception:                                   # pragma: no cover
        well_covered = lambda _q: False                 # noqa: E731
    markers, samples, own = [], [], None
    for city, d in (per_city or {}).items():
        if not d:
            continue
        if d.get("own") is not None and own is None:
            own = {"parent": d.get("parent"), "markers": list(d["own"])}
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
            cities[c] = _vote(markers, [s for s in ss if s["city"] == c], own,
                              present_freq)["verdict"]
        weeks.append({"start": wk.isoformat(), "answer": _vote(markers, ss, own, present_freq),
                      "n_samples": len(ss), "n_well": sum(s["well"] for s in ss),
                      "cities": cities})
    return {"answer": _vote(markers, samples, own, present_freq),
            "n_cities": len({s["city"] for s in samples}), "weeks": weeks}


def sudden(nat: Optional[dict], quiet_weeks: int = 2, min_cities: int = 2) -> Optional[dict]:
    """A lineage that appears suddenly (2026-10-23): seen in the newest week,
    in at least `min_cities` cities, after `quiet_weeks` weeks in which all
    cities pooled said not seen, each on at least 2 well-covered samples (a
    trustworthy "not seen"). PJ.2.1, 5 Sep – 5 Oct: not seen on 2–5
    well-covered samples a week, then seen on 3–4 Oct in Lugano (64 % of
    counted reads), Basel, Zürich and Genève, each city with a different set
    of its mutations. A real arrival and a lab batch look the same here; the
    next samples tell them apart. None when it doesn't apply."""
    if not nat or len(nat.get("weeks") or []) < quiet_weeks + 1:
        return None
    wk = nat["weeks"]
    last, before = wk[-1], wk[-1 - quiet_weeks:-1]
    if last["answer"]["verdict"] != "present":
        return None
    if not all(w["answer"]["verdict"] == "absent" and w["n_well"] >= 2 for w in before):
        return None
    cities = [c for c, v in last["cities"].items() if v == "present"]
    if len(cities) < min_cities:
        return None
    return {"week": last["start"], "cities": cities,
            "quiet": [(w["start"], w["n_well"]) for w in before]}


def grade(per_city: Dict[str, dict], cfg: Optional[dict] = None) -> Optional[dict]:
    """How sure a found lineage is, from all cities pooled (2026-10-23).

    per_city: {city: check_in_data result}. Returns None without samples, else
      cities       cities whose own check says seen
      pooled_seen  all cities pooled say seen at check.confirm_pooled_freq (5 %)
      weeks, n_weeks  weeks of the all-cities line that say seen, of all weeks
      missing, n_measured  its mutations clearly missing pooled (< 1 %), of
                   those measured
      confirmed    pooled_seen, and in >= check.confirm_min_cities cities or
                   >= check.confirm_min_weeks weeks, and missing <=
                   check.confirm_max_missing of n_measured
    bench_scorecard (6 cities, 7 Jun – 5 Sep and 8 Mar – 6 Jun): confirmed held
    no gone or not-in-Europe lineage; found-and-not-confirmed ("early") held
    the one not-in-Europe false alarm (RF.5.1, Zürich) and real first
    sightings (TK.1, PQ.16.1.1, RV.2.1.1). Finding stays per city: pooling
    only grades (finding by pooling added up mutations from different
    cities, XFG.6.3.3)."""
    from process.cooc import _check_cfg
    c = _check_cfg(cfg)
    conf = national(per_city, c["confirm_pooled_freq"])
    if not conf:
        return None
    a = conf["answer"]
    cities = sorted(k for k, d in (per_city or {}).items() if (d or {}).get("state") == "present")
    weeks = sum(w["answer"]["verdict"] == "present" for w in conf["weeks"])
    missing = sorted(m for m, x in (a.get("markers") or {}).items() if x.get("status") == "absent")
    n_meas = a.get("n_measured", 0)
    pooled_seen = a["verdict"] == "present"
    ok_missing = len(missing) <= int(c["confirm_max_missing"] * n_meas)
    confirmed = bool(pooled_seen and ok_missing and (len(cities) >= c["confirm_min_cities"]
                                                     or weeks >= c["confirm_min_weeks"]))
    return {"cities": cities, "pooled_seen": pooled_seen, "weeks": weeks, "n_weeks": len(conf["weeks"]),
            "missing": missing, "n_measured": n_meas, "ok_missing": ok_missing,
            "confirmed": confirmed}
