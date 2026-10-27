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


def national(per_city: Dict[str, dict], present_freq: Optional[float] = None,
             only: Optional[Dict[str, set]] = None) -> Optional[dict]:
    """The pooled answer for the chosen dates and one per week, or None
    without samples. present_freq: the "seen" share (default check.present_freq).
    only: {city: {ISO date}} — pool only these samples (other cities left out;
    grade(), 2026-10-27: the samples where its coocs were seen).

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
        if only is not None and city not in only:
            continue
        for dt, cells in (d.get("samples") or {}).items():
            day = str(dt)[:10]
            if only is not None and day not in only[city]:
                continue
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
    """Not shown since 2026-10-27 (the marker vote; the weekly grade,
    weekly(), shows an arrival). A lineage that appears suddenly (2026-10-23): seen in the newest week,
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


def city_seen(d: Optional[dict]) -> bool:
    """Does this city count as seen (2026-10-26)? With pair counts, only when a
    sample had its own coocs (process.cooc.pair_evidence: the solid green
    cell); single mutations alone (light green) don't count. Without pair
    counts (older results), the marker vote as before. Pair counts looked for
    but none read (no pair of its own mutations, or no reads at them;
    2026-10-27): not seen — single mutations never count."""
    d = d or {}
    pev = d.get("pairs") or {}
    if pev.get("measured") or pev.get("positive"):
        return bool(pev.get("positive"))
    if d.get("pairs_checked"):
        return False
    return d.get("state") == "present"


def _pair_weeks(per_city: Dict[str, dict]) -> Optional[set]:
    """Mondays of the well-covered samples that counted, over all cities;
    None when no city has pair counts."""
    if not any(((d or {}).get("pairs") or {}).get("measured") for d in (per_city or {}).values()):
        return None
    return {_monday(x) for d in (per_city or {}).values()
            for x in (((d or {}).get("pairs") or {}).get("positive_well") or [])}


def counted_samples(per_city: Dict[str, dict]) -> Optional[Dict[str, set]]:
    """{city: {ISO dates}} of the samples where its own coocs were seen
    (process.cooc.pair_evidence "positive"); None when no city has pair
    counts (results before 2026-10-26) or none was seen."""
    pcs = {c: (d or {}).get("pairs") or {} for c, d in (per_city or {}).items()}
    if not any(p.get("measured") for p in pcs.values()):
        return None
    out = {c: {str(x)[:10] for x in p.get("positive") or []} for c, p in pcs.items()}
    out = {c: ds for c, ds in out.items() if ds}
    return out or None


def grade(per_city: Dict[str, dict], cfg: Optional[dict] = None) -> Optional[dict]:
    """How sure a found lineage is, from all cities pooled (2026-10-23).

    per_city: {city: check_in_data result}. Returns None without samples, else
      cities       cities seen (city_seen: their own coocs, or the marker vote
                   for results without pair counts), with at least
                   check.thin_well_min well-covered samples (2026-10-26)
      cities_thin  cities whose check says seen on thinner samples (shown,
                   not counted)
      pooled_seen  all cities pooled say seen at check.confirm_pooled_freq (5 %);
                   with pair counts, pooled over the samples where its coocs
                   were seen (n_counted, counted: {city: dates}; 2026-10-27)
      weeks, n_weeks  weeks with a well-covered sample that counted (pair
                   counts), else weeks of the all-cities line that say seen; of
                   all weeks
      missing, n_measured  its mutations clearly missing pooled (< 1 %; the same
                   samples as pooled_seen), of
                   those measured
      missing_inherited  of those, the ones it inherits (not its own markers),
                   when the result names its own markers
      confirmed    (2026-10-27: cities and weeks only; pooled_seen and
                   ok_missing count only with check.confirm_use_pooled /
                   confirm_use_missing on, both off by default)
                   pooled_seen, and in >= check.confirm_min_cities cities AND
                   >= check.confirm_min_weeks weeks, and missing <=
                   check.confirm_max_missing of n_measured, and never fewer than
                   check.confirm_missing_allow (1) allowed (2026-10-27)
    Until 2026-10-26 cities OR weeks, and thin cities counted: RV.2.1.1
    (Sep-Oct, 1 of 4 weeks) was confirmed while PJ.2.1 (1 of 4 weeks) was
    early; PJ.2.1's "seen in BS, GR, TI" rested on Basel (0 well covered)
    and Chur (1).
    bench_scorecard (6 cities, 7 Jun - 5 Sep and 8 Mar - 6 Jun, 23 Oct, with
    "or"): confirmed held no gone or not-in-Europe lineage; found-and-not-
    confirmed ("early") held the one not-in-Europe false alarm (RF.5.1,
    Zurich) and real first sightings (TK.1, PQ.16.1.1, RV.2.1.1). Finding
    stays per city: pooling only grades (finding by pooling added up
    mutations from different cities, XFG.6.3.3)."""
    from process.cooc import _check_cfg
    c = _check_cfg(cfg)
    conf = national(per_city, c["confirm_pooled_freq"])
    if not conf:
        return None
    # its mutations judged on the samples where its coocs were seen
    # (2026-10-27): pooled over every sample of every city, a rare lineage's
    # shares are diluted by the samples without it (PJ.2.1, Sep-Oct: coocs on 6
    # of 48 samples, "4 of 18 absent" on all 48)
    only = counted_samples(per_city)
    own_conf = national(per_city, c["confirm_pooled_freq"], only=only) if only else None
    a = (own_conf or conf)["answer"]
    n_counted = sum(len(ds) for ds in only.values()) if only and own_conf else None
    well_min = int(c.get("thin_well_min", 2))
    seen = {k: d for k, d in (per_city or {}).items() if city_seen(d)}
    def thin(d):
        # seen through specific co-occurrences only (2026-10-26): thin when
        # every sample that showed them was thin (Zürich PJ.2.1: 3 Oct well
        # covered counts; a city with only a thin positive sample doesn't)
        pev = d.get("pairs") or {}
        if pev.get("positive"):                  # the same samples the cell shows
            return not pev.get("positive_well")
        return d.get("n_well") is not None and d["n_well"] < well_min
    cities = sorted(k for k, d in seen.items() if not thin(d))
    cities_thin = sorted(k for k, d in seen.items() if thin(d))
    # weeks: with pair counts, the weeks with a well-covered sample that
    # counted (2026-10-26); before, the all-cities marker vote per week
    pw = _pair_weeks(per_city)
    weeks = (len(pw) if pw is not None
             else sum(w["answer"]["verdict"] == "present" for w in conf["weeks"]))
    missing = sorted(m for m, x in (a.get("markers") or {}).items() if x.get("status") == "absent")
    own = next((set(d["own"]) for d in (per_city or {}).values() if d and d.get("own") is not None), None)
    inherited = [m for m in missing if own is not None and m not in own]
    n_meas = a.get("n_measured", 0)
    pooled_seen = a["verdict"] == "present"
    ok_missing = len(missing) <= max(int(c.get("confirm_missing_allow", 1)),
                                     int(c["confirm_max_missing"] * n_meas))
    confirmed = bool((pooled_seen or not c.get("confirm_use_pooled", False))
                     and (ok_missing or not c.get("confirm_use_missing", False))
                     and len(cities) >= c["confirm_min_cities"]
                     and weeks >= c["confirm_min_weeks"])
    return {"cities": cities, "cities_thin": cities_thin, "pooled_seen": pooled_seen,
            "weeks": weeks, "n_weeks": len(conf["weeks"]), "missing": missing,
            "missing_inherited": inherited, "n_measured": n_meas, "ok_missing": ok_missing,
            "confirmed": confirmed, "n_counted": n_counted,
            "counted": {k: sorted(v) for k, v in (only or {}).items()} if n_counted else None}


def week_ends(d0: str, d1: str, days: int = 28) -> List[str]:
    """Sundays up to d1 whose `days`-day window starts on or after d0 (ISO).
    Samples after the last Sunday count from the next week on."""
    a, b = date.fromisoformat(str(d0)[:10]), date.fromisoformat(str(d1)[:10])
    x = b - timedelta(days=(b.weekday() + 1) % 7)          # last Sunday <= d1
    out = []
    while x - timedelta(days=days - 1) >= a:
        out.append(x.isoformat())
        x -= timedelta(days=7)
    return out[::-1]


def weekly(raw: Dict[str, dict], d0: str, d1: str, cfg: Optional[dict] = None) -> dict:
    """The grade on the 4 weeks up to each Sunday (2026-10-27), for + Add.

    raw: {city: {markers, per_date, dates, quality, own}} — the inputs of
    process.variant_explorer.check_in_data, over the chosen dates d0..d1.
    Each Sunday S with a full window (S - add_days + 1 >= d0): every city's
    check on the samples of that window, then grade().

    Returns {weeks: [{end, status, grade}], in_a_row, need, ready, short}
      status     confirmed / watch (seen in a city, not confirmed) / none
      in_a_row   confirmed weeks counted back from the newest
      need       check.add_weeks
      short      fewer weeks than `need` fit in the dates (a window under
                 add_days + 7 x (need - 1) days): + Add can't be judged
      ready      in_a_row >= need
    bench_scorecard 27 Oct: see check.add_weeks."""
    from process.cooc import _check_cfg
    from process.variant_explorer import check_in_data
    c = _check_cfg(cfg)
    need, days = int(c.get("add_weeks", 3)), int(c.get("add_days", 28))
    weeks = []
    for end in week_ends(d0, d1, days):
        start = (date.fromisoformat(end) - timedelta(days=days - 1)).isoformat()
        inside = lambda d: start <= str(d)[:10] <= end          # noqa: E731
        per = {}
        for city, r in (raw or {}).items():
            if not r or not r.get("markers"):
                continue
            pd = {d: x for d, x in (r.get("per_date") or {}).items() if inside(d)}
            ds = [d for d in (r.get("dates") or []) if inside(d)]
            if not ds and not pd:
                continue
            per[city] = check_in_data(r["markers"], pd, ds, quality=r.get("quality"),
                                      own=r.get("own"))
        g = grade(per, cfg) if per else None
        status = ("confirmed" if g and g.get("confirmed") else
                  "watch" if any(city_seen(x) for x in per.values()) else "none")
        weeks.append({"end": end, "status": status, "grade": g})
    n = 0
    for w in reversed(weeks):
        if w["status"] != "confirmed":
            break
        n += 1
    return {"weeks": weeks, "in_a_row": n, "need": need,
            "short": len(weeks) < need, "ready": n >= need}


def row_status(wk: Optional[dict]) -> Optional[dict]:
    """The row's status from its weekly grades (2026-10-27): the decision
    looks at the last weeks, whatever range was chosen; a longer range only
    adds history (cells show the whole range).

      add       the last check.add_weeks weeks confirmed in a row
      watch     the newest week confirmed (not yet in a row), or found
                there and never confirmed in the chosen dates
      fading    confirmed in an earlier week of the chosen dates, the newest
                week not confirmed
      earlier   never confirmed, not found in the newest week

    Returns {kind, in_a_row, need, last_confirmed, last_found} or None
    without weekly grades (dates shorter than one window)."""
    weeks = (wk or {}).get("weeks") or []
    if not weeks:
        return None
    need = int((wk or {}).get("need", 3))
    n = int((wk or {}).get("in_a_row", 0))
    conf = [w["end"] for w in weeks if w["status"] == "confirmed"]
    found = [w["end"] for w in weeks if w["status"] != "none"]
    newest = weeks[-1]["status"]
    if n >= need:
        kind = "add"
    elif newest == "confirmed":
        kind = "watch"
    elif conf:
        kind = "fading"
    elif newest == "watch":
        kind = "watch"
    else:
        kind = "earlier"
    return {"kind": kind, "in_a_row": n, "need": need,
            "last_confirmed": conf[-1] if conf else None,
            "last_found": found[-1] if found else None}


def check_start(d0: str, d1: str, cfg: Optional[dict] = None) -> str:
    """Start of the dates the found lineages' check reads (2026-10-27): the
    chosen start, or check.weekly_days before the end if that is earlier, so
    + Add (3 weekly windows of 28 days) can be judged on "last month". Only
    the cheap check reads further back; the deep scan keeps the chosen dates."""
    from process.cooc import _check_cfg
    c = _check_cfg(cfg)
    a = date.fromisoformat(str(d0)[:10])
    b = date.fromisoformat(str(d1)[:10]) - timedelta(days=int(c.get("weekly_days", 49)) - 1)
    return min(a, b).isoformat()


def trim_check(fc: Optional[dict], d0: str, d1: str) -> Optional[dict]:
    """A lineages_check result cut to the chosen dates (2026-10-27): the cells,
    the grade and the family rule see the chosen dates only; the weekly grade
    gets the untrimmed result."""
    if not fc:
        return fc
    a, b = str(d0)[:10], str(d1)[:10]
    inside = lambda d: a <= str(d)[:10] <= b                     # noqa: E731
    out = dict(fc)
    out["dates"] = [d for d in fc.get("dates") or [] if inside(d)]
    out["per_variant"] = {
        v: dict(pv, per_date={d: x for d, x in (pv.get("per_date") or {}).items() if inside(d)})
        for v, pv in (fc.get("per_variant") or {}).items()}
    return out
