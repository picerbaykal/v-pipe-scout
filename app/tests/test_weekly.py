"""+ Add needs the grade confirmed in the 4 weeks up to each of the last 3
Sundays, in a row (process.national.weekly, 2026-10-27). bench_scorecard 27
Oct, 6 cities, Mar-Oct: offered TA.1 and BA.3.2.2 only; one month confirmed
(TK.1, RV.2.1.1, PY.1.1.1) is not enough."""
import importlib.util
import sys
from datetime import date, timedelta
from pathlib import Path

sys.path.insert(0, "/app_shared")

from process.national import week_ends, weekly      # noqa: E402

MK = ["1000A", "1100C", "5000G", "5100T"]
PAIRS = ["1000A+1100C", "5000G+5100T"]
WELL = {"median": 5000, "share_100": 0.9}


def _days(d0, d1, every=3):
    a, b = date.fromisoformat(d0), date.fromisoformat(d1)
    out = []
    while a <= b:
        out.append(a.isoformat())
        a += timedelta(days=every)
    return out


def _city(dates, seen_from="0000"):
    """samples on `dates`; its coocs (and 30 % markers) from `seen_from` on"""
    pd = {}
    for d in dates:
        on = d >= seen_from
        h = 300 if on else 0
        c = {m: [1000, h, h, h] for m in MK}
        for p in PAIRS:
            c[p] = [1000, h, h, h]
        pd[d] = c
    return {"markers": MK, "per_date": pd, "dates": dates,
            "quality": {d: WELL for d in dates}, "own": None}


D0, D1 = "2026-08-10", "2026-10-05"


def test_week_ends_need_a_full_window():
    assert week_ends(D0, D1) == ["2026-09-06", "2026-09-13", "2026-09-20",
                                 "2026-09-27", "2026-10-04"]
    assert week_ends("2026-09-05", "2026-10-05") == ["2026-10-04"]


def test_seen_all_along_in_two_cities_is_ready():
    ds = _days(D0, D1)
    w = weekly({"GE": _city(ds), "ZH": _city(ds)}, D0, D1)
    assert [x["status"] for x in w["weeks"]] == ["confirmed"] * 5
    assert w["in_a_row"] == 5 and w["ready"] and not w["short"]


def test_one_month_is_not_enough():
    # coocs only from 21 Sep: the last window confirms, the ones before don't
    ds = _days(D0, D1)
    w = weekly({"GE": _city(ds, "2026-09-21"), "ZH": _city(ds, "2026-09-21")}, D0, D1)
    assert w["weeks"][-1]["status"] == "confirmed"
    assert w["in_a_row"] < 3 and not w["ready"]


def test_one_city_is_watch_not_ready():
    ds = _days(D0, D1)
    w = weekly({"ZH": _city(ds), "GE": _city(ds, "2099")}, D0, D1)
    assert {x["status"] for x in w["weeks"]} == {"watch"} and not w["ready"]


def test_a_gap_starts_the_count_again():
    ds = [d for d in _days(D0, D1, 1) if not ("2026-08-24" <= d <= "2026-09-20")]
    w = weekly({"GE": _city(ds), "ZH": _city(ds)}, D0, D1)
    # no samples 24 Aug - 20 Sep: confirmed, watch, none, watch, confirmed
    assert [x["status"] for x in w["weeks"]] == ["confirmed", "watch", "none", "watch", "confirmed"]
    assert w["in_a_row"] == 1 and not w["ready"]


def test_last_month_is_too_short():
    ds = _days("2026-09-05", D1)
    w = weekly({"GE": _city(ds), "ZH": _city(ds)}, "2026-09-05", D1)
    assert w["short"] and not w["ready"]


def _vt():
    path = Path(__file__).resolve().parents[1] / "components" / "variants_table.py"
    spec = importlib.util.spec_from_file_location("vt_weekly", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_button_only_when_ready():
    vt = _vt()
    from process.variant_explorer import check_in_data
    C = ["Genève (GE)", "Zürich (ZH)"]
    ds = _days(D0, D1)
    pc = {c: {"days": 3, "check": check_in_data(MK, _city(ds)["per_date"], ds,
                                                quality={d: WELL for d in ds})} for c in C}
    ready = weekly({c: _city(ds) for c in C}, D0, D1)
    f = {"per_city": pc, "addable": True, "add_ready": True, "weekly": ready}
    assert "data-add" in vt._finding_overview("TA.1", f, C, [])[1]
    late = weekly({c: _city(ds, "2026-09-21") for c in C}, D0, D1)
    f = dict(f, add_ready=False, weekly=late)
    ov = vt._finding_overview("TK.1", f, C, [])[1]
    assert "data-add" not in ov and "Watch" in ov and f"· {late['in_a_row']} of 3 weeks" in ov
    # dates too short: said once above the table, not on every row
    short = weekly({c: _city(_days("2026-09-05", D1)) for c in C}, "2026-09-05", D1)
    note = vt._add_note({"TK.1": dict(f, weekly=short), "TA.1": dict(f, weekly=short)})
    assert note.count("+ Add needs 3 weeks in a row") == 1 and "choose 42 days or more" in note
    assert vt._add_note({"TA.1": dict(f, weekly=ready)}) == ""


def test_last_month_check_reads_seven_weeks():
    # 2026-10-27: the found lineages' check reads check.weekly_days (49) back,
    # so 3 weekly windows fit on "last month"; the cells keep the chosen dates
    from process.national import check_start, trim_check
    d0 = check_start("2026-09-05", "2026-10-05")
    assert d0 == "2026-08-18"
    assert week_ends(d0, "2026-10-05")[-3:] == ["2026-09-20", "2026-09-27", "2026-10-04"]
    assert check_start("2026-06-07", "2026-10-05") == "2026-06-07"
    fc = {"dates": ["2026-08-20", "2026-09-10"],
          "per_variant": {"V": {"markers": ["1A"], "per_date": {"2026-08-20": {}, "2026-09-10": {}}}}}
    t = trim_check(fc, "2026-09-05", "2026-10-05")
    assert t["dates"] == ["2026-09-10"] and list(t["per_variant"]["V"]["per_date"]) == ["2026-09-10"]
    assert fc["dates"] == ["2026-08-20", "2026-09-10"]          # the original is kept
