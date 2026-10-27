"""Investigate a variant gives the tree's answer (2026-10-27): the same worker
result through both paths — the tree's found row (check_in_data on the
trimmed check, grade, weekly on the 7-week check) and process.variant_explorer
.judge — gives the same city answers, the same all-cities grade, the same weeks."""
import sys
from datetime import date, timedelta

sys.path.insert(0, "/app_shared")
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))

from process.national import check_start, grade, trim_check, weekly   # noqa: E402
from process.variant_explorer import check_in_data, judge              # noqa: E402

MK = ["1000A", "1100C", "5000G", "5100T"]
PAIRS = ["1000A+1100C", "5000G+5100T"]
WELL = {"median": 5000, "share_100": 0.9}
THIN = {"median": 40, "share_100": 0.2}


def _days(d0, d1, every=3):
    a, b = date.fromisoformat(d0), date.fromisoformat(d1)
    out = []
    while a <= b:
        out.append(a.isoformat())
        a += timedelta(days=every)
    return out


def _check(dates, seen_from, thin=()):
    pd = {}
    for d in dates:
        h = 300 if d >= seen_from else 0
        c = {m: [1000, h, h, h] for m in MK}
        for p in PAIRS:
            c[p] = [1000, h, h, h]
        pd[d] = c
    return {"dates": dates, "per_variant": {"X.1": {"markers": MK, "per_date": pd}}}, \
        {d: (THIN if d in thin else WELL) for d in dates}


D0, D1 = "2026-09-05", "2026-10-05"


def _both(cities):
    dc = check_start(D0, D1)
    ext = {c: _check(_days(dc, D1), s, t) for c, (s, t) in cities.items()}
    # the tree (abundance_cooc _fcheck / _grade / _raw + weekly)
    tree_pc = {}
    for c, (fc, q) in ext.items():
        tr = trim_check(fc, D0, D1)
        pv = tr["per_variant"]["X.1"]
        tree_pc[c] = check_in_data(pv["markers"], pv["per_date"], tr["dates"], quality=q)
    tree_raw = {c: {"markers": fc["per_variant"]["X.1"]["markers"],
                    "per_date": fc["per_variant"]["X.1"]["per_date"],
                    "dates": fc["dates"], "quality": q, "own": None}
                for c, (fc, q) in ext.items()}
    tree = {"per_city": tree_pc, "grade": grade(tree_pc), "weekly": weekly(tree_raw, dc, D1)}
    return tree, judge(tree_raw, D0, D1, dc)


def _same(a, b):
    assert set(a["per_city"]) == set(b["per_city"])
    for c in a["per_city"]:
        x, y = a["per_city"][c], b["per_city"][c]
        for k in ("state", "state_markers", "timeline", "pairs", "n_present", "n_measured"):
            assert x.get(k) == y.get(k), (c, k)
    assert a["grade"] == b["grade"]
    assert a["weekly"] == b["weekly"]


def test_seen_all_along():
    _same(*_both({"GE": ("0000", ()), "ZH": ("0000", ())}))


def test_seen_late_one_city_thin():
    _same(*_both({"GE": ("2026-09-21", ()), "ZH": ("2026-09-28", ("2026-09-28",)),
                  "BS": ("9999", ())}))


def test_judge_cuts_cells_to_the_chosen_dates():
    _t, j = _both({"GE": ("0000", ())})
    assert all(D0 <= d <= D1 for d, _m in j["per_city"]["GE"]["timeline"])
    assert j["weekly"]["weeks"]                     # the weeks see 7 weeks back


def test_view_draws_a_tree_row_with_open_details():
    from components.variants_table import investigate_view
    _t, j = _both({"Genève (GE)": ("0000", ()), "Zürich (ZH)": ("0000", ())})
    cities = ["Genève (GE)", "Zürich (ZH)"]
    for in_panel in (False, True):
        h = investigate_view("X.1", j, cities, cities[0], in_panel=in_panel, dates=(D0, D1))
        assert "class='vt'" in h and "X.1" in h
        assert "class='dt g-det-inv'" in h          # details open
