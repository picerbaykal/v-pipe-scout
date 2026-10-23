"""Found lineages graded on all cities pooled (process.national.grade,
2026-10-23): "New · confirmed" when all cities pooled say seen, in >= 2 cities
or >= 2 weeks, with at most 10 % of its mutations clearly missing; otherwise
"New · early". A seen cell in a city with fewer than 2 well-covered samples is
striped. bench_scorecard 23 Oct: confirmed held no gone or not-in-Europe
lineage in summer or spring."""
import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, "/app_shared")

from process.national import grade                   # noqa: E402
from process.variant_explorer import check_in_data    # noqa: E402

MK = ["m1", "m2", "m3", "m4"]
WELL = {"median": 5000, "share_100": 0.9}
THIN = {"median": 8, "share_100": 0.05}   # counts as 80 reads a sample


def _city(days, share=0.3, q=WELL, zero=()):
    """one city: {date: share}; markers in `zero` at 0 %"""
    h = int(share * 1000)
    pd = {d: {m: [1000, 0, 0, 0] if m in zero else [1000, h, h, h] for m in MK} for d in days}
    return check_in_data(MK, pd, list(pd), quality={d: q for d in days})


def test_two_cities_confirm():
    g = grade({"GE": _city(["2026-09-15"]), "ZH": _city(["2026-09-16"])})
    assert g["confirmed"] and g["cities"] == ["GE", "ZH"] and g["weeks"] == 1


def test_one_city_one_week_is_early():
    g = grade({"ZH": _city(["2026-09-15"]), "GE": _city(["2026-09-16"], share=0.0)})
    assert g["pooled_seen"] and g["cities"] == ["ZH"] and not g["confirmed"]


def test_one_city_two_weeks_confirms():
    g = grade({"ZH": _city(["2026-09-15", "2026-09-22"])})
    assert g["weeks"] == 2 and g["confirmed"]


def test_mutations_missing_pooled_keep_it_early():
    # PJ.2.1, Sep–Oct: 4 of 18 at 0.0–0.2 % pooled; here 1 of 4 at 0 %
    g = grade({"GE": _city(["2026-09-15"], zero=("m4",)), "ZH": _city(["2026-09-16"], zero=("m4",))})
    assert g["missing"] == ["m4"] and not g["ok_missing"] and not g["confirmed"]


def _vt():
    path = Path(__file__).resolve().parents[1] / "components" / "variants_table.py"
    spec = importlib.util.spec_from_file_location("vt_grade", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_row_says_confirmed_or_early_and_why():
    vt = _vt()
    C = ["Genève (GE)", "Zürich (ZH)"]
    two = {"per_city": {c: {"days": 3, "check": _city([d])} for c, d in
                        zip(C, ("2026-09-15", "2026-09-16"))}, "addable": True}
    kind, ov, _tip, dot, _n, dk = vt._finding_overview("TA.1", two, C, [])
    assert kind == "new" and "New · confirmed" in ov and "seen in GE, ZH, 1 of 1 week" in ov
    assert dk == "filled"
    one = {"per_city": {C[1]: {"days": 3, "check": _city(["2026-09-16"])},
                        C[0]: {"days": 3, "check": _city(["2026-09-15"], share=0.0)}},
           "addable": True}
    kind, ov, _tip, _dot, _n, dk = vt._finding_overview("RF.5.1", one, C, [])
    assert "New · early" in ov and "seen in ZH only, 1 of 1 week" in ov and dk == "ring"


def test_seen_cell_in_a_thin_city_is_striped():
    vt = _vt()
    thin = _city(["2026-09-15", "2026-09-16", "2026-09-17"], q=THIN)
    well = _city(["2026-09-15", "2026-09-16"])
    assert thin["state"] == "present" and thin["n_well"] == 0
    assert vt._scale(thin)[0] == vt.S_YES + vt._LOWQ_CELL
    assert vt._scale(well)[0] == vt.S_YES
    assert "striped: 0 of 3 samples well covered" in vt._cell_tip("TA.1", "Basel", thin)
