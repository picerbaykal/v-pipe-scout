"""Found lineages graded on all cities pooled (process.national.grade,
2026-10-23): "New · confirmed" when all cities pooled say seen, in >= 2 well-covered
cities and >= 2 weeks (2026-10-26; was "or", thin cities counted), with at most 10 % of its mutations clearly missing; otherwise
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


W1, W2 = ["2026-09-15", "2026-09-16"], ["2026-09-22", "2026-09-23"]


def test_two_cities_two_weeks_confirm():
    g = grade({"GE": _city(W1 + W2), "ZH": _city(W1 + W2)})
    assert g["confirmed"] and g["cities"] == ["GE", "ZH"] and g["weeks"] == 2


def test_two_cities_one_week_is_early():
    # RV.2.1.1, Sep-Oct: GE and ZH, 1 of 4 weeks (confirmed until 26 Oct)
    g = grade({"GE": _city(W1), "ZH": _city(W1)})
    assert g["cities"] == ["GE", "ZH"] and g["weeks"] == 1 and not g["confirmed"]


def test_one_city_two_weeks_is_early():
    g = grade({"ZH": _city(W1 + W2), "GE": _city(W1 + W2, share=0.0)})
    assert g["pooled_seen"] and g["cities"] == ["ZH"] and g["weeks"] == 2 and not g["confirmed"]


def test_thin_city_is_shown_not_counted():
    # PJ.2.1, Sep-Oct: "seen in BS, GR, TI", Basel 0 and Chur 1 well covered
    g = grade({"TI": _city(W1 + W2), "BS": _city(W1 + W2, q=THIN)})
    assert g["cities"] == ["TI"] and g["cities_thin"] == ["BS"] and not g["confirmed"]


def _own(d, own, parent="P"):
    d["own"], d["parent"] = list(own), parent
    return d


def test_mutations_missing_pooled_keep_it_early(monkeypatch):
    # with the absent check on (off by default since 2026-10-27)
    import process.cooc as pc
    monkeypatch.setitem(pc.CHECK_DEFAULTS, "confirm_use_missing", True)
    # PJ.2.1, Sep-Oct: 4 of 19 absent on its samples, all inherited from PJ.2
    g = grade({c: _own(_city(W1 + W2, zero=("m3", "m4")), ["m1", "m2"]) for c in ("GE", "ZH")})
    assert g["missing"] == ["m3", "m4"] and not g["ok_missing"] and not g["confirmed"]
    assert g["missing_inherited"] == ["m3", "m4"]


def test_absent_mutations_dont_block_by_default():
    # 2026-10-27: confirmed = cities and weeks; an absent mutation is reported only
    import process.cooc as pc
    g = grade({c: _city(W1 + W2, zero=("m4",)) for c in ("GE", "ZH")},
              cfg={"confirm_missing_allow": 0})
    assert g["missing"] == ["m4"] and not g["ok_missing"] and g["confirmed"]
    assert pc.CHECK_DEFAULTS["confirm_use_missing"] is False


def test_one_absent_mutation_is_always_allowed(monkeypatch):
    import process.cooc as pc
    monkeypatch.setitem(pc.CHECK_DEFAULTS, "confirm_use_missing", True)
    # 2026-10-27: TA.1's 14808C at 0.0 % in Sep-Oct (7 mutations: 10 % allowed none)
    g = grade({c: _city(W1 + W2, zero=("m4",)) for c in ("GE", "ZH")})
    assert g["missing"] == ["m4"] and g["ok_missing"] and g["confirmed"]


def _vt():
    path = Path(__file__).resolve().parents[1] / "components" / "variants_table.py"
    spec = importlib.util.spec_from_file_location("vt_grade", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_row_says_confirmed_or_early_and_why():
    vt = _vt()
    C = ["Genève (GE)", "Zürich (ZH)"]
    two = {"per_city": {c: {"days": 3, "check": _city(W1 + W2)} for c in C}, "addable": True}
    # confirmed this month, but + Add needs weeks in a row (2026-10-27): Watch
    kind, ov, tip, dot, _n, dk = vt._finding_overview("TA.1", two, C, [])
    assert kind == "new" and "Watch" in ov and "re-run for weeks" in ov
    assert "seen in GE, ZH, 2 of 2 weeks" in tip and dk == "ring"
    ready = dict(two, add_ready=True, weekly={"in_a_row": 3, "need": 3, "weeks": []})
    kind, ov, _tip, _dot, _n, dk = vt._finding_overview("TA.1", ready, C, [])
    assert "Add to panel" in ov and "3 weeks in a row" in ov and "data-add" in ov and dk == "filled"
    wk1 = {"per_city": {c: {"days": 3, "check": _city(W1)} for c in C}, "addable": True}
    kind, ov, tip, _dot, _n, dk = vt._finding_overview("RV.2.1.1", wk1, C, [])
    assert "Watch" in ov and "seen in GE, ZH, 1 of 1 week" in tip and dk == "ring"
    assert "and 2 weeks" in tip
    one = {"per_city": {C[1]: {"days": 3, "check": _city(W1 + W2)},
                        C[0]: {"days": 3, "check": _city(W1 + W2, share=0.0)}},
           "addable": True}
    kind, ov, _tip, _dot, _n, dk = vt._finding_overview("RF.5.1", one, C, [])
    assert "Watch" in ov and "· 1 city" in ov and dk == "ring"


def test_row_names_thin_cities_and_inherited_missing(monkeypatch):
    # 4 mutations, 1 absent: none allowed here, to read the wording
    import process.cooc as pc
    monkeypatch.setitem(pc.CHECK_DEFAULTS, "confirm_missing_allow", 0)
    monkeypatch.setitem(pc.CHECK_DEFAULTS, "confirm_use_missing", True)
    vt = _vt()
    C = ["Basel (BS)", "Lugano (TI)"]
    own = ["m1", "m2", "m3"]
    z = ("m4",)
    f = {"per_city": {C[0]: {"days": 3, "check": _own(_city(W1 + W2, q=THIN, zero=z), own, "PJ.2")},
                      C[1]: {"days": 3, "check": _own(_city(W1 + W2, zero=z), own, "PJ.2")}},
         "addable": True}
    kind, ov, tip, _dot, _n, dk = vt._finding_overview("PJ.2.1", f, C, [])
    assert "Watch" in ov and "· 1 of 4 mutations absent" in ov
    assert "seen in TI (+ BS thin)" in tip
    assert "1 of its 4 mutations missing, all inherited from PJ.2" in tip


def test_seen_cell_in_a_thin_city_is_striped():
    vt = _vt()
    thin = _city(["2026-09-15", "2026-09-16", "2026-09-17"], q=THIN)
    well = _city(["2026-09-15", "2026-09-16"])
    assert thin["state"] == "present" and thin["n_well"] == 0
    assert vt._scale(thin)[0] == vt.S_YES + vt._LOWQ_CELL
    assert vt._scale(well)[0] == vt.S_YES
    assert "striped: 0 of 3 samples well covered" in vt._cell_tip("TA.1", "Basel", thin)


def test_together_hover_names_the_cooccurring_mutations():
    # TA.1, Zürich: 23157G's reads reach 23148C and 23401T (W-ASAP amplicon 75,
    # Jaccard 1.0); 6196A's reach only XFG's 6183G (amplicon 20, .97)
    vt = _vt()
    from collections import Counter
    sigs = {"TA.1": {"23148C", "23157G", "23401T", "23005A", "6196A", "6183G", "21249A"},
            "P": {"23148C", "23005A", "6183G"}}
    sigs.update({f"L{i}": {"23005A", "6183G"} for i in range(600)})
    par = {"TA.1": "P"}
    vt._pango_sig = lambda: (sigs, par, Counter(m for v in sigs.values() for m in v))
    tip = vt._together_tip("TA.1", {"marker": "23157G", "link": 0.99, "link_n": 4012})
    assert tip.startswith("99 % of 4,012 reads carrying 23157G")
    assert "23148C (in 2 lineages)" in tip and "23401T (new in TA.1, in 1 lineage)" in tip
    assert "1 common one (in 602 lineages: 23005A)" in tip
    assert "judged on its share alone" in vt._together_tip("TA.1", {"marker": "21249A"})
