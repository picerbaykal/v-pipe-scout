"""The "covered by" rule (process.family.family_rows, 2026-10-23): X itself
is a share of READS compared with the 5 % limit; a city with fewer than 2
well-covered samples can't tell; the lineage below is read at the markers its
own check sees. XFG.14.8.1 next to TA.1 (5 Sep – 5 Oct) stayed "New"."""
import sys

sys.path.insert(0, "/app_shared")

from process.family import family_rows        # noqa: E402

XM = ["100A", "200C"]                 # X's markers: TA.1 carries both
DM = ["100A", "200C", "300G", "400T"]  # D's markers; 400T weak (not seen)


def _chk(n_well, detail):
    return {"state": "present", "n_well": n_well, "basis": "weighted",
            "quality": {"2026-09-10": {"median": 5000, "share_100": 0.9}},
            "detail": [{"marker": m, "status": s} for m, s in detail.items()]}


def _case(n_well_x, share_d, exclude=None):
    pd_x = {"2026-09-10": {m: [1000, 200, 200, 200] for m in XM}}          # X 20 %
    pd_d = {"2026-09-10": {**{m: [1000, int(share_d * 1000), 0, 0] for m in DM[:3]},
                           "400T": [1000, 20, 0, 0]}}                      # 400T at 2 %
    fc = {"below_all": {"X": ["D"], "D": []}, "sigs": {"X": XM, "D": DM},
          "per_variant": {"X": {"markers": XM, "per_date": pd_x},
                          "D": {"markers": DM, "per_date": pd_d}}}
    findings = {
        "X": {"status": "confirmed", "per_city": {"GE": {"check": _chk(n_well_x, {m: "present" for m in XM})}}},
        "D": {"status": "confirmed", "per_city": {"GE": {"check": _chk(7, {
            "100A": "present", "200C": "present", "300G": "present", "400T": "absent"})}}}}
    family_rows(findings, {"GE": {"findings_check": fc}}, {}, [], exclude=exclude)
    return findings["X"]


def test_x_whose_markers_the_lineage_below_carries_is_covered():
    x = _case(7, 0.19)
    assert x["status"] == "family" and x["descendants"] == ["D"]
    val, _parts = x["itself"]["GE"]
    assert val < 0.05          # a share of reads, not the left-over fraction


def test_x_with_a_real_share_of_its_own_stays():
    x = _case(7, 0.05)         # D explains 5 of X's 20 %: 15 % of reads left
    assert x["status"] != "family" and x["itself"]["GE"][0] > 0.05


def test_a_thin_city_cannot_tell():
    # its share can't be measured there; with the lineage below seen in the
    # same (only) city it is covered by it (2026-10-27, thin-only rule)
    x = _case(0, 0.05)
    assert x["itself"]["GE"] == (None, []) and x["status"] == "family"


def test_a_folded_lineage_below_covers_nothing():
    # 2026-10-27: XFG.6.3.3 "Covered by SW.2.1" while SW.2.1 was folded (coocs
    # on 1 sample, 4 of 5 mutations missing)
    assert _case(7, 0.19, exclude={"D"})["status"] != "family"
    assert _case(7, 0.19, exclude={"X"})["status"] != "family"


def test_thin_only_parent_is_covered_by_its_child_seen_there():
    # 2026-10-27: PJ.2 seen only in thin Basel and Chur, PJ.2.1 seen in both
    findings = {"X": {"status": "confirmed", "per_city": {}}, "D": {"status": "confirmed", "per_city": {}}}
    sc = {}
    for loc in ("BS", "GR"):
        fc = {"below_all": {"X": ["D"], "D": []}, "sigs": {"X": XM, "D": DM},
              "per_variant": {"X": {"markers": XM, "per_date": {}}, "D": {"markers": DM, "per_date": {}}}}
        sc[loc] = {"findings_check": fc}
        findings["X"]["per_city"][loc] = {"check": _chk(0, {m: "present" for m in XM})}
        findings["D"]["per_city"][loc] = {"check": _chk(0, {m: "present" for m in DM})}
    family_rows(findings, sc, {}, [])
    assert findings["X"]["status"] == "family" and findings["X"]["descendants"] == ["D"]
    # a well-covered city without the child: not covered
    findings["X"]["status"] = "confirmed"
    findings["X"]["per_city"]["ZH"] = {"check": _chk(7, {m: "present" for m in XM})}
    sc["ZH"] = {"findings_check": {"below_all": {"X": ["D"], "D": []}, "sigs": {"X": XM, "D": DM},
                                   "per_variant": {"X": {"markers": XM, "per_date": {}}}}}
    family_rows(findings, sc, {}, [])
    assert findings["X"]["status"] != "family"


def test_thin_city_checked_by_the_cross_check_uses_the_tree_from_other_cities():
    # PJ.2 in Basel and Chur came from the cross-check: no "below" list there;
    # Zürich's scan has PJ.2 > PJ.2.1
    findings = {"X": {"status": "confirmed", "per_city": {}}, "D": {"status": "confirmed", "per_city": {}}}
    sc = {"ZH": {"findings_check": {"below_all": {"X": ["D"], "D": []}, "per_variant": {}}}}
    for loc in ("BS", "GR"):
        sc[loc] = {"findings_check": {"below_all": {}, "per_variant": {}}}
        findings["X"]["per_city"][loc] = {"check": _chk(0, {m: "present" for m in XM})}
        findings["D"]["per_city"][loc] = {"check": _chk(0, {m: "present" for m in DM})}
    family_rows(findings, sc, {}, [])
    assert findings["X"]["status"] == "family" and findings["X"]["descendants"] == ["D"]
