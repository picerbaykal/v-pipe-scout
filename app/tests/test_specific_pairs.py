"""process.cooc.specific_pairs (2026-10-26): pairs of a lineage's mutations
one read can carry, specific by the ★ marker rule applied to the pair."""
import types

import process.recombinants as rec
from process.cooc import specific_pairs


class _RC:
    def rec_root(self, v):
        return v

    def tolerated(self, l, v, m, own):
        return False


def _setup(monkeypatch):
    monkeypatch.setattr(rec, "get_recombinants", lambda *a, **k: _RC())


def test_common_mutations_can_make_a_specific_pair(monkeypatch):
    # TA.1-like: 23148C and 23401T each in many lineages, together only in TA.1
    _setup(monkeypatch)
    lin = {"TA.1": {"23148C", "23401T", "23005A", "6196A"},
           "P": {"23005A"}}
    lin.update({f"A{i}": {"23148C", "23005A"} for i in range(30)})
    lin.update({f"B{i}": {"23401T", "23005A"} for i in range(30)})
    parent = {"TA.1": "P", "P": "", **{k: "" for k in lin if k[0] in "AB"}}
    out = specific_pairs({"TA.1": lin["TA.1"]}, lin, parent, {"pair_span": 300, "out_other_max": 5})
    pairs = {d["pair"]: d for d in out["TA.1"]}
    assert "23148C+23401T" in pairs and pairs["23148C+23401T"]["n_out"] == 0
    assert pairs["23148C+23401T"]["markers"] == 0
    assert "23005A+23148C" not in pairs          # 30 lineages carry both
    assert not any("6196A" in p for p in pairs)  # too far from the others


def test_pair_carried_by_another_panel_variant_is_not_specific(monkeypatch):
    _setup(monkeypatch)
    lin = {"X": {"100A", "200C"}, "Y": {"100A", "200C", "300G"}}
    parent = {"X": "", "Y": ""}
    out = specific_pairs({"X": lin["X"], "Y": lin["Y"]}, lin, parent,
                         {"pair_span": 300, "out_other_max": 5})
    assert out["X"] == []
    assert "200C+300G" in {d["pair"] for d in out["Y"]}


def test_pairs_are_tallied_like_markers():
    from process.cooc import accumulate_check_stats
    rows = [{"date": "2026-10-03", "count": 70, "[23148]": "C", "[23157]": "G", "[23401]": "T"},
            {"date": "2026-10-03", "count": 20, "[23148]": "C", "[23157]": "G", "[23401]": "G"},
            {"date": "2026-10-03", "count": 10, "[23148]": "N", "[23157]": "T", "[23401]": "T"},
            {"date": "2026-10-03", "count": 900, "[23148]": "T", "[23157]": "T", "[23401]": "G"}]
    sig = {"TA.1": {"23148C", "23157G", "23401T"}}
    stats = {}
    accumulate_check_stats(rows, [23148, 23157, 23401], sig, {"TA.1": ["23157G"]}, stats,
                           pairs={"TA.1": ["23148C+23401T"]})
    day = stats["TA.1"]["2026-10-03"]
    assert day["23148C+23401T"] == [990, 70, 70, 70]   # 10 reads miss 23148: not counted
    assert day["23157G"][:2] == [1000, 90]             # the marker as before


def test_no_pairs_changes_nothing():
    from process.cooc import accumulate_check_stats
    rows = [{"date": "2026-10-03", "count": 5, "[100]": "A", "[200]": "C"}]
    a, b = {}, {}
    accumulate_check_stats(rows, [100, 200], {"X": {"100A", "200C"}}, {"X": ["100A"]}, a)
    accumulate_check_stats(rows, [100, 200], {"X": {"100A", "200C"}}, {"X": ["100A"]}, b,
                           pairs={})
    assert a == b and "100A+200C" not in a["X"]["2026-10-03"]


ZH_PJ = {  # PJ.2.1, Zürich, Sep-Oct (diag_pairs, 26 Oct): [cov, hit, link_n, link_ok]
    "2026-09-20": {"26985A+27259C": [37494, 14023, 14023, 14023],
                   "22879A+22899C": [42178, 0, 0, 0]},
    "2026-09-29": {"22879A+22899C": [19792, 0, 0, 0], "24077A+24079A": [19167, 0, 0, 0]},
    "2026-10-03": {"22879A+22899C": [10097, 7815, 7815, 7815],
                   "24077A+24079A": [14805, 8824, 8824, 8824],
                   "23429T+23585A": [3207, 1392, 1392, 1392],
                   "24999G+25067T": [23237, 0, 0, 0]},
}
WELL = {"median": 5000, "share_100": 0.9}


def test_two_stretches_on_one_sample_make_it_positive():
    from process.cooc import pair_evidence
    ev = pair_evidence(ZH_PJ, {d: WELL for d in ZH_PJ})
    assert ev["positive"] == ["2026-10-03"] and ev["positive_well"] == ["2026-10-03"]
    assert ev["n_stretches"] == 3                    # around 22879, 23429 and 24077
    assert [p[0] for p in ev["per_date"]["2026-10-03"]][0] == "22879A+22899C"
    # 20 Sep: one stretch alone (26985A+27259C at 37 %) is not enough
    assert "2026-09-20" not in ev["positive"]


def test_no_pair_counts_gives_none():
    from process.cooc import pair_evidence
    assert pair_evidence({"2026-10-03": {"100A": [500, 10, 0, 0]}}) is None


def test_pairs_turn_a_hidden_city_into_seen():
    # the marker vote says not seen (most single markers 0 %); the pairs say seen
    from process.variant_explorer import check_in_data
    mk = ["m1", "m2", "m3", "m4"]
    pd = {d: {m: [1000, 0, 0, 0] for m in mk} for d in ZH_PJ}
    for d, cells in ZH_PJ.items():
        pd[d].update(cells)
    r = check_in_data(mk, pd, list(pd), quality={d: WELL for d in pd})
    assert r["state_markers"] == "absent" and r["state"] == "present"
    assert r["pairs"]["positive"] == ["2026-10-03"]


def test_cell_and_grade_use_the_pairs():
    import importlib.util
    from pathlib import Path
    from process.national import grade
    from process.variant_explorer import check_in_data
    spec = importlib.util.spec_from_file_location(
        "vt_pairs", Path(__file__).resolve().parents[1] / "components" / "variants_table.py")
    vt = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(vt)
    mk = ["m1", "m2", "m3", "m4"]

    def city(q):
        pd = {d: {m: [1000, 0, 0, 0] for m in mk} for d in ZH_PJ}
        for d, cells in ZH_PJ.items():
            pd[d].update(cells)
        return check_in_data(mk, pd, list(pd), quality={d: q for d in pd})
    zh = city(WELL)
    style, label, *_ = vt._scale(zh, found=True)
    assert label == "3\non 1 sample"
    assert "22879A+22899C 77 % of 10,097" in vt._cell_tip("PJ.2.1", "Zürich", zh)
    thin = city({"median": 8, "share_100": 0.05})
    g = grade({"ZH": zh, "BS": thin})
    assert g["cities"] == ["ZH"] and g["cities_thin"] == ["BS"]


def test_every_cell_with_pair_counts_speaks_one_measure():
    import importlib.util
    from pathlib import Path
    from process.variant_explorer import check_in_data
    spec = importlib.util.spec_from_file_location(
        "vt_pairs2", Path(__file__).resolve().parents[1] / "components" / "variants_table.py")
    vt = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(vt)
    mk = ["m1", "m2", "m3"]
    dates = ["2026-09-20", "2026-10-03"]
    # markers seen, pairs read but on one stretch only -> light green, old count
    pd = {d: {**{m: [1000, 300, 300, 300] for m in mk}, "100A+200C": [500, 0, 0, 0]}
          for d in dates}
    r = check_in_data(mk, pd, dates, quality={d: WELL for d in dates})
    sty, lab, *_ = vt._scale(r, found=True)
    assert sty == vt.S_SINGLE and lab == "3:0"
    assert "single mutations only" in vt._cell_tip("ST.1", "Zürich", r)
    # markers not seen, pairs read, none -> 0 and the samples measured
    pd = {d: {**{m: [1000, 0, 0, 0] for m in mk}, "100A+200C": [500, 0, 0, 0]} for d in dates}
    r = check_in_data(mk, pd, dates, quality={d: WELL for d in dates})
    sty, lab, *_ = vt._scale(r, found=True)
    assert lab == "0\n2 samples" and "not seen" in vt._cell_tip("X", "Zürich", r)
