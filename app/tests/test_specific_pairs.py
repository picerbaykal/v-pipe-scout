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
