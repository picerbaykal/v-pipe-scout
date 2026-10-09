"""Unit tests for the co-occurrence check (process.cooc): specific_markers
(which markers) and accumulate_check_stats + check_verdicts (the marker vote).
Pure functions, no LAPIS."""
import sys

sys.path.insert(0, "/app_shared")

from process.cooc import (CHECK_DEFAULTS, accumulate_check_stats,  # noqa: E402
                          check_verdicts, specific_markers)
from process.recombinants import override_parents  # noqa: E402

CFG = dict(CHECK_DEFAULTS)

# Tiny tree: KP -> LF, KP -> NB, KP -> OTHER; XV and XW are recombinant roots
PARENT = {"KP": "", "LF": "KP", "XV": "", "XV.1": "XV", "XV.2": "XV",
          "NB": "KP", "XW": "", "OTHER": "KP"}
LINEAGES = {
    "KP":    {"100A", "200C"},
    "LF":    {"100A", "200C", "300G"},
    "XV":    {"100A", "200C", "300G", "400T", "500A"},
    "XV.1":  {"100A", "200C", "300G", "400T", "500A", "600C"},
    "XV.2":  {"100A", "200C", "300G", "400T", "500A"},
    "NB":    {"100A", "200C", "700G", "800T"},
    "XW":    {"100A", "400T"},          # recombinant that inherited 400T
    "OTHER": {"100A", "800T"},          # competitor carrying 800T
}
PANEL = {v: LINEAGES[v] for v in ("XV", "NB", "KP")}


def test_own_sublineages_do_not_count_against_a_variant():
    m = specific_markers(PANEL, LINEAGES, PARENT, CFG)
    assert "400T" in m["XV"] and "500A" in m["XV"]
    assert "700G" in m["NB"]
    assert m["KP"] == []                    # ancestor: nothing of its own


def test_strict_threshold_drops_markers_carried_elsewhere():
    # XW made from unrelated parents: its 400T counts against XV (since
    # 2026-10-02 only recombinants made FROM the variant are tolerated)
    with override_parents({"XW": ["KP", "NB"]}):
        m = specific_markers(PANEL, LINEAGES, PARENT, dict(CFG, out_other_max=0))
    assert "300G" not in m["XV"]            # parent LF carries it
    assert "800T" not in m["NB"]            # OTHER carries it
    assert "400T" not in m["XV"]            # XW carries it, not made from XV


def test_recombinant_made_from_the_variant_is_tolerated():
    # XW got its 400T from XV.1, so seeing 400T still means "XV material"
    with override_parents({"XW": ["XV.1", "NB"]}):
        m = specific_markers(PANEL, LINEAGES, PARENT, dict(CFG, out_other_max=0))
    assert "400T" in m["XV"]


SIG = {"XV": {"300G", "400T", "500A", "600C", "700G"}}


def _tally(rows_by_date, markers):
    stats = {}
    for rows in rows_by_date:
        accumulate_check_stats(rows, [300, 400, 500, 600, 700], SIG,
                               {"XV": markers}, stats)
    return stats.get("XV", {})


def _row(d, n, **bases):
    r = {"date": d, "count": n}
    r.update({f"[{k[1:]}]": v for k, v in bases.items()})
    return r


def test_present_when_markers_present_on_variant_reads():
    rows = [[_row("a", 900, p400="T", p300="G"), _row("a", 100, p400="C", p300="G"),
             _row("a", 500, p500="A", p300="G")]]
    out = check_verdicts(["400T", "500A"], _tally(rows, ["400T", "500A"]), CFG)
    assert out["verdict"] == "present"
    assert (out["n_present"], out["n_measured"]) == (2, 2)
    assert out["markers"]["400T"]["status"] == "present"


def test_absent_and_one_convergent_marker_does_not_flip_it():
    rows = [[_row("a", 900, p400="T", p300="G"), _row("a", 1100, p400="C", p300="G"),
             _row("a", 2000, p500="G", p600="T", p700="A", p300="G")]]
    mk = ["400T", "500A", "600C", "700G"]
    out = check_verdicts(mk, _tally(rows, mk), CFG)
    assert (out["n_present"], out["n_measured"]) == (1, 4)
    assert out["verdict"] == "absent"


def test_real_split_is_mixed():
    rows = [[_row("a", 2000, p400="T", p300="G"), _row("a", 2000, p500="G", p300="G")]]
    out = check_verdicts(["400T", "500A"], _tally(rows, ["400T", "500A"]), CFG)
    assert out["verdict"] == "mixed"


def test_marker_on_non_variant_reads_is_not_present():
    # 400T seen, but neighbour 300 is reference -> not on a variant haplotype
    rows = [[_row("a", 500, p400="T", p300="A"), _row("a", 500, p400="C", p300="A")]]
    out = check_verdicts(["400T"], _tally(rows, ["400T"]), CFG)
    assert out["markers"]["400T"]["status"] == "unmeasured"
    assert out["verdict"] == "not_covered"


def test_counts_are_pooled_over_dates():
    rows = [[_row("a", 60, p400="T", p300="G")], [_row("b", 60, p400="T", p300="G")]]
    out = check_verdicts(["400T"], _tally(rows, ["400T"]), dict(CFG, min_present=1))
    assert out["markers"]["400T"]["cov"] == 120
    assert out["verdict"] == "present"


def test_too_few_reads_or_no_markers_not_covered():
    assert check_verdicts([], {}, CFG)["verdict"] == "not_covered"
    rows = [[_row("a", 30, p400="T", p300="G")]]
    assert check_verdicts(["400T"], _tally(rows, ["400T"]), CFG)["verdict"] == "not_covered"


def test_positions_outside_batch_ignored():
    stats = {}
    accumulate_check_stats([_row("a", 900, p400="T", p500="A")], [500], SIG,
                           {"XV": ["400T"]}, stats)
    assert stats == {}

# ── at least min_present markers for "present" (2026-10-09) ───────────────
# One marker alone can come from an unrelated lineage that carries it (XFZ
# carries LP.8's 1954A), so one present marker is "mixed", never "present".

def _cells(**markers):
    """{marker: [cov, hit, link_n, link_ok]} with a link that passes."""
    return {m: [cov, hit, hit, hit] for m, (cov, hit) in markers.items()}


def test_one_present_marker_is_mixed_not_present():
    out = check_verdicts(["400T"], {"d1": _cells(**{"400T": (200, 100)})}, CFG)
    assert (out["n_present"], out["n_measured"]) == (1, 1)
    assert out["verdict"] == "mixed"


def test_two_present_markers_are_present():
    pd_ = {"d1": _cells(**{"400T": (200, 100), "500A": (200, 90)})}
    assert check_verdicts(["400T", "500A"], pd_, CFG)["verdict"] == "present"


def test_one_marker_of_an_unrelated_lineage_does_not_make_it_present():
    # LP.8 case: 1954A high (another lineage carries it), its other markers at 0
    pd_ = {"d1": _cells(**{"1954A": (500, 200), "100T": (500, 0), "200C": (500, 0)})}
    out = check_verdicts(["1954A", "100T", "200C"], pd_, CFG)
    assert out["verdict"] == "mixed"          # 1 of 3: not present, not absent


def test_day_marks_vote_per_marker_instead_of_pooling():
    from process.variant_explorer import check_in_data
    mk = ["1954A", "100T", "200C"]
    # pooled, this day was "present": 300 of 1,500 reads = 20 % >= 5 %
    pd_ = {"2026-09-01": _cells(**{"1954A": (500, 300), "100T": (500, 0), "200C": (500, 0)}),
           "2026-09-05": _cells(**{"1954A": (500, 0), "100T": (500, 0), "200C": (500, 0)}),
           "2026-09-09": _cells(**{"1954A": (500, 60), "100T": (500, 50), "200C": (500, 40)}),
           "2026-09-12": _cells(**{"1954A": (50, 5)})}
    r = check_in_data(mk, pd_, sorted(pd_))
    assert dict(r["timeline"]) == {"2026-09-01": "mixed", "2026-09-05": "absent",
                                   "2026-09-09": "present", "2026-09-12": "uncovered"}


def test_day_marks_use_the_cells_rule_and_say_when_a_day_is_too_thin():
    from process.variant_explorer import check_in_data
    mk = ["400T", "500A", "600C"]
    pd_ = {
        # one marker decided (absent), the others under 100 reads -> thin
        "2026-09-05": {"400T": [500, 0, 0, 0], "500A": [40, 10, 10, 10]},
        # two decided, they disagree -> mixed
        "2026-09-09": {"400T": [500, 100, 100, 100], "500A": [500, 0, 0, 0]},
    }
    r = check_in_data(mk, pd_, sorted(pd_))
    assert dict(r["timeline"]) == {"2026-09-05": "thin", "2026-09-09": "mixed"}
    # the hover says why, marker by marker
    assert "500A 40 reads (under 100)" in r["day_detail"]["2026-09-05"]
    assert "400T 20 % of 500 reads ✓" in r["day_detail"]["2026-09-09"]


def test_a_day_uses_the_link_test_of_the_whole_window():
    # XFG in Chur: each day has only ~12 reads for the link test, the window has
    # plenty -> the days can be present, as the cell is
    from process.variant_explorer import check_in_data
    mk = ["400T", "500A"]
    pd_ = {d: {"400T": [300, 30, 12, 12], "500A": [300, 25, 12, 12]}
           for d in ("2026-09-01", "2026-09-05", "2026-09-09")}
    r = check_in_data(mk, pd_, sorted(pd_))
    assert [m for _d, m in r["timeline"]] == ["present"] * 3
    # a marker whose reads don't look like the lineage stays undecided
    bad = {d: {"400T": [300, 30, 12, 0], "500A": [300, 25, 12, 12]} for d in pd_}
    assert [m for _d, m in check_in_data(mk, bad, sorted(bad))["timeline"]] == ["thin"] * 3


def test_a_one_marker_lineage_can_be_absent_on_a_day():
    from process.variant_explorer import check_in_data
    r = check_in_data(["400T"], {"2026-09-01": {"400T": [500, 0, 0, 0]}}, ["2026-09-01"])
    assert r["timeline"] == [["2026-09-01", "absent"]]
