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


# ── a sample counts at most cap_reads reads per marker (2026-10-13) ───────
def test_one_huge_sample_cannot_outvote_the_others():
    # Lugano, XFG's 8350C: 51 % of 176 reads on 20 Aug, 0.06 % of 724,267 on 4 Sep
    pd_ = {"2026-08-20": {"8350C": [176, 90, 90, 90]},
           "2026-09-04": {"8350C": [724267, 460, 460, 460]}}
    capped = check_verdicts(["8350C"], pd_, dict(CFG, cap_reads=1000))["markers"]["8350C"]
    summed = check_verdicts(["8350C"], pd_, dict(CFG, cap_reads=0))["markers"]["8350C"]
    assert summed["freq"] < 0.001 and summed["status"] == "absent"
    assert 0.07 < capped["freq"] < 0.09 and capped["status"] == "present"
    assert (capped["cov"], capped["cov_raw"]) == (1176, 724443)


def test_marker_with_no_neighbour_in_reach_is_judged_on_its_share():
    # NB.1.8.1's 2095A: nearest own mutation 695 bases away, so no read can
    # take the link test (link_n 0 with 120 carriers) -> share alone
    pd_ = {"d1": {"2095A": [1000, 120, 0, 0], "823T": [1000, 360, 300, 290]}}
    out = check_verdicts(["2095A", "823T"], pd_, CFG)
    assert out["no_link"] == ["2095A"]
    assert out["markers"]["2095A"]["status"] == "present"
    assert out["verdict"] == "present"
    # some reads do reach a neighbour, too few: still undecided
    pd_ = {"d1": {"2095A": [1000, 120, 5, 5], "823T": [1000, 360, 300, 290]}}
    assert check_verdicts(["2095A", "823T"], pd_, CFG)["markers"]["2095A"]["status"] == "unmeasured"
    # too few carriers to tell: undecided
    pd_ = {"d1": {"2095A": [1000, 10, 0, 0]}}
    assert check_verdicts(["2095A"], pd_, CFG)["no_link"] == []


def test_day_squares_use_the_windows_no_neighbour_markers():
    from process.variant_explorer import check_in_data
    mk = ["2095A", "823T"]
    pd_ = {"2026-09-01": {"2095A": [1000, 120, 0, 0], "823T": [1000, 360, 300, 290]},
           # this day alone has only 12 carriers, the window decides
           "2026-09-05": {"2095A": [100, 12, 0, 0], "823T": [1000, 360, 300, 290]}}
    r = check_in_data(mk, pd_, sorted(pd_))
    assert dict(r["timeline"]) == {"2026-09-01": "present", "2026-09-05": "present"}
    assert any("share only" in x["why"] for x in r["detail"])


# ── lineages below a found lineage: how much of it is itself (14 Oct) ──

def test_lineages_below_all_uses_the_tree_and_recombinants():
    from process.cooc import lineages_below_all
    parents = {"ANC": "", "KID": "ANC", "GRAND": "KID", "XR": "", "XR.1": "XR",
               "OTHER": ""}
    markers = {"ANC": ["1000A", "5000T"], "KID": ["1100G"], "GRAND": ["7000C"],
               "XR": ["9100A"], "XR.1": ["9200G"], "OTHER": ["9000C"]}
    sigs = {"ANC": {"1000A", "5000T"},
            "KID": {"1000A", "1100G"},                    # lost 5000T: still below ANC
            "GRAND": {"1000A", "1100G", "7000C"},
            "XR": {"1100G", "9100A"},                     # recombinant with KID's part
            "XR.1": {"1100G", "9100A", "9200G"},
            "OTHER": {"9000C", "1000A"}}                  # carries 1000A, not in ANC's line
    b = lineages_below_all(list(parents), markers, sigs, parents)
    assert set(b["ANC"]) == {"KID", "GRAND"}               # by the tree only
    assert set(b["KID"]) == {"GRAND", "XR", "XR.1"}        # recombinants carrying 1100G
    assert b["XR"] == ["XR.1"] and b["OTHER"] == [] and b["GRAND"] == []


def test_nearest_below_skips_lineages_not_taken():
    from process.cooc import nearest_below
    ball = {"X": ["A", "A1", "A2", "B"], "A": ["A1", "A2"], "A1": ["A2"]}
    assert nearest_below("X", ball, lambda d: True) == ["A", "B"]
    # A mixed: the walk goes on below it, to A1 (A2 is inside A1)
    assert nearest_below("X", ball, lambda d: d != "A") == ["A1", "B"]
    assert nearest_below("X", ball, lambda d: d == "A2") == ["A2"]


def test_kept_markers_drop_markers_lost_below_with_reads():
    from process.cooc import kept_markers
    ball = {"JN1": ["A", "B"]}
    sigs = {"JN1": {"27810C", "3565C", "18894T"},
            "A": {"27810C", "18894T"},                    # lost 3565C, has reads
            "B": {"27810C", "3565C"}}                     # lost 18894T, absent here
    mk = ["27810C", "3565C", "18894T"]
    assert kept_markers(mk, "JN1", ball, sigs, lambda e: e == "A") == ["27810C", "18894T"]
    assert kept_markers(mk, "JN1", ball, sigs, lambda e: True) == ["27810C"]


def test_lineage_share_and_jaccard():
    from process.cooc import family_jaccard, lineage_share
    pd_ = {"2026-09-01": {"1000A": [1000, 960, 0, 0], "5000T": [1000, 950, 0, 0],
                          "1100G": [1000, 950, 0, 0], "9999C": [50, 50, 0, 0]}}
    cfg = dict(CFG, cap_reads=0)
    sx = lineage_share(["1000A", "5000T"], pd_, cfg=cfg)
    sd = lineage_share(["1100G", "9999C"], pd_, cfg=cfg)  # 9999C: under 100 reads, left out
    assert abs(sx - 0.955) < 1e-9 and abs(sd - 0.95) < 1e-9
    assert abs(family_jaccard(sx, sd) - 0.95 / 0.955) < 1e-9   # ANC itself ≈ 0.5 %
    assert lineage_share(["9999C"], pd_, cfg=cfg) is None
    assert family_jaccard(None, 0.5) is None and family_jaccard(0.2, 0.3) == 1.0
    assert lineage_share(["1000A"], pd_, dates={"2026-09-02"}, cfg=cfg) is None


# ── every sample votes, weighted by its coverage (2026-10-16) ─────────────
def test_sample_cap_follows_median_depth():
    from process.cooc import sample_cap
    c = dict(CFG, cap_reads=1000, cov_factor=10)
    assert sample_cap(c, None, "2026-09-04") == 1000           # no coverage data
    assert sample_cap(c, {}, "2026-09-04") == 1000             # date not measured
    assert sample_cap(c, {"2026-09-04": 2.0}, "2026-09-04") == 20
    assert sample_cap(c, {"2026-09-04": 0.0}, "2026-09-04") == 1   # never 0
    assert sample_cap(c, {"2026-08-24": 1356.0}, "2026-08-24") == 1000


def test_a_jackpot_sample_counts_as_few_reads():
    # Lugano 4 Sep: median depth 2, 724,266 reads at 8350C, 0 % carry it
    pd_ = {"2026-08-20": {"8350C": [176, 90, 90, 90]},
           "2026-09-04": {"8350C": [724266, 0, 0, 0]}}
    c = dict(CFG, cap_reads=1000, cov_factor=10)
    flat = check_verdicts(["8350C"], pd_, c)["markers"]["8350C"]
    wtd = check_verdicts(["8350C"], pd_, c,
                         depth={"2026-08-20": 22.0, "2026-09-04": 2.0})["markers"]["8350C"]
    assert flat["cov"] == 1176                     # 176 + 1,000
    assert wtd["cov"] == 196                       # 176 (under 220) + 20
    assert abs(wtd["freq"] - 90 / 196) < 1e-9 and wtd["cov_raw"] == 724442


# Lugano, 7 Jun – 5 Sep 2026, XFG's two ★ markers (bench_samples.py):
# {date: (median depth, well covered, (reads, % carrying) at 4214C, at 8350C)}
_LUGANO = {
    "2026-06-12": (28, False, (61, 21), (383, 1)),
    "2026-06-28": (2498, True, None, (9, 22)),
    "2026-07-07": (37, False, (13, 54), (112, 8)),
    "2026-07-11": (54, False, (242235, 98), (49, 6)),
    "2026-08-20": (22, False, (23, 83), (176, 51)),
    "2026-08-24": (1356, True, (71030, 0), (1, 0)),
    "2026-09-04": (2, False, (171707, 0), (724266, 0)),
}


def _lugano():
    per_date, quality = {}, {}
    for d, (md, well, a, b) in _LUGANO.items():
        cells = {}
        for m, x in (("4214C", a), ("8350C", b)):
            if x:
                hit = round(x[0] * x[1] / 100)
                cells[m] = [x[0], hit, hit, hit]          # carriers look like XFG nearby
        per_date[d] = cells
        quality[d] = {"median": float(md), "share_100": 0.6 if well else 0.1}
    return per_date, quality


def test_two_well_covered_jackpots_no_longer_decide_alone():
    # until 16 Oct the 2 well-covered samples decided: 4214C 0 % of 71,030 ->
    # XFG "absent" in Lugano, while every other sample shows it
    from process.variant_explorer import check_in_data
    per_date, quality = _lugano()
    r = check_in_data(["4214C", "8350C"], per_date, sorted(per_date), quality=quality)
    assert r["basis"] == "weighted" and r["n_well"] == 2
    assert (r["state"], r["n_present"], r["n_measured"]) == ("present", 2, 2)


def test_without_coverage_data_every_sample_is_capped_flat():
    from process.variant_explorer import check_in_data
    per_date, _q = _lugano()
    r = check_in_data(["4214C", "8350C"], per_date, sorted(per_date), quality=None)
    assert r["basis"] == "all"


def test_lineage_share_uses_the_same_weights():
    from process.cooc import lineage_share
    pd_ = {"2026-08-20": {"8350C": [176, 90, 0, 0]},
           "2026-09-04": {"8350C": [724266, 0, 0, 0]}}
    c = dict(CFG, cap_reads=1000, cov_factor=10)
    flat = lineage_share(["8350C"], pd_, cfg=c)
    wtd = lineage_share(["8350C"], pd_, cfg=c, depth={"2026-08-20": 22.0, "2026-09-04": 2.0})
    assert abs(flat - 90 / 1176) < 1e-9 and abs(wtd - 90 / 196) < 1e-9
