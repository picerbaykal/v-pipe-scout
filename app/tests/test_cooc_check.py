"""Unit tests for the co-occurrence check (process.cooc): specific_markers
(which markers) and accumulate_check_stats + check_verdicts (the marker vote).
Pure functions, no LAPIS."""
import sys

sys.path.insert(0, "/app_shared")

from process.cooc import (CHECK_DEFAULTS, accumulate_check_stats,  # noqa: E402
                          check_verdicts, specific_markers)

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
    m = specific_markers(PANEL, LINEAGES, PARENT, dict(CFG, out_other_max=0))
    assert "300G" not in m["XV"]            # parent LF carries it
    assert "800T" not in m["NB"]            # OTHER carries it
    assert "400T" in m["XV"]                # XW is another recombinant (tolerated)


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


def test_confirmed_when_markers_present_on_variant_reads():
    rows = [[_row("a", 900, p400="T", p300="G"), _row("a", 100, p400="C", p300="G"),
             _row("a", 500, p500="A", p300="G")]]
    out = check_verdicts(["400T", "500A"], _tally(rows, ["400T", "500A"]), CFG)
    assert out["verdict"] == "confirmed"
    assert (out["n_present"], out["n_measured"]) == (2, 2)
    assert out["markers"]["400T"]["status"] == "present"


def test_not_found_and_one_convergent_marker_does_not_flip_it():
    rows = [[_row("a", 900, p400="T", p300="G"), _row("a", 1100, p400="C", p300="G"),
             _row("a", 2000, p500="G", p600="T", p700="A", p300="G")]]
    mk = ["400T", "500A", "600C", "700G"]
    out = check_verdicts(mk, _tally(rows, mk), CFG)
    assert (out["n_present"], out["n_measured"]) == (1, 4)
    assert out["verdict"] == "not_found"


def test_real_split_is_inconsistent():
    rows = [[_row("a", 2000, p400="T", p300="G"), _row("a", 2000, p500="G", p300="G")]]
    out = check_verdicts(["400T", "500A"], _tally(rows, ["400T", "500A"]), CFG)
    assert out["verdict"] == "inconsistent"


def test_marker_on_non_variant_reads_is_not_present():
    # 400T seen, but neighbour 300 is reference -> not on a variant haplotype
    rows = [[_row("a", 500, p400="T", p300="A"), _row("a", 500, p400="C", p300="A")]]
    out = check_verdicts(["400T"], _tally(rows, ["400T"]), CFG)
    assert out["markers"]["400T"]["status"] == "unmeasured"
    assert out["verdict"] == "cant_confirm"


def test_counts_are_pooled_over_dates():
    rows = [[_row("a", 60, p400="T", p300="G")], [_row("b", 60, p400="T", p300="G")]]
    out = check_verdicts(["400T"], _tally(rows, ["400T"]), CFG)
    assert out["markers"]["400T"]["cov"] == 120
    assert out["verdict"] == "confirmed"


def test_too_few_reads_or_no_markers_cant_confirm():
    assert check_verdicts([], {}, CFG)["verdict"] == "cant_confirm"
    rows = [[_row("a", 30, p400="T", p300="G")]]
    assert check_verdicts(["400T"], _tally(rows, ["400T"]), CFG)["verdict"] == "cant_confirm"


def test_positions_outside_batch_ignored():
    stats = {}
    accumulate_check_stats([_row("a", 900, p400="T", p500="A")], [500], SIG,
                           {"XV": ["400T"]}, stats)
    assert stats == {}