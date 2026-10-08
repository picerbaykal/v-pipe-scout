"""'Panel variant + 1 change' (near-panel) — per read in process.cooc, own
band in process.completeness_composition."""
import sys

sys.path.insert(0, "/app_shared")

from process.completeness_composition import compute_completeness_composition  # noqa: E402
from process.cooc import near_panel_label  # noqa: E402

XFG = {"22792T", "22893G", "22896G", "23021G"}
NB = {"22865T", "22896A", "22930A"}
BA32 = {"22896C", "22899A"}


def test_other_base_at_own_position_is_one_change():
    # XFG read with 22896C instead of 22896G (the Zürich case)
    read = {"22792T", "22893G", "22896C", "23021G"}
    assert near_panel_label(read, set(), {"XFG": XFG, "NB.1.8.1": NB}) == "XFG + 22896C"


def test_result_does_not_depend_on_other_panel_members():
    read = {"22792T", "22893G", "22896C", "23021G"}
    with_ba32 = near_panel_label(read, set(), {"XFG": XFG, "BA.3.2": BA32})
    without = near_panel_label(read, set(), {"XFG": XFG})
    assert with_ba32 == without == "XFG + 22896C"


def test_missing_mutation_is_one_change():
    read = {"22792T", "22893G", "22896G"}
    assert near_panel_label(read, {"23021G"}, {"XFG": XFG}) == "XFG − 23021G"


def test_two_changes_is_not_near():
    read = {"22792T", "22893G", "22896C", "23021G", "22930A"}
    assert near_panel_label(read, set(), {"XFG": XFG, "NB.1.8.1": NB}) == ""


def test_needs_two_shared_mutations():
    assert near_panel_label({"22792T", "99999A"}, set(), {"XFG": XFG}) == ""


def _res(ups):
    return {"dates": ["d1"], "matched_counts": [900], "unexplained_counts": [100],
            "unexplained_patterns": ups}


def test_near_is_its_own_band_not_green():
    ups = [{"date": "d1", "count": 60, "near_count": 50, "near_label": "XFG + 22896C",
            "confirmed_present": ["22792T", "22896C"]}]
    r = compute_completeness_composition(_res(ups), {}, min_reads=1,
                                         panel_union=XFG | BA32)[0]
    assert r["explained"] == 900 and r["near"] == 50 and r["noise"] == 50
    assert abs(r["explained_pct"] + r["near_pct"] + r["noise_pct"] - 1) < 1e-9
    assert r["near_top"] == [("XFG + 22896C", 50)]


def test_legacy_result_keeps_old_rule():
    ups = [{"date": "d1", "count": 60, "confirmed_present": ["22792T", "22896C"]}]
    r = compute_completeness_composition(_res(ups), {}, min_reads=1,
                                         panel_union=XFG | BA32)[0]
    assert r["near"] == 60


def test_change_that_is_a_scanner_marker_goes_to_addable_whoever_is_in_panel():
    # XFG read + 9999A, where 9999A is a marker of a scanner finding. An extinct
    # control carrying 9999A in the panel must not turn it into near/green.
    scan = {"resolved_clade": [{"member_blocks": [{"discriminating": ["9999A", "9998C"]}]}]}
    ups = [{"date": "d1", "count": 60, "near_count": 60, "near_label": "XFG + 9999A",
            "confirmed_present": ["22792T", "22893G", "9999A"]}]
    for union in (XFG, XFG | {"9999A"}):                 # without / with the control
        r = compute_completeness_composition(_res(ups), scan, min_reads=1,
                                             panel_union=union)[0]
        assert r["addable"] == 60 and r["near"] == 0


def test_novel_patterns_have_no_band_their_reads_stay_grey():
    # 2026-10-06: no novel band. Reads of a novel pattern are unexplained and not
    # attributed (grey), and the bands still add up to 100 %.
    scan = {"novel": {"top_patterns": [{"mutations": ["21588T", "21653C", "100A"]}]}}
    ups = [{"date": "d1", "count": 90, "near_count": 0, "near_label": "",
            "confirmed_present": ["21588T", "21653C", "100A"]}]
    r = compute_completeness_composition(_res(ups), scan, min_reads=1,
                                         panel_union={"21618T"})[0]
    assert "novel" not in r and "novel_pct" not in r
    assert r["noise"] == 100 and r["addable"] == 0 and r["near"] == 0   # all unexplained → grey
    assert abs(r["explained_pct"] + r["near_pct"] + r["addable_pct"] + r["noise_pct"] - 1) < 1e-9


def test_lone_mutation_counts_for_a_finding_only_if_star():
    scan = {"resolved_clade": [{"member_blocks": [
        {"discriminating": ["4214C", "4184A"], "mut_star": {"4214C": True, "4184A": False}}]}]}
    star = [{"date": "d1", "count": 50, "confirmed_present": ["4321T", "4214C"],
             "near_count": 0, "near_label": ""}]
    shared = [{"date": "d1", "count": 50, "confirmed_present": ["4321T", "4184A"],
               "near_count": 0, "near_label": ""}]
    assert compute_completeness_composition(_res(star), scan, min_reads=1,
                                            panel_union={"4321T"})[0]["addable"] == 50
    assert compute_completeness_composition(_res(shared), scan, min_reads=1,
                                            panel_union={"4321T"})[0]["addable"] == 0