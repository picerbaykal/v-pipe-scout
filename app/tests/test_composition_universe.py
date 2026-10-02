"""Deep scan: a finding found on positions the graph didn't read still colours
the graph, matched on the positions the graph did read."""
from process.completeness_composition import compute_completeness_composition


def _cooc():
    return {"dates": ["2026-08-01"], "matched_counts": [5000], "unexplained_counts": [1000],
            "unexplained_patterns": [{"date": "2026-08-01", "count": 300,
                                      "confirmed_present": ["100A", "200C", "300G"],
                                      "confirmed_absent": ["400T"]}]}


def test_finding_with_unseen_positions_still_matches():
    # 3 of its 6 discriminating mutations were never read by phase 1
    scan = {"resolved_clade": [{"member_blocks": [{"discriminating":
            ["100A", "200C", "300G", "900A", "901C", "902G"]}]}]}
    r = compute_completeness_composition(_cooc(), scan)[0]
    assert r["addable"] == 300


def test_finding_with_under_two_seen_mutations_is_dropped():
    scan = {"resolved_clade": [{"member_blocks": [{"discriminating": ["100A", "900A", "901C"]}]}]}
    r = compute_completeness_composition(_cooc(), scan)[0]
    assert r["addable"] == 0 and r["noise"] == 1000