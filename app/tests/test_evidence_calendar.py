"""Step 7 (2026-10-02): one measure for every lineage cell, a day calendar in
Evidence for every row, clues on novel rows."""
import importlib.util
from pathlib import Path

# app/tests/components/ (a test package) shadows app/components/ on the test
# path, so load the module from its file — tidy in the test cleanup (step 8)
_spec = importlib.util.spec_from_file_location(
    "variants_table_under_test",
    Path(__file__).resolve().parents[1] / "components" / "variants_table.py")
vt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vt)

from process.scanner import _novel_clues, _novel_groups

DAY = {"2026-08-01": 10_000, "2026-08-05": 10_000, "2026-08-09": 10_000}


def test_phrase_present_latest_and_gone():
    tl = [["2026-08-01", "absent"], ["2026-08-05", "present"], ["2026-08-09", "present"]]
    assert vt._cal_phrase(tl) == "present 2 of 3 days · in the latest sample"
    tl = [["2026-08-01", "present"], ["2026-08-05", "absent"], ["2026-08-09", "uncovered"]]
    assert vt._cal_phrase(tl) == "present 1 of 2 days · last 1 Aug"
    assert vt._cal_phrase([["2026-08-01", "uncovered"]]).startswith("not covered")
    assert vt._cal_phrase([["2026-08-01", "absent"]]) == "absent all 1 day"


def test_strip_shows_the_last_covered_samples_only():
    tl = [[f"2026-08-{i:02d}", m] for i, m in enumerate(
        ["present"] * 6 + ["uncovered", "absent", "uncovered", "absent", "absent", "present",
                           "absent"], start=1)]
    cal = vt._day_cal(tl, vt.GREEN)
    assert cal.count("<i style='background:#d1d5db;'") == 4      # 4 absent of the last 5
    assert "uncovered" not in cal and "too few" not in cal


def test_recent_vote():
    from process.variant_explorer import check_in_data
    old = {"1A": [500, 400, 400, 400]}
    gone = {"1A": [500, 0, 0, 0]}
    pdd = {f"2026-07-{i:02d}": old for i in range(1, 11)}
    pdd.update({f"2026-08-{i:02d}": gone for i in range(1, 6)})
    whole = check_in_data(["1A"], pdd, list(pdd))
    now = check_in_data(["1A"], pdd, list(pdd), recent=5)
    assert whole["state"] == "present" and now["state"] == "absent"
    assert now["recent_dates"] == [f"2026-08-{i:02d}" for i in range(1, 6)]


def test_found_cell_uses_the_check():
    chk = {"n_markers": 3, "n_measured": 3, "n_present": 3, "state": "present",
           "timeline": [["2026-08-01", "present"]]}
    f = {"per_city": {"Zürich (ZH)": {"days": 4, "stars": ["1A"], "regions": [(1, 300)],
                                      "check": chk}}}
    cells = vt._finding_cells("PY.1.1.1", f, ["Zürich (ZH)"], "Zürich (ZH)")
    assert ">3/3<" in cells
    ev = vt._finding_evidence(f, "Zürich (ZH)")
    assert "present 1 of 1 day" in ev and "class='cal'" in ev


def test_older_results_still_show_days():
    f = {"per_city": {"Zürich (ZH)": {"days": 4, "stars": ["1A"], "regions": [(1, 300)]}}}
    assert ">4<" in vt._finding_cells("PY.1.1.1", f, ["Zürich (ZH)"], "Zürich (ZH)")


def test_novel_timeline_marks():
    pats = [{"mutations": ["30000A", "30010C"], "count": 100, "date": "2026-08-01"},
            {"mutations": ["30000A", "30010C"], "count": 5, "date": "2026-08-05"}]
    cov = {"2026-08-01": {"30000": 500, "30010": 500}, "2026-08-05": {"30000": 500, "30010": 500},
           "2026-08-09": {"30000": 500, "30010": 40}}
    g = _novel_groups(pats, DAY, sorted(DAY), cov)[0]
    assert [m for _d, m in g["timeline"]] == ["day", "seen", "uncovered"]


def test_clue_real_mutation_plus_misread():
    # 11851T (RF.6's) on 10,000 reads covering 11929; 11929A on 100 of them
    pats = [{"mutations": ["11851T", "11929A"], "count": 100, "date": "2026-08-01"},
            {"mutations": ["11851T", "12071A"], "count": 80, "date": "2026-08-01"}]
    obs = [(frozenset({"11851T", "11929A"}), frozenset(), 100, "2026-08-01"),
           (frozenset({"11851T"}), frozenset({"11929A"}), 9_900, "2026-08-01")]
    groups = [{"mutations": ["11851T", "11929A"]}]

    class _T:
        def dominant_clade(self, cs, f):
            return None
    _novel_clues(groups, pats, obs, {"11851T": {"RF.6"}}, _T())
    c = groups[0]["clue"]
    assert c["anchor"] == "11851T" and c["lineage"] == "RF.6" and c["partners"] == 2
    assert c["n_lineages"] == 1
    assert abs(c["together"] - 0.01) < 1e-6
    ev = vt._novel_evidence(["11851T", "11929A"], {"timeline": [], "clue": c})
    assert "RF.6" in ev and "likely a real mutation + misread" in ev