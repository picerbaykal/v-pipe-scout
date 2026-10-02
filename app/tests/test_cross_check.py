"""Phase 3, cross-check: a lineage named in another city gets the ★ marker
check here, shown with a dashed border (components.variants_table)."""
import importlib.util
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "variants_table_xcheck",
    Path(__file__).resolve().parents[1] / "components" / "variants_table.py")
vt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vt)

CHK = {"n_markers": 3, "n_measured": 1, "n_present": 1, "state": "present",
       "timeline": [["2026-08-03", "absent"], ["2026-08-10", "present"]]}
F = {"per_city": {
    "Zürich (ZH)": {"days": 3, "stars": ["25461G"], "regions": [(1, 300)], "check": CHK},
    "Basel (BS)": {"days": 0, "check": CHK, "xcheck": True, "found_in": ["Zürich (ZH)"]}}}


def test_cross_checked_cell_is_dashed_with_the_vote():
    cells = vt._finding_cells("PQ.16.1.1", F, ["Basel (BS)"], "Basel (BS)")
    assert ">1/1<" in cells and "dashed #b91c1c" in cells and "found in Zürich" in cells


def test_cross_checked_evidence_says_so():
    ev = vt._finding_evidence(F, "Basel (BS)")
    assert "cross-check" in ev and "present 1 of 2 days" in ev