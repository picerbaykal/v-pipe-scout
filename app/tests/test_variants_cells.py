"""The Lineages tree's city cells (components.variants_table, 2026-10-20):
the label counts every ★ marker (✓ present · ~ no vote · ✗ absent), the
colour shows the answer only (three colours), "no ★ marker" is said once
across the row, and a panel lineage absent in the chosen city gets one
warning in its row. Pure functions, no Streamlit, no LAPIS."""
import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, "/app_shared")


def _load_variants_table():
    """Load app/components/variants_table.py by its path: tests/components/
    is a package too, and `import components` would find that one first."""
    path = Path(__file__).resolve().parents[1] / "components" / "variants_table.py"
    spec = importlib.util.spec_from_file_location("variants_table_under_test", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


vt = _load_variants_table()


def _d(n_markers, n_measured, n_present, **kw):
    return dict(n_markers=n_markers, n_measured=n_measured, n_present=n_present, **kw)


XFG_TI = _d(2, 1, 1)          # Lugano: 8350C present, 4214C at 2 % (no vote)
BA32_TI = _d(24, 24, 0)       # Lugano: 24 markers absent


def test_label_counts_every_marker():
    assert vt._counts_label(XFG_TI) == "1✓ 1~"       # was "1/1": 4214C hidden
    assert vt._counts_label(BA32_TI) == "24✗"
    assert vt._counts_label(_d(4, 3, 2)) == "2✓ 1~ 1✗"


def test_colour_is_the_answer_only():
    assert vt._scale(_d(3, 3, 3))[0] == vt.S_YES
    assert vt._scale(_d(3, 3, 3), found=True)[0] == vt.S_YES   # was red outside the panel
    assert vt._scale(BA32_TI)[0] == vt.S_NO                     # was orange in the panel
    assert vt._scale(_d(3, 3, 0), found=True)[0] == vt.S_NO


def test_tie_not_called_and_nothing_decided_share_one_style():
    assert vt._scale(_d(2, 2, 1))[0] == vt.S_OPEN               # tie
    assert vt._scale(XFG_TI)[0] == vt.S_OPEN                    # 1 present, needs 2
    st, lab = vt._scale(_d(3, 0, 0))[:2]                         # nothing decided
    assert st == vt.S_OPEN and lab == "3~"


def test_no_light_shade_any_more():
    # weak evidence (no coverage data) keeps the full colour; the hover says it
    d = _d(2, 2, 2, basis="all", state="present")
    assert vt._scale(d)[0] == vt.S_YES
    assert "weak" in vt._cell_tip("XFG", "Lugano", d)


def test_no_marker_said_once_across_the_cities():
    cities = ["Basel (BS)", "Lugano (TI)"]
    pc = {c: _d(0, 0, 0) for c in cities}
    html = vt._panel_cells("KP.2", pc, cities, "Lugano (TI)")
    assert html.count("<td") == 1 and f"colspan='{len(cities)}'" in html
    assert "no ★ marker" in html


def test_panel_lineage_present_nowhere_gets_one_warning():
    # LP.8: absent in 4 cities, all markers without a vote in 2 (Chur, Bern)
    pc = {"BS": _d(4, 4, 0), "GR": _d(4, 0, 0), "GE": _d(4, 4, 0),
          "BE": _d(4, 0, 0), "TI": _d(4, 2, 0), "ZH": _d(4, 4, 0)}
    w = vt._panel_warn(pc, in_panel=True)
    assert "in your panel" in w
    assert "present in none of 6 cities: absent in 4, not called in 2" in w
    assert vt._panel_warn(pc, in_panel=False) == ""            # not in the panel
    # present in one city: no warning, whichever city is selected
    pc["ZH"] = _d(4, 4, 3)
    assert vt._panel_warn(pc, in_panel=True) == ""


def test_legends_describe_the_new_cells():
    for text in (vt.SHORT_LEGEND, vt.legend_full_html()):
        assert "light = weak" not in text
        assert "✓" in text and "~" in text and "✗" in text
