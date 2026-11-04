"""A recurring unassigned pattern that clinical says is a tree lineage's
(2026-11-04): 25409A 25455T 25528T was ST.1 in 90 of 90 clinical sequences,
in no signature (ST.1's has 25528T only) — a branch of ST.1 with no name yet.
The pattern's row says who carries it, and a "No name yet" row hangs under
ST.1. Pure functions, no Streamlit, no LAPIS."""
import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, "/app_shared")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from process.clinical import carrier, carrier_words  # noqa: E402

CFG = {"branch_min_sequences": 20, "branch_min_share": 0.9}
PARENT = {"ST.1": "XFG.6.5.1", "XFG.6.5.1": "XFG.6.5", "XFG.6.5": "XFG.6", "XFG.6": "XFG",
          "XFG": "", "ST.1.1": "ST.1", "PJ.2.1": "KP.3", "KP.3": ""}


def _vt():
    path = Path(__file__).resolve().parents[1] / "components" / "variants_table.py"
    spec = importlib.util.spec_from_file_location("variants_table_branch", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_deepest_tree_row_that_holds_most_sequences():
    # all ST.1: ST.1, not its ancestor XFG
    assert carrier({"ST.1": 90}, PARENT, {"XFG", "ST.1"}, CFG) == ("ST.1", 90, 90)
    # a child counts for its parent row
    assert carrier({"ST.1": 50, "ST.1.1": 40}, PARENT, {"XFG", "ST.1"}, CFG) == ("ST.1", 90, 90)
    # ST.1 not a row: the nearest row above it
    assert carrier({"ST.1": 90}, PARENT, {"XFG"}, CFG) == ("XFG", 90, 90)


def test_no_carrier_when_mixed_or_too_few():
    assert carrier({"ST.1": 50, "PJ.2.1": 40}, PARENT, {"ST.1", "PJ.2.1"}, CFG) is None
    assert carrier({"ST.1": 19}, PARENT, {"ST.1"}, CFG) is None
    assert carrier({}, PARENT, {"ST.1"}, CFG) is None


def test_words():
    assert carrier_words({"ST.1": 90}, PARENT, {"ST.1"}, CFG) == (
        "clinical: ST.1 (90 of 90 sequences)", "ST.1")
    w, x = carrier_words({"ST.1": 50, "PJ.2.1": 40}, PARENT, {"ST.1", "PJ.2.1"}, CFG)
    assert x is None and w == "clinical: mixed · most ST.1 (50 of 90)"
    assert carrier_words({"ST.1": 7}, PARENT, {"ST.1"}, CFG)[0] == (
        "clinical: 7 sequences, too few (most ST.1, 7)")
    assert carrier_words({}, PARENT, {"ST.1"}, CFG) == ("clinical: no sequences with all of them", None)
    assert carrier_words(None, PARENT, {"ST.1"}, CFG) == ("", None)     # LAPIS not reached


def _build(vt, sig, counts):
    cities = ["Genève (GE)", "Lugano (TI)", "Zürich (ZH)"]
    base = dict(depth=1, guides=[], last=True, has_children=False)
    rows = [dict(base, node="ST.1", kind="panel")]
    ver = {c: {"ST.1": dict(n_markers=5, n_measured=5, n_present=3)} for c in cities}
    muts = ["25409A", "25455T", "25528T"]
    novel = [(muts, {"Genève (GE)": 2, "Lugano (TI)": 3, "Zürich (ZH)": 2})]
    world = {"2026-09-15": {"ST.1": 100, "": 1000}, "2026-09-22": {"ST.1": 30, "": 900}}
    return vt.build(cities, cities[1], rows, ver, {}, ["ST.1"], novel=novel,
                    clinical={"World": world}, pango_parent=PARENT,
                    pattern_clinical=lambda m: counts, signature_of=lambda l: sig)


def test_branch_row_under_its_carrier():
    vt = _vt()
    html = _build(vt, {"25528T", "232C"}, {"ST.1": 90})
    assert "No name yet" in html and "class='branch g-br-ST_1 hid'" in html   # folded
    assert "◇ 1 branch, no name yet" in html and "data-grp='br-ST_1'" in html
    assert "+25409A +25455T" in html and "+25528T" not in html     # only what ST.1 lacks
    assert "clinical: 90 of 130 ST.1 sequences (69 %), no other lineage" in html
    assert "clinical: ST.1 (90 of 90 sequences) → under ST.1 in the tree" in html
    assert "no + Add" in html


def test_no_branch_row_when_mixed_or_all_in_the_signature():
    vt = _vt()
    mixed = _build(vt, {"25528T"}, {"ST.1": 50, "PJ.2.1": 40})
    assert "No name yet" not in mixed and "clinical: mixed" in mixed and "◇" not in mixed
    own = _build(vt, {"25409A", "25455T", "25528T"}, {"ST.1": 90})
    assert "No name yet" not in own
    off = _build(vt, {"25528T"}, None)                             # LAPIS not reached
    assert "No name yet" not in off and "clinical:" not in off


def test_a_lineage_that_cannot_be_checked_still_shows_its_share():
    # 2026-11-04: KP.3 held 9 % in Zürich with 0 ★ markers and the row showed
    # no number at all, so "can't be checked" and the share were never read
    # together. The three branches that can't judge a lineage now carry it.
    vt = _vt()
    cities = ["Genève (GE)", "Zürich (ZH)"]
    tr = {"Zürich (ZH)": 0.092}
    st = {"Zürich (ZH)": {"mean": 0.092, "last": 0.091, "last_date": "2026-10-16",
                          "peak": 0.097, "peak_date": "2026-10-11"}}

    def d(mk, meas, pres, **kw):
        return dict(n_markers=mk, n_measured=meas, n_present=pres, **kw)

    no_marker = {c: d(0, 0, 0) for c in cities}
    html, tip, _dot, _ink, _k = vt._panel_overview(no_marker, cities, None, tr, st)
    assert "can't be checked" in html and "ZH 9.2 % over 2 weeks" in html
    assert "deconvolution share" in tip

    no_pairs = {c: d(3, 0, 0, pairs_checked=True, pairs={}) for c in cities}
    html, _t, _d2, _i, _k2 = vt._panel_overview(no_pairs, cities, None, tr, st)
    assert "can't be checked with coocs" in html and "ZH 9.2 %" in html

    absent = {c: d(4, 4, 0, pairs={"measured": 3, "positive": 0}) for c in cities}
    html, _t, _d3, _i2, _k3 = vt._panel_overview(absent, cities, None, tr, st)
    assert "not seen in any city" in html and "ZH 9.2 %" in html

    # without a trial there is no share line, as before
    html, _t, _d4, _i3, _k4 = vt._panel_overview(no_marker, cities)
    assert "can't be checked" in html and "%" not in html
