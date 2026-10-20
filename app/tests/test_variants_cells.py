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
    # present : absent (couldn't answer), like a score
    assert vt._counts_label(XFG_TI) == "1:0 (1)"     # was "1/1": 4214C hidden
    assert vt._counts_label(_d(39, 24, 0)) == "0:24 (15)"   # BA.3.2 Lugano (was 0/24)
    assert vt._counts_label(_d(4, 4, 3)) == "3:1"
    assert vt._counts_label(_d(0, 0, 0)) == "–"


def test_colour_is_the_answer_only():
    assert vt._scale(_d(3, 3, 3))[0] == vt.S_YES
    assert vt._scale(_d(3, 3, 3), found=True)[0] == vt.S_YES   # was red outside the panel
    assert vt._scale(BA32_TI)[0] == vt.S_NO                     # was orange in the panel
    assert vt._scale(_d(3, 3, 0), found=True)[0] == vt.S_NO


def test_tie_not_called_and_nothing_decided_share_one_style():
    assert vt._scale(_d(2, 2, 1))[0] == vt.S_OPEN               # tie
    assert vt._scale(XFG_TI)[0] == vt.S_OPEN                    # 1 present, needs 2
    st, lab = vt._scale(_d(3, 0, 0))[:2]                         # nothing decided
    assert st == vt.S_OPEN and lab == "0:0 (3)"


def test_undecided_cells_are_not_orange():
    # orange is kept for one meaning: a panel lineage no city calls present
    assert vt.AMBER not in vt.S_OPEN and vt.AMBER_T not in vt.S_OPEN


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
    assert "no mutation of its own" in html


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


# ── kept from test_variants_table_cells.py (9–13 Oct; its light colours and
#    trend arrows are gone, these three checks still apply) ────────────────
def _marks(p, m, n, marks):
    d = _d(n, m, p)
    d["timeline"] = [[f"2026-08-{10 + i:02d}", x] for i, x in enumerate(marks)]
    d["state"] = vt._state(d)
    return d


def test_hover_says_on_how_few_samples_it_was_clear():
    # what the light shade said until 20 Oct, now in the hover
    d = _marks(2, 2, 3, ["absent", "thin", "present", "mixed"])
    assert "only 1 sample on its own (12 Aug)" in vt._clear_note(d)
    assert "only 1 sample on its own" in vt._cell_tip("X", "Basel", d)


def test_old_results_use_only_their_well_covered_samples():
    from process.sample_qc import PROBES, summarize
    good = summarize({p: 9000 for p in PROBES})
    bad = summarize({p: 1 for p in PROBES})
    tl = [[f"2026-08-{10 + i:02d}", m] for i, m in enumerate(["present", "absent", "present", "absent"])]
    d = {"n_present": 2, "n_measured": 2, "n_markers": 2, "state": "present", "timeline": tl,
         "basis": "well", "quality": {x: (good if i in (1, 3) else bad) for i, (x, _m) in enumerate(tl)}}
    assert vt._basis_samples(d)[0] == [("2026-08-11", "absent"), ("2026-08-13", "absent")]


def test_deep_scan_evidence_says_when_samples_are_low_coverage():
    from process.sample_qc import PROBES, summarize
    good = summarize({p: 9000 for p in PROBES})
    bad = summarize({p: (200000 if i < 15 else 1) for i, p in enumerate(PROBES)})
    q = {"2026-08-01": bad, "2026-08-05": bad, "2026-08-09": good}
    d = {"ev_dates": ["2026-08-01", "2026-08-05"], "check": {"quality": q}}
    assert vt._ev_quality(d) == ", all low coverage"
    d["ev_dates"].append("2026-08-09")
    assert vt._ev_quality(d) == " (1 well covered)"
    assert vt._ev_quality({"ev_dates": ["2026-08-01"], "check": {}}) == ""


# ── rows by importance: the overview column (2026-10-20) ─────────────────
CITIES = ["Basel (BS)", "Chur (GR)", "Genève (GE)", "Laupen (BE)", "Lugano (TI)", "Zürich (ZH)"]


def _sc(p, a, u=0):
    """a check result from a score p:a (u)"""
    return _d(p + a + u, p + a, p)


def _find(per_city, **kw):
    """per_city: {code: (days, (p, a, u), xcheck)}"""
    pc = {}
    for c in CITIES:
        code = c.split("(")[1][:2]
        if code in per_city:
            days, sc, x = per_city[code]
            pc[c] = {"days": days, "check": _sc(*sc), "xcheck": x}
    return dict({"status": "confirmed", "per_city": pc, "addable": False}, **kw)


def test_panel_overview_seen_not_seen_cant_check():
    nb = {c: s for c, s in zip(CITIES, [_sc(3, 1), _sc(1, 0, 3), _sc(2, 2), _sc(4, 0),
                                         _sc(3, 0, 1), _sc(3, 0, 1)])}
    html, _tip, dot, _ink, _k = vt._panel_overview(nb, CITIES)
    assert "seen in 4 of 6" in html and dot == vt.GREEN
    lp8 = {c: s for c, s in zip(CITIES, [_sc(0, 4), _sc(0, 0, 4), _sc(0, 4), _sc(0, 0, 4),
                                          _sc(0, 2, 2), _sc(0, 4)])}
    html, tip, dot, _ink, _k = vt._panel_overview(lp8, CITIES)
    assert "not seen in any city" in html and dot == vt.AMBER
    assert "absent in 4, not called in 2" in tip
    kp2 = {c: _d(0, 0, 0) for c in CITIES}
    html, _tip, dot, _ink, _k = vt._panel_overview(kp2, CITIES)
    assert "can't be checked" in html and dot != vt.AMBER       # amber means one thing only


def test_suggestion_new_possible_folded():
    # RV.2.1.1: seen in BE and ZH, but suggested on 1 sample only -> possible
    rv = _find({"BE": (1, (4, 0, 2), False), "ZH": (0, (5, 1), True), "GE": (0, (1, 3, 2), True)})
    assert vt._finding_status("RV.2.1.1", rv, CITIES) == (
        "possible", "seen in BE, ZH · suggested on 1 sample only")
    # RF.5.4: suggested 5x in ZH, where its one mutation couldn't answer -> stays
    rf = _find({"ZH": (5, (0, 0, 1), False), "BS": (0, (0, 1), True), "GE": (0, (0, 1), True)})
    assert vt._finding_status("RF.5.4", rf, CITIES) == (
        "possible", "couldn't check in ZH · suggested 5× in ZH")
    # ST.1: not seen where suggested (BS), but a tie in TI -> stays
    st = _find({"BS": (2, (0, 4, 1), False), "TI": (0, (1, 1, 3), True), "GE": (0, (0, 5), True)})
    assert vt._finding_status("ST.1", st, CITIES) == ("possible", "unclear in TI · suggested 2× in BS")
    # PJ.2.1: not seen where suggested, not seen or unclear anywhere -> folded
    pj = _find({"BE": (1, (0, 13, 6), False), "GE": (0, (3, 15, 1), True), "ZH": (0, (0, 19), True)})
    assert vt._finding_status("PJ.2.1", pj, CITIES)[0] == "folded"
    assert vt._finding_status("PJ.2.1", pj, CITIES, has_children=True)[0] == "possible"
    # + Add -> new; family -> covered
    assert vt._finding_status("PQ.16.1.1", _find({}, addable=True), CITIES)[0] == "new"
    assert vt._finding_status("SW.2", _find({}, status="family"), CITIES)[0] == "covered"


def test_table_has_an_overview_column_and_folds_unsupported_suggestions():
    base = dict(depth=1, guides=[], last=True, has_children=False)
    rows = [dict(base, node="LP.8", kind="panel"), dict(base, node="PJ.2.1", kind="finding"),
            dict(base, node="RV.2.1.1", kind="finding")]
    ver = {c: {"LP.8": _sc(0, 4)} for c in CITIES}
    findings = {
        "PJ.2.1": _find({"BE": (1, (0, 13, 6), False), "GE": (0, (3, 15, 1), True)}),
        "RV.2.1.1": _find({"BE": (1, (4, 0, 2), False), "ZH": (0, (5, 1), True)})}
    html = vt.build(CITIES, CITIES[3], rows, ver, findings, ["LP.8"])
    assert "Overview · all cities" in html
    assert "not seen in any city" in html                      # LP.8, amber
    assert "Possible" in html and "suggested on 1 sample only" in html   # RV.2.1.1
    assert "1 suggestion not supported by its own mutations" in html and "PJ.2.1" in html
    assert html.count("class='nm'") == 2                       # PJ.2.1 has no row


def test_legends_describe_the_new_cells():
    for text in (vt.SHORT_LEGEND, vt.legend_full_html()):
        assert "light = weak" not in text
        assert "couldn't answer" in text
