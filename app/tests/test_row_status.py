"""The row's status from its weekly grades (process.national.row_status,
2026-10-27): add / watch / fading / earlier, judged on the last weeks whatever
range is chosen; the cells show the whole range. Weeks patterns from the
summer scorecard (7 Jun - 5 Oct)."""
import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, "/app_shared")

from process.national import row_status      # noqa: E402

ENDS = ["2026-07-05", "2026-07-12", "2026-07-19", "2026-07-26", "2026-08-02", "2026-08-09",
        "2026-08-16", "2026-08-23", "2026-08-30", "2026-09-06", "2026-09-13", "2026-09-20",
        "2026-09-27", "2026-10-04"]
S = {"C": "confirmed", "F": "watch", ".": "none"}


def _wk(pattern):
    weeks = [{"end": e, "status": S[x], "grade": {"cities": []}} for e, x in zip(ENDS, pattern)]
    n = 0
    for w in reversed(weeks):
        if w["status"] != "confirmed":
            break
        n += 1
    return {"weeks": weeks, "in_a_row": n, "need": 3, "ready": n >= 3, "short": False}


def test_add_watch_fading_earlier():
    assert row_status(_wk("......FFCCCCCC"))["kind"] == "add"          # TA.1
    assert row_status(_wk(".......FFFFFCC"))["kind"] == "watch"        # TK.1, 2 of 3
    assert row_status(_wk(".......FFFFFFC"))["kind"] == "watch"        # PJ.2.1
    assert row_status(_wk("FFFF..FFFFFFFF"))["kind"] == "watch"        # never confirmed, found now
    r = row_status(_wk("FCFFFF..FFFFF."))                              # PY.1.1.1
    assert r["kind"] == "fading" and r["last_confirmed"] == "2026-07-12"
    r = row_status(_wk("FFFF....FFFF.."))                              # XFZ.4.1.2
    assert r["kind"] == "earlier" and r["last_found"] == "2026-09-20"
    assert row_status({"weeks": []}) is None


def _vt():
    path = Path(__file__).resolve().parents[1] / "components" / "variants_table.py"
    spec = importlib.util.spec_from_file_location("vt_status", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_weeks_strip_keeps_the_newest_12_and_sums_up_the_rest():
    vt = _vt()
    html = vt._weeks_strip(_wk("FCFFFF..FFFFF."))
    assert html.count("<i ") == 12 and "⋯" in html
    assert "2 older weeks: 1 confirmed, last confirmed 12 Jul" in html


def test_rows_say_fading_or_earlier_only():
    vt = _vt()
    C = ["Genève (GE)", "Zürich (ZH)"]
    seen = {"pairs": {"positive": ["2026-07-10"], "positive_well": ["2026-07-10"],
                      "measured": ["2026-07-10"]}, "state": "present"}
    f = {"status": "confirmed", "per_city": {C[1]: {"days": 2, "check": seen}},
         "weekly": _wk("FCFFFF..FFFFF.")}
    ov = vt._finding_overview("PY.1.1.1", f, C, [])[1]
    assert "<b>Fading</b> · last confirmed 12 Jul" in ov and "data-add" not in ov
    ov = vt._finding_overview("XFZ.4.1.2", dict(f, weekly=_wk("FFFF....FFFF..")), C, [])[1]
    assert "<b>Earlier only</b> · last seen 20 Sep" in ov


def test_panel_row_says_not_confirmed_since():
    vt = _vt()
    C = ["Genève (GE)", "Zürich (ZH)"]
    d = {"pairs": {"positive": ["2026-07-10"], "positive_well": ["2026-07-10"],
                   "measured": ["2026-07-10"]}, "state": "present", "n_markers": 3,
         "n_present": 3, "n_measured": 3}
    ov = vt._panel_overview({c: d for c in C}, C, _wk("FCFFFF..FFFFF."))[0]
    assert "not confirmed since 12 Jul" in ov
    ov = vt._panel_overview({c: d for c in C}, C, _wk("FFFF....FFFFFC"))[0]
    assert "· not confirmed" not in ov


def test_watch_says_how_many_weeks_in_a_row():
    vt = _vt()
    C = ["Genève (GE)", "Zürich (ZH)"]
    seen = {"pairs": {"positive": ["2026-09-30"], "positive_well": ["2026-09-30"],
                      "measured": ["2026-09-30"]}, "state": "present"}
    f = {"status": "confirmed", "per_city": {C[1]: {"days": 2, "check": seen}},
         "weekly": _wk(".......FFFFFCC")}
    assert "<b>Watch</b> <span style='color:#4b5563'>· 2 of 3 weeks</span>" in vt._finding_overview("TK.1", f, C, [])[1]
