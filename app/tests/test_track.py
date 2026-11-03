"""+ Add needs the weeks in a row AND a share worth tracking in a trial
deconvolution (process.track, 2026-11-03)."""
import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, "/app_shared")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from process.track import size, trial_shares   # noqa: E402

C = ["Basel (BS)", "Lugano (TI)", "Zürich (ZH)"]


def _res(city, variant, vals):
    return {city: {variant: {"timeseriesSummary": [
        {"date": f"2026-10-{d:02d}", "proportion": p} for d, p in vals]}}}


def test_shares_are_the_last_two_weeks_in_well_covered_cities():
    r = {"Lugano (TI)": _res("Lugano (TI)", "ST.1", [(1, 0.0), (5, 0.30), (12, 0.50)]),
         "Chur (GR)": _res("Chur (GR)", "ST.1", [(12, 0.40)])}
    sh = trial_shares(r, ["ST.1"], ["Lugano (TI)"], cfg={"last_days": 8})
    assert set(sh["ST.1"]) == {"Lugano (TI)"}                  # Chur thin: left out
    assert abs(sh["ST.1"]["Lugano (TI)"] - 0.40) < 1e-9        # 5 and 12 Oct, not 1 Oct


def test_size_words():
    assert size({"TI": 0.43})[0] == "track"
    assert size({"TI": 0.017})[0] == "small"
    assert size({"TI": 0.002})[0] == "negligible"
    assert size({}) is None


def _vt():
    path = Path(__file__).resolve().parents[1] / "components" / "variants_table.py"
    spec = importlib.util.spec_from_file_location("vt_track", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _f(trial):
    seen = {"pairs": {"positive": ["2026-10-01", "2026-10-08"],
                      "positive_well": ["2026-10-01", "2026-10-08"],
                      "measured": ["2026-10-01", "2026-10-08"]}, "state": "present",
            "pairs_checked": True}
    wk = {"weeks": [{"end": "2026-10-04", "status": "confirmed"},
                    {"end": "2026-10-11", "status": "confirmed"},
                    {"end": "2026-10-18", "status": "confirmed"}],
          "in_a_row": 3, "need": 3, "ready": True}
    return {"per_city": {c: {"days": 2, "check": seen} for c in C}, "weekly": wk,
            "add_ready": True, "trial": trial}


def test_add_needs_a_share_worth_tracking():
    vt = _vt()
    ov = vt._finding_overview("PJ.2.1", _f({"Lugano (TI)": 0.23}), C, [])[1]
    assert "data-add" in ov and "up to 23 % (TI)" in ov
    ov = vt._finding_overview("SW.2.1", _f({"Lugano (TI)": 0.002}), C, [])[1]
    assert "data-add" not in ov and "Watch" in ov and "negligible: 0.2 %" in ov
    ov = vt._finding_overview("SB.1.1.1", _f({"Lugano (TI)": 0.017}), C, [])[1]
    assert "data-add" not in ov and "small: up to 1.7 %" in ov
    # no trial yet: + Add as before
    assert "data-add" in vt._finding_overview("TA.1", _f(None), C, [])[1]


def test_panel_row_says_under_one_percent():
    vt = _vt()
    d = {"pairs": {"positive": ["2026-10-01"], "positive_well": ["2026-10-01"],
                   "measured": ["2026-10-01"]}, "state": "present", "n_markers": 3,
         "n_present": 3, "n_measured": 3, "pairs_checked": True}
    ov = vt._panel_overview({c: d for c in C}, C, None, {"Lugano (TI)": 0.003})[0]
    assert "under 1 % in every city" in ov
