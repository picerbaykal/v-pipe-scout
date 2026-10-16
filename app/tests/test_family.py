"""process.family: the Lineages table's decisions without Streamlit (2026-10-15)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from process.family import family_rows, plus_add  # noqa: E402

D = "2026-09-01"


def _pd(shares):
    return {D: {m: [1000, int(1000 * v), 0, 0] for m, v in shares.items()}}


def test_plus_add_needs_present_and_more_than_one_sample():
    present = {"check": {"state": "present"}}
    assert plus_add({"ZH": {**present, "days": 2}}, "confirmed", 2)
    assert not plus_add({"ZH": {**present, "days": 1}}, "1 day", 2)
    # one day in each of two cities is enough
    assert plus_add({"ZH": {**present, "days": 1}, "BS": {"check": {}, "days": 1}}, "1 day", 2)
    assert not plus_add({"ZH": {"check": {"state": "mixed"}, "days": 5}}, "confirmed", 2)
    assert not plus_add({"ZH": {**present, "days": 5}}, "family", 2)
    assert not plus_add({"ZH": {**present, "days": 5}}, "named", 2)


def test_family_rows_skips_mixed_children_and_lost_markers():
    # BA.2.86 > JN.1 > LF.7 > XFG > (XFG.14 mixed, XFG.1.1 present); NB.1.8.1 below JN.1
    pv = {
        "BA.2.86": {"markers": ["1A", "2A"], "per_date": _pd({"1A": 1.0, "2A": 1.0})},
        # 3565C lost again below JN.1 (XFG lacks it): JN.1 is read at 27810C only
        "JN.1": {"markers": ["27810C", "3565C"], "per_date": _pd({"27810C": .997, "3565C": .58})},
        "LF.7": {"markers": ["5A"], "per_date": _pd({"5A": .6})},
        "XFG": {"markers": ["6A"], "per_date": _pd({"6A": .45})},
        "XFG.14": {"markers": ["7A"], "per_date": _pd({"7A": .1})},
        "XFG.1.1": {"markers": ["8A"], "per_date": _pd({"8A": .05})},
        "NB.1.8.1": {"markers": ["9A"], "per_date": _pd({"9A": .38})},
    }
    sig = {"BA.2.86": ["1A", "2A"], "JN.1": ["1A", "2A", "27810C", "3565C"],
           "LF.7": ["1A", "2A", "27810C", "3565C", "5A"],
           "XFG": ["1A", "2A", "27810C", "5A", "6A"],
           "XFG.14": ["1A", "2A", "27810C", "5A", "6A", "7A"],
           "XFG.1.1": ["1A", "2A", "27810C", "5A", "6A", "8A"],
           "NB.1.8.1": ["1A", "2A", "27810C", "9A"]}
    ball = {"BA.2.86": ["JN.1", "LF.7", "XFG", "XFG.14", "XFG.1.1", "NB.1.8.1"],
            "JN.1": ["LF.7", "XFG", "XFG.14", "XFG.1.1", "NB.1.8.1"],
            "LF.7": ["XFG", "XFG.14", "XFG.1.1"], "XFG": ["XFG.14", "XFG.1.1"],
            "XFG.14": [], "XFG.1.1": [], "NB.1.8.1": []}
    states = {"BA.2.86": "present", "JN.1": "present", "LF.7": "present", "XFG": "present",
              "XFG.14": "mixed", "XFG.1.1": "present", "NB.1.8.1": "present"}
    findings = {v: {"status": "confirmed", "addable": True,
                    "per_city": {"ZH": {"check": {"state": s, "basis": "all"}}}}
                for v, s in states.items()}
    scan = {"ZH": {"findings_check": {"per_variant": pv, "below_all": ball, "sigs": sig}}}
    family_rows(findings, scan, {"ZH": {}}, [])
    assert findings["BA.2.86"]["status"] == "family"
    assert findings["JN.1"]["status"] == "family"
    # XFG: only XFG.1.1 (5 % of 45 %) is present below it -> itself ~ 89 %, stays
    assert findings["XFG"]["status"] == "confirmed"
    assert abs(findings["XFG"]["itself"]["ZH"][0] - (1 - .05 / .45)) < 1e-6
