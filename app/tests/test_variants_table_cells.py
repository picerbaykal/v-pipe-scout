"""The Lineages table's cell colours (components.variants_table, 2026-10-09):
one fixed colour per answer of the ★ check, light when present in the sum but
clearly present in fewer than evidence.min_days samples on their own."""
import os
import sys

sys.path.insert(0, "/app_shared")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from components import variants_table as vt  # noqa: E402


def _d(p, m, n, marks=None):
    d = {"n_present": p, "n_measured": m, "n_markers": n}
    if marks is not None:
        d["timeline"] = [[f"2026-08-{10 + i:02d}", x] for i, x in enumerate(marks)]
    d["state"] = vt._state(d)
    return d


def test_present_in_two_samples_is_solid():
    st = vt._scale(_d(2, 2, 2, ["present", "thin", "present"]))[0]
    assert st.startswith(f"background:{vt.GREEN}")


def test_present_in_one_sample_only_is_light():
    d = _d(2, 2, 3, ["absent", "thin", "present", "mixed"])
    assert vt._scale(d)[0] == vt.S_PRESENT_LIGHT
    assert vt._scale(d, found=True)[0] == vt.S_FOUND_LIGHT
    assert "only 1 sample on its own (12 Aug)" in vt._clear_note(d)


def test_mixed_and_absent_keep_their_colours():
    assert vt._scale(_d(2, 5, 5, ["present"]))[0] == vt.S_MIXED
    assert vt._scale(_d(0, 4, 4, ["absent"]))[0] == vt.S_NOTFOUND
    assert vt._scale(_d(0, 4, 4, ["absent"]), found=True)[0] == vt.S_PALE


def test_one_present_marker_is_never_present():
    assert vt._state(_d(1, 1, 1)) == "mixed"


# ── trend arrows and the well-covered basis (2026-10-13) ──────────────────
def _tl(marks, basis="all", well=()):
    from process.sample_qc import PROBES, summarize
    good = summarize({p: 9000 for p in PROBES})
    bad = summarize({p: 1 for p in PROBES})
    tl = [[f"2026-08-{10 + i:02d}", m] for i, m in enumerate(marks)]
    return {"n_present": 2, "n_measured": 2, "n_markers": 2, "state": "present",
            "timeline": tl, "basis": basis,
            "quality": {d: (good if i in well else bad) for i, (d, _m) in enumerate(tl)}}


def test_trend_new_gone_and_steady():
    assert vt._trend(_tl(["absent", "absent", "present", "present"]))[0] == "↑"
    assert vt._trend(_tl(["present", "present", "absent", "absent"]))[0] == "↓"
    assert vt._trend(_tl(["present", "present", "present"]))[0] == "→"
    assert vt._trend(_tl(["absent", "present"]))[0] == ""        # one sample is no trend


def test_trend_uses_only_well_covered_samples_when_they_decided():
    d = _tl(["present", "absent", "present", "absent"], basis="well", well=(1, 3))
    assert vt._basis_samples(d)[0] == [("2026-08-11", "absent"), ("2026-08-13", "absent")]
    assert vt._trend(d)[0] == ""                                  # absent throughout


def test_weak_basis_is_light():
    d = {"n_present": 0, "n_measured": 2, "n_markers": 2, "state": "absent", "basis": "all"}
    assert vt._scale(d)[0] == vt.S_NOTFOUND_LIGHT
    d["basis"] = "well"
    assert vt._scale(d)[0] == vt.S_NOTFOUND


def test_deep_scan_evidence_says_when_samples_are_low_coverage():
    from components.variants_table import _ev_quality
    from process.sample_qc import PROBES, summarize
    good = summarize({p: 9000 for p in PROBES})
    bad = summarize({p: (200000 if i < 15 else 1) for i, p in enumerate(PROBES)})
    q = {"2026-08-01": bad, "2026-08-05": bad, "2026-08-09": good}
    d = {"ev_dates": ["2026-08-01", "2026-08-05"], "check": {"quality": q}}
    assert _ev_quality(d) == ", all low coverage"
    d["ev_dates"].append("2026-08-09")
    assert _ev_quality(d) == " (1 well covered)"
    assert _ev_quality({"ev_dates": ["2026-08-01"], "check": {}}) == ""
