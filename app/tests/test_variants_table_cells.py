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
