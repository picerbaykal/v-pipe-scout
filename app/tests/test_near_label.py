"""'± 1 change' label when several panel variants are one change away
(2026-10-09): pango decides, else the dominant variant, else alphabetical —
never the panel's order."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from process.cooc import near_panel_label  # noqa: E402

# the read carries A and B (both variants have them) and C (neither has it)
SIGS = {"XFG": {"100A", "200C", "300G"}, "NB.1.8.1": {"100A", "200C", "400T"}}
READ, ABSENT = {"100A", "200C", "999T"}, set()


def _label(sigs, **kw):
    return near_panel_label(READ, ABSENT, sigs, **kw)


def test_panel_order_no_longer_decides():
    flipped = dict(reversed(list(SIGS.items())))
    assert _label(SIGS) == _label(flipped) == "NB.1.8.1 + 999T"      # alphabetical


def test_pango_decides_when_one_variant_has_a_sublineage_with_the_change():
    stats = {}
    got = _label(SIGS, in_tree=lambda v, sign, m: v == "XFG" and m == "999T", stats=stats)
    assert got == "XFG + 999T" and stats == {"ties": 1, "pango": 1}


def test_dominant_variant_when_pango_does_not_decide():
    stats = {}
    got = _label(SIGS, in_tree=lambda *a: False, support={"XFG": 900, "NB.1.8.1": 40},
                 stats=stats)
    assert got == "XFG + 999T" and stats == {"ties": 1, "dominant": 1}


def test_equal_support_falls_back_to_alphabetical_and_is_counted():
    stats = {}
    got = _label(SIGS, support={"XFG": 10, "NB.1.8.1": 10}, stats=stats)
    assert got == "NB.1.8.1 + 999T" and stats == {"ties": 1, "unresolved": 1}


def test_more_shared_mutations_still_win_first():
    sigs = dict(SIGS, XFG={"100A", "200C", "300G", "777A"})
    read = {"100A", "200C", "777A", "999T"}
    assert near_panel_label(read, set(), sigs) == "XFG + 999T"   # 3 shared beats 2


def test_missing_mutation_label_unchanged():
    sigs = {"XFG": {"100A", "200C", "300G"}}
    assert near_panel_label({"100A", "200C"}, {"300G"}, sigs) == "XFG − 300G"


def test_annotation_uses_the_days_dominant_variant():
    import pandas as pd
    from process.cooc import annotate_cooc_dataframe
    amp = {100: ["A"], 200: ["C"], 300: ["G"], 400: ["T"], 999: ["T"]}
    pos = [100, 200, 300, 400, 999]
    rows = [  # 800 reads only XFG explains, 50 only NB.1.8.1, then the tied read
        {"date": "d1", "count": 800, "[100]": "A", "[200]": "C", "[300]": "G", "[400]": "N", "[999]": "N"},
        {"date": "d1", "count": 50, "[100]": "A", "[200]": "C", "[300]": "N", "[400]": "T", "[999]": "N"},
        {"date": "d1", "count": 30, "[100]": "A", "[200]": "C", "[300]": "N", "[400]": "N", "[999]": "T"},
    ]
    stats = {}
    out = annotate_cooc_dataframe(pd.DataFrame(rows), pos, amp, SIGS, stats=stats)
    assert out.loc[out["near"] != "", "near"].tolist() == ["XFG + 999T"]
    assert stats == {"ties": 1, "dominant": 1}
