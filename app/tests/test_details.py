"""The row's details, reworked (2026-10-27): links, one column per week
(wastewater, Europe / World clinical, each city's samples), its coocs seen per
city, its mutations x city; clinical counts from process.clinical."""
import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, "/app_shared")

import test_counted_samples as cs                  # noqa: E402
from process.clinical import family, lineage_weeks  # noqa: E402


def _vt():
    path = Path(__file__).resolve().parents[1] / "components" / "variants_table.py"
    spec = importlib.util.spec_from_file_location("vt_details", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


C = ["Genève (GE)", "Zürich (ZH)"]


def _per():
    return {c: cs._city(cs.D, cs.POS[c[-3:-1]], gone=("9000A",)) for c in C}


def test_clinical_counts_include_the_lineage_below():
    parent = {"TK.1": "RV.2.1.1", "RV.2.1.1": "RV.2.1", "TK.1.1": "TK.1"}
    assert family("RV.2.1.1", parent) == {"RV.2.1.1", "TK.1", "TK.1.1"}
    weekly = {"2026-09-07": {"": 50, "TK.1": 3, "RV.2.1.1": 1, "XFG": 20}}
    assert lineage_weeks(weekly, family("RV.2.1.1", parent), ["2026-09-07", "2026-09-14"]) == {
        "2026-09-07": (4, 50), "2026-09-14": (0, 0)}
    assert lineage_weeks(None, {"X"}, ["2026-09-07"]) is None


def test_details_have_weeks_coocs_and_mutations():
    vt = _vt()
    clin = {"Europe": {"2026-09-07": {"": 52, "X": 6}}, "World": {"2026-09-07": {"": 350, "X": 58}}}
    html = vt._detail_row("g", "X", _per(), C, C[0], 8, found=True, clin=clin, parent={})
    # weeks: 7 Sep to 28 Sep, wastewater + 2 clinical rows + 2 cities
    for w in ("7 Sep", "14 Sep", "21 Sep", "28 Sep"):
        assert f"<th>{w}</th>" in html
    assert "wastewater, all cities" in html and "Europe clinical" in html and "World clinical" in html
    assert "6/52" in html and "58/350" in html
    assert "Its coocs seen" in html and "1000A + 1100C" in html
    assert "X and below / all sequences" in html and "well covered" in html
    assert ">clinical</th>" in html and html.count("together, clinical") == 1      # the top line only
    assert "Its mutations, share of counted reads" in html and "9000A" in html
    assert "nucMutations=" in html and "/Europe/" in html
    # samples: coocs on 12 Sep in Genève
    assert "Genève, 12 Sep: coocs" in html


def test_without_clinical_data_the_rows_are_left_out():
    vt = _vt()
    html = vt._detail_row("g", "X", _per(), C, C[0], 8)
    assert "wastewater, all cities" in html and "Europe clinical" not in html


def test_coocs_table_counts_every_pair_and_marks_common_partners():
    vt = _vt()
    d = {"pairs": {"per_date_all": {"2026-09-12": [["24469A+24604G", 0.66, 300, False],
                                                    ["6183G+6196A", 0.98, 900, True]],
                                     "2026-09-19": [["6183G+6196A", 0.90, 500, True]]}}}
    html = vt._cooc_table({C[0]: d}, [C[0]])
    # 2026-10-27: one row per stretch; a stretch whose pairs all have a
    # backbone mutation is marked; pairs in the hover
    assert "6,183–6,196 · 1 pair" in html and "2 · 98 %" in html
    assert "24,469–24,604 · 1 pair · with a backbone mutation" in html
    assert "24469A + 24604G (backbone mutation)" in html


def test_pairs_closer_than_a_read_are_one_stretch():
    vt = _vt()
    d = {"pairs": {"per_date_all": {"2026-09-12": [["22792T+22896G", 1.0, 900, True],
                                                    ["22893G+23021G", 1.0, 900, True],
                                                    ["23021G+23039G", 1.0, 900, True],
                                                    ["8350C+8393A", 0.99, 900, True]]}}}
    html = vt._cooc_table({C[0]: d}, [C[0]])
    assert "22,792–23,039 · 3 pairs" in html and "8,350–8,393 · 1 pair" in html


def test_couldnt_answer_says_why():
    vt = _vt()
    base = {"status": "unmeasured", "freq": 0.3, "why": ""}
    assert vt._why_short(dict(base, cov=65, cov_raw=142546)) == "<100"
    assert vt._why_short(dict(base, cov=500, cov_raw=500, freq=0.02)) == "1–5 %"
    assert vt._why_short(dict(base, cov=500, cov_raw=500, link_n=8)) == "link"
