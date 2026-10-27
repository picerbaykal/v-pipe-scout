"""Its mutations are judged on the samples where its coocs were seen
(process.national.grade, 2026-10-27). Pooled over every sample of every city,
a rare lineage's shares are diluted by the samples without it: PJ.2.1, Sep-Oct,
coocs on 6 of 48 samples, "4 of 18 absent" on all 48."""
import importlib.util
import sys
from datetime import date, timedelta
from pathlib import Path

sys.path.insert(0, "/app_shared")

from process.national import counted_samples, grade      # noqa: E402
from process.variant_explorer import check_in_data       # noqa: E402

MK = ["1000A", "1100C", "5000G", "5100T", "9000A"]
PAIRS = ["1000A+1100C", "5000G+5100T"]
WELL = {"median": 5000, "share_100": 0.9}


def _days(n, start="2026-09-07"):
    a = date.fromisoformat(start)
    return [(a + timedelta(days=i)).isoformat() for i in range(n)]


def _city(days, pos, gone=()):
    """its coocs and 10 % of reads on `pos` only; `gone` at 0 % everywhere"""
    pd = {}
    for d in days:
        h = 100 if d in pos else 0
        c = {m: [1000, 0 if m in gone else h, 0 if m in gone else h, 0 if m in gone else h]
             for m in MK}
        for p in PAIRS:
            c[p] = [1000, h, h, h]
        pd[d] = c
    return check_in_data(MK, pd, days, quality={d: WELL for d in days})


D = _days(24)
POS = {"GE": [D[5], D[12]], "ZH": [D[13], D[20]]}


def test_counted_samples_are_the_ones_with_coocs():
    per = {c: _city(D, ds) for c, ds in POS.items()}
    assert counted_samples(per) == {c: set(ds) for c, ds in POS.items()}


def test_a_rare_lineage_is_not_absent_by_dilution():
    # 10 % on 2 of 24 samples per city = 0.8 % pooled over all 48: "absent"
    # before; on its 4 samples it is 10 %
    per = {c: _city(D, ds) for c, ds in POS.items()}
    g = grade(per)
    assert g["n_counted"] == 4 and g["missing"] == [] and g["pooled_seen"]


def test_a_mutation_absent_on_its_own_samples_still_counts():
    per = {c: _city(D, ds, gone=("9000A",)) for c, ds in POS.items()}
    g = grade(per)
    assert g["missing"] == ["9000A"] and g["n_measured"] == 5


def test_the_row_says_where():
    path = Path(__file__).resolve().parents[1] / "components" / "variants_table.py"
    spec = importlib.util.spec_from_file_location("vt_counted", path)
    vt = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(vt)
    C = ["Genève (GE)", "Zürich (ZH)"]
    gone = ("9000A", "5100T")       # 2 of 5 absent on its samples: over 10 %
    f = {"status": "confirmed", "per_city": {
        c: {"days": 2, "check": _city(D, POS[c[-3:-1]], gone=gone)} for c in C}}
    _kind, ov, *_ = vt._finding_overview("X", f, C, [])
    assert "2 of 5 mutations absent on its 4 samples" in ov


def test_details_say_where_and_when_for_panel_and_found_rows():
    path = Path(__file__).resolve().parents[1] / "components" / "variants_table.py"
    spec = importlib.util.spec_from_file_location("vt_counted2", path)
    vt = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(vt)
    C = ["Genève (GE)", "Zürich (ZH)"]
    per = {c: _city(D, POS[c[-3:-1]], gone=("9000A",)) for c in C}
    html = vt._detail_row("g", "X", per, C, C[0], 8)
    assert "its 4 samples with coocs, oldest first: 12 Sep GE · 19 Sep GE · 20 Sep ZH · 27 Sep ZH" in html
    assert ("its mutations on its 4 samples with coocs, all cities together: 5 measured, "
            "1 absent: 9000A") in html
    assert vt._detail_row("g", "X", per, C, C[0], 8, found=True).count("its 4 samples") == 2
