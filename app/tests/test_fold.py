"""Found lineages without coocs of their own fold under the tree
(variants_table._finding_status, 2026-10-27); coocs on one sample with most
of its mutations missing too (BA.3.2 Genève July: 32 of 32; XFW.1.1.1: 7 of 9)."""
import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, "/app_shared")

from process.variant_explorer import check_in_data     # noqa: E402

MK = ["1000A", "1100C", "5000G", "5100T"]
PAIRS = ["1000A+1100C", "5000G+5100T"]
WELL = {"median": 5000, "share_100": 0.9}
C = ["Genève (GE)", "Zürich (ZH)"]


def _vt():
    path = Path(__file__).resolve().parents[1] / "components" / "variants_table.py"
    spec = importlib.util.spec_from_file_location("vt_fold", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _check(days, pos=(), zero=()):
    """samples on `days`; coocs on the days in `pos`; markers in `zero` at 0 %"""
    pd = {}
    for d in days:
        h = 300 if d in pos else 0
        c = {m: [1000, 0 if m in zero else h, 0 if m in zero else h, 0 if m in zero else h]
             for m in MK}
        for p in PAIRS:
            c[p] = [1000, h, h, h]
        pd[d] = c
    return check_in_data(MK, pd, list(days), quality={d: WELL for d in days})


def _f(**per):
    return {"status": "confirmed", "per_city": {c: {"days": 2, "check": r} for c, r in per.items()}}


D = ["2026-09-15", "2026-09-22", "2026-09-29"]


def test_coocs_in_a_city_is_new():
    vt = _vt()
    f = _f(**{C[1]: _check(D, pos=D[1:]), C[0]: _check(D)})
    assert vt._finding_status("X", f, C) == ("new", "")


def test_one_sample_with_most_mutations_missing_folds():
    vt = _vt()
    # coocs on 1 sample only; 3 of its 4 mutations at 0 % pooled
    f = _f(**{C[1]: _check(D, pos=D[2:], zero=("1100C", "5000G", "5100T")), C[0]: _check(D)})
    kind, why = vt._finding_status("X", f, C)
    assert kind == "folded" and why.startswith("coocs on 1 sample in ZH")


def test_no_coocs_folds_with_the_cell_words():
    vt = _vt()
    f = _f(**{C[1]: _check(D), C[0]: _check(D)})
    kind, why = vt._finding_status("X", f, C)
    assert kind == "folded" and "not seen there" in why
    assert vt._finding_status("X", f, C, has_children=True)[0] == "suggested"


def test_half_absent_on_one_sample_folds():
    vt = _vt()
    # XFG.6.3.3, 27 Oct: 2 of 4 absent on one Zürich sample -> folded
    f = _f(**{C[1]: _check(D, pos=D[2:], zero=("5000G", "5100T")), C[0]: _check(D)})
    kind, why = vt._finding_status("X", f, C)
    assert kind == "folded" and why.endswith("2 of 4 mutations absent")
