"""Novel groups: a position with several new bases is flagged as a likely
sequencing error (scanner._novel_groups)."""
from process.scanner import _novel_groups, _multi_alt_positions

DAY = {"2026-08-12": 100_000, "2026-08-20": 100_000}


def _p(muts, n, d="2026-08-12"):
    return {"mutations": muts, "count": n, "date": d}


def test_multi_alt_flags_error_hotspot():
    pats = [_p(["17857A", "17859C"], 3000), _p(["17857G", "17859C"], 2000),
            _p(["17857C", "17859C"], 1500), _p(["17857A", "17859C"], 3000, "2026-08-20")]
    assert _multi_alt_positions(pats, DAY) == {17857: ["A", "C", "G"]}
    g = _novel_groups(pats, DAY)
    assert all(x["likely_error"] for x in g)
    assert g[0]["error_positions"] == {"17857": ["A", "C", "G"]}


def test_single_base_is_not_an_error():
    pats = [_p(["23139A", "23148C"], 3000), _p(["23139A", "23148C"], 2500, "2026-08-20")]
    g = _novel_groups(pats, DAY)
    assert g[0]["likely_error"] is False and g[0]["days"] == ["2026-08-12", "2026-08-20"]


def test_a_stray_second_base_below_the_day_rule_does_not_count():
    # 5 reads of 29279C never pass the day rule (>= 20 reads and >= 0.5 %)
    pats = [_p(["29272T", "29279T"], 3000), _p(["29272T", "29279C"], 5)]
    assert _multi_alt_positions(pats, DAY) == {}


def test_errors_do_not_crowd_out_real_groups():
    pats = [_p([f"{17000 + i}A", "17859C"], 5000) for i in range(25)]
    pats += [_p([f"{17000 + i}G", "17859C"], 5000) for i in range(25)]
    pats += [_p(["23139A", "23148C"], 100)]          # small but real
    g = _novel_groups(pats, DAY, top=20, top_errors=10)
    assert g[0]["mutations"] == ["23139A", "23148C"] and not g[0]["likely_error"]
    assert sum(x["likely_error"] for x in g) == 10