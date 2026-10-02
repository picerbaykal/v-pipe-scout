"""One-day findings: what happened after their day (scanner._one_day_after)."""
from process.scanner import _one_day_after

DATES = ["2026-08-12", "2026-08-14", "2026-08-20", "2026-08-22"]


def _finding(day, evidence=None, stars=("18585T",)):
    return {"counted_days": [day],
            "evidence_days": evidence or {day: 5000},
            "member_blocks": [{"discriminating": list(stars),
                               "mut_star": {m: True for m in stars}}]}


def _cov(n, dates=DATES, pos="18585"):
    return {d: {pos: n} for d in dates}


def test_latest_sample():
    a = _one_day_after(_finding("2026-08-22"), DATES, _cov(5000))
    assert a["state"] == "latest" and a["later"] == 0


def test_not_seen_since_when_later_samples_covered_it():
    a = _one_day_after(_finding("2026-08-12"), DATES, _cov(5000))
    assert a["state"] == "not_seen" and a["later"] == 3 and a["covered"] == 3


def test_not_covered_since():
    a = _one_day_after(_finding("2026-08-12"), DATES, _cov(40))     # < 100 reads
    assert a["state"] == "not_covered" and a["covered"] == 0


def test_seen_again_below_the_day_rule():
    f = _finding("2026-08-12", {"2026-08-12": 5000, "2026-08-20": 8})
    a = _one_day_after(f, DATES, _cov(5000))
    assert a["state"] == "seen_again" and a["seen"] == 1


def test_without_coverage_data_falls_back_to_not_seen():
    a = _one_day_after(_finding("2026-08-12"), DATES, None)
    assert a["state"] == "not_seen" and a["covered"] is None