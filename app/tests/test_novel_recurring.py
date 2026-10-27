"""Unassigned 'recurring' needs well-covered evidence days in >= 2 weeks (2026-10-27)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from components.variants_table import _novel_weight

M = ["C100T", "G200A"]


def _info(days):
    return {tuple(M): {"ZH": {"timeline": [(d, "day") for d in days] + [("2026-09-02", "seen")]}}}


def test_two_days_same_week_not_recurring():
    w = {"ZH": {"2026-09-01": True, "2026-09-03": True}}
    assert _novel_weight(_info(["2026-09-01", "2026-09-03"]), w, M) == (1, 1)


def test_two_weeks_recurring():
    w = {"ZH": {"2026-09-01": True, "2026-09-10": True}}
    assert _novel_weight(_info(["2026-09-01", "2026-09-10"]), w, M) == (2, 1)


def test_thin_days_not_counted():
    w = {"ZH": {"2026-09-01": True, "2026-09-10": False}}
    assert _novel_weight(_info(["2026-09-01", "2026-09-10"]), w, M) == (1, 1)


def test_seen_marks_not_counted():
    w = {"ZH": {"2026-09-02": True}}
    assert _novel_weight(_info([]), w, M) == (0, 0)
