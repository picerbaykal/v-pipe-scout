"""Worker: reads covering each position (cooc.coverage_from_rows)."""
from cooc import coverage_from_rows


def test_counts_reads_with_a_base_not_N_or_deletion():
    rows = [{"date": "2026-08-12", "count": 10, "[100]": "A", "[200]": "N"},
            {"date": "2026-08-12", "count": 5, "[100]": "-", "[200]": "C"},
            {"date": "2026-08-12", "count": 2, "[100]": "T", "[200]": None}]
    assert coverage_from_rows(rows, [100, 200, 300]) == {"100": 12, "200": 5}


def test_empty():
    assert coverage_from_rows([], [100]) == {}