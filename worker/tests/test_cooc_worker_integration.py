"""Integration test for the worker-side cooc pipeline.

Marked skip_in_ci because it hits the real WASAP endpoint and takes
several seconds to complete.
"""

from datetime import datetime

import pytest

import sys
sys.path.insert(0, "/app_shared")  # if running from worker perspective


@pytest.mark.skip_in_ci
def test_run_cooc_panel_completeness_end_to_end():
    """
    Compute panel completeness for Lugano over a 2-week window.

    Verifies the full pipeline: amp_dict building, amplicon batching,
    LAPIS queries, annotation, aggregation.
    """
    # Import here so pytest collection doesn't fail on non-worker envs
    from cooc import run_cooc_panel_completeness

    result = run_cooc_panel_completeness(
        location="Lugano (TI)",
        start_date=datetime(2026, 3, 9),
        end_date=datetime(2026, 4, 22),
        variants=["KP.2", "NB.1.8.1", "XFG"],
    )

    assert result["location"] == "Lugano (TI)"
    assert isinstance(result["dates"], list)
    assert isinstance(result["matched_counts"], list)
    assert isinstance(result["unexplained_counts"], list)
    assert isinstance(result["completeness"], list)

    # All lists aligned
    n = len(result["dates"])
    assert len(result["matched_counts"]) == n
    assert len(result["unexplained_counts"]) == n
    assert len(result["completeness"]) == n

    # Expect at least a couple of dates in that window
    assert n >= 2

    # Completeness is a valid fraction
    for c in result["completeness"]:
        if c == c:  # not NaN
            assert 0.0 <= c <= 1.0

    print(f"\n=== Panel completeness for Lugano ===")
    print(f"Dates: {result['dates']}")
    print(f"Matched counts: {result['matched_counts']}")
    print(f"Unexplained counts: {result['unexplained_counts']}")
    print(f"Completeness: {[round(c, 3) for c in result['completeness']]}")