"""
Panel completeness signal integration tests.

Tests that the co-occurrence completeness metric correctly signals when
the circulating population is not explained by the panel.

Empirical findings from Lugano (TI), verified against WASAP LAPIS:

  2026-06-17: Panel (KP.2, NB.1.8.1, KP.3, XFG) explains ~99.9% of reads.
              XFG carries 897A which dominates this date — the panel is adequate.

  2026-07-07: Panel (KP.2, NB.1.8.1, KP.3, XFG, LP.8) explains only ~3% of
              reads. 832T is at 21.4% on this date and is not in any panel
              variant's signature — the metric correctly flags the gap.

These tests run against the live WASAP endpoint and are marked skip_in_ci.
Configuration used: discriminating_positions_only=false, per-amplicon batching.
"""
from datetime import datetime
import pytest


@pytest.mark.skip_in_ci
def test_completeness_high_when_panel_adequate():
    """
    On 2026-06-17, XFG (which carries 897A) dominates Lugano reads.
    A panel including XFG should show high completeness — the metric
    confirms the panel explains what is circulating.

    Empirically verified: ~99.9% completeness, 1.3M informative reads.
    """
    import sys
    sys.path.insert(0, "/app_shared")
    from cooc import run_cooc_panel_completeness

    result = run_cooc_panel_completeness(
        location="Lugano (TI)",
        start_date=datetime(2026, 6, 17),
        end_date=datetime(2026, 6, 17),
        variants=["KP.2", "NB.1.8.1", "KP.3", "XFG"],
    )

    assert result["dates"], "Expected data for 2026-06-17"
    c = result["completeness"][0]
    informative = result["matched_counts"][0] + result["unexplained_counts"][0]

    print(f"\n2026-06-17: completeness={c:.4f}  informative={informative:,}")

    assert informative > 100_000, (
        f"Expected many informative reads on a well-covered date, got {informative:,}. "
        "Check that discriminating_positions_only=false in cooc_config.yaml."
    )
    assert c > 0.95, (
        f"Expected high completeness when panel includes the dominant variant, "
        f"got {c:.4f}. Panel: KP.2, NB.1.8.1, KP.3, XFG."
    )


@pytest.mark.skip_in_ci
def test_completeness_low_when_panel_missing_dominant_variant():
    """
    On 2026-07-07, 832T is at 21.4% of reads in Lugano. No panel variant
    carries 832T as a signature mutation. The completeness metric should
    signal that the panel cannot explain what is circulating.

    Empirically verified: ~3% completeness, 390k informative reads.
    This is the core surveillance value of the metric — it flags stale panels.
    """
    import sys
    sys.path.insert(0, "/app_shared")
    from cooc import run_cooc_panel_completeness

    result = run_cooc_panel_completeness(
        location="Lugano (TI)",
        start_date=datetime(2026, 7, 7),
        end_date=datetime(2026, 7, 7),
        variants=["KP.2", "NB.1.8.1", "KP.3", "XFG", "LP.8"],
    )

    assert result["dates"], "Expected data for 2026-07-07"
    c = result["completeness"][0]
    informative = result["matched_counts"][0] + result["unexplained_counts"][0]

    print(f"\n2026-07-07: completeness={c:.4f}  informative={informative:,}")

    assert informative > 100_000, (
        f"Expected many informative reads on a well-covered date, got {informative:,}. "
        "Check that discriminating_positions_only=false in cooc_config.yaml."
    )
    assert c < 0.1, (
        f"Expected low completeness when dominant variant (832T carrier) is "
        f"missing from panel, got {c:.4f}. The metric should flag this gap."
    )


@pytest.mark.skip_in_ci
def test_completeness_varies_between_dates():
    """
    The completeness metric produces meaningfully different values across
    dates, reflecting genuine changes in the circulating population.

    2026-06-17: panel-adequate date (~99.9%)
    2026-07-07: panel-inadequate date (~3%)

    The gap between them should be large — this is not noise, it is the
    metric correctly tracking variant turnover.
    """
    import sys
    sys.path.insert(0, "/app_shared")
    from cooc import run_cooc_panel_completeness

    variants = ["KP.2", "NB.1.8.1", "KP.3", "XFG", "LP.8"]

    r1 = run_cooc_panel_completeness(
        location="Lugano (TI)",
        start_date=datetime(2026, 6, 17),
        end_date=datetime(2026, 6, 17),
        variants=variants,
    )
    r2 = run_cooc_panel_completeness(
        location="Lugano (TI)",
        start_date=datetime(2026, 7, 7),
        end_date=datetime(2026, 7, 7),
        variants=variants,
    )

    c1 = r1["completeness"][0] if r1["dates"] else None
    c2 = r2["completeness"][0] if r2["dates"] else None

    print(f"\n2026-06-17 completeness: {c1:.4f}")
    print(f"2026-07-07 completeness: {c2:.4f}")
    print(f"Gap: {abs(c1 - c2):.4f}")

    assert c1 is not None and c2 is not None
    assert c1 > 0.95, f"2026-06-17 should be panel-adequate, got {c1:.4f}"
    assert c2 < 0.10, f"2026-07-07 should be panel-inadequate, got {c2:.4f}"
    assert c1 - c2 > 0.85, (
        f"Expected large completeness gap between dates "
        f"(2026-06-17={c1:.4f}, 2026-07-07={c2:.4f}). "
        "Same panel, different circulating variants — the metric should reflect this."
    )