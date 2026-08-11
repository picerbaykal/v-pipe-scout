"""Integration test for cooc queries against live LAPIS.

Marked skip_in_ci because it hits the real WASAP endpoint.
"""

import asyncio
from datetime import datetime

import pytest

from api.wiseloculus import WiseLoculusLapis
from utils.config import get_wiseloculus_url


@pytest.mark.skip_in_ci
def test_get_cooccurrence_returns_data():
    """
    Query cooc at positions 8350 & 8380 for Lugano over a short window.

    These positions are 30bp apart on amplicon 28 — well within any read
    length — so we expect real read-level combinations, not N-dominated rows.
    """
    async def _run():
        client = WiseLoculusLapis(get_wiseloculus_url())
        return await client.get_cooccurrence(
            locationName="Lugano (TI)",
            date_range=(datetime(2026, 3, 9), datetime(2026, 4, 22)),
            positions=[8350, 8380],
        )

    df = asyncio.run(_run())
    assert len(df) > 0, "Expected at least one row"
    assert "date" in df.columns
    assert "count" in df.columns
    assert "[8350]" in df.columns
    assert "[8380]" in df.columns
    assert df["date"].nunique() > 0
    # We should have some non-N combinations (since positions are 30bp apart)
    non_n = df[(df["[8350]"] != "N") & (df["[8380]"] != "N")]
    assert len(non_n) > 0, "Expected some rows with both positions called"
    print(f"\nFetched {len(df)} rows across {df['date'].nunique()} dates")
    print(f"Non-N combinations: {len(non_n)}")
    print(df.head(10))


@pytest.mark.skip_in_ci
def test_get_cooccurrence_empty_range_returns_empty():
    """A date range with no sampling dates returns an empty DataFrame."""
    async def _run():
        client = WiseLoculusLapis(get_wiseloculus_url())
        return await client.get_cooccurrence(
            locationName="Lugano (TI)",
            date_range=(datetime(2020, 1, 1), datetime(2020, 1, 2)),
            positions=[8350, 8380],
        )

    df = asyncio.run(_run())
    assert len(df) == 0