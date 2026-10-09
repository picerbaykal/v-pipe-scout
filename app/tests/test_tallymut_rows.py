"""Deconvolution input from /component/nucleotideMutationsOverTime (2026-10-08).

The rows handed to LolliPop must keep absent mutations (count 0, coverage > 0):
they are the only evidence against a panel variant. Coverage-0 cells (not
sequenced that day) are left out.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from api.wiseloculus import WiseLoculusLapis  # noqa: E402

W = WiseLoculusLapis


def _answer(names, dates, grid):
    return {"data": {"mutations": names,
                     "dateRanges": [{"dateFrom": d, "dateTo": d} for d in dates],
                     "data": grid}}


def test_absent_mutation_is_kept_uncovered_is_dropped():
    ans = _answer(["22896T", "23021G"], ["2026-03-20", "2026-03-24"],
                  [[{"count": 6300, "coverage": 9000}, {"count": 0, "coverage": 0}],
                   [{"count": 0, "coverage": 8500}, {"count": 3, "coverage": 40}]])
    rows = W._component_to_rows(ans, ["22896T", "23021G"])
    got = {(r["pos"], r["date"]): (r["var"], r["cov"]) for r in rows}
    assert got == {("22896T", "2026-03-20"): (6300, 9000),
                   ("23021G", "2026-03-20"): (0, 8500),      # absent: kept
                   ("23021G", "2026-03-24"): (3, 40)}        # low coverage: kept (step 2)
    assert ("22896T", "2026-03-24") not in got                # coverage 0: dropped


def test_names_with_reference_base_map_back():
    ans = _answer(["G22896T"], ["2026-03-20"], [[{"count": 1, "coverage": 10}]])
    rows = W._component_to_rows(ans, ["22896T"])
    assert rows == [{"date": "2026-03-20", "pos": "22896T", "cov": 10, "var": 1}]


def test_unrequested_mutation_is_ignored():
    ans = _answer(["241T"], ["2026-03-20"], [[{"count": 5, "coverage": 10}]])
    assert W._component_to_rows(ans, ["22896T"]) == []


def test_bare_mutation():
    assert W._bare_mutation("C241T") == "241T"
    assert W._bare_mutation("241T") == "241T"
    assert W._bare_mutation(" 241- ") == "241-"


class _Loader:
    _reconstructed_signatures = {"BA.3.2": True}

    def get_signature(self, v):
        return {"XFG": ["22896T", "8350C", "11288-"], "NB.1.8.1": ["22896T", "823T"],
                "BA.3.2": ["999A"]}[v]


def test_panel_mutations_union_sorted_no_deletions():
    assert W._panel_mutations(["XFG", "NB.1.8.1"], _Loader()) == ["823T", "8350C", "22896T"]


def test_reconstructed_node_uses_cowwid_signature():
    muts = W._panel_mutations(["BA.3.2"], _Loader(), cowwid_variants={"BA.3.2": {"1234G"}})
    assert muts == ["1234G"]


# ── 2026-10-09: smaller and safer LAPIS requests ──────────────────────────

def test_mutations_every_variant_carries_are_not_fetched():
    loader = _Loader()
    loader._reconstructed_signatures = {}
    got = W._panel_mutations(["XFG", "NB.1.8.1"], loader, informative_only=True)
    assert got == ["823T", "8350C"]                 # 22896T is in both: LolliPop drops it
    assert "22896T" in W._panel_mutations(["XFG", "NB.1.8.1"], loader)


def test_gateway_timeout_is_retried_then_succeeds(monkeypatch):
    import asyncio
    from api.exceptions import APIError
    w = W.__new__(W)
    w._TALLYMUT_RETRY_WAITS = (0, 0)
    calls = []

    async def fake(session, loc, muts, dates):
        calls.append(1)
        if len(calls) < 3:
            raise APIError(f"nucleotideMutationsOverTime failed: status 504", status_code=504)
        return {"ok": True}
    w._fetch_tallymut_block = fake
    out = asyncio.run(w._fetch_tallymut_block_retry(None, "Lugano (TI)", ["1A"], ["2026-01-01"],
                                                    asyncio.Semaphore(2)))
    assert out == {"ok": True} and len(calls) == 3


def test_other_errors_are_not_retried():
    import asyncio
    import pytest
    from api.exceptions import APIError
    w = W.__new__(W)
    w._TALLYMUT_RETRY_WAITS = (0, 0)
    calls = []

    async def fake(session, loc, muts, dates):
        calls.append(1)
        raise APIError("bad request: status 400", status_code=400)
    w._fetch_tallymut_block = fake
    with pytest.raises(APIError):
        asyncio.run(w._fetch_tallymut_block_retry(None, "Lugano (TI)", ["1A"], ["2026-01-01"],
                                                  asyncio.Semaphore(2)))
    assert len(calls) == 1
