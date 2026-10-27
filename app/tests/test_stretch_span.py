"""Seen pairs with a gap of at most check.stretch_span (400) bases count as one
stretch (2026-10-27): reads in LAPIS span up to ~390 bases (diag_span.py,
Zürich Sep 2026: median 330, max 370 on a 10-base grid), so two pairs 350
bases apart can be one read, one jackpot."""
import sys

sys.path.insert(0, "/app_shared")

from process.cooc import pair_evidence     # noqa: E402

WELL = {"median": 5000, "share_100": 0.9}


def _sample(*pairs):
    return {"2026-09-15": {p: [1000, 300, 300, 300] for p in pairs}}


def test_two_pairs_350_bases_apart_are_one_stretch():
    r = pair_evidence(_sample("1000A+1010C", "1360G+1370T"), {"2026-09-15": WELL})
    assert r["positive"] == [] and r["n_stretches"] == 0


def test_two_pairs_800_bases_apart_are_two_stretches():
    r = pair_evidence(_sample("1000A+1010C", "1810G+1820T"), {"2026-09-15": WELL})
    assert r["positive"] == ["2026-09-15"]
