"""Which seen pairs count as one stretch for a sample.

With the amplicon inserts (2026-11-10) a stretch IS an amplicon: a read comes
off one, so every pair inside one could share a read and no pair across two
ever can. TA.1's 23,148-23,604 was one stretch under the distance rule and is
amplicons 75 and 76 (ArticV542inserts.bed: 75 = 23110-23452, 76 = 23259-23609),
which W-ASAP reports separately at 5 and 6 mutations.

Without a BED, the fallback: pairs with a gap of at most check.stretch_span
(400) bases count as one stretch (2026-10-27): reads in LAPIS span up to ~390
bases (diag_span.py, Zürich Sep 2026: median 330, max 370 on a 10-base grid),
so two pairs 350 bases apart can be one read, one jackpot."""
import sys

sys.path.insert(0, "/app_shared")

import process.amplicons as amplicons      # noqa: E402
from process.cooc import pair_evidence     # noqa: E402

WELL = {"median": 5000, "share_100": 0.9}
DAY = "2026-09-15"

# the three inserts around TA.1's merged stretch, 1-based inclusive
INSERTS = [
    ("SARS-CoV-2_INSERT_74", 22775, 23119),
    ("SARS-CoV-2_INSERT_75", 23110, 23452),
    ("SARS-CoV-2_INSERT_76", 23259, 23609),
]


def _sample(*pairs):
    return {DAY: {p: [1000, 300, 300, 300] for p in pairs}}


def _run(monkeypatch, inserts, *pairs):
    # pair_evidence imports get_inserts at call time, so patching the module
    # attribute decides which geometry it uses
    monkeypatch.setattr(amplicons, "get_inserts", lambda *a, **k: list(inserts))
    return pair_evidence(_sample(*pairs), {DAY: WELL})


# ── with the amplicon inserts ─────────────────────────────────────────────

def test_ta1_pairs_on_amplicons_75_and_76_are_two_stretches(monkeypatch):
    # 22 bases apart, so the distance rule chained them into one; they sit on
    # different amplicons, so they are independent evidence
    r = _run(monkeypatch, INSERTS, "23148C+23401T", "23423T+23604G")
    assert r["positive"] == [DAY]
    assert r["n_stretches"] == 2


def test_two_pairs_on_one_amplicon_are_one_stretch(monkeypatch):
    # NB.1.8.1's pair and another inside amplicon 74's insert
    r = _run(monkeypatch, INSERTS, "22930A+23039G", "22786C+22813T")
    assert r["positive"] == [] and r["n_stretches"] == 0


def test_a_pair_in_an_overlap_goes_to_the_lower_amplicon(monkeypatch):
    # 23401+23423 is inside both 75 and 76: it joins 75, deterministically, so
    # with another 75 pair it is still one stretch
    r = _run(monkeypatch, INSERTS, "23148C+23157G", "23401T+23423T")
    assert r["positive"] == []


def test_pairs_on_three_amplicons_are_three_stretches(monkeypatch):
    r = _run(monkeypatch, INSERTS, "22930A+23039G", "23148C+23401T", "23423T+23604G")
    assert r["positive"] == [DAY]
    assert r["n_stretches"] == 3


# ── fallback without a BED: unchanged behaviour ───────────────────────────

def test_without_bed_two_pairs_350_bases_apart_are_one_stretch(monkeypatch):
    r = _run(monkeypatch, [], "1000A+1010C", "1360G+1370T")
    assert r["positive"] == [] and r["n_stretches"] == 0


def test_without_bed_two_pairs_800_bases_apart_are_two_stretches(monkeypatch):
    r = _run(monkeypatch, [], "1000A+1010C", "1810G+1820T")
    assert r["positive"] == [DAY]


def test_without_bed_ta1_pairs_are_chained_into_one(monkeypatch):
    # documents the old behaviour the amplicon rule corrects
    r = _run(monkeypatch, [], "23148C+23401T", "23423T+23604G")
    assert r["positive"] == []
