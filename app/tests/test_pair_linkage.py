"""Pairs take the linkage test markers have had since 2026-09 (2026-11-10).

NB.1.8.1, Zürich 11 and 15 Oct: its pair 22930A+23039G was on 13,870 reads,
and every one of them showed C at 22896 and 22995, where NB.1.8.1 has A and
T. Those reads belong to another lineage; until now a pair's link_n and
link_ok were set to its hit count, so it passed the test by construction and
the reads were counted as NB.1.8.1's."""
import sys

sys.path.insert(0, "/app_shared")

import process.amplicons as amplicons                          # noqa: E402
from process.cooc import accumulate_check_stats, pair_evidence  # noqa: E402

DAY = "2026-10-11"
WELL = {"median": 5000, "share_100": 0.9}

# a lineage with two pairs far apart, each with a third signature position
# beside it that one read also reaches
SIG = {"L": {"1000A", "1010C", "1020G", "5000T", "5010A", "5020C"}}
PAIRS = {"L": ["1000A+1010C", "5000T+5010A"]}
POS = [1000, 1010, 1020, 5000, 5010, 5020]


def _row(n, **bases):
    r = {"count": n, "date": DAY}
    r.update({f"[{p}]": "N" for p in POS})
    r.update({f"[{k[1:]}]": v for k, v in bases.items()})
    return r


def _verdict(monkeypatch, rows):
    monkeypatch.setattr(amplicons, "get_inserts", lambda *a, **k: [])  # distance fallback
    stats = {}
    accumulate_check_stats(rows, POS, SIG, {}, stats, pairs=PAIRS)
    return stats, pair_evidence(stats["L"], {DAY: WELL})


def test_reads_matching_the_lineage_nearby_count(monkeypatch):
    rows = [_row(500, p1000="A", p1010="C", p1020="G"),
            _row(500, p5000="T", p5010="A", p5020="C")]
    stats, r = _verdict(monkeypatch, rows)
    assert stats["L"][DAY]["1000A+1010C"] == [500, 500, 500, 500]
    assert r["positive"] == [DAY]


def test_reads_contradicting_the_lineage_nearby_do_not(monkeypatch):
    # the NB.1.8.1 case: both mutations of the pair present, the lineage's
    # own neighbour showing another base on every read
    rows = [_row(500, p1000="A", p1010="C", p1020="T"),
            _row(500, p5000="T", p5010="A", p5020="G")]
    stats, r = _verdict(monkeypatch, rows)
    assert stats["L"][DAY]["1000A+1010C"] == [500, 500, 500, 0]
    assert r["positive"] == []


def test_one_contradicted_pair_leaves_one_stretch(monkeypatch):
    rows = [_row(500, p1000="A", p1010="C", p1020="G"),
            _row(500, p5000="T", p5010="A", p5020="G")]
    _stats, r = _verdict(monkeypatch, rows)
    assert r["positive"] == []


def test_a_few_contradicting_reads_are_tolerated(monkeypatch):
    # link_min 0.8: errors and genomes off their consensus don't sink a pair
    rows = [_row(450, p1000="A", p1010="C", p1020="G"),
            _row(50, p1000="A", p1010="C", p1020="T"),
            _row(500, p5000="T", p5010="A", p5020="C")]
    _stats, r = _verdict(monkeypatch, rows)
    assert r["positive"] == [DAY]


def test_no_neighbour_reached_is_judged_on_share_alone(monkeypatch):
    # PJ.2.1-like: the reads cover only the pair -- nothing to test, as a
    # marker with no neighbour in reach
    rows = [_row(500, p1000="A", p1010="C"),
            _row(500, p5000="T", p5010="A")]
    stats, r = _verdict(monkeypatch, rows)
    assert stats["L"][DAY]["1000A+1010C"] == [500, 500, 0, 0]
    assert r["positive"] == [DAY]


def test_results_saved_before_still_read_as_linked(monkeypatch):
    monkeypatch.setattr(amplicons, "get_inserts", lambda *a, **k: [])
    old = {DAY: {"1000A+1010C": [500, 500, 500, 500], "5000T+5010A": [500, 500, 500, 500]}}
    assert pair_evidence(old, {DAY: WELL})["positive"] == [DAY]
