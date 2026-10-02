"""Error hotspots out of novel patterns (scanner._drop_hotspots, 2026-10-02).

A position where the city's data shows >= 2 different new bases is an error
hotspot. A mutation there that no pango lineage carries is left out of novel
patterns; a lineage's own mutation is always kept."""
from process.scanner import _drop_hotspots, _multi_alt_positions, _novel_groups

# 11851T and 22896G: lineage mutations; 11878T, 11933A: noise; 30000A/30010C: new
INDEX = {"11851T": {"RF.6"}, "22896G": {"XFG.7"}, "1A": {"L"}, "2C": {"L"}}
DAY = {"2026-08-12": 10_000, "2026-08-20": 10_000}


def _cand(muts):
    acc = None
    for m in muts:
        acc = set(INDEX.get(m, set())) if acc is None else acc & INDEX.get(m, set())
    return acc or set()


def _p(muts, n=100, d="2026-08-12"):
    return {"mutations": muts, "count": n, "date": d}


def test_real_mutation_plus_noise_is_not_novel():
    kept, hot = _drop_hotspots([_p(["11851T", "11878T"])], [11878], DAY, INDEX, _cand)
    assert kept == [] and hot["patterns"] == 1 and hot["reads"] == 100


def test_lineage_mutation_at_a_hotspot_is_kept():
    pats = [_p(["22896G", "30000A", "30010C"])]
    kept, hot = _drop_hotspots(pats, [22896], DAY, INDEX, _cand)
    assert kept[0]["mutations"] == ["22896G", "30000A", "30010C"] and hot["patterns"] == 0


def test_noise_removed_rest_still_novel():
    kept, _ = _drop_hotspots([_p(["11933A", "30000A", "30010C"])], [11933], DAY, INDEX, _cand)
    assert kept[0]["mutations"] == ["30000A", "30010C"]


def test_rest_carried_by_a_lineage_is_not_novel():
    kept, hot = _drop_hotspots([_p(["1A", "2C", "11933A"])], [11933], DAY, INDEX, _cand)
    assert kept == [] and hot["patterns"] == 1


def test_clean_positions_untouched():
    pats = [_p(["30000A", "30010C"])]
    assert _drop_hotspots(pats, [11878], DAY, INDEX, _cand)[0] == pats


def test_without_data_hotspots_come_from_the_reads():
    # 17857 shows A, C and G on days passing the day rule -> hotspot
    pats = [_p([f"17857{b}", "30000A"], 200) for b in "ACG"]
    assert set(_multi_alt_positions(pats, DAY)) == {17857}
    kept, hot = _drop_hotspots(pats, None, DAY, INDEX, _cand)
    assert kept == [] and hot["patterns"] == 3


def test_groups_count_days():
    pats = [_p(["30000A", "30010C"], 100, "2026-08-12"), _p(["30000A", "30010C"], 100, "2026-08-20")]
    g = _novel_groups(pats, DAY)
    assert g[0]["days"] == ["2026-08-12", "2026-08-20"] and "likely_error" not in g[0]