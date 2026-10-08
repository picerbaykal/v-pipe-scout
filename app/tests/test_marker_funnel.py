"""Investigate a variant: the ★ marker steps (marker_funnel) give the same
markers as the check, and the per-marker numbers flag near misses."""
from process.recombinants import override_parents
from process.variant_explorer import check_in_data, marker_funnel, variant_markers


class _L:
    RAW = {"A": {"parent": ""}, "A.1": {"parent": "A"}, "B": {"parent": ""},
           **{f"C{i}": {"parent": ""} for i in range(7)}}
    SIG = {"A": {"100A", "200C"}, "A.1": {"100A", "200C", "150G", "900T", "1500C"},
           "B": {"100A", "300C"}, **{f"C{i}": {"900T", f"{5000 + i}G"} for i in range(7)}}

    def get_raw_data(self):
        return self.RAW

    def get_signature(self, v):
        return self.SIG.get(v, set())


def test_funnel_matches_the_check():
    with override_parents({}):
        f = marker_funnel("A.1", _L(), ["B"])
        assert sorted(f["markers"]) == sorted(variant_markers("A.1", _L(), ["B"]))
    assert f["by_panel"] == {"100A": ["B"]}                # step 2
    assert "900T" not in f["markers"]                      # 7 outsiders > 5
    assert [r["mutation"] for r in f["near"]] == ["900T"]  # 1–5 over: near miss


def test_detail_flags_a_near_miss():
    pd = {"2026-08-01": {"150G": [500, 22, 22, 22]}}       # 4.4 %: present needs 5 %
    d = check_in_data(["150G"], pd, ["2026-08-01"])["detail"][0]
    assert d["why"].startswith("between") and "present needs" in d["near"]
    