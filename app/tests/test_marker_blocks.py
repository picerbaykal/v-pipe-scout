"""Investigate a variant → Signal over time: groups of a lineage's mutations a
read can cover, kept when they hold a ★ marker (process.variant_explorer)."""
from process.recombinants import override_parents
from process.variant_explorer import marker_blocks


class _L:
    RAW = {"A": {"parent": ""}, "A.1": {"parent": "A"}, "B": {"parent": ""}}
    SIG = {"A": {"100A", "200C"}, "A.1": {"100A", "200C", "150G", "900T", "1500C"},
           "B": {"100A", "300C"}}

    def get_raw_data(self):
        return self.RAW

    def get_signature(self, v):
        return self.SIG.get(v, set())


def test_groups_within_a_read_with_markers():
    with override_parents({}):
        bl = marker_blocks("A.1", _L(), ["B"])
    assert [b["discriminating"] for b in bl] == [["100A", "150G", "200C"], ["900T"], ["1500C"]]
    assert bl[0]["mut_star"] == {"100A": False, "150G": True, "200C": True}   # B carries 100A