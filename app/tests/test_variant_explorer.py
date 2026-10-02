"""Investigate a variant: ★ markers with the panel check's rule, and the
per-city answer from read counts."""
from process.recombinants import override_parents
from process.variant_explorer import check_in_data, investigate_variant


class _Loader:
    # A -> A.1 -> A.1.1 ; B unrelated
    RAW = {"A": {"parent": ""}, "A.1": {"parent": "A"}, "A.1.1": {"parent": "A.1"},
           "B": {"parent": ""}}
    SIG = {"A": {"1A", "2C"}, "A.1": {"1A", "2C", "10G", "11T"},
           "A.1.1": {"1A", "2C", "10G", "11T", "20A"}, "B": {"1A", "30C"}}

    def get_raw_data(self):
        return self.RAW

    def get_signature(self, v):
        return self.SIG.get(v, set())


def test_markers_and_reason():
    with override_parents({}):
        r = investigate_variant("A.1", _Loader(), panel=["B"])
    assert r["found"] and r["detectable"]
    # 1A: carried by panel variant B -> not a marker; 2C: one outsider (A) <= 5
    assert set(r["markers"]) == {"2C", "10G", "11T"}
    assert r["children"] == ["A.1.1"]


def test_unknown_name():
    assert investigate_variant("ZZ.9", _Loader())["found"] is False


def _cells(cov, hit):
    return {"10G": [cov, hit, hit, hit], "11T": [cov, hit, hit, hit]}


def test_present_absent_uncovered_days():
    pd = {"2026-08-01": _cells(500, 300), "2026-08-05": _cells(500, 0),
          "2026-08-09": _cells(20, 10)}
    r = check_in_data(["10G", "11T"], pd, ["2026-08-01", "2026-08-05", "2026-08-09",
                                           "2026-08-12"])
    assert [m for _, m in r["timeline"]] == ["present", "absent", "uncovered", "uncovered"]
    assert r["state"] == "present"                   # pooled 600/1040 reads


def test_old_variant_is_absent_not_silent():
    pd = {d: _cells(800, 0) for d in ("2026-07-01", "2026-07-08")}
    assert check_in_data(["10G", "11T"], pd, list(pd))["state"] == "absent"


def test_no_marker():
    assert check_in_data([], {}, ["2026-08-01"])["state"] == "no_marker"