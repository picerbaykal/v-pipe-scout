"""The own-marker rule (process.cooc.check_verdicts, check.own_marker,
2026-10-23): "seen" also needs one of the mutations the lineage gained itself,
not its parent's, when one of them is measured not seen. SB.1.1, 5 Sep – 5 Oct,
Zürich: 2 seen : 1 not seen, but both seen markers are SB.1's (21414C) and
XFG.6.2.2's (18585T); its own 1501A is at 0.8 % of 2,587 reads."""
import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, "/app_shared")

from process.cooc import check_verdicts, own_markers    # noqa: E402
from process.national import national                  # noqa: E402
from process.variant_explorer import check_in_data     # noqa: E402

MK = ["1501A", "21414C", "18585T"]
OWN = {"parent": "SB.1", "markers": ["1501A"]}


def _cells(own_hit, own_cov=2600):
    """one sample: 21414C and 18585T at 10 %, 1501A at own_hit of own_cov"""
    return {"2026-09-30": {"1501A": [own_cov, own_hit, own_hit, own_hit],
                           "21414C": [2000, 200, 200, 200],
                           "18585T": [2000, 200, 200, 200]}}


def test_own_mutation_not_seen_makes_it_unclear():
    off = check_verdicts(MK, _cells(20), dict(cap_reads=0))
    on = check_verdicts(MK, _cells(20), dict(cap_reads=0), own=OWN)
    assert off["verdict"] == "present" and not off["own_missing"]
    assert on["verdict"] == "mixed" and on["own_missing"]
    assert on["own"] == ["1501A"] and on["parent"] == "SB.1"


def test_own_mutation_seen_or_unmeasured_keeps_seen():
    assert check_verdicts(MK, _cells(200), dict(cap_reads=0), own=OWN)["verdict"] == "present"
    # 50 reads: too few to say; the vote stays as it is
    assert check_verdicts(MK, _cells(1, 50), dict(cap_reads=0), own=OWN)["verdict"] == "present"
    # 1–5 %: no vote either way
    assert check_verdicts(MK, _cells(60), dict(cap_reads=0), own=OWN)["verdict"] == "present"


def test_setting_off_and_no_own_markers():
    assert check_verdicts(MK, _cells(20), dict(cap_reads=0, own_marker=False),
                          own=OWN)["verdict"] == "present"
    # no own marker among its markers (PQ.16.1.1): the rule can't apply
    none = {"parent": "PQ.16.1", "markers": []}
    assert check_verdicts(MK, _cells(20), dict(cap_reads=0), own=none)["verdict"] == "present"


def test_own_markers_from_the_tree():
    sigs = {"SB.1": {"21414C", "18585T"}, "SB.1.1": {"21414C", "18585T", "1501A", "22865T"},
            "XFG": {"4214C"}}
    par = {"SB.1": "XFG.6.2.2", "SB.1.1": "SB.1", "XFG": ""}
    assert own_markers("SB.1.1", MK, sigs, par) == OWN
    assert own_markers("SB.1.1", None, sigs, par)["markers"] == ["1501A", "22865T"]
    assert own_markers("XFG", ["4214C"], sigs, par) is None       # recombinant root


def test_city_and_all_cities_use_it():
    q = {"2026-09-30": {"median": 5000, "share_100": 0.9}}
    r = check_in_data(MK, _cells(20), ["2026-09-30"], quality=q, own=OWN)
    assert r["state"] == "mixed" and r["own_missing"] and r["parent"] == "SB.1"
    assert r["timeline"] == [["2026-09-30", "mixed"]]
    assert national({"ZH": r})["answer"]["verdict"] == "mixed"
    assert check_in_data(MK, _cells(20), ["2026-09-30"], quality=q)["state"] == "present"


def _vt():
    path = Path(__file__).resolve().parents[1] / "components" / "variants_table.py"
    spec = importlib.util.spec_from_file_location("vt_own", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_row_says_why():
    vt = _vt()
    q = {"2026-09-30": {"median": 5000, "share_100": 0.9}}
    r = check_in_data(MK, _cells(20), ["2026-09-30"], quality=q, own=OWN)
    r.pop("pairs_checked", None)        # a result from before pair counts
    f = {"status": "confirmed", "addable": False,
         "per_city": {"Zürich (ZH)": {"days": 3, "check": r}}}
    kind, why = vt._finding_status("SB.1.1", f, ["Zürich (ZH)"])
    # no coocs of its own: folded under the tree, the reason in cell words (2026-10-27)
    assert kind == "folded" and why == "suggested 3× in ZH · only SB.1's mutations"
    tip = vt._cell_tip("SB.1.1", "Zürich (ZH)", r)
    assert "the mutation new in SB.1.1 (1501A) not seen" in tip
