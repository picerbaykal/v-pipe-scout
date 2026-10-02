"""Option C: a carrier under another recombinant doesn't count against a
variant's marker only if that recombinant inherited the mutation from it
(process.recombinants, process.cooc.specific_markers). Uses the real
app/data/alias_key.json with a small hand-made tree."""
from process.recombinants import Recombinants, _CACHE
from process.cooc import specific_markers

PARENT = {
    "B.1.1": "", "B.1.1.7": "B.1.1",
    "LF.7": "", "LP.8.1": "", "LP.8.1.2": "LP.8.1",
    "XFG": "", "XFG.3": "XFG", "XFG.3.3": "XFG.3", "XFG.3.3.1": "XFG.3.3", "XFG.3.4.1": "XFG.3",
    "XFV": "", "XFV.1": "XFV",
    "XDV": "", "NB.1.8.1": "XDV", "SV.4": "NB.1.8.1",
    "XHD": "",
}
SIGS = {
    "B.1.1.7": {"3267T", "28048T"},
    "XFG": {"8350C", "4214C"}, "XFG.3": {"8350C", "4214C"}, "XFG.3.3": {"8350C", "4214C"},
    "XFG.3.3.1": {"8350C", "4214C"}, "XFG.3.4.1": {"8350C", "4214C"},
    "XFV": {"8350C", "1954A"}, "XFV.1": {"8350C", "1954A"},
    "LP.8.1": {"1954A"}, "LP.8.1.2": {"1954A"},
    "NB.1.8.1": {"29164C", "823T"}, "SV.4": {"29164C", "823T", "3267T"},
    "XHD": {"29164C", "8350C"},
}


def _rc():
    return Recombinants(PARENT, SIGS)


def test_recombinant_made_from_the_variant_is_tolerated():
    rc = _rc()
    assert rc.made_from("XFV", "XFG", "8350C")             # XFV = LP.8.1 × XFG.3.3.1
    assert rc.tolerated("XFV.1", "XFG", "8350C", {"XFG"})


def test_unrelated_circulating_recombinant_counts():
    rc = _rc()                                           # SV.4 sits under XDV
    assert not rc.tolerated("SV.4", "B.1.1.7", "3267T")


def test_parent_must_carry_the_mutation():
    rc = _rc()                                           # XFV got 1954A from LP.8.1, not XFG
    assert not rc.made_from("XFV", "XFG", "1954A")
    assert rc.made_from("XFV", "LP.8.1", "1954A")


def test_unaliased_parent_name():
    rc = _rc()                                           # XHD = XFG.3.4.1 × XDV.1.5.1.1.8.1
    assert rc.made_from("XHD", "NB.1.8.1", "29164C")


def test_specific_markers_option_c():
    _CACHE.clear()
    vs = {"B.1.1.7": SIGS["B.1.1.7"], "XFG": SIGS["XFG"]}
    out = specific_markers(vs, SIGS, PARENT, cfg={"out_other_max": 0})
    assert "8350C" in out["XFG"]        # XFV and XHD carry it, both made from XFG
    assert "4214C" in out["XFG"]
    assert "3267T" not in out["B.1.1.7"]                # SV.4 counts now
    assert "28048T" in out["B.1.1.7"]


def test_unrelated_recombinant_drops_a_marker():
    _CACHE.clear()
    parent = dict(PARENT, XEC="")                       # XEC = KS.1.1 × KP.3.3
    sigs = dict(SIGS, XEC={"8350C"})
    out = specific_markers({"XFG": sigs["XFG"]}, sigs, parent, cfg={"out_other_max": 0})
    assert "8350C" not in out["XFG"] and "4214C" in out["XFG"]