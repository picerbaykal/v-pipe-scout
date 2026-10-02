"""Which recombinants inherited a mutation from a given variant (2026-10-02).

The ★ marker rule and the scanner tolerate carriers of a mutation that sit
under a recombinant MADE FROM the variant: XFV = LP.8.1 × XFG.3.3.1 carries
XFG's 8350C because it got that piece of genome from XFG, so seeing 8350C still
means "XFG material". Before, every carrier under any other recombinant was
tolerated (up to 60) — including unrelated, circulating ones: SV.4 (an NB.1.8.1
descendant under XDV) carries B.1.1.7's 3267T and showed up at 52 % in Zürich,
making an extinct variant look present.

Rule: recombinant R inherited mutation m from variant V when one of R's parents
  - is V or a descendant of V, and carries m, or
  - is itself under a recombinant that inherited m from V (followed up to 6
    levels: XGA = XFJ.4.1 × PY.1 × XFG.6.2).
Parents come from pango-designation's alias_key.json (app/data/alias_key.json);
names are compared un-aliased (NB.1.8.1 = XDV.1.5.1.1.8.1).

Without alias_key.json every recombinant counts as inherited — the old
behaviour — and a warning is logged.
"""
import json
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Set

logger = logging.getLogger(__name__)

_REC_ROOT_RE = re.compile(r"^X[A-Z]+$")
_ALIAS_PATH = Path(__file__).resolve().parents[1] / "data" / "alias_key.json"

# Tests use small made-up trees whose recombinant names (XG, XY …) also exist in
# the real alias key; they set their own parents here (see override_parents).
ALIAS_OVERRIDE: Optional[dict] = None


class Recombinants:
    def __init__(self, parent_map: Dict[str, str],
                 lineage_signatures: Dict[str, Set[str]],
                 alias_path: Optional[Path] = None):
        self.parent = parent_map
        self.sigs = lineage_signatures
        self.kids: Dict[str, List[str]] = {}
        for c, p in parent_map.items():
            if p:
                self.kids.setdefault(p, []).append(c)
        raw = {}
        try:
            raw = (dict(ALIAS_OVERRIDE) if ALIAS_OVERRIDE is not None
                   else json.loads(Path(alias_path or _ALIAS_PATH).read_text()))
        except Exception as e:
            logger.warning(f"[recombinants] alias_key.json not readable ({e}): "
                           "every recombinant carrier is treated as inherited")
        self.available = bool(raw)
        self.parents_of = {k: list(dict.fromkeys(x.rstrip("*") for x in v if x))
                           for k, v in raw.items() if k.startswith("X") and isinstance(v, list)}
        self.alias = {k: v for k, v in raw.items()
                      if isinstance(v, str) and v and not k.startswith("X")}
        self._full: Dict[str, str] = {}
        self._canon = {}
        for lin in set(parent_map) | set(lineage_signatures):
            self._canon.setdefault(self.full(lin), lin)
        self._fam: Dict[str, Set[str]] = {}
        self._rroot: Dict[str, Optional[str]] = {}
        self._memo: Dict[tuple, bool] = {}

    # ── names ──────────────────────────────────────────────────────────────
    def full(self, name: str) -> str:
        """'NB.1.8.1' -> 'XDV.1.5.1.1.8.1'."""
        if name in self._full:
            return self._full[name]
        x = name.rstrip("*")
        for _ in range(10):
            head, _, rest = x.partition(".")
            if head in self.alias:
                x = self.alias[head] + ("." + rest if rest else "")
            else:
                break
        self._full[name] = x
        return x

    def canon(self, name: str) -> str:
        """A name as alias_key writes it -> the lineage name in the pango file."""
        name = name.rstrip("*")
        if name in self.parent or name in self.sigs:
            return name
        return self._canon.get(self.full(name), name)

    # ── tree ───────────────────────────────────────────────────────────────
    def family(self, v: str) -> Set[str]:
        if v not in self._fam:
            fam, stack = {v}, [v]
            while stack:
                for c in self.kids.get(stack.pop(), []):
                    if c not in fam:
                        fam.add(c)
                        stack.append(c)
            self._fam[v] = fam
        return self._fam[v]

    def rec_root(self, lin: str) -> Optional[str]:
        """Topmost X.. ancestor-or-self (recombinant roots have an empty parent)."""
        if lin not in self._rroot:
            x, found, seen = lin, None, set()
            while x and x not in seen:
                seen.add(x)
                if _REC_ROOT_RE.match(x):
                    found = x
                x = self.parent.get(x, "")
            self._rroot[lin] = found
        return self._rroot[lin]

    def in_family(self, lin: str, v: str) -> bool:
        if lin in self.family(v):
            return True
        lf, vf = self.full(lin), self.full(v)
        return lf == vf or lf.startswith(vf + ".")

    # ── the rule ───────────────────────────────────────────────────────────
    def made_from(self, root: str, v: str, m: Optional[str] = None) -> bool:
        """Did recombinant `root` get mutation m from variant v? With m None:
        is v (an ancestor of) one of its parents (read-level use, where a read
        carries several mutations)."""
        if not self.available:
            return True
        key = (root, v, m)
        if key not in self._memo:
            self._memo[key] = self._made_from(root, v, m, 0, set())
        return self._memo[key]

    def _made_from(self, root, v, m, depth, seen) -> bool:
        if depth > 6 or root in seen:
            return False
        seen.add(root)
        for p in self.parents_of.get(root, []):
            pc = self.canon(p)
            if m is not None and pc in self.sigs and m not in self.sigs[pc]:
                continue                      # this parent didn't pass m on
            if self.in_family(pc, v):
                return True
            pr = self.rec_root(pc) if pc in self.parent else (pc if pc in self.parents_of else None)
            if pr and pr != root and self._made_from(pr, v, m, depth + 1, seen):
                return True
        return False

    def tolerated(self, lin: str, v: str, m: Optional[str] = None,
                  own_roots: Optional[Set[Optional[str]]] = None) -> bool:
        """Carrier `lin` (outside v's family) doesn't count against v: it sits
        under a recombinant (other than v's own) that inherited m from v."""
        r = self.rec_root(lin)
        if not r or r in (own_roots or set()):
            return False
        return self.made_from(r, v, m)


_CACHE: list = []


class override_parents:
    """Context manager for tests: `with override_parents({"XY1": ["XG.1", "J"]}):`
    uses these recombinant parents instead of alias_key.json."""

    def __init__(self, parents: dict):
        self.parents = parents

    def __enter__(self):
        global ALIAS_OVERRIDE
        self._old, ALIAS_OVERRIDE = ALIAS_OVERRIDE, dict(self.parents)
        _CACHE.clear()
        return self

    def __exit__(self, *exc):
        global ALIAS_OVERRIDE
        ALIAS_OVERRIDE = self._old
        _CACHE.clear()


def get_recombinants(parent_map: Dict[str, str],
                     lineage_signatures: Dict[str, Set[str]]) -> Recombinants:
    """One instance per (parent map, signatures) object — the worker keeps both
    for its lifetime."""
    if _CACHE and _CACHE[0][0] is parent_map and _CACHE[0][1] is lineage_signatures:
        return _CACHE[0][2]
    rc = Recombinants(parent_map, lineage_signatures)
    _CACHE[:] = [(parent_map, lineage_signatures, rc)]
    return rc