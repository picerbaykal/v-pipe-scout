from __future__ import annotations

import json
import re
import logging
import shutil
import urllib.request
from pathlib import Path
from typing import cast
import logging

PANGO_DATA_DIR = Path(__file__).parent.parent / "data"
PANGO_SUMMARY_DEFAULT = PANGO_DATA_DIR / "pango_summary.json"
PANGO_SUMMARY_CACHE = Path("/app/.cache/pango/pango_summary.json")  # Docker runtime path

def _cache_is_valid(path: Path) -> bool:
    """A cache is valid only if it has no ORPHANED lineages — lineages with an
    empty parent whose naming-parent exists in the file. Orphans (introduced by
    an old merge step) break clade traversal and corrupt scanner resolution, so
    a cache containing them must be rejected in favour of the clean checked-in
    file. Cheap check: scan parents once."""
    try:
        with path.open("r", encoding="utf-8") as f:
            d = json.load(f)
        for lin, entry in d.items():
            if entry.get("parent"):
                continue
            if "." in lin and ".".join(lin.split(".")[:-1]) in d:
                return False  # orphaned lineage → corrupted cache
        return True
    except Exception:
        return False


def get_pango_summary_path() -> Path:
    """
    Return the path to the best available pango_summary.json.

    Prefers the runtime cache ONLY if it passes validation (no orphaned
    lineages). A corrupted cache (orphans from an old merge) is rejected and the
    clean checked-in default is used instead. This prevents a stale/corrupted
    cache — which can survive volume/container resets — from breaking clade-level
    scanner detection.
    """
    if PANGO_SUMMARY_CACHE.exists() and _cache_is_valid(PANGO_SUMMARY_CACHE):
        return PANGO_SUMMARY_CACHE
    return PANGO_SUMMARY_DEFAULT

class PangoLoader:
    """
       Load and process lineage mutation signatures from pango_summary.json.

       Provides access to:
       - full mutation signatures
       - private mutations
       - designation dates
    """

    path: Path
    raw_data: dict[str, dict[str, object]]
    _signatures: dict[str, set[str]]
    _private_mutations: dict[str, set[str]]
    _designation_dates: dict[str, str | None]

    def __init__(self, pango_summary_path: str | Path) -> None:
        # convert input to Path
        self.path = Path(pango_summary_path)

        # check file exists
        if not self.path.exists():
            raise FileNotFoundError(f"Pango summary not found: {self.path}")

        # load JSON
        with self.path.open("r", encoding="utf-8") as handle:
            self.raw_data = cast(
                dict[str, dict[str, object]],
                json.load(handle),
            )

        # initialize processed containers — signatures built lazily on first access
        self._signatures = {}
        self._private_mutations = {}
        self._designation_dates = {}
        self._reconstructed_signatures: set[str] = set()
        # pre-process only designation dates (lightweight, needed for sorting/display)
        # signatures are built lazily in get_signature() to avoid loading all 6,007
        # lineages into memory when only 2-29 are needed for a given panel
        self._process_dates()


    @staticmethod
    def _normalize_substitution(mutation: str) -> str:
        """
        Convert a pango_summary substitution to scan-compatible format.

        Parameters
        ----------
        mutation : str
            Mutation in pango_summary format, for example "C241T".

        Returns
        -------
        str
            Mutation in scan-compatible format, for example "241T".
        """
        return mutation[1:]


    def _process_dates(self) -> None:
        """Lightweight pass — only extract designation dates (needed for display/sorting)."""
        for lineage, entry in self.raw_data.items():
            designation_date = entry.get("designationDate")
            self._designation_dates[lineage] = (
                designation_date if isinstance(designation_date, str) else None
            )

    def _process_lineage(self, lineage: str) -> None:
        """Process one lineage on demand — called lazily from get_signature()."""
        if lineage in self._signatures:
            return  # already processed
        entry = self.raw_data.get(lineage, {})
        substitutions = entry.get("nucSubstitutions", [])
        if not isinstance(substitutions, list):
            substitutions = []
        private_substitutions = entry.get("nucSubstitutionsNew", [])
        if not isinstance(private_substitutions, list):
            private_substitutions = []
        signature: set[str] = set()
        for mutation in substitutions:
            if isinstance(mutation, str) and len(mutation) > 1:
                signature.add(self._normalize_substitution(mutation))
        private_mutations: set[str] = set()
        for mutation in private_substitutions:
            if isinstance(mutation, str) and len(mutation) > 1:
                private_mutations.add(self._normalize_substitution(mutation))
        # handle deletions
        nuc_deletions = entry.get("nucDeletions", [])
        if not isinstance(nuc_deletions, list):
            nuc_deletions = []
        for deletion in nuc_deletions:
            if isinstance(deletion, str) and "-" in deletion:
                try:
                    start, end = deletion.split("-")
                    for pos in range(int(start), int(end) + 1):
                        signature.add(f"{pos}-")
                except (ValueError, IndexError):
                    pass
        nuc_deletions_new = entry.get("nucDeletionsNew", [])
        if not isinstance(nuc_deletions_new, list):
            nuc_deletions_new = []
        for deletion in nuc_deletions_new:
            if isinstance(deletion, str) and "-" in deletion:
                try:
                    start, end = deletion.split("-")
                    for pos in range(int(start), int(end) + 1):
                        private_mutations.add(f"{pos}-")
                except (ValueError, IndexError):
                    pass
        self._signatures[lineage] = signature
        self._private_mutations[lineage] = private_mutations

    def _process(self) -> None:
        """
        Process raw lineage data into cleaned internal lookup tables.
        """
        for lineage, entry in self.raw_data.items():
            substitutions = entry.get("nucSubstitutions", [])
            if not isinstance(substitutions, list):
                substitutions = []

            private_substitutions = entry.get("nucSubstitutionsNew", [])
            if not isinstance(private_substitutions, list):
                private_substitutions = []

            designation_date = entry.get("designationDate")

            signature: set[str] = set()
            for mutation in substitutions:
                if isinstance(mutation, str) and len(mutation) > 1:
                    signature.add(self._normalize_substitution(mutation))

            private_mutations: set[str] = set()
            for mutation in private_substitutions:
                if isinstance(mutation, str) and len(mutation) > 1:
                    private_mutations.add(self._normalize_substitution(mutation))

            # Add nucleotide deletions — pango_summary stores as ranges
            # e.g. "11288-11296" → expand to "11288-", "11289-", ..., "11296-"
            # matching the tallymut pos format used by LolliPop
            nuc_deletions = entry.get("nucDeletions", [])
            if not isinstance(nuc_deletions, list):
                nuc_deletions = []

            for deletion in nuc_deletions:
                if isinstance(deletion, str) and '-' in deletion:
                    try:
                        start, end = deletion.split('-')
                        for pos in range(int(start), int(end) + 1):
                            signature.add(f"{pos}-")
                    except (ValueError, IndexError):
                        pass

            nuc_deletions_new = entry.get("nucDeletionsNew", [])
            if not isinstance(nuc_deletions_new, list):
                nuc_deletions_new = []

            for deletion in nuc_deletions_new:
                if isinstance(deletion, str) and '-' in deletion:
                    try:
                        start, end = deletion.split('-')
                        for pos in range(int(start), int(end) + 1):
                            private_mutations.add(f"{pos}-")
                    except (ValueError, IndexError):
                        pass


            self._signatures[lineage] = signature
            self._private_mutations[lineage] = private_mutations
            self._designation_dates[lineage] = (
                designation_date if isinstance(designation_date, str) else None
            )


    def _fill_empty_node_signatures(self) -> None:
        """
        Reconstruct signatures for lineage nodes that have zero mutations
        because no sequences are directly labeled with that name.

        This happens for intermediate nodes like BA.3.2, whose sequences are
        all labeled BA.3.2.1, BA.3.2.2, etc. The pango_summary centroid method
        produces an empty signature for such nodes because it averages over zero
        sequences. The same pattern occurred historically with Omicron (before
        BA.1/BA.2 split) and Delta (B.1.617.1 vs B.1.617.2).

        Fix: for any variant with an empty signature that has children in the
        tree, set its signature to the intersection of its children's signatures.
        This is equivalent to what CovSpectrum computes for "BA.3.2*" — the
        mutations common to all descendants. We process children before parents
        (topological order) so nested empty nodes are handled correctly.
        """
        # Build children map from raw_data's "parent" field
        children: dict[str, list[str]] = {v: [] for v in self.raw_data}
        for lineage, entry in self.raw_data.items():
            parent = entry.get("parent")
            if parent and isinstance(parent, str) and parent.strip():
                parent = parent.strip()
                if parent in children:
                    children[parent].append(lineage)

        # Topological sort: process leaves first so nested empty parents resolve
        # correctly. Use post-order DFS from all roots.
        visited: set[str] = set()
        topo_order: list[str] = []

        def _dfs(node: str) -> None:
            if node in visited:
                return
            visited.add(node)
            for child in children.get(node, []):
                _dfs(child)
            topo_order.append(node)

        for lineage in self.raw_data:
            _dfs(lineage)

        # Now process in topo order (children before parents)
        filled: list[str] = []
        for lineage in topo_order:
            if self._signatures.get(lineage):  # non-empty → nothing to do
                continue
            child_list = children.get(lineage, [])
            if not child_list:
                continue  # genuine leaf with zero mutations — leave as-is

            # Collect the (possibly already-filled) signatures of children
            child_sigs = [
                self._signatures[c]
                for c in child_list
                if self._signatures.get(c)  # skip children that also have zero sig
            ]
            if not child_sigs:
                continue  # all children are also empty — can't reconstruct

            # Intersection = mutations shared by ALL children (common core)
            intersection = child_sigs[0].intersection(*child_sigs[1:])
            if intersection:
                self._signatures[lineage] = intersection
                filled.append(lineage)
                self._reconstructed_signatures.add(lineage)

        if filled:
            logging.getLogger(__name__).info(
                "Reconstructed signatures for %d empty-node lineage(s) from "
                "children intersection: %s",
                len(filled),
                filled,
            )

    def get_signature(self, lineage: str) -> set[str]:
        # lazy: process this lineage if not yet done
        if lineage not in self._signatures:
            if lineage in self.raw_data:
                self._process_lineage(lineage)
            else:
                # lineage absent from file — reconstruct from children via prefix match
                prefix = lineage + "."
                # ensure all candidate children are processed
                for name in self.raw_data:
                    if name.startswith(prefix) and name not in self._signatures:
                        self._process_lineage(name)
                child_sigs = [
                    self._signatures[name]
                    for name in self._signatures
                    if name.startswith(prefix) and self._signatures[name]
                ]
                if not child_sigs:
                    self._signatures[lineage] = set()
                    return set()
                intersection = child_sigs[0].intersection(*child_sigs[1:])
                self._signatures[lineage] = intersection
                if intersection:
                    self._reconstructed_signatures.add(lineage)
                return intersection
        return self._signatures[lineage]

    def get_private_mutations(self, lineage: str) -> set[str]:
        return self._private_mutations.get(lineage, set())

    def is_reconstructed(self, lineage: str) -> bool:
        """
        Return True if this lineage had zero mutations in pango_summary and its
        signature was reconstructed from children's intersection (e.g. BA.3.2).
        """
        return lineage in self._reconstructed_signatures

    def is_known(self, lineage: str) -> bool:
        """Return True if we can compute a non-empty signature for this lineage,
        either directly from pango_summary or reconstructed from sub-lineages."""
        return bool(self.get_signature(lineage))

    def get_raw_data(self) -> dict:
        """Return the raw pango_summary dict for parent graph construction."""
        return self.raw_data

PANGO_SUMMARY_URL = (
    "https://raw.githubusercontent.com/corneliusroemer/pango-sequences"
    "/refs/heads/main/data/pango-consensus-sequences_summary.json"
)
FREYJA_BARCODES_URL = (
    "https://raw.githubusercontent.com/andersen-lab/Freyja"
    "/main/freyja/data/usher_barcodes.feather"
)

def _get_pango_source() -> str:
    """Read pango.source from config.yaml, default to cornelius."""
    try:
        import yaml
        config_path = Path(__file__).parent.parent / "config.yaml"
        with open(config_path) as f:
            cfg = yaml.safe_load(f)
        return cfg.get("pango", {}).get("source", "cornelius")
    except Exception:
        return "cornelius"


NEXTCLADE_TREE_URL = (
    "https://raw.githubusercontent.com/nextstrain/nextclade_data"
    "/refs/heads/master/data/nextstrain/sars-cov-2/wuhan-hu-1/orfs/tree.json"
)

_TREE_MUT = re.compile(r"^([ACGTN-])(\d+)([ACGTN-])$")
_RECOMB = re.compile(r"^X[A-Z]+$")


def _nuc_ranges(positions) -> list:
    out: list[list[int]] = []
    for pos in sorted(positions):
        if out and pos == out[-1][1] + 1:
            out[-1][1] = pos
        else:
            out.append([pos, pos])
    return [f"{a}-{b}" for a, b in out]


def build_summary_from_tree(tree: dict, existing_dates: dict | None = None) -> dict:
    """Convert a Nextclade Auspice v2 tree.json into a pango_summary-compatible
    dict. Each node carries only BRANCH mutations, so a lineage's full signature
    is the accumulation root -> clade-root (shallowest node with that
    Nextclade_pango) with reversions/deletions applied in order. designationDate
    is absent from the tree and is carried over from `existing_dates`. Recombinant
    roots (X[A-Z]+) get an empty parent, matching the corneliusroemer schema."""
    from collections import defaultdict
    ref = tree["root_sequence"]["nuc"]

    def ref_base(pos: int) -> str:
        return ref[pos - 1]

    clade: dict[str, dict] = {}

    def walk(node, subs, dels, parent_lin):
        subs = dict(subs); dels = set(dels)
        branch_subs = []; branch_dels = set()
        for m in node.get("branch_attrs", {}).get("mutations", {}).get("nuc", []):
            mm = _TREE_MUT.match(m)
            if not mm:
                continue
            _, pos, alt = mm.group(1), int(mm.group(2)), mm.group(3)
            if alt == "N":
                continue
            if alt == "-":
                dels.add(pos); subs.pop(pos, None); branch_dels.add(pos)
            elif alt == ref_base(pos):
                subs.pop(pos, None); dels.discard(pos)
            else:
                subs[pos] = alt; dels.discard(pos); branch_subs.append((pos, alt))
        pango = node.get("node_attrs", {}).get("Nextclade_pango", {}).get("value")
        if pango and pango not in clade:
            clade[pango] = {
                "subs": dict(subs), "dels": set(dels),
                "new_subs": list(branch_subs), "new_dels": set(branch_dels),
                "parent": parent_lin,
                "nsClade": node.get("node_attrs", {}).get("clade_nextstrain", {}).get("value", ""),
                "alias": node.get("node_attrs", {}).get("partiallyAliased", {}).get("value", ""),
            }
        for child in node.get("children", []):
            walk(child, subs, dels, pango if pango else parent_lin)

    walk(tree["tree"], {}, set(), "")
    lineages = set(clade)
    existing_dates = existing_dates or {}

    def parent_of(lin: str, tree_parent: str) -> str:
        if _RECOMB.match(lin):
            return ""
        return tree_parent if (tree_parent in lineages and tree_parent) else ""

    def _pos(mstr: str) -> int:
        return int(re.match(r"^[ACGTN-](\d+)", mstr).group(1))

    children: dict[str, list] = defaultdict(list)
    out: dict[str, dict] = {}
    for lin, c in clade.items():
        parent = parent_of(lin, c["parent"])
        if parent:
            children[parent].append(lin)
        out[lin] = {
            "lineage": lin,
            "unaliased": c.get("alias") or lin,
            "parent": parent,
            "children": [],
            "nextstrainClade": c.get("nsClade", ""),
            "nucSubstitutions": sorted(
                (f"{ref_base(p)}{p}{a}" for p, a in c["subs"].items()), key=_pos),
            "nucSubstitutionsNew": sorted(
                (f"{ref_base(p)}{p}{a}" for p, a in c["new_subs"]), key=_pos),
            "nucDeletions": _nuc_ranges(c["dels"]),
            "nucDeletionsNew": _nuc_ranges(c["new_dels"]),
            "designationDate": existing_dates.get(lin),
        }
    for parent, kids in children.items():
        if parent in out:
            out[parent]["children"] = sorted(kids)
    return out


def download_pango_summary(local_path: str | Path) -> dict:
    """
    Fetch the Nextclade reference tree and write it to local_path as a
    pango_summary-compatible file. Replaces the corneliusroemer download
    (upstream frozen since 2025-06, no PJ.2+) and the old UShER/Freyja merge
    (which introduced orphaned lineages that quarantined the cache).
    designationDate is carried over from the file already at local_path.

    Returns:
        {success, new_variants, old_variants, added, error}
    """
    local = Path(local_path)

    old_variants: set[str] = set()
    existing_dates: dict[str, str] = {}
    if local.exists():
        try:
            with local.open("r", encoding="utf-8") as f:
                cur = json.load(f)
            old_variants = set(cur)
            existing_dates = {
                lin: entry.get("designationDate")
                for lin, entry in cur.items()
                if isinstance(entry.get("designationDate"), str)
            }
        except Exception:
            pass

    logging.info("pango_loader: fetching Nextclade reference tree -> pango_summary")
    try:
        with urllib.request.urlopen(NEXTCLADE_TREE_URL, timeout=60) as resp:
            tree = json.loads(resp.read())

        new_data = build_summary_from_tree(tree, existing_dates)
        new_variants = set(new_data)

        local.parent.mkdir(parents=True, exist_ok=True)
        tmp = local.with_suffix(".tmp")
        with tmp.open("w", encoding="utf-8") as f:
            json.dump(new_data, f)
        shutil.move(str(tmp), str(local))

        return {
            "success": True,
            "new_variants": len(new_variants),
            "old_variants": len(old_variants),
            "added": sorted(new_variants - old_variants),
            "error": None,
        }

    except Exception as exc:
        return {
            "success": False,
            "new_variants": 0,
            "old_variants": len(old_variants),
            "added": [],
            "error": str(exc),
        }