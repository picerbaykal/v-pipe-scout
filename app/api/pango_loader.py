from __future__ import annotations

import json
import logging
import shutil
import urllib.request
from pathlib import Path
from typing import cast
import logging

PANGO_DATA_DIR = Path(__file__).parent.parent / "data"
PANGO_SUMMARY_DEFAULT = PANGO_DATA_DIR / "pango_summary.json"
PANGO_SUMMARY_CACHE = Path("/app/.cache/pango/pango_summary.json")  # Docker runtime path

def get_pango_summary_path() -> Path:
    """
    Return the path to the best available pango_summary.json.

    Two copies may exist:
    - A checked-in default at app/data/pango_summary.json — always
      present, never overwritten at runtime. Used for local dev,
      CI, and as the fallback on a fresh deployment before any
      update check has run.
    - A runtime copy at /app/.cache/pango/pango_summary.json —
      written by download_pango_summary() on first startup and
      refreshed whenever upstream has a newer version. Lives in a
      Docker volume so it persists across container restarts.

    Prefers the runtime copy when it exists (i.e. after at least
    one successful download), falls back to the checked-in default
    otherwise.
    """
    if PANGO_SUMMARY_CACHE.exists():
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

        # initialize processed containers
        self._signatures = {}
        self._private_mutations = {}
        self._designation_dates = {}
        self._reconstructed_signatures: set[str] = set()  # empty-node parents filled from children
        self._process()
        self._fill_empty_node_signatures()


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
        if lineage in self._signatures:
            return self._signatures[lineage]
        # Lineage absent from pango_summary (e.g. BA.3.2) — try to reconstruct
        # from sub-lineages present in pango_summary via prefix matching.
        # Equivalent to querying CovSpectrum with "BA.3.2*".
        prefix = lineage + "."
        child_sigs = [
            sig for name, sig in self._signatures.items()
            if name.startswith(prefix) and sig
        ]
        if not child_sigs:
            self._signatures[lineage] = set()
            return set()
        intersection = child_sigs[0].intersection(*child_sigs[1:])
        self._signatures[lineage] = intersection
        if intersection:
            self._reconstructed_signatures.add(lineage)
        return intersection

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

def download_pango_summary(local_path: str | Path) -> dict:
    """
    Download the latest pango_summary.json and overwrite local_path.

    Returns:
        {
            "success": bool,
            "new_variants": int,   # variants in new file
            "old_variants": int,   # variants in old file (0 if didn't exist)
            "added": list[str],    # newly added variant names
            "error": str | None,
        }
    """

    local = Path(local_path)

    old_variants: set[str] = set()
    if local.exists():
        try:
            with local.open("r", encoding="utf-8") as f:
                old_data = json.load(f)
            old_variants = set(old_data.keys())
        except Exception:
            pass

    try:
        with urllib.request.urlopen(PANGO_SUMMARY_URL, timeout=30) as resp:
            raw = resp.read()
            remote_etag = resp.headers.get("ETag")

        # validate JSON before overwriting
        new_data = json.loads(raw)
        new_variants = set(new_data.keys())

        # atomic write via temp file
        local.parent.mkdir(parents=True, exist_ok=True)
        tmp = local.with_suffix(".tmp")
        tmp.write_bytes(raw)
        shutil.move(str(tmp), str(local))

        # save ETag sidecar so next check knows current version
        if remote_etag:
            local.with_suffix(".etag").write_text(remote_etag)

        added = sorted(new_variants - old_variants)
        return {
            "success": True,
            "new_variants": len(new_variants),
            "old_variants": len(old_variants),
            "added": added,
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



