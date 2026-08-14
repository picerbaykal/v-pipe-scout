"""Cached signature loaders for the synchronous scanner.

The scanner classifies unexplained patterns against known lineage signatures.
Loading all ~5000 pango signatures takes a few seconds, so we cache them
with st.cache_resource — they load once per Streamlit session and stay
cached across reruns.
"""

import re
from typing import Dict, Set

import streamlit as st


@st.cache_resource(show_spinner=False)
def get_all_lineage_signatures() -> Dict[str, Set[str]]:
    """All pango lineage signatures, keyed by lineage name.

    Format: {lineage_name: set of "{pos}{alt}" strings}
    Only lineages with >= 2 substitution mutations are included.
    Cached for the session — loads once (~a few seconds), then instant.
    """
    from api.pango_loader import PangoLoader, get_pango_summary_path
    pl = PangoLoader(get_pango_summary_path())
    raw = pl.get_raw_data()
    result: Dict[str, Set[str]] = {}
    for lineage in raw:
        try:
            sig = {m for m in pl.get_signature(lineage)
                   if re.match(r"^\d+[ACGT]$", m)}
            if len(sig) >= 2:
                result[lineage] = sig
        except Exception:
            continue
    return result


@st.cache_resource(show_spinner=False)
def get_panel_parent_map() -> Dict[str, str]:
    """Map each pango lineage to its parent lineage. Cached for the session."""
    from api.pango_loader import PangoLoader, get_pango_summary_path
    pl = PangoLoader(get_pango_summary_path())
    raw = pl.get_raw_data()
    return {lineage: data.get("parent", "") for lineage, data in raw.items()}


@st.cache_resource(show_spinner=False)
def get_cowwid_signatures() -> Dict[str, Set[str]]:
    """All cowwid surveillance variant signatures, keyed by variant name.

    Format: {variant_name: set of "{pos}{alt}" strings}
    Cached for the session.
    """
    from api.signatures import get_variant_list
    return {
        v.name: {m[1:] for m in v.signature_mutations if len(m) > 1}
        for v in get_variant_list().variants
    }