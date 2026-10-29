"""“ⓘ How to read this page”: one popover button above the result tabs
(2026-10-13). What the page shows, in its own words: every colour, label and
square of the Lineages table, its details and the panel check graph
(variants_table.legend_full_html).

2026-10-28: the "Rules and thresholds" tab is gone (it described the
algorithm, not the page, and still listed the old marker vote); the rules,
their settings and how they were tested are in RULES.md.

The Lineages table keeps a one-line legend on top (variants_table.SHORT_LEGEND)
so the colours stay in view while reading the tree."""
from __future__ import annotations

import streamlit as st

_GRAPH = ("<div style='font-size:13px;color:#4b5563;margin-top:10px;'>"
          "<b style='font-size:12px;color:#6b7280;text-transform:uppercase;'>Panel check graph</b>"
          " &nbsp;green = explained by the panel · light green = panel ± 1 change · "
          "red = a new lineage the reads suggest · grey = not explained (too broad to name, "
          "new combinations, errors)</div>")


def render_how_to_read() -> None:
    """The ⓘ popover button. Call it in the row of the result tabs."""
    with st.popover("ⓘ How to read this page", use_container_width=True):
        from components.variants_table import legend_full_html
        st.markdown(legend_full_html() + _GRAPH, unsafe_allow_html=True)
