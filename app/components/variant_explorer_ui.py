"""components/variant_explorer_ui.py

UI for the on-demand variant explorer: a search box + a fact-sheet + a focused
connections view. Renders the result of process.variant_explorer.investigate_variant.
Pure rendering — all logic lives in process/variant_explorer.py.
"""
from __future__ import annotations

import streamlit as st

from process.variant_explorer import investigate_variant


def render_variant_explorer(pango_loader, panel=None, cooc_result=None,
                            key_prefix="acooc_explorer"):
    """Render the 'Investigate a variant' section.

    Args:
        pango_loader: provides get_raw_data() / get_signature().
        panel: current panel variant names (for the in-panel flag).
        cooc_result: optional cooc result for the active location (found-in-data).
    """
    st.markdown("#### 🔎 Investigate a variant")
    st.caption(
        "Look up any variant's co-occurrence detectability and how it connects to "
        "others — including blind spots the scanner can't surface."
    )

    query = st.text_input(
        "Variant", placeholder="e.g. KP.2, XFG, NB.1.8.1",
        key=f"{key_prefix}_query", label_visibility="collapsed",
    )
    if not query:
        return
    variant = query.strip()

    r = investigate_variant(variant, pango_loader, panel=panel,
                            cooc_result=cooc_result)
    if not r.get("found"):
        st.info(r.get("reason", "Not found."))
        return

    # ── header ──
    if r["detectable"]:
        badge = ("<span style='background:#e6f4ef;color:#0f6e56;font-size:11px;"
                 "font-weight:600;padding:1px 8px;border-radius:10px;'>✓ detectable</span>")
    else:
        badge = ("<span style='background:#f3f4f6;color:#6b7280;font-size:11px;"
                 "font-weight:600;padding:1px 8px;border-radius:10px;'>· blind spot</span>")
    recomb = ("<span style='background:#f3e8f1;color:#8a4e82;font-size:11px;"
              "padding:1px 8px;border-radius:10px;margin-left:4px;'>recombinant</span>"
              if r["is_recombinant"] else "")
    panel_b = ("<span style='background:#e6f0f9;color:#185fa5;font-size:11px;"
               "padding:1px 8px;border-radius:10px;margin-left:4px;'>in panel</span>"
               if r["in_panel"] else "")
    st.markdown(
        f"<div style='font-size:15px;font-weight:600;margin-top:4px;'>{r['name']} "
        f"{badge}{recomb}{panel_b}</div>", unsafe_allow_html=True)

    # ── fact sheet ──
    _found = ("yes — %d reads" % r["found_reads"]) if r["found_in_data"] else "not distinctly"
    st.markdown(
        f"<div style='font-size:12.5px;line-height:1.7;margin:6px 0;'>"
        f"<b>Detectability:</b> {r['reason']}<br>"
        f"<b>In panel:</b> {'yes' if r['in_panel'] else 'no'} · "
        f"<b>Found in data:</b> {_found}"
        + ("" if r["detectable"] else
           " · <b>Quantify with:</b> deconvolution (not co-occurrence)")
        + "</div>",
        unsafe_allow_html=True,
    )

    # ── connections ──
    lines = []
    if r["parent"]:
        lines.append(f"<span style='color:#6b7280;'>Belongs to:</span> <b>{r['parent']}</b>")
    if r["oscillates_with"]:
        osc = ", ".join(r["oscillates_with"][:6])
        lines.append(
            f"<span style='color:#6b7280;'>Near-identical to (can't be told apart):</span> "
            f"<span style='background:#fdf4e6;color:#ba7517;padding:0 6px;"
            f"border-radius:8px;'>{osc}</span>")
    if r["children"]:
        lines.append(f"<span style='color:#6b7280;'>Descendants:</span> "
                     + ", ".join(r["children"][:6]))
    if not r["detectable"]:
        if r["nearest_detectable"]:
            lines.append(
                f"<span style='color:#6b7280;'>Detectable variant in this branch "
                f"(trackable instead):</span> "
                f"<span style='color:#0f6e56;font-weight:600;'>{r['nearest_detectable']}</span>")
        else:
            lines.append("<span style='color:#6b7280;'>No detectable variant in "
                         "this branch — the whole cluster is co-occurrence-blind.</span>")
    if lines:
        st.markdown(
            "<div style='border:0.5px solid #e5e7eb;border-radius:8px;padding:8px 12px;"
            "font-size:12.5px;line-height:1.9;'>" + "<br>".join(lines) + "</div>",
            unsafe_allow_html=True,
        )