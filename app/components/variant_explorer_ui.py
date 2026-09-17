"""components/variant_explorer_ui.py

Trimmed on-demand variant explorer: pick any variant and see WHY it is or isn't
co-occurrence-detectable, plus its immediate connections (parent, near-identical
siblings). Uses a searchable selectbox (no free-typing / typos). The heavier
"nearest detectable relative" and per-city "found in data" were dropped to keep
it simple and avoid the subtle bugs those introduced.
"""
from __future__ import annotations

import streamlit as st

from process.variant_explorer import investigate_variant


def render_variant_explorer(pango_loader, panel=None, options=None,
                            disabled=False, key_prefix="acooc_explorer"):
    """Render the trimmed 'Investigate a variant' section.

    Args:
        pango_loader: provides get_raw_data() / get_signature().
        panel: current panel variant names (for the in-panel flag).
        options: list of selectable lineage names (defaults to all known).
        disabled: when True (e.g. a scan is still running) the picker is disabled
                  to avoid rerun races with the scanner autorefresh.
    """
    st.markdown("#### 🔎 Investigate a variant")
    st.caption("Pick any variant to see why co-occurrence can or can't detect it, "
               "and how it relates to others — including blind spots.")

    if disabled:
        st.info("Available once the scan completes.")
        return

    if options is None:
        options = sorted(pango_loader.get_raw_data().keys())

    variant = st.selectbox(
        "Variant", options=[""] + options, index=0,
        format_func=lambda v: "Search a lineage…" if v == "" else v,
        key=f"{key_prefix}_pick", label_visibility="collapsed",
    )
    if not variant:
        return

    r = investigate_variant(variant, pango_loader, panel=panel)
    if not r.get("found"):
        st.info(r.get("reason", "Not found."))
        return

    # header
    if r["detectable"]:
        badge = ("<span style='background:#e6f4ef;color:#0f6e56;font-size:11px;"
                 "font-weight:600;padding:1px 8px;border-radius:10px;'>&#10003; detectable</span>")
    else:
        badge = ("<span style='background:#f3f4f6;color:#6b7280;font-size:11px;"
                 "font-weight:600;padding:1px 8px;border-radius:10px;'>&#183; blind spot</span>")
    recomb = ("<span style='background:#f3e8f1;color:#8a4e82;font-size:11px;"
              "padding:1px 8px;border-radius:10px;margin-left:4px;'>recombinant</span>"
              if r["is_recombinant"] else "")
    panel_b = ("<span style='background:#e6f0f9;color:#185fa5;font-size:11px;"
               "padding:1px 8px;border-radius:10px;margin-left:4px;'>in panel</span>"
               if r["in_panel"] else "")
    st.markdown(
        f"<div style='font-size:15px;font-weight:600;margin-top:4px;'>{r['name']} "
        f"{badge}{recomb}{panel_b}</div>", unsafe_allow_html=True)

    # detectability + connections (trimmed — no nearest-detectable, no found-in-data)
    lines = [f"<b>Detectability:</b> {r['reason']}"]
    if not r["detectable"]:
        lines.append("<b>Quantify with:</b> deconvolution (not co-occurrence)")
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
    st.markdown(
        "<div style='border:0.5px solid #e5e7eb;border-radius:8px;padding:8px 12px;"
        "font-size:12.5px;line-height:1.9;margin-top:4px;'>" + "<br>".join(lines)
        + "</div>",
        unsafe_allow_html=True,
    )
