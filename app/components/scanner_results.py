"""Scanner results component.

Renders the three-bucket output of run_cooc_scanner_lapis:
  - Missing from panel (cowwid tracked variants not selected)
  - Emerging sublineages (descendants of panel variants)
  - Possibly new (no known lineage explains the pattern)

Each bucket shows actionable information: which variants to add,
which sublineages to watch, and the read evidence behind each finding.
"""

import streamlit as st


def render_scanner_results(
    result: dict,
    selected_variants: list,
    on_add_variant=None,
) -> None:
    """
    Render the scanner classification results.

    Args:
        result: Output of run_cooc_scanner_lapis — dict with keys
            missing_from_panel, emerging_sublineage, possibly_new,
            total_unexplained_reads, summary.
        selected_variants: Currently selected panel variants (for context).
        on_add_variant: Optional callback(variant_name) when user clicks
            "Add to panel". If None, shows the variant name only.
    """
    if not result:
        st.caption("No scanner results available.")
        return

    total = result.get("total_unexplained_reads", 0)
    summary = result.get("summary", "")

    if total == 0:
        st.success("No unexplained patterns found — panel looks complete.")
        return

    st.caption(f"{total:,} unexplained reads classified · {summary}")
    st.markdown("---")

    # ── Bucket 1: missing from panel ─────────────────────────────────────────
    missing = result.get("missing_from_panel", [])
    if missing:
        st.markdown("#### 🔴 Missing from panel")
        st.caption(
            "These officially tracked (cowwid) variants are not in your panel "
            "but explain unexplained reads. Add them to improve completeness."
        )
        for item in missing:
            variant = item["variant"]
            reads = item["total_reads"]
            n_patterns = item["pattern_count"]
            col1, col2 = st.columns([4, 1])
            with col1:
                st.markdown(
                    f"<div style='background:#fef2f2; border:1px solid #fecaca; "
                    f"border-radius:6px; padding:8px 12px; margin:4px 0;'>"
                    f"<span style='font-weight:600; color:#dc2626;'>{variant}</span>"
                    f"<span style='color:#6b7280; font-size:0.82rem; margin-left:8px;'>"
                    f"{reads:,} reads · {n_patterns} pattern(s)</span>"
                    f"</div>",
                    unsafe_allow_html=True,
                )
            with col2:
                if on_add_variant:
                    if st.button(
                        "＋ Add",
                        key=f"scanner_add_{variant}",
                        use_container_width=True,
                        type="primary",
                    ):
                        on_add_variant(variant)
        st.markdown("")

    # ── Bucket 2: emerging sublineages ───────────────────────────────────────
    emerging = result.get("emerging_sublineage", [])
    if emerging:
        st.markdown("#### 🟡 Emerging sublineages")
        st.caption(
            "Descendants of your panel variants with rising co-occurrence signal. "
            "Not yet on the official surveillance list."
        )
        for item in emerging[:10]:  # cap display at 10
            lineage = item["lineage"]
            parent = item["parent"]
            reads = item["total_reads"]
            n_patterns = item["pattern_count"]
            col1, col2 = st.columns([4, 1])
            with col1:
                st.markdown(
                    f"<div style='background:#fffbeb; border:1px solid #fde68a; "
                    f"border-radius:6px; padding:8px 12px; margin:4px 0;'>"
                    f"<span style='font-weight:600; color:#92400e;'>{lineage}</span>"
                    f"<span style='color:#6b7280; font-size:0.82rem; margin-left:8px;'>"
                    f"sublineage of {parent} · {reads:,} reads · {n_patterns} pattern(s)</span>"
                    f"</div>",
                    unsafe_allow_html=True,
                )
            with col2:
                if on_add_variant:
                    if st.button(
                        "＋ Add",
                        key=f"scanner_add_{lineage}",
                        use_container_width=True,
                    ):
                        on_add_variant(lineage)
        if len(emerging) > 10:
            st.caption(f"… and {len(emerging) - 10} more sublineages.")
        st.markdown("")

    # ── Bucket 3: possibly new ───────────────────────────────────────────────
    possibly_new = result.get("possibly_new", {})
    pn_reads = possibly_new.get("total_reads", 0)
    pn_patterns = possibly_new.get("pattern_count", 0)
    top_patterns = possibly_new.get("top_patterns", [])

    if pn_reads > 0:
        st.markdown("#### 🔵 Possibly new")
        st.caption(
            f"{pn_reads:,} reads in {pn_patterns} pattern(s) match no known lineage. "
            "Could be a novel variant, recombinant, or sequencing artifact."
        )
        if top_patterns:
            with st.expander("Show top patterns", expanded=False):
                for pat in top_patterns[:5]:
                    mutations = ", ".join(pat["mutations"])
                    st.markdown(
                        f"<div style='background:#eff6ff; border:1px solid #bfdbfe; "
                        f"border-radius:6px; padding:6px 10px; margin:3px 0; "
                        f"font-size:0.82rem;'>"
                        f"<span style='color:#1d4ed8; font-family:monospace;'>{mutations}</span>"
                        f"<span style='color:#6b7280; margin-left:8px;'>"
                        f"{pat['count']:,} reads · {pat['date']}</span>"
                        f"</div>",
                        unsafe_allow_html=True,
                    )