"""Scanner results component.

Renders the three-bucket output of run_cooc_scanner_lapis:
  - Missing from panel (cowwid tracked variants not selected)
  - Emerging sublineages (descendants of panel variants)
  - Possibly new (no known lineage explains the pattern)

Each bucket is collapsible. Missing from panel is expanded by default
since it's the most actionable. Each missing variant has a heatmap
showing its mutation signal over time.
"""

import streamlit as st


def render_scanner_progress(location: str, task_id: str, celery_app, redis_client) -> None:
    """
    Render progress bar while the scanner task is running.
    Follows the same pattern as render_location_progress.
    """
    import json

    try:
        task = celery_app.AsyncResult(task_id)
        if task.ready():
            return  # caller handles completion

        # try to get progress from Redis
        progress_key = f"task_progress:{task_id}"
        progress_data = redis_client.get(progress_key)
        if progress_data:
            progress = json.loads(progress_data)
            current = progress.get("current", 0)
            total = progress.get("total", 3)
            status = progress.get("status", "Running scanner…")
            st.progress(current / total if total > 0 else 0, text=status)
        else:
            st.progress(0, text="Scanner starting…")

    except Exception:
        st.info("Scanner running…")


def render_scanner_results(
    result: dict,
    selected_variants: list,
    client=None,
    location: str = "",
    date_range: tuple = None,
    on_add_variant=None,
) -> None:
    """
    Render the scanner classification results with expandable sections.

    Args:
        result: Output of run_cooc_scanner_lapis.
        selected_variants: Currently selected panel variants.
        client: WiseLoculusLapis instance for heatmap queries.
        location: Location name for heatmap queries.
        date_range: (start_datetime, end_datetime) for heatmap queries.
        on_add_variant: Optional callback(variant_name) for Add button.
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

    missing = result.get("missing_from_panel", [])
    emerging = result.get("emerging_sublineage", [])
    possibly_new = result.get("possibly_new", {})
    pn_reads = possibly_new.get("total_reads", 0)

    # ── Bucket 1: missing from panel ─────────────────────────────────────────
    with st.expander(
        f"🔴 Missing from panel ({len(missing)} variant{'s' if len(missing) != 1 else ''})",
        expanded=False,
    ):
        if not missing:
            st.caption("No missing tracked variants detected.")
        else:
            st.caption(
                "Officially tracked variants not in your panel that explain unexplained reads. "
                "Add them to improve completeness."
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
                        f"<span style='color:#6b7280; font-size:0.82rem; margin-left:8px;' "
                        f"title='{reads:,} reads carry patterns matching this sublineage "
                        f"(may overlap with related sublineages) · "
                        f"{n_patterns} distinct co-occurrence pattern(s)'>"
                        f"{reads:,} reads · {n_patterns} pattern(s)</span>"
                        f"</div>",
                        unsafe_allow_html=True,
                    )
                with col2:
                    if on_add_variant:
                        if st.button(
                            "＋ Add",
                                key=f"scanner_add_{location}_{variant.replace(' ', '_').replace('(', '').replace(')', '')}",
                                use_container_width=True,
                            type="primary",
                        ):
                            on_add_variant(variant)

                # heatmap expander per variant
                if client and item.get("observed_mutations") and date_range:
                    from components.scanner_heatmap import render_scanner_heatmap
                    with st.expander(f"Signal over time — {location} — {variant}", expanded=False):
                        render_scanner_heatmap(
                            variant=variant,
                            mutations=item["observed_mutations"],
                            client=client,
                            location=location,
                            date_range=date_range,
                            max_mutations=20,
                        )


    # ── Bucket 2: emerging sublineages ───────────────────────────────────────
    with st.expander(
        f"🟡 Emerging sublineages ({len(emerging)} found)",
        expanded=False,
    ):
        if not emerging:
            st.caption("No emerging sublineages detected.")
        else:
            st.caption(
                "Descendants of your panel variants with rising co-occurrence signal. "
                "Not yet on the official surveillance list. "
                "Sorted by read count — note that closely related sublineages sharing "
                "the same mutations may have overlapping counts."
            )
            from collections import defaultdict
            by_parent = defaultdict(list)
            for item in emerging:
                by_parent[item["parent"]].append(item)

            # sort groups by total reads (strongest signal first)
            sorted_groups = sorted(
                by_parent.items(),
                key=lambda x: -sum(i["total_reads"] for i in x[1])
            )

            for parent_variant, items in sorted_groups:
                # sort sublineages within group by reads descending
                items_sorted = sorted(items, key=lambda x: -x["total_reads"])
                group_reads = sum(i["total_reads"] for i in items_sorted)

                with st.expander(
                        f"**{parent_variant}** sublineages — {len(items_sorted)} found · {group_reads:,} reads",
                        expanded=False,
                ):
                    for item in items_sorted[:10]:
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
                                        key=f"scanner_add_{location}_{lineage.replace(' ', '_').replace('(', '').replace(')', '')}",
                                        use_container_width=True,
                                ):
                                    on_add_variant(lineage)

                        if client and item.get("observed_mutations") and date_range:
                            from components.scanner_heatmap import render_scanner_heatmap
                            with st.expander(f"Signal over time — {location} — {lineage}", expanded=False):
                                render_scanner_heatmap(
                                    variant=lineage,
                                    mutations=item["observed_mutations"],
                                    client=client,
                                    location=location,
                                    date_range=date_range,
                                    max_mutations=20,
                            )
                if len(items) > 10:
                    st.caption(f"… and {len(items_sorted) - 10} more {parent_variant} sublineages.")
    # ── Bucket 3: possibly new ───────────────────────────────────────────────
    with st.expander(
        f"🔵 Possibly new ({pn_reads:,} reads)",
        expanded=False,
    ):
        if pn_reads == 0:
            st.caption("No unexplained patterns that don't match any known lineage.")
        else:
            pn_patterns = possibly_new.get("pattern_count", 0)
            top_patterns = possibly_new.get("top_patterns", [])
            st.caption(
                f"{pn_reads:,} reads in {pn_patterns} pattern(s) match no known lineage. "
                "Could be a novel variant, recombinant, or sequencing artifact."
            )
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