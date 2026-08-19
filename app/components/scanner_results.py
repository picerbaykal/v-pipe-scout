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
            return

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


def _jaccard(sig_a: set, sig_b: set) -> float:
    """Jaccard similarity between two signature sets."""
    if not sig_a or not sig_b:
        return 0.0
    return len(sig_a & sig_b) / len(sig_a | sig_b)


def _similarity_warning(
    candidate: str,
    candidate_sig: set,
    panel_variants: list,
    panel_sigs: dict,
) -> tuple[str, str, float]:
    """
    Check how similar a candidate is to all current panel variants.

    Returns (level, message, max_jaccard) where level is:
        "blocked"  — Jaccard >= 0.90 or one contains the other (pvt=0)
        "caution"  — Jaccard 0.85–0.90
        "safe"     — Jaccard < 0.85
    """
    max_j = 0.0
    worst_variant = ""
    worst_pvt_candidate = 0
    worst_pvt_panel = 0

    for pv in panel_variants:
        sig_p = panel_sigs.get(pv, set())
        if not sig_p or not candidate_sig:
            continue
        j = _jaccard(candidate_sig, sig_p)
        if j > max_j:
            max_j = j
            worst_variant = pv
            worst_pvt_candidate = len(candidate_sig - sig_p)
            worst_pvt_panel = len(sig_p - candidate_sig)

    if not worst_variant:
        return "safe", "", 0.0

    # containment: one signature is fully inside the other
    contained = (worst_pvt_candidate == 0 or worst_pvt_panel == 0)

    if contained or max_j >= 0.90:
        if contained:
            msg = (
                f"⚠️ **{candidate}** fully overlaps with **{worst_variant}** "
                f"already in your panel — adding both is redundant. "
                f"Choose one or the other."
            )
        else:
            msg = (
                f"⚠️ **{candidate}** is very similar to **{worst_variant}** "
                f"(Jaccard {max_j:.2f}, {worst_pvt_candidate} private mutations). "
                f"Adding both will destabilize deconvolution — choose one."
            )
        return "blocked", msg, max_j

    if max_j >= 0.85:
        msg = (
            f"⚠️ **{candidate}** is moderately similar to **{worst_variant}** "
            f"(Jaccard {max_j:.2f}). Adding both may widen confidence intervals. "
            f"{candidate} has {worst_pvt_candidate} private mutations."
        )
        return "caution", msg, max_j

    return "safe", "", max_j


def _get_panel_sigs(panel_variants: list) -> dict:
    """Load pango signatures for all panel variants. Cached in session state."""
    import re
    cache_key = "scanner_panel_sigs_cache"
    cached = st.session_state.get(cache_key, {})
    missing = [v for v in panel_variants if v not in cached]
    if missing:
        try:
            from api.pango_loader import PangoLoader, get_pango_summary_path
            pl = PangoLoader(get_pango_summary_path())
            for v in missing:
                try:
                    s = pl.get_signature(v)
                    cached[v] = {m for m in s if re.match(r'^\d+[ACGT]$', m)}
                except Exception:
                    cached[v] = set()
            st.session_state[cache_key] = cached
        except Exception:
            pass
    return cached


def _add_button_with_jaccard(
    candidate: str,
    candidate_sig: set,
    panel_variants: list,
    panel_sigs: dict,
    location: str,
    on_add_variant,
    key_suffix: str = "",
) -> None:
    """
    Render an Add button with Jaccard-based guidance.

    - blocked  (>=0.90 or containment): show warning, no Add button.
    - caution  (0.85-0.90): show warning + Add button (user decides).
    - safe     (<0.85): show Add button, optionally show similarity.
    """
    if not on_add_variant:
        return

    level, msg, max_j = _similarity_warning(
        candidate, candidate_sig, panel_variants, panel_sigs
    )
    safe_key = f"scanner_add_{location}_{candidate.replace(' ', '_').replace('(', '').replace(')', '')}{key_suffix}"

    if level == "blocked":
        st.warning(msg)
    elif level == "caution":
        st.warning(msg)
        if st.button("＋ Add anyway", key=safe_key, use_container_width=True):
            on_add_variant(candidate)
    else:
        # safe — show similarity as subtle info if non-trivial
        if max_j >= 0.70:
            st.caption(f"Similarity to nearest panel variant: {max_j:.2f}")
        if st.button("＋ Add", key=safe_key, use_container_width=True, type="primary"):
            on_add_variant(candidate)


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

    # load panel signatures once for Jaccard checks
    panel_sigs = _get_panel_sigs(selected_variants)

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
                import re
                variant = item["variant"]
                reads = item["total_reads"]
                n_patterns = item["pattern_count"]

                # get candidate signature for Jaccard check
                try:
                    from api.pango_loader import PangoLoader, get_pango_summary_path
                    pl = PangoLoader(get_pango_summary_path())
                    raw_sig = pl.get_signature(variant)
                    candidate_sig = {m for m in raw_sig if re.match(r'^\d+[ACGT]$', m)}
                except Exception:
                    candidate_sig = set()

                col1, col2 = st.columns([4, 1])
                with col1:
                    st.markdown(
                        f"<div style='background:#fef2f2; border:1px solid #fecaca; "
                        f"border-radius:6px; padding:8px 12px; margin:4px 0;'>"
                        f"<span style='font-weight:600; color:#dc2626;'>{variant}</span>"
                        f"<span style='color:#6b7280; font-size:0.82rem; margin-left:8px;' "
                        f"title='{reads:,} reads carry patterns matching this variant "
                        f"(may overlap with related variants) · "
                        f"{n_patterns} distinct co-occurrence pattern(s)'>"
                        f"{reads:,} reads · {n_patterns} pattern(s)</span>"
                        f"</div>",
                        unsafe_allow_html=True,
                    )
                with col2:
                    _add_button_with_jaccard(
                        candidate=variant,
                        candidate_sig=candidate_sig,
                        panel_variants=selected_variants,
                        panel_sigs=panel_sigs,
                        location=location,
                        on_add_variant=on_add_variant,
                    )

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

            sorted_groups = sorted(
                by_parent.items(),
                key=lambda x: -sum(i["total_reads"] for i in x[1])
            )

            for parent_variant, items in sorted_groups:
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
                        obs_muts = item.get("observed_mutations", [])

                        st.markdown(
                            f"<div style='background:#fffbeb; border:1px solid #fde68a; "
                            f"border-radius:6px; padding:8px 12px; margin:4px 0;'>"
                            f"<span style='font-weight:600; color:#92400e;'>{lineage}</span>"
                            f"<span style='color:#6b7280; font-size:0.82rem; margin-left:8px;'>"
                            f"sublineage of {parent} · {reads:,} reads · {n_patterns} pattern(s)</span>"
                            f"</div>",
                            unsafe_allow_html=True,
                        )
                        # surveillance-only message — no Add button
                        n_private = len(obs_muts)
                        st.caption(
                            f"📋 Surveillance finding — counted within **{parent}**'s proportion. "
                            f"{n_private} private mutation(s) observed: "
                            f"{', '.join(obs_muts[:5])}{'…' if n_private > 5 else ''}. "
                            f"Not recommended to add separately (Jaccard ≥ 0.90 with parent)."
                        )

                        if client and obs_muts and date_range:
                            from components.scanner_heatmap import render_scanner_heatmap
                            with st.expander(f"Signal over time — {location} — {lineage}", expanded=False):
                                render_scanner_heatmap(
                                    variant=lineage,
                                    mutations=obs_muts,
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