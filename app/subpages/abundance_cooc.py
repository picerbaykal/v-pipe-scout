"""
Abundance & Co-occurrence tab.

Estimates variant abundance over time using LAPIS-sourced wastewater
mutation data, with LolliPop deconvolution (including bootstrap confidence
intervals) and pango-lineage-based sublineage enrichment.

Co-occurrence features (panel-completeness signal, emerging-sublineage and
de-novo variant detection) are planned but not yet implemented — they
depend on a SaneQL/LAPIS endpoint not yet deployed on the WASAP SILO
instance. The tree view's co-occurrence controls are present but disabled
until that lands.
"""

import streamlit as st
import os
from celery import Celery
import redis
from api.pango_loader import PangoLoader, get_pango_summary_path
import logging
logger = logging.getLogger(__name__)

from api.wiseloculus import WiseLoculusLapis
from utils.config import get_wiseloculus_url

from components.abundance_cooc_tree import render_panel_tree
from components.jaccard_heatmap import render_jaccard_heatmap

# Initialize Celery
celery_app = Celery(
    'tasks',
    broker=os.environ.get('CELERY_BROKER_URL', 'redis://redis:6379/0'),
    backend=os.environ.get('CELERY_RESULT_BACKEND', 'redis://redis:6379/0')
)

# Initialize Redis client for checking task status
redis_client = redis.Redis(
    host=os.environ.get('REDIS_HOST', 'redis'),
    port=int(os.environ.get('REDIS_PORT', 6379)),
    password=os.environ.get('REDIS_PASSWORD', 'defaultpassword123'),
    db=0
)

wise_server_ip = get_wiseloculus_url()
wiseLoculus = WiseLoculusLapis(wise_server_ip)


@st.cache_resource
def cached_get_pango_loader() -> PangoLoader:
    """
    Return a cached PangoLoader instance, loaded once per app session.
    Uses get_pango_summary_path() to prefer the runtime copy when available.
    """
    return PangoLoader(get_pango_summary_path())

@st.cache_data
def cached_get_variant_names() -> list:
    """
    Return cached list of curated variant names from cowwid.
    Cached per session to avoid repeated GitHub API calls.
    """
    from api.signatures import get_variant_names
    return get_variant_names()

def _render_completeness(result: dict) -> None:
    """Plot per-date panel completeness with an informative-read count below."""
    import pandas as pd
    import plotly.graph_objects as go

    if not result.get("dates"):
        st.warning("No co-occurrence data for this location and range.")
        return

    df = pd.DataFrame({
        "date": pd.to_datetime(result["dates"]),
        "matched": result["matched_counts"],
        "unexplained": result["unexplained_counts"],
        "completeness": result["completeness"],
    })
    df["informative"] = df["matched"] + df["unexplained"]

    fig = go.Figure(go.Scatter(
        x=df["date"], y=df["completeness"],
        mode="lines+markers",
        marker=dict(size=[6 + min(8, n / 5000) for n in df["informative"]]),
        line=dict(width=1.5),
        customdata=df[["informative"]],
        hovertemplate="%{x|%Y-%m-%d}<br>completeness %{y:.1%}"
                      "<br>%{customdata[0]:,} informative reads<extra></extra>",
    ))
    fig.update_yaxes(range=[-0.05, 1.05], tickformat=".0%", title="completeness")
    fig.update_layout(height=200, margin=dict(t=10, b=30, l=50, r=20),
                      template="plotly_white")
    st.plotly_chart(fig, use_container_width=True)
    st.caption(
        f"Marker size scales with informative reads "
        f"({df['informative'].min():,}–{df['informative'].max():,} across dates). "
        "Sparse dates give unstable values."
    )



def app():
    st.session_state.setdefault("location_results", {})
    st.session_state.setdefault("acooc_location_tasks", {})
    st.session_state.setdefault("acooc_cooc_tasks", {})  # add
    st.session_state.setdefault("acooc_cooc_results", {})  # add

    # ── Header ───────────────────────────────────────────────────────────────
    st.title("Abundance & Co-occurrence")
    st.subheader(
        "Estimate the proportion of variants circulating in wastewater over time, "
        "using live LAPIS mutation data and LolliPop deconvolution with bootstrap "
        "confidence intervals."
    )
    st.write(
        "Build a custom variant panel — curated variants of interest, high-signal "
        "sublineages, or any pango lineage — and run on-demand deconvolution across "
        "one or more sampling locations. Co-occurrence features (panel-completeness "
        "check, emerging-sublineage and de-novo variant detection) are planned but "
        "not yet available — they depend on a LAPIS/SaneQL endpoint not yet deployed "
        "on this instance."
    )

    st.markdown("---")

    # ── Two-column layout ─────────────────────────────────────────────────────
    col_controls, col_results = st.columns([1, 3])

    with col_controls:

        # ── Variant Panel ────────────────────────────────────────────────────────
        st.subheader("Variant Panel")

        # Surveillance panel - cowwid curated variants
        curated_variants = cached_get_variant_names()

        col_select_all, col_clear = st.columns(2)
        with col_select_all:
            if st.button("Select all", key="acooc_select_all", use_container_width=True):
                st.session_state["acooc_variant_multiselect"] = curated_variants
        with col_clear:
            if st.button("Clear", key="acooc_clear_variants", use_container_width=True):
                st.session_state["acooc_variant_multiselect"] = []

        selected_variants = st.multiselect(
            "Surveillance variants",
            options=curated_variants,
            help="Curated variants from the Swiss wastewater surveillance panel (cowwid).",
            key="acooc_variant_multiselect",
        )

        # Manual add - for v0, user can add any lineage they already know about
        # In v1, this will be replaced/augmented by the scanner's + button on the tree
        with st.expander("➕ Add lineage manually", expanded=False):
            pango_loader = cached_get_pango_loader()
            available_lineages = sorted(pango_loader.get_raw_data().keys())
            extra_options = [v for v in available_lineages if v not in selected_variants]

            extra_variants = st.multiselect(
                "Search pango lineage",
                options=extra_options,
                placeholder="e.g. KP.2.3",
                help="Add any pango lineage to the panel. Strategy D enrichment applied automatically.",
                key="acooc_extra_variants",
            )

        # Combined panel
        all_selected_variants = selected_variants + extra_variants

        if len(all_selected_variants) < 2:
            st.warning("Select at least 2 variants to run deconvolution.")
        else:
            st.caption(f"{len(all_selected_variants)} variants in panel")

        st.markdown("---")

        # ── Variant Tree ──────────────────────────────────────────────────────
        st.subheader("Variant Tree")
        render_panel_tree(
            selected_variants=all_selected_variants,
            yaml_variants=curated_variants,
            pango_loader=cached_get_pango_loader(),
        )

        st.markdown("---")

        # ── Locations ─────────────────────────────────────────────
        st.subheader("Locations")

        available_locations = wiseLoculus.fetch_locations()

        selected_locations = st.multiselect(
            "Select sampling locations",
            options=available_locations,
            placeholder="Select one or more locations...",
            key="acooc_location_multiselect",
        )

        if not selected_locations:
            st.warning("Select at least one location.")

        st.markdown("---")

        # ── Date Range ────────────────────────────────────────────────────────
        st.subheader("Date Range")

        default_start, default_end, min_date, max_date = wiseLoculus.get_cached_date_range_with_bounds("abundance_cooc")

        col_start, col_end = st.columns(2)

        with col_start:
            start_date = st.date_input(
                "Start Date",
                value=default_start,
                min_value=min_date,
                max_value=max_date,
                key="acooc_start_date",
            )
        with col_end:
            end_date = st.date_input(
                "End Date",
                value=default_end,
                min_value=min_date,
                max_value=max_date,
                key="acooc_end_date",
            )

        if end_date <= start_date:
            st.warning("End date must be after start date.")

        # ── LolliPop Parameters ───────────────────────────────────────────────
        with st.expander("⚙️ LolliPop Parameters", expanded=False):
            bootstrap_options = {
                "Rapid (Fast)": 50,
                "Standard": 100,
                "Reliable (Slower)": 300,
            }
            selected_bootstrap = st.radio(
                "Bootstrap Iterations",
                options=list(bootstrap_options.keys()),
                index=1,
                help="More iterations = tighter confidence intervals but slower runtime.",
                key="acooc_bootstrap",
            )
            bootstraps = bootstrap_options[selected_bootstrap]
            st.caption(f"Selected: {bootstraps} bootstrap iterations")

            bandwidth_options = {
                "Narrow": 10,
                "Medium": 20,
                "Wide": 30,
            }
            selected_bandwidth = st.radio(
                "Bandwidth (Gaussian Kernel Smoothing)",
                options=list(bandwidth_options.keys()),
                index=0,
                help="Narrow preserves short-term variation, Wide smooths long-term trends.",
                key="acooc_bandwidth",
            )
            bandwidth = bandwidth_options[selected_bandwidth]
            st.caption(f"Selected: Bandwidth = {bandwidth}")

        st.markdown("---")

        # ── Run Button ────────────────────────────────────────────────────────
        can_run = (
                len(all_selected_variants) >= 2
                and len(selected_locations) >= 1
                and end_date > start_date
        )

        if st.button(
                "▶ Run",
                type="primary",
                disabled=not can_run,
                key="acooc_run_button",
                help="Select at least 2 variants and 1 location to run.",
        ):
            st.session_state["acooc_trigger_run"] = True

    with col_results:

        # ── Task submission ───────────────────────────────────────────────────
        if st.session_state.get("acooc_trigger_run"):
            st.session_state["acooc_trigger_run"] = False

            location_tasks = {}
            cooc_tasks = {}
            for loc in selected_locations:
                task = celery_app.send_task(
                    "tasks.run_deconvolve_lapis",
                    kwargs={
                        "location": loc,
                        "start_date": start_date.isoformat(),
                        "end_date": end_date.isoformat(),
                        "variants": all_selected_variants,
                        "bootstraps": bootstraps,
                        "bandwidth": bandwidth,
                    }
                )
                location_tasks[loc] = task.id

                cooc_task = celery_app.send_task(
                    "tasks.run_cooc_completeness_lapis",
                    kwargs={
                        "location": loc,
                        "start_date": start_date.isoformat(),
                        "end_date": end_date.isoformat(),
                        "variants": all_selected_variants,
                    }
                )

                cooc_tasks[loc] = cooc_task.id
                logger.info(f"Submitted task {task.id} for {loc}")

            st.session_state["acooc_location_tasks"] = location_tasks
            st.session_state["location_results"] = {}
            st.session_state["acooc_cooc_tasks"] = cooc_tasks
            st.session_state["acooc_cooc_results"] = {}
            st.rerun()

        # ── Results ───────────────────────────────────────────────────────────
        from components.multi_location_results import (
            render_single_location_result,
            render_location_progress,
        )

        location_tasks = st.session_state.get("acooc_location_tasks", {})
        location_results = st.session_state.get("location_results", {})

        if not location_tasks:
            st.info("Results will appear here after running deconvolution.")
        else:
            location_tasks = st.session_state.get("acooc_location_tasks", {})

            # Check task states
            cooc_task_ids = list(st.session_state.get("acooc_cooc_tasks", {}).values())
            all_task_ids = list(location_tasks.values()) + cooc_task_ids

            states = {
                task_id: celery_app.AsyncResult(task_id).state
                for task_id in all_task_ids
            }

            any_running = any(
                s in ("PENDING", "STARTED", "RETRY")
                for s in states.values()
            )

            n_results = (
                    len(st.session_state.location_results)
                    + len(st.session_state.get("acooc_cooc_results", {}))
            )
            any_just_completed = (
                    any(s == "SUCCESS" for s in states.values())
                    and n_results < len(all_task_ids)
            )

            if any_running or any_just_completed:
                from streamlit_autorefresh import st_autorefresh
                st_autorefresh(interval=2000, key="acooc_autorefresh")

                # ── Custom tab rendering with cooc/scanner placeholders ──────
            location_names = list(location_tasks.keys())
            tabs = st.tabs([f"📍 {loc}" for loc in location_names])

            for location, tab in zip(location_names, tabs):
                with tab:
                    task_id = location_tasks[location]

                    # 1. Cooc panel completeness
                    st.markdown("#### Panel completeness (co-occurrence)")
                    cooc_tasks = st.session_state.get("acooc_cooc_tasks", {})
                    cooc_results = st.session_state.get("acooc_cooc_results", {})

                    if location in cooc_results:
                        _render_completeness(cooc_results[location])
                    elif location in cooc_tasks:
                        cooc_task = celery_app.AsyncResult(cooc_tasks[location])
                        if cooc_task.ready():
                            try:
                                cooc_results[location] = cooc_task.get()
                                st.rerun()
                            except Exception as e:
                                st.error(f"Co-occurrence failed: {e}")
                        else:
                            st.info("Computing panel completeness…")
                    else:
                        st.caption("Run to compute panel completeness.")

                    # 2. Deconvolution (real result)
                    st.markdown("#### Variant deconvolution")
                    if location in st.session_state.location_results:
                        render_single_location_result(
                            location, st.session_state.location_results[location]
                        )
                    else:
                        render_location_progress(
                            location, task_id, celery_app, redis_client
                        )

                    # 3. Scanner (placeholder button)
                    st.markdown("#### Scanner")
                    scan_col1, scan_col2 = st.columns([1, 5])
                    with scan_col1:
                        scan_clicked = st.button(
                            "▶ Run scanner",
                            key=f"acooc_scan_{location}",
                        )
                    with scan_col2:
                        if scan_clicked:
                            st.warning(
                                "Scanner requires the enhanced `/aggregated` endpoint, "
                                "not yet deployed on WASAP LAPIS."
                            )
                        else:
                            st.caption(
                                "Detects mutation combinations not explained by the selected panel."
                            )

        # ── Signature similarity ──────────────────────────────────────────────
        if len(all_selected_variants) >= 2:
            st.markdown("---")
            render_jaccard_heatmap(
                variants=all_selected_variants,
                pango_loader=cached_get_pango_loader(),
            )



if __name__ == "__main__":
    app()