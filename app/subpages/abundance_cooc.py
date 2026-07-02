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


def app():
    st.session_state.setdefault("location_results", {})
    st.session_state.setdefault("acooc_location_tasks", {})

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

        selected_variants = st.multiselect(
            "Surveillance variants",
            options=curated_variants,
            placeholder="Start typing to search for a lineage...",
            help="Curated variants from the Swiss wastewater surveillance panel (cowwid).",
            key="accoc_variant_multiselect"
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
                logger.info(f"Submitted task {task.id} for {loc}")

            st.session_state["acooc_location_tasks"] = location_tasks
            st.session_state["location_results"] = {}
            st.rerun()

        # ── Results ───────────────────────────────────────────────────────────
        from components.multi_location_results import render_location_results_tabs

        location_tasks = st.session_state.get("acooc_location_tasks", {})
        location_results = st.session_state.get("location_results", {})

        if not location_tasks:
            st.info("Results will appear here after running deconvolution.")
        else:
            location_tasks = st.session_state.get("acooc_location_tasks", {})

            # Check task states
            states = {
                task_id: celery_app.AsyncResult(task_id).state
                for task_id in location_tasks.values()
            }

            any_running = any(
                s in ("PENDING", "STARTED", "RETRY")
                for s in states.values()
            )

            any_just_completed = any(
                s == "SUCCESS"
                for s in states.values()
            ) and len(st.session_state.location_results) < len(location_tasks)

            if any_running or any_just_completed:
                from streamlit_autorefresh import st_autorefresh
                st_autorefresh(interval=2000, key="acooc_autorefresh")

            render_location_results_tabs(
                location_tasks=location_tasks,
                location_results=st.session_state.location_results,
                celery_app=celery_app,
                redis_client=redis_client,
            )


if __name__ == "__main__":
    app()