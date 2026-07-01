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

from api.wiseloculus import WiseLoculusLapis
from utils.config import get_wiseloculus_url
from datetime import date

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


def app():
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

    # ── Variant Panel ────────────────────────────────────────────────────────
    st.subheader("Variant Panel")
    pango_loader = cached_get_pango_loader()
    available_lineages = sorted(pango_loader.get_raw_data().keys())

    selected_variants = st.multiselect(
        "Select variants of interest",
        options=available_lineages,
        placeholder="Start typing to search for a lineage...",
        help="Select any pango lineage. Requires at least 2 variants to run deconvolution.",
        key="accoc_variant_multiselect"
    )

    if len(selected_variants) < 2:
        st.warning("Select at least 2 variants to run deconvolution.")

    st.markdown("---")

    # ── Location & Date Range ─────────────────────────────────────────────
    st.subheader("Sampling Locations & Date Range")

    available_locations = wiseLoculus.fetch_locations()

    selected_locations = st.multiselect(
        "Select sampling locations",
        options=available_locations,
        default=st.session_state.get("acooc_selected_locations", []),
        placeholder="Select one or more locations...",
        key="acooc_location_multiselect",
    )

    if not selected_locations:
        st.warning("Select at least one location.")

    default_start, default_end, min_date, max_date = wiseLoculus.get_cached_date_range_with_bounds("abundance_cooc")

    default_start, default_end, min_date, max_date = wiseLoculus.get_cached_date_range_with_bounds("abundance_cooc")

    col1, col2 = st.columns(2)

    with col1:
        start_date = st.date_input(
            "Start Date",
            value=default_start,
            min_value=min_date,
            max_value=max_date,
            key="acooc_start_date",
        )
    with col2:
        end_date = st.date_input(
            "End Date",
            value=default_end,
            min_value=min_date,
            max_value=max_date,
            key="acooc_end_date",
        )

    if end_date <= start_date:
        st.warning("End date must be after start date.")
    else:
        date_range = (start_date, end_date)

        if len(date_range) != 2:
            st.warning("Please select a start and end date.")



if __name__ == "__main__":
    app()