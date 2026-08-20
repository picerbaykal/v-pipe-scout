"""
Abundance & Co-occurrence tab.

Estimates variant abundance over time using LAPIS-sourced wastewater
mutation data, with LolliPop deconvolution (including bootstrap confidence
intervals) and pango-lineage-based sublineage enrichment.

Co-occurrence features (panel-completeness check, emerging-sublineage and
de-novo variant detection) are live — powered by the LAPIS /aggregated
[position] bracket syntax (PR #1768, deployed 2026-08-19 on WASAP 0.8.5).
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
from components.scanner_results import render_scanner_results

from datetime import datetime
import pandas as pd

celery_app = Celery(
    'tasks',
    broker=os.environ.get('CELERY_BROKER_URL', 'redis://redis:6379/0'),
    backend=os.environ.get('CELERY_RESULT_BACKEND', 'redis://redis:6379/0')
)

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
    return PangoLoader(get_pango_summary_path())


@st.cache_data
def cached_get_variant_names() -> list:
    from api.signatures import get_variant_names
    return get_variant_names()


def _render_completeness(result: dict) -> None:
    """Plot per-date panel completeness with read-weighted smoothed line."""
    import numpy as np
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
    }).sort_values("date").reset_index(drop=True)
    df["informative"] = df["matched"] + df["unexplained"]

    WINDOW = 5
    comp = pd.to_numeric(df["completeness"], errors="coerce")
    w = df["informative"].astype(float).where(comp.notna(), 0.0)
    cw = (comp.fillna(0.0) * w)
    num = cw.rolling(WINDOW, center=True, min_periods=1).sum()
    den = w.rolling(WINDOW, center=True, min_periods=1).sum()
    df["smoothed"] = np.where(den > 0, num / den, np.nan)

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df["date"], y=comp, mode="markers",
        marker=dict(size=[6 + min(8, n / 5000) for n in df["informative"]],
                    color="rgba(120,120,120,0.35)"),
        name="per-date",
        customdata=df[["informative"]],
        hovertemplate="%{x|%Y-%m-%d}<br>completeness %{y:.1%}"
                      "<br>%{customdata[0]:,} informative reads<extra></extra>",
    ))
    fig.add_trace(go.Scatter(
        x=df["date"], y=df["smoothed"], mode="lines",
        line=dict(width=2.5, color="#4C6EF5", shape="spline"),
        name="weighted smoothed",
        hovertemplate="%{x|%Y-%m-%d}<br>smoothed %{y:.1%}<extra></extra>",
    ))
    fig.update_yaxes(range=[-0.05, 1.05], tickformat=".0%", title="completeness")
    fig.update_layout(
        height=220, margin=dict(t=10, b=30, l=50, r=20),
        template="plotly_white",
        legend=dict(orientation="h", yanchor="bottom", y=1.0, x=0),
        showlegend=True,
    )
    st.plotly_chart(fig, use_container_width=True)
    st.caption(
        f"Line = read-weighted smoothed (window {WINDOW}). Points = per-date, "
        f"size scales with informative reads "
        f"({df['informative'].min():,}–{df['informative'].max():,} across dates). "
        "Weighting lets high-read dates anchor the curve and down-weights sparse ones."
    )


def _step_label(n: int, label: str, done: bool = False, active: bool = False) -> None:
    """Render a numbered step header in the left column."""
    if done:
        icon = "✓"
        color = "#3B6D11"
        bg = "#EAF3DE"
    elif active:
        icon = str(n)
        color = "#fff"
        bg = "#E24B4A"
    else:
        icon = str(n)
        color = "var(--text-muted)"
        bg = "var(--surface-1)"
    st.markdown(
        f"<div style='display:flex;align-items:center;gap:8px;margin-bottom:4px;'>"
        f"<span style='width:20px;height:20px;border-radius:50%;background:{bg};"
        f"color:{color};font-size:11px;font-weight:500;display:flex;align-items:center;"
        f"justify-content:center;flex:none;border:0.5px solid var(--border);'>{icon}</span>"
        f"<span style='font-size:13px;font-weight:500;"
        f"color:{'var(--text-muted)' if not done and not active else 'var(--text-primary)'};"
        f"'>{label}</span></div>",
        unsafe_allow_html=True,
    )


def app():
    # Apply any pending variant additions from the scanner BEFORE widgets render
    for key in list(st.session_state.keys()):
        if key.startswith("acooc_add_variant_pending_"):
            v = st.session_state.pop(key)
            current = st.session_state.get("acooc_variant_multiselect", [])
            if v not in current:
                st.session_state["acooc_variant_multiselect"] = current + [v]
            break

    st.session_state.setdefault("location_results", {})
    st.session_state.setdefault("acooc_location_tasks", {})
    st.session_state.setdefault("acooc_cooc_tasks", {})
    st.session_state.setdefault("acooc_cooc_results", {})
    st.session_state.setdefault("acooc_scanner_results", {})
    st.session_state.setdefault("acooc_scanner_panels", {})

    # ── Header ───────────────────────────────────────────────────────────────
    st.title("Abundance & Co-occurrence")
    st.subheader(
        "Estimate the proportion of variants circulating in wastewater over time, "
        "using live LAPIS mutation data and LolliPop deconvolution with bootstrap "
        "confidence intervals."
    )
    st.caption(
        "Build a custom variant panel and run on-demand deconvolution across one or "
        "more sampling locations. The panel-completeness check and scanner guide you "
        "to a complete panel iteratively — follow the numbered steps."
    )
    st.markdown("---")

    # ── Two-column layout ─────────────────────────────────────────────────────
    col_controls, col_results = st.columns([1, 3])

    # derive state flags for step indicators
    has_variants = len(st.session_state.get("acooc_variant_multiselect", [])) >= 2
    has_locations = len(st.session_state.get("acooc_location_multiselect", [])) >= 1
    has_run = bool(st.session_state.get("acooc_location_tasks"))
    has_completeness = bool(st.session_state.get("acooc_cooc_results"))
    has_scanner = bool(st.session_state.get("acooc_scanner_results"))

    with col_controls:
        curated_variants = cached_get_variant_names()

        # ── Step 1: Variant panel ─────────────────────────────────────────────
        _step_label(1, "Variant panel", done=has_variants)

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
            label_visibility="collapsed",
        )

        with st.expander("Add lineage manually", expanded=False):
            pango_loader = cached_get_pango_loader()
            available_lineages = sorted(pango_loader.get_raw_data().keys())
            extra_options = [v for v in available_lineages if v not in selected_variants]
            extra_variants = st.multiselect(
                "Search pango lineage",
                options=extra_options,
                placeholder="e.g. KP.2.3",
                help="Add any pango lineage. Use the scanner (step 6) to discover missing variants automatically.",
                key="acooc_extra_variants",
            )

        all_selected_variants = selected_variants + extra_variants
        if len(all_selected_variants) >= 2:
            st.caption(f"{len(all_selected_variants)} variants in panel")
        else:
            st.caption("Select at least 2 variants.")

        st.markdown("---")

        # ── Step 2: Variant tree ──────────────────────────────────────────────
        _step_label(2, "Variant tree", done=has_variants)
        render_panel_tree(
            selected_variants=all_selected_variants,
            yaml_variants=curated_variants,
            pango_loader=cached_get_pango_loader(),
            scanner_added_for=st.session_state.get("acooc_scanner_added_for", {}),
        )

        st.markdown("---")

        # ── Step 3: Locations ─────────────────────────────────────────────────
        _step_label(3, "Locations", done=has_locations)
        available_locations = wiseLoculus.fetch_locations()
        col_loc_all, col_loc_clear = st.columns(2)
        with col_loc_all:
            if st.button("Select all", key="acooc_loc_select_all", use_container_width=True):
                st.session_state["acooc_location_multiselect"] = available_locations
        with col_loc_clear:
            if st.button("Clear", key="acooc_loc_clear", use_container_width=True):
                st.session_state["acooc_location_multiselect"] = []
        selected_locations = st.multiselect(
            "Select sampling locations",
            options=available_locations,
            placeholder="Select one or more locations...",
            key="acooc_location_multiselect",
            label_visibility="collapsed",
        )

        st.markdown("---")

        # ── Step 3b: Date range + LolliPop params ─────────────────────────────
        _step_label(3, "Date range", done=has_locations)
        default_start, default_end, min_date, max_date = wiseLoculus.get_cached_date_range_with_bounds("abundance_cooc")
        col_start, col_end = st.columns(2)
        with col_start:
            start_date = st.date_input(
                "Start", value=default_start, min_value=min_date,
                max_value=max_date, key="acooc_start_date",
            )
        with col_end:
            end_date = st.date_input(
                "End", value=default_end, min_value=min_date,
                max_value=max_date, key="acooc_end_date",
            )
        if end_date <= start_date:
            st.warning("End date must be after start date.")

        with st.expander("LolliPop parameters", expanded=False):
            bootstrap_options = {"Rapid": 50, "Standard": 100, "Reliable": 300}
            selected_bootstrap = st.radio(
                "Bootstrap iterations", options=list(bootstrap_options.keys()),
                index=1, key="acooc_bootstrap",
            )
            bootstraps = bootstrap_options[selected_bootstrap]
            bandwidth_options = {"Narrow": 10, "Medium": 20, "Wide": 30}
            selected_bandwidth = st.radio(
                "Bandwidth", options=list(bandwidth_options.keys()),
                index=0, key="acooc_bandwidth",
            )
            bandwidth = bandwidth_options[selected_bandwidth]

        st.markdown("---")

        # ── Step 4: Run analysis ──────────────────────────────────────────────
        can_run = (
            len(all_selected_variants) >= 2
            and len(selected_locations) >= 1
            and end_date > start_date
        )
        _step_label(4, "Run analysis", done=has_run, active=not has_run and can_run)
        if st.button(
            "▶ Run analysis",
            type="primary",
            disabled=not can_run,
            key="acooc_run_button",
            use_container_width=True,
            help="Runs deconvolution + panel-completeness check for all selected locations.",
        ):
            st.session_state["acooc_trigger_run"] = True


        # scanner state — read by the right column city tabs
        scanner_results = st.session_state.get("acooc_scanner_results", {})
        scanner_panels = st.session_state.get("acooc_scanner_panels", {})

    # ── Right column: all outputs ─────────────────────────────────────────────
    with col_results:

        # task submission
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
                logger.info(f"Submitted task {cooc_task.id} for {loc}")

            st.session_state["acooc_location_tasks"] = location_tasks
            st.session_state["location_results"] = {}
            st.session_state["acooc_cooc_tasks"] = cooc_tasks
            existing_cooc = st.session_state.get("acooc_cooc_results", {})
            for loc in cooc_tasks:
                existing_cooc.pop(loc, None)
            st.session_state["acooc_cooc_results"] = existing_cooc
            st.rerun()

        from components.multi_location_results import (
            render_single_location_result,
            render_location_progress,
        )

        location_tasks = st.session_state.get("acooc_location_tasks", {})

        if not location_tasks:
            st.info("Complete steps 1–5 on the left to see results here.")
        else:
            cooc_task_ids = list(st.session_state.get("acooc_cooc_tasks", {}).values())
            scanner_task_ids = list(st.session_state.get("acooc_scanner_tasks", {}).values())
            all_task_ids = list(location_tasks.values()) + cooc_task_ids + scanner_task_ids
            states = {tid: celery_app.AsyncResult(tid).state for tid in all_task_ids}
            n_results = (
                len(st.session_state.location_results)
                + len(st.session_state.get("acooc_cooc_results", {}))
            )
            deconv_or_cooc_running = any(
                celery_app.AsyncResult(tid).state in ("PENDING", "STARTED", "RETRY")
                for tid in all_task_ids
            )
            any_just_completed = (
                any(s == "SUCCESS" for s in states.values())
                and n_results < len(all_task_ids)
            )
            if deconv_or_cooc_running or any_just_completed:
                from streamlit_autorefresh import st_autorefresh
                st_autorefresh(interval=2000, key="acooc_autorefresh")

            # collect completed scanner task results
            scanner_tasks_map = st.session_state.get("acooc_scanner_tasks", {})
            for _loc, _tid in list(scanner_tasks_map.items()):
                if _loc not in scanner_results:
                    _t = celery_app.AsyncResult(_tid)
                    if _t.ready():
                        try:
                            scanner_results[_loc] = _t.get()
                            st.session_state["acooc_scanner_results"] = scanner_results
                        except Exception as _e:
                            logger.error(f"Scanner task failed for {_loc}: {_e}")

            location_names = list(location_tasks.keys())

            tab_labels = ["📊 Summary"] + [f"📍 {loc}" for loc in location_names]
            all_tabs = st.tabs(tab_labels)
            summary_tab = all_tabs[0]
            city_tabs = all_tabs[1:]

            added_for = st.session_state.get("acooc_scanner_added_for", {})

            for location, tab in zip(location_names, city_tabs):
                with tab:
                    task_id = location_tasks[location]
                    cooc_tasks_map = st.session_state.get("acooc_cooc_tasks", {})
                    cooc_results = st.session_state.get("acooc_cooc_results", {})

                    # ── completeness ──────────────────────────────────────────
                    st.markdown("#### Panel completeness (co-occurrence)")
                    if location in cooc_results:
                        _render_completeness(cooc_results[location])
                    elif location in cooc_tasks_map:
                        cooc_task = celery_app.AsyncResult(cooc_tasks_map[location])
                        if cooc_task.ready():
                            try:
                                cooc_results[location] = cooc_task.get()
                                st.session_state["acooc_cooc_results"] = cooc_results
                                _render_completeness(cooc_results[location])
                            except Exception as e:
                                st.error(f"Co-occurrence failed: {e}")
                        else:
                            st.info("Computing panel completeness…")
                    else:
                        st.caption("Run analysis to compute panel completeness.")

                    # ── Jaccard collapsed ─────────────────────────────────────
                    if len(all_selected_variants) >= 2:
                        with st.expander("Signature similarity (Jaccard)", expanded=False):
                            render_jaccard_heatmap(
                                variants=all_selected_variants,
                                pango_loader=cached_get_pango_loader(),
                            )

                    # ── deconvolution ─────────────────────────────────────────
                    st.markdown("#### Variant deconvolution")
                    if location in st.session_state.location_results:
                        render_single_location_result(
                            location, st.session_state.location_results[location]
                        )
                    else:
                        render_location_progress(
                            location, task_id, celery_app, redis_client
                        )

                    # ── scanner ───────────────────────────────────────────────
                    st.markdown("---")
                    _scan_col1, _scan_col2 = st.columns([3, 1])
                    with _scan_col1:
                        _scan_badge = ""
                        if location in cooc_results:
                            _matched = sum(cooc_results[location].get("matched_counts", []))
                            _unexplained = sum(cooc_results[location].get("unexplained_counts", []))
                            _total = _matched + _unexplained
                            if _total > 0:
                                _pct = _matched / _total
                                _scan_badge = f"{'✓ complete' if _pct >= 0.95 else '⚠ panel incomplete'}"
                        st.markdown(f"#### Scanner {_scan_badge}")
                    with _scan_col2:
                        _can_scan = location in cooc_results and bool(
                            cooc_results[location].get("unexplained_patterns")
                        )
                        if _can_scan:
                            _stale = set(scanner_panels.get(location, [])) != set(all_selected_variants)
                            _btn_label = "▶ Re-scan" if location in scanner_results and not _stale else "▶ Run scanner"
                            if st.button(_btn_label, key=f"acooc_scan_{location}", use_container_width=True):
                                _scan_task = celery_app.send_task(
                                    "tasks.run_cooc_scanner_lapis",
                                    kwargs={
                                        "location": location,
                                        "start_date": start_date.isoformat(),
                                        "end_date": end_date.isoformat(),
                                        "variants": all_selected_variants,
                                        "unexplained_patterns": cooc_results[location]["unexplained_patterns"],
                                    }
                                )
                                _sct = st.session_state.get("acooc_scanner_tasks", {})
                                _sct[location] = _scan_task.id
                                scanner_panels[location] = list(all_selected_variants)
                                st.session_state["acooc_scanner_tasks"] = _sct
                                st.session_state["acooc_scanner_panels"] = scanner_panels
                                scanner_results.pop(location, None)
                                st.session_state["acooc_scanner_results"] = scanner_results
                                # st.rerun() - autorefresh picks up the new task
                        else:
                            st.caption("Available after completeness.")

                    if location in scanner_results:
                        _stale = set(scanner_panels.get(location, [])) != set(all_selected_variants)
                        if _stale:
                            st.warning("Panel changed since last scan — re-run scanner.")

                        _loc_ref = location  # capture for closure
                        def _add_variant(v, _loc=_loc_ref):
                            st.session_state[f"acooc_add_variant_pending_{_loc}"] = v
                            _af = st.session_state.get("acooc_scanner_added_for", {})
                            _af[v] = _loc
                            st.session_state["acooc_scanner_added_for"] = _af
                            # no st.rerun() — Streamlit's natural button rerun handles it

                        render_scanner_results(
                            scanner_results[location],
                            selected_variants=all_selected_variants,
                            client=wiseLoculus,
                            location=location,
                            date_range=(
                                datetime.combine(start_date, datetime.min.time()),
                                datetime.combine(end_date, datetime.min.time()),
                            ),
                            on_add_variant=_add_variant,
                        )

                        # ── re-run nudge ──────────────────────────────────────
                        _loc_added = [
                            v for v, l in added_for.items()
                            if l == location and v in all_selected_variants
                        ]
                        if _loc_added:
                            _short = location.split("(")[0].strip()
                            with st.container():
                                st.markdown(
                                    f"<div style='background:#EFF6FF;border:0.5px solid #BFDBFE;"
                                    f"border-radius:8px;padding:8px 12px;margin-top:8px;'>"
                                    f"<span style='font-size:12px;font-weight:500;color:#1E40AF;'>"
                                    f"{', '.join(f'+{v}' for v in _loc_added)} added to panel</span>"
                                    f"</div>",
                                    unsafe_allow_html=True,
                                )
                                if st.button(
                                    f"▶ Re-run {_short}",
                                    key=f"acooc_rerun_{location}",
                                    type="primary",
                                    use_container_width=True,
                                ):
                                    _rt = celery_app.send_task(
                                        "tasks.run_deconvolve_lapis",
                                        kwargs={
                                            "location": location,
                                            "start_date": start_date.isoformat(),
                                            "end_date": end_date.isoformat(),
                                            "variants": all_selected_variants,
                                            "bootstraps": bootstraps,
                                            "bandwidth": bandwidth,
                                        }
                                    )
                                    _rct = celery_app.send_task(
                                        "tasks.run_cooc_completeness_lapis",
                                        kwargs={
                                            "location": location,
                                            "start_date": start_date.isoformat(),
                                            "end_date": end_date.isoformat(),
                                            "variants": all_selected_variants,
                                        }
                                    )
                                    st.session_state["acooc_location_tasks"][location] = _rt.id
                                    st.session_state["acooc_cooc_tasks"][location] = _rct.id
                                    st.session_state["location_results"].pop(location, None)
                                    st.session_state["acooc_cooc_results"].pop(location, None)
                                    st.session_state["acooc_scanner_results"].pop(location, None)
                                    st.session_state["acooc_scanner_added_for"] = {
                                        v: l for v, l in added_for.items() if l != location
                                    }
                                    logger.info(f"Re-run submitted for {location} only")
                                    st.session_state["acooc_force_rerun"] = True
                    elif location in scanner_tasks_map:
                        _st = celery_app.AsyncResult(scanner_tasks_map[location])
                        if _st.state in ("PENDING", "STARTED", "RETRY"):
                            with st.status("Scanner running…", expanded=True) as _scan_status:
                                st.write(f"Classifying unexplained patterns for {location}…")
                                st.write("Checking against known lineages and panel variants…")
                            _scan_status.update(label="Scanner running…", state="running")
                        elif _st.state == "FAILURE":
                            st.error(f"Scanner failed — check worker logs.")
                    else:
                        st.caption("Run scanner to classify unexplained patterns.")

            # ── Summary tab (last) ─────────────────────────────────────────────
            with summary_tab:
                st.markdown("#### Location status")
                cooc_results = st.session_state.get("acooc_cooc_results", {})

                # status table
                _cols = st.columns([2, 2, 3, 1])
                _cols[0].caption("Location")
                _cols[1].caption("Completeness")
                _cols[2].caption("Status")
                _cols[3].caption("")
                st.markdown("<hr style='margin:4px 0;opacity:.3;'>", unsafe_allow_html=True)

                for loc in location_names:
                    _c1, _c2, _c3, _c4 = st.columns([2, 2, 3, 1])
                    _c1.write(f"📍 {loc.split('(')[0].strip()}")
                    if loc in cooc_results:
                        _m = sum(cooc_results[loc].get("matched_counts", []))
                        _u = sum(cooc_results[loc].get("unexplained_counts", []))
                        _t = _m + _u
                        _pct = int(_m / _t * 100) if _t > 0 else 0
                        _c2.write(f"{_pct}%")
                        if _pct >= 95:
                            _c3.success("✓ complete")
                        elif loc in scanner_results:
                            _n_miss = len(scanner_results[loc].get("missing_from_panel", []))
                            _c3.error(f"⚠ {_n_miss} missing")
                        else:
                            _c3.warning("⚠ run scanner")
                        _c4.write("")
                    else:
                        _c2.caption("—")
                        _c3.caption("not run")
                        _c4.write("")

                st.markdown("---")

                # multi-city report
                _rep_col1, _rep_col2 = st.columns([3, 1])
                with _rep_col1:
                    st.markdown("#### Multi-city report")
                    st.caption("All deconvolution curves in one view · download-ready")
                with _rep_col2:
                    if st.button(
                        "Generate report",
                        key="acooc_gen_report",
                        use_container_width=True,
                    ):
                        st.session_state["acooc_show_report"] = True

                if st.session_state.get("acooc_show_report"):
                    _completed = [loc for loc in location_names if loc in st.session_state.location_results]
                    if not _completed:
                        st.info("Run analysis for at least one location to generate the report.")
                    else:
                        import plotly.graph_objects as go
                        from plotly.subplots import make_subplots
                        import math
                        _ncols = 2
                        _nrows = math.ceil(len(_completed) / _ncols)
                        _fig = make_subplots(
                            rows=_nrows, cols=_ncols,
                            subplot_titles=_completed,
                            vertical_spacing=0.12,
                            horizontal_spacing=0.08,
                        )
                        for _i, _loc in enumerate(_completed):
                            _row = _i // _ncols + 1
                            _col = _i % _ncols + 1
                            _res = st.session_state.location_results[_loc]
                            for _variant, _data in _res.items():
                                if _variant == "undetermined":
                                    continue
                                _ts = _data.get("timeseriesSummary", [])
                                if not _ts:
                                    continue
                                _dates = [e.get("date") for e in _ts]
                                _props = [e.get("proportion", 0) for e in _ts]
                                _fig.add_trace(
                                    go.Scatter(x=_dates, y=_props, name=_variant,
                                               mode="lines+markers",
                                               marker=dict(size=4),
                                               showlegend=(_i == 0)),
                                    row=_row, col=_col,
                                )
                        _fig.update_layout(
                            height=300 * _nrows,
                            template="plotly_white",
                            title_text="Variant Proportion Estimates — All Locations",
                        )
                        st.plotly_chart(_fig, use_container_width=True)
                        # download button
                        import io
                        try:
                            _buf = io.BytesIO()
                            _fig.write_image(_buf, format="pdf")
                            st.download_button(
                                "⬇ Download PDF",
                                data=_buf.getvalue(),
                                file_name="vpipe_scout_report.pdf",
                                mime="application/pdf",
                                key="acooc_download_pdf",
                            )
                        except Exception:
                            st.caption("Install kaleido for PDF export: pip install kaleido")


if __name__ == "__main__":
    app()