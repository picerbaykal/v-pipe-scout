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
            st.rerun()
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
            help="Runs deconvolution, completeness check, and scanner for all locations.",
        ):
            st.session_state["acooc_trigger_run"] = True


        # scanner state — read by result collection in right column
        scanner_results = st.session_state.get("acooc_scanner_results", {})
        scanner_panels = st.session_state.get("acooc_scanner_panels", {})

        # scanner status in left column (auto-runs, no button needed)
        _stasks_left = st.session_state.get("acooc_scanner_tasks", {})
        if _stasks_left:
            st.markdown("---")
            _step_label(6, "Scanner", done=bool(st.session_state.get("acooc_scanner_results")), active=False)
            for _sloc, _stid in _stasks_left.items():
                _sstate = celery_app.AsyncResult(_stid).state
                if _sstate in ("PENDING", "STARTED", "RETRY"):
                    st.caption(f"⟳ Scanning {_sloc.split('(')[0].strip()}…")
                elif _sstate == "SUCCESS":
                    _sr = st.session_state.get("acooc_scanner_results", {})
                    if _sloc in _sr:
                        _nm = len(_sr[_sloc].get("missing_from_panel", []))
                        st.caption(f"✓ {_sloc.split('(')[0].strip()} — {_nm} missing" if _nm else f"✓ {_sloc.split('(')[0].strip()} — panel ok")
                elif _sstate == "FAILURE":
                    st.caption(f"✗ {_sloc.split('(')[0].strip()} failed")
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

            # autorefresh only while tasks are running — stops when all done
            # using 3s interval to give render cycles enough time to complete
            _all_tids = list(location_tasks.values()) + cooc_task_ids + scanner_task_ids
            _running_tids = [
                tid for tid in _all_tids
                if tid and celery_app.AsyncResult(tid).state in ("PENDING", "STARTED", "RETRY")
            ]
            if _running_tids:
                from streamlit_autorefresh import st_autorefresh
                st_autorefresh(interval=3000, key="acooc_autorefresh")

            # collect completed results
            scanner_tasks_map = st.session_state.get("acooc_scanner_tasks", {})
            for _loc, _tid in list(scanner_tasks_map.items()):
                if _loc not in scanner_results:
                    _t = celery_app.AsyncResult(_tid)
                    if _t.ready():
                        try:
                            scanner_results[_loc] = _t.get()
                            st.session_state["acooc_scanner_results"] = scanner_results
                            logger.info(f"Scanner results collected for {_loc}")
                        except Exception as _e:
                            logger.error(f"Scanner task failed for {_loc}: {_e}")

            # collect completed cooc results + auto-submit scanner
            _cooc_res = st.session_state.get("acooc_cooc_results", {})
            _scanner_tasks = st.session_state.get("acooc_scanner_tasks", {})
            for _loc, _tid in list(st.session_state.get("acooc_cooc_tasks", {}).items()):
                if _loc not in _cooc_res:
                    _t = celery_app.AsyncResult(_tid)
                    if _t.ready():
                        try:
                            _cooc_res[_loc] = _t.get()
                            st.session_state["acooc_cooc_results"] = _cooc_res
                        except Exception:
                            pass
                # auto-submit scanner when completeness is ready and scanner not yet run
                if (_loc in _cooc_res
                        and _cooc_res[_loc].get("unexplained_patterns")
                        and _loc not in _scanner_tasks
                        and _loc not in scanner_results):
                    _stask = celery_app.send_task(
                        "tasks.run_cooc_scanner_lapis",
                        kwargs={
                            "location": _loc,
                            "start_date": start_date.isoformat(),
                            "end_date": end_date.isoformat(),
                            "variants": all_selected_variants,
                            "unexplained_patterns": _cooc_res[_loc]["unexplained_patterns"],
                        }
                    )
                    _scanner_tasks[_loc] = _stask.id
                    scanner_panels[_loc] = list(all_selected_variants)
                    logger.info(f"Auto-submitted scanner for {_loc}")
            st.session_state["acooc_scanner_tasks"] = _scanner_tasks
            st.session_state["acooc_scanner_panels"] = scanner_panels

            location_names = list(location_tasks.keys())

            # ── Progress header ───────────────────────────────────────────────
            _cooc_res2 = st.session_state.get("acooc_cooc_results", {})
            _scan_res2 = st.session_state.get("acooc_scanner_results", {})
            _scan_tasks2 = st.session_state.get("acooc_scanner_tasks", {})

            def _city_status(loc):
                """Return (pct, status_label, step_states) for a city."""
                _deconv_done = loc in st.session_state.get("location_results", {})
                _cooc_done = loc in _cooc_res2
                _scan_done = loc in _scan_res2
                _scan_running = (loc in _scan_tasks2 and
                    celery_app.AsyncResult(_scan_tasks2[loc]).state in ("PENDING","STARTED","RETRY"))
                _deconv_running = (loc in location_tasks and
                    celery_app.AsyncResult(location_tasks[loc]).state in ("PENDING","STARTED","RETRY"))
                _cooc_running = (loc in st.session_state.get("acooc_cooc_tasks",{}) and
                    celery_app.AsyncResult(st.session_state["acooc_cooc_tasks"][loc]).state in ("PENDING","STARTED","RETRY"))
                _pct = None
                if _cooc_done:
                    _m = sum(_cooc_res2[loc].get("matched_counts",[]))
                    _u = sum(_cooc_res2[loc].get("unexplained_counts",[]))
                    _t = _m + _u
                    _pct = int(_m/_t*100) if _t>0 else 0
                _n_miss = len(_scan_res2.get(loc,{}).get("missing_from_panel",[])) if _scan_done else 0
                return _pct, _n_miss, _scan_done, _scan_running, _deconv_done, _cooc_done, _deconv_running, _cooc_running

            def _tid_running(tid):
                return bool(tid) and celery_app.AsyncResult(tid).state in ("PENDING","STARTED","RETRY")

            _n_complete = sum(1 for loc in location_names
                if loc in st.session_state.get("location_results", {})
                and loc in _cooc_res2
                and loc in _scan_res2)
            _n_running = sum(1 for loc in location_names
                if _tid_running(location_tasks.get(loc,""))
                or _tid_running(st.session_state.get("acooc_cooc_tasks",{}).get(loc,""))
                or _tid_running(st.session_state.get("acooc_scanner_tasks",{}).get(loc,"")))
            _n_pending = len(location_names) - _n_complete - _n_running

            with st.container():
                _ph1, _ph2 = st.columns([3,1])
                with _ph1:
                    st.markdown(
                        f"<div style='display:flex;gap:16px;align-items:center;padding:4px 0;'>"
                        f"<span style='font-size:13px;font-weight:500;'>Analysis progress</span>"
                        f"<span style='font-size:12px;color:#3B6D11;font-weight:500;'>✓ {_n_complete} complete</span>"
                        f"<span style='font-size:12px;color:#185FA5;'>⟳ {_n_running} running</span>"
                        f"<span style='font-size:12px;color:#898781;'>○ {_n_pending} pending</span>"
                        f"</div>",
                        unsafe_allow_html=True,
                    )
                with _ph2:
                    # download report button
                    _completed_locs = [loc for loc in location_names if loc in st.session_state.get("location_results",{})]
                    if _completed_locs:
                        if st.button("⬇ Download report", key="acooc_dl_report", use_container_width=True):
                            st.session_state["acooc_show_report"] = True

            # ── Compact progress rows ─────────────────────────────────────────
            st.markdown(
                "<div style='display:flex;gap:10px;margin-bottom:4px;font-size:10px;color:#898781;'>"
                "<span style='display:flex;align-items:center;gap:3px;'>"
                "<span style='width:12px;height:3px;border-radius:2px;background:#3B6D11;display:inline-block;'></span>done</span>"
                "<span style='display:flex;align-items:center;gap:3px;'>"
                "<span style='width:12px;height:3px;border-radius:2px;background:#93C5FD;display:inline-block;'></span>running</span>"
                "<span style='display:flex;align-items:center;gap:3px;'>"
                "<span style='width:12px;height:3px;border-radius:2px;background:#E5E3DC;display:inline-block;'></span>pending</span>"
                "<span style='margin-left:6px;'>deconv · completeness · scanner</span></div>",
                unsafe_allow_html=True,
            )
            def _dot_color(done, running):
                return "#3B6D11" if done else ("#93C5FD" if running else "#E5E3DC")

            for _loc in location_names:
                _pct, _n_miss, _scan_done, _scan_running, _deconv_done, _cooc_done, _deconv_running, _cooc_running = _city_status(_loc)
                _c1 = _dot_color(_deconv_done, _deconv_running)
                _c2 = _dot_color(_cooc_done, _cooc_running)
                _c3 = _dot_color(_scan_done, _scan_running)
                _pct_str = f"{_pct}%" if _pct is not None else "—"
                if _scan_done:
                    _badge = f"⚠ {_n_miss} missing" if _n_miss > 0 else ("✓ complete" if _pct is not None and _pct >= 95 else "✓ scanned")
                    _badge_bg, _badge_col = ("#FCEBEB","#A32D2D") if _n_miss > 0 else ("#EAF3DE","#3B6D11")
                elif _scan_running:
                    _badge, _badge_bg, _badge_col = "scanning…","#EFF6FF","#1E40AF"
                elif _cooc_running or _deconv_running:
                    _badge, _badge_bg, _badge_col = "running…","#EFF6FF","#1E40AF"
                else:
                    _badge, _badge_bg, _badge_col = "queued","#F1EFE8","#898781"
                st.markdown(
                    f"<div style='display:flex;align-items:center;gap:8px;padding:4px 0;"
                    f"border-bottom:0.5px solid rgba(0,0,0,.06);'>"
                    f"<span style='font-size:11px;font-weight:500;width:90px;flex:none;"
                    f"overflow:hidden;text-overflow:ellipsis;white-space:nowrap;'>{_loc.split('(')[0].strip()}</span>"
                    f"<div style='display:flex;gap:2px;flex:none;'>"
                    f"<div style='width:16px;height:4px;border-radius:2px;background:{_c1};'></div>"
                    f"<div style='width:16px;height:4px;border-radius:2px;background:{_c2};'></div>"
                    f"<div style='width:16px;height:4px;border-radius:2px;background:{_c3};'></div>"
                    f"</div>"
                    f"<span style='font-size:10px;color:#898781;width:28px;flex:none;text-align:right;'>{_pct_str}</span>"
                    f"<span style='font-size:10px;padding:1px 5px;border-radius:5px;"
                    f"background:{_badge_bg};color:{_badge_col};flex:none;'>{_badge}</span>"
                    f"</div>",
                    unsafe_allow_html=True,
                )

            st.markdown("<hr style='margin:8px 0 6px;opacity:.15;'>", unsafe_allow_html=True)

            # ── Tab strip ────────────────────────────────────────────────────
            _city_options = [f"📍 {loc}" for loc in location_names]
            if ("acooc_selected_city" not in st.session_state or
                    st.session_state.get("acooc_selected_city") not in _city_options):
                st.session_state["acooc_selected_city"] = _city_options[0] if _city_options else ""

            # HTML visual tab strip
            _tab_html = (
                "<div style='display:flex;gap:0;background:#fff;border:0.5px solid rgba(0,0,0,.08);"
                "border-radius:10px;overflow:hidden;margin-bottom:4px;'>"
            )
            for _ti, _loc in enumerate(location_names):
                _p2, _nm2, _sd2, _sr2, _dd2, _cd2, _dr2, _cr2 = _city_status(_loc)
                _dc1 = _dot_color(_dd2, _dr2)
                _dc2 = _dot_color(_cd2, _cr2)
                _dc3 = _dot_color(_sd2, _sr2)
                _is_on = st.session_state.get("acooc_selected_city","") == f"📍 {_loc}"
                _ts = "border-bottom:2px solid #E24B4A;font-weight:500;background:#faf9f5;" if _is_on else "border-bottom:2px solid transparent;"
                _br = "" if _ti == len(location_names)-1 else "border-right:0.5px solid rgba(0,0,0,.06);"
                _sn = _loc.split("(")[0].strip()
                _tc = "#0b0b0b" if _is_on else "#898781"
                _tab_html += (
                    f"<div style='padding:6px 10px;font-size:11px;flex:1;text-align:center;color:{_tc};{_ts}{_br}'>"
                    f"<div style='display:flex;gap:2px;justify-content:center;margin-bottom:2px;'>"
                    f"<div style='width:10px;height:3px;border-radius:2px;background:{_dc1};'></div>"
                    f"<div style='width:10px;height:3px;border-radius:2px;background:{_dc2};'></div>"
                    f"<div style='width:10px;height:3px;border-radius:2px;background:{_dc3};'></div>"
                    f"</div>{_sn}</div>"
                )
            _tab_html += "</div>"
            st.markdown(_tab_html, unsafe_allow_html=True)

            # invisible radio drives actual selection (hidden via CSS)
            _selected = st.radio(
                "City", options=_city_options, horizontal=True,
                label_visibility="collapsed", key="acooc_selected_city",
            )
            st.markdown(
                "<style>div[data-testid='stRadio']{display:none!important;}</style>",
                unsafe_allow_html=True,
            )
            st.markdown("<hr style='margin:2px 0 10px;opacity:.1;'>", unsafe_allow_html=True)

            def _city_tab_content(location, task_id):
                """Shared content for both active and idle city tab fragments."""
                _cooc_tasks_map = st.session_state.get("acooc_cooc_tasks", {})
                _cooc_results = st.session_state.get("acooc_cooc_results", {})
                _scanner_results = st.session_state.get("acooc_scanner_results", {})
                _scanner_panels = st.session_state.get("acooc_scanner_panels", {})
                _added_for = st.session_state.get("acooc_scanner_added_for", {})
                _scanner_tasks_map = st.session_state.get("acooc_scanner_tasks", {})

                # ── completeness ──────────────────────────────────────────────
                st.markdown("#### Panel completeness (co-occurrence)")
                if location in _cooc_results:
                    _render_completeness(_cooc_results[location])
                elif location in _cooc_tasks_map:
                    _cooc_task = celery_app.AsyncResult(_cooc_tasks_map[location])
                    if _cooc_task.ready():
                        try:
                            _cooc_results[location] = _cooc_task.get()
                            st.session_state["acooc_cooc_results"] = _cooc_results
                            _render_completeness(_cooc_results[location])
                        except Exception as _e:
                            st.error(f"Co-occurrence failed: {_e}")
                    else:
                        st.info("Computing panel completeness…")
                else:
                    st.caption("Run analysis to compute panel completeness.")

                # ── Jaccard collapsed ─────────────────────────────────────────
                if len(all_selected_variants) >= 2:
                    with st.expander("Signature similarity (Jaccard)", expanded=False):
                        render_jaccard_heatmap(
                            variants=all_selected_variants,
                            pango_loader=cached_get_pango_loader(),
                        )

                # ── deconvolution ─────────────────────────────────────────────
                st.markdown("#### Variant deconvolution")
                if location in st.session_state.location_results:
                    render_single_location_result(
                        location, st.session_state.location_results[location]
                    )
                else:
                    render_location_progress(
                        location, task_id, celery_app, redis_client
                    )

                # ── scanner (auto-runs after completeness) ────────────────────
                st.markdown("---")
                _scan_badge = ""
                if location in _cooc_results:
                    _mm = sum(_cooc_results[location].get("matched_counts", []))
                    _uu = sum(_cooc_results[location].get("unexplained_counts", []))
                    _tt = _mm + _uu
                    if _tt > 0:
                        _pp = _mm / _tt
                        _scan_badge = "✓ complete" if _pp >= 0.95 else "⚠ panel incomplete"
                _scan_col_h, _scan_col_note = st.columns([2, 2])
                with _scan_col_h:
                    st.markdown(f"#### Scanner {_scan_badge}")
                with _scan_col_note:
                    if location not in _scanner_results and location not in _scanner_tasks_map:
                        st.caption("auto-runs after completeness")
                    elif location in _scanner_tasks_map and location not in _scanner_results:
                        st.caption("⟳ scanning…")

                if location in _scanner_results:
                    _stale2 = set(_scanner_panels.get(location, [])) != set(all_selected_variants)
                    if _stale2:
                        st.warning("Panel changed — re-run analysis to refresh scanner.")

                    # scanner is read-only — no add button
                    render_scanner_results(
                        _scanner_results[location],
                        selected_variants=all_selected_variants,
                        client=wiseLoculus,
                        location=location,
                        date_range=(
                            datetime.combine(start_date, datetime.min.time()),
                            datetime.combine(end_date, datetime.min.time()),
                        ),
                        on_add_variant=None,
                    )

                elif location in _scanner_tasks_map:
                    _stt = celery_app.AsyncResult(_scanner_tasks_map[location])
                    if _stt.state in ("PENDING", "STARTED", "RETRY"):
                        st.info("⟳ Scanner running — auto-runs after completeness…")
                    elif _stt.state == "FAILURE":
                        st.error("Scanner failed — check worker logs.")
                else:
                    st.caption("Scanner will run automatically after completeness.")

                # ── Per-city panel additions ──────────────────────────────────
                st.markdown("---")
                _short2 = location.split("(")[0].strip()
                st.markdown(f"#### Refine panel for {_short2}")

                # get scanner suggestions for this city
                _suggestions = []
                if location in _scanner_results:
                    _suggestions = [
                        item["variant"]
                        for item in _scanner_results[location].get("missing_from_panel", [])
                        if item["variant"] not in all_selected_variants
                    ][:5]  # top 5 suggestions

                _city_extras = st.session_state.get("acooc_city_extras", {})
                _current_extras = _city_extras.get(location, [])

                st.caption(
                    f"Add variants for {_short2} only — these extend the base panel ({len(all_selected_variants)} variants) "
                    f"without affecting other cities."
                )

                # clickable suggestion pills
                if _suggestions:
                    st.markdown(
                        "<span style='font-size:11px;color:#898781;'>Scanner suggestions:</span>",
                        unsafe_allow_html=True,
                    )
                    _pill_cols = st.columns(min(5, len(_suggestions)))
                    for _pi, _sug in enumerate(_suggestions):
                        with _pill_cols[_pi % 5]:
                            if _sug not in _current_extras:
                                if st.button(f"+ {_sug}", key=f"acooc_suggest_{location}_{_sug}", use_container_width=True):
                                    _city_extras[location] = _current_extras + [_sug]
                                    st.session_state["acooc_city_extras"] = _city_extras
                            else:
                                st.markdown(
                                    f"<div style='font-size:11px;padding:4px 8px;background:#EAF3DE;"
                                    f"color:#3B6D11;border-radius:6px;text-align:center;'>✓ {_sug}</div>",
                                    unsafe_allow_html=True,
                                )

                # multiselect for local additions (includes suggestions + any pango lineage)
                _pango_loader_local = cached_get_pango_loader()
                _all_lineages = sorted(_pango_loader_local.get_raw_data().keys())
                _extra_options = [v for v in _all_lineages if v not in all_selected_variants]
                _new_extras = st.multiselect(
                    "Additional variants",
                    options=_extra_options,
                    default=_current_extras,
                    placeholder="Search any pango lineage…",
                    key=f"acooc_extras_ms_{location}",
                    label_visibility="collapsed",
                )
                # sync multiselect back to session state
                if _new_extras != _current_extras:
                    _city_extras[location] = _new_extras
                    st.session_state["acooc_city_extras"] = _city_extras
                    _current_extras = _new_extras

                # re-run button
                _combined_variants = all_selected_variants + [v for v in _current_extras if v not in all_selected_variants]
                _rerun_col1, _rerun_col2 = st.columns([3, 1])
                with _rerun_col1:
                    if _current_extras:
                        st.caption(f"Base ({len(all_selected_variants)}) + {len(_current_extras)} local = {len(_combined_variants)} variants total")
                    else:
                        st.caption("No local additions — will use base panel only")
                with _rerun_col2:
                    if st.button(f"▶ Re-run {_short2}", key=f"acooc_rerun_{location}", type="primary", use_container_width=True):
                        _rt = celery_app.send_task(
                            "tasks.run_deconvolve_lapis",
                            kwargs={
                                "location": location,
                                "start_date": start_date.isoformat(),
                                "end_date": end_date.isoformat(),
                                "variants": _combined_variants,
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
                                "variants": _combined_variants,
                            }
                        )
                        st.session_state["acooc_location_tasks"][location] = _rt.id
                        st.session_state["acooc_cooc_tasks"][location] = _rct.id
                        st.session_state["location_results"].pop(location, None)
                        st.session_state["acooc_cooc_results"].pop(location, None)
                        st.session_state["acooc_scanner_results"].pop(location, None)
                        # clear scanner task so it re-runs after new completeness
                        _sct2 = st.session_state.get("acooc_scanner_tasks", {})
                        _sct2.pop(location, None)
                        st.session_state["acooc_scanner_tasks"] = _sct2
                        logger.info(f"Re-run {location} with {len(_combined_variants)} variants ({len(_current_extras)} local extras)")

            _active_loc = _selected.replace("📍 ", "")
            if _active_loc in location_tasks:
                _city_tab_content(_active_loc, location_tasks[_active_loc])


            # ── Download report (triggered by button in progress header) ───────
            if st.session_state.get("acooc_show_report"):
                _completed_locs2 = [loc for loc in location_names if loc in st.session_state.get("location_results",{})]
                if _completed_locs2:
                    import plotly.graph_objects as go
                    from plotly.subplots import make_subplots
                    import math
                    _ncols = 2
                    _nrows = math.ceil(len(_completed_locs2)/_ncols)
                    _fig = make_subplots(rows=_nrows, cols=_ncols, subplot_titles=_completed_locs2,
                        vertical_spacing=0.12, horizontal_spacing=0.08)
                    for _i, _loc in enumerate(_completed_locs2):
                        _row = _i//_ncols+1; _col = _i%_ncols+1
                        _res = st.session_state.location_results[_loc]
                        for _variant, _data in _res.items():
                            if _variant == "undetermined": continue
                            _ts = _data.get("timeseriesSummary",[])
                            if not _ts: continue
                            _fig.add_trace(go.Scatter(
                                x=[e.get("date") for e in _ts],
                                y=[e.get("proportion",0) for e in _ts],
                                name=_variant, mode="lines+markers",
                                marker=dict(size=4), showlegend=(_i==0)),
                                row=_row, col=_col)
                    _fig.update_layout(height=300*_nrows, template="plotly_white",
                        title_text="Variant Proportion Estimates — All Locations")
                    st.plotly_chart(_fig, use_container_width=True)
                    import io
                    try:
                        _buf = io.BytesIO()
                        _fig.write_image(_buf, format="pdf")
                        st.download_button("⬇ Download PDF", data=_buf.getvalue(),
                            file_name="vpipe_scout_report.pdf", mime="application/pdf",
                            key="acooc_download_pdf")
                    except Exception:
                        st.caption("Install kaleido for PDF export.")
                if st.button("Close report", key="acooc_close_report"):
                    st.session_state["acooc_show_report"] = False


if __name__ == "__main__":
    app()