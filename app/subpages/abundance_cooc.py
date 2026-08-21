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


@st.cache_data(ttl=3600)
def cached_fetch_locations() -> list:
    """Cache locations for an hour so transient LAPIS failures or reruns don't
    empty the options list (which would silently drop the user's selection)."""
    return wiseLoculus.fetch_locations()


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
    # Apply any pending variant additions from the scanner BEFORE widgets render.
    # Scanner-found variants may not be in the curated cowwid list, so they go
    # into the manual-add multiselect (acooc_extra_variants), which accepts any
    # pango lineage. Curated ones could go either way; extra_variants is safe.
    _pending_changed = False
    # removals first (for "track instead" swaps: remove parent, add sublineage)
    for key in list(st.session_state.keys()):
        if key.startswith("acooc_remove_variant_pending_"):
            v = st.session_state.pop(key)
            _cm = st.session_state.get("acooc_variant_multiselect", [])
            if v in _cm:
                st.session_state["acooc_variant_multiselect"] = [x for x in _cm if x != v]
            _ce = st.session_state.get("acooc_extra_variants", [])
            if v in _ce:
                st.session_state["acooc_extra_variants"] = [x for x in _ce if x != v]
            _pending_changed = True
    # then additions
    for key in list(st.session_state.keys()):
        if key.startswith("acooc_add_variant_pending_"):
            v = st.session_state.pop(key)
            _curated = cached_get_variant_names()
            if v in _curated:
                _cur = st.session_state.get("acooc_variant_multiselect", [])
                if v not in _cur:
                    st.session_state["acooc_variant_multiselect"] = _cur + [v]
            else:
                _cur = st.session_state.get("acooc_extra_variants", [])
                if v not in _cur:
                    st.session_state["acooc_extra_variants"] = _cur + [v]
            _pending_changed = True
    if _pending_changed:
        st.rerun()

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
            # options must ALWAYS include whatever is currently in session state,
            # otherwise Streamlit silently drops a value not in options — which is
            # exactly what happens when the scanner adds a variant: it lands in
            # session state, but if it's filtered out of options the widget drops
            # it and the panel count doesn't grow (breaking the Run button).
            _current_extras = st.session_state.get("acooc_extra_variants", [])
            extra_options = sorted(set(
                [v for v in available_lineages if v not in selected_variants]
                + list(_current_extras)
            ))
            extra_variants = st.multiselect(
                "Search pango lineage",
                options=extra_options,
                placeholder="e.g. KP.2.3",
                help="Add any pango lineage. Scanner suggestions below highlight missing variants automatically.",
                key="acooc_extra_variants",
            )

        # dedup: a scanner variant routed to extras might also be curated
        all_selected_variants = list(dict.fromkeys(selected_variants + extra_variants))
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
        )

        st.markdown("---")

        # ── Step 3: Locations ─────────────────────────────────────────────────
        _step_label(3, "Locations", done=has_locations)
        available_locations = cached_fetch_locations()
        # guard: if the fetch transiently returned empty, fall back to whatever
        # the user already selected so their choice isn't silently dropped
        if not available_locations:
            available_locations = st.session_state.get("acooc_location_multiselect", [])
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
        # Persist a stable copy of the selection. On a mid-rerun (e.g. right
        # after adding a scanner variant), the widget can transiently return
        # empty if its options briefly mismatch; recover from the stable copy
        # so can_run / re-run don't see locs=0 and disable the Run button.
        if selected_locations:
            st.session_state["acooc_locations_stable"] = list(selected_locations)
        elif st.session_state.get("acooc_locations_stable"):
            selected_locations = st.session_state["acooc_locations_stable"]

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

        # scanner state
        scanner_results = st.session_state.get("acooc_scanner_results", {})
        scanner_panels = st.session_state.get("acooc_scanner_panels", {})

        # is anything currently running? (guards Run button)
        # A task counts as "running" only if it's in a live state AND its result
        # hasn't been collected yet. Celery reports expired/unknown task IDs as
        # PENDING forever, so we cross-check against the results dicts to avoid
        # a stale PENDING keeping the button disabled after everything finished.
        def _any_running():
            _lr = st.session_state.get("location_results", {})
            _cr = st.session_state.get("acooc_cooc_results", {})
            _sr = st.session_state.get("acooc_scanner_results", {})
            _checks = [
                (st.session_state.get("acooc_location_tasks", {}), _lr),
                (st.session_state.get("acooc_cooc_tasks", {}), _cr),
                (st.session_state.get("acooc_scanner_tasks", {}), _sr),
            ]
            for _tasks, _results in _checks:
                for _loc, _t in _tasks.items():
                    if _loc in _results:
                        continue  # result already collected — not running
                    if not _t:
                        continue
                    _state = celery_app.AsyncResult(_t).state
                    # PENDING/STARTED/RETRY without a collected result = running.
                    # The results cross-check above already skipped tasks whose
                    # results we have, so a stale expired-PENDING that already
                    # produced a result won't reach here. This keeps the Run
                    # button disabled while work is genuinely in flight.
                    if _state in ("PENDING", "STARTED", "RETRY"):
                        return True
            return False
        _busy = _any_running()

        # ── Step 5: Run analysis ──────────────────────────────────────────────
        can_run = (
            len(all_selected_variants) >= 2
            and len(selected_locations) >= 1
            and end_date > start_date
            and not _busy
        )
        _step_label(5, "Run analysis", done=has_run, active=not has_run and can_run)
        if st.button(
            "▶ Run analysis",
            type="primary",
            disabled=not can_run,
            key="acooc_run_button",
            use_container_width=True,
            help="Runs deconvolution, completeness, and scanner for all locations."
                 if not _busy else "Wait for the current analysis to finish.",
        ):
            st.session_state["acooc_trigger_run"] = True
        st.caption("runs all cities with base panel" if not _busy else "⟳ analysis in progress…")
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
            st.session_state["acooc_cooc_results"] = {}
            # clear scanner state too so it re-scans fresh (was showing stale 6/6)
            st.session_state["acooc_scanner_results"] = {}
            st.session_state["acooc_scanner_tasks"] = {}
            st.session_state["acooc_scanner_panels"] = {}
            # reset bucket expand flags so scanner starts collapsed on a new run
            st.session_state["acooc_exp_missing"] = False
            st.session_state["acooc_exp_sub"] = False

        from components.multi_location_results import (
            render_single_location_result,
            render_location_progress,
        )

        location_tasks = st.session_state.get("acooc_location_tasks", {})

        if not location_tasks:
            st.info("Complete steps 1–5 on the left to see results here.")
        else:
            # collect completed results — track if anything new arrives this cycle
            _new_collected = False

            # deconvolution results
            _loc_res = st.session_state.get("location_results", {})
            for _loc, _tid in list(location_tasks.items()):
                if _loc not in _loc_res:
                    _t = celery_app.AsyncResult(_tid)
                    if _t.ready():
                        try:
                            _loc_res[_loc] = _t.get()
                            st.session_state["location_results"] = _loc_res
                            _new_collected = True
                        except Exception:
                            pass

            # scanner results
            scanner_tasks_map = st.session_state.get("acooc_scanner_tasks", {})
            for _loc, _tid in list(scanner_tasks_map.items()):
                if _loc not in scanner_results:
                    _t = celery_app.AsyncResult(_tid)
                    if _t.ready():
                        try:
                            scanner_results[_loc] = _t.get()
                            st.session_state["acooc_scanner_results"] = scanner_results
                            _new_collected = True
                            logger.info(f"Scanner results collected for {_loc}")
                        except Exception as _e:
                            logger.error(f"Scanner task failed for {_loc}: {_e}")

            # cooc results + auto-submit scanner
            _cooc_res = st.session_state.get("acooc_cooc_results", {})
            _scanner_tasks = st.session_state.get("acooc_scanner_tasks", {})
            for _loc, _tid in list(st.session_state.get("acooc_cooc_tasks", {}).items()):
                if _loc not in _cooc_res:
                    _t = celery_app.AsyncResult(_tid)
                    if _t.ready():
                        try:
                            _cooc_res[_loc] = _t.get()
                            st.session_state["acooc_cooc_results"] = _cooc_res
                            _new_collected = True
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

            # ── Autorefresh decision — AFTER collection + scanner submission ───
            # Base it on "is there outstanding work?" rather than raw task state,
            # so newly-submitted scanner tasks keep the refresh alive and the bars
            # update to green without needing a manual click.
            _lr_now = st.session_state.get("location_results", {})
            _cr_now = st.session_state.get("acooc_cooc_results", {})
            _sr_now = st.session_state.get("acooc_scanner_results", {})
            _outstanding = False
            for _ln in location_names:
                # deconv or completeness not yet collected → running
                if _ln not in _lr_now or _ln not in _cr_now:
                    _outstanding = True
                    break
                # completeness done but scanner not yet done → running
                if (_cr_now.get(_ln, {}).get("unexplained_patterns")
                        and _ln not in _sr_now):
                    _outstanding = True
                    break
            if _outstanding:
                from streamlit_autorefresh import st_autorefresh
                st_autorefresh(interval=3000, key="acooc_autorefresh")
            elif _new_collected:
                # final result(s) just arrived and nothing is left running —
                # force one full-page rerun so the left column (Run button) and
                # the progress bars reflect the completed state without a click.
                st.rerun()

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

            def _dot_color(done, running):
                return "#3B6D11" if done else ("#93C5FD" if running else "#E5E3DC")

            # counts for the three stages
            _lr = st.session_state.get("location_results", {})
            _n_deconv = sum(1 for l in location_names if l in _lr)
            _n_cooc = sum(1 for l in location_names if l in _cooc_res2)
            _n_scan = sum(1 for l in location_names if l in _scan_res2)
            _n_tot = len(location_names)

            with st.container():
                _ph1, _ph2 = st.columns([3, 1])
                with _ph1:
                    st.markdown("<div style='font-size:13px;font-weight:500;'>Analysis progress</div>",
                                unsafe_allow_html=True)
                with _ph2:
                    _completed_locs = [loc for loc in location_names if loc in _lr]
                    if _completed_locs:
                        if st.button("⬇ Download report", key="acooc_dl_report", use_container_width=True):
                            st.session_state["acooc_show_report"] = True

            def _agg_bar(name, n_done, color):
                _frac = n_done / _n_tot if _n_tot else 0
                st.markdown(
                    f"<div style='margin-bottom:8px;'>"
                    f"<div style='display:flex;justify-content:space-between;font-size:11px;margin-bottom:3px;'>"
                    f"<span style='font-weight:500;'>{name}</span>"
                    f"<span style='color:#898781;'>{n_done} / {_n_tot} cities</span></div>"
                    f"<div style='height:8px;background:#F1EFE8;border-radius:4px;overflow:hidden;'>"
                    f"<div style='height:100%;border-radius:4px;width:{_frac*100:.0f}%;background:{color};'></div>"
                    f"</div></div>",
                    unsafe_allow_html=True,
                )

            _agg_bar("Deconvolution", _n_deconv, "#185FA5")
            _agg_bar("Completeness", _n_cooc, "#3B6D11")
            _agg_bar("Scanner", _n_scan, "#EF9F27")

            # per-city checklist
            _chk_html = "<div style='display:flex;flex-wrap:wrap;gap:6px;margin-top:6px;padding-top:8px;border-top:0.5px solid rgba(0,0,0,.06);'>"
            for _loc in location_names:
                _pct, _n_miss, _sd, _sr, _dd, _cd, _dr, _cr = _city_status(_loc)
                # a city is "done" when deconv is done (primary deliverable)
                if _dd:
                    _icon, _ic = "✓", "#3B6D11"
                elif _dr or _cr or _sr:
                    _icon, _ic = "⟳", "#93C5FD"
                else:
                    _icon, _ic = "○", "#898781"
                _chk_html += (
                    f"<span style='display:flex;align-items:center;gap:5px;font-size:11px;"
                    f"padding:3px 9px;border-radius:20px;background:#faf9f5;"
                    f"border:0.5px solid rgba(0,0,0,.08);'>"
                    f"<span style='width:14px;height:14px;border-radius:50%;background:{_ic};"
                    f"color:#fff;display:flex;align-items:center;justify-content:center;"
                    f"font-size:9px;flex:none;'>{_icon}</span>"
                    f"{_loc.split('(')[0].strip()}</span>"
                )
            _chk_html += "</div>"
            st.markdown(_chk_html, unsafe_allow_html=True)

            st.markdown("<hr style='margin:10px 0 8px;opacity:.15;'>", unsafe_allow_html=True)

            # ── Tab strip ────────────────────────────────────────────────────
            _city_options = [f"📍 {loc}" for loc in location_names]
            if ("acooc_selected_city" not in st.session_state or
                    st.session_state.get("acooc_selected_city") not in _city_options):
                st.session_state["acooc_selected_city"] = _city_options[0] if _city_options else ""

            # Tab strip using actual buttons (guaranteed clickable, no grey-out)
            _tab_cols = st.columns(len(location_names))
            for _ti, _loc in enumerate(location_names):
                _p2, _nm2, _sd2, _sr2, _dd2, _cd2, _dr2, _cr2 = _city_status(_loc)
                _dc0 = _dot_color(_dd2, _dr2)  # deconvolution
                _dc1 = _dot_color(_cd2, _cr2)  # completeness
                _dc2 = _dot_color(_sd2, _sr2)  # scanner
                _is_on = st.session_state.get("acooc_selected_city","") == f"📍 {_loc}"
                _sn = _loc.split("(")[0].strip()
                _pct_t = _p2 if _p2 is not None else None
                _tab_label = f"{_sn} · {_pct_t}%" if _pct_t is not None else _sn
                with _tab_cols[_ti]:
                    if st.button(
                        _tab_label,
                        key=f"acooc_tab_{_loc}",
                        use_container_width=True,
                        type="primary" if _is_on else "secondary",
                    ):
                        st.session_state["acooc_selected_city"] = f"📍 {_loc}"

            _selected = st.session_state.get("acooc_selected_city", _city_options[0] if _city_options else "")
            st.markdown("<hr style='margin:6px 0 10px;opacity:.15;'>", unsafe_allow_html=True)

            def _city_tab_content(location, task_id):
                """Shared content for both active and idle city tab fragments."""
                _cooc_tasks_map = st.session_state.get("acooc_cooc_tasks", {})
                _cooc_results = st.session_state.get("acooc_cooc_results", {})
                _scanner_results = st.session_state.get("acooc_scanner_results", {})
                _scanner_panels = st.session_state.get("acooc_scanner_panels", {})
                _added_for = st.session_state.get("acooc_scanner_added_for", {})
                _scanner_tasks_map = st.session_state.get("acooc_scanner_tasks", {})

                # ── deconvolution (primary output) ────────────────────────────
                st.markdown("#### Variant deconvolution")
                st.caption("Primary output — estimated variant proportions over time.")
                if location in st.session_state.location_results:
                    render_single_location_result(
                        location, st.session_state.location_results[location]
                    )
                else:
                    render_location_progress(
                        location, task_id, celery_app, redis_client
                    )

                st.markdown("---")

                # ── completeness (guiding indicator) ──────────────────────────
                st.markdown("#### Panel completeness")
                st.caption("Guiding indicator — how much of the co-occurrence signal your panel explains.")
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
                    st.caption("Runs alongside deconvolution.")

                # ── Jaccard collapsed ─────────────────────────────────────────
                if len(all_selected_variants) >= 2:
                    with st.expander("Signature similarity (Jaccard)", expanded=False):
                        render_jaccard_heatmap(
                            variants=all_selected_variants,
                            pango_loader=cached_get_pango_loader(),
                        )

                # scanner results shown in a combined summary below all city tabs
                # (not per-tab) — see the "Scanner findings" section after the tabs

            _active_loc = _selected.replace("📍 ", "")
            if _active_loc in location_tasks:
                _city_tab_content(_active_loc, location_tasks[_active_loc])

            # ── Scanner (one section, aggregated across all cities) ────────────
            _scan_res_all = st.session_state.get("acooc_scanner_results", {})
            if _scan_res_all:
                from collections import defaultdict as _ddict
                st.markdown("---")
                st.markdown("### Scanner")
                st.caption(
                    "Diagnostic across all cities — expand a category to see findings "
                    "and add variants. Adding applies to the whole panel; re-run to apply."
                )

                # coverage caption (instant, no scan needed): panel ∩ OT vs all OT
                _ot_set = set(cached_get_variant_names())
                _n_panel_ot = sum(1 for _v in all_selected_variants if _v in _ot_set)
                st.markdown(
                    f"<div style='border:0.5px solid #BFD9F2;background:#EFF6FF;border-radius:8px;"
                    f"padding:8px 12px;margin:2px 0 10px;font-size:12px;color:#1E3A5F;'>"
                    f"<b>Panel coverage:</b> {_n_panel_ot} of {len(_ot_set)} officially tracked "
                    f"variants selected. The scanner ranks the missing ones by how much "
                    f"co-occurrence signal they actually have in your samples.</div>",
                    unsafe_allow_html=True,
                )

                def _chip_html(cities):
                    return "".join(
                        f"<span style='display:inline-block;font-size:10px;padding:1px 6px;"
                        f"border-radius:4px;background:#F1EFE8;color:#5F5E5A;margin:1px 2px 1px 0;'>"
                        f"{c.split('(')[0].strip()}</span>"
                        for c in cities
                    )

                # ---- Bucket 1 (reframed): Missing variants WITH SIGNAL ----
                # Aggregate each missing variant across cities (sum reads, collect
                # cities), then group by the scanner's stable cluster_key so
                # near-identical OT lineages matched by the same reads collapse to
                # ONE representative row. cluster_key is identical across cities
                # (see scanner._cluster_missing_ot), so grouping is consistent;
                # the representative is picked here by summed reads across cities.
                def _human_reads(_n):
                    if _n >= 1_000_000:
                        return f"{_n/1_000_000:.1f}M".replace(".0M", "M")
                    if _n >= 1_000:
                        return f"{_n/1_000:.0f}K"
                    return str(_n)

                _agg_miss = {}
                for _loc, _res in _scan_res_all.items():
                    for _item in _res.get("missing_from_panel", []):
                        _v = _item["variant"]
                        if _v in all_selected_variants:
                            continue
                        _slot = _agg_miss.setdefault(_v, {
                            "variant": _v, "reads": 0, "cities": [],
                            "cluster_key": _item.get("cluster_key", _v),
                        })
                        _slot["reads"] += int(_item.get("total_reads", 0))
                        _slot["cities"].append(_loc)

                _by_cluster = _ddict(list)
                for _rec in _agg_miss.values():
                    _by_cluster[_rec["cluster_key"]].append(_rec)

                _clusters = []
                for _members in _by_cluster.values():
                    # representative = most summed reads, tie-break most specific
                    # (deepest pango), then name — stable across cities.
                    _rep = max(_members, key=lambda r: (r["reads"], len(r["variant"].split(".")), r["variant"]))
                    _cl_cities = []
                    for _m in _members:
                        for _c in _m["cities"]:
                            if _c not in _cl_cities:
                                _cl_cities.append(_c)
                    _clusters.append({
                        "rep": _rep, "members": _members,
                        "cities": _cl_cities, "reads": _rep["reads"],
                        "size": len(_members),
                    })
                _clusters.sort(key=lambda c: -c["reads"])

                with st.expander(f"🔴 Missing variants with signal ({len(_clusters)} distinct)", expanded=st.session_state.get("acooc_exp_missing", False)):
                    if not _clusters:
                        st.caption("No missing tracked variants with signal.")
                    else:
                        st.caption("Officially tracked variants not in your panel that show real co-occurrence signal. Ranked by supporting reads.")
                        for _cl in _clusters:
                            _rep = _cl["rep"]
                            _badge = f"{_human_reads(_cl['reads'])} reads"
                            if _cl["size"] == 1:
                                _m1, _m2 = st.columns([4, 1])
                                with _m1:
                                    st.markdown(
                                        f"<div style='font-size:13px;font-weight:600;color:#dc2626;'>{_rep['variant']} "
                                        f"<span style='font-size:10px;font-weight:400;color:#92400E;'>· {_badge}</span></div>"
                                        f"<div style='margin-top:2px;'>{_chip_html(_cl['cities'])}</div>",
                                        unsafe_allow_html=True,
                                    )
                                with _m2:
                                    if st.button("＋ Add", key=f"acooc_addmiss_{_rep['variant']}", use_container_width=True):
                                        st.session_state[f"acooc_add_variant_pending_{_rep['variant']}"] = _rep["variant"]
                                        st.session_state["acooc_exp_missing"] = True
                                        st.rerun()
                            else:
                                _open_key = f"acooc_cluster_open_{_rep['variant']}"
                                _names = ", ".join(_m["variant"] for _m in _cl["members"])
                                _n_others = _cl["size"] - 1
                                st.markdown(
                                    f"<div style='border:0.5px solid #FCA5A5;background:#FEF2F2;"
                                    f"border-radius:8px;padding:8px 11px;margin:6px 0;'>"
                                    f"<div style='font-size:13px;font-weight:600;color:#991B1B;'>"
                                    f"{_rep['variant']} <span style='font-weight:400;color:#92400E;'>"
                                    f"+ {_n_others} near-identical lineage{'s' if _n_others != 1 else ''}"
                                    f"</span></div>"
                                    f"<div style='font-size:11px;color:#6b7280;margin-top:1px;'>"
                                    f"best match · {_badge}</div>"
                                    f"<div style='margin-top:3px;'>{_chip_html(_cl['cities'])}</div>"
                                    f"</div>",
                                    unsafe_allow_html=True,
                                )
                                st.caption(
                                    f"⚠ {_names} all matched by the same ~{_human_reads(_cl['reads'])} reads "
                                    f"(Jaccard >0.9). Counted once — adding more than one would destabilize deconvolution."
                                )
                                _a1, _a2 = st.columns([2, 3])
                                with _a1:
                                    if st.button(f"＋ Add {_rep['variant']}", key=f"acooc_addclust_{_rep['variant']}", use_container_width=True):
                                        st.session_state[f"acooc_add_variant_pending_{_rep['variant']}"] = _rep["variant"]
                                        st.session_state["acooc_exp_missing"] = True
                                        st.rerun()
                                with _a2:
                                    _lbl = "▾ hide cluster" if st.session_state.get(_open_key) else f"▸ show all {_cl['size']} lineages"
                                    if st.button(_lbl, key=f"acooc_clustbtn_{_rep['variant']}", use_container_width=True):
                                        st.session_state[_open_key] = not st.session_state.get(_open_key, False)
                                        st.session_state["acooc_exp_missing"] = True
                                        st.rerun()
                                if st.session_state.get(_open_key):
                                    _rows = "".join(
                                        f"<div style='padding:1px 0;'>• <b>{_m['variant']}</b> · {_human_reads(_m['reads'])} reads</div>"
                                        for _m in sorted(_cl["members"], key=lambda r: -r["reads"])
                                    )
                                    st.markdown(
                                        f"<div style='font-size:11px;color:#5F5E5A;margin:2px 0 6px 4px;'>{_rows}</div>",
                                        unsafe_allow_html=True,
                                    )

                # ---- Bucket 2: Emerging sublineages (aggregated) ----
                _agg_sub = {}  # lineage -> {parent, reads, cities, obs_muts}
                for _loc, _res in _scan_res_all.items():
                    for _item in _res.get("emerging_sublineage", []):
                        _lin = _item["lineage"]
                        if _lin not in _agg_sub:
                            _agg_sub[_lin] = {
                                "parent": _item["parent"],
                                "reads": _item["total_reads"],
                                "cities": [],
                                "obs": _item.get("observed_mutations", []),
                            }
                        _agg_sub[_lin]["cities"].append(_loc)
                _sub_list = sorted(_agg_sub.items(), key=lambda x: -x[1]["reads"])

                with st.expander(f"🟡 Emerging sublineages ({len(_sub_list)} found)", expanded=st.session_state.get("acooc_exp_sub", False)):
                    if not _sub_list:
                        st.caption("No emerging sublineages detected.")
                    else:
                        st.caption("Descendants of panel variants with rising signal. Highly similar to their parent — see options.")
                        for _lin, _d in _sub_list:
                            _parent = _d["parent"]
                            _in_panel = _lin in all_selected_variants
                            st.markdown(
                                f"<div style='border:0.5px solid #FDE68A;background:#FFFBEB;"
                                f"border-radius:8px;padding:8px 11px;margin:6px 0;'>"
                                f"<div style='font-size:13px;font-weight:600;color:#92400E;'>{_lin}</div>"
                                f"<div style='font-size:11px;color:#6b7280;margin-top:1px;'>"
                                f"sublineage of <b>{_parent}</b> · {_d['reads']:,} reads</div>"
                                f"<div style='margin-top:3px;'>{_chip_html(_d['cities'])}</div>"
                                f"</div>",
                                unsafe_allow_html=True,
                            )
                            _parent_in = _parent in all_selected_variants
                            if _in_panel:
                                st.caption(f"✓ {_lin} already in panel.")
                            else:
                                st.caption(
                                    f"⚠ Very similar to parent {_parent} — adding both may destabilize "
                                    f"deconvolution. 'Track instead' swaps {_parent} → {_lin}."
                                )
                                _b1, _b2, _b3 = st.columns([2, 2, 3])
                                with _b1:
                                    if _parent_in and st.button(f"⇄ Track instead of {_parent}", key=f"acooc_swap_{_lin}", use_container_width=True):
                                        st.session_state[f"acooc_remove_variant_pending_{_parent}"] = _parent
                                        st.session_state[f"acooc_add_variant_pending_{_lin}"] = _lin
                                        st.session_state["acooc_exp_sub"] = True
                                        st.rerun()
                                with _b2:
                                    if st.button("＋ Add anyway", key=f"acooc_addsub_{_lin}", use_container_width=True):
                                        st.session_state[f"acooc_add_variant_pending_{_lin}"] = _lin
                                        st.session_state["acooc_exp_sub"] = True
                                        st.rerun()

                # ---- Bucket 3: Possibly new (aggregated reads) ----
                _pn_total = sum(
                    _res.get("possibly_new", {}).get("total_reads", 0)
                    for _res in _scan_res_all.values()
                )
                with st.expander(f"🔵 Possibly new ({_pn_total:,} reads)", expanded=False):
                    if _pn_total == 0:
                        st.caption("No unexplained patterns without a known lineage.")
                    else:
                        st.caption(
                            f"{_pn_total:,} reads across all cities match no known lineage. "
                            "Could be a novel variant, recombinant, or artifact. No action — monitor."
                        )
                        for _loc, _res in _scan_res_all.items():
                            _pn = _res.get("possibly_new", {})
                            if _pn.get("total_reads", 0) > 0:
                                st.markdown(
                                    f"<span style='font-size:11px;'><b>{_loc.split('(')[0].strip()}</b>: "
                                    f"{_pn['total_reads']:,} reads, {_pn.get('pattern_count',0)} pattern(s)</span>",
                                    unsafe_allow_html=True,
                                )


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