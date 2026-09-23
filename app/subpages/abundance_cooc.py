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
from components.variant_explorer_ui import render_variant_explorer

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
    """Officially tracked variants — the curated reporting panel. Read from the
    local config (app/config/officially_tracked.yaml), NOT fetched from cowwid,
    so it reflects exactly the set you report on. Falls back to the cowwid list
    only if the config is missing."""
    import os
    import yaml
    _cfg = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                        "config", "officially_tracked.yaml")
    try:
        with open(_cfg) as _f:
            _data = yaml.safe_load(_f) or {}
        _vars = _data.get("officially_tracked") or []
        if _vars:
            return list(_vars)
    except Exception:
        pass
    from api.signatures import get_variant_names
    return get_variant_names()


@st.cache_data(ttl=3600)
def cached_fetch_locations() -> list:
    """Cache locations for an hour so transient LAPIS failures or reruns don't
    empty the options list (which would silently drop the user's selection)."""
    return wiseLoculus.fetch_locations()


def _render_composition(cooc_result: dict, scanner_result: dict, key: str = "",
                        show_legend: bool = True, show_caption: bool = True) -> None:
    """Normalized (0-100%) per-date composition: explained / addable / novel /
    noise. Green height = completeness; other bands = what the gap is made of.
    Empty/low-read dates dropped. Built from existing outputs — scanner untouched."""
    try:
        from process.completeness_composition import compute_completeness_composition
    except Exception:
        return
    rows = compute_completeness_composition(cooc_result, scanner_result or {})
    if not rows:
        st.caption("Not enough co-occurrence reads to show composition.")
        return
    import plotly.graph_objects as go
    dates = [r["date"] for r in rows]
    layers = [
        ("explained by panel", "explained_pct", "#0F6E56", "rgba(15,110,86,0.85)"),
        ("addable (not in panel)", "addable_pct", "#dc2626", "rgba(220,38,38,0.55)"),
        ("novel (investigate)", "novel_pct", "#2563eb", "rgba(37,99,235,0.45)"),
        ("unresolved / noise", "noise_pct", "#9ca3af", "rgba(156,163,175,0.40)"),
    ]
    fig = go.Figure()
    for name, _fld, line_c, fill_c in layers:
        fig.add_trace(go.Scatter(
            x=dates, y=[r[_fld] for r in rows], name=name, mode="lines",
            stackgroup="one", line=dict(width=0.5, color=line_c), fillcolor=fill_c,
            hovertemplate="%{x|%Y-%m-%d}<br>" + name + " %{y:.0%}<extra></extra>"))
    fig.update_layout(
        height=234, margin=dict(t=18, b=26, l=46, r=12),
        template="plotly_white",
        legend=dict(orientation="h", yanchor="bottom", y=1.0, x=0),
        showlegend=show_legend,
    )
    fig.update_yaxes(
        range=[0, 1], autorange=False, fixedrange=True,
        tickmode="array", tickvals=[0, 0.25, 0.5, 0.75, 1.0],
        ticktext=["0%", "25%", "50%", "75%", "100%"], title="share")
    st.plotly_chart(fig, use_container_width=True, key=f"comp_stack_{key}")
    if show_caption:
        st.caption(
            "Green = explained by your panel (its height = completeness) · "
            "red = addable (scanner found it) · blue = novel · grey = noise.")


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
    # Scanner-added variants live in their OWN plain session key (NOT a widget
    # key), so adding one never mutates the multiselect widgets' state — which
    # was wiping the OT panel on rerun. They're merged into the panel at build
    # time (see all_selected_variants below). Removals just drop from this list
    # or from the widget lists as appropriate.
    st.session_state.setdefault("acooc_scanner_added", [])
    _removal_changed = False
    for key in list(st.session_state.keys()):
        if key.startswith("acooc_remove_variant_pending_"):
            v = st.session_state.pop(key)
            st.session_state["acooc_scanner_added"] = [
                x for x in st.session_state.get("acooc_scanner_added", []) if x != v]
            _panel = st.session_state.get("acooc_panel", [])
            if v in _panel:
                st.session_state["acooc_panel"] = [x for x in _panel if x != v]
            _removal_changed = True
    for key in list(st.session_state.keys()):
        if key.startswith("acooc_add_variant_pending_"):
            v = st.session_state.pop(key)
            # add directly to the unified panel so it shows IN the multiselect.
            # Done here at the top, BEFORE the widget renders, so the widget picks
            # it up. (Also track in acooc_scanner_added for the "from scanner"
            # provenance / heatmap wiring.)
            _panel = st.session_state.get("acooc_panel", [])
            if v not in _panel:
                st.session_state["acooc_panel"] = _panel + [v]
            _sa = st.session_state.setdefault("acooc_scanner_added", [])
            if v not in _sa:
                _sa.append(v)
    if _removal_changed:
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
    has_variants = len(list(dict.fromkeys(
        st.session_state.get("acooc_panel", [])
        + st.session_state.get("acooc_scanner_added", [])))) >= 2
    has_locations = len(st.session_state.get("acooc_location_multiselect", [])) >= 1
    has_run = bool(st.session_state.get("acooc_location_tasks"))
    has_completeness = bool(st.session_state.get("acooc_cooc_results"))
    has_scanner = bool(st.session_state.get("acooc_scanner_results"))

    with col_controls:
        curated_variants = cached_get_variant_names()

        # ── Step 1: Variant panel ─────────────────────────────────────────────
        _step_label(1, "Variant panel", done=has_variants)

        pango_loader = cached_get_pango_loader()
        available_lineages = sorted(pango_loader.get_raw_data().keys())

        col_add_ot, col_clear = st.columns(2)
        with col_add_ot:
            if st.button("+ Add surveillance panel",
                         key="acooc_add_all_ot", use_container_width=True,
                         help="Add the officially-tracked surveillance panel (7 variants)"):
                _cur = st.session_state.get("acooc_panel", [])
                st.session_state["acooc_panel"] = list(dict.fromkeys(_cur + curated_variants))
                st.rerun()
        with col_clear:
            if st.button("Clear", key="acooc_clear_variants", use_container_width=True):
                st.session_state["acooc_panel"] = []
                st.session_state["acooc_scanner_added"] = []
                st.rerun()

        # ONE multiselect for the whole panel — any pango lineage. Scanner-added
        # variants are merged in (below) so everything lives in one list the user
        # can edit. OT variants are the curated set; the "Add all OT" button is a
        # convenience to load them. options always include current selections +
        # scanner-added so Streamlit never silently drops a value.
        _scanner_added = st.session_state.get("acooc_scanner_added", [])
        _cur_panel = st.session_state.get("acooc_panel", [])
        _panel_options = sorted(set(available_lineages) | set(_cur_panel)
                                | set(_scanner_added) | set(curated_variants))
        selected_variants = st.multiselect(
            "Variant panel",
            options=_panel_options,
            default=None,
            help="Add any pango lineage. Use “Add surveillance panel” to load the "
                 "officially-tracked variants; the tree marks them “OT”. Scanner "
                 "findings you add appear here too.",
            key="acooc_panel",
            placeholder="Search any pango lineage (e.g. KP.2.3, XFG)…",
        )

        # scanner-added variants are written into acooc_panel directly (above),
        # so the multiselect selection IS the full panel — single source of truth.
        all_selected_variants = list(dict.fromkeys(list(selected_variants)))
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
            st.warning("Could not fetch locations from the API; using your current selection.")
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
        # quick date-range presets — anchored on the latest available data date
        # (max_date), clamped to the available window [min_date, max_date].
        import datetime as _dt
        def _acooc_set_range(_days):
            _end = max_date
            _start = (min_date if _days is None
                      else max(min_date, _end - _dt.timedelta(days=_days)))
            st.session_state["acooc_start_date"] = _start
            st.session_state["acooc_end_date"] = _end
            st.rerun()
        _pcol1, _pcol3, _pcolall = st.columns(3)
        with _pcol1:
            if st.button("Last month", key="acooc_preset_1m", use_container_width=True):
                _acooc_set_range(30)
        with _pcol3:
            if st.button("Last 3 months", key="acooc_preset_3m", use_container_width=True):
                _acooc_set_range(90)
        with _pcolall:
            if st.button("All", key="acooc_preset_all", use_container_width=True):
                _acooc_set_range(None)
        # seed once, then let the widgets read from session_state (the preset
        # buttons set it). No value= arg → no "default value but also set via
        # Session State" warning.
        st.session_state.setdefault("acooc_start_date", default_start)
        st.session_state.setdefault("acooc_end_date", default_end)
        col_start, col_end = st.columns(2)
        with col_start:
            start_date = st.date_input(
                "Start", min_value=min_date,
                max_value=max_date, key="acooc_start_date",
            )
        with col_end:
            end_date = st.date_input(
                "End", min_value=min_date,
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
        _busy = _any_running() or bool(st.session_state.get("acooc_trigger_run"))

        # ── Step 5: Run analysis ──────────────────────────────────────────────
        can_run = (
            len(all_selected_variants) >= 2
            and len(selected_locations) >= 1
            and end_date > start_date
            and not _busy
        )
        # ── timing estimate (benchmarked model) ───────────────────────────
        if end_date > start_date and selected_locations:
            _days  = (end_date - start_date).days
            # wastewater ~2 samples/week
            _n_dates = max(1, int(_days * 2 / 7))
            _n_locs  = len(selected_locations)
            # amp_dict positions from cutoff (2024-01-01 default ≈ 2141)
            _positions = 2141
            _per_date  = 0.30 + 0.00022 * _positions      # s / date / location
            _cooc_s    = _per_date * _n_dates * _n_locs * 1.2   # safety factor
            _deconv_s  = 90 * _n_locs                       # ~1.5 min/loc deconv
            _total_s   = _cooc_s + _deconv_s
            _total_min = _total_s / 60

            if _total_s <= 45:
                st.caption(
                    f"⏱ Estimated run time: ~{_total_s:.0f}s "
                    f"({_n_dates} sampling dates × {_n_locs} location(s))"
                )
            elif _total_s <= 120:
                st.info(
                    f"⏱ ~{_total_min:.1f} min — {_n_dates} dates × {_n_locs} location(s). "
                    "Runs in the background; you can keep working."
                )
            elif _total_s <= 300:
                st.warning(
                    f"⏱ ~{_total_min:.0f} min — this is a large query "
                    f"({_n_dates} dates × {_n_locs} locations). "
                    "Consider fewer locations or a shorter date range for faster results."
                )
            else:
                st.error(
                    f"⏱ ~{_total_min:.0f} min — very large query "
                    f"({_n_dates} dates × {_n_locs} locations). "
                    "Strongly consider narrowing the date range or selecting fewer locations."
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
            st.session_state["acooc_ran_panel"] = sorted(all_selected_variants)
            st.session_state["acooc_ran_locations"] = sorted(selected_locations)
            st.session_state["acooc_ran_dates"] = (start_date.isoformat(), end_date.isoformat())
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
            # Smart "re-run needed" warning. A re-run is needed only when the
            # existing results become stale/incomplete:
            #  - variants changed (deconv+scanner are computed for that panel)
            #  - date range changed (results are for that window)
            #  - a NEW city was added (it has no results yet)
            # Removing a city does NOT need a re-run — the remaining cities'
            # results are still valid (each city is computed independently).
            _ran = st.session_state.get("acooc_ran_panel")
            _ran_locs = st.session_state.get("acooc_ran_locations", [])
            _ran_dates = st.session_state.get("acooc_ran_dates")
            _now_dates = (start_date.isoformat(), end_date.isoformat())
            _reasons = []
            if _ran is not None and sorted(all_selected_variants) != _ran:
                _reasons.append("variant panel changed")
            if _ran_dates is not None and _now_dates != _ran_dates:
                _reasons.append("date range changed")
            _added_cities = [c for c in selected_locations if c not in _ran_locs]
            if _added_cities:
                _reasons.append(f"new location(s) added ({', '.join(_added_cities)})")
            if _reasons:
                st.warning("⚠ Re-run needed — " + "; ".join(_reasons)
                           + ". Results below reflect the previous run.")
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
                _n_miss = len(_scan_res2.get(loc,{}).get("resolved_clade",[])) if _scan_done else 0
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
                        # build the deconvolution zip inline so the button
                        # downloads directly (no intermediate report view)
                        import io as _io, csv as _csv, zipfile as _zip
                        _zbuf = _io.BytesIO(); _n_dl = 0
                        with _zip.ZipFile(_zbuf, "w", _zip.ZIP_DEFLATED) as _zf:
                            for _loc in _completed_locs:
                                _res = st.session_state.location_results[_loc]
                                # deconv result is wrapped by location name; unwrap
                                if isinstance(_res, dict) and _loc in _res and isinstance(_res[_loc], dict):
                                    _res = _res[_loc]
                                _rows = []
                                for _variant, _data in _res.items():
                                    if _variant == "undetermined":
                                        continue
                                    for _e in (_data.get("timeseriesSummary", []) or []):
                                        _rows.append({
                                            "location": _loc, "variant": _variant,
                                            "date": _e.get("date", ""),
                                            "proportion": _e.get("proportion", ""),
                                            "proportion_lower": _e.get("proportionLower", _e.get("ci_lower", "")),
                                            "proportion_upper": _e.get("proportionUpper", _e.get("ci_upper", "")),
                                        })
                                if not _rows:
                                    continue
                                _sio = _io.StringIO()
                                _w = _csv.DictWriter(_sio, fieldnames=["location","variant","date","proportion","proportion_lower","proportion_upper"])
                                _w.writeheader(); _w.writerows(_rows)
                                _safe = _loc.split("(")[0].strip().replace(" ", "_")
                                _zf.writestr(f"deconvolution_{_safe}.csv", _sio.getvalue())
                                _n_dl += 1
                        if _n_dl:
                            st.download_button(
                                "⬇ Download deconvolution (CSV)",
                                data=_zbuf.getvalue(),
                                file_name="vpipe_scout_deconvolution.zip",
                                mime="application/zip",
                                key="acooc_download_zip",
                                use_container_width=True)

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
            _agg_bar("Scanner", _n_scan, "#EF9F27")

            # per-city checklist
            _chk_html = "<div style='display:flex;flex-wrap:wrap;gap:6px;margin-top:6px;padding-top:8px;border-top:0.5px solid rgba(0,0,0,.06);'>"
            for _loc in location_names:
                _pct, _n_miss, _sd, _sr, _dd, _cd, _dr, _cr = _city_status(_loc)
                # a city is "done" only when ALL stages are done — deconv AND
                # completeness AND scanner. Showing green on deconv alone is
                # misleading: the scanner findings / addable band aren't ready yet.
                if _dd and _cd and _sd:
                    _icon, _ic = "✓", "#3B6D11"
                elif _dr or _cr or _sr or _dd or _cd:
                    # any stage running, or an earlier stage done but later ones
                    # still pending → in progress
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

            # ── Section switcher (top-level tabs) ─────────────────────────────
            # Buttons persist selection across autorefresh (st.tabs ghosted +
            # reset to the first tab on each rerun). Styled as underline tabs (via
            # CSS on their container) so they read as navigation, not as another
            # row of city buttons.
            st.session_state.setdefault("acooc_section", "Deconvolution results")
            _sections = ["Deconvolution results", "Co-occurrence results", "Investigate a variant"]
            _scols = st.columns(len(_sections))
            for _si, _snm in enumerate(_sections):
                with _scols[_si]:
                    _active = st.session_state["acooc_section"] == _snm
                    if st.button(_snm, key=f"acooc_sec_{_si}",
                                 use_container_width=True,
                                 type="primary" if _active else "secondary"):
                        st.session_state["acooc_section"] = _snm
                        st.rerun()
            _active_section = st.session_state["acooc_section"]
            st.markdown("<hr style='margin:2px 0 10px;opacity:.12;'>", unsafe_allow_html=True)

            if _active_section == "Deconvolution results":
              # ── City selector (radio — lighter than the section tabs above) ───
              _city_options = [f"📍 {loc}" for loc in location_names]
              if ("acooc_selected_city" not in st.session_state or
                      st.session_state.get("acooc_selected_city") not in _city_options):
                  st.session_state["acooc_selected_city"] = _city_options[0] if _city_options else ""
              if len(location_names) > 1:
                  _city_labels = [loc.split("(")[0].strip() for loc in location_names]
                  _cur_idx = 0
                  _cur_city = st.session_state.get("acooc_selected_city", "")
                  for _ci, _loc in enumerate(location_names):
                      if f"📍 {_loc}" == _cur_city:
                          _cur_idx = _ci
                  _pick = st.radio("City", options=list(range(len(location_names))),
                                   format_func=lambda i: _city_labels[i],
                                   index=_cur_idx, horizontal=True,
                                   key="acooc_city_radio", label_visibility="collapsed")
                  st.session_state["acooc_selected_city"] = f"📍 {location_names[_pick]}"

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

                  # ── co-occurrence check per deconvolution variant ─────────────
                  # Annotate each deconvolution result with whether co-occurrence
                  # corroborates it (confirmed / oscillating / can't-confirm).
                  if location in st.session_state.location_results:
                      _dec = st.session_state.location_results[location]
                      # deconv result is wrapped by location name:
                      # {loc: {variant: {timeseriesSummary...}}}. Unwrap to the
                      # inner variant dict. Handle both wrapped and flat shapes.
                      if isinstance(_dec, dict) and location in _dec and isinstance(_dec[location], dict):
                          _dec = _dec[location]
                      # Corroboration is about PANEL variants and depends only on
                      # deconvolution + the co-occurrence completeness pipeline
                      # (panel_confirmations) — NOT the scanner (which finds
                      # non-panel variants for the separate discovery section).
                      _found_nodes = set()
                      _cooc_res = _cooc_results.get(location)
                      if _cooc_res is None and location in _cooc_tasks_map:
                          _t = celery_app.AsyncResult(_cooc_tasks_map[location])
                          if _t.ready():
                              try:
                                  _cooc_res = _t.get()
                                  _cooc_results[location] = _cooc_res
                                  st.session_state["acooc_cooc_results"] = _cooc_results
                              except Exception:
                                  _cooc_res = None
                      # panel variant confirmed if its discriminating haplotype was
                      # observed co-occurring (panel_confirmations from cooc pipeline)
                      _panel_conf = (_cooc_res or {}).get("panel_confirmations", {}) or {}
                      _MIN_CONFIRM_READS = 100
                      for _pv, _reads in _panel_conf.items():
                          if _reads and _reads >= _MIN_CONFIRM_READS:
                              _found_nodes.add(_pv)
                      # only render once the cooc completeness result is available,
                      # so verdicts don't flip as the scanner streams in
                      _cooc_ready = _cooc_res is not None
                      # header always shows so the section doesn't pop in late
                      st.markdown(
                          "<div style='font-weight:600;font-size:13px;"
                          "margin:6px 0 2px;'>Co-occurrence check</div>",
                          unsafe_allow_html=True)
                      st.caption(
                          "Does each panel variant's deconvolution abundance have "
                          "distinctive mutations co-occurring on reads to back it up? "
                          "Confirmed = yes; oscillating = swaps with a near-identical "
                          "relative (trust the sum); can't confirm = no distinctive "
                          "haplotype (blind spot); ▲ not found in WW = distinctive "
                          "haplotype was co-covered but never co-occurred in your data.")
                      if not _cooc_ready:
                          st.info("⏳ Waiting for the co-occurrence scan to finish…")
                      if _cooc_ready and all_selected_variants:
                          try:
                              from process.variant_annotation import (
                                  VariantIndex, oscillating_pairs, annotate_variant)
                              _sigs = st.session_state.get("acooc_all_sigs_cache")
                              if not _sigs:
                                  _pl0 = cached_get_pango_loader()
                                  _sigs = {lin: _pl0.get_signature(lin)
                                           for lin in _pl0.raw_data}
                                  st.session_state["acooc_all_sigs_cache"] = _sigs
                              if _sigs:
                                  _pl = cached_get_pango_loader()
                                  _pmap = {l: _pl.get_raw_data().get(l, {}).get("parent", "")
                                           for l in _sigs}
                                  _idx = VariantIndex(_sigs, _pmap)
                                  _osc = oscillating_pairs(all_selected_variants, _sigs)
                                  _dec_vars = [v for v in _dec.keys()
                                               if v != "undetermined" and v in _sigs]
                                  # verdict thresholds — calibrated on real
                                  # data (Alpha -> not_found; current variants
                                  # -> confirmed). Override in cooc_config.yaml.
                                  try:
                                      from utils.config import get_cooc_setting as _gcs
                                  except Exception:
                                      _gcs = None
                                  def _ppcfg(_k, _d):
                                      try:
                                          _val = _gcs(_k, _d) if _gcs else _d
                                          return _d if _val is None else _val
                                      except Exception:
                                          return _d
                                  _PP_CONF_FREQ = float(_ppcfg("presence.confirm_freq", 0.01))
                                  _PP_CONF_MINP = int(_ppcfg("presence.confirm_min_present", 10))
                                  _PP_ABS_COV = int(_ppcfg("presence.absent_min_cov", 3000))
                                  _PP_ABS_MAXP = int(_ppcfg("presence.absent_max_present", 2))
                                  _PP_CON_FLOOR = float(_ppcfg("presence.con_floor", 0.5))
                                  _pp = (_cooc_res or {}).get("panel_presence", {}) or {}
                                  _rows = []
                                  for _v in _dec_vars:
                                      _pi = _pp.get(_v, {}) or {}
                                      _pres = int(_pi.get("present", 0) or 0)
                                      _cov = int(_pi.get("co_covered", 0) or 0)
                                      _nd = int(_pi.get("distinctive", 0) or 0)
                                      _con_p = int(_pi.get("con_present", 0) or 0)
                                      _con_t = int(_pi.get("con_testable", 0) or 0)
                                      _frac = (_pres / _cov) if _cov else 0.0
                                      _sib = _osc.get(_v, [])
                                      if _frac >= _PP_CONF_FREQ and _pres >= _PP_CONF_MINP:
                                          _status = "confirmed"
                                          _reason = f"co-occurs in {_frac*100:.1f}% of {_cov:,} covered reads"
                                      elif _nd < 2:
                                          _status = "cant_confirm"
                                          _reason = ("no distinctive haplotype (blind spot)"
                                                     + (f" — near-identical to {', '.join(_sib)}; "
                                                        "trust the sum" if _sib else ""))
                                      elif _cov < _PP_ABS_COV:
                                          # co-occurrence blind (distinctive muts
                                          # don't co-occur) -> constellation fallback
                                          if _con_t < 2:
                                              _status = "cant_confirm"
                                              _reason = "too few readable distinctive positions — no data"
                                          elif (_con_p / _con_t) < _PP_CON_FLOOR:
                                              _status = "not_found"
                                              _reason = f"constellation absent — {_con_p}/{_con_t} distinctive mutations present"
                                          else:
                                              _status = "detectable"
                                              _reason = f"{_con_p}/{_con_t} distinctive mutations present but not co-occurring — needs external check"
                                      elif _pres < _PP_ABS_MAXP:
                                          _status = "not_found"
                                          _reason = f"looked at {_cov:,} reads; haplotype not co-occurring — not found"
                                      elif _sib:
                                          _status = "oscillating"
                                          _reason = f"oscillates with {', '.join(_sib)} — trust the sum"
                                      else:
                                          _status = "detectable"
                                          _reason = f"weak co-occurrence ({_frac*100:.2f}%)"
                                      _ts = _dec.get(_v, {}).get("timeseriesSummary", [])
                                      _ab = ([e.get("proportion", 0) for e in _ts]
                                             if _ts else [])
                                      _abmean = (sum(_ab) / len(_ab)) if _ab else 0.0
                                      _rows.append((_v, _status, _reason, _abmean,
                                                    _cov, _frac, _nd, _con_p, _con_t))
                                  if _rows:
                                      # confirmed first, then oscillating, then blind
                                      _order = {"confirmed": 0, "oscillating": 1,
                                                "detectable": 2, "cant_confirm": 3,
                                                "not_found": 4}
                                      _rows.sort(key=lambda r: (_order.get(r[1], 3),
                                                                -r[3]))
                                      _meta = {
                                          "confirmed":   ("#0f6e56", "#e6f4ef", "✓ confirmed"),
                                          "oscillating": ("#ba7517", "#fdf4e6", "⚠ oscillating"),
                                          "detectable":  ("#6b7280", "#f3f4f6", "· detectable"),
                                          "cant_confirm":("#6b7280", "#f3f4f6", "· can't confirm"),
                                          "not_found":   ("#b45309", "#fff7ed", "▲ not found in WW"),
                                      }
                                      _th = ("padding:4px 8px;font-weight:600;"
                                             "color:#6b7280;text-align:left;")
                                      _html = (
                                          "<table style='width:100%;border-collapse:collapse;"
                                          "font-size:12px;'><thead><tr style='border-bottom:"
                                          "1px solid #e5e7eb;'>"
                                          f"<th style='{_th}'>Variant</th>"
                                          f"<th style='{_th}'>Abund.</th>"
                                          f"<th style='{_th}' title='distinctive mutations "
                                          "(unique within the panel)'>Distinct.</th>"
                                          f"<th style='{_th}' title='haplotype frequency "
                                          "(present / co-covered reads); blank when the "
                                          "distinctive mutations are too spread out to "
                                          "co-occur'>Co-occurrence</th>"
                                          f"<th style='{_th}' title='distinctive mutations "
                                          "present / testable individually'>Constellation</th>"
                                          f"<th style='{_th}'>Verdict</th>"
                                          "</tr></thead><tbody>")
                                      for _row in _rows:
                                          (_v, _s, _r, _ab, _cov, _frac,
                                           _nd, _con_p, _con_t) = _row
                                          _fg, _bg, _lbl = _meta.get(
                                              _s, ("#6b7280", "#f3f4f6", _s))
                                          if _cov >= _PP_ABS_COV:
                                              _cooc_cell = (
                                                  f"{_frac*100:.0f}% <span style='color:#9ca3af;'>"
                                                  f"({_cov:,})</span>")
                                          else:
                                              _cooc_cell = "<span style='color:#c9c7bf;'>—</span>"
                                          _con_cell = (f"{_con_p}/{_con_t}" if _con_t
                                                       else "<span style='color:#c9c7bf;'>—</span>")
                                          _rt = _r.replace("'", "&#39;")
                                          _td = "padding:5px 8px;border-bottom:0.5px solid #f0f0f0;"
                                          _html += (
                                              f"<tr title='{_rt}'>"
                                              f"<td style='{_td}font-weight:600;'>{_v}</td>"
                                              f"<td style='{_td}color:#6b7280;'>{_ab*100:.0f}%</td>"
                                              f"<td style='{_td}color:#9ca3af;'>{_nd}</td>"
                                              f"<td style='{_td}'>{_cooc_cell}</td>"
                                              f"<td style='{_td}'>{_con_cell}</td>"
                                              f"<td style='{_td}'><span style='background:{_bg};"
                                              f"color:{_fg};font-size:11px;font-weight:600;"
                                              f"padding:1px 8px;border-radius:10px;"
                                              f"white-space:nowrap;'>{_lbl}</span></td></tr>")
                                      _html += "</tbody></table>"
                                      st.markdown(_html, unsafe_allow_html=True)
                          except Exception as _e:
                              st.caption(f"(co-occurrence check unavailable: {_e})")

                  # ── Jaccard (signature similarity) ────────────────────────────
                  if len(all_selected_variants) >= 2:
                      st.markdown("---")
                      st.markdown("<div style='font-weight:600;font-size:13px;'>"
                                  "🧬 Signature similarity (Jaccard)</div>",
                                  unsafe_allow_html=True)
                      st.caption("How much each pair of panel variants shares mutations "
                                 "— high similarity means deconvolution may struggle to "
                                 "tell them apart.")
                      with st.expander("Show similarity heatmap", expanded=False):
                          render_jaccard_heatmap(
                              variants=all_selected_variants,
                              pango_loader=cached_get_pango_loader(),
                          )

                  # scanner results shown in a combined summary below all city tabs
                  # (not per-tab) — see the "Scanner findings" section after the tabs

              _active_loc = _selected.replace("📍 ", "")
              if _active_loc in location_tasks:
                  _city_tab_content(_active_loc, location_tasks[_active_loc])


            if _active_section == "Co-occurrence results":
              # panel-changed warning: scanner findings reflect the panel that was
              # run, so flag when the current selection differs.
              _ran_sc = st.session_state.get("acooc_ran_panel")
              _sc_stale = (
                  (_ran_sc is not None and sorted(all_selected_variants) != _ran_sc)
                  or (st.session_state.get("acooc_ran_dates") is not None
                      and (start_date.isoformat(), end_date.isoformat())
                      != st.session_state.get("acooc_ran_dates"))
                  or any(c not in st.session_state.get("acooc_ran_locations", [])
                         for c in selected_locations)
              )
              if _sc_stale:
                  st.warning("⚠ Re-run needed — the panel, dates, or locations "
                             "changed since these results were computed.")
              # ── Panel completeness (at top of scanner — shows what each city's
              #    signal is made of, before the findings that fill the gaps) ─────
              _cr_all = st.session_state.get("acooc_cooc_results", {})
              _sr_all = st.session_state.get("acooc_scanner_results", {})
              _ready_locs = [l for l in location_names
                             if _cr_all.get(l) is not None and _sr_all.get(l) is not None]
              if _ready_locs:
                  st.markdown("#### Panel completeness")
                  st.caption("How much of each city's co-occurrence signal your panel "
                             "explains (green) vs the rest. Green = explained · red = "
                             "addable (scanner found it) · blue = novel · grey = noise.")
                  # Consistent layout: reserve a slot for EVERY selected city
                  # so plots are the same size and labelled from the start,
                  # regardless of which city finishes first. Cities still
                  # computing show a placeholder in their slot. Legend on the
                  # first ready plot only.
                  _grid_locs = list(location_names)
                  _legend_used = [False]
                  def _one_city(_lc):
                      st.markdown(f"<div style='font-size:12px;font-weight:600;'>"
                                  f"{_lc}</div>", unsafe_allow_html=True)
                      if _cr_all.get(_lc) is not None and _sr_all.get(_lc) is not None:
                          _render_composition(_cr_all[_lc], _sr_all[_lc], key=_lc,
                                              show_legend=not _legend_used[0],
                                              show_caption=False)
                          _legend_used[0] = True
                      else:
                          st.caption("\u23f3 computing\u2026")
                  if len(_grid_locs) == 1:
                      _one_city(_grid_locs[0])
                  else:
                      for _i in range(0, len(_grid_locs), 2):
                          _cols = st.columns(2)
                          for _j, _lc in enumerate(_grid_locs[_i:_i+2]):
                              with _cols[_j]:
                                  _one_city(_lc)
                  st.markdown("---")

              # ── Scanner (one section, aggregated across all cities) ────────────
              _scan_res_all = st.session_state.get("acooc_scanner_results", {})
              # Show progress only while a scanner task is actively running. Once
              # no task is running, show whatever results were collected (avoids
              # hiding findings forever on a location-matching mismatch).
              # scanner is "running" if any location with a completeness result
              # doesn't yet have a collected scanner result (same reliable signal
              # the autorefresh uses — the task-map .ready() check was unreliable).
              # Show scanner findings only when NOTHING is still running (_outstanding
              # is the same signal that drives the autorefresh — it's True while any
              # deconv/completeness/scanner task is uncollected). The previous
              # _scan_running check relied on cooc results being present, but they
              # aren't collected yet at this point (cr_keys empty), so it never
              # triggered and partial findings leaked through.
              _scan_running = _outstanding
              st.markdown("---")
              st.markdown("### Scanner")
              st.caption(
                  "Variants circulating that your panel doesn't cover — across all "
                  "cities. Adding applies to the whole panel; re-run to apply.")
              if _scan_running:
                  st.info("🔍 Scanning for variants not in your panel… "
                          "findings and Add buttons appear when the scan completes.")
                  _scan_res_all = {}  # suppress partial rendering below

              if _scan_res_all:
                  from collections import defaultdict as _ddict
                  # (reads-threshold slider removed — the heatmap now defaults to
                  # regions with a discriminating mutation and offers a per-finding
                  # "show all regions" toggle instead.)

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

                  # ══ Option C scanner rendering ══════════════════════════════
                  # signatures for heatmap classification (lazy, cached in session)
                  _all_sigs = st.session_state.get("acooc_all_sigs_cache")
                  if _all_sigs is None:
                      try:
                          from api.pango_loader import PangoLoader, get_pango_summary_path
                          _pl = PangoLoader(get_pango_summary_path())
                          _all_sigs = {lin: _pl.get_signature(lin) for lin in _pl.raw_data}
                          st.session_state["acooc_all_sigs_cache"] = _all_sigs
                      except Exception:
                          _all_sigs = {}
                  # Aggregate the new clade-based scanner output across cities.
                  # Categories: not-in-panel clades, sublineage clades,
                  # unresolved, novel. Honest clade labels, member counts.
                  def _human_reads(_n):
                      if _n >= 1_000_000:
                          return f"{_n/1_000_000:.1f}M".replace(".0M", "M")
                      if _n >= 1_000:
                          return f"{_n/1_000:.0f}K"
                      return str(_n)

                  # collect clade findings across cities, keyed by node
                  _agg_new = {}      # node -> {reads, cities, member_count, members, designation, muts}
                  _agg_sub = {}
                  _agg_unres = {}    # frozenset(fp) -> {reads, cand, anc}
                  _agg_mnh = {}      # node -> matched-but-no-haplotype
                  _novel_total = 0
                  _novel_pats = 0

                  def _ingest_clade(_c, _loc):
                      # Route one finding (top-level OR a promoted sub-finding) into
                      # the new/sub bucket and aggregate its reads across cities.
                      _bucket = _agg_new if _c.get("relationship") == "new_lineage" else _agg_sub
                      _node = _c["node"]
                      _slot = _bucket.setdefault(_node, {
                          "node": _node, "reads": 0, "signal_reads": 0,
                          "cities": [],
                          "member_count": _c.get("member_count", 1),
                          "members": _c.get("members", []),
                          "designation": _c.get("designation", ""),
                          "panel_ancestor": _c.get("panel_ancestor", ""),
                          "muts": _c.get("observed_mutations", []),
                          "member_blocks": _c.get("member_blocks", []),
                          "shared_mutations": _c.get("shared_mutations", []),
                          "associated": _c.get("associated_members", []),
                          "confidence": _c.get("confidence", "weak"),
                          "verdict": _c.get("verdict", ""),
                          "trend": _c.get("trend", "flat"),
                          "trend_series": list(_c.get("trend_series", [])),
                          "peak_date": _c.get("peak_date", ""),
                          "trend_by_city": {},
                          "per_city": {},
                      })
                      _slot["reads"] += int(_c.get("total_reads", 0))
                      # signal_reads is the strongest discriminating region's
                      # reads; across cities take the MAX (summing would
                      # double-count the same discriminating reads and can
                      # exceed the total). Capped at total as a safety net.
                      _slot["signal_reads"] = max(
                          _slot.get("signal_reads", 0),
                          int(_c.get("signal_reads", 0)))
                      if _loc not in _slot["cities"]:
                          _slot["cities"].append(_loc)
                      # keep each city's own series so trend can be shown
                      # per-city or aggregated (summed element-wise)
                      _slot["trend_by_city"][_loc] = list(_c.get("trend_series", []))
                      # per-city discriminating (star) mutations: a mutation is
                      # discriminating if <= _STAR_MAX lineages carry it. Kept per
                      # city so the card shows which city has which. Regex-free
                      # leading-digit position sort.
                      _STAR_MAX = 30
                      def _pcpos(_m):
                          _d = ""
                          for _ch in _m:
                              if _ch.isdigit():
                                  _d += _ch
                              else:
                                  break
                          return int(_d) if _d else 0
                      _pc_car = {}
                      for _blk in _c.get("member_blocks", []):
                          _pc_car.update(_blk.get("mut_carriers", {}))
                      _pc_star = sorted(
                          {_m for _blk in _c.get("member_blocks", [])
                           for _m in _blk.get("discriminating", [])
                           if _pc_car.get(_m, 999) <= _STAR_MAX},
                          key=_pcpos)
                      _slot["per_city"][_loc] = {
                          "star": _pc_star,
                          "car": {_m: _pc_car.get(_m) for _m in _pc_star},
                      }

                  for _loc, _res in _scan_res_all.items():
                      for _c in _res.get("resolved_clade", []):
                          _ingest_clade(_c, _loc)
                          # A confirmed sub-finding (e.g. LF.7 nested under JN.1) is
                          # its own real finding — surface it as a card too, else a
                          # tiny parent (e.g. a 4k-read JN.1) hides a huge descendant
                          # (a 460k-read LF.7). Only promote sub-findings that are
                          # actually confirmed (have discriminating blocks); their
                          # reads are disjoint from the parent's, so no double count.
                          for _sf in _c.get("sub_findings", []) or []:
                              if _sf.get("member_blocks"):
                                  _ingest_clade(_sf, _loc)
                      for _u in _res.get("unresolved", []):
                          _k = tuple(_u["fingerprint"])
                          _s = _agg_unres.setdefault(_k, {
                              "fp": _u["fingerprint"], "reads": 0,
                              "cand": _u.get("candidate_count", 0),
                              "anc": _u.get("common_ancestor", ""),
                          })
                          _s["reads"] += int(_u.get("total_reads", 0))
                      for _mnh in _res.get("matched_no_haplotype", []):
                          _k = _mnh["node"]
                          _s = _agg_mnh.setdefault(_k, {
                              "node": _mnh["node"], "reads": 0,
                              "member_count": _mnh.get("member_count", 1),
                              "muts": _mnh.get("observed_mutations", []),
                          })
                          _s["reads"] += int(_mnh.get("total_reads", 0))
                      _nv = _res.get("novel", {})
                      _novel_total += int(_nv.get("total_reads", 0))
                      _novel_pats += int(_nv.get("pattern_count", 0))

                  def _clade_label(_slot):
                      return f"{_slot['node']} clade" if _slot["member_count"] > 1 else _slot["node"]

                  def _render_finding(_slot, _accent, _bg, _border, _is_sub=False):
                      _label = _clade_label(_slot)
                      _sig_reads = _slot.get("signal_reads", 0)
                      _tot_reads = _slot["reads"]
                      # headline = discriminating co-occurrence reads (the real
                      # signal); the broad fingerprint total is shown as context.
                      _reads = _sig_reads if _sig_reads > 0 else _tot_reads
                      _mc = _slot["member_count"]
                      _desig = _slot["designation"]
                      _members = _slot["members"]
                      _member_txt = ""
                      if _mc > 1:
                          _samp = ", ".join(_members[:4])
                          _more = f" +{_mc-4}" if _mc > 4 else ""
                          _member_txt = f" · {_mc} members: {_samp}{_more}"
                      _desig_txt = f" \u00b7 designated {_desig}" if _desig else ""
                      # ── per-city discriminating-mutation rows (no verdict /
                      # trend / sparkline / read-count headline). Show, per city,
                      # the discriminating (star) mutations found there; cities
                      # that share the same set are grouped on one row; a city
                      # with none is flagged backbone-only.
                      _per_city = _slot.get("per_city", {})
                      st.markdown(
                          f"<div style='background:{_bg};border:1px solid {_border};"
                          f"border-radius:6px;padding:8px 12px;margin:4px 0;'>"
                          f"<span style='font-weight:600;color:{_accent};'>{_label}</span>"
                          f"<span style='color:#6b7280;font-size:0.8rem;margin-left:8px;'>"
                          f"{_member_txt.lstrip(' \u00b7') if _member_txt else ''}"
                          f"{_desig_txt}</span></div>",
                          unsafe_allow_html=True,
                      )
                      _CAP_MUTS = 6
                      _groups = {}
                      for _cty in _slot.get("cities", []):
                          _pc = _per_city.get(_cty, {})
                          _key = tuple(_pc.get("star", []))
                          _groups.setdefault(_key, []).append(_cty)
                      for _star_key, _cts in sorted(_groups.items(),
                                                    key=lambda kv: -len(kv[0])):
                          _tags = "".join(
                              f"<span style='display:inline-block;font-size:10px;"
                              f"padding:1px 7px;border-radius:10px;background:{_bg};"
                              f"border:0.5px solid {_border};color:{_accent};"
                              f"margin:0 3px 2px 0;'>{_c2.split('(')[0].strip()}</span>"
                              for _c2 in _cts)
                          if _star_key:
                              # Option 1: show only the COUNT of discriminating
                              # mutations — the names aren't actionable on the card;
                              # the heatmap shows which ones over time.
                              _n = len(_star_key)
                              st.markdown(
                                  f"<div style='margin:2px 0 2px 4px;font-size:12px;"
                                  f"color:#374151;'>{_tags}"
                                  f"<span style='margin-left:4px;'>{_n} discriminating "
                                  f"mutation{'s' if _n != 1 else ''}</span></div>",
                                  unsafe_allow_html=True,
                              )
                          else:
                              st.markdown(
                                  f"<div style='margin:2px 0 2px 4px;font-size:12px;"
                                  f"color:#9ca3af;'>{_tags}"
                                  f"<span style='margin-left:4px;'>\u26a0 backbone only "
                                  f"\u2014 no discriminating mutation here</span></div>",
                                  unsafe_allow_html=True,
                              )
                      _assoc = _slot.get("associated", [])
                      if _assoc:
                          st.caption(
                              "⚠ includes recombinants sharing these mutations: "
                              + ", ".join(_assoc[:6])
                              + " — check the heatmap to see which is driving the signal."
                          )
                      if _is_sub:
                          _anc = _slot.get("panel_ancestor", "")
                          st.caption(f"↳ sublineage of {_anc} — signal partly counted within its proportion; adding refines it.")
                      _addkey = f"acooc_addc_{_slot['node']}"
                      if _scan_running:
                          st.caption(f"⏳ finishing scan… Add enables shortly")
                      else:
                          if st.button(f"+ Add {_slot['node']}", key=_addkey):
                              st.session_state[f"acooc_add_variant_pending_{_slot['node']}"] = _slot["node"]
                              st.rerun()  # process immediately (no autorefresh running once scan done)
                      # drill-down heatmap — let the user pick which city when the
                      # finding spans several (the finding's reads are summed across
                      # cities, but a heatmap shows one city's signal at a time).
                      if _slot["muts"] and wiseLoculus and _slot["cities"]:
                          _cities = _slot["cities"]
                          with st.expander(f"Signal over time — {_label}", expanded=False):
                              if len(_cities) > 1:
                                  _hloc = st.selectbox(
                                      "City",
                                      _cities,
                                      key=f"acooc_hmcity_{_slot['node']}",
                                      help="This finding was seen in several cities; "
                                           "pick which city's signal to display.",
                                  )
                              else:
                                  _hloc = _cities[0]
                              st.caption(f"Showing {_hloc}")
                              from components.scanner_heatmap import render_clade_heatmap
                              render_clade_heatmap(
                                  clade_node=_slot["node"],
                                  shared_mutations=_slot.get("shared_mutations", []),
                                  member_blocks=_slot.get("member_blocks", []),
                                  client=wiseLoculus,
                                  location=_hloc,
                                  date_range=(start_date, end_date),
                              )

                  # ---- Not in panel (addable) — merges new-lineage clades AND
                  #      sublineages of the panel. Both are "consider adding"; the
                  #      relationship (unrelated vs sublineage) is shown per-row.
                  _new_list = sorted(_agg_new.values(), key=lambda x: -x["reads"])
                  _sub_list = sorted(_agg_sub.values(), key=lambda x: -x["reads"])
                  _addable_n = len(_new_list) + len(_sub_list)
                  with st.expander(
                      f"🔴 Not in your panel — {_addable_n} finding(s)",
                      expanded=st.session_state.get("acooc_exp_missing", False),
                  ):
                      if not _addable_n:
                          st.caption("Nothing circulating that your panel doesn't cover.")
                      else:
                          st.caption(
                              "Co-occurrence signal your panel doesn't explain, mapped "
                              "to the tightest pango clade the mutations support. "
                              "Consider adding these to the panel."
                          )
                          # unrelated new lineages first (bigger gaps), then sublineages
                          for _slot in _new_list:
                              _render_finding(_slot, "#dc2626", "#fef2f2", "#fecaca")
                          for _slot in _sub_list:
                              _render_finding(_slot, "#dc2626", "#fef2f2", "#fecaca",
                                              _is_sub=True)

                  # ---- Matched a lineage but no discriminating co-occurrence ----
                  _mnh_list = sorted(_agg_mnh.values(), key=lambda x: -x["reads"])
                  if _mnh_list:
                      _mnh_reads = sum(m["reads"] for m in _mnh_list)
                      with st.expander(
                          f"🟣 Matched but not co-occurrence-confirmed — "
                          f"{len(_mnh_list)} lineage(s) · {_mnh_reads:,} reads",
                          expanded=False,
                      ):
                          st.caption(
                              "These lineages match some observed mutations, but their "
                              "distinguishing mutations don't co-occur on reads — so "
                              "co-occurrence can't confirm them (they may still be present; "
                              "deconvolution is the tool to quantify them)."
                          )
                          for _m in _mnh_list[:15]:
                              _lbl = f"{_m['node']} clade" if _m["member_count"] > 1 else _m["node"]
                              _muts = ", ".join(_m["muts"][:5])
                              st.markdown(
                                  f"<div style='background:#faf5ff;border:1px solid #e9d5ff;"
                                  f"border-radius:6px;padding:6px 10px;margin:3px 0;font-size:0.82rem;'>"
                                  f"<span style='font-weight:600;color:#7c3aed;'>{_lbl}</span>"
                                  f"<span style='color:#6b7280;margin-left:8px;'>"
                                  f"{_m['reads']:,} reads · matched: {_muts}</span></div>",
                                  unsafe_allow_html=True,
                              )

                  # ---- Unresolved ----
                  _unres_list = sorted(_agg_unres.values(), key=lambda x: -x["reads"])
                  if _unres_list:
                      _ur_reads = sum(u["reads"] for u in _unres_list)
                      with st.expander(
                          f"⚪ Unresolved — {len(_unres_list)} pattern(s) · {_ur_reads:,} reads",
                          expanded=False,
                      ):
                          st.caption(
                              "Co-occurring mutations match many lineages across "
                              "unrelated clades — too broad to name."
                          )
                          for _u in _unres_list[:15]:
                              _fp = ", ".join(_u["fp"][:5])
                              _anc = f" · nearest ancestor {_u['anc']}" if _u["anc"] else ""
                              st.markdown(
                                  f"<div style='background:#f8f9fa;border:1px solid #e5e7eb;"
                                  f"border-radius:6px;padding:6px 10px;margin:3px 0;font-size:0.82rem;'>"
                                  f"<span style='font-family:monospace;color:#374151;'>{_fp}</span>"
                                  f"<span style='color:#6b7280;margin-left:8px;'>"
                                  f"{_u['cand']} lineages · {_u['reads']:,} reads{_anc}</span></div>",
                                  unsafe_allow_html=True,
                              )

                  # ---- Novel ----
                  with st.expander(
                      f"🔵 Novel — no pango match ({_novel_total:,} reads)",
                      expanded=False,
                  ):
                      if _novel_total == 0:
                          st.caption("No unexplained patterns without a pango match.")
                      else:
                          st.caption(
                              f"{_novel_total:,} reads in {_novel_pats} pattern(s) match "
                              "no known lineage. Could be novel, recombinant, or artifact. "
                              "A rising pattern is the most worth investigating."
                          )
                          for _loc, _res in _scan_res_all.items():
                              _pn = _res.get("novel", {}) or _res.get("possibly_new", {})
                              for _pi, _pat in enumerate(_pn.get("top_patterns", [])[:5]):
                                  _pmuts = _pat.get("mutations", [])
                                  _muts = ", ".join(_pmuts)
                                  st.markdown(
                                      f"<div style='background:#eff6ff;border:1px solid #bfdbfe;"
                                      f"border-radius:6px;padding:6px 10px;margin:3px 0;font-size:0.8rem;'>"
                                      f"<span style='color:#1d4ed8;font-family:monospace;'>{_muts}</span>"
                                      f"<span style='color:#6b7280;margin-left:8px;'>"
                                      f"{_pat.get('count',0):,} · {_loc.split('(')[0].strip()}</span></div>",
                                      unsafe_allow_html=True,
                                  )
                                  # trend heatmap for this novel pattern — is it rising?
                                  if _pmuts and wiseLoculus and len(_pmuts) >= 2:
                                      with st.expander(
                                          f"Trend over time — {_loc.split('(')[0].strip()} — "
                                          f"{_muts[:40]}",
                                          expanded=False,
                                      ):
                                          from components.scanner_heatmap import render_clade_heatmap
                                          # carrier count per mutation (how many
                                          # lineages carry it) so the novel heatmap
                                          # can mark discriminating (★) rows instead
                                          # of labelling everything backbone. Uses
                                          # the signatures already cached on the page.
                                          _nov_sigs = st.session_state.get("acooc_all_sigs_cache") or {}
                                          _nov_car = ({_m: sum(1 for _s in _nov_sigs.values() if _m in _s)
                                                       for _m in _pmuts} if _nov_sigs else {})
                                          render_clade_heatmap(
                                              clade_node=f"novel_{_loc}_{_pi}",
                                              shared_mutations=[],
                                              member_blocks=[{
                                                  "member": "novel pattern",
                                                  "discriminating": _pmuts,
                                                  "member_count": 1,
                                                  "reads": _pat.get("count", 0),
                                                  "mut_carriers": _nov_car,
                                              }],
                                              client=wiseLoculus,
                                              location=_loc,
                                              date_range=(start_date, end_date),
                                          )



            if _active_section == "Investigate a variant":
              # ── Investigate a variant (on-demand explorer) ─────────────────────
              st.markdown("---")
              render_variant_explorer(
                  pango_loader=cached_get_pango_loader(),
                  panel=all_selected_variants,
                  disabled=_outstanding,  # avoid rerun races during a scan
              )



if __name__ == "__main__":
    app()