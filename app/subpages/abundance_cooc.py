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
import json
import time
from celery import Celery
import redis
from api.pango_loader import PangoLoader, get_pango_summary_path
import logging
logger = logging.getLogger(__name__)

from api.wiseloculus import WiseLoculusLapis
from utils.config import get_wiseloculus_url, get_cooc_setting

from components.abundance_cooc_tree import render_panel_tree, recombinant_parents, tree_rows
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
                        show_legend: bool = True, panel_union=None) -> None:
    """Normalized (0-100%) per-date composition: explained / + 1 change /
    addable / not attributed (2026-10-06: no novel band — see
    process.completeness_composition). Green height = completeness; other bands = what the gap is made of.
    Empty/low-read dates dropped. Built from existing outputs — scanner untouched."""
    try:
        from process.completeness_composition import compute_completeness_composition
    except Exception:
        return
    rows = compute_completeness_composition(cooc_result, scanner_result or {},
                                            panel_union=panel_union)
    if not rows:
        st.caption("Not enough co-occurrence reads to show composition.")
        return
    import plotly.graph_objects as go
    dates = [r["date"] for r in rows]
    layers = [
        ("explained by panel", "explained_pct", "#0F6E56", "rgba(15,110,86,0.85)"),
        ("panel variant + 1 change", "near_pct", "#4ade80", "rgba(134,239,172,0.85)"),
        ("found, not in panel", "addable_pct", "#dc2626", "rgba(220,38,38,0.55)"),
        ("unexplained, not attributed", "noise_pct", "#9ca3af", "rgba(156,163,175,0.40)"),
    ]
    fig = go.Figure()
    _near_txt = ["<br>".join(f"{_l} ({_n:,} reads)" for _l, _n in r.get("near_top", []))
                 for r in rows]
    for name, _fld, line_c, fill_c in layers:
        _is_near = _fld == "near_pct"
        fig.add_trace(go.Scatter(
            x=dates, y=[r.get(_fld, 0) for r in rows], name=name, mode="lines",
            stackgroup="one", line=dict(width=0.5, color=line_c), fillcolor=fill_c,
            customdata=_near_txt if _is_near else None,
            hovertemplate=("%{x|%Y-%m-%d}<br>" + name + " %{y:.1%}"
                           + ("<br>%{customdata}" if _is_near else "")
                           + "<extra></extra>")))
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


def _panel_verdicts(cooc_res: dict, variants: list) -> dict:
    """Per panel variant in ONE city: {variant: {"state", "reason"}} from the
    worker's panel_check (reads only, never the abundance).

    state: present / absent / mixed (the marker vote of
    process.cooc.check_verdicts), no_marker (nothing specific to test),
    no_data (no check data, or no marker measurable)."""
    from process.cooc import check_verdicts
    from process.variant_explorer import check_in_data
    pc = (cooc_res or {}).get("panel_check", {}) or {}
    _dates = [str(d)[:10] for d in ((cooc_res or {}).get("dates") or [])]
    _RECENT = int(get_cooc_setting("check.recent_samples", default=5))
    sym = {"present": "✓", "absent": "✗", "unmeasured": "?"}
    out = {}
    for v in variants:
        ci = pc.get(v)
        if ci is None:
            out[v] = {"state": "no_data", "reason": "no check data — re-run the scan",
                      "n_present": 0, "n_measured": 0, "n_markers": 0}
            continue
        # the vote uses the last RECENT covered samples ("is it there now?",
        # 2026-10-02); the timeline covers the window
        _cid = check_in_data(ci.get("markers") or [], ci.get("per_date") or {}, _dates,
                             recent=_RECENT)
        _pdr = {d: (ci.get("per_date") or {}).get(d, {}) for d in _cid["recent_dates"]}
        res = check_verdicts(ci.get("markers") or [], _pdr)
        nd, np_, nm = res["n_markers"], res["n_present"], res["n_measured"]
        mk = [f"{k} {(m['freq'] or 0) * 100:.0f}% {sym[m['status']]}"
              if m["cov"] else f"{k} no reads ?" for k, m in res["markers"].items()]
        mtxt = " · ".join(mk[:6]) + (f" · +{len(mk) - 6} more" if len(mk) > 6 else "")
        if nd == 0:
            state, reason = "no_marker", ("every mutation is shared with another panel "
                                          "variant or common outside its family")
        elif nm == 0:
            state, reason = "no_data", f"too few reads on its markers — {mtxt}"
        else:
            state = res["verdict"] if res["verdict"] in ("present", "absent", "mixed") else "no_data"
            reason = f"{np_} of {nm} measurable markers present — {mtxt}"
        out[v] = {"state": state, "reason": reason,
                  "n_present": np_, "n_measured": nm, "n_markers": nd,
                  # one mark per sampling day, for the Evidence calendar
                  "timeline": _cid["timeline"], "recent_dates": _cid["recent_dates"]}
    return out


def _deconv_mean(location_result: dict, location: str, variant: str):
    """Mean deconvolution proportion of a variant in one city, or None."""
    d = location_result
    if isinstance(d, dict) and location in d and isinstance(d[location], dict):
        d = d[location]
    ts = ((d or {}).get(variant) or {}).get("timeseriesSummary", []) or []
    vals = [e.get("proportion", 0) or 0 for e in ts]
    return (sum(vals) / len(vals)) if vals else None


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
                         help="Add the officially-tracked surveillance panel "
                              f"({len(curated_variants)} variants)"):
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
        st.caption("Your panel (structural). Per-city results are coloured on the "
                   "tree in Co-occurrence results.")
        # structural tree only — always black, updates live with the selection
        # (no verdicts here, so it never flickers and needs no run).
        render_panel_tree(
            selected_variants=all_selected_variants,
            yaml_variants=curated_variants,
            pango_loader=cached_get_pango_loader(),
            variant_status={},
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
                format_func=lambda k: f"{k} · {bootstrap_options[k]} bootstraps",
                index=1, key="acooc_bootstrap",
            )
            bootstraps = bootstrap_options[selected_bootstrap]
            bandwidth_options = {"Narrow": 10, "Medium": 20, "Wide": 30}
            selected_bandwidth = st.radio(
                "Bandwidth", options=list(bandwidth_options.keys()),
                format_func=lambda k: f"{k} · bandwidth {bandwidth_options[k]}",
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
            st.session_state["acooc_failed"] = {}
            st.session_state["acooc_scanner_tasks"] = {}
            st.session_state["acooc_scanner_panels"] = {}
            st.session_state["acooc_xcheck_tasks"] = {}
            st.session_state["acooc_xcheck_results"] = {}
            # start times of the two progress phases (for "time left")
            st.session_state["acooc_phase_t0"] = {0: time.time()}
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
            # quick look-up works without a run: ★ markers and relatives come
            # from the pango tree; "Check in data" needs a run's cities
            st.markdown("---")
            render_variant_explorer(pango_loader=cached_get_pango_loader(),
                                    panel=all_selected_variants)
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
            # failed tasks: {stage: {city: message}}. A failed task is "ready" but
            # has no result; without recording it the poller saw it as a fresh
            # result on every tick and reran the whole page every 3 s forever
            # (flickering, clicks lost).
            _failed = st.session_state.setdefault("acooc_failed", {})

            def _fail(stage, loc, err):
                _failed.setdefault(stage, {})[loc] = str(err)[:300]
                logger.error(f"{stage} task failed for {loc}: {err}")

            # deconvolution results
            _loc_res = st.session_state.get("location_results", {})
            for _loc, _tid in list(location_tasks.items()):
                if _loc not in _loc_res and _loc not in _failed.get("deconvolution", {}):
                    _t = celery_app.AsyncResult(_tid)
                    if _t.ready():
                        try:
                            _loc_res[_loc] = _t.get()
                            st.session_state["location_results"] = _loc_res
                        except Exception as _e:
                            _fail("deconvolution", _loc, _e)
                        _new_collected = True

            # scanner results
            scanner_tasks_map = st.session_state.get("acooc_scanner_tasks", {})
            for _loc, _tid in list(scanner_tasks_map.items()):
                if _loc not in scanner_results and _loc not in _failed.get("scanner", {}):
                    _t = celery_app.AsyncResult(_tid)
                    if _t.ready():
                        try:
                            scanner_results[_loc] = _t.get()
                            st.session_state["acooc_scanner_results"] = scanner_results
                            logger.info(f"Scanner results collected for {_loc}")
                        except Exception as _e:
                            _fail("scanner", _loc, _e)
                        _new_collected = True

            # cooc results + auto-submit scanner
            _cooc_res = st.session_state.get("acooc_cooc_results", {})
            _scanner_tasks = st.session_state.get("acooc_scanner_tasks", {})
            try:
                _ref_locs = [str(x) for x in (cached_fetch_locations() or [])]
            except Exception:
                _ref_locs = list(location_tasks)
            for _loc, _tid in list(st.session_state.get("acooc_cooc_tasks", {}).items()):
                if _loc not in _cooc_res and _loc not in _failed.get("completeness", {}):
                    _t = celery_app.AsyncResult(_tid)
                    if _t.ready():
                        try:
                            _cooc_res[_loc] = _t.get()
                            st.session_state["acooc_cooc_results"] = _cooc_res
                        except Exception as _e:
                            _fail("completeness", _loc, _e)
                        _new_collected = True
                # phase 2 — the deep scan — once this city's completeness (phase 1)
                # is ready: the worker reads today's positions plus those where
                # the data shows a mutation and runs the scanner on them
                # (scope.data_positions; the graph stays on phase 1)
                if (_loc in _cooc_res
                        and _loc not in _scanner_tasks
                        and _loc not in scanner_results):
                    _stask = celery_app.send_task(
                        "tasks.run_cooc_deep_scan_lapis",
                        kwargs={
                            "location": _loc,
                            "start_date": start_date.isoformat(),
                            "end_date": end_date.isoformat(),
                            "variants": all_selected_variants,
                            # error hotspots come from every available
                            # location, not only the run's cities
                            "reference_locations": _ref_locs,
                        }
                    )
                    _scanner_tasks[_loc] = _stask.id
                    scanner_panels[_loc] = list(all_selected_variants)
                    logger.info(f"Auto-submitted scanner for {_loc}")
            st.session_state["acooc_scanner_tasks"] = _scanner_tasks
            st.session_state["acooc_scanner_panels"] = scanner_panels

            location_names = list(location_tasks.keys())

            # ── phase 3: cross-check (2026-10-02) ─────────────────────────────
            # A lineage the deep scan named in SOME city is checked by its ★
            # markers in EVERY city (the panel's check). A city's scan names a
            # lineage only from reads with >= 2 mutations beyond the panel, so a
            # young sublineage of a panel variant (PQ.16.1.1 under NB.1.8.1 in
            # Basel: one extra mutation per read) stayed "–" although its marker
            # was on 44 % of the reads.
            _xt = st.session_state.setdefault("acooc_xcheck_tasks", {})
            _xr = st.session_state.setdefault("acooc_xcheck_results", {})
            for _ln, _tid in list(_xt.items()):
                if _ln not in _xr and _ln not in _failed.get("cross-check", {}):
                    _t = celery_app.AsyncResult(_tid)
                    if _t.ready():
                        try:
                            _xr[_ln] = _t.get()
                        except Exception as _e:
                            _fail("cross-check", _ln, _e)
                        _new_collected = True
            _scan_all = st.session_state.get("acooc_scanner_results", {})
            _scans_done = all(_l in _scan_all or _l in _failed.get("scanner", {})
                              for _l in location_names)
            if _scans_done and not _xt and _scan_all:
                _found_in = {}
                for _l, _r in _scan_all.items():
                    for _k in ("resolved_clade", "one_day"):
                        for _c in (_r.get(_k) or []):
                            if _c.get("node"):
                                _found_in.setdefault(_c["node"], set()).add(_l)
                _ran_pan = set(st.session_state.get("acooc_ran_panel") or all_selected_variants)
                _xmax = int(get_cooc_setting("lists.cross_check_max", default=40))
                for _l in location_names:
                    _miss_all = sorted(n for n, cs in _found_in.items()
                                       if _l not in cs and n not in _ran_pan)
                    _miss = _miss_all[:_xmax]
                    if len(_miss_all) > _xmax:
                        logger.warning(f"Cross-check {_l}: checking {_xmax} of "
                                       f"{len(_miss_all)} lineages (lists.cross_check_max)")
                    if not _miss or _l in _failed.get("scanner", {}):
                        _xr[_l] = {"dates": [], "per_variant": {}, "found_in": {}}
                        _xt[_l] = ""
                        continue
                    _xtask = celery_app.send_task("tasks.run_cooc_lineages_check_lapis", kwargs={
                        "location": _l, "start_date": start_date.isoformat(),
                        "end_date": end_date.isoformat(), "variants": _miss,
                        "panel": sorted(_ran_pan)})
                    _xt[_l] = _xtask.id
                st.session_state["acooc_xcheck_found_in"] = {n: sorted(cs) for n, cs in _found_in.items()}
                logger.info(f"Cross-check submitted: {sum(1 for t in _xt.values() if t)} cities")

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
                if ((_ln not in _lr_now and _ln not in _failed.get("deconvolution", {}))
                        or (_ln not in _cr_now and _ln not in _failed.get("completeness", {}))):
                    _outstanding = True
                    break
                # completeness done but scanner not yet done → running
                if (_ln in _cr_now
                        and _ln not in _sr_now and _ln not in _failed.get("scanner", {})):
                    _outstanding = True
                    break
            # phase 3 (cross-check) not yet in → running
            if not _outstanding and _sr_now:
                _xt_now = st.session_state.get("acooc_xcheck_tasks", {})
                _xr_now = st.session_state.get("acooc_xcheck_results", {})
                if not _xt_now or any(_l not in _xr_now and _l not in _failed.get("cross-check", {})
                                      for _l in location_names):
                    _outstanding = True
            if not _outstanding and _new_collected:
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

            # Two phases (2026-10-02). Phase 1 = deconvolution + completeness:
            # the graph and abundances. Phase 2 = the deep scan: Variants table
            # and the red band. Each bar: % of the work (from the
            # worker's own progress, per city) and the time left at the pace so
            # far. Drawn inside a fragment that reruns every 3 s while work is
            # outstanding, so the bars move without redrawing the charts; the
            # whole page reruns only when a task has finished.
            _deep_on = get_cooc_setting("scope.data_positions", default=True)
            _PHASES = [
                ("1 · Abundance and completeness graph",
                 [("deconvolution", "acooc_location_tasks", "location_results"),
                  ("completeness", "acooc_cooc_tasks", "acooc_cooc_results")], "#475569"),
                ("2 · " + ("Deep scan — panel + positions with a mutation in the data"
                           if _deep_on else "Scanner — panel positions only"),
                 [("scanner", "acooc_scanner_tasks", "acooc_scanner_results")], "#475569"),
                ("3 · Cross-check — lineages found in any city, checked in every city",
                 [("cross-check", "acooc_xcheck_tasks", "acooc_xcheck_results")], "#475569"),
            ]

            def _task_frac(tid):
                """0-1 for a running task, from task_progress:{id} in redis."""
                try:
                    _raw = redis_client.get(f"task_progress:{tid}")
                    if not _raw:
                        return 0.0, ""
                    _d = json.loads(_raw)
                    _tot = max(1, int(_d.get("total") or 4))
                    _cur = max(1, int(_d.get("current") or 1))
                    _f = (_cur - 1 + float(_d.get("frac") or 0)) / _tot
                    return min(0.99, max(0.0, _f)), str(_d.get("status") or "")
                except Exception:
                    return 0.0, ""

            def _left(sec):
                if sec < 60:
                    return "< 1 min left"
                return f"~{int(round(sec / 60))} min left"

            def _phase_bars(live):
                import time as _time
                _ss = st.session_state
                _fl = _ss.get("acooc_failed", {})
                _t0s = _ss.setdefault("acooc_phase_t0", {})
                _ready = False
                _html = ""
                for _pi, (_title, _stages, _col) in enumerate(_PHASES):
                    _fs, _now_txt, _n_bad = [], "", 0
                    for _stg, _tk, _rk in _stages:
                        _tasks, _done = _ss.get(_tk, {}), _ss.get(_rk, {})
                        _bad = _fl.get(_stg, {})
                        _n_bad += len(_bad)
                        for _l in location_names:
                            if _l in _done or _l in _bad:
                                _fs.append(1.0)
                            elif _l in _tasks:
                                if live and celery_app.AsyncResult(_tasks[_l]).ready():
                                    _ready = True
                                _f, _msg = _task_frac(_tasks[_l]) if live else (0.0, "")
                                _fs.append(_f)
                                if _msg and not _now_txt:
                                    _now_txt = f"{_l.split('(')[0].strip()}: {_msg}"
                            else:
                                _fs.append(0.0)
                    _frac = sum(_fs) / len(_fs) if _fs else 0.0
                    _started = any(f > 0 for f in _fs) or any(
                        _ss.get(_tk) for _, _tk, _ in _stages)
                    if _started and _pi not in _t0s:
                        _t0s[_pi] = _time.time()
                    _right = f"{_frac * 100:.0f}%"
                    if _frac >= 0.999:
                        _right = "done"
                    elif live and _started and _frac > 0.03 and _pi in _t0s:
                        _el = _time.time() - _t0s[_pi]
                        _right += " · " + _left(_el * (1 - _frac) / _frac)
                    elif not _started:
                        _right = ("waiting for phase 1", "waiting for the deep scan")[_pi - 1] if _pi else "—"
                    if _n_bad:
                        _right += f" · {_n_bad} failed"
                    _sub = (f"<div style='font-size:10.5px;color:#898781;margin-top:2px;"
                            f"white-space:nowrap;overflow:hidden;text-overflow:ellipsis;'>"
                            f"{_now_txt}</div>" if (live and _now_txt and _frac < 0.999) else "")
                    _html += (
                        f"<div style='margin-bottom:8px;'>"
                        f"<div style='display:flex;justify-content:space-between;font-size:11px;"
                        f"margin-bottom:3px;'><span style='font-weight:500;'>{_title}</span>"
                        f"<span style='color:#898781;'>{_right}</span></div>"
                        f"<div style='height:8px;background:#F1EFE8;border-radius:4px;overflow:hidden;'>"
                        f"<div style='height:100%;border-radius:4px;width:{_frac*100:.0f}%;"
                        f"background:{_col};'></div></div>{_sub}</div>")
                st.markdown(_html, unsafe_allow_html=True)
                return _ready

            if _outstanding:
                @st.fragment(run_every=3)
                def _poll_tasks():
                    if _phase_bars(live=True):
                        st.rerun(scope="app")
                _poll_tasks()
            else:
                _phase_bars(live=False)
            _fl_all = st.session_state.get("acooc_failed", {})
            if any(_fl_all.values()):
                st.error("Failed: " + " · ".join(
                    f"{_stg} in {', '.join(c.split('(')[0].strip() for c in _cs)}"
                    for _stg, _cs in _fl_all.items() if _cs)
                    + " — the worker may have restarted or run out of memory; re-run to retry.")

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
            _sections = ["Deconvolution results", "Co-occurrence results", "Signal over time",
                         "Investigate a variant"]
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
              # City radio removed — the deconv plots for all cities are shown
              # together in a small-multiples grid below.

              _selected = st.session_state.get("acooc_selected_city", _city_options[0] if _city_options else "")
              st.markdown("<hr style='margin:6px 0 10px;opacity:.15;'>", unsafe_allow_html=True)

              def _city_tab_content(location, task_id, show_deconv=True,
                                    show_similarity=True):
                  """One city's deconvolution plot (the co-occurrence check and
                  per-city tree moved to Co-occurrence results)."""

                  # ── deconvolution (primary output) ────────────────────────────
                  # In the multi-city grid the deconv plots are drawn once by
                  # render_location_grid; here only when show_deconv is set.
                  if show_deconv:
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

                  # ── Jaccard (signature similarity) — city-independent, so it's
                  #    shown once above the tabs (show_similarity=False here). ──
                  if show_similarity and len(all_selected_variants) >= 2:
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

              # ── signature-similarity matrix: city-independent, shown ONCE
              #    above the tabs (not repeated per city). ──
              if len(all_selected_variants) >= 2:
                  st.markdown("<div style='font-weight:600;font-size:13px;'>"
                              "🧬 Signature similarity (Jaccard)</div>",
                              unsafe_allow_html=True)
                  st.caption("How much each pair of panel variants shares mutations "
                             "— high similarity means deconvolution may struggle to "
                             "tell them apart. (Same for every city.)")
                  with st.expander("Show similarity heatmap", expanded=False):
                      render_jaccard_heatmap(
                          variants=all_selected_variants,
                          pango_loader=cached_get_pango_loader(),
                      )
              # ── per-city TABS: each tab is one city's full-size deconvolution
              #    plot (with confidence bands). ──
              _dtabs = st.tabs([loc for loc in location_names])
              for _dtab, _loc in zip(_dtabs, location_names):
                  with _dtab:
                      if _loc in location_tasks:
                          _city_tab_content(_loc, location_tasks[_loc],
                                            show_deconv=True, show_similarity=False)
                      else:
                          st.caption("Not started.")


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
              # phase 1 is enough for the graph and the tree (panel variants and
              # their +1 changes); the deep scan adds findings when it ends
              _ready_locs = [l for l in location_names if _cr_all.get(l) is not None]
              # last sampling date with data per city — the window may run past it
              _win_end = str((st.session_state.get("acooc_ran_dates")
                              or (None, end_date.isoformat()))[1])[:10]
              _data_until = {l: max(str(d)[:10] for d in (_cr_all[l].get("dates") or []))
                             for l in location_names
                             if _cr_all.get(l) and _cr_all[l].get("dates")}
              _until_short = {l: d for l, d in _data_until.items() if d < _win_end}

              def _dshort(d):
                  try:
                      from datetime import date as _dt
                      _x = _dt.fromisoformat(d)
                      return f"{_x.day} {_x.strftime('%b')}"
                  except Exception:
                      return d
              # the graph needs only phase 1 (completeness); the deep scan colours
              # its red band in when it arrives
              _graph_locs = [l for l in location_names if _cr_all.get(l) is not None]
              if _graph_locs:
                  st.markdown("#### Panel completeness")
                  if _until_short:
                      st.caption("📅 Data until " + " · ".join(
                          f"{l.split('(')[0].strip()} {_dshort(d)}"
                          for l, d in sorted(_until_short.items()))
                          + f" — the window runs to {_dshort(_win_end)}; later dates have "
                          "no samples yet.")
                  st.caption("How much of each city's co-occurrence signal your panel "
                             "explains (green) vs the rest. Light green = a panel variant with one "
                             "change (hover for which) · red = a lineage the deep scan named · "
                             "grey = unexplained, not attributed (novel patterns: Variants table).")
                  # One shared legend ABOVE the grid; every plot has
                  # show_legend=False and the same fixed height, so all plot areas
                  # are identical. (Previously the legend went on the first plot
                  # only and ate into its fixed height, making it shorter.)
                  _comp_leg = [
                      ("#0F6E56", "explained by panel"),
                      ("#4ade80", "panel variant + 1 change"),
                      ("#dc2626", "found, not in panel"),
                      ("#9ca3af", "unexplained, not attributed"),
                  ]
                  # panel union (city-independent) for the near-panel split —
                  # a read that is a panel variant + <2 stray mutations counts
                  # toward completeness, not the grey gap.
                  _pl_pu = cached_get_pango_loader()
                  _panel_union = {}   # {variant: signature}
                  # use the panel that was RUN (frozen with the cached results),
                  # not the live selection — so editing the panel without
                  # re-running doesn't silently change the completeness graph.
                  for _puv in (st.session_state.get("acooc_ran_panel")
                               or all_selected_variants):
                      for _pum in (_pl_pu.get_signature(_puv) or []):
                          if _pum and _pum[-1] in "ACGT":
                              _panel_union.setdefault(_puv, set()).add(_pum)
                  st.markdown(
                      "<div style='display:flex;gap:14px;flex-wrap:wrap;"
                      "font-size:11.5px;color:#6b7280;margin:2px 0 8px;'>"
                      + "".join(
                          f"<span><span style='display:inline-block;width:11px;"
                          f"height:11px;border-radius:2px;background:{_c};"
                          f"vertical-align:-1px;margin-right:5px;'></span>{_t}</span>"
                          for _c, _t in _comp_leg)
                      + "</div>", unsafe_allow_html=True)
                  _grid_locs = list(location_names)
                  def _one_city(_lc):
                      _du = (f" <span style='font-weight:400;color:#6b7280;'>· data until "
                             f"{_dshort(_until_short[_lc])}</span>" if _lc in _until_short else "")
                      if _cr_all.get(_lc) is not None and _sr_all.get(_lc) is None:
                          _du += (" <span style='font-weight:400;color:#6b7280;'>· deep scan "
                                  "running — red appears when it ends</span>")
                      st.markdown(f"<div style='font-size:12px;font-weight:600;'>"
                                  f"{_lc}{_du}</div>", unsafe_allow_html=True)
                      if _cr_all.get(_lc) is not None:
                          _render_composition(_cr_all[_lc], _sr_all.get(_lc) or {}, key=_lc,
                                              show_legend=False,
                                              panel_union=_panel_union)
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

              # ── Variants: one city selector; the tree (left) and one table for
              #    everything co-occurrence says (right); then the chosen city's
              #    signal over time. Uses the panel that was RUN, so the colours
              #    always match the results. (2026-09-30) ─────────────────────
              _run_panel = st.session_state.get("acooc_ran_panel") or all_selected_variants
              _verdicts_by_city = {
                  _l: _panel_verdicts(_cr_all[_l], _run_panel)
                  for _l in location_names if _cr_all.get(_l) is not None}
              # scanner findings only when NOTHING is still running
              # (_outstanding drives the autorefresh), so partial results
              # never leak into the table
              _scan_running = _outstanding
              _scan_res_all = ({} if _scan_running
                               else st.session_state.get("acooc_scanner_results", {}) or {})
              _agg_new, _agg_sub, _agg_unres, _agg_mnh, _agg_one = {}, {}, {}, {}, {}
              _novel_total, _novel_pats = 0, 0

              def _human_reads(_n):
                  if _n >= 1_000_000:
                      return f"{_n/1_000_000:.1f}M".replace(".0M", "M")
                  if _n >= 1_000:
                      return f"{_n/1_000:.0f}K"
                  return str(_n)

              if _scan_res_all:
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
                  # collect clade findings across cities, keyed by node
                  _agg_new = {}      # node -> {reads, cities, member_count, members, designation, muts}
                  _agg_sub = {}
                  _agg_unres = {}    # frozenset(fp) -> {reads, cand, anc}
                  _agg_mnh = {}      # node -> matched-but-no-haplotype
                  _agg_one = {}      # node -> ★ evidence on 1 day only, in every city
                  _novel_total = 0
                  _novel_pats = 0

                  from process.variant_explorer import check_in_data as _cid

                  def _fcheck(_loc, _node):
                      """The panel check for a found lineage in one city (the deep
                      scan's findings_check): same measure as the panel cells."""
                      _fc = (_scan_res_all.get(_loc) or {}).get("findings_check") or {}
                      _pv = (_fc.get("per_variant") or {}).get(_node)
                      if not _pv:
                          return None
                      return _cid(_pv.get("markers") or [], _pv.get("per_date") or {},
                                  _fc.get("dates") or [],
                                  recent=int(get_cooc_setting("check.recent_samples", default=5)))

                  def _xcheck(_loc, _node):
                      """Phase 3: the ★ marker check of a lineage named in another
                      city, for this city (None when it wasn't checked here)."""
                      _x = (st.session_state.get("acooc_xcheck_results", {}) or {}).get(_loc) or {}
                      _pv = (_x.get("per_variant") or {}).get(_node)
                      if not _pv or not _pv.get("markers"):
                          return None
                      return _cid(_pv.get("markers") or [], _pv.get("per_date") or {},
                                  _x.get("dates") or [],
                                  recent=int(get_cooc_setting("check.recent_samples", default=5)))

                  def _ingest_clade(_c, _loc, _bucket=None):
                      # Route one finding (top-level OR a promoted sub-finding) into
                      # the new/sub bucket and aggregate its reads across cities.
                      if _bucket is None:
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
                      # per-city ★ mutations (the scanner's mut_star). Kept per
                      # city so the card shows which city has which. Regex-free
                      # leading-digit position sort.
                      def _pcpos(_m):
                          _d = ""
                          for _ch in _m:
                              if _ch.isdigit():
                                  _d += _ch
                              else:
                                  break
                          return int(_d) if _d else 0
                      _pc_car, _pc_flag, _pc_out = {}, {}, {}
                      for _blk in _c.get("member_blocks", []):
                          _pc_car.update(_blk.get("mut_carriers", {}))
                          _pc_flag.update(_blk.get("mut_star", {}) or {})
                          _pc_out.update(_blk.get("mut_outside", {}) or {})
                      # ★ = rare outside the finding's own family, decided by the
                      # scanner (mut_star)
                      _pc_star = sorted(
                          {_m for _blk in _c.get("member_blocks", [])
                           for _m in _blk.get("discriminating", [])
                           if _pc_flag.get(_m, False)},
                          key=_pcpos)
                      _slot["per_city"][_loc] = {
                          "star": _pc_star,
                          "car": {_m: _pc_out.get(_m, _pc_car.get(_m)) for _m in _pc_star},
                          # days with ★ evidence in this city (shown on the chip)
                          "days": len(_c.get("counted_days", []) or []),
                          # facts for the table (scanner _add_notes)
                          "evidence_days": dict(_c.get("evidence_days", {}) or {}),
                          "regions": list(_c.get("evidence_regions", []) or []),
                          # one-day findings: what happened after that day
                          "after": _c.get("one_day_after"),
                          # the panel check on its ★ markers (2026-10-02)
                          "check": _fcheck(_loc, _node),
                          # this city's own blocks, for this city's heatmap
                          "blocks": _c.get("member_blocks", []) or [],
                      }
                      _slot.setdefault("descendants", set()).update(
                          _c.get("confirmed_descendants", []) or [])

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
                          _s.setdefault("days", {})[_loc] = len(_u.get("days", []) or [])
                      for _mnh in _res.get("matched_no_haplotype", []):
                          _k = _mnh["node"]
                          _s = _agg_mnh.setdefault(_k, {
                              "node": _mnh["node"], "reads": 0,
                              "member_count": _mnh.get("member_count", 1),
                              "muts": _mnh.get("observed_mutations", []),
                              "relationship": _mnh.get("relationship", ""),
                              "cities": [],
                          })
                          _s["reads"] += int(_mnh.get("total_reads", 0))
                          if _loc not in _s["cities"]:
                              _s["cities"].append(_loc)
                      _nv = _res.get("novel", {})
                      _novel_total += int(_nv.get("total_reads", 0))
                      _novel_pats += int(_nv.get("pattern_count", 0))
                  # ★ evidence on 1 day only: a node confirmed in another city
                  # shows this city as a "1 day" chip on its card; otherwise it
                  # goes to "Seen on 1 day only" (a possible jackpot — watch it)
                  for _loc, _res in _scan_res_all.items():
                      for _c in _res.get("one_day", []) or []:
                          _b = (_agg_new if _c["node"] in _agg_new
                                else _agg_sub if _c["node"] in _agg_sub else _agg_one)
                          _ingest_clade(_c, _loc, _b)

              if _ready_locs:
                  from components import variants_table as _vt
                  st.markdown("---")
                  st.markdown("### Variants")
                  st.caption("Evidence from reads only, one column per city. The tree is "
                             "coloured by the city you click in the header (shaded). Hover "
                             "anything for details; + Add puts a finding in your panel — "
                             "re-run to apply.")
                  # the chosen city (click a city code in the table header)
                  if st.session_state.get("acooc_tree_city") not in _ready_locs:
                      st.session_state["acooc_tree_city"] = _ready_locs[0]
                  _tcity = st.session_state["acooc_tree_city"]
                  _cities_all = [l for l in location_names
                                 if l in _verdicts_by_city or l in _scan_res_all]

                  # ---- findings not in the panel ----
                  _conf = sorted(list(_agg_new.values()) + list(_agg_sub.values()),
                                 key=lambda x: -x["reads"])
                  _broader = [x for x in _conf if x.get("descendants")]
                  _conf = [x for x in _conf if not x.get("descendants")]
                  _ED = int(get_cooc_setting("evidence.min_days", default=2))
                  _conf_ok = [x for x in _conf if any(
                      (p.get("days") or 0) >= _ED for p in x["per_city"].values())]
                  _one = ([x for x in _conf if x not in _conf_ok]
                          + sorted(_agg_one.values(), key=lambda x: -x["reads"]))
                  _shown = {x["node"] for x in _conf + _broader + _one}
                  # a node confirmed in one city but only named in another is
                  # one row (its city cells say which)
                  _named_all = [x for x in sorted(_agg_mnh.values(), key=lambda x: -x["reads"])
                                if x["node"] not in _shown]
                  _named = _named_all[:int(get_cooc_setting("lists.named_only_max", default=25))]

                  # findings as the view wants them: node -> {status, per_city, addable}
                  def _f(_slot, _status):
                      if _status == "named":
                          _pc = {c: {"days": 0} for c in _slot.get("cities", [])}
                      else:
                          _pc = {c: {"days": p.get("days", 0), "stars": p.get("star", []),
                                     "regions": p.get("regions", []), "after": p.get("after"),
                                     "check": p.get("check")}
                                 for c, p in (_slot.get("per_city") or {}).items()}
                          # phase 3: the other cities, by the cross-check
                          _fin = (st.session_state.get("acooc_xcheck_found_in") or {}).get(
                              _slot.get("node"), sorted(_pc))
                          for _xc in location_names:
                              if _xc in _pc:
                                  continue
                              _chk = _xcheck(_xc, _slot.get("node"))
                              if _chk is not None:
                                  _pc[_xc] = {"days": 0, "check": _chk, "xcheck": True,
                                              "found_in": [c for c in _fin if c != _xc]}
                      _ok = any((d.get("days") or 0) >= _ED for d in _pc.values())
                      # + Add (2026-10-02): the panel check would confirm it in at
                      # least one city, and it is more than one sample (>= min_days
                      # days in a city, or 1 day in >= 2 cities)
                      _present = any((d.get("check") or {}).get("state") == "present"
                                     for d in _pc.values())
                      _samples = _ok or sum(1 for d in _pc.values()
                                            if (d.get("days") or 0) >= 1) >= 2
                      return {"status": _status if _status != "confirmed" or _ok else "1 day",
                              "per_city": _pc,
                              "addable": _status != "named" and _present and _samples}

                  _findings = {}
                  if not _scan_running:
                      for x in _named:
                          _findings[x["node"]] = _f(x, "named")
                      for x in _one:
                          _findings[x["node"]] = _f(x, "1 day")
                      for x in _conf_ok + _broader:
                          _findings[x["node"]] = _f(x, "confirmed")

                  # ---- unnamed signal: novel patterns grouped by mutations ----
                  _nov = {}   # mutations -> {city: days}
                  _nov_info = {}   # mutations -> {city: {timeline, clue}}
                  # error-hotspot patterns the scanner left out (_drop_hotspots)
                  _hot_p = sum(int(((_r.get("novel") or {}).get("hotspot") or {}).get("patterns", 0))
                               for _r in _scan_res_all.values())
                  _hot_r = sum(int(((_r.get("novel") or {}).get("hotspot") or {}).get("reads", 0))
                               for _r in _scan_res_all.values())
                  for _loc, _res in _scan_res_all.items():
                      for _g in (_res.get("novel", {}) or {}).get("groups", []) or []:
                          if _g.get("likely_error"):      # results from before 2026-10-02
                              continue
                          _k = tuple(_g["mutations"])
                          _nov.setdefault(_k, {})[_loc] = len(_g.get("days", []))
                          _nov_info.setdefault(_k, {})[_loc] = {
                              "timeline": _g.get("timeline"), "clue": _g.get("clue")}
                  _nov_order = lambda kv: (-max(kv[1].values(), default=0),
                                           -sum(kv[1].values()), kv[0])
                  _nov_list = sorted([(list(k), v) for k, v in _nov.items()], key=_nov_order)
                  _nov_listed = sum(
                      _g.get("reads", 0) for _res in _scan_res_all.values()
                      for _g in (_res.get("novel", {}) or {}).get("groups", []) or [])
                  _nov_rest = max(0, _novel_total - _nov_listed)
                  _broad = [(u["fp"], u.get("days", {}), u["cand"], u["anc"])
                            for u in sorted(_agg_unres.values(), key=lambda x: -x["reads"])[
                                :int(get_cooc_setting("lists.broad_max", default=25))]]

                  if _scan_running:
                      st.markdown(
                          "<style>@keyframes acoocspin{to{transform:rotate(360deg)}}"
                          "@keyframes acoocpulse{50%{box-shadow:0 0 0 4px rgba(71,85,105,.15)}}</style>"
                          "<div style='display:flex;align-items:center;gap:10px;padding:9px 14px;"
                          "margin:4px 0 10px;border-radius:8px;background:#f8fafc;"
                          "border:1px solid #cbd5e1;border-left:4px solid #475569;"
                          "animation:acoocpulse 2s ease-in-out infinite;'>"
                          "<span style='width:14px;height:14px;border-radius:50%;flex:none;"
                          "border:2px solid #cbd5e1;border-top-color:#475569;"
                          "animation:acoocspin .9s linear infinite;'></span>"
                          "<span style='font-size:13px;color:#1f2937;'><b>Deep scan running</b> — "
                          "lineages not in your panel and novel patterns appear here when it ends."
                          "</span></div>", unsafe_allow_html=True)
                  if "acooc_recomb_parents" not in st.session_state:
                      st.session_state["acooc_recomb_parents"] = recombinant_parents()
                  _rows = tree_rows(_run_panel, curated_variants, cached_get_pango_loader(),
                                    findings=list(_findings),
                                    recomb_parents=st.session_state["acooc_recomb_parents"])
                  # ── panel variant ± 1 change (the amber band, per change) ──
                  from process.near_changes import near_changes, where_in_tree
                  from process.scanner import (EVIDENCE_MIN_DAYS, EVIDENCE_MIN_READS,
                                               EVIDENCE_MIN_SHARE)
                  _near = near_changes({c: _cr_all[c] for c in _cities_all if _cr_all.get(c)},
                                       _run_panel, EVIDENCE_MIN_READS, EVIDENCE_MIN_SHARE)
                  if _near:
                      _pl_raw = cached_get_pango_loader().get_raw_data()
                      if "acooc_pango_children" not in st.session_state:
                          _kids = {}
                          for _l, _e in _pl_raw.items():
                              if _e.get("parent"):
                                  _kids.setdefault(_e["parent"], []).append(_l)
                          st.session_state["acooc_pango_children"] = _kids
                      _sigs_n = st.session_state.get("acooc_all_sigs_cache")
                      if _sigs_n is None:
                          _ld = cached_get_pango_loader()
                          _sigs_n = {l: _ld.get_signature(l) for l in _pl_raw}
                          st.session_state["acooc_all_sigs_cache"] = _sigs_n
                      for _v, _chs in _near.items():
                          for _ch in _chs:
                              _ch["where"] = where_in_tree(
                                  _v, _ch["sign"], _ch["mut"], _sigs_n,
                                  st.session_state["acooc_pango_children"])
                  _view = _vt.build(_cities_all, _tcity, _rows, _verdicts_by_city, _findings,
                                    current_panel=set(all_selected_variants),
                                    ot=curated_variants, novel=_nov_list, broad=_broad,
                                    novel_hotspot=((_hot_p, _human_reads(_hot_r))
                                                   if _hot_p else None),
                                    novel_info=_nov_info,
                                    novel_rest=_human_reads(_nov_rest) if _nov_rest else None,
                                    broad_total=len(_agg_unres),
                                    named_hidden=(0 if _scan_running
                                                  else len(_named_all) - len(_named)),
                                    near=_near, near_min_days=EVIDENCE_MIN_DAYS,
                                    data_until={l: _dshort(d) for l, d in _until_short.items()})
                  _click = _vt.render(_view, key="acooc_variants_view")
                  if _click and _click.get("t") != st.session_state.get("acooc_vv_last_click"):
                      st.session_state["acooc_vv_last_click"] = _click.get("t")
                      if _click.get("add"):
                          st.session_state[f"acooc_add_variant_pending_{_click['add']}"] = _click["add"]
                      elif _click.get("city") in _ready_locs:
                          st.session_state["acooc_tree_city"] = _click["city"]
                      st.rerun()

                  st.caption("📈 Per-week heatmaps of any lineage or novel pattern: "
                             "**Signal over time**.")


            if _active_section == "Investigate a variant":
              # ── Investigate a variant (on-demand explorer) ─────────────────────
              st.markdown("---")
              # "Check in data" uses the run's cities, window and panel, so its
              # answer matches the results above
              _rd = st.session_state.get("acooc_ran_dates") or (start_date.isoformat(),
                                                                end_date.isoformat())
              render_variant_explorer(
                  pango_loader=cached_get_pango_loader(),
                  panel=st.session_state.get("acooc_ran_panel") or all_selected_variants,
                  # the look-up always works; "Check in data" waits for the run
                  cities=None if _outstanding else list(location_names), start_date=_rd[0], end_date=_rd[1],
                  celery_app=celery_app,
              )

            if _active_section == "Signal over time":
              # ── per-week heatmaps (2026-10-02; were at the bottom of the
              #    co-occurrence results) ─────────────────────────────────────────
              st.markdown("---")
              from components.signal_over_time_ui import render_signal_over_time
              _rd2 = st.session_state.get("acooc_ran_dates") or (start_date.isoformat(),
                                                                 end_date.isoformat())
              render_signal_over_time(
                  pango_loader=cached_get_pango_loader(), client=wiseLoculus,
                  cities=list(location_names), start_date=_rd2[0], end_date=_rd2[1],
                  panel=st.session_state.get("acooc_ran_panel") or all_selected_variants,
                  # the scan's groups and evidence days once it is done; before
                  # that, every lineage uses its ★ marker groups
                  scanner_results=({} if _outstanding else
                                   st.session_state.get("acooc_scanner_results", {}) or {}),
                  default_city=st.session_state.get("acooc_tree_city"),
                  scanning=_outstanding,
              )



if __name__ == "__main__":
    app()