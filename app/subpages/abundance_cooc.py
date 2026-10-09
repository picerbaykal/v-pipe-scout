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


def _render_covvfit(r: dict) -> None:
    """CovvFit's result (2026-10-09): the figure, the growth advantages between
    panel variants (each pair once, the faster one first), the rows against
    covvfit's "other" apart, and the downloads."""
    import base64 as _b64
    import math as _math
    import pandas as _pd
    st.caption(
        f"Fitted on each sample's own estimate (deconvolution without smoothing), "
        f"{r.get('date_min', '')} – {r.get('date_max', '')}, predicted {r.get('horizon', '')} "
        f"days ahead (grey area). One growth advantage per variant, shared by "
        f"{len(r.get('locations') or [])} location(s).")
    for _l, _why in (r.get("skipped") or {}).items():
        st.warning(f"Left out: {_l} ({_why})")
    if r.get("figure_png"):
        st.image(_b64.b64decode(r["figure_png"]), width="stretch")

    rows = r.get("pairwise") or []
    main = [x for x in rows if not x.get("involves_other") and x["estimate"] > 0]
    if main:
        def _reading(x):
            return ("grows faster" if x["lower"] > 0 else "no clear difference")
        df = _pd.DataFrame([{
            "Variant": x["variant"], "against": x["reference"],
            "advantage per week": round(x["estimate"], 3),
            "95 % range": f"{x['lower']:.3f} – {x['upper']:.3f}",
            "odds × per week": round(_math.exp(x["estimate"]), 2),
            "reading": _reading(x)} for x in sorted(main, key=lambda x: -x["estimate"])])
        st.markdown("<div style='font-size:13px;font-weight:600;margin-top:6px;'>"
                    "Between your panel variants</div>", unsafe_allow_html=True)
        st.caption("Advantage per week on the log-odds scale: 0.5 means the odds of the "
                   "first variant against the second multiply by e^0.5 ≈ 1.65 each week. "
                   "\"No clear difference\" when the 95 % range includes 0.")
        st.dataframe(df, width="stretch", hide_index=True)
    oth = [x for x in rows if x.get("involves_other") and x["variant"] != "other"]
    if oth:
        with st.expander("Against \"other\" (the share your panel doesn't explain)",
                         expanded=False):
            st.caption("covvfit's \"other\" is 1 minus your panel variants: LolliPop's "
                       "undetermined. When the panel explains almost everything, it is "
                       "near 0 and these advantages can't be estimated (very wide ranges).")
            st.dataframe(_pd.DataFrame([{
                "Variant": x["variant"], "advantage per week vs other": round(x["estimate"], 3),
                "95 % range": f"{x['lower']:.3f} – {x['upper']:.3f}"} for x in oth]),
                width="stretch", hide_index=True)

    c1, c2, c3, c4 = st.columns(4)
    with c1:
        if r.get("figure_png"):
            st.download_button("⬇ Figure (PNG)", _b64.b64decode(r["figure_png"]),
                               file_name="covvfit_figure.png", mime="image/png",
                               key="acooc_cvf_png", width="stretch")
    with c2:
        if r.get("figure_pdf"):
            st.download_button("⬇ Figure (PDF)", _b64.b64decode(r["figure_pdf"]),
                               file_name="covvfit_figure.pdf", mime="application/pdf",
                               key="acooc_cvf_pdf", width="stretch")
    with c3:
        if r.get("pairwise_csv"):
            st.download_button("⬇ Advantages (TSV)", r["pairwise_csv"],
                               file_name="covvfit_pairwise_fitnesses.tsv",
                               mime="text/tab-separated-values",
                               key="acooc_cvf_pw", width="stretch")
    with c4:
        if r.get("predictions_csv"):
            st.download_button("⬇ Predictions (TSV)", r["predictions_csv"],
                               file_name="covvfit_predictions.tsv",
                               mime="text/tab-separated-values",
                               key="acooc_cvf_pred", width="stretch")


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
    # the kinds of task a run can start: (stage, tasks key, results key).
    # Stage names are the keys of acooc_failed; "stopped" there = Stop pressed.
    # Each is {location: task id}, except CovvFit: one task for all locations,
    # under the key _ALL (2026-10-09).
    _RUN_STAGES = [("completeness", "acooc_cooc_tasks", "acooc_cooc_results"),
                   ("deconvolution", "acooc_location_tasks", "location_results"),
                   ("scanner", "acooc_scanner_tasks", "acooc_scanner_results"),
                   ("cross-check", "acooc_xcheck_tasks", "acooc_xcheck_results"),
                   ("covvfit", "acooc_covvfit_tasks", "acooc_covvfit_results")]
    _STOPPED = "stopped"
    _ALL = "all locations"

    # ── Header ───────────────────────────────────────────────────────────────
    st.title("Abundance & Co-occurrence")
    st.subheader(
        "Estimate the proportion of variants circulating in wastewater over time, "
        "using live LAPIS mutation data and LolliPop deconvolution with bootstrap "
        "confidence intervals."
    )
    st.caption(
        "Build a custom variant panel and run on-demand deconvolution across one or "
        "more sampling locations. The panel check and the deep scan guide you "
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
    # a run = the completeness tasks of ▶ Run (2026-10-08: deconvolution and
    # the deep scan have their own buttons)
    has_run = bool(st.session_state.get("acooc_cooc_tasks"))
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
                   "tree in the Lineages tab.")
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
            _ss = st.session_state
            _fl = _ss.get("acooc_failed", {})
            _checks = [(_stg, _ss.get(_tk, {}), _ss.get(_rk, {}))
                       for _stg, _tk, _rk in _RUN_STAGES]
            for _stg, _tasks, _results in _checks:
                for _loc, _t in _tasks.items():
                    if _loc in _results or _loc in _fl.get(_stg, {}):
                        continue  # collected, failed or stopped — not running
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
            # ▶ Run = completeness only (deconvolution has its own button;
            # it took ~1.5 min per location)
            _total_s   = _cooc_s
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
            help="Runs the panel check for all locations. Deconvolution and "
                 "the deep scan have their own buttons above the results."
                 if not _busy else "Wait for the current analysis to finish, or stop it.",
        ):
            st.session_state["acooc_trigger_run"] = True
        st.caption("panel check for all locations" if not _busy
                   else "⟳ analysis in progress… (Stop is above the results)")
    # ── Right column: all outputs ─────────────────────────────────────────────
    with col_results:

        # task submission
        if st.session_state.get("acooc_trigger_run"):
            st.session_state["acooc_trigger_run"] = False
            # ▶ Run = completeness only (2026-10-08). Deconvolution, the deep
            # scan (+ cross-check) and CovvFit start from their buttons, with
            # the cities, panel and dates saved here (acooc_ran_*).
            cooc_tasks = {}
            for loc in selected_locations:
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

            st.session_state["acooc_location_tasks"] = {}
            st.session_state["acooc_deconv_settings"] = None
            st.session_state["acooc_want_scan"] = False
            st.session_state["acooc_stopped"] = False
            st.session_state["acooc_section"] = "Panel check"
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
            st.session_state["acooc_covvfit_tasks"] = {}
            st.session_state["acooc_covvfit_results"] = {}
            # start times of the two progress phases (for "time left")
            st.session_state["acooc_phase_t0"] = {"completeness": time.time()}
            # reset bucket expand flags so scanner starts collapsed on a new run
            st.session_state["acooc_exp_missing"] = False
            st.session_state["acooc_exp_sub"] = False

        from components.multi_location_results import (
            render_single_location_result,
            render_location_progress,
        )

        location_tasks = st.session_state.get("acooc_location_tasks", {})

        if not st.session_state.get("acooc_cooc_tasks"):
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
                           + ". Results below, and the buttons in the tabs, "
                           "use the previous run.")
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

            # ── the run's saved settings: every button below uses them, never
            #    the sidebar as it is now (2026-10-08). Before, the deep scan
            #    and cross-check read the sidebar when they started, so a panel
            #    or date change during a run gave a city a scan for other
            #    settings than its graph. ──────────────────────────────────────
            location_names = list(st.session_state.get("acooc_cooc_tasks", {}).keys())
            _run_panel = list(st.session_state.get("acooc_ran_panel") or all_selected_variants)
            _run_d0, _run_d1 = (st.session_state.get("acooc_ran_dates")
                                or (start_date.isoformat(), end_date.isoformat()))

            def _stage_open(stage, tkey, rkey):
                """Cities whose task of this stage was started and is not yet
                collected, failed or stopped."""
                _t = st.session_state.get(tkey, {}) or {}
                _r = st.session_state.get(rkey, {}) or {}
                _f = _failed.get(stage, {})
                return [l for l, tid in _t.items() if tid and l not in _r and l not in _f]

            # ── Buttons (2026-10-08): each sits in the tab where its result
            #    appears (Deconvolution, CovvFit: Deconvolution tab; deep scan:
            #    Lineages tab); ■ Stop in the progress box. ─────────────────────
            def _btn_deconv():
                _open = _stage_open(*_RUN_STAGES[1])
                _label = "Run deconvolution" if not location_tasks else "Re-run deconvolution"
                if st.button(_label, key="acooc_btn_deconv", use_container_width=True,
                             disabled=bool(_open),
                             help=("Running…" if _open else
                                   f"LolliPop per location with the settings of step 4 "
                                   f"({bootstraps} bootstraps, bandwidth {bandwidth}).")):
                    _new = {}
                    for _l in location_names:
                        _new[_l] = celery_app.send_task(
                            "tasks.run_deconvolve_lapis",
                            kwargs={"location": _l, "start_date": _run_d0,
                                    "end_date": _run_d1, "variants": _run_panel,
                                    "bootstraps": bootstraps, "bandwidth": bandwidth}).id
                    st.session_state["acooc_location_tasks"] = _new
                    st.session_state["location_results"] = {}
                    st.session_state["acooc_deconv_settings"] = (bootstraps, bandwidth)
                    _failed.pop("deconvolution", None)
                    st.session_state.setdefault("acooc_phase_t0", {})["deconvolution"] = time.time()
                    logger.info(f"Deconvolution submitted for {len(_new)} location(s)")
                    st.rerun()

            def _variant_colors():
                """{variant: '#rrggbb'}: the deconvolution plots' colours (one per
                variant, the same in every location), as hex for covvfit."""
                from components.multi_location_results import build_variant_color_map
                cm = build_variant_color_map(st.session_state.get("location_results", {}) or {},
                                             location_names)
                out = {}
                for v, c in cm.items():
                    c = str(c).strip()
                    if c.startswith("rgb"):
                        try:
                            r, g, b = [int(float(x)) for x in
                                       c[c.index("(") + 1:c.index(")")].split(",")[:3]]
                            c = f"#{r:02x}{g:02x}{b:02x}"
                        except Exception:
                            continue
                    out[v] = c
                return out

            def _btn_covvfit():
                """CovvFit (2026-10-09): growth advantages of the panel variants
                from the locations whose deconvolution is in. Settings in a
                closed expander under the button."""
                _open = _stage_open(*_RUN_STAGES[4])
                _lr = st.session_state.get("location_results", {}) or {}
                _locs = [l for l in location_names if l in _lr]
                _ready = bool(_locs) and not _stage_open(*_RUN_STAGES[1])
                _done = _ALL in (st.session_state.get("acooc_covvfit_results") or {})
                _h = int(st.session_state.get("acooc_covvfit_horizon", 60))
                if st.button("Re-run CovvFit" if _done else "Run CovvFit",
                             key="acooc_btn_covvfit", use_container_width=True,
                             disabled=bool(_open) or not _ready,
                             help=("Running…" if _open else
                                   "Available once the deconvolution is done." if not _ready else
                                   f"Growth advantage of each panel variant, from the "
                                   f"{len(_locs)} location(s) with a deconvolution; "
                                   f"predicts {_h} days ahead.")):
                    _tid = celery_app.send_task(
                        "tasks.run_covvfit_lapis",
                        kwargs={"locations": _locs, "start_date": _run_d0,
                                "end_date": _run_d1, "variants": _run_panel,
                                "horizon": _h, "colors": _variant_colors()}).id
                    st.session_state["acooc_covvfit_tasks"] = {_ALL: _tid}
                    st.session_state["acooc_covvfit_results"] = {}
                    _failed.pop("covvfit", None)
                    st.session_state.setdefault("acooc_phase_t0", {}).pop("covvfit", None)
                    logger.info(f"CovvFit submitted for {len(_locs)} location(s)")
                    st.rerun()
                with st.expander("CovvFit settings", expanded=False):
                    st.slider("Days to predict after the last date", 7, 120, 60, step=1,
                              key="acooc_covvfit_horizon",
                              help="covvfit's --horizon; the fit itself uses the run's dates.")

            def _btn_scan():
                _scan_req = (bool(st.session_state.get("acooc_want_scan"))
                             or bool(st.session_state.get("acooc_scanner_tasks")))
                _scan_stopped = [l for l, m in _failed.get("scanner", {}).items()
                                 if m == _STOPPED]
                _xc_stopped = any(m == _STOPPED for m in _failed.get("cross-check", {}).values())
                _resume = bool(_scan_stopped) or _xc_stopped
                _wait_graph = _cooc_outstanding and not _scan_req
                _done = (bool(st.session_state.get("acooc_xcheck_tasks"))
                         and not _stage_open(*_RUN_STAGES[2])
                         and not _stage_open(*_RUN_STAGES[3]))
                if st.button("Resume deep scan" if _resume else "Run deep scan",
                             key="acooc_btn_scan", use_container_width=True,
                             disabled=_wait_graph or (_scan_req and not _resume),
                             help=("Waits for the panel check of every location."
                                   if _wait_graph else
                                   ("Done for these results; press ▶ Run for new ones."
                                    if _done else "Running: every location, then the cross-check.")
                                   if _scan_req and not _resume else
                                   "Looks for variants outside your panel in every location; "
                                   "the cross-check follows.")):
                    # after a Stop: the stopped locations scan again; the
                    # cross-check (it needs every scan) starts over
                    _st = st.session_state.get("acooc_scanner_tasks", {})
                    for _l in _scan_stopped:
                        _st.pop(_l, None)
                        _failed.get("scanner", {}).pop(_l, None)
                    _failed.pop("cross-check", None)
                    st.session_state["acooc_xcheck_tasks"] = {}
                    st.session_state["acooc_xcheck_results"] = {}
                    st.session_state["acooc_want_scan"] = True
                    st.session_state["acooc_stopped"] = False
                    st.rerun()

            def _btn_stop():
                if st.button("■ Stop", key="acooc_btn_stop", use_container_width=True,
                             help="Cancels every task of this run that hasn't finished. "
                                  "Results already in stay."):
                    _ids = []
                    for _stg, _tk, _rk in _RUN_STAGES:
                        for _l in _stage_open(_stg, _tk, _rk):
                            _ids.append(st.session_state[_tk][_l])
                            _failed.setdefault(_stg, {})[_l] = _STOPPED
                    if _ids:
                        # terminate=True also ends tasks already running
                        # (the worker runs Celery's default prefork pool)
                        celery_app.control.revoke(_ids, terminate=True)
                    st.session_state["acooc_want_scan"] = False
                    st.session_state["acooc_stopped"] = True
                    logger.info(f"Stop: revoked {len(_ids)} task(s)")
                    st.rerun()

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

            # CovvFit result (one task for all locations, 2026-10-09)
            _cvr_all = st.session_state.setdefault("acooc_covvfit_results", {})
            for _k, _tid in list((st.session_state.get("acooc_covvfit_tasks") or {}).items()):
                if _tid and _k not in _cvr_all and _k not in _failed.get("covvfit", {}):
                    _t = celery_app.AsyncResult(_tid)
                    if _t.ready():
                        try:
                            _cvr_all[_k] = _t.get()
                        except Exception as _e:
                            _fail("covvfit", _k, _e)
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
                _ref_locs = list(location_names)
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
            # phase 2 — the deep scan: the worker reads today's positions plus
            # those where the data shows a mutation and runs the scanner on them
            # (scope.data_positions; the graph stays on phase 1). Only after
            # "Run deep scan" was pressed, which is possible once completeness
            # is in for every location (2026-10-08: was per location as soon as
            # its own completeness was in). All locations start together, with
            # the run's saved settings. The task doesn't use the completeness
            # result; the wait is for the workflow (look at the graph, then decide).
            if st.session_state.get("acooc_want_scan"):
                for _loc in location_names:
                    if (_loc in _cooc_res
                            and _loc not in _scanner_tasks
                            and _loc not in scanner_results):
                        _stask = celery_app.send_task(
                            "tasks.run_cooc_deep_scan_lapis",
                            kwargs={
                                "location": _loc,
                                "start_date": _run_d0,
                                "end_date": _run_d1,
                                "variants": _run_panel,
                                # error hotspots come from every available
                                # location, not only the run's cities
                                "reference_locations": _ref_locs,
                            }
                        )
                        _scanner_tasks[_loc] = _stask.id
                        scanner_panels[_loc] = list(_run_panel)
                        logger.info(f"Submitted deep scan for {_loc}")
            st.session_state["acooc_scanner_tasks"] = _scanner_tasks
            st.session_state["acooc_scanner_panels"] = scanner_panels

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
            if (st.session_state.get("acooc_want_scan")
                    and _scans_done and not _xt and _scan_all):
                _found_in = {}
                for _l, _r in _scan_all.items():
                    for _k in ("resolved_clade", "one_day"):
                        for _c in (_r.get(_k) or []):
                            if _c.get("node"):
                                _found_in.setdefault(_c["node"], set()).add(_l)
                _ran_pan = set(_run_panel)
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
                        "location": _l, "start_date": _run_d0,
                        "end_date": _run_d1, "variants": _miss,
                        "panel": sorted(_ran_pan)})
                    _xt[_l] = _xtask.id
                st.session_state["acooc_xcheck_found_in"] = {n: sorted(cs) for n, cs in _found_in.items()}
                logger.info(f"Cross-check submitted: {sum(1 for t in _xt.values() if t)} cities")

            # ── Autorefresh decision — AFTER collection + scanner submission ───
            # Base it on "is there outstanding work?" rather than raw task state,
            # so newly-submitted scanner tasks keep the refresh alive and the bars
            # update to green without needing a manual click.
            # Three kinds of outstanding work (2026-10-08), so a long
            # deconvolution doesn't hold back the scan's results and the reverse:
            #   completeness — the graph; deconvolution — the LolliPop plots;
            #   scan — deep scan + cross-check, once "Deep scan" was pressed.
            _cr_now = st.session_state.get("acooc_cooc_results", {})
            _sr_now = st.session_state.get("acooc_scanner_results", {})
            _want_scan = bool(st.session_state.get("acooc_want_scan"))
            _cooc_outstanding = any(
                _ln not in _cr_now and _ln not in _failed.get("completeness", {})
                for _ln in location_names)
            _deconv_outstanding = bool(_stage_open(
                "deconvolution", "acooc_location_tasks", "location_results"))
            _scan_outstanding = False
            if _want_scan:
                for _ln in location_names:
                    # requested (completeness is in everywhere) but the scan
                    # not yet done → running
                    if _ln in _failed.get("completeness", {}):
                        continue
                    if _ln not in _sr_now and _ln not in _failed.get("scanner", {}):
                        _scan_outstanding = True
                        break
                # phase 3 (cross-check) not yet in → running
                if not _scan_outstanding and _sr_now:
                    _xt_now = st.session_state.get("acooc_xcheck_tasks", {})
                    _xr_now = st.session_state.get("acooc_xcheck_results", {})
                    if not _xt_now or any(_l not in _xr_now and _l not in _failed.get("cross-check", {})
                                          for _l in location_names):
                        _scan_outstanding = True
            _covvfit_outstanding = bool(_stage_open(*_RUN_STAGES[4]))
            _outstanding = (_cooc_outstanding or _deconv_outstanding or _scan_outstanding
                            or _covvfit_outstanding)
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

            def _deconv_download():
                """Zip of the deconvolution CSVs (were in the progress header)."""
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

            # Progress box (2026-10-08): one summary line for every kind of task
            # that was started, a bar only for what is running now (% of the
            # work from the worker's own progress, time left at the pace so
            # far). Drawn inside a fragment that reruns every 3 s while work is
            # outstanding, so the bars move without redrawing the charts; the
            # whole page reruns only when a task has finished.
            _deep_on = get_cooc_setting("scope.data_positions", default=True)
            _ws = bool(st.session_state.get("acooc_want_scan"))
            # (key, short name, bar title, [(stage, tasks key, results key)]).
            # The deep scan and its cross-check are ONE bar (2026-10-08): the
            # tree updates only when both are done, so "done" must mean both.
            _PHASES = [
                ("completeness", "Panel check", "Panel check", [_RUN_STAGES[0]]),
                ("deconvolution", "Deconvolution", "Deconvolution", [_RUN_STAGES[1]]),
                ("covvfit", "CovvFit", "CovvFit — growth advantages", [_RUN_STAGES[4]]),
                ("scanner", "Deep scan",
                 "Deep scan — then the cross-check of what it found"
                 + ("" if _deep_on else " (panel positions only)"),
                 [_RUN_STAGES[2], _RUN_STAGES[3]]),
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
                """The progress box: a summary line, then one bar per kind of
                task that is running, split in one segment per location (each
                with its name and state under it). For the deep scan a
                location's segment is half its scan, half its cross-check."""
                import time as _time
                _ss = st.session_state
                _fl = _ss.get("acooc_failed", {})
                _t0s = _ss.setdefault("acooc_phase_t0", {})
                _ready = [False]
                _summary, _bars = [], ""

                def _one(stage, tk, rk, l):
                    """(fraction 0-1, state, status message) of one location's task:
                    state = done / failed / stopped / running / none."""
                    _tasks, _done = _ss.get(tk, {}), _ss.get(rk, {})
                    _bad = _fl.get(stage, {})
                    if l in _bad:
                        return 1.0, ("stopped" if _bad[l] == _STOPPED else "failed"), ""
                    if l in _done:
                        return 1.0, "done", ""
                    if _tasks.get(l):
                        if live and celery_app.AsyncResult(_tasks[l]).ready():
                            _ready[0] = True
                        _f, _msg = _task_frac(_tasks[l]) if live else (0.0, "")
                        return _f, "running", _msg
                    return 0.0, "none", ""

                _MARKS = {"done": "✓", "stopped": "■", "failed": "✗"}
                for _pk, _short, _title, _stages in _PHASES:
                    _has_tasks = any(_ss.get(_tk) for _, _tk, _ in _stages)
                    if not (_has_tasks or (_pk == "scanner" and _ws)):
                        continue
                    _segs, _now_txt = [], ""
                    for _l in ([_ALL] if _pk == "covvfit" else location_names):
                        _nm = _l.split("(")[0].strip()
                        _f, _st, _msg = _one(*_stages[0], _l)
                        if len(_stages) == 2 and _st == "done":
                            # deep scan done here → its cross-check
                            _f2, _st2, _msg2 = _one(*_stages[1], _l)
                            _f = 0.5 + 0.5 * _f2
                            _st, _msg = _st2, _msg2
                            _lab = (_MARKS.get(_st2) if _st2 in _MARKS else
                                    "cross-check" if _st2 == "running" else "scanned")
                        elif len(_stages) == 2:
                            _f = 1.0 if _st in ("failed", "stopped") else 0.5 * _f
                            _lab = (_MARKS.get(_st) if _st in _MARKS else
                                    f"{_f * 200:.0f} %" if _st == "running" else "starting")
                        else:
                            _lab = (_MARKS.get(_st) if _st in _MARKS else
                                    f"{_f * 100:.0f} %" if _st == "running" else "waiting")
                        if _msg and not _now_txt:
                            _now_txt = f"{_nm}: {_msg}"
                        _segs.append((_nm, _f, _st, _lab))
                    _frac = sum(x[1] for x in _segs) / len(_segs) if _segs else 0.0
                    _n_bad = sum(1 for x in _segs if x[2] == "failed")
                    _n_stop = sum(1 for x in _segs if x[2] == "stopped")
                    _extra = ((f" · {_n_bad} failed" if _n_bad else "")
                              + (f" · {_n_stop} stopped" if _n_stop else ""))
                    if _frac >= 0.999:
                        _summary.append(f"{_short} {'done' if not (_n_bad or _n_stop) else 'ended'}{_extra}")
                        continue
                    _summary.append(f"{_short} running{_extra}")
                    if _pk not in _t0s:
                        _t0s[_pk] = _time.time()
                    _right = f"{_frac * 100:.0f}%"
                    if live and _frac > 0.03:
                        _el = _time.time() - _t0s[_pk]
                        _right += " · " + _left(_el * (1 - _frac) / _frac)
                    # one segment per location; at < 110 px each (many
                    # locations) they wrap onto a second row
                    _track = "".join(
                        f"<div style='flex:1 1 0;min-width:110px;'>"
                        f"<div style='height:8px;background:#D3D1C7;border-radius:4px;overflow:hidden;'>"
                        f"<div style='height:100%;width:{f * 100:.0f}%;background:#475569;"
                        f"border-radius:4px;'></div></div>"
                        f"<div style='font-size:12px;color:#5F5E5A;margin-top:3px;white-space:nowrap;"
                        f"overflow:hidden;text-overflow:ellipsis;'>{nm} {lab}</div></div>"
                        for nm, f, _, lab in _segs)
                    _sub = (f"<div style='font-size:12px;color:#898781;margin-top:2px;"
                            f"white-space:nowrap;overflow:hidden;text-overflow:ellipsis;'>"
                            f"{_now_txt}</div>" if (live and _now_txt) else "")
                    _bars += (
                        f"<div style='margin-top:10px;'>"
                        f"<div style='display:flex;justify-content:space-between;font-size:13px;"
                        f"margin-bottom:4px;'><span style='font-weight:500;'>{_title}</span>"
                        f"<span style='color:#898781;font-size:12px;'>{_right}</span></div>"
                        f"<div style='display:flex;flex-wrap:wrap;gap:6px 4px;'>{_track}</div>{_sub}</div>")
                st.markdown(f"<div style='font-size:13px;color:#5F5E5A;'>{' · '.join(_summary)}</div>"
                            + _bars, unsafe_allow_html=True)
                return _ready[0]

            _pbox = st.container(border=True)
            with _pbox:
                _pb1, _pb2 = st.columns([5, 1])
                with _pb1:
                    if _outstanding:
                        @st.fragment(run_every=3)
                        def _poll_tasks():
                            if _phase_bars(live=True):
                                st.rerun(scope="app")
                        _poll_tasks()
                    else:
                        _phase_bars(live=False)
                with _pb2:
                    if _outstanding:
                        _btn_stop()
            with _pbox:
              _fl_all = st.session_state.get("acooc_failed", {})
              _err = {_stg: [c for c, m in _cs.items() if m != _STOPPED]
                      for _stg, _cs in _fl_all.items()}
              _stp = {_stg: [c for c, m in _cs.items() if m == _STOPPED]
                      for _stg, _cs in _fl_all.items()}

              _SHOWN = {"completeness": "panel check", "scanner": "deep scan", "covvfit": "CovvFit"}

              def _by_stage(d):
                  return " · ".join(f"{_SHOWN.get(_stg, _stg)} in "
                                    f"{', '.join(c.split('(')[0].strip() for c in _cs)}"
                                    for _stg, _cs in d.items() if _cs)
              if any(_err.values()):
                  st.error("Failed: " + _by_stage(_err)
                           + " — the worker may have restarted or run out of memory; re-run to retry.")
              if any(_stp.values()):
                  st.info("Stopped: " + _by_stage(_stp)
                          + ". Results already in are kept; press a button again to restart that part.")

            st.markdown("<hr style='margin:10px 0 8px;opacity:.15;'>", unsafe_allow_html=True)

            # ── Section switcher (top-level tabs) ─────────────────────────────
            # Buttons persist selection across autorefresh (st.tabs ghosted +
            # reset to the first tab on each rerun). Grouped (2026-10-08) in the
            # order results arrive: panel and abundance (is the panel good
            # enough for a deconvolution? then the deconvolution), the deep scan
            # (lineages, signal over time), and the look-up. A mark on the label
            # tells the state of the tab's task: ✓ done, ⟳ running, ■ stopped,
            # ○ not started.
            _GROUPS = [("Panel and abundance", ["Panel check", "Deconvolution"]),
                       ("Deep scan", ["Lineages", "Signal over time"]),
                       ("Look-up", ["Investigate a variant"])]
            _all_secs = [_x for _, _xs in _GROUPS for _x in _xs]
            if st.session_state.get("acooc_section") not in _all_secs:
                st.session_state["acooc_section"] = "Panel check"

            def _mark(running, stages):
                if running:
                    return " ⟳"
                if any(m == _STOPPED for _sg in stages for m in _failed.get(_sg[0], {}).values()):
                    return " ■"
                if any(st.session_state.get(_sg[1]) for _sg in stages):
                    return " ✓"
                return " ○"
            _MARK = {"Panel check": _mark(_cooc_outstanding, _RUN_STAGES[0:1]),
                     "Deconvolution": _mark(_deconv_outstanding or _covvfit_outstanding,
                                            _RUN_STAGES[1:2]),
                     "Lineages": _mark(_scan_outstanding, _RUN_STAGES[2:4])}
            _spec = []
            for _gi, (_, _xs) in enumerate(_GROUPS):
                if _gi:
                    _spec.append(0.06)                 # separator
                _spec += [1.15 if _x == "Investigate a variant" else 1 for _x in _xs]
            _SEP = ("<div style='border-left:1px solid rgba(0,0,0,.18);height:{h}px;"
                    "width:0;margin:0 auto;'></div>")
            _lab_cols, _btn_cols = st.columns(_spec), st.columns(_spec)
            _ci = 0
            for _gi, (_gname, _xs) in enumerate(_GROUPS):
                if _gi:
                    _lab_cols[_ci].markdown(_SEP.format(h=16), unsafe_allow_html=True)
                    _btn_cols[_ci].markdown(_SEP.format(h=38), unsafe_allow_html=True)
                    _ci += 1
                _lab_cols[_ci].markdown(f"<div style='font-size:11px;color:#898781;'>{_gname}</div>",
                                        unsafe_allow_html=True)
                for _snm in _xs:
                    with _btn_cols[_ci]:
                        _active = st.session_state["acooc_section"] == _snm
                        if st.button(_snm + _MARK.get(_snm, ""), key=f"acooc_sec_{_snm}",
                                     use_container_width=True,
                                     type="primary" if _active else "secondary"):
                            st.session_state["acooc_section"] = _snm
                            st.rerun()
                    _ci += 1
            _active_section = st.session_state["acooc_section"]
            st.markdown("<hr style='margin:2px 0 10px;opacity:.12;'>", unsafe_allow_html=True)

            if _active_section == "Deconvolution":
              st.markdown("#### Deconvolution")
              st.caption("Variant shares over time per location, with LolliPop and the "
                         "settings of step 4.")
              _dc1, _dc2, _dc3 = st.columns(3)
              with _dc1:
                  _btn_deconv()
              with _dc2:
                  _btn_covvfit()
              with _dc3:
                  _deconv_download()
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
                          # one colour per variant in every location, the same
                          # colours CovvFit gets (2026-10-09; was by position)
                          from components.multi_location_results import (
                              build_variant_color_map, create_variant_plot)
                          _rd = st.session_state.location_results[location]
                          _vd = _rd.get(location) if isinstance(_rd, dict) and location in _rd else _rd
                          _fig = create_variant_plot(
                              _vd, location, color_map=build_variant_color_map(
                                  st.session_state.location_results, location_names))
                          if _fig is not None:
                              st.plotly_chart(_fig, use_container_width=True,
                                              key=f"deconv_plot_{location}")
                          else:
                              st.warning("No plottable deconvolution data.")
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
                          variants=_run_panel,
                          pango_loader=cached_get_pango_loader(),
                      )
              # deconvolution has its own button (2026-10-08)
              _dset = st.session_state.get("acooc_deconv_settings")
              if not location_tasks:
                  st.info("Deconvolution hasn't run for these results yet: press "
                          "Run deconvolution.")
              elif _dset:
                  st.caption(f"LolliPop: {_dset[0]} bootstraps, bandwidth {_dset[1]}.")
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

              # ── CovvFit: growth advantages (2026-10-09) ────────────────────
              _cvr = (st.session_state.get("acooc_covvfit_results") or {}).get(_ALL)
              _cvf = _failed.get("covvfit", {}).get(_ALL)
              if _cvr or _covvfit_outstanding or _cvf:
                  st.markdown("---")
                  st.markdown("#### Growth advantages (CovvFit)")
                  if _covvfit_outstanding:
                      st.caption("CovvFit running (progress above).")
                  elif _cvf == _STOPPED:
                      st.caption("CovvFit was stopped. Press Re-run CovvFit to start it again.")
                  elif _cvf:
                      st.error(f"CovvFit failed: {_cvf}")
                  if _cvr and not _covvfit_outstanding:
                      _render_covvfit(_cvr)


            if _active_section in ("Panel check", "Lineages"):
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
              if _graph_locs and _active_section == "Panel check":
                  st.markdown("#### Is your panel good enough for deconvolution?")
                  st.caption("Share of each location's reads your panel variants explain. "
                             "Mostly green = the deconvolution can split the signal "
                             "between them.")
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
                      if (_cr_all.get(_lc) is not None and _sr_all.get(_lc) is None
                              and st.session_state.get("acooc_want_scan")):
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
                  st.caption("Grey = reads your panel doesn't explain. Red appears after "
                             "the deep scan (Lineages tab).")
                  _n1, _n2, _n3 = st.columns(3)
                  with _n1:
                      if st.button("Go to deconvolution →", key="acooc_go_deconv",
                                   use_container_width=True):
                          st.session_state["acooc_section"] = "Deconvolution"
                          st.rerun()
                  with _n2:
                      if st.button("Find what's missing (deep scan) →", key="acooc_go_scan",
                                   use_container_width=True):
                          st.session_state["acooc_section"] = "Lineages"
                          st.rerun()
                  with _n3:
                      st.button("Suggest a panel · release 2", key="acooc_suggest_panel",
                                use_container_width=True, disabled=True,
                                help="Coming in release 2: a panel proposed from the data.")

              if _active_section == "Lineages":
                st.markdown("#### Lineages")
                st.caption("Which lineages are really there. Your panel rows have their ★ "
                           "check now; the deep scan adds what's outside your panel, then the "
                           "cross-check checks them in every location.")
                _ls1, _ls2 = st.columns([1, 2])
                with _ls1:
                    _btn_scan()
                if _cooc_outstanding and not st.session_state.get("acooc_want_scan"):
                    with _ls2:
                        st.caption("Available once the panel check is in for "
                                   "every location.")
                elif _scan_outstanding:
                    with _ls2:
                        st.caption("Deep scan running: the lineages it finds appear "
                                   "when the cross-check ends (progress above).")
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
                # (2026-10-08: only the scan's own work counts — a deconvolution
                # still running doesn't hide the findings)
                _scan_running = _scan_outstanding
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
                        """The ★ check for a found lineage in one city (the deep
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
                            # the ★ check on its markers (2026-10-02)
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
                        # + Add (2026-10-02): the ★ check would confirm it in at
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
                  cities=None if (_cooc_outstanding or _scan_outstanding) else list(location_names), start_date=_rd[0], end_date=_rd[1],
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
                  scanner_results=({} if _scan_outstanding else
                                   st.session_state.get("acooc_scanner_results", {}) or {}),
                  default_city=st.session_state.get("acooc_tree_city"),
                  scanning=_scan_outstanding,
              )



if __name__ == "__main__":
    app()