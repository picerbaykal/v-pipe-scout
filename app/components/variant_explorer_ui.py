"""components/variant_explorer_ui.py

'Investigate a variant': pick any lineage, see whether co-occurrence can tell
it apart (★ markers, the same rule as the panel check) and how it relates to
others, then "Check in data": the worker reads its ★ markers in the run's
cities and window, one line per city (2026-10-02). The heatmaps live in their
own section, "Signal over time" (components/signal_over_time_ui.py).
"""
from __future__ import annotations

import streamlit as st

from process.variant_explorer import check_in_data, investigate_variant
_STATE = {
    "present":     ("#dcfce7", "#166534", "present"),
    "absent":      ("#f3f4f6", "#4b5563", "absent"),
    "mixed":       ("#fef3c7", "#92400e", "mixed"),
    "not_covered": ("#ffffff", "#6b7280", "not covered"),
    "no_marker":   ("#f3f4f6", "#6b7280", "no ★ marker"),
}


def _pill(bg, fg, txt, border=None):
    b = f"border:1px solid {border};" if border else ""
    return (f"<span style='background:{bg};color:{fg};{b}font-size:11px;font-weight:600;"
            f"padding:1px 8px;border-radius:10px;'>{txt}</span>")


def _check_section(variant, panel, cities, start_date, end_date, celery_app, key_prefix):
    if not (cities and start_date and end_date and celery_app):
        st.caption("Run the analysis (or wait for it to finish) to check it in the "
                   "data of your cities.")
        return
    k = f"{key_prefix}_check"
    cur = st.session_state.get(k)
    win = (str(start_date)[:10], str(end_date)[:10])
    sig = (variant, tuple(sorted(cities)), win, tuple(sorted(panel or [])))
    if cur and cur.get("sig") != sig:
        cur = None
    n_city = len(cities)
    if cur is None:
        if st.button(f"Check {variant} in data · {n_city} cit{'y' if n_city == 1 else 'ies'}",
                     key=f"{key_prefix}_go", type="primary"):
            tasks = {}
            for c in cities:
                t = celery_app.send_task("tasks.run_cooc_variant_check_lapis", kwargs={
                    "location": c, "start_date": win[0], "end_date": win[1],
                    "variant": variant, "panel": list(panel or [])})
                tasks[c] = t.id
            st.session_state[k] = {"sig": sig, "tasks": tasks, "res": {}, "err": {}}
            st.rerun()
        return

    def _collect():
        for c, tid in cur.get("tasks", {}).items():
            if c in cur["res"] or c in cur["err"]:
                continue
            t = celery_app.AsyncResult(tid)
            if t.ready():
                try:
                    r = t.get()
                    cur["res"][c] = check_in_data(r.get("markers") or [],
                                                  r.get("per_date") or {},
                                                  r.get("dates") or [])
                except Exception as e:
                    cur["err"][c] = str(e)[:200]
        return [c for c in cur.get("tasks", {}) if c not in cur["res"] and c not in cur["err"]]

    def _draw():
        rows = []
        for c in cities:
            name = c.split("(")[0].strip()
            if c in cur["err"]:
                right = "<span style='color:#dc2626;'>failed — try again</span>"
            elif c in cur["res"]:
                r = cur["res"][c]
                bg, fg, word = _STATE[r["state"]]
                _low = r["n_markers"] - r["n_measured"]
                cnt = (f"<span style='color:#6b7280;'>{r['n_present']} of {r['n_measured']} "
                       f"measurable ★ present"
                       + (f" ({_low} with too few reads)" if _low > 0 else "")
                       + "</span>" if r["n_measured"] else "")
                right = (_pill(bg, fg, word, "#d1d5db" if bg == "#ffffff" else None)
                         + f" {cnt}")
            else:
                right = "<span style='color:#2563eb;'>⟳ reading…</span>"
            rows.append(f"<div style='display:flex;align-items:center;gap:10px;"
                        f"padding:3px 0;'><span style='width:90px;font-weight:500;'>{name}</span>"
                        f"<span>{right}</span></div>")
        st.markdown("<div style='font-size:12.5px;margin-top:6px;'>" + "".join(rows)
                    + "</div><div style='font-size:11px;color:#6b7280;margin-top:2px;'>"
                    "Pooled over the run's window · present / absent by the panel check "
                    "(★ marker ≥ 5 % / < 1 % of ≥ 100 reads) · its signal per week: "
                    "<b>Signal over time</b></div>", unsafe_allow_html=True)

    if _collect():
        # poll only while something is still running; the whole page reruns
        # once when the last city is in, which draws the final, static view
        @st.fragment(run_every=2)
        def _poll():
            if not _collect():
                st.rerun(scope="app")
            _draw()
        _poll()
    else:
        _draw()
        if st.button("Check again", key=f"{key_prefix}_again"):
            st.session_state.pop(k, None)
            st.rerun()


def render_variant_explorer(pango_loader, panel=None, options=None,
                            disabled=False, key_prefix="acooc_explorer",
                            cities=None, start_date=None, end_date=None,
                            celery_app=None):
    """Render 'Investigate a variant'.

    Args:
        pango_loader: provides get_raw_data() / get_signature().
        panel: the run's panel (★ markers are computed as if the variant were added).
        options: selectable lineage names (defaults to all known).
        disabled: True while a run is still going (avoids rerun races).
        cities, start_date, end_date, celery_app: enable "Check in data" for
            these cities and window.
    """
    st.markdown("#### 🔎 Investigate a variant")
    st.caption("Quick look-up of any lineage, in the panel or not: can reads tell it "
               "apart (★ markers), how it relates to others — and, after a run, is it "
               "in your cities' data.")

    if disabled:
        st.info("Available once the run completes.")
        return

    if options is None:
        options = sorted(pango_loader.get_raw_data().keys())

    # index=None + placeholder: the box starts empty, so you can type at once
    # (an "" option shown as text had to be deleted first)
    variant = st.selectbox(
        "Variant", options=options, index=None, placeholder="Search a lineage…",
        key=f"{key_prefix}_pick", label_visibility="collapsed",
    )
    if not variant:
        return

    r = investigate_variant(variant, pango_loader, panel=panel)
    if not r.get("found"):
        st.info(r.get("reason", "Not found."))
        return

    badge = (_pill("#e6f4ef", "#0f6e56", f"★ {len(r['markers'])} markers") if r["detectable"]
             else _pill("#f3f4f6", "#6b7280", "no ★ marker"))
    recomb = (" " + _pill("#f3e8f1", "#8a4e82", "recombinant")) if r["is_recombinant"] else ""
    panel_b = (" " + _pill("#e6f0f9", "#185fa5", "in panel")) if r["in_panel"] else ""
    st.markdown(f"<div style='font-size:15px;font-weight:600;margin-top:4px;'>{r['name']} "
                f"{badge}{recomb}{panel_b}</div>", unsafe_allow_html=True)

    lines = [f"<b>Markers:</b> {r['reason']}"]
    if r["parent"]:
        lines.append(f"<span style='color:#6b7280;'>Belongs to:</span> <b>{r['parent']}</b>")
    if r["oscillates_with"]:
        lines.append(
            "<span style='color:#6b7280;'>Near-identical to (can't be told apart):</span> "
            f"<span style='background:#fdf4e6;color:#ba7517;padding:0 6px;border-radius:8px;'>"
            f"{', '.join(r['oscillates_with'][:6])}</span>")
    if r["children"]:
        lines.append("<span style='color:#6b7280;'>Descendants:</span> "
                     + ", ".join(r["children"][:6]))
    st.markdown("<div style='border:0.5px solid #e5e7eb;border-radius:8px;padding:8px 12px;"
                "font-size:12.5px;line-height:1.9;margin-top:4px;'>" + "<br>".join(lines)
                + "</div>", unsafe_allow_html=True)

    if r["detectable"]:
        _check_section(variant, panel, cities, start_date, end_date, celery_app, key_prefix)