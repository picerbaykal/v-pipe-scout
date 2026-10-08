"""components/variant_explorer_ui.py

'Investigate a variant': pick any lineage, see whether co-occurrence can tell
it apart (★ markers, the same rule as the panel check) and how it relates to
others, then "Check in data": the worker reads its ★ markers in the run's
cities and window, one line per city (2026-10-02). The heatmaps live in their
own section, "Signal over time" (components/signal_over_time_ui.py).
"""
from __future__ import annotations

import streamlit as st

from process.variant_explorer import check_in_data, investigate_variant, marker_funnel
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


def _mono(muts, n=None):
    ms = list(muts)[:n] if n else list(muts)
    more = f" … +{len(muts) - len(ms)}" if n and len(muts) > len(ms) else ""
    return ("<span style='font-family:ui-monospace,Menlo,monospace;font-size:12px'>"
            + " ".join(ms) + "</span>" + more)


def _funnel(variant, pango_loader, panel):
    """'How the ★ markers were chosen': the rule step by step, with what each
    step removed and the near misses (marker_funnel)."""
    f = marker_funnel(variant, pango_loader, panel)
    lim = f["settings"]["out_other_max"]
    n1 = len(f["signature"])
    n2 = len(f["rows"])
    n3 = len(f["markers"])
    with st.expander(f"How the ★ markers were chosen · {n1} → {n2} → {n3}", expanded=False):
        box = "border-left:3px solid #e5e7eb;padding:2px 10px;margin:6px 0;font-size:12.5px;line-height:1.7"
        st.markdown(
            f"<div style='{box}'><b>1 · Signature</b> — every substitution of {variant} "
            f"(deletions left out): <b>{n1}</b><br>{_mono(f['signature'], 40)}</div>",
            unsafe_allow_html=True)
        by_who = {}
        for m, who in f["by_panel"].items():
            by_who.setdefault(", ".join(who), []).append(m)
        drop2 = "".join(f"<br>✗ also in <b>{w}</b> ({len(ms)}): {_mono(ms, 20)}"
                        for w, ms in sorted(by_who.items(), key=lambda kv: -len(kv[1])))
        st.markdown(
            f"<div style='{box}'><b>2 · Not in another panel variant</b> — a mutation another "
            f"variant of your panel carries can't point to {variant} (its own descendants "
            f"excepted): <b>{n1} → {n2}</b>{drop2 or '<br>nothing dropped'}</div>",
            unsafe_allow_html=True)
        rows = f["rows"]
        kept = [r for r in rows if r["outside"] <= lim]
        bad = [r for r in rows if r["outside"] > lim]
        forg = [r for r in rows if r["forgiven"]]

        def _r(r, mark):
            t = f"{mark} {_mono([r['mutation']])} · {r['outside']} outside"
            if r["forgiven"]:
                t += f" (+{len(r['forgiven'])} forgiven: {', '.join(r['forgiven'][:3])})"
            if r["outside"] and mark != "★":
                t += f" · e.g. {', '.join(r['examples'][:3])}"
            return t
        # kept markers on one line ("150G·0 200C·1" = outsiders), the forgiven
        # and the near misses spelled out
        lines = []
        if kept:
            lines.append(f"★ {len(kept)} kept: " + _mono(
                [f"{r['mutation']}·{r['outside']}" for r in kept], 30)
                + " <span style='color:#6b7280'>(·n = outsiders)</span>")
        lines += [_r(r, "★") + " — kept thanks to forgiven recombinants"
                  for r in kept if r["forgiven"]][:5]
        lines += [_r(r, "⚠") for r in f["near"]]
        n_other_bad = len(bad) - len(f["near"])
        if n_other_bad:
            lines.append(f"✗ {n_other_bad} more with more than {lim + 5} outsiders")
        st.markdown(
            f"<div style='{box}'><b>3 · Few carriers outside the family</b> — lineages outside "
            f"{variant} + its descendants that carry it; a recombinant <i>made from</i> "
            f"{variant} that inherited it is forgiven (option C); kept if ≤ <b>{lim}</b>: "
            f"<b>{n2} → {n3}</b><br>" + "<br>".join(lines or ["nothing left to test"])
            + (f"<br><span style='color:#b45309'>⚠ near miss = rejected by 1–5 outsiders "
               f"over the limit</span>" if f["near"] else "")
            + "</div>", unsafe_allow_html=True)
        st.markdown(
            f"<div style='font-size:11.5px;color:#6b7280;margin-top:4px'>Settings: "
            f"<code>markers.out_other_max = {lim}</code> · panel: "
            f"{', '.join(f['settings']['panel']) or '—'} · {len(forg)} mutation(s) kept "
            f"thanks to forgiven recombinants</div>", unsafe_allow_html=True)


def _detail_table(r):
    """Per-marker numbers behind one city's verdict."""
    rows = []
    for x in r.get("detail") or []:
        f = x["freq"]
        num = (f"{f * 100:.1f} % of {x['cov']:,}" if f is not None and x["cov"]
               else f"{x['cov']:,} reads")
        link = (f"link {x['link'] * 100:.0f} % of {x['link_n']:,}" if x["link"] is not None
                else "")
        flag = (f" <span style='color:#b45309'>⚠ {x['near']}</span>" if x["near"] else "")
        ok = "✓" if x["why"] == "present" else "·"
        rows.append(f"<tr><td style='padding-right:12px;font-family:ui-monospace,Menlo,monospace'>"
                    f"{ok} {x['marker']}</td><td style='padding-right:12px'>{num}</td>"
                    f"<td style='padding-right:12px;color:#6b7280'>{link}</td>"
                    f"<td>{x['why']}{flag}</td></tr>")
    return "<table style='font-size:12px;margin:2px 0 8px 12px'>" + "".join(rows) + "</table>"


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
        if cur.get("res"):
            from process.cooc import _check_cfg
            c = _check_cfg()
            with st.expander("Reads at the markers, per city", expanded=False):
                st.markdown(
                    f"<div style='font-size:11.5px;color:#6b7280'>The same ★ markers everywhere; "
                    f"what differs per city is the reads at them. Pooled over the run's window · "
                    f"present ≥ {c['present_freq'] * 100:g} % · absent &lt; "
                    f"{c['absent_freq'] * 100:g} % · ≥ {c['min_cov']} reads · link ≥ "
                    f"{c['link_min'] * 100:g} % of ≥ {c['min_link']} reads · verdict: present "
                    f"if ≥ {c['confirm_share'] * 100:g} % of measurable markers are present, "
                    f"absent if ≤ {c['notfound_share'] * 100:g} % · ⚠ = just missed a rule"
                    f"</div>", unsafe_allow_html=True)
                for city in cities:
                    if city in cur["res"]:
                        st.markdown(f"<div style='font-size:12.5px;font-weight:600;"
                                    f"margin-top:6px'>{city.split('(')[0].strip()}</div>"
                                    + _detail_table(cur["res"][city]), unsafe_allow_html=True)
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
    try:
        _funnel(variant, pango_loader, panel)
    except Exception as e:                          # never break the look-up
        st.caption(f"(marker steps unavailable: {e})")

    if r["detectable"]:
        _check_section(variant, panel, cities, start_date, end_date, celery_app, key_prefix)