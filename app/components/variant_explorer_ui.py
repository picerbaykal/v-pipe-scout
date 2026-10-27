"""components/variant_explorer_ui.py

'Investigate a variant': pick any lineage, see whether co-occurrence can tell
it apart (★ markers, the same rule as the ★ check of panel variants) and how it relates to
others, then "Check in data": the worker reads it in the run's cities and window, and
the page shows it as a row of the tree, with the tree's rules (2026-10-27;
until then one present / absent / mixed line per city, the marker vote). The heatmaps live in their
own section, "Signal over time" (components/signal_over_time_ui.py).
"""
from __future__ import annotations

import streamlit as st

from process.variant_explorer import investigate_variant, marker_funnel


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


def _check_section(variant, panel, cities, start_date, end_date, celery_app, key_prefix,
                   pango_loader=None, quality=None, clinical=None, pango_parent=None):
    """Check in data (2026-10-27): the worker reads the lineage in each city,
    7 weeks back like the tree's check; the page judges it with the tree's
    rules (process.variant_explorer.judge) and draws it as a tree row
    (components.variants_table.investigate_view)."""
    if not (cities and start_date and end_date and celery_app):
        st.caption("Run the analysis (or wait for it to finish) to check it in the "
                   "data of your cities.")
        return
    from process.national import check_start
    k = f"{key_prefix}_check"
    cur = st.session_state.get(k)
    win = (str(start_date)[:10], str(end_date)[:10])
    d_check = check_start(*win)
    sig = ("v2", variant, tuple(sorted(cities)), win, tuple(sorted(panel or [])))
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
                    "variant": variant, "panel": list(panel or []),
                    "check_start_date": d_check})
                tasks[c] = t.id
            st.session_state[k] = {"sig": sig, "tasks": tasks, "raw": {}, "err": {}}
            st.rerun()
        return

    def _collect():
        for c, tid in cur.get("tasks", {}).items():
            if c in cur["raw"] or c in cur["err"]:
                continue
            t = celery_app.AsyncResult(tid)
            if t.ready():
                try:
                    cur["raw"][c] = t.get()
                except Exception as e:
                    cur["err"][c] = str(e)[:200]
        return [c for c in cur.get("tasks", {}) if c not in cur["raw"] and c not in cur["err"]]

    def _progress():
        rows = []
        for c in cities:
            name = c.split("(")[0].strip()
            right = ("<span style='color:#dc2626;'>failed — try again</span>" if c in cur["err"]
                     else "<span style='color:#166534;'>read</span>" if c in cur["raw"]
                     else "<span style='color:#2563eb;'>⟳ reading…</span>")
            rows.append(f"<div style='display:flex;gap:10px;padding:2px 0;'>"
                        f"<span style='width:90px;font-weight:500;'>{name}</span>{right}</div>")
        st.markdown("<div style='font-size:12.5px;margin-top:6px;'>" + "".join(rows) + "</div>",
                    unsafe_allow_html=True)

    if _collect():
        # poll only while something is still running; the whole page reruns
        # once when the last city is in, which draws the final view
        @st.fragment(run_every=2)
        def _poll():
            if not _collect():
                st.rerun(scope="app")
            _progress()
        _poll()
        return

    # every city in: judge once, with the tree's rules
    if "judged" not in cur:
        from process.variant_explorer import judge, own_of
        try:
            own = own_of(variant, pango_loader) if pango_loader is not None else None
        except Exception:
            own = None
        raw = {}
        for c, r in cur["raw"].items():
            if r and r.get("markers"):
                raw[c] = {"markers": r["markers"], "per_date": r.get("per_date") or {},
                          "dates": list(r.get("dates") or []),
                          "quality": (quality or {}).get(c), "own": own}
        cur["judged"] = judge(raw, win[0], win[1], d_check)
    from components import variants_table as _vt
    shown = [c for c in cities if c in cur["raw"]]
    sel = cur.get("sel") if cur.get("sel") in shown else (shown[0] if shown else None)
    if cur["err"]:
        st.caption("Failed: " + ", ".join(c.split("(")[0].strip() for c in cur["err"])
                   + " — Check again to retry.")
    if shown:
        view = _vt.investigate_view(variant, cur["judged"], shown, sel,
                                    in_panel=variant in (panel or []),
                                    current_panel=panel or [], clinical=clinical,
                                    pango_parent=pango_parent, dates=win)
        click = _vt.render(view, key=f"{key_prefix}_view")
        if click and click.get("t") != cur.get("last_click"):
            cur["last_click"] = click.get("t")
            if click.get("add"):
                st.session_state[f"acooc_add_variant_pending_{click['add']}"] = click["add"]
            elif click.get("city") in shown:
                cur["sel"] = click["city"]
            st.rerun()
    if st.button("Check again", key=f"{key_prefix}_again"):
        st.session_state.pop(k, None)
        st.rerun()


def render_variant_explorer(pango_loader, panel=None, options=None,
                            key_prefix="acooc_explorer",
                            cities=None, start_date=None, end_date=None,
                            celery_app=None, quality=None, clinical=None,
                            pango_parent=None):
    """Render 'Investigate a variant'.

    Args:
        pango_loader: provides get_raw_data() / get_signature().
        panel: the run's panel (★ markers are computed as if the variant were added).
        options: selectable lineage names (defaults to all known).
        cities, start_date, end_date, celery_app: enable "Check in data" for
            these cities and window.
        quality: {city: sample_qc} of the run, clinical / pango_parent: the
            details' clinical rows — so the answer is the tree's (2026-10-27).
    """
    st.markdown("#### 🔎 Investigate a variant")
    st.caption("Quick look-up of any lineage, in the panel or not: can reads tell it "
               "apart (★ markers), how it relates to others — and, after a run, is it "
               "in your cities' data.")

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
        _check_section(variant, panel, cities, start_date, end_date, celery_app, key_prefix,
                       pango_loader=pango_loader, quality=quality, clinical=clinical,
                       pango_parent=pango_parent)