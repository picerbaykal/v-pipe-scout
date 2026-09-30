"""
components/variants_table.py

The "Variants" view of the Co-occurrence results page: the panel tree and
the per-city evidence in ONE table (2026-09-30).

  Variant   the tree, coloured by the evidence in the chosen city; the
            recombinant families (XEC, XFG, XDV → NB.1.8.1 …) under their
            own root, not under B; "+ Add" next to a confirmed finding
  cities    one cell per city — click a city code to colour the tree by it
  Evidence  markers / regions / ⚠ — the facts to judge

Below the tree: unnamed signal — novel patterns (recurring, then 1-day, in
expandable groups) and patterns too broad to name, each with per-city days.

Rendered as a small bidirectional component (plain HTML + the Streamlit
component message protocol, no build step): hover works inside it, and a
click returns {"add": node, "t": ms} or {"city": name, "t": ms} to Python.
"""
from __future__ import annotations

import html
import re
from pathlib import Path

# one palette for cells, tree dots and legend
GREEN = "#16a34a"          # panel variant present
AMBER = "#f59e0b"          # panel variant not found (fill / dot)
AMBER_T = "#b45309"        # … as text
VIOLET = "#7c3aed"         # markers disagree (panel) · named only (finding)
GREY = "#9ca3af"           # no ★ marker / no data
RED = "#dc2626"            # found, not in panel
BLUE = "#2563eb"           # novel
SLATE = "#64748b"          # too broad to name
TRACKED = "#185FA5"

_PANEL = {  # check state -> (fill, text-on-fill, symbol, label, name colour)
    "confirmed": (GREEN, "#fff", "✓", "present", GREEN),
    "not_found": (AMBER, "#422006", "✗", "not found", AMBER_T),
    "inconsistent": (VIOLET, "#fff", "≠", "markers disagree", VIOLET),
    "no_marker": ("#e5e7eb", "#4b5563", "·", "no ★ marker", "#6b7280"),
    "no_data": (None, GREY, "?", "too few reads", "#6b7280"),
}

CSS = """
<style>
.vt { border-collapse:collapse; font-size:15px; width:100%; }
.vt th { text-align:left; font-size:12px; letter-spacing:.03em; text-transform:uppercase; color:#6b7280;
         font-weight:600; padding:0 8px 6px; border-bottom:1px solid #e6e8ec; white-space:nowrap; }
.vt th.cc { width:34px; padding:0 2px 6px; text-align:center; letter-spacing:0; }
.vt th.cc button { font:inherit; font-size:12px; font-weight:700; color:#6b7280; background:none; border:0;
                   border-radius:5px; padding:2px 5px; cursor:pointer; }
.vt th.cc button:hover { background:#eef0f3; color:#31333f; }
.vt th.cc.sel button { background:#31333f; color:#fff; }
.vt td { padding:0 8px; height:34px; border-bottom:1px solid #eef0f2; white-space:nowrap; vertical-align:middle; }
.vt td.cc { padding:0 2px; width:34px; text-align:center; }
.vt td.cc.sel { background:rgba(49,51,63,.07); }
.vt tr.sp td { height:24px; border-bottom:none; }
.vt td.tr { position:relative; padding-right:12px; }
.vt .v { position:absolute; top:0; bottom:0; border-left:1.5px solid #d3d1c7; }
.vt .v.half { bottom:50%; }
.vt .v.down { top:50%; }
.vt .h { position:absolute; top:50%; height:0; border-top:1.5px solid #d3d1c7; }
.vt .dot { position:absolute; top:50%; border-radius:50%; box-sizing:border-box; }
.vt .nm { font-weight:600; }
.vt .spl { color:#9a988f; font-size:12.5px; display:inline-block; max-width:280px; overflow:hidden;
           text-overflow:ellipsis; vertical-align:middle; }
.vt .par { color:#8a4e82; font-size:12.5px; margin-left:8px; }
.vt .grp-root { color:#6b7280; font-size:12.5px; font-weight:700; letter-spacing:.02em; }
.vt .chip { font-size:11px; padding:0 7px; border-radius:8px; margin-left:6px; vertical-align:1px;
            background:#f1efe8; color:#5f5e5a; }
.vt .add { font:inherit; font-size:12.5px; font-weight:600; margin-left:9px; padding:0 9px; line-height:20px;
           border-radius:10px; border:1.5px solid #dc2626; color:#dc2626; background:#fff; cursor:pointer;
           vertical-align:1px; }
.vt .add:hover { background:#dc2626; color:#fff; }
.vt .added { font-size:12.5px; margin-left:9px; color:#6b7280; }
.vt .c { display:inline-flex; align-items:center; justify-content:center; width:28px; height:22px;
         border-radius:5px; font-size:13px; font-weight:700; box-sizing:border-box; color:#c4c7cd; cursor:default; }
.vt td.ev { white-space:normal; line-height:1.3; padding-left:12px; font-size:14.5px; }
.vt .dim { color:#8b8f97; }
.vt tr.sec td { height:auto; padding:16px 8px 6px; font-size:12px; font-weight:700; color:#6b7280;
                letter-spacing:.03em; text-transform:uppercase; border-bottom:1px solid #e6e8ec; }
.vt tr.grp td { height:30px; cursor:pointer; color:#31333f; font-size:14px; background:#f7f8fa; }
.vt tr.grp td b { font-weight:600; }
.vt tr.grp .arr { display:inline-block; width:14px; color:#6b7280; }
.vt tr.hid { display:none; }
.vt tr.note td { height:auto; padding:6px 8px; color:#6b7280; font-size:13.5px; white-space:normal; }
.vt .mono { font-family:ui-monospace,Menlo,monospace; font-size:13px; }
.lg { font-size:13.5px; color:#4b5563; margin:12px 0 2px; line-height:1.9; }
.lg .row { display:flex; flex-wrap:wrap; gap:2px 18px; align-items:center; }
.lg .t { color:#6b7280; font-weight:600; font-size:12px; text-transform:uppercase; letter-spacing:.03em;
         min-width:52px; }
.lg span.i { display:inline-flex; align-items:center; gap:6px; }
.lg .c { display:inline-flex; align-items:center; justify-content:center; width:24px; height:19px;
         border-radius:4px; font-size:12px; font-weight:700; box-sizing:border-box; }
.lg .d { display:inline-block; width:11px; height:11px; border-radius:50%; box-sizing:border-box; }
</style>
"""


def _e(s) -> str:
    return html.escape(str(s), quote=True)


def city_code(city: str) -> str:
    m = re.search(r"\(([A-Z]{2})\)", city or "")
    return m.group(1) if m else (city or "?")[:2].upper()


def city_name(city: str) -> str:
    return (city or "").split("(")[0].strip()


def _fill(bg, fg):
    return f"background:{bg};color:{fg};"


def _dashed(color):
    return f"background:#fff;border:1.5px dashed {color};color:{color};"


# cell styles, shared with the legend
S_PRESENT = _fill(GREEN, "#fff")
S_NOTFOUND = _fill(AMBER, "#422006")
S_DISAGREE = _fill(VIOLET, "#fff")
S_NOMARK = _fill("#e5e7eb", "#4b5563")
S_NODATA = _dashed(GREY)
S_FOUND = _fill(RED, "#fff")
S_ONEDAY = _dashed(RED)
S_NAMED = _fill("#ede9fe", VIOLET)
S_NOVEL = _fill(BLUE, "#fff")
S_NOVEL1 = _dashed(BLUE)
S_BROAD = _fill("#e2e8f0", "#334155")


def _cell(label, style, tip, sel) -> str:
    return (f"<td class='cc{' sel' if sel else ''}'><span class='c' style='{style}' "
            f"data-tip='{_e(tip)}'>{_e(label)}</span></td>")


def _empty(tip, sel):
    return _cell("–", "", tip, sel)


def _range(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return "–"
    lo, hi = min(vals), max(vals)
    return f"{lo}" if lo == hi else f"{lo}–{hi}"


def _pl(n, word):
    return word if n == 1 else word + "s"


# ── per-row pieces ─────────────────────────────────────────────────────────

def _panel_cells(v, per_city, cities, sel):
    out = []
    for c in cities:
        d = per_city.get(c)
        if not d:
            out.append(_empty(f"{city_name(c)}: no result", c == sel))
            continue
        fill, fg, sym, lab, _ = _PANEL.get(d["state"], _PANEL["no_data"])
        style = _fill(fill, fg) if fill else S_NODATA
        out.append(_cell(sym, style, f"{v} · {city_name(c)}: {lab} — {d.get('reason', '')}",
                         c == sel))
    return "".join(out)


def _panel_evidence(per_city):
    n_mk = max((d.get("n_markers", 0) for d in per_city.values()), default=0)
    if n_mk == 0:
        return "<span class='dim'>no ★ marker — can't be checked on its own</span>"
    meas = [d.get("n_present") for d in per_city.values() if d.get("n_measured")]
    if not meas:
        return f"<span class='dim'>too few reads on its {n_mk} {_pl(n_mk, 'marker')}</span>"
    return f"{_range(meas)} of {n_mk} {_pl(n_mk, 'marker')} present"


def _panel_tint(per_city):
    states = [d["state"] for d in per_city.values()]
    if "confirmed" in states:
        return GREEN
    if "not_found" in states:
        return AMBER
    if "inconsistent" in states:
        return VIOLET
    return None


def _finding_cells(v, f, cities, sel):
    pc = f.get("per_city", {}) or {}
    out = []
    for c in cities:
        d = pc.get(c)
        if d is None:
            out.append(_empty(f"{v} · {city_name(c)}: not seen", c == sel))
            continue
        n = int(d.get("days", 0) or 0)
        stars, regs = d.get("stars") or [], d.get("regions") or []
        if n == 0:
            out.append(_cell("0", S_NAMED, f"{v} · {city_name(c)}: named only — its mutations "
                             "are shared with related lineages, no ★ marker on reads", c == sel))
            continue
        tip = (f"{v} · {city_name(c)}: evidence on {n} {_pl(n, 'day')} · "
               + (f"{len(stars)} ★ {_pl(len(stars), 'marker')} ({', '.join(stars[:5])})"
                  if stars else "no ★ marker (combinations only)")
               + (f" · {len(regs)} {_pl(len(regs), 'region')}: "
                  + ", ".join(f"{lo:,}–{hi:,}" for lo, hi in regs[:5]) if regs else ""))
        out.append(_cell(n, S_FOUND if n >= 2 else S_ONEDAY, tip, c == sel))
    return "".join(out)


def _finding_evidence(f):
    pc = f.get("per_city", {}) or {}
    ev = [d for d in pc.values() if int(d.get("days", 0) or 0) >= 1]
    if f.get("status") == "named" or not ev:
        return "<span class='dim'>named only — shared mutations, no ★ marker</span>"
    conf = f.get("status") == "confirmed"
    n_star = [len(d.get("stars") or []) for d in ev]
    n_reg = [len(d.get("regions") or []) for d in ev]
    weak = conf and (max(n_star) == 0 or max(n_reg) <= 1)
    star = (f"{_range(n_star)} ★ {_pl(max(n_star), 'marker')}" if max(n_star)
            else "combinations only")
    txt = f"{star} · {_range(n_reg)} {_pl(max(n_reg), 'region')}"
    return f"<b>{txt} ⚠</b>" if weak else txt


def _tree_cell(r, color, dot_kind, name_html, tip):
    X0, IND = 10, 18
    d = r["depth"]
    parts = []
    for k, g in enumerate(r["guides"]):
        if g:
            parts.append(f"<span class='v' style='left:{X0 + k * IND}px'></span>")
    if d >= 1:
        x = X0 + (d - 1) * IND
        parts.append(f"<span class='v{' half' if r['last'] else ''}' style='left:{x}px'></span>")
        parts.append(f"<span class='h' style='left:{x}px;width:{IND - 4}px'></span>")
    if r["has_children"]:
        parts.append(f"<span class='v down' style='left:{X0 + d * IND}px'></span>")
    rad = 3 if dot_kind == "spine" else 6
    dot = {"spine": "background:#d3d1c7;",
           "filled": f"background:{color};",
           "dashed": f"background:#fff;border:1.8px dashed {color};",
           "ring": f"background:#fff;border:2px solid {color};"}[dot_kind]
    parts.append(f"<span class='dot' style='left:{X0 + d * IND - rad + 1}px;width:{2 * rad}px;"
                 f"height:{2 * rad}px;margin-top:-{rad}px;{dot}' data-tip='{_e(tip)}'></span>")
    pad = X0 + d * IND + rad + 8
    return (f"<td class='tr' style='padding-left:{pad}px'>{''.join(parts)}"
            f"<span data-tip='{_e(tip)}'>{name_html}</span></td>")


def _group(gid, label, n_cols, rows, open_=False):
    """A clickable row that shows / hides the rows under it."""
    if not rows:
        return []
    arr = "▾" if open_ else "▸"
    head = (f"<tr class='grp' data-grp='{gid}'><td colspan='{n_cols}'>"
            f"<span class='arr'>{arr}</span>{label}</td></tr>")
    body = [r.replace("<tr", f"<tr class='g-{gid}{'' if open_ else ' hid'}'", 1) for r in rows]
    return [head] + body


def _pattern_row(muts, days, cities, sel, style_full, style_one, color, evidence, tip_word):
    cells = []
    for c in cities:
        n = days.get(c)
        if n is None:
            cells.append(_empty(f"{city_name(c)}: not seen", sel == c))
        elif n == 0:
            cells.append(_cell("0", S_NOMARK, f"{city_name(c)}: seen, but on no day above "
                               "the evidence rule (20 reads and 0.5% of the day)", sel == c))
        else:
            cells.append(_cell(n, style_full if n >= 2 else style_one,
                               f"{city_name(c)}: {tip_word} on {n} {_pl(n, 'day')} — "
                               f"{' '.join(muts)}", sel == c))
    m = " ".join(muts[:4]) + (" …" if len(muts) > 4 else "")
    return (f"<tr><td style='padding-left:14px'><span class='mono' style='color:{color}' "
            f"data-tip='{_e(' '.join(muts))}'>{_e(m)}</span></td>{''.join(cells)}"
            f"<td class='ev'>{evidence}</td></tr>")


# ── the whole view ─────────────────────────────────────────────────────────

def build(cities, sel, rows, verdicts, findings, current_panel, ot=(),
          novel=(), broad=(), novel_rest=None) -> str:
    """cities: city names in column order; sel: the chosen city.
    rows: components.abundance_cooc_tree.tree_rows(...).
    verdicts: {city: {variant: {state, reason, n_present, n_measured, n_markers}}}.
    findings: {node: {status: confirmed|1 day|named, per_city: {city: {days,
      stars, regions}}, addable: bool}}.
    current_panel: variants in the panel now (an added finding shows "added").
    novel: [(mutations, {city: days})]; broad: [(mutations, {city: days},
      n_lineages, ancestor)]; novel_rest: reads text of novel patterns not listed."""
    ot = set(ot)
    n_c = len(cities)
    n_cols = n_c + 2
    head = ("<tr><th>Variant</th>"
            + "".join(f"<th class='cc{' sel' if c == sel else ''}'><button data-city='{_e(c)}' "
                      f"data-tip='{_e(city_name(c))} — click to colour the tree by this city'>"
                      f"{_e(city_code(c))}</button></th>" for c in cities)
            + "<th style='padding-left:12px'>Evidence</th></tr>")
    body = []
    sname = city_name(sel)
    for r in rows:
        v, kind = r["node"], r["kind"]
        chips = " <span class='chip'>OT</span>" if v in ot and kind not in ("spine", "group") else ""
        if r.get("parents"):
            chips += (f"<span class='par' data-tip='{_e(r.get('parents_tip', ''))}'>"
                      f"{_e(r['parents'])}</span>")
        if kind == "group":
            body.append(f"<tr class='sp'>{_tree_cell(r, GREY, 'spine', f'''<span class='grp-root'>{_e(r['label'])}</span>''', r.get('tip', r['label']))}"
                        f"<td colspan='{n_c + 1}'></td></tr>")
            continue
        if kind == "spine":
            name = f"<span class='spl'>{_e(r['label'])}</span>{chips}"
            body.append(f"<tr class='sp'>{_tree_cell(r, GREY, 'spine', name, r['label'])}"
                        f"<td colspan='{n_c + 1}'></td></tr>")
            continue
        if kind in ("panel", "panel_ot"):
            pc = {c: verdicts[c][v] for c in cities if v in verdicts.get(c, {})}
            st_ = (pc.get(sel) or {}).get("state", "no_data")
            fill, _fg, _sym, lab, ncol = _PANEL.get(st_, _PANEL["no_data"])
            tip = f"{v} · {sname}: {lab} — {(pc.get(sel) or {}).get('reason', 'no result')}"
            name = f"<span class='nm' style='color:{ncol}'>{_e(v)}</span>{chips}"
            dcol = GREY if st_ in ("no_marker", "no_data") else fill
            tc = _tree_cell(r, dcol, "dashed" if st_ == "no_data" else "filled", name, tip)
            tint = _panel_tint(pc)
            cells, ev = _panel_cells(v, pc, cities, sel), _panel_evidence(pc)
        elif kind == "finding" and v in findings:
            f = findings[v]
            n_sel = int(((f.get("per_city") or {}).get(sel) or {}).get("days", 0) or 0)
            named = f.get("status") == "named"
            col = VIOLET if named else RED
            dk = ("ring" if named or n_sel == 0 else "filled" if n_sel >= 2 else "dashed")
            tip = (f"{v}: not in your panel · {sname}: "
                   + ("named only" if named else f"evidence on {n_sel} {_pl(n_sel, 'day')}"
                      if n_sel else "not seen here"))
            btn = ""
            if f.get("addable"):
                btn = ("<span class='added'>added · re-run</span>" if v in current_panel
                       else f"<button class='add' data-add='{_e(v)}'>+ Add</button>")
            name = (f"<span class='nm' style='color:{col}'>{_e(v)}</span>{chips}{btn}")
            tc = _tree_cell(r, col, dk, name, tip)
            tint = RED if f.get("status") == "confirmed" else None
            cells, ev = _finding_cells(v, f, cities, sel), _finding_evidence(f)
        else:  # tracked, not selected
            name = f"<span style='color:{TRACKED}'>{_e(v)}</span>{chips}"
            tc = _tree_cell(r, TRACKED, "ring", name, f"{v}: officially tracked, not in your panel")
            tint, cells, ev = None, f"<td colspan='{n_c}'></td>", "<span class='dim'>not in the run</span>"
        bg = f" style='background:{tint}14'" if tint else ""
        body.append(f"<tr{bg}>{tc}{cells}<td class='ev'>{ev}</td></tr>")

    # ── unnamed signal ──
    rec = [(m, d) for m, d in novel if max(d.values(), default=0) >= 2]
    one = [(m, d) for m, d in novel if max(d.values(), default=0) == 1]
    low = [(m, d) for m, d in novel if max(d.values(), default=0) == 0]
    if novel or broad or novel_rest:
        body.append(f"<tr class='sec'><td colspan='{n_cols}'>Unnamed signal</td></tr>")
        nov_ev = "<span class='dim'>novel — no pango lineage has this combination</span>"
        body += _group("nrec", f"<b style='color:{BLUE}'>Novel, recurring ({len(rec)})</b> "
                       "<span class='dim'>· on ≥ 2 days in a city</span>", n_cols,
                       [_pattern_row(m, d, cities, sel, S_NOVEL, S_NOVEL1, BLUE, nov_ev, "seen")
                        for m, d in rec], open_=True)
        body += _group("none", f"<b style='color:{BLUE}'>Novel, 1 day ({len(one)})</b>", n_cols,
                       [_pattern_row(m, d, cities, sel, S_NOVEL, S_NOVEL1, BLUE, nov_ev, "seen")
                        for m, d in one])
        body += _group("nlow", f"<b style='color:{BLUE}'>Novel, below the day rule ({len(low)})</b> "
                       "<span class='dim'>· never 20 reads and 0.5% of a day</span>", n_cols,
                       [_pattern_row(m, d, cities, sel, S_NOVEL, S_NOVEL1, BLUE, nov_ev, "seen")
                        for m, d in low])
        body += _group("broad", f"<b style='color:{SLATE}'>Too broad to name ({len(broad)})</b> "
                       "<span class='dim'>· the combination fits many unrelated lineages</span>",
                       n_cols,
                       [_pattern_row(m, d, cities, sel, S_BROAD, _dashed(SLATE), SLATE,
                                     f"<span class='dim'>fits {n} lineages"
                                     + (f" · mostly {_e(a)}" if a else "") + "</span>", "seen")
                        for m, d, n, a in broad])
        if novel_rest:
            body.append(f"<tr class='note'><td colspan='{n_cols}'>+ more small novel patterns "
                        f"({novel_rest} reads) — not listed</td></tr>")

    return (CSS + f"<table class='vt'><thead>{head}</thead><tbody>{''.join(body)}</tbody></table>"
            + LEGEND)


def _li(style, sym, text):
    return f"<span class='i'><i class='c' style='{style}'>{sym}</i>{text}</span>"


def _ld(style, text):
    return f"<span class='i'><i class='d' style='{style}'></i>{text}</span>"


LEGEND = (
    "<div class='lg'>"
    "<div class='row'><span class='t'>Panel</span>"
    + _li(S_PRESENT, "✓", "present")
    + _li(S_NOTFOUND, "✗", "not found")
    + _li(S_DISAGREE, "≠", "markers disagree")
    + _li(S_NOMARK, "·", "no ★ marker")
    + _li(S_NODATA, "?", "too few reads")
    + "</div><div class='row'><span class='t'>Found</span>"
    + _li(S_FOUND, "5", "not in panel · days with evidence")
    + _li(S_ONEDAY, "1", "1 day only")
    + _li(S_NAMED, "0", "named only (shared mutations)")
    + _li(S_NOVEL, "3", "novel · days")
    + _li(S_BROAD, "2", "too broad to name · days")
    + "</div><div class='row'><span class='t'>Tree</span>"
    + _ld(f"background:{RED}", "found in this city")
    + _ld(f"border:1.8px dashed {RED}", "1 day here")
    + _ld(f"border:2px solid {RED}", "not seen here")
    + _ld(f"border:2px solid {TRACKED}", "tracked, not selected")
    + "<span class='i'>panel dots = the cell colours</span>"
    + "<span class='i'>⚠ weak fact to judge · hover anything for details · "
      "click a city code to colour the tree</span></div></div>")


# ── the structural tree (left column): same rows and look, no evidence ──
_TREE_CSS = """
<style>
body { margin:0; padding:2px 2px 4px; background:transparent;
       font-family:"Source Sans Pro","Source Sans 3",-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
       color:#31333f; }
.vt { font-size:14px; }
.vt td { height:28px; border-bottom:none; }
.vt tr.sp td { height:21px; }
.vt .spl { font-size:11.5px; max-width:210px; }
.vt .par { font-size:11.5px; }
.vt .grp-root { font-size:11.5px; }
.lg { font-size:12px; margin-top:6px; line-height:1.6; }
#tip { position:fixed; z-index:10; display:none; max-width:300px; padding:5px 8px; border-radius:6px;
       background:#31333f; color:#fff; font-size:12px; pointer-events:none; }
</style>
"""

_TIP_JS = """
<div id="tip"></div>
<script>
var tip = document.getElementById("tip");
document.addEventListener("mousemove", function (e) {
  var t = e.target.closest ? e.target.closest("[data-tip]") : null;
  if (!t) { tip.style.display = "none"; return; }
  tip.textContent = t.getAttribute("data-tip"); tip.style.display = "block";
  var w = tip.offsetWidth, h = tip.offsetHeight, x = e.clientX + 12, y = e.clientY + 14;
  if (x + w > window.innerWidth - 4) x = Math.max(4, e.clientX - w - 12);
  if (y + h > window.innerHeight - 4) y = Math.max(4, e.clientY - h - 12);
  tip.style.left = x + "px"; tip.style.top = y + "px";
});
document.addEventListener("mouseleave", function () { tip.style.display = "none"; });
</script>
"""


def build_tree(rows, ot=()) -> tuple:
    """The panel tree alone (left column): (html, height in px)."""
    ot = set(ot)
    body, h = [], 0
    for r in rows:
        v, kind = r["node"], r["kind"]
        chips = " <span class='chip'>OT</span>" if v in ot and kind not in ("spine", "group") else ""
        if r.get("parents"):
            chips += (f"<span class='par' data-tip='{_e(r.get('parents_tip', ''))}'>"
                      f"{_e(r['parents'])}</span>")
        if kind == "group":
            name = f"<span class='grp-root'>{_e(r['label'])}</span>"
            body.append(f"<tr class='sp'>{_tree_cell(r, GREY, 'spine', name, r.get('tip', ''))}</tr>")
            h += 21
        elif kind == "spine":
            name = f"<span class='spl'>{_e(r['label'])}</span>{chips}"
            body.append(f"<tr class='sp'>{_tree_cell(r, GREY, 'spine', name, r['label'])}</tr>")
            h += 21
        elif kind in ("panel", "panel_ot", "finding"):
            name = f"<span class='nm'>{_e(v)}</span>{chips}"
            body.append(f"<tr>{_tree_cell(r, '#1f2430', 'filled', name, f'{v}: in your panel')}</tr>")
            h += 28
        else:
            name = f"<span style='color:{TRACKED}'>{_e(v)}</span>{chips}"
            body.append(f"<tr>{_tree_cell(r, TRACKED, 'ring', name, f'{v}: officially tracked, not selected')}</tr>")
            h += 28
    legend = ("<div class='lg'><span class='i'><i class='d' style='background:#1f2430'></i>in your panel</span> · "
              f"<span class='i'><i class='d' style='border:2px solid {TRACKED}'></i>tracked, not selected</span> · "
              "<span class='chip' style='margin:0'>OT</span> officially tracked</div>")
    html_ = (CSS + _TREE_CSS + f"<table class='vt'><tbody>{''.join(body)}</tbody></table>"
             + legend + _TIP_JS)
    return html_, h + 60


# ── Streamlit component ────────────────────────────────────────────────────
_FRONTEND = Path(__file__).parent / "variants_table_frontend"
_COMPONENT = None


def render(view_html: str, key: str):
    """Show the view; returns the last click, {"add": node, "t": ms} or
    {"city": name, "t": ms} (it stays the value on later reruns — compare t)."""
    global _COMPONENT
    if _COMPONENT is None:
        import streamlit.components.v1 as components
        _COMPONENT = components.declare_component("variants_table", path=str(_FRONTEND))
    return _COMPONENT(html=view_html, key=key, default=None)