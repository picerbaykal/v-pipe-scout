"""
components/variants_table.py

The "Variants" view of the Co-occurrence results page: the panel tree and
the per-city evidence in ONE table (2026-09-30).

  Variant   the tree, coloured by the evidence in the chosen city; the
            recombinant families (XEC, XFG, XDV → NB.1.8.1 …) under their
            own root, not under B; "+ Add" next to a confirmed finding
  cities    one cell per city — click a city code to colour the tree by it
  Evidence  for the chosen city: what was counted (markers voted, days,
            ★ markers, genome regions) and, in words, when a fact is weak

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
VIOLET = "#7c3aed"         # mixed markers (panel) · named only (finding)
GREY = "#9ca3af"           # no ★ marker / no data
RED = "#dc2626"            # found, not in panel
BLUE = "#2563eb"           # novel
SLATE = "#64748b"          # too broad to name
TRACKED = "#185FA5"

_PANEL = {  # check state -> (fill, text-on-fill, symbol, label, name colour)
    "confirmed": (GREEN, "#fff", "✓", "present", GREEN),
    "not_found": (AMBER, "#422006", "✗", "not found", AMBER_T),
    "inconsistent": (VIOLET, "#fff", "±", "mixed", VIOLET),
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
.vt .near-tag { font-size:12px; font-weight:600; margin-left:9px; padding:0 8px; line-height:19px;
                border-radius:10px; background:#dcfce7; color:#166534; cursor:pointer; vertical-align:1px; }
.vt .near-tag:hover { background:#bbf7d0; }
.vt .near-tag.off { background:#f3f4f6; color:#9ca3af; font-weight:500; }
.vt tr.near td { height:28px; }
.vt .nl { color:#166534; font-size:13.5px; }
.vt .nl .mono { font-size:13px; }
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
.vt .cal { display:inline-flex; gap:3px; vertical-align:-1px; margin-left:6px; }
.vt .cal i, .lg .cal i { display:inline-block; width:11px; height:11px; border-radius:2px; box-sizing:border-box; }
.vt .calw { font-size:12.5px; margin-left:5px; color:#6b7280; }
.lg { font-size:13.5px; color:#4b5563; margin:12px 0 2px; line-height:1.9; }
.lg .lgt { cursor:pointer; font-weight:600; color:#31333f; }
.lg .hid { display:none; }
.lg .row { display:flex; flex-wrap:wrap; gap:2px 18px; align-items:center; }
.lg .t { color:#6b7280; font-weight:600; font-size:12px; text-transform:uppercase; letter-spacing:.03em;
         min-width:52px; }
.lg span.i { display:inline-flex; align-items:center; gap:6px; }
.lg .c { font-style:normal; display:inline-flex; align-items:center; justify-content:center; width:24px; height:19px;
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
LIGHT_GREEN = "#4ade80"    # panel variant ± 1 change (as in the completeness graph)
S_NEAR = _fill("#bbf7d0", "#14532d")
S_NEAR1 = _dashed("#16a34a")


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
    if n == 1:
        return word
    if word.endswith("y") and word[-2:-1] not in "aeiou":
        return word[:-1] + "ies"           # city → cities, but day → days
    return word + "s"


# ── per-row pieces ─────────────────────────────────────────────────────────

# ── panel cells: a colour scale instead of cut-offs (2026-10-01) ──────────
# hue   = share of the MEASURED markers that are present: orange (0 %) →
#         pale grey (50 %) → green (100 %)
# depth = share of ALL its markers that could be measured: faint when few
#         were (2 of 5), full when all were
_SC_LO, _SC_MID, _SC_HI = (245, 158, 11), (209, 213, 219), (22, 163, 74)


def _hex(rgb):
    return "#%02x%02x%02x" % tuple(int(round(x)) for x in rgb)


def _mix(a, b, t):
    return tuple(a[i] + (b[i] - a[i]) * t for i in range(3))


_SC_GREY, _SC_RED = (229, 231, 235), (220, 38, 38)


def _scale(d, found=False):
    """(cell style, label, dot colour, name colour, share present or None)
    for one city's check result. found: a lineage not in the panel — grey
    (absent) → red (present) instead of orange → green; absent is no warning
    for a lineage you don't track."""
    n_mk, p, a, u = _vote_counts(d)
    if n_mk == 0:
        return S_NOMARK, "·", GREY, "#6b7280", None
    m = p + a
    if m == 0:
        return S_NODATA, "?", GREY, "#6b7280", None
    share = p / m
    if found:
        hue = _mix(_SC_GREY, _SC_RED, share)
        depth = 0.3 + 0.7 * (m / n_mk)
        bg = _mix((255, 255, 255), hue, depth)
        fg = "#fff" if share >= 0.7 and depth > 0.8 else "#111827"
        return _fill(_hex(bg), fg), f"{p}/{m}", _hex(bg), RED, share
    hue = _mix(_SC_LO, _SC_MID, share * 2) if share <= 0.5 else _mix(_SC_MID, _SC_HI, share * 2 - 1)
    depth = 0.3 + 0.7 * (m / n_mk)                 # 30 % strength at the least
    bg = _mix((255, 255, 255), hue, depth)
    fg = "#111827"                                 # always dark: one ink for all cells
    ink = "#15803d" if share >= 0.7 else AMBER_T if share <= 0.3 else "#4b5563"
    return _fill(_hex(bg), fg), f"{p}/{m}", _hex(bg), ink, share


def _panel_cells(v, per_city, cities, sel):
    out = []
    for c in cities:
        d = per_city.get(c)
        if not d:
            out.append(_empty(f"{city_name(c)}: no result", c == sel))
            continue
        style, label, _dot, _ink, share = _scale(d)
        out.append(_cell(label, style, f"{v} · {city_name(c)}: "
                         + (_breakdown(d, v, html=False) if d.get("n_markers") else _vote_text(d))
                         + _marker_details(d, v), c == sel))
    return "".join(out)


def _vote_counts(d):
    n_mk, n_meas, n_p = d.get("n_markers", 0), d.get("n_measured", 0), d.get("n_present", 0)
    n_meas = min(n_meas, n_mk)
    n_p = min(n_p, n_meas)
    return n_mk, n_p, n_meas - n_p, n_mk - n_meas


def _vote_text(d):
    """'1 present · 3 absent · 1 too few reads (of 5 ★ markers)' — exactly what
    the verdict voted on: present ÷ (present + absent); ≥ 75 % present,
    ≤ 25 % not found, in between mixed."""
    n_mk, p, a, u = _vote_counts(d)
    if n_mk == 0:
        return "no ★ marker — it can't be told apart from related lineages on its own"
    parts = [f"{p} present", f"{a} absent"] + ([f"{u} too few reads"] if u else [])
    return " · ".join(parts) + f" (of {n_mk} ★ {_pl(n_mk, 'marker')})"


def _marker_list(d):
    r = d.get("reason", "")
    return f" — {r.split(' — ', 1)[1]}" if " — " in r else ""


def _thr():
    """The check's marker thresholds, from the config (so the words follow it)."""
    try:
        from process.cooc import _check_cfg
        c = _check_cfg()
    except Exception:
        c = {"min_cov": 100, "present_freq": 0.05, "absent_freq": 0.01}
    return (int(c["min_cov"]), float(c["present_freq"]) * 100, float(c["absent_freq"]) * 100)


def _pct(x):
    return f"{x:g}"


def _breakdown(d, v, html=True):
    """'2 present (≥ 5 % of ≥ 100 reads) · 1 at 1–5 % · 1 on reads unlike X ·
    1 under 100 reads' — every marker in one of five parts, each part named
    by its rule. Uses the per-marker groups when the result has them."""
    mc, pf, af = _thr()
    mk = d.get("markers")
    if mk is None:                                   # results from before the groups
        n_mk, p, a, u = _vote_counts(d)
        cnt = {"present": p, "absent": a, "low": u}
    else:
        cnt = {}
        for x in mk:
            cnt[x["group"]] = cnt.get(x["group"], 0) + 1
    b = (lambda t: f"<b>{t}</b>") if html else (lambda t: t)
    parts = [(b(f"{cnt.get('present', 0)} present") if cnt.get("present") else "0 present")
             + f" (≥ {_pct(pf)} % of ≥ {mc} reads)"]
    if cnt.get("absent"):
        parts.append(b(f"{cnt['absent']} absent") + f" (< {_pct(af)} % of ≥ {mc} reads)")
    if cnt.get("between"):
        parts.append(f"{cnt['between']} at {_pct(af)}–{_pct(pf)} %")
    if cnt.get("nofit"):
        parts.append(f"{cnt['nofit']} on reads unlike {v}")
    if cnt.get("low"):
        parts.append(f"{cnt['low']} under {mc} reads")
    return " · ".join(parts)


def _marker_details(d, v):
    """Hover text: every marker with its numbers and what they mean."""
    mc, pf, af = _thr()
    mk = d.get("markers")
    if not mk:
        return _marker_list(d)
    why = {"present": "present", "absent": "absent",
           "between": f"at {_pct(af)}–{_pct(pf)} %: too much for noise, too little to count",
           "nofit": f"≥ {_pct(pf)} %, but on reads that don't look like {v} at its other "
                    "positions (another lineage carries it)",
           "low": f"under {mc} reads"}
    order = ["present", "absent", "between", "nofit", "low"]
    out = []
    for x in sorted(mk, key=lambda x: order.index(x["group"])):
        f = x.get("freq")
        num = (f"{f * 100:.1f} % of {x['cov']:,} reads" if f is not None and x["cov"]
               else "no reads")
        out.append(f"{x['marker']} {num} — {why[x['group']]}")
    more = f" · +{len(out) - 10} more" if len(out) > 10 else ""
    return " — " + " · ".join(out[:10]) + more


# ── the Evidence calendar (2026-10-02): one square per sampling day of the
#    chosen city, for every row ─────────────────────────────────────────────
def _day_marks(top):
    """mark -> (square style, meaning). top = colour of a present day."""
    return {
        # lineages (★ markers pooled per day, process.variant_explorer.check_in_data)
        "present": (f"background:{top};", "★ markers present (≥ 5 % of ≥ 100 reads)"),
        "weak": (f"background:{top}55;", "★ markers between absent and present"),
        "absent": ("background:#d1d5db;", "★ markers absent (< 1 % of ≥ 100 reads)"),
        # novel combinations (scanner timeline)
        "day": (f"background:{top};", "passes the day rule (≥ 20 reads, ≥ 0.5 % of the day)"),
        "seen": (f"background:{top}55;", "seen, below the day rule"),
        "gone": ("background:#d1d5db;", "≥ 100 reads at its positions, not there"),
        "uncovered": ("background:#fff;border:1.5px dashed #9ca3af;",
                      "under 100 reads at its positions — couldn't show"),
        "unknown": ("background:#e5e7eb;", "no coverage data (re-run)"),
    }


def _recent_n():
    try:
        from process.cooc import _check_cfg
        return int(_check_cfg().get("recent_samples", 5))
    except Exception:
        return 5


def _day_cal(tl, top):
    """The last N covered samples only (check.recent_samples): what the cell
    voted on — "is it there now?". Days with too few reads are skipped."""
    n = _recent_n()
    last = [(d, m) for d, m in (tl or []) if m not in ("uncovered", "unknown")][-n:]
    if not last:
        return ""
    mk = _day_marks(top)
    sq = "".join(f"<i style='{mk.get(m, mk['unknown'])[0]}' "
                 f"data-tip='{_e(_short_date(d) + ': ' + mk.get(m, mk['unknown'])[1])}'></i>"
                 for d, m in last)
    pad = "".join("<i style='background:transparent'></i>" for _ in range(n - len(last)))
    return f"<span class='cal' style='margin-left:0;margin-right:8px'>{pad}{sq}</span>"


def _cal_phrase(tl, word="present"):
    """'present on 4 of 9 covered days · in the latest sample' — what the
    calendar says, in a few words."""
    if not tl:
        return ""
    on = {"present", "day"}
    cov = [(d, m) for d, m in tl if m not in ("uncovered", "unknown")]
    hits = [d for d, m in cov if m in on]
    if not cov:
        return "not covered — under 100 reads at its positions"
    if not hits:
        weak = sum(1 for _d, m in cov if m in ("weak", "seen"))
        return (f"absent all {len(cov)} {_pl(len(cov), 'day')}" if not weak
                else f"weak {weak} of {len(cov)} {_pl(len(cov), 'day')}, never present")
    after = sum(1 for d, _m in cov if d > hits[-1])
    tail = "in the latest sample" if not after else f"last {_short_date(hits[-1])}"
    return f"{word} {len(hits)} of {len(cov)} {_pl(len(cov), 'day')} · {tail}"


def _panel_evidence(per_city, sel, v=""):
    d = per_city.get(sel)
    if not d:
        return "<span class='dim'>no result in this city</span>"
    n_mk, p, a, u = _vote_counts(d)
    if n_mk == 0:
        return f"<span class='dim'>{_vote_text(d)}</span>"
    tip = (f"{v} · the cell votes on the last {_recent_n()} samples with enough reads: "
           f"present / (present + absent) = {p}/{p + a} — "
           + _breakdown(d, v, html=False)
           + ". Markers that are neither present nor absent don't count either way"
           + _marker_details(d, v))
    tl = d.get("timeline")
    if not tl:                                      # results from before 2026-10-02
        return f"<span data-tip='{_e(tip)}'>{_breakdown(d, v)}</span>"
    return f"{_day_cal(tl, GREEN)}<span data-tip='{_e(tip)}'>{_e(_cal_phrase(tl))}</span>"


def _panel_tint(per_city, sel):
    """Row tint: the chosen city's cell colour (none when not measured)."""
    d = per_city.get(sel)
    if not d:
        return None
    _style, _label, dot, _ink, share = _scale(d)
    return None if share is None else dot


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
        chk = d.get("check")
        if d.get("xcheck") and chk:
            # phase 3: not named by this city's scan, checked because it was
            # named elsewhere — same colours, dashed border
            style, label, _dot, _ink, _share = _scale(chk, found=True)
            style += "border:1.5px dashed #b91c1c;"
            fi = ", ".join(city_name(x) for x in d.get("found_in") or [])
            out.append(_cell(label, style, f"{v} · {city_name(c)}, last {_recent_n()} samples "
                             f"with enough reads: {_vote_text(chk)} · not named by this city's "
                             f"scan (its reads here may carry only one mutation beyond your "
                             f"panel); checked because it was found in {fi}", c == sel))
            continue
        if chk and n > 0:
            # the same measure as the panel: ★ markers present / measurable
            style, label, _dot, _ink, _share = _scale(chk, found=True)
            out.append(_cell(label, style, f"{v} · {city_name(c)}, last {_recent_n()} samples "
                             f"with enough reads: {_vote_text(chk)} · the scanner counted "
                             f"{n} {_pl(n, 'day')} of evidence in the window", c == sel))
            continue
        if n == 0:
            out.append(_cell("0", S_NAMED, f"{v} · {city_name(c)}: named only — reads point to "
                             f"{v}, but they also fit related lineages, so no day has evidence "
                             f"specific to {v}", c == sel))
            continue
        tip = (f"{v} · {city_name(c)}: evidence on {n} {_pl(n, 'day')} · "
               + (f"{len(stars)} ★ {_pl(len(stars), 'marker')} ({', '.join(stars[:5])})"
                  if stars else "no ★ marker (combinations only)")
               + (f" · {len(regs)} genome {_pl(len(regs), 'region')}: "
                  + ", ".join(f"{lo:,}–{hi:,}" for lo, hi in regs[:5]) if regs else ""))
        if n == 1 and d.get("after"):
            tip += " · 1 day " + _after_text(d["after"], html=False)
        out.append(_cell(n, S_FOUND if n >= 2 else S_ONEDAY, tip, c == sel))
    return "".join(out)


def _after_text(a, html=True):
    """A one-day finding: what happened after its day, in words."""
    if not a:
        return ""
    st, k, cov, seen = a.get("state"), a.get("later", 0), a.get("covered"), a.get("seen", 0)
    day = _short_date(a.get("day", ""))
    if st == "latest":
        t = f"on {day}, the latest sample — too early to tell, watch it"
    elif st == "seen_again":
        t = (f"on {day}; back on {seen} of {k} later {_pl(k, 'sample')}, but too weak "
             "to count as a second day")
    elif st == "not_covered":
        t = (f"on {day}; {k} later {_pl(k, 'sample')}, none with ≥ 100 reads at its "
             "positions — could not have shown up")
    elif cov is None:
        t = f"on {day}; not seen in {k} later {_pl(k, 'sample')}"
    else:
        t = (f"on {day}; not seen since — {cov} later {_pl(cov, 'sample')} covered its "
             "positions without it (likely a one-sample artefact)")
    if html and st in ("latest", "seen_again"):
        return f"<b>{t}</b>"
    return t


_CAL = {  # one-day finding, per sample from its day on: square style, meaning
    "day": (f"background:{RED};", "the day it was seen"),
    "seen": ("background:#fca5a5;", "back, but too weak to count as a day"),
    "gone": ("background:#9ca3af;", "covered (≥ 100 reads), not there"),
    "uncovered": ("background:#fff;border:1.5px dashed #9ca3af;", "too few reads at its positions"),
    "unknown": ("background:#e5e7eb;", "no coverage data (re-run)"),
}
_CAL_WORD = {"latest": "latest sample · watch", "seen_again": "back, weak",
             "not_covered": "not covered since", "not_seen": "gone since"}


def _after_cal(a):
    """A one-day finding as a tiny calendar: one square per sample from its
    day to the last sample, plus one or two words."""
    tl = (a or {}).get("timeline") or []
    if not tl:
        return ""
    sq = "".join(f"<i style='{_CAL.get(m, _CAL['unknown'])[0]}' "
                 f"data-tip='{_e(_short_date(d) + ': ' + _CAL.get(m, _CAL['unknown'])[1])}'></i>"
                 for d, m in tl)
    w = _CAL_WORD.get(a.get("state"), "")
    wst = " style='color:#b91c1c;font-weight:600'" if a.get("state") in ("latest", "seen_again") else ""
    return (f"<span class='cal' data-tip='{_e(_after_text(a, html=False))}'>{sq}</span>"
            f"<span class='calw'{wst}>{_e(w)}</span>")


def _short_date(d):
    """'2026-08-22' -> '22 Aug'."""
    try:
        from datetime import date
        x = date.fromisoformat(str(d)[:10])
        return f"{x.day} {x.strftime('%b')}"
    except Exception:
        return str(d)


_REGION_TIP = ("genome regions = separate stretches of the genome where evidence reads "
               "were seen (each read covers ~250 bases; overlapping reads merge into one "
               "region). More regions = independent pieces of the genome agree, so a "
               "PCR artefact or one convergent mutation is unlikely.")


def _finding_evidence(f, sel):
    pc = f.get("per_city", {}) or {}
    d = pc.get(sel)
    n_cities = sum(1 for x in pc.values() if int(x.get("days", 0) or 0) >= 1)
    if d is None:
        where = (f" · evidence in {n_cities} other {_pl(n_cities, 'city')}" if n_cities else "")
        return f"<span class='dim'>not seen in this city{where}</span>"
    n = int(d.get("days", 0) or 0)
    if d.get("xcheck"):
        tl = (d.get("check") or {}).get("timeline")
        fi = ", ".join(city_name(x) for x in d.get("found_in") or [])
        why = (f"not named by this city's scan; checked because it was found in {fi} — "
               "its ★ markers one by one, as for panel variants")
        return (f"{_day_cal(tl, RED)}<span data-tip='{_e(why)}'>{_e(_cal_phrase(tl))}"
                f" <span class='dim'>· cross-check</span></span>")
    if n == 0:
        return ("<span class='dim'>named only — reads point to it but also fit related "
                "lineages; no day with specific evidence</span>")
    stars, regs = d.get("stars") or [], d.get("regions") or []
    facts = (f"scanner: {n} {_pl(n, 'day')} · "
             + (f"{len(stars)} ★ {_pl(len(stars), 'marker')}" if stars else "combinations only")
             + f" · {len(regs)} genome {_pl(len(regs), 'region')}"
             + (" · 1 day " + _after_text(d["after"], html=False)
                if n == 1 and d.get("after") else ""))
    tl = (d.get("check") or {}).get("timeline")
    if tl:
        txt = (f"{_day_cal(tl, RED)}<span data-tip='{_e(facts + '. ' + _REGION_TIP)}'>"
               f"{_e(_cal_phrase(tl))}</span>")
    else:                                        # results from before 2026-10-02
        txt = (f"{n} {_pl(n, 'day')} · "
               + (f"{len(stars)} ★ {_pl(len(stars), 'marker')}" if stars else "combinations only")
               + f" · <span data-tip='{_e(_REGION_TIP)}'>{len(regs)} genome "
               f"{_pl(len(regs), 'region')}</span>")
        if n == 1 and d.get("after"):
            txt += _after_cal(d["after"]) or (" — " + _after_text(d["after"]))
    weak = []
    if not stars:
        weak.append("no ★ marker")
    if len(regs) <= 1:
        weak.append("one region only")
    if weak and n >= 2:
        txt += f" — <b>weak: {', '.join(weak)}</b>"
    return txt


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

def _near_rows(r, v, changes, cities, sel, gid, n_cols):
    """The folded "± 1 change" rows under panel variant v (hidden until its
    tag is clicked)."""
    d = r["depth"] + 1
    guides = (r["guides"] + [not r["last"]]) if r["depth"] >= 1 else []
    out = []
    for i, c in enumerate(changes):
        last = (i == len(changes) - 1) and not r["has_children"]
        pr = {"depth": d, "guides": guides, "last": last, "has_children": False}
        gain = c["sign"] == "+"
        name = (f"<span class='nl'>{_e(v)} {'+' if gain else '−'} "
                f"<span class='mono'>{_e(c['mut'])}</span></span>")
        what = (f"reads = {v} plus {c['mut']}" if gain
                else f"reads = {v} but the reference base at {c['mut'][:-1]} (no {c['mut']})")
        tc = _tree_cell(pr, LIGHT_GREEN, "dashed" if not gain else "ring", name,
                        f"{c['label']}: {what}. One mutation only — a hint to watch, not a finding.")
        cells = []
        for city in cities:
            n = c["per_city"].get(city)
            reads = c.get("reads", {}).get(city, 0)
            if not n:
                cells.append(_empty(f"{city_name(city)}: "
                                    + (f"{reads:,} reads, no day above the day rule" if reads
                                       else "not seen"), city == sel))
                continue
            cells.append(_cell(n, S_NEAR if n >= 2 else S_NEAR1,
                               f"{c['label']} · {city_name(city)}: {n} {_pl(n, 'day')} "
                               f"({reads:,} reads)", city == sel))
        n_sel = c["per_city"].get(sel, 0)
        w = c.get("where")
        if w:
            more = (f" and {w['n_roots'] - len(w['roots'])} more" if w["n_roots"] > len(w["roots"])
                    else "")
            place = (f"{'gained' if gain else 'lost'} in {', '.join(w['roots'])}{more} "
                     f"({w['n_lineages']} {_pl(w['n_lineages'], 'lineage')} with it)")
        else:
            place = (f"no designated sublineage of {v} {'carries' if gain else 'lacks'} it"
                     + (" — not designated yet?" if gain else
                        " — reversion, or reads at a coverage edge"))
        ev = (f"{n_sel} {_pl(n_sel, 'day')} here · " if n_sel else
              "<span class='dim'>not in this city · </span>") + place
        if not gain:
            ev = f"<span class='dim'>{ev}</span>"
        out.append(f"<tr class='near g-{gid} hid'>{tc}{''.join(cells)}<td class='ev'>{ev}</td></tr>")
    return out


def build(cities, sel, rows, verdicts, findings, current_panel, ot=(),
          novel=(), broad=(), novel_rest=None, near=None, near_min_days=2,
          novel_hotspot=None, data_until=None, novel_info=None) -> str:
    """cities: city names in column order; sel: the chosen city.
    rows: components.abundance_cooc_tree.tree_rows(...).
    verdicts: {city: {variant: {state, reason, n_present, n_measured, n_markers}}}.
    findings: {node: {status: confirmed|1 day|named, per_city: {city: {days,
      stars, regions}}, addable: bool}}.
    current_panel: variants in the panel now (an added finding shows "added").
    novel: [(mutations, {city: days})]; broad: [(mutations, {city: days},
      n_lineages, ancestor)]; novel_rest: reads text of novel patterns not listed.
    novel_hotspot: (patterns, reads text) left out at error-hotspot positions
      (scanner._drop_hotspots), shown as one grey line.
    novel_info: {tuple(mutations): {city: {timeline, clue}}} — the Evidence
      calendar and clues of each novel row (scanner._novel_groups / _novel_clues).
    near: {panel variant: [change, …]} from process.near_changes (each with
      "where"); a variant with some gets a tag "◆ N changes" (gains counting
      on >= near_min_days days) that unfolds them."""
    near = near or {}
    data_until = data_until or {}
    ot = set(ot)
    n_c = len(cities)
    n_cols = n_c + 2
    head = ("<tr><th>Variant</th>"
            + "".join(f"<th class='cc{' sel' if c == sel else ''}'><button data-city='{_e(c)}' "
                      f"data-tip='{_e(city_name(c))}"
                      + (f" · data until {_e(data_until[c])}" if c in data_until else "")
                      + f" — click to colour the tree by this city'>"
                      f"{_e(city_code(c))}</button></th>" for c in cities)
            + f"<th style='padding-left:12px'>Evidence · {_e(city_name(sel))}"
            + (f" <span style='font-weight:400;text-transform:none;letter-spacing:0'>· data until "
               f"{_e(data_until[sel])}</span>" if sel in data_until else "")
            + "</th></tr>")
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
            dsel = pc.get(sel) or {}
            _st, _lb, dcol, ncol, share = _scale(dsel) if dsel else (None, "", GREY, "#6b7280", None)
            tip = (f"{v} · {sname}: " + ((_breakdown(dsel, v, html=False)
                                          if dsel.get("n_markers") else _vote_text(dsel))
                                         if dsel else "no result"))
            name = f"<span class='nm' style='color:{ncol}'>{_e(v)}</span>{chips}"
            no_meas = share is None and dsel.get("n_markers", 0) > 0
            tc = _tree_cell(r, dcol, "dashed" if (no_meas or not dsel) else "filled", name, tip)
            tint = _panel_tint(pc, sel)
            cells, ev = _panel_cells(v, pc, cities, sel), _panel_evidence(pc, sel, v)
        elif kind == "finding" and v in findings:
            f = findings[v]
            _dsel = (f.get("per_city") or {}).get(sel) or {}
            n_sel = int(_dsel.get("days", 0) or 0)
            named = f.get("status") == "named"
            col = VIOLET if named else RED
            dk = ("ring" if named or n_sel == 0 else "filled" if n_sel >= 2 else "dashed")
            if _dsel.get("xcheck"):
                dk = "dashed" if (_dsel.get("check") or {}).get("state") == "present" else "ring"
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
            cells, ev = _finding_cells(v, f, cities, sel), _finding_evidence(f, sel)
        else:  # tracked, not selected
            name = f"<span style='color:{TRACKED}'>{_e(v)}</span>{chips}"
            tc = _tree_cell(r, TRACKED, "ring", name, f"{v}: officially tracked, not in your panel")
            tint, cells, ev = None, f"<td colspan='{n_c}'></td>", "<span class='dim'>not in the run</span>"
        bg = f" style='background:{tint}{'22' if kind in ('panel', 'panel_ot') else '14'}'" if tint else ""
        nch = near.get(v) if kind in ("panel", "panel_ot") else None
        if nch:
            gid = f"near{len(body)}"
            # the tag follows the chosen city, like the Evidence column
            nch = sorted(nch, key=lambda c: -c["per_city"].get(sel, 0))
            here = [c for c in nch if c["per_city"].get(sel)]
            n_tag = sum(1 for c in here if c["sign"] == "+"
                        and c["per_city"][sel] >= near_min_days)
            if here:
                tip = (f"{v} ± 1 change in {sname} — reads that are {v} with one mutation "
                       "more (+) or one less (−), on days with ≥ 20 such reads and ≥ 0.5 % "
                       "of the day: " + "; ".join(f"{c['label']} ({c['per_city'][sel]} d)"
                                                  for c in here[:8]))
            else:
                tip = (f"No ± 1 change on any day in {sname}; in other cities: "
                       + "; ".join(f"{c['label']}" for c in nch[:8]))
            label = (f"◆ {n_tag} {_pl(n_tag, 'change')}" if n_tag
                     else f"◆ {len(here)} weak {_pl(len(here), 'change')}" if here
                     else "◆ elsewhere")
            cls = "near-tag" + ("" if here else " off")
            tag = (f"<span class='{cls}' data-grp='{gid}' data-tip='{_e(tip)}'>"
                   f"{label} <span class='arr'>▸</span></span>")
            tc = tc.replace("</span></td>", f"</span>{tag}</td>", 1) if tc.endswith("</span></td>") else tc
        body.append(f"<tr{bg}>{tc}{cells}<td class='ev'>{ev}</td></tr>")
        if nch:
            body += _near_rows(r, v, nch[:8], cities, sel, gid, n_cols)

    # ── unnamed signal ──
    rec = [(m, d) for m, d in novel if max(d.values(), default=0) >= 2]
    one = [(m, d) for m, d in novel if max(d.values(), default=0) == 1]
    low = [(m, d) for m, d in novel if max(d.values(), default=0) == 0]
    if novel or broad or novel_rest or novel_hotspot:
        body.append(f"<tr class='sec'><td colspan='{n_cols}'>Unnamed signal</td></tr>")
        novel_info = novel_info or {}

        def nov_ev_for(m, d):
            info = (novel_info.get(tuple(m)) or {}).get(sel)
            if not info:
                n_other = sum(1 for c, x in d.items() if c != sel)
                return (f"<span class='dim'>not in this city"
                        + (f" · in {n_other} other {_pl(n_other, 'city')}" if n_other else "")
                        + "</span>")
            return _novel_evidence(m, info)
        body += _group("nrec", f"<b style='color:{BLUE}'>Novel, recurring ({len(rec)})</b> "
                       "<span class='dim'>· on ≥ 2 days in a city</span>", n_cols,
                       [_pattern_row(m, d, cities, sel, S_NOVEL, S_NOVEL1, BLUE, nov_ev_for(m, d), "seen")
                        for m, d in rec], open_=True)
        body += _group("none", f"<b style='color:{BLUE}'>Novel, 1 day ({len(one)})</b>", n_cols,
                       [_pattern_row(m, d, cities, sel, S_NOVEL, S_NOVEL1, BLUE, nov_ev_for(m, d), "seen")
                        for m, d in one])
        body += _group("nlow", f"<b style='color:{BLUE}'>Novel, below the day rule ({len(low)})</b> "
                       "<span class='dim'>· never 20 reads and 0.5% of a day</span>", n_cols,
                       [_pattern_row(m, d, cities, sel, S_NOVEL, S_NOVEL1, BLUE, nov_ev_for(m, d), "seen")
                        for m, d in low])
        body += _group("broad", f"<b style='color:{SLATE}'>Too broad to name ({len(broad)})</b> "
                       "<span class='dim'>· the combination fits many unrelated lineages</span>",
                       n_cols,
                       [_pattern_row(m, d, cities, sel, S_BROAD, _dashed(SLATE), SLATE,
                                     f"<span class='dim'>fits {n} lineages"
                                     + (f" · mostly {_e(a)}" if a else "") + "</span>", "seen")
                        for m, d, n, a in broad])
        if novel_hotspot:
            _tip = ("A position where the city's data shows two or more different new bases "
                    "(each on ≥ 1 % of ≥ 100 reads) is an error hotspot — a virus has one. "
                    "A mutation there that no pango lineage carries is left out of novel "
                    "patterns; a pattern left with one mutation, or one a lineage carries, "
                    "is not novel. Lineage mutations are always kept.")
            body.append(f"<tr class='note'><td colspan='{n_cols}'><span data-tip='{_e(_tip)}'>"
                        f"{novel_hotspot[0]:,} patterns ({novel_hotspot[1]} reads) at error-hotspot "
                        "positions — left out ⓘ</span></td></tr>")
        if novel_rest:
            body.append(f"<tr class='note'><td colspan='{n_cols}'>+ more small novel patterns "
                        f"({novel_rest} reads, summed over cities and days) — not listed</td></tr>")

    return (CSS + f"<table class='vt'><thead>{head}</thead><tbody>{''.join(body)}</tbody></table>"
            + LEGEND)


def _novel_evidence(muts, info):
    """A novel row in the chosen city: its calendar, then the clues that it may
    be a real mutation + a steady misread rather than a new variant."""
    tl = info.get("timeline") or []
    cl = info.get("clue") or {}
    words = []
    a, lin, n_p, tog = cl.get("anchor"), cl.get("lineage"), cl.get("partners"), cl.get("together")
    n_l = cl.get("n_lineages") or 0
    if a and n_l == 1 and lin:
        words.append(f"{a} is {lin}'s")
    elif a and n_l > 1:
        words.append(f"{a} is in {n_l} lineages" + (f", mostly under {lin}" if lin else ""))
    if a and n_p and n_p >= 2:
        words.append(f"{a} seen with {n_p} partners")
    if tog is not None and len(muts) >= 2:
        words.append(f"together on {tog * 100:.0f} % of covering reads"
                     if tog >= 0.01 else "together on < 1 % of covering reads")
    pf = _thr()[1] / 100
    hint = ""
    if n_l and tog is not None and tog < pf:
        hint = " → likely a real mutation + misread"
    elif tog is not None and tog >= 0.5:
        hint = " → the mutations travel together"
    tip = ("Clues, nothing hidden: a real new combination is on most reads covering its "
           "positions and keeps the same partners; a real lineage mutation plus a steady "
           "misread is on few of them (around the error rate) and changes partners. "
           f"Hint when together < {_pct(pf * 100)} % (the check's present share).")
    phrase = _cal_phrase(tl, word="passes") if tl else ""
    clue = " · ".join(words)
    return (f"{_day_cal(tl, BLUE)}<span class='dim'>{_e(phrase)}</span>"
            + (f"<br><span data-tip='{_e(tip)}'>{_e(clue)}<b>{_e(hint)}</b></span>" if clue else ""))


def _li(style, sym, text):
    return f"<span class='i'><i class='c' style='{style}'>{sym}</i>{text}</span>"


def _ld(style, text):
    return f"<span class='i'><i class='d' style='{style}'></i>{text}</span>"


_LEGEND_ROWS = (
    ""
    "<div class='row'><span class='t'>Panel</span>"
    "<span class='i'>markers present / measured:"
    + "".join(f"<i class='c' style='width:auto;padding:0 6px;{_scale(d)[0]}'>{_scale(d)[1]}</i>"
              for d in ({"n_markers": 4, "n_measured": 4, "n_present": 0},
                        {"n_markers": 4, "n_measured": 4, "n_present": 2},
                        {"n_markers": 4, "n_measured": 4, "n_present": 4},
                        {"n_markers": 5, "n_measured": 2, "n_present": 2}))
    + " the last: only 2 of 5 markers measured, so faint</span>"
    + _li(S_NOMARK, "·", "no ★ marker")
    + _li(S_NODATA, "?", "too few reads")
    + "</div><div class='row'><span class='t'>Found</span>"
    + "<span class='i'>not in your panel, same check:"
    + "".join(f"<i class='c' style='width:auto;padding:0 6px;{_scale(d, found=True)[0]}'>"
              f"{_scale(d, found=True)[1]}</i>"
              for d in ({"n_markers": 3, "n_measured": 3, "n_present": 0},
                        {"n_markers": 3, "n_measured": 3, "n_present": 3}))
    + "</span>"
    + "<span class='i'><i class='c' style='width:auto;padding:0 6px;"
    + _scale({"n_markers": 2, "n_measured": 1, "n_present": 1}, found=True)[0]
    + "border:1.5px dashed #b91c1c;'>1/1</i>dashed = not named by this city's scan, "
      "checked because found in another city</span>"
    + _li(S_FOUND, "5", "older results: days with evidence")
    + _li(S_NAMED, "0", "named only (shared mutations)")

    + _li(S_NOVEL, "3", "novel · days")
    + _li(S_BROAD, "2", "too broad to name · days")
    + _li(S_NEAR, "3", "◆ panel variant ± 1 change · days (a hint)")
    + "</div><div class='row'><span class='t'>Days</span>"
    + "<span class='i'>Evidence, one square per sampling day in the chosen city:"
    + "<span class='cal'>" + "".join(
        f"<i style='{_day_marks(GREEN)[k][0]}'></i>" for k in ("present", "weak", "absent", "uncovered"))
    + "</span> present · weak · absent · under 100 reads (red for found, blue for novel)</span>"
    + "</div><div class='row'><span class='t'>Tree</span>"
    + _ld(f"background:{RED}", "found in this city")
    + _ld(f"border:1.8px dashed {RED}", "1 day here")
    + _ld(f"border:2px solid {RED}", "not seen here")
    + _ld(f"border:2px solid {TRACKED}", "tracked, not selected")
    + "<span class='i'>panel dots = the cell colours</span>"
    + "<span class='i'>weak = no ★ marker or one genome region only · hover anything for details · "
      "click a city code to colour the tree</span></div>")

LEGEND = ("<div class='lg'><div class='row'><span class='lgt' data-grp='legend'>"
          "<span class='arr'>▸</span> How to read the colours</span>"
          "<span class='dim'>· or hover anything to see what it means</span></div>"
          "<div class='g-legend hid'>" + _LEGEND_ROWS + "</div></div>")


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


# Click behaviour, sent with every render (index.html runs it once per
# change): "+ Add", a city code, and unfolding a group / "± 1 change" tag.
CLICK_JS = """
window.vtClick = function (e, send, sendHeight) {
  var b = e.target.closest("[data-add]");
  if (b) {
    b.disabled = true; b.textContent = "adding…";
    send("streamlit:setComponentValue",
         {value: {add: b.getAttribute("data-add"), t: Date.now()}, dataType: "json"});
    return;
  }
  var c = e.target.closest("[data-city]");
  if (c) {
    send("streamlit:setComponentValue",
         {value: {city: c.getAttribute("data-city"), t: Date.now()}, dataType: "json"});
    return;
  }
  var g = e.target.closest("[data-grp]");
  if (g) {
    var id = g.getAttribute("data-grp"), open = false;
    document.querySelectorAll(".g-" + id).forEach(function (r) {
      open = r.classList.toggle("hid") === false;
    });
    var arr = g.querySelector(".arr");
    if (arr) arr.textContent = open ? "▾" : "▸";
    sendHeight();
  }
};
"""


def render(view_html: str, key: str):
    """Show the view; returns the last click, {"add": node, "t": ms} or
    {"city": name, "t": ms} (it stays the value on later reruns — compare t)."""
    global _COMPONENT
    if _COMPONENT is None:
        import streamlit.components.v1 as components
        _COMPONENT = components.declare_component("variants_table", path=str(_FRONTEND))
    return _COMPONENT(html=view_html, js=CLICK_JS, key=key, default=None)