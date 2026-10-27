"""
components/variants_table.py

The "Variants" view of the Co-occurrence results page: the panel tree and
the per-city evidence in ONE table (2026-09-30).

  Variant   the tree, coloured by the evidence in the chosen city; the
            recombinant families (XEC, XFG, XDV → NB.1.8.1 …) under their
            own root, not under B; "+ Add" next to a confirmed finding
  cities    one cell per city — click a city code, or any cell of its column,
            to colour the tree by it
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
from urllib.parse import quote_plus
from pathlib import Path

# one palette for cells, tree dots and legend
GREEN = "#16a34a"          # panel variant present
AMBER = "#f59e0b"          # panel variant not found (fill / dot)
AMBER_T = "#b45309"        # … as text
VIOLET = "#7c3aed"         # no longer drawn (2026-10-14: one colour = one meaning)
GREY = "#9ca3af"           # no ★ marker / no data
RED = "#dc2626"            # found, not in panel
BLUE = "#2563eb"           # novel
SLATE = "#64748b"          # too broad to name
TRACKED = "#1f2430"        # officially tracked, not selected: hollow ring in the
                           #   panel ink (filled = in your panel); blue means novel

CSS = """
<style>
/* narrow screens (2026-10-22): the table scrolls sideways inside its box
   instead of being cut at the right, and the lineage column stays put */
.vtw { overflow-x:auto; max-width:100%; -webkit-overflow-scrolling:touch; }
.vt { border-collapse:collapse; font-size:15px; width:100%; }
.vt td.tr, .vt th:first-child { position:sticky; left:0; z-index:2; background:#fff; }
.vt tr.dt td, .vt tr.sec td, .vt tr.grp td, .vt tr.note td { position:static; }
@media (max-width: 1100px) {
  .vt { font-size:14px; }
  .vt td.ov { min-width:160px; max-width:220px; }
  .vt th.cc, .vt td.cc { width:54px; }
  .vt td.cc .c { width:50px; font-size:11.5px; }
  .vt .spl { max-width:160px; }
}
.vt th { text-align:left; font-size:12px; letter-spacing:.03em; text-transform:uppercase; color:#6b7280;
         font-weight:600; padding:0 8px 6px; border-bottom:1px solid #e6e8ec; white-space:nowrap; }
.vt th.cc { width:66px; padding:0 2px 6px; text-align:center; letter-spacing:0; }
.vt th.cc button { font:inherit; font-size:12px; font-weight:700; color:#6b7280; background:none; border:0;
                   border-radius:5px; padding:2px 5px; cursor:pointer; }
.vt th.cc button:hover { background:#eef0f3; color:#31333f; }
.vt th.cc.sel button { background:#31333f; color:#fff; }
.vt td { padding:0 8px; height:34px; border-bottom:1px solid #eef0f2; white-space:nowrap; vertical-align:middle; }
.vt td.cc { padding:0 2px; width:66px; text-align:center; cursor:pointer; }
.vt td.nomk { color:#8b8f97; font-size:13px; padding-left:8px; }
.vt td.ov { white-space:normal; padding:2px 8px; min-width:250px; max-width:300px; }
.vt .ovb { display:inline-block; font-size:12.5px; line-height:1.3; padding:2px 7px; border-radius:5px; }
.vt .ovw { display:flex; justify-content:space-between; align-items:center; gap:10px; }
.vt table.wkt { border-collapse:separate; border-spacing:2px 2px; font-size:12px; margin:4px 0 2px; width:auto; }
.vt table.wkt th { font-weight:400; color:#6b7280; padding:1px 5px; text-align:center;
                   text-transform:none; letter-spacing:0; border:0 !important; background:none; }
.vt table.wkt td, .vt tr.dt table.wkt td { padding:2px 6px; text-align:center; border-radius:3px; border:0 !important;
                   height:auto; line-height:1.35; vertical-align:middle; }
.vt table.wkt td.rl, .vt tr.dt table.wkt td.rl { text-align:left; padding-right:10px; white-space:nowrap; }
.vt table.wkt td.rr, .vt tr.dt table.wkt td.rr { text-align:left; color:#6b7280; padding-left:10px; white-space:nowrap; }
.vt table.wkt tr.sep0:first-of-type td, .vt table.wkt tr.sep0 td { }
.vt table.wkt i.sq { display:inline-block; width:10px; height:10px; border-radius:2px; margin:0 1px; vertical-align:-1px; }
.vt .dim i.sq { display:inline-block; width:10px; height:10px; border-radius:2px; margin:0 2px 0 6px; vertical-align:-1px; }
.vt .dh { font-size:12.5px; font-weight:600; margin:10px 0 2px; }
.vt .ovw .ovb { white-space:nowrap; }
.vt .wks { display:inline-flex; gap:2px; align-items:center; white-space:nowrap; flex:none; }
.vt .wks i { display:inline-block; width:12px; height:14px; border-radius:2px; }
.vt .wks .wkm { font-size:12px; color:#9ca3af; margin-right:2px; }
.vt .ovb b { font-weight:600; }
.vt .ovb .add { margin-left:6px; }
.vt td.cc .c { cursor:pointer; }
.vt td.cc:not(.sel):hover { background:rgba(49,51,63,.04); }
.vt td.cc.sel { background:rgba(49,51,63,.07); }
.vt tr.sp td { height:24px; border-bottom:none; }
.vt td.tr { padding-right:12px; }   /* sticky (top of CSS) also anchors the tree lines */
.vt .v { position:absolute; top:0; bottom:0; border-left:1.5px solid #d3d1c7; }
.vt .v.half { bottom:50%; }
.vt .v.down { top:50%; }
.vt .h { position:absolute; top:50%; height:0; border-top:1.5px solid #d3d1c7; }
.vt .dot { position:absolute; top:50%; border-radius:50%; box-sizing:border-box; }
.vt .nm { font-weight:600; }
.vt .spl { color:#9a988f; font-size:12.5px; display:inline-block; max-width:280px; overflow:hidden;
           text-overflow:ellipsis; vertical-align:middle; }
.vt .par { color:#6b7280; font-size:12.5px; margin-left:8px; }
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
.vt .c { display:inline-flex; align-items:center; justify-content:center; min-width:28px; padding:0 3px; height:22px;
         border-radius:5px; font-size:12.5px; font-weight:700; box-sizing:border-box; color:#c4c7cd; cursor:default;
         white-space:nowrap; font-variant-numeric:tabular-nums; }
.vt td.cc .c { width:62px; padding:0 2px; }   /* one width for every city cell (2026-10-20) */
.vt .c .cu, .lg .c .cu { font-weight:400; font-size:10.5px; margin-left:2px; opacity:.8; }
.vt td.ev { white-space:normal; line-height:1.3; padding-left:12px; font-size:14.5px; }
.vt .dim { color:#8b8f97; }
.vt .tr-a { font-size:11px; margin-left:2px; }
.vt .det { color:#2563eb; font-size:12.5px; cursor:pointer; margin-left:6px; white-space:nowrap; }
.vt tr.dt td { height:auto; padding:8px 12px 10px; background:#f8fafc; white-space:normal;
               font-size:13px; line-height:1.55; border-left:3px solid #e2e8f0; }
.vt .dl { display:flex; align-items:center; gap:10px; margin:2px 0; }
.vt .dl b, .vt .dl .dn { min-width:84px; }
.vt .dl b { font-weight:600; }
.vt .gap { color:#9ca3af; font-size:11px; margin:0 3px; }
.vt table.mk { border-collapse:collapse; margin-top:6px; }
.vt table.mk td, .vt table.mk th { height:auto; padding:1px 14px 1px 0; border:0; font-size:12.5px;
                                   text-transform:none; letter-spacing:0; background:none; }
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
.lg.lgs { margin:0 0 6px; padding:4px 8px; border:1px solid #e5e7eb; border-radius:8px; background:#f8fafc; font-size:12.5px; line-height:1.7; }
.lg.lgs .row { gap:2px 14px; }
.lg .cal { display:inline-flex; }
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
S_NOMARK = _fill("#e5e7eb", "#4b5563")
S_NODATA = _dashed(GREY)
S_FOUND = _fill(RED, "#fff")
S_ONEDAY = _dashed(RED)
S_NAMED = _fill("#f3f4f6", "#6b7280")   # not confirmed: grey, neither yes nor no
S_NOVEL = _fill(BLUE, "#fff")
S_NOVEL1 = _dashed(BLUE)
S_BROAD = _fill("#e2e8f0", "#334155")
LIGHT_GREEN = "#4ade80"    # panel variant ± 1 change (as in the completeness graph)
S_NEAR = _fill("#bbf7d0", "#14532d")
S_NEAR1 = _dashed("#16a34a")


def _lab_html(label) -> str:
    """'0:24 (15)' with the bracket smaller and lighter: the score stands out.
    '3\non 2 samples' (specific co-occurrences, 2026-10-26): the second line
    small, under the number."""
    if "\n" in str(label):
        top, sub = str(label).split("\n", 1)
        return (f"{_e(top)}<span class='cu' style='display:block;font-size:10px;"
                f"line-height:1.1'>{_e(sub)}</span>")
    return re.sub(r" \((\d+)\)$", r"<span class='cu'>(\1)</span>", _e(label))


def _cell(label, style, tip, sel) -> str:
    return (f"<td class='cc{' sel' if sel else ''}'><span class='c' style='{style}' "
            f"data-tip='{_e(tip).replace(chr(10), '&#10;')}'>{_lab_html(label)}</span></td>")


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
    if word.endswith(("ch", "sh", "s", "x")):
        return word + "es"                 # stretch → stretches
    return word + "s"


# ── per-row pieces ─────────────────────────────────────────────────────────

# ── cells: one fixed colour per answer of the ★ check (2026-10-09) ─────────
# The check answers present / mixed / absent (process.cooc.check_verdicts);
# the cell shows that answer, the count says how many markers voted.
#                 in your panel                  not in your panel
#   present       green: fine                    red: missing from your panel
#   mixed         grey                           grey
#   absent        orange: in the panel, not in   pale: nothing to do
#                 the wastewater
# (Until 2026-10-09 a gradient by the share present, on two different scales
# for panel and found rows, so the same 2/5 had two colours.)
# Present in the sum, but clearly present in fewer than evidence.min_days
# samples on their own (2026-10-09): a light fill with a thin solid edge of
# the same hue (XFZ in Lugano: 2/2 from one sample, 24 Aug).
# (2026-10-20) The cell shows the check's answer only, in three colours, and
# its label counts every marker: "2✓ 1~ 1✗" (present · no vote · absent), so
# no marker drops out of view (the old "1/1" for XFG in Lugano hid 4214C, in
# between). Whether the lineage is in your panel is said once in the row
# ("in your panel · absent in 6 of 6 cities"), not by the cell's colour; the
# light shades are gone (their two meanings go into the hover in words).
#   S_YES   present                     green
#   S_NO    absent                      grey
#   S_OPEN  tie, or not called (fewer   white, dashed grey edge (orange is
#           than min_present present),  kept for one meaning only: a lineage
#           or nothing decided          of your panel no city calls present)
S_YES = _fill(GREEN, "#fff")
S_NO = _fill("#f3f4f6", "#4b5563") + "box-shadow:inset 0 0 0 1px #d1d5db;"
S_OPEN = "background:#fff;color:#374151;border:1.5px dashed #9ca3af;"
# seen by its single mutations, but no sample with specific co-occurrences
# (2026-10-26): light, outlined green, weaker than the filled green
S_SINGLE = "background:#f0fdf4;color:#15803d;box-shadow:inset 0 0 0 1.5px #86efac;"
# kept for the legend of older pieces and the deep scan's own cells
S_PRESENT_LIGHT = "background:#fff;color:#14532d;box-shadow:inset 0 0 0 1.5px #16a34a;"
S_NOTFOUND_LIGHT = "background:#fff;color:#78350f;box-shadow:inset 0 0 0 1.5px #f59e0b;"
S_FOUND_LIGHT = "background:#fff;color:#7f1d1d;box-shadow:inset 0 0 0 1.5px #dc2626;"
MIXED_BG, PALE_BG, PALE_FG = "#b4bac3", "#f3f4f6", "#6b7280"
S_MIXED = _fill(MIXED_BG, "#111827")
S_PALE = _fill(PALE_BG, PALE_FG) + "box-shadow:inset 0 0 0 1px #e5e7eb;"
S_FEW = _fill(PALE_BG, "#374151") + "box-shadow:inset 0 0 0 1px #9ca3af;"   # too few markers
D_PALE = f"background:{PALE_BG};box-shadow:inset 0 0 0 1px #d1d5db;"   # a pale day square


def _state(d):
    """The check's answer for one cell: the result's own state, or (for
    results without one, e.g. the legend) the vote from its counts.
    "few" (2026-10-16): "mixed" with no marker absent — every measured marker
    present, but fewer than check.min_present measured: too little data, not
    markers that disagree."""
    st = d.get("state")
    n_mk, p, a, u = _vote_counts(d)
    if st not in ("present", "absent", "mixed"):
        if p + a == 0:
            return "not_covered"
        # the same majority vote as process.cooc.check_verdicts (2026-10-20)
        if p > a and p >= _min_present():
            st = "present"
        else:
            st = "absent" if a > p else "mixed"
    if st == "mixed" and a == 0 and 0 < p < _min_present():
        return "few"
    return st


def _clear_samples(d):
    """Samples where the lineage is clearly present on its own (the
    calendar's present squares); None for results without a calendar."""
    tl = d.get("timeline")
    if tl is None:
        return None
    return [x for x, mk in tl if mk == "present"]


def _few_clear(d):
    """Present in the sum, but clear in fewer than evidence.min_days samples."""
    c = _clear_samples(d)
    return c is not None and len(c) < _num()["ed"]


def _well_note(d):
    """' — 1 of 12 samples well covered: read with care' for cell hovers."""
    n, k = d.get("n_samples"), d.get("n_well")
    if not n or k is None:
        return ""
    return (f" — {k} of {n} {_pl(n, 'sample')} well covered"
            + (": read with care" if k * 2 < n else ""))


def _clear_note(d):
    """Hover words for a light cell: on which sample(s) it was clear."""
    if _state(d) != "present" or not _few_clear(d):
        return ""
    c = _clear_samples(d) or []
    when = (f"only {len(c)} {_pl(len(c), 'sample')} on its own ("
            + ", ".join(_short_date(x) for x in c) + ")") if c else "no sample on its own"
    return (f" — clearly present in {when}; the sum of the samples says present, the "
            "single samples were too thin, mixed or absent")


def _counts_label(d):
    """'2:1 (1)': markers present : absent, and in brackets those that
    couldn't answer (too few reads, between 1 and 5 %, or reads unlike the
    lineage) — like a score, nothing left out (2026-10-20)."""
    n_mk, p, a, u = _vote_counts(d)
    if n_mk == 0:
        return "–"
    return f"{p}:{a}" + (f" ({u})" if u else "")


def _thin(d):
    """A city whose samples in the chosen dates are mostly thin: fewer than
    check.thin_well_min well covered (None when the result has no coverage
    data)."""
    nw = d.get("n_well")
    return nw is not None and nw < int(_check_cfg_safe().get("thin_well_min", 2))


def _scale(d, found=False):
    """(cell style, label, dot colour, name colour, share present or None)
    for one city's check result. found: a lineage not in the panel.
    The cell (style, label) shows the answer only (2026-10-20); the dot and
    name colours of the tree are unchanged."""
    n_mk, p, a, u = _vote_counts(d)
    pev = d.get("pairs") or {}
    if pev.get("positive"):
        # specific co-occurrences (2026-10-26): how many stretches of genome
        # carried the lineage's own combinations, on how many samples
        # one short label for every cell with pair counts (2026-10-26): samples
        # carrying its own combinations on 2+ stretches / samples where pairs
        # were read. The stretches are in the hover.
        k, n = len(pev["positive"]), len(pev.get("measured") or pev["positive"])
        lab = f"{k}/{n}"
        thin = not pev.get("positive_well")
        sty = S_YES + _LOWQ_CELL if thin else S_YES
        col = RED if found else GREEN
        return (sty, lab, (("#fca5a5" if found else "#86efac") if thin else col),
                (RED if found else "#15803d"), None)
    if pev and pev.get("measured"):
        # pairs counted, no sample with 2 stretches (2026-10-26): one measure
        # in every cell. Seen by single mutations only -> light green with the
        # old count; otherwise 0 and on how many samples the pairs were read
        k = len(pev["measured"])
        lab = f"0/{k}"
        if d.get("state_markers") == "present":
            return (S_SINGLE, lab, ("#fca5a5" if found else "#86efac"),
                    (RED if found else "#15803d"), None)
        if d.get("state_markers") != "absent":
            # single mutations undecided (a tie, too few reads): no answer,
            # not "not seen" (XFG in Basel, Sep-Oct: 1:0 (1), dominant there)
            return S_OPEN, lab, GREY, (RED if found else "#6b7280"), None
        if not found:
            return S_NO, lab, AMBER, AMBER_T, 0.0
        return S_NO, lab, PALE_BG, RED, 0.0
    if n_mk == 0:
        return S_NOMARK, "·", GREY, "#6b7280", None
    lab = _counts_label(d)
    m = p + a
    if m == 0:
        return S_OPEN, lab, GREY, "#6b7280", None
    share, st = p / m, _state(d)
    weak = d.get("basis") == "all"          # no well-covered sample decided it
    if st == "present":
        col = RED if found else GREEN
        # seen in a city with fewer than check.thin_well_min well-covered
        # samples: striped, as a low-coverage square (2026-10-23, display only)
        sty = S_YES + _LOWQ_CELL if _thin(d) else S_YES
        if weak or _few_clear(d):
            return (sty, lab, ("#fca5a5" if found else "#86efac"),
                    (RED if found else "#15803d"), share)
        return sty, lab, col, (RED if found else "#15803d"), share
    if st == "absent" and not found:
        return S_NO, lab, ("#fcd9a3" if weak else AMBER), AMBER_T, share
    if st == "absent":
        return S_NO, lab, PALE_BG, RED, share
    if st == "few":
        return S_OPEN, lab, PALE_BG, (RED if found else "#6b7280"), share
    return S_OPEN, lab, MIXED_BG, (RED if found else "#4b5563"), share


def _basis_samples(d):
    """(samples the vote used, as [(date, mark)], 'weighted' / 'well' /
    'all'): every sample, weighted by coverage (2026-10-16); the well-covered
    ones for results saved before; all when there's no coverage data."""
    tl = [(x, m) for x, m in (d.get("timeline") or []) if m != "unknown"]
    if d.get("basis") == "weighted":
        return tl, "weighted"
    if d.get("basis") != "well":
        return tl, "all"
    try:
        from process.sample_qc import well_covered
    except Exception:
        return tl, "all"
    q = d.get("quality") or {}
    return [(x, m) for x, m in tl if well_covered(q.get(x))], "well"


def _summary(d):
    """One line for the selected city: the answer, and what it rests on (the
    trend arrows were removed 2026-10-16: the details calendar shows the
    samples)."""
    n_mk, p, a, u = _vote_counts(d)
    if n_mk == 0:
        return _vote_text(d)
    st = _state(d)
    word = _WORD.get(st, "too few reads")
    samples, basis = _basis_samples(d)
    clear = [x for x, m in samples if m == "present"]
    k = len(samples)
    if basis == "weighted":
        nw = d.get("n_well")
        txt = (f"{word} on all {k} {_pl(k, 'sample')}, weighted by coverage"
               + (f" ({nw} well covered)" if nw is not None else ""))
    elif basis == "well":
        txt = f"{word} on the {k} well-covered {_pl(k, 'sample')}"
    else:
        txt = f"{word} · weak: no well-covered sample, all {k} samples counted"
    if clear:
        txt += (f" · last {_short_date(clear[-1])}" if basis == "well"
                else f" · clearly present in {len(clear)}, last {_short_date(clear[-1])}")
    return txt


# hover symbols (2026-10-13): short tips instead of sentences; the legend
# explains them once
SYM = {"present": "✓", "absent": "✗", "between": "~", "low": "?", "nofit": "≠"}


def _marker_syms(d, n_max=6):
    """'4214C 50 % ✓ · 8350C 0 % ✗ · 1148T 40 reads ?' from a check's detail."""
    out = []
    for x in d.get("detail") or []:
        reads = x.get("cov_raw", x.get("cov", 0)) or 0
        f, st, why = x.get("freq"), x.get("status"), x.get("why") or ""
        m = x["marker"]
        if not reads:
            out.append(f"{m} no reads ?")
        elif why.startswith("under"):                 # too few reads at the marker
            cnt = x.get("cov", reads) or 0
            out.append(f"{m} {reads:,} reads ?" if cnt >= reads else
                       f"{m} {reads:,} reads, count as {cnt:,} ?")   # low-coverage samples
        elif "cover a neighbour" in why:              # share fine, too few reads to test the link
            out.append(f"{m} {100 * (f or 0):.0f} % ({x.get('link_n', 0)} linked reads) ?")
        elif "look like" in why:                      # its reads look like another lineage
            out.append(f"{m} {100 * (f or 0):.0f} % ≠")
        else:
            sym = ("✓" if st == "present" or why.startswith("present")
                   else "✗" if why.startswith("absent") else "~")
            out.append(f"{m} {100 * (f or 0):.0f} % {sym}"
                       + (" (share only)" if "share only" in why else ""))
    more = f" · +{len(out) - n_max}" if len(out) > n_max else ""
    return " · ".join(out[:n_max]) + more


def _pairs_text(pev, limit=3):
    """'3 Oct (well covered): 22879A+22899C 77 % of 10,098, 23854A+24077A 60 %
    of 14,802 · 25 Sep (thin): …' — the specific co-occurrences behind a cell."""
    q = pev.get("positive_well") or []
    out = []
    for dt in sorted(pev.get("positive") or [], reverse=True)[:limit]:
        prs = ", ".join(f"{x[0]} {100 * x[1]:.0f} % of {x[2]:,}"
                        + ("" if len(x) < 4 or x[3] else " (backbone mutation)")
                        for x in pev["per_date"].get(dt, []))
        out.append(f"{_short_date(dt)} ({'well covered' if dt in q else 'thin'}): {prs}")
    return " · ".join(out)


def _cell_tip(v, city, d, extra=""):
    """'XFG · Basel — mixed (1 ✓ 1 ✗ of 2 ★) · 2 well-covered samples · 4214C
    50 % ✓ · 8350C 0 % ✗' (+ why it was checked)."""
    pev = d.get("pairs") or {}
    if pev.get("positive") or pev.get("measured"):
        # one line per point (2026-10-26; the tip keeps line breaks)
        nk, sp, sa, su = _vote_counts(d)
        seen_m = [x["marker"] for x in (d.get("detail") or []) if x.get("status") == "present"]
        single = (f"Single mutations alone: {sp} seen"
                  + (f" ({', '.join(seen_m[:4])}{' …' if len(seen_m) > 4 else ''})" if seen_m else "")
                  + f", {sa} not seen" + (f", {su} couldn't answer" if su else "")
                  + f" — of {nk}") if nk else "No single mutation of its own"
        n_read = len(pev.get("measured") or [])
        if pev.get("positive"):
            k, n = len(pev["positive"]), pev["n_stretches"]
            lines = [f"{v} · {city} — seen on {k} of {n_read} {_pl(n_read, 'sample')}",
                     "Its own combinations on the same reads (pair · share · reads):"]
            q = pev.get("positive_well") or []
            for dt in sorted(pev["positive"], reverse=True)[:4]:
                prs = " · ".join(f"{x[0]} {100 * x[1]:.0f} % of {x[2]:,}"
                                 + ("" if len(x) < 4 or x[3] else " (backbone mutation)")
                                 for x in pev["per_date"].get(dt, []))
                lines.append(f"  {_short_date(dt)} ({'well covered' if dt in q else 'thin'}): {prs}")
            lines.append(f"{n} separate {'cooc' if n == 1 else 'coocs'} in all; a sample counts with "
                         "2 or more, at least one without a backbone mutation")
            if not q:
                lines.append("Striped: only on thin samples, weak on its own")
            lines.append(single)
        elif d.get("state_markers") == "present":
            lines = [f"{v} · {city} — weaker: seen by single mutations only", single,
                     f"Its own combinations on the same reads: none on 2 or more stretches of one "
                     f"sample ({n_read} {_pl(n_read, 'sample')} read)",
                     "A single mutation can come from another lineage's reads: check the next "
                     "samples"]
        elif d.get("state_markers") != "absent":
            lines = [f"{v} · {city} — undecided",
                     f"Its own combinations on the same reads: none on 2 or more stretches of one "
                     f"sample ({n_read} {_pl(n_read, 'sample')} read)", single,
                     "Too little to say seen or not seen"]
        else:
            lines = [f"{v} · {city} — not seen",
                     f"Its own combinations on the same reads: none on 2 or more stretches of one "
                     f"sample ({n_read} {_pl(n_read, 'sample')} read)", single]
        return "\n".join(lines) + extra
    n_mk, p, a, u = _vote_counts(d)
    if n_mk == 0:
        return f"{v} · {city} — no mutation of its own"
    # plain words (2026-10-21): the old ✓ ✗ ? ★ symbols are gone from the page
    word = {"present": "seen", "absent": "not seen"}.get(_state(d), "undecided")
    counts = (f"{p} seen, {a} not seen" + (f", {u} couldn't answer" if u else "")
              + f", of {n_mk} own {_pl(n_mk, 'mutation')}")
    samples, basis = _basis_samples(d)
    base = (f"{len(samples)} {_pl(len(samples), 'sample')} weighted by coverage (each ≤ "
            f"{_cov_factor():g} × its median depth, ≤ {_cap_reads():,} reads)" if basis == "weighted"
            else f"{len(samples)} well-covered {_pl(len(samples), 'sample')}" if basis == "well"
            else "weak: no coverage data (all counted, ≤ "
                 f"{_cap_reads():,} reads each)")
    if _state(d) == "few":
        base += f" · only {p} seen, {_min_present()} needed"
    ms = " · ".join(f"{x['marker']}: {_answer_text(x)}" for x in (d.get("detail") or [])[:6])
    if len(d.get("detail") or []) > 6:
        ms += f" · +{len(d['detail']) - 6} more"
    # what the light shades said until 2026-10-20 is in words: "weak: no
    # coverage data" in base, and on how few samples it was clear
    if _state(d) == "present" and _thin(d):
        base += (f" · striped: {d.get('n_well')} of {d.get('n_samples')} samples well covered, "
                 "so this city alone is weak evidence")
    if d.get("own_missing"):
        own, n_o = ", ".join(d.get("own") or []), len(d.get("own") or [])
        base += (f" · {'the mutation' if n_o == 1 else 'the mutations'} new in {v} ({own}) "
                 f"not seen: the seen ones {d.get('parent') or 'its parent'} carries too")
    return (f"{v} · {city} — {word} ({counts}) · {base}" + _clear_note(d)
            + (f" · {ms}" if ms else "") + extra)


def _det_link(gid):
    return f"<span class='det' data-grp='{gid}'>details <span class='arr'>▸</span></span>"


_SHORT = {"present": "✓ present", "absent": "✗ absent", "mixed": "◐ mixed",
          "thin": "? too few reads", "uncovered": "? no reads at its markers"}


def _square(x, m, mk, detail, quality):
    """'27 Aug · ✓ well covered (depth 9,212 · 84 % ≥ 100×) · ✗ absent —
    4214C 0 % of 12,480 reads ✗ · …'"""
    style, _meaning = mk.get(m, mk["unknown"])
    extra, _qtxt = _qual_bits(quality, x)
    q = (quality or {}).get(x)
    head = ""
    if q:
        num = f"(depth {q.get('median', 0):,.0f} · {100 * q.get('share_100', 0):.0f} % ≥ 100×)"
        head = ("▧ low coverage " if extra is None else "✓ well covered ") + num + " · "
    if extra is None:                  # low coverage: stripes over its colour
        extra = _LOWQ_EMPTY if m in ("thin", "uncovered") else _LOWQ
    tip = (_short_date(x) + " · " + head + _SHORT.get(m, m)
           + (" — " + detail[x] if (detail or {}).get(x) else ""))
    return f"<i style='{style}{extra}' data-tip='{_e(tip)}'></i>"


def _cal_line(d, top, low=None):
    """The newest check.calendar_samples samples, and before them, after a
    small gap, any older sample where it was clearly present (2026-10-13)."""
    tl = [(x, m) for x, m in (d.get("timeline") or []) if m != "unknown"]
    if not tl:
        return ""
    n = _cal_max()
    mk = _day_marks(top, low)
    det, q = d.get("day_detail") or {}, d.get("quality") or {}
    try:
        from process.sample_qc import well_covered
    except Exception:
        well_covered = lambda _q: False                      # noqa: E731
    # older samples that matter: where it was present, and the well-covered
    # ones (the cell may rest on them)
    older = [(x, m) for x, m in tl[:-n] if m == "present" or well_covered(q.get(x))]
    head = "".join(_square(x, m, mk, det, q) for x, m in older)
    if older:
        hidden = len(tl) - n - len(older)
        head += (f"<span class='gap' data-tip='{hidden} older {_pl(hidden, 'sample')} "
                 f"not shown (low coverage, without it)'>⋯</span>")
    return ("<span class='cal' style='margin-left:0'>" + head
            + "".join(_square(x, m, mk, det, q) for x, m in tl[-n:]) + "</span>")


# Outside links (2026-10-13): the lineage's wastewater page on GenSpectrum
# (mutation heatmaps, chosen city) and its ★ markers in clinical sequences on
# CoV-Spectrum (which lineages carry them, alone or together).
GENSPECTRUM_WW = ("https://genspectrum.org/swiss-wastewater/covid?locationName={city}"
                  "&samplingDate=Most+recent+30+days&granularity=day&analysisMode=variant"
                  "&sequenceType=nucleotide&signatureType=computed&variant={lineage}*"
                  "&minProportion=0.8&minCount=15&minJaccard=0.75&timeFrame=6months")
COVSPECTRUM_MUTS = ("https://cov-spectrum.org/explore/World/AllSamples/Past6M/variants"
                    "?nucMutations={muts}")


def _ext(url, text, tip=""):
    t = f" data-tip='{_e(tip)}'" if tip else ""
    return (f"<a href='{_e(url)}' target='_blank' rel='noopener'{t} "
            f"style='color:#2563eb;text-decoration:none'>{_e(text)} ↗</a>")


def _ww_link(v, city):
    return GENSPECTRUM_WW.format(city=quote_plus(city), lineage=quote_plus(v))


def _clin_link(muts):
    return COVSPECTRUM_MUTS.format(muts=quote_plus(",".join(muts)))


# ── the details' words (2026-10-20): numbers in plain words ──────────────
def _check_cfg_safe():
    try:
        from process.cooc import _check_cfg
        return _check_cfg()
    except Exception:
        from process.cooc import CHECK_DEFAULTS
        return dict(CHECK_DEFAULTS)


def _city_line(d):
    """'2:0 · 9 samples, 4 well covered · seen on 1 alone, 28 Aug'"""
    parts = [_counts_label(d)]
    n, k = d.get("n_samples"), d.get("n_well")
    if n:
        parts.append(f"{n} {_pl(n, 'sample')}" + (f", {k} well covered" if k is not None else ""))
    clear = _clear_samples(d) or []
    if clear:
        parts.append(f"seen on {len(clear)} alone, last {_short_date(clear[-1])}")
    if d.get("own_missing"):
        parts.append(f"only mutations {d.get('parent') or 'its parent'} carries too seen")
    return " · ".join(parts)


def _share_text(x):
    """'44.1 % of 65 counted reads (142,546 in total)' — the share is computed
    on the counted reads; until 2026-10-20 it was printed next to the total."""
    cov, raw, f = x.get("cov", 0) or 0, x.get("cov_raw", x.get("cov", 0)) or 0, x.get("freq")
    if not raw or f is None:
        return "no reads"
    if cov < raw:
        return f"{100 * f:.1f} % of {cov:,} counted reads ({raw:,} in total)"
    return f"{100 * f:.1f} % of {raw:,} reads"


def _answer_text(x):
    """seen / not seen / couldn't answer, with the reason in numbers."""
    c = _check_cfg_safe()
    st, why, f = x.get("status"), x.get("why") or "", x.get("freq")
    pf, af = c["present_freq"] * 100, c["absent_freq"] * 100
    cov, ln = x.get("cov", 0) or 0, x.get("link_n", 0) or 0
    if st == "present":
        return (f"seen (≥ {pf:g} %), on its share alone: no neighbour in reach"
                if "share only" in why else f"seen (≥ {pf:g} %)")
    if st == "absent":
        return f"not seen (< {af:g} %)"
    if not cov:
        return "couldn't answer: no reads"
    if cov < c["min_cov"]:
        return f"couldn't answer: {cov:,} counted reads, needs {c['min_cov']}"
    if f is not None and f < c["present_freq"]:
        return f"couldn't answer: between {af:g} and {pf:g} %"
    if ln < c["min_link"]:
        return f"couldn't answer: only {ln} of its reads reach a neighbour, needs {c['min_link']}"
    link = x.get("link")
    return (f"couldn't answer: its reads look like the lineage nearby on "
            f"{100 * (link or 0):.0f} %, needs {c['link_min'] * 100:g} %")


def _together_text(x):
    """'96 % of 1,240 reads' — the co-occurrence at the mutation's neighbours."""
    if x.get("link") is not None and x.get("link_n"):
        return f"{100 * x['link']:.0f} % of {x['link_n']:,} reads"
    if "share only" in (x.get("why") or ""):
        return "no neighbour in reach"
    return "–"


# ── Which mutations "reads together" is about (2026-10-26) ──────────────────
# The check counts, for a marker, the reads that also reach another position of
# the lineage and match it there (process.cooc.accumulate_check_stats). Those
# neighbours are the lineage's own mutations within a read's reach: the same
# combination W-ASAP shows per amplicon (TA.1, Zürich, Sep–Oct: 6196A with
# 6183G, its amplicon 20 at Jaccard .97; 23157G with 23148C and 23401T, its
# amplicon 75 at 1.0). The hover names them, with how many lineages carry each:
# the marker is the specific part, a neighbour carried by thousands of lineages
# only places it on the right background.
NB_SPAN = 300      # bases one read covers (process.variant_explorer.marker_blocks)
NB_RARE = 500      # neighbours carried by fewer lineages are named one by one
_NB_CACHE: dict = {}


def _pango_sig():
    """(signatures, parent, carriers) from the page's pango loader, or None."""
    # the tree file the page uses (api.pango_loader.tree_info), loaded once per
    # file; not imported from the page module, which would run the page
    try:
        from api.pango_loader import PangoLoader, tree_info
        from process.variant_explorer import _indices
        t = tree_info()
        if _NB_CACHE.get("tree") != t["id"]:
            _NB_CACHE.clear()
            _NB_CACHE["tree"], _NB_CACHE["loader"] = t["id"], PangoLoader(t["file"])
        _raw, sigs, parent, _t = _indices(_NB_CACHE["loader"])
    except Exception:
        return None
    if _NB_CACHE.get("sigs") is not sigs:
        from collections import Counter
        _NB_CACHE["sigs"] = sigs
        _NB_CACHE["car"] = Counter(m for v in sigs.values() for m in v)
    return sigs, parent, _NB_CACHE["car"]


def _neighbours(v, m):
    """[(mutation, inherited from the parent, lineages carrying it)] of `v`
    within NB_SPAN bases of marker `m`, nearest first; None without pango data."""
    ps = _pango_sig()
    if not ps:
        return None
    sigs, parent, car = ps
    pos = lambda t: int("".join(ch for ch in t if ch.isdigit()) or 0)   # noqa: E731
    p, psig = pos(m), sigs.get(parent.get(v, ""), set())
    nb = [t for t in sigs.get(v, set()) if t != m and 0 < abs(pos(t) - p) <= NB_SPAN]
    return [(t, t in psig, car.get(t, 0)) for t in sorted(nb, key=lambda t: abs(pos(t) - p))]


def _together_tip(v, x):
    """The hover of 'reads together': which of the lineage's mutations the
    marker's reads were checked against."""
    m = x.get("marker", "")
    nb = _neighbours(v, m)
    if nb is None:
        return ""
    if not nb:
        return (f"no other {v} mutation within a read's reach ({NB_SPAN} bases) of {m}: "
                "judged on its share alone")
    rare = sorted([t for t in nb if t[2] < NB_RARE], key=lambda t: t[2])
    common = [t for t in nb if t[2] >= NB_RARE]
    named = ", ".join(f"{t} ({'' if inh else 'new in ' + v + ', '}in {n:,} "
                      f"{_pl(n, 'lineage')})" for t, inh, n in rare)
    lo, hi = (min(t[2] for t in common), max(t[2] for t in common)) if common else (0, 0)
    more = (f"{len(common)} common {_pl(len(common), 'one')} (in "
            f"{lo:,}{'' if lo == hi else f'–{hi:,}'} lineages: "
            f"{', '.join(t[0] for t in common[:4])}{' …' if len(common) > 4 else ''})"
            if common else "")
    head = (f"{100 * x['link']:.0f} % of {x['link_n']:,} reads carrying {m} also carry {v}'s "
            "mutations they reach. " if x.get("link") is not None and x.get("link_n") else "")
    return (head + f"{v}'s mutations within a read's reach ({NB_SPAN} bases): "
            + "; ".join(t for t in (named, more) if t)
            + f". {m} is the specific part; a common one only shows it sits on the "
            "lineage's background, not on another branch's")


# ── All cities, one square per week (2026-10-22): on top of the city lines in
#    a row's details. Each square is the cells' vote on every sample of every
#    city that week (process.national); striped = fewer than 2 well-covered
#    samples that week. ──────────────────────────────────────────────────────
_NAT_WORD = {"present": "seen", "absent": "not seen", "mixed": "tie",
             "not_covered": "too few reads"}


def _nat_line(per_city, top, low=None):
    """'All 6 cities  ■■■□□  6:0 (1) · seen' — or "" without data."""
    try:
        from process.national import national
        nat = national(per_city)
    except Exception:
        return ""
    if not nat or not nat["weeks"]:
        return ""
    mk = _day_marks(top, low)
    sq = []
    for w in nat["weeks"]:
        a = w["answer"]
        m = {"present": "present", "absent": "absent", "mixed": "mixed"}.get(a["verdict"], "thin")
        style, _ = mk[m]
        extra = "" if w["n_well"] >= 2 else (_LOWQ_EMPTY if m == "thin" else _LOWQ)
        cities = " · ".join(f"{city_name(c)} {_NAT_WORD.get(v, v)}"
                            for c, v in w["cities"].items())
        tip = (f"week of {_short_date(w['start'])} · {_counts_label(a)} · "
               f"{_NAT_WORD.get(a['verdict'], '')} · {w['n_samples']} "
               f"{_pl(w['n_samples'], 'sample')}, {w['n_well']} well covered — {cities}")
        sq.append(f"<i style='{style}{extra}' data-tip='{_e(tip)}'></i>")
    a, n_c = nat["answer"], nat["n_cities"]
    tip = ("all cities pooled, every sample in the chosen dates: its own mutations seen : "
           "not seen (couldn't answer), the same rule as the cells")
    return (f"<div class='dl'><b class='dn'>All {n_c} {_pl(n_c, 'city')}</b>"
            f"<span class='cal' style='margin-left:0'>{''.join(sq)}</span>"
            f"<span data-tip='{_e(tip)}'><b>{_e(_counts_label(a))}</b> · "
            f"{_e(_NAT_WORD.get(a['verdict'], ''))}</span></div>")


_CLIN_WORLD = "https://cov-spectrum.org/explore/World/AllSamples/Past6M/variants?nucMutations={muts}"
_CLIN_EUROPE = "https://cov-spectrum.org/explore/Europe/AllSamples/Past6M/variants?nucMutations={muts}"


def _clin_links(muts, what="its own mutations together"):
    """'… World ↗ · Europe ↗' on CoV-Spectrum (2026-10-27)"""
    q = quote_plus(",".join(muts))
    return (f"<span class='dim'>{_e(what)}, clinical:</span> "
            + _ext(_CLIN_WORLD.format(muts=q), "World", "clinical sequences carrying all of them, worldwide")
            + " · " + _ext(_CLIN_EUROPE.format(muts=q), "Europe", "the same, Europe only"))


def _sample_kind(d, day):
    """one sample's square, in the cells' words (2026-10-27)"""
    pev = (d or {}).get("pairs") or {}
    pos = {str(x)[:10] for x in pev.get("positive") or []}
    posw = {str(x)[:10] for x in pev.get("positive_well") or []}
    mark = dict((str(a)[:10], b) for a, b in ((d or {}).get("timeline") or [])).get(day)
    if day in pos:
        return "coocs" if day in posw else "coocs, thin sample"
    if mark == "present":
        return "single mutations only"
    if mark == "absent":
        return "not seen"
    return "undecided" if mark in ("mixed", "thin") else "too few reads"


def _kind_style():
    """sample square styles (a function: _LOWQ is defined further down)"""
    return {"coocs": f"background:{GREEN}",
               "coocs, thin sample": f"background:{GREEN};" + _LOWQ,
               "single mutations only": "background:#c0dd97",
               "not seen": D_PALE,
               "undecided": "background:#fff;box-shadow:inset 0 0 0 1px #9ca3af;border-radius:2px;"
                            "outline:1px dashed #9ca3af;outline-offset:-1px",
               "too few reads": "background:#fff;box-shadow:inset 0 0 0 1px #e5e7eb"}


def _week_block(v, per_city, cities, clin=None, parent=None):
    """One column per week, oldest left (2026-10-27): wastewater (samples with
    coocs / samples, all cities), Europe and World clinical (sequences of the
    lineage and its descendants / all), then each city's samples."""
    from datetime import date as _d, timedelta as _td
    days = sorted({str(a)[:10] for d in per_city.values() for a, _m in ((d or {}).get("timeline") or [])})
    if not days:
        return ""
    def mon(x):
        y = _d.fromisoformat(x)
        return (y - _td(days=y.weekday())).isoformat()
    weeks, w = [], _d.fromisoformat(mon(days[0]))
    while w.isoformat() <= mon(days[-1]):
        weeks.append(w.isoformat())
        w += _td(days=7)
    head = "".join(f"<th>{_e(_short_date(x))}</th>" for x in weeks)
    rows = []
    # wastewater, all cities
    cells = []
    for wk in weeks:
        n = k = 0
        for c in cities:
            for a, _m in ((per_city.get(c) or {}).get("timeline") or []):
                if mon(str(a)[:10]) == wk:
                    n += 1
                    k += _sample_kind(per_city[c], str(a)[:10]).startswith("coocs")
        shade = ("#97c459" if k >= 3 else "#c0dd97" if k == 2 else "#eaf3de" if k == 1 else "")
        cells.append(f"<td style='{'background:' + shade if shade else 'color:#9ca3af'}' "
                     f"data-tip='week of {_short_date(wk)}: {k} of {n} samples with its coocs, all cities'>"
                     f"{k}/{n}</td>" if n else "<td></td>")
    rows.append(f"<tr><td class='rl'>wastewater, all cities</td>{''.join(cells)}"
                "<td class='rr'>samples with coocs</td></tr>")
    # clinical
    if clin:
        try:
            from process.clinical import family, lineage_weeks
            fam = family(v, parent or {})
        except Exception:
            fam = {v}
        for reg in ("Europe", "World"):
            lw = lineage_weeks((clin or {}).get(reg), fam, weeks)
            if lw is None:
                continue
            cells = []
            for i, wk in enumerate(weeks):
                k, n = lw[wk]
                late = i >= len(weeks) - 2
                thin = n < 20
                sty = (("background:#c0dd97" if k and n and k / n >= .1 else "background:#eaf3de" if k else "")
                       + (";" + _LOWQ_EMPTY if (thin or late) else ""))
                tip = (f"week of {_short_date(wk)}: {k} of {n} clinical sequences ({reg}) are {v} or "
                       "below it" + (" · fewer than 20 sequences" if thin else "")
                       + (" · still filling in" if late else ""))
                cells.append(f"<td style='{sty}' data-tip='{_e(tip)}'>{k}/{n}</td>")
            rows.append(f"<tr><td class='rl'>{reg} clinical</td>{''.join(cells)}"
                        f"<td class='rr'>{_e(v)} and below / all sequences</td></tr>")
    # cities
    for c in cities:
        d = per_city.get(c)
        if not d or not d.get("timeline"):
            continue
        cells = []
        for wk in weeks:
            sq = []
            for a, _m in d.get("timeline") or []:
                day = str(a)[:10]
                if mon(day) != wk:
                    continue
                kind = _sample_kind(d, day)
                pairs = [x[0] for x in (((d.get("pairs") or {}).get("per_date") or {}).get(day) or [])]
                tip = f"{city_name(c)}, {_short_date(day)}: {kind}" + (f" · {', '.join(pairs[:4])}" if pairs else "")
                sq.append(f"<i class='sq' style='{_kind_style().get(kind, '')}' data-tip='{_e(tip)}'></i>")
            cells.append(f"<td>{''.join(sq)}</td>")
        pev = d.get("pairs") or {}
        lab = (f"{len(pev.get('positive') or [])}/{len(pev.get('measured') or d.get('timeline'))}"
               f" · {d.get('n_well') or 0} well covered")
        rows.append(f"<tr class='sep0'><td class='rl'>{_e(city_name(c))}</td>{''.join(cells)}"
                    f"<td class='rr'>{_e(lab)}</td></tr>")
    legend = ("<div class='dim' style='font-size:11px;margin-top:3px'>one square = one sample, hover for "
              f"date and coocs · <i class='sq' style='{_kind_style()['coocs']}'></i>coocs "
              f"<i class='sq' style='{_kind_style()['coocs, thin sample']}'></i>on a thin sample "
              f"<i class='sq' style='{_kind_style()['single mutations only']}'></i>single mutations only "
              f"<i class='sq' style='{_kind_style()['not seen']}'></i>not seen "
              f"<i class='sq' style='{_kind_style()['undecided']}'></i>undecided"
              + (" · striped clinical week = under 20 sequences or still filling in" if clin else "")
              + "</div>")
    return (f"<table class='wkt'><tr><th></th>{head}<th></th></tr>{''.join(rows)}</table>" + legend)


def _pair_pos(k):
    ps = []
    for m in str(k).split("+"):
        mm = re.match(r"^(\d+)", m)
        if mm:
            ps.append(int(mm.group(1)))
    return (min(ps), max(ps)) if ps else (0, 0)


def _cooc_table(per_city, cities):
    """Its coocs seen, one row per stretch of genome (2026-10-27): pairs with a
    gap of at most check.stretch_span bases are one stretch (they count once
    for a sample); cell = samples where a pair of the stretch was seen · best
    share; the pairs are in the hover."""
    seen, backbone = {}, set()                 # pair -> city -> {date: share}
    for c in cities:
        pev = (per_city.get(c) or {}).get("pairs") or {}
        # every pair seen on its counted samples; results saved before only
        # kept the best-read pair of each stretch
        pd = pev.get("per_date_all") or pev.get("per_date") or {}
        for day, lst in pd.items():
            for x in lst or []:
                k, share = x[0], float(x[1] or 0)
                if len(x) > 3 and not x[3]:
                    backbone.add(k)
                d = seen.setdefault(k, {}).setdefault(c, {})
                d[day] = max(d.get(day, 0.0), share)
    if not seen:
        return ""
    span = int((_check_cfg_safe() or {}).get("stretch_span", 400))
    groups, last = [], None
    for k in sorted(seen, key=_pair_pos):
        lo, hi = _pair_pos(k)
        if groups and lo - last <= span:
            groups[-1].append(k)
            last = max(last, hi)
        else:
            groups.append([k])
            last = hi
    def stat(g, c):
        days = {}
        for k in g:
            for day, sh in (seen[k].get(c) or {}).items():
                days[day] = max(days.get(day, 0.0), sh)
        return len(days), (max(days.values()) if days else 0.0)
    order = sorted(groups, key=lambda g: (-sum(stat(g, c)[0] for c in cities), _pair_pos(g[0])))
    shown, more = order[:8], max(0, len(order) - 8)
    head = "".join(f"<th>{_e(city_code(c))}</th>" for c in cities)
    rows = []
    n_bb = int((_check_cfg_safe() or {}).get("pair_common", 2000))
    for g in shown:
        cells = []
        for c in cities:
            n, b = stat(g, c)
            sty = ("background:#97c459" if n >= 3 else "background:#c0dd97" if n == 2 else
                   "background:#eaf3de" if n == 1 else "color:#d1d5db")
            cells.append(f"<td style='{sty}'>{f'{n} · {100 * b:.0f} %' if n else '–'}</td>")
        lo, hi = _pair_pos(g[0])[0], max(_pair_pos(k)[1] for k in g)
        weak = all(k in backbone for k in g)
        tip = ("pairs seen in this stretch: " + "; ".join(
            k.replace("+", " + ") + (" (backbone mutation)" if k in backbone else "") for k in g)
            + (f" · backbone mutation: inherited from a deep ancestor and carried by more than "
               f"{n_bb:,} lineages, so weaker evidence; a sample also needs a stretch without one"
               if any(k in backbone for k in g) else ""))
        lab = (f"{lo:,}–{hi:,} · {len(g)} {_pl(len(g), 'pair')}"
               + (" · with a backbone mutation" if weak else ""))
        best = max(g, key=lambda k: sum(len(v) for v in seen[k].values()))
        q = quote_plus(best.replace("+", ","))
        rows.append(f"<tr><td class='rl' style='{'color:#9ca3af' if weak else ''}' data-tip='{_e(tip)}'>"
                    f"{_e(lab)}</td>{''.join(cells)}<td class='rr'>"
                    + _ext(_CLIN_WORLD.format(muts=q), "World", f"clinical sequences carrying {best.replace('+', ' and ')}, worldwide")
                    + " " + _ext(_CLIN_EUROPE.format(muts=q), "Europe", "the same, Europe only")
                    + "</td></tr>")
    return ("<div class='dh'>Its coocs seen, by stretch of genome <span class='dim'>(its own mutations "
            "together on the same reads; pairs closer than a read count once; cell = samples where "
            "seen · best share; hover a stretch for its pairs)</span></div>"
            f"<table class='wkt'><tr><th style='text-align:left'>stretch</th>{head}"
            "<th style='text-align:left'>clinical</th></tr>"
            f"{''.join(rows)}</table>"
            + (f"<div class='dim' style='font-size:11px'>+ {more} more {_pl(more, 'stretch')}</div>"
               if more else ""))


def _why_short(x):
    """why a mutation couldn't be answered, in a few characters (2026-10-27)"""
    a = _answer_text(x)
    if "needs 100" in a or "counted reads" in a or not x.get("cov"):
        return "<100"
    if "between" in a:
        return "1–5 %"
    if "neighbour" in a or "look like" in a:
        return "link"
    return "?"


def _mut_grid(v, per_city, cities):
    """Its mutations × city, share of counted reads (2026-10-27); its own
    first, the ones it inherits tagged; 8 shown, the rest behind 'show all'."""
    det = {c: {x.get("marker"): x for x in ((per_city.get(c) or {}).get("detail") or [])} for c in cities}
    marks = list(dict.fromkeys(m for c in cities for m in det[c]))
    if not marks:
        return ""
    any_d = next((d for d in per_city.values() if d and d.get("own") is not None), {}) or {}
    own, par = set(any_d.get("own") or []), any_d.get("parent") or ""
    pos = {m: i for i, m in enumerate(marks)}
    marks.sort(key=lambda m: (0 if (not own or m in own) else 1, pos[m]))
    head = "".join(f"<th>{_e(city_code(c))}</th>" for c in cities)

    def row(m):
        tag = ("" if not own else "<span class='dim'> own</span>" if m in own
               else f"<span class='dim'> {_e(par or 'parent')}'s</span>")
        cells = []
        for c in cities:
            x = det[c].get(m)
            if not x:
                cells.append("<td></td>")
                continue
            st, f = x.get("status"), x.get("freq")
            txt = (_why_short(x) if st not in ("present", "absent") or f is None
                   else f"{100 * f:.0f} %")
            sty = ("background:#97c459" if st == "present" and (f or 0) >= .2 else
                   "background:#c0dd97" if st == "present" else
                   "color:#9ca3af" if st == "absent" else "color:#9ca3af;font-style:italic")
            cells.append(f"<td style='{sty}' data-tip='{_e(city_name(c))}: {_e(_share_text(x))} · "
                         f"{_e(_answer_text(x))}'>{txt}</td>")
        lk = _ext(_clin_link([m]), m, "which lineages carry it in clinical sequences (CoV-Spectrum)")
        return f"<tr><td class='rl mono'>{lk}{tag}</td>{''.join(cells)}</tr>"
    first, rest = marks[:8], marks[8:]
    t1 = f"<table class='wkt'><tr><th></th>{head}</tr>{''.join(row(m) for m in first)}</table>"
    t2 = (f"<details><summary class='dim' style='font-size:11.5px;cursor:pointer'>show all "
          f"{len(marks)}</summary><table class='wkt'><tr><th></th>{head}</tr>"
          f"{''.join(row(m) for m in rest)}</table></details>" if rest else "")
    return ("<div class='dh'>Its mutations, share of counted reads <span class='dim'>(couldn't "
            "answer: &lt;100 = too few counted reads · 1–5 % = in between · link = its reads don't "
            "look like the lineage nearby; hover for reads)</span></div>" + t1 + t2)


def _detail_row(gid, v, per_city, cities, sel, n_cols, found=False, clin=None, parent=None):
    """The row's details (2026-10-27, reworked as mocked): links, one column
    per week (wastewater, Europe / World clinical, each city's samples), its
    coocs seen per city, its mutations × city."""
    ds = per_city.get(sel) or {}
    own = [m for m in (next((d.get("own") for d in per_city.values() if d and d.get("own")), None)
                       or [x.get("marker") for x in ds.get("detail") or []]) if m]
    links = ("<div style='font-size:12.5px;margin-bottom:6px'>"
             + (_clin_links(own[:12]) + " · " if own else "")
             + _ext(_ww_link(v, sel), f"wastewater, {city_name(sel)} (GenSpectrum)",
                    "mutation heatmaps of the lineage and its sublineages")
             + "<span class='dim'> · each mutation per week: Signal over time (Look-up)</span></div>")
    al, ct = _absent_line(per_city), _counted_text(per_city, cities)
    allc = (("<div style='font-size:12.5px;margin:6px 0;white-space:normal;max-width:900px'><b>All cities</b> · "
             + " · ".join(_e(x) for x in (ct, al) if x) + "</div>") if (al or ct) else "")
    body = (links + _week_block(v, per_city, cities, clin, parent) + allc
            + _cooc_table(per_city, cities) + _mut_grid(v, per_city, cities))
    return f"<tr class='dt g-{gid} hid'><td colspan='{n_cols}'>{body}</td></tr>"


def _panel_cells(v, per_city, cities, sel):
    res = [d for d in per_city.values() if d]
    if res and all(_vote_counts(d)[0] == 0 for d in res):
        # no ★ marker: the same in every city, so said once (2026-10-20)
        return (f"<td class='nomk' colspan='{len(cities)}'>no mutation of its own — all are "
                "shared with other lineages in the panel</td>")
    out = []
    for c in cities:
        d = per_city.get(c)
        if not d:
            out.append(_empty(f"{city_name(c)}: no result", c == sel))
            continue
        style, label, _dot, _ink, share = _scale(d)
        out.append(_cell(label, style,
                         _cell_tip(v, city_name(c), d), c == sel))
    return "".join(out)


def _vote_counts(d):
    n_mk, n_meas, n_p = d.get("n_markers", 0), d.get("n_measured", 0), d.get("n_present", 0)
    n_meas = min(n_meas, n_mk)
    n_p = min(n_p, n_meas)
    return n_mk, n_p, n_meas - n_p, n_mk - n_meas


def _vote_text(d):
    """'1 present · 3 absent · 1 too few reads (of 5 ★ markers)' — exactly what
    the verdict voted on: more present than absent (and >= check.min_present
    present) -> present, more absent -> absent, a tie -> mixed (majority,
    2026-10-20; until then >= 75 % present, <= 25 % absent)."""
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


def _num():
    """Every number shown in a text, from the config (2026-10-08): mc reads to
    measure a marker, pf / af present / absent % (check.*); er / es / ed the
    evidence day rule and days to confirm (evidence.*)."""
    mc, pf, af = _thr()
    try:
        from process.scanner import EVIDENCE_MIN_DAYS, EVIDENCE_MIN_READS, EVIDENCE_MIN_SHARE
    except Exception:
        EVIDENCE_MIN_READS, EVIDENCE_MIN_SHARE, EVIDENCE_MIN_DAYS = 20, 0.005, 2
    return {"mc": mc, "pf": _pct(pf), "af": _pct(af), "mp": _min_present(),
            "er": int(EVIDENCE_MIN_READS),
            "es": _pct(float(EVIDENCE_MIN_SHARE) * 100), "ed": int(EVIDENCE_MIN_DAYS)}


_WORD = {"present": "present", "absent": "absent", "mixed": "mixed",
         "few": "too few markers measured"}


def _cov_factor():
    """A sample counts at most this × its median depth (check.cov_factor)."""
    try:
        from process.cooc import _check_cfg
        return float(_check_cfg().get("cov_factor", 10))
    except Exception:
        return 10.0


def _cap_reads():
    """Reads a sample counts at most per marker (check.cap_reads)."""
    try:
        from process.cooc import _check_cfg
        return int(_check_cfg().get("cap_reads", 1000))
    except Exception:
        return 1000


def _min_present():
    """Markers that must be present for "present" (check.min_present)."""
    try:
        from process.cooc import _check_cfg
        return int(_check_cfg().get("min_present", 2))
    except Exception:
        return 2


def _list_max(key, default):
    """lists.* caps from the config."""
    try:
        from utils.config import get_cooc_setting
        return int(get_cooc_setting(f"lists.{key}", default) or default)
    except Exception:
        return default


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
def _day_marks(top, low=None):
    """mark -> (square style, meaning). top = colour of a present day, low =
    of an absent one: the same colours as the cells (orange for a panel
    variant, pale for a found lineage)."""
    n = _num()
    return {
        # lineages: the ★ check's vote on that day's reads
        # (process.variant_explorer.check_in_data)
        "present": (f"background:{top};", f"★ markers present (≥ {n['mp']} markers at ≥ "
                                          f"{n['pf']} % of ≥ {n['mc']} reads)"),
        "mixed": (f"background:{MIXED_BG};", "★ markers disagree"),
        "weak": (f"background:{MIXED_BG};", "★ markers disagree"),
        "thin": ("background:#fff;border:1.5px dotted #6b7280;",
                 "too few reads to tell from this sample alone (the cell adds up the "
                 "reads of all samples in the chosen dates)"),
        "absent": (f"background:{low};" if low else D_PALE, f"★ markers absent (< {n['af']} % of ≥ {n['mc']} reads)"),
        # novel combinations (scanner timeline)
        "day": (f"background:{top};", f"passes the evidence rule (≥ {n['er']} reads, ≥ {n['es']} % of the sample's reads)"),
        "seen": (f"background:{top}55;", "seen, below the evidence rule"),
        "gone": ("background:#d1d5db;", f"≥ {n['mc']} reads at its positions, not there"),
        "uncovered": ("background:#fff;border:1.5px dotted #6b7280;",
                      f"under {n['mc']} reads at its markers — too few reads to tell"),
        "unknown": ("background:#e5e7eb;", "no coverage data (re-run)"),
    }


def _cal_max():
    """Squares shown at most (check.calendar_samples): the newest ones."""
    try:
        from process.cooc import _check_cfg
        return int(_check_cfg().get("calendar_samples", 12))
    except Exception:
        return 12


# a low-coverage sample (process.sample_qc, 2026-10-13): its square keeps its
# colour, crossed by thin diagonal lines — counted, but it weighs little
_LOWQ = ("background-image:repeating-linear-gradient(135deg,rgba(255,255,255,.85) 0 1.5px,"
         "transparent 1.5px 3.5px);")
# the same stripes on a city cell, lighter so its count stays readable
_LOWQ_CELL = ("background-image:repeating-linear-gradient(135deg,rgba(255,255,255,.3) 0 2px,"
              "transparent 2px 6px);text-shadow:0 0 2px #14532d,0 0 1px #14532d;")
_LOWQ_EMPTY = ("background-image:repeating-linear-gradient(135deg,#9ca3af 0 1px,"
               "transparent 1px 3.5px);")


def _qual_bits(quality, d):
    """(extra square style, hover prefix) for one sample's coverage quality."""
    try:
        from process.sample_qc import quality_text, well_covered
    except Exception:
        return "", ""
    q = (quality or {}).get(d)
    ok = well_covered(q)
    if ok is None:
        return "", ""
    return ("" if ok else None), quality_text(q)


def _day_cal(tl, top, low=None, detail=None, quality=None):
    """One square per sample in the chosen dates, oldest left, newest right;
    the newest check.calendar_samples when there are more (the text after
    the squares counts all of them).
    A sample with too few reads is a dotted square."""
    tl = [(d, m) for d, m in (tl or []) if m != "unknown"]
    if not tl:
        return ""
    n = _cal_max()
    last = tl[-n:]          # the text after the squares counts all samples
    mk = _day_marks(top, low)
    detail = detail or {}
    sq = []
    for d, m in last:
        style, meaning = mk.get(m, mk["unknown"])
        extra, qtxt = _qual_bits(quality, d)
        if extra is None:              # low coverage: stripes over its colour
            extra = _LOWQ_EMPTY if m in ("thin", "uncovered") else _LOWQ
            qtxt += (" — the sample is there, but each part of the genome comes from "
                     "few virus molecules, so its shares are snapshots, not the city's mix")
        tip = (_short_date(d) + ": " + (qtxt + " · " if qtxt else "") + meaning
               + (" — " + detail[d] if detail.get(d) else ""))
        sq.append(f"<i style='{style}{extra}' data-tip='{_e(tip)}'></i>")
    sq = "".join(sq)
    pad = "".join("<i style='background:transparent'></i>" for _ in range(n - len(last)))
    return f"<span class='cal' style='margin-left:0;margin-right:8px'>{pad}{sq}</span>"


def _cal_phrase(tl, word="present"):
    """'present on 4 of 9 covered days · in the latest sample' — what the
    calendar says, in a few words."""
    if not tl:
        return ""
    on = {"present", "day"}
    cov = [(d, m) for d, m in tl if m != "unknown"]
    hits = [d for d, m in cov if m in on]
    if not cov or all(m in ("uncovered", "thin") for _d, m in cov):
        return f"too few reads in all {len(cov)} {_pl(len(cov), 'sample')}"
    if not hits:
        mid = sum(1 for _d, m in cov if m in ("mixed", "weak", "seen"))
        thin = sum(1 for _d, m in cov if m in ("thin", "uncovered"))
        if not mid and not thin:
            return f"absent all {len(cov)} {_pl(len(cov), 'sample')}"
        if word != "present":
            return f"weak {mid} of {len(cov)} {_pl(len(cov), 'sample')}, never {word}"
        ab = len(cov) - mid - thin
        parts = ([f"{ab} absent"] if ab else []) + ([f"{mid} mixed"] if mid else []) \
            + ([f"{thin} too few reads"] if thin else [])
        return f"never present in {len(cov)} {_pl(len(cov), 'sample')}: " + " · ".join(parts)
    after = sum(1 for d, _m in cov if d > hits[-1])
    tail = "in the latest sample" if not after else f"last {_short_date(hits[-1])}"
    return f"{word} {len(hits)} of {len(cov)} {_pl(len(cov), 'sample')} · {tail}"


def _star_phrase(tl):
    """The ★ check's days, named as such so they can't be mixed up with the
    deep scan's evidence days."""
    return "★ markers: " + _cal_phrase(tl)


def _panel_warn(per_city, in_panel):
    """'in your panel · present in none of 6 cities: absent in 4, not called
    in 2 · ' for a lineage of the panel that no city calls present — the same
    whichever city is selected (2026-10-20: the warning the orange cells
    gave, now said once, with the counts)."""
    if not in_panel:
        return ""
    res = [x for x in per_city.values() if x and _vote_counts(x)[0]]
    if not res or any(_state(x) == "present" for x in res):
        return ""
    n_abs = sum(1 for x in res if _state(x) == "absent")
    rest = len(res) - n_abs
    parts = ([f"absent in {n_abs}"] if n_abs else []) + ([f"not called in {rest}"] if rest else [])
    return (f"<b style='color:{AMBER_T}'>in your panel</b> · present in none of {len(res)} "
            f"{_pl(len(res), 'city')}: {', '.join(parts)} · ")


# ── the overview column (2026-10-20): the whole lineage over all cities, in
#    plain words, the same whichever city is selected. The per-city words are
#    in the details. Returns (html, tip, dot colour, name colour, dot kind). ──
def _ovb(text, style="color:#4b5563;"):
    return f"<span class='ovb' style='{style}'>{text}</span>"


def _city_seen(d):
    """process.national.city_seen, with the cell's own vote for results that
    carry counts but no state (2026-10-26)."""
    pev = (d or {}).get("pairs") or {}
    if pev.get("measured") or pev.get("positive"):
        return bool(pev.get("positive"))
    return _state(d or {}) == "present"


def _panel_overview(per_city, cities, wk=None):
    res = [per_city[c] for c in cities if per_city.get(c)]
    if not res:
        return (_ovb("<b>In panel</b> · no result"), "", GREY, "#6b7280", "dashed")
    if all(_vote_counts(x)[0] == 0 for x in res):
        return (_ovb("<b>In panel</b> · can't be checked", "color:#6b7280;"),
                "no mutation of its own: all its mutations are shared with other lineages in "
                "the panel, so its share in Deconvolution can't be checked here",
                GREY, "#6b7280", "filled")
    if _panel_warn(per_city, True):
        return (_ovb("<b>In panel</b> · not seen in any city",
                     f"background:#fffbeb;color:{AMBER_T};"),
                _panel_warn(per_city, True).replace("<b style='color:" + AMBER_T + "'>", "")
                .replace("</b>", "").rstrip(" ·"),
                AMBER, AMBER_T, "filled")
    meas = [c for c in cities if per_city.get(c) and _vote_counts(per_city[c])[0]]
    seen = [c for c in meas if _city_seen(per_city[c])]
    # a city seen only on thin samples says so (2026-10-26)
    def _tag(c):
        pev = per_city[c].get("pairs") or {}
        return city_code(c) + (" thin" if pev.get("positive") and not pev.get("positive_well")
                               else "")
    where = f" ({', '.join(_tag(c) for c in seen)})" if 0 < len(seen) <= 3 else ""
    # not confirmed in the newest weeks (2026-10-27): information only
    note, tip = "", ""
    try:
        from process.national import row_status
        rs = row_status(wk)
    except Exception:
        rs = None
    if rs and wk["weeks"][-1]["status"] != "confirmed":
        note = (f" · not confirmed since {_short_date(rs['last_confirmed'])}"
                if rs.get("last_confirmed") else " · not confirmed")
        tip = (f"not confirmed in the 4 weeks to {_short_date(wk['weeks'][-1]['end'])}"
               + (f" (last confirmed {_short_date(rs['last_confirmed'])})" if rs.get("last_confirmed") else "")
               + " · information only")
    where_all = ", ".join(_tag(c) for c in seen)
    tip = (f"In panel: seen in {len(seen)} of {len(meas)} cities"
           + (f" ({where_all})" if where_all else "") + ("\n" + tip if tip else ""))
    return (_ovw(_ovb(f"<b style='color:#1f2328'>In panel</b> · {len(seen)} of "
                      f"{len(meas)}{note}"), _weeks_strip(wk)),
            tip, GREEN if seen else "#d1d5db", "#15803d" if seen else "#6b7280", "filled")


def _absent_line(chk):
    """The same line for every row (2026-10-27): its mutations on the samples
    where its coocs were seen, all cities together (process.national.grade):
    'its mutations on its 13 samples with coocs: 2 measured, 0 absent'. A
    lineage carries all its mutations; several absent on its own samples
    means the reads are probably a close relative (PJ.2.1, Sep-Oct: 4 of 19,
    all inherited from PJ.2)."""
    try:
        from process.national import grade
        g = grade({c: d for c, d in (chk or {}).items() if d})
    except Exception:
        g = None
    if not g or not g.get("n_measured"):
        return ""
    n = g.get("n_counted")
    where = (f"on its {n} {_pl(n, 'sample')} with coocs" if n
             else "on all samples (no sample with coocs)")
    miss, inh = g.get("missing") or [], g.get("missing_inherited") or []
    par = next(((d or {}).get("parent") for d in (chk or {}).values() if (d or {}).get("parent")), "")
    txt = f"its mutations {where}, all cities together: {g['n_measured']} measured, {len(miss)} absent"
    if miss:
        txt += (": " + ", ".join(miss[:8]) + (" …" if len(miss) > 8 else "")
                + (f" (all inherited from {par})" if inh and len(inh) == len(miss) and par else ""))
    return txt + f" (absent = under {(_check_cfg_safe() or {}).get('absent_freq', 0.01) * 100:g} % of the reads)"


def _cell_word(d):
    """A city's cell in words (2026-10-27), as the legend names the colours."""
    d = d or {}
    pev = d.get("pairs") or {}
    if pev.get("positive"):
        return "coocs" if pev.get("positive_well") else "coocs on thin samples"
    sm = d.get("state_markers") or d.get("state") or (_state(d) if d else None)
    if d.get("own_missing"):
        return f"only {d.get('parent') or 'its parent'}'s mutations"
    if pev.get("measured") and sm == "present":
        return "single mutations only"
    if sm == "absent":
        return "not seen"
    return "undecided" if pev.get("measured") or sm in ("mixed", "present") else "too few reads"


def _finding_status(v, f, cities, has_children=False):
    """'new' / 'covered' / 'named' / 'folded' / 'suggested' for a lineage the
    reads suggested (a deep scan find), and the reason for folded / suggested
    (2026-10-27).

    new        its own coocs in at least one city (the cell's rule,
               process.national.city_seen): Watch, or Add to panel
    folded     no coocs of its own anywhere — or coocs on a single sample with
               half or more of its measured mutations clearly missing pooled
               (BA.3.2 Genève July 32 of 32, XFW.1.1.1 7 of 9): listed in one
               closed line under the tree, with the reason
    suggested  the same, but other rows hang below it in the tree: stays in
               place, grey, so no branch loses its children
    Until 2026-10-27 "new" followed + Add (the marker vote) and the rest was
    "possible" or, when its own mutations were not seen wherever the reads
    suggested it, "folded"."""
    st = f.get("status")
    if st == "family":
        return "covered", ""
    if st == "named":
        return "named", ""
    pc = f.get("per_city") or {}
    chk = {c: (pc.get(c) or {}).get("check") for c in cities if (pc.get(c) or {}).get("check")}
    seen = [c for c in cities if _city_seen(chk.get(c) or {})]
    n_pos = sum(len(((chk.get(c) or {}).get("pairs") or {}).get("positive") or []) for c in seen)
    if seen:
        g = None
        try:
            from process.national import grade
            g = grade(chk)
        except Exception:
            g = None
        g = g or {}
        miss, n_meas = len(g.get("missing") or []), int(g.get("n_measured") or 0)
        # half or more (2026-10-27; was "more than half": XFG.6.3.3, 2 of 4, on
        # the same Zürich sample that folded its child SW.2.1, 4 of 5)
        if n_pos == 1 and n_meas and 2 * miss >= n_meas:
            why = (f"coocs on 1 sample in {city_code(seen[0])} · {miss} of {n_meas} "
                   "mutations absent")
            return ("suggested" if has_children else "folded"), why
        return "new", ""
    sugg = {c: int((pc.get(c) or {}).get("days") or 0) for c in cities
            if not (pc.get(c) or {}).get("xcheck") and int((pc.get(c) or {}).get("days") or 0) > 0}
    top = max(sugg, key=sugg.get) if sugg else None
    parts = []
    if top:
        n = sugg[top]
        parts.append((f"suggested {n}× in {city_code(top)}" if n > 1
                      else f"suggested on 1 sample in {city_code(top)}")
                     + f" · {_cell_word(chk.get(top))}{' there' if _cell_word(chk.get(top)) == 'not seen' else ''}")
    # one more city that isn't a plain "not seen", if any
    other = [c for c in cities if c != top and c in chk
             and _cell_word(chk[c]) not in ("not seen", "too few reads")]
    if other:
        parts.append(f"{_cell_word(chk[other[0]])} in {city_code(other[0])}")
    if not parts:
        parts.append("its own mutations not seen anywhere" if chk else "no check result")
    return ("suggested" if has_children else "folded"), " · ".join(parts)


def _weeks_text(wk):
    """the weekly grade in words, for the hover (process.national.weekly)"""
    need, n = wk.get("need", 3), wk.get("in_a_row", 0)
    days = (_check_cfg_safe() or {}).get("add_days", 28)
    tail = " · ".join(f"{_short_date(w['end'])} {w['status']}" for w in (wk.get("weeks") or [])[-need - 2:])
    return (f"+ Add when {need} weeks in a row are confirmed, each on the {days} days up to that "
            f"Sunday: {n} of {need}" + (f" ({tail})" if tail else "")
            + (f"; these dates fit {len(wk.get('weeks') or [])}, choose {days + 7 * (need - 1)} days "
               "or more" if wk.get("short") else ""))


def _short_reason(f, cities):
    """What keeps a found lineage from + Add, in a few words (2026-10-27), the
    first that applies: thin samples only, 1 city, 1 week, then the weeks in
    a row (mutations absent / not clear on its samples only when
    check.confirm_use_missing / confirm_use_pooled are on; off since
    2026-10-27)."""
    pc = f.get("per_city") or {}
    chk = {c: (pc.get(c) or {}).get("check") for c in cities if (pc.get(c) or {}).get("check")}
    g = None
    # the newest 4-week window (2026-10-27): the status is about the last
    # weeks, not the whole chosen range
    wks = (f.get("weekly") or {}).get("weeks") or []
    if wks and wks[-1].get("grade"):
        g = wks[-1]["grade"]
    if g is None:
        try:
            from process.national import grade
            g = grade(chk) if chk else None
        except Exception:
            g = None
    g = g or f.get("grade") or {}
    c = _check_cfg_safe() or {}
    if g and c.get("confirm_use_missing", False) and not g.get("ok_missing", True):
        # its mutations on the samples with its coocs, all cities together;
        # absent = under check.absent_freq (1 %) (2026-10-27)
        on = (f" on its {g['n_counted']} {_pl(g['n_counted'], 'sample')}" if g.get("n_counted")
              else " (all samples)")
        return f"{len(g.get('missing') or [])} of {g.get('n_measured', 0)} mutations absent{on}"
    well, thin = g.get("cities") or [], g.get("cities_thin") or []
    if not well and thin:
        return "thin samples only"
    if len(well) < c.get("confirm_min_cities", 2):
        return "1 city" if well else "not seen"
    if g and g.get("weeks", 0) < c.get("confirm_min_weeks", 2):
        return "1 week" if g.get("weeks") else "no well-covered week"
    if g and c.get("confirm_use_pooled", False) and not g.get("pooled_seen", True):
        return ("not clear on its samples" if g.get("n_counted")
                else "not clear with cities pooled")
    wk = f.get("weekly")
    if not wk:
        return "re-run for weeks"
    return f"{wk.get('in_a_row', 0)} of {wk.get('need', 3)} weeks"


def _counted_text(chk, cities):
    """'its 6 samples with coocs: 13 Sep GE · 3 Oct BS, GR, TI, ZH · 4 Oct ZH'
    (2026-10-27): when they were, newest last, so a run of recent samples
    (an arrival) and scattered single ones read differently.
    chk: {city: check_in_data result}."""
    by_day = {}
    for c in cities:
        pev = ((chk or {}).get(c) or {}).get("pairs") or {}
        for d in pev.get("positive") or []:
            by_day.setdefault(str(d)[:10], []).append(city_code(c))
    if not by_day:
        return ""
    n = sum(len(v) for v in by_day.values())
    return (f"its {n} {_pl(n, 'sample')} with coocs, oldest first: "
            + " · ".join(f"{_short_date(d)} {', '.join(sorted(cs))}" for d, cs in sorted(by_day.items())))


def _short_hover(f, cities, rs, wk, why):
    """Three short lines for a found row's hover (2026-10-27): the status and
    its weeks, where it was seen in the chosen dates, and absent mutations
    if any; everything else is in the details."""
    weeks = (wk or {}).get("weeks") or []
    conf = [w["end"] for w in weeks if w["status"] == "confirmed"]
    if rs and rs["kind"] == "add":
        last = [_short_date(e) for e in conf[-rs["in_a_row"]:]]
        l1 = f"Add to panel: confirmed {rs['in_a_row']} weeks in a row ({', '.join(last)})"
    elif rs and rs["kind"] == "fading":
        l1 = f"Fading: last confirmed {_short_date(rs['last_confirmed'])}, not since"
    elif rs and rs["kind"] == "earlier":
        l1 = (f"Earlier only: last seen {_short_date(rs['last_found'])}, never confirmed"
              if rs.get("last_found") else "Earlier only: never confirmed")
    elif rs and weeks and weeks[-1]["status"] == "confirmed":
        last = [_short_date(e) for e in conf[-rs["in_a_row"]:]]
        l1 = (f"Watch: confirmed {rs['in_a_row']} of the {rs['need']} weeks needed "
              f"({', '.join(last)})")
    else:
        l1 = (f"Watch: not confirmed in the 4 weeks to {_short_date(weeks[-1]['end'])} · "
              f"{_short_reason(f, cities)}" if weeks else f"Watch: {_short_reason(f, cities)}")
    lines = [l1, f"in the chosen dates: {why}" if why else ""]
    pc = f.get("per_city") or {}
    chk = {c: (pc.get(c) or {}).get("check") for c in cities if (pc.get(c) or {}).get("check")}
    try:
        from process.national import grade
        g = grade(chk) if chk else None
    except Exception:
        g = None
    miss = (g or {}).get("missing") or []
    if miss:
        lines.append(f"{len(miss)} of its {g.get('n_measured', 0)} mutations absent on its "
                     f"samples: {', '.join(miss[:4])}{' …' if len(miss) > 4 else ''} → details")
    return "\n".join(x for x in lines if x)


def _ovw(label, strip):
    """status left, weeks squares right-aligned in the overview cell, so the
    squares line up like a column (2026-10-27)"""
    return f"<span class='ovw'>{label}{strip}</span>" if strip else label


def _weeks_strip(wk, keep=12):
    """The weekly grades as small squares, newest right (2026-10-27): each =
    the 4 weeks up to that Sunday; solid confirmed, light found, empty not
    found. More than `keep` weeks: '⋯' first, its hover sums up the older."""
    weeks = (wk or {}).get("weeks") or []
    if len(weeks) < 2:
        return ""
    sty = {"confirmed": f"background:{GREEN}", "watch": "background:#c0dd97",
           "none": "background:#fff;box-shadow:inset 0 0 0 1px #d1d5db"}
    word = {"confirmed": "confirmed", "watch": "found, not confirmed", "none": "not found"}
    days = (_check_cfg_safe() or {}).get("add_days", 28)
    old, shown = weeks[:-keep] if len(weeks) > keep else [], weeks[-keep:]
    out = []
    if old:
        oc = [w for w in old if w["status"] == "confirmed"]
        tip = (f"{len(old)} older {_pl(len(old), 'week')}: {len(oc)} confirmed"
               + (f", last confirmed {_short_date(oc[-1]['end'])}" if oc else ""))
        out.append(f"<span class='wkm' data-tip='{_e(tip)}'>⋯</span>")
    for w in shown:
        g = w.get("grade") or {}
        where = ", ".join(city_code(c) for c in g.get("cities") or [])
        tip = (f"{days} days to {_short_date(w['end'])}: {word.get(w['status'], '')}"
               + (f" · seen in {where}" if where else ""))
        out.append(f"<i style='{sty.get(w['status'], '')}' data-tip='{_e(tip)}'></i>")
    return f"<span class='wks'>{''.join(out)}</span>"


def _frames_note(findings, panel_weekly, dates):
    """'cells: chosen dates … · status: last weeks, up to …' (2026-10-27)"""
    if not dates:
        return ""
    ends = [w["end"] for wk in list((panel_weekly or {}).values())
            + [f.get("weekly") for f in (findings or {}).values()] if wk
            for w in (wk.get("weeks") or [])[-1:]]
    last = max(ends) if ends else None
    return (f"<div class='dim' style='font-size:12.5px;margin:4px 0 2px'>cells: the chosen dates, "
            f"{_e(_short_date(dates[0]))} – {_e(_short_date(dates[1]))} · status: the last weeks"
            + (f", up to {_e(_short_date(last))}" if last else "")
            + " (each square: the 4 weeks up to a Sunday)</div>")


def _add_note(findings):
    """One line above the table when the dates are too short to judge + Add
    (2026-10-27), instead of the same words on every row."""
    wks = [f.get("weekly") for f in (findings or {}).values()
           if f.get("status") not in ("named", "family") and f.get("weekly")]
    if not wks or not any(w.get("short") for w in wks):
        return ""
    c = _check_cfg_safe() or {}
    need, days = int(c.get("add_weeks", 3)), int(c.get("add_days", 28))
    fit = max(len(w.get("weeks") or []) for w in wks)
    return (f"<div class='dim' style='font-size:12.5px;margin:4px 0 6px'>Found lineages · + Add needs "
            f"{need} weeks in a row confirmed: these dates allow {fit} "
            f"{_pl(fit, 'week')}, choose {days + 7 * (need - 1)} days or more to judge</div>")


def _new_grade(f, cities):
    """(confirmed, reason, hover) for a "New" row, from all cities pooled
    (process.national.grade). Confirmed: >= 2 well-covered cities AND >= 2
    weeks (2026-10-27; the pooled and absent checks are switches, off);
    otherwise early, with the first condition it misses. Absent mutations
    are always named in the hover."""
    pc = f.get("per_city") or {}
    chk = {c: (pc.get(c) or {}).get("check") for c in cities if (pc.get(c) or {}).get("check")}
    g = None
    if chk:                     # recomputed: a grade saved before 26 Oct used the old rule
        try:
            from process.national import grade
            g = grade(chk)
        except Exception:
            g = None
    if g is None:
        g = f.get("grade")
    seen = [c for c in cities if _city_seen(chk.get(c) or {})]
    if not g:
        codes = ", ".join(city_code(c) for c in seen)
        return False, (f"seen in {codes}" if seen else "seen"), "no pooled answer (no samples)"
    if g and "cities_thin" not in g:                 # a grade saved before 26 Oct
        g = dict(g, cities_thin=[c for c in seen if _thin(chk.get(c) or {})])
        g["cities"] = [c for c in seen if c not in g["cities_thin"]]
    c = _check_cfg_safe()
    well, thin = [x for x in seen if x in g["cities"]], [x for x in seen if x in g["cities_thin"]]
    n_c, wk = len(well), f"{g['weeks']} of {g['n_weeks']} {_pl(g['n_weeks'], 'week')}"
    codes = ", ".join(city_code(x) for x in well)
    thin_txt = (f" (+ {', '.join(city_code(x) for x in thin)} thin)" if thin else "")
    where = (f"seen in {n_c} {_pl(n_c, 'city')}" if n_c > 3 else (f"seen in {codes}" if well else
             ("seen only in thin samples" if thin else ""))) + (thin_txt if well else "")
    use_p, use_m = c.get("confirm_use_pooled", False), c.get("confirm_use_missing", False)
    rule = (f"confirmed needs {c.get('confirm_min_cities', 2)} cities with its coocs on a "
            f"well-covered sample and {c.get('confirm_min_weeks', 2)} weeks"
            + (", its mutations seen together on its samples" if use_p else "")
            + (f", at most {c.get('confirm_max_missing', 0.1) * 100:g} % of its mutations absent "
               f"(at least {c.get('confirm_missing_allow', 1)} allowed)" if use_m else ""))
    # absent mutations: always said, a check only with confirm_use_missing
    k, inh = len(g.get("missing") or []), g.get("missing_inherited") or []
    par = next(((chk.get(x) or {}).get("parent") for x in cities if (chk.get(x) or {}).get("parent")), "")
    inh_txt = ("" if not inh else (f", all inherited from {par}" if len(inh) == k and par else
               f", {len(inh)} inherited{' from ' + par if par else ''}"))
    on = (f"on its {g['n_counted']} {_pl(g['n_counted'], 'sample')} with coocs"
          if g.get("n_counted") else "on all samples")
    absent = (f" · note: {k} of its {g.get('n_measured', 0)} mutations absent {on}{inh_txt}: "
              f"{', '.join((g.get('missing') or [])[:8])} (under "
              f"{c.get('absent_freq', 0.01) * 100:g} % of the reads). A virus of this lineage "
              "carries all of them; on a few samples one can read 0 % because that sample's "
              "amplicon was copied from another virus, several suggest a close relative. "
              "Check before + Add." if k else "")
    if g["confirmed"]:
        return True, f"{where}, {wk}" if where else wk, (
            f"confirmed: in {n_c} well-covered {_pl(n_c, 'city')} and {wk}" + absent)
    early = "early: found, not confirmed yet · check the next samples · " + rule + absent
    if use_m and not g["ok_missing"]:
        return False, (f"{where} · " if where else "") + (
            f"{k} of its {g['n_measured']} mutations missing{inh_txt}"), early
    if use_p and not g["pooled_seen"]:
        return False, f"{where} · not clear on its samples", early
    if n_c < c.get("confirm_min_cities", 2):
        return False, (f"{where} only, {wk}" if well else f"{where}, {wk}"), early
    return False, f"{where}, {wk}", early


def _finding_overview(v, f, cities, current_panel, has_children=False):
    kind, why = _finding_status(v, f, cities, has_children)
    pc = f.get("per_city") or {}
    if kind == "covered":
        desc = f.get("descendants") or []
        by = ", ".join(desc[:2]) + (f" +{len(desc) - 2}" if len(desc) > 2 else "")
        return (kind, _ovb(f"<b>Covered</b> by {_e(by)}", "color:#9ca3af;"),
                f"its own mutations are on the reads of the lineages below it ({', '.join(desc[:6])}); "
                "add those instead", "#d1d5db", "#9ca3af", "ring")
    if kind == "named":
        return (kind, _ovb("not confirmed: reads also fit related lineages", "color:#9ca3af;"),
                "", GREY, "#9ca3af", "ring")
    if kind == "new":
        # a badge and the one thing holding it back (2026-10-27); the
        # sentence goes to the hover
        conf, why, tip = _new_grade(f, cities)
        wk = f.get("weekly")
        try:
            from process.national import row_status
            rs = row_status(wk)
        except Exception:
            rs = None
        hover = _short_hover(f, cities, rs, wk, why)
        strip = _weeks_strip(wk)
        if f.get("add_ready") or (rs and rs["kind"] == "add"):
            btn = ("<span class='added'>added · re-run needed</span>" if v in current_panel
                   else f"<button class='add' data-add='{_e(v)}'>+ Add</button>")
            n = (wk or {}).get("in_a_row", 0)
            return (kind, _ovw(_ovb(f"<b>Add</b> · {n} weeks in a row{btn}",
                                    "background:#fef2f2;color:#b91c1c;"), strip),
                    hover, RED, RED, "filled")
        if rs and rs["kind"] == "fading":
            return (kind, _ovw(_ovb(f"<b>Fading</b> · last confirmed {_e(_short_date(rs['last_confirmed']))}",
                                    "color:#6b7280;box-shadow:inset 0 0 0 1px #e5e7eb;"), strip),
                    hover, "#d1d5db", "#6b7280", "ring")
        if rs and rs["kind"] == "earlier":
            return (kind, _ovw(_ovb(f"<b>Earlier only</b> · last seen "
                                    f"{_e(_short_date(rs['last_found'])) if rs.get('last_found') else '–'}",
                                    "color:#6b7280;box-shadow:inset 0 0 0 1px #e5e7eb;"), strip),
                    hover, "#d1d5db", "#6b7280", "ring")
        # confirmed in the newest weeks, not yet in a row: say how far it got
        reason = (f"{rs['in_a_row']} of {rs['need']} weeks"
                  if rs and (wk or {}).get("weeks") and wk["weeks"][-1]["status"] == "confirmed"
                  else _short_reason(f, cities))
        return (kind, _ovw(_ovb(f"<b>Watch</b> <span style='color:#4b5563'>· {_e(reason)}</span>",
                                "background:#eff6ff;color:#1d4ed8;"), strip), hover, RED, RED, "ring")
    if kind in ("folded", "suggested"):
        # no coocs of its own (2026-10-27): grey, the reason in the cells' words
        return (kind, _ovb(_e(why), "color:#6b7280;"),
                "the deep scan suggested this name, but no city has coocs of its own mutations "
                "(or only one sample, with half or more of its mutations absent): not counted; if it is "
                "real, the next batch moves it up", "#d1d5db", "#9ca3af", "ring")


def _panel_evidence(per_city, sel, v="", gid="", in_panel=False):
    """One line for the selected city (2026-10-13: the per-sample squares
    moved into the row's details)."""
    d = per_city.get(sel)
    if not d:
        return "<span class='dim'>no result in this city</span>"
    n_mk, p, a, u = _vote_counts(d)
    if n_mk == 0:
        return ""                     # said once across the city cells (2026-10-20)
    tip = _cell_tip(v, city_name(sel), d)
    return (_panel_warn(per_city, in_panel)
            + f"<span data-tip='{_e(tip)}'>{_e(_summary(d))}</span>"
            + (_det_link(gid) if gid and d.get("timeline") else ""))


def _panel_tint(per_city, sel):
    """Row tint: the chosen city's cell colour (none when not measured)."""
    d = per_city.get(sel)
    if not d:
        return None
    _style, _label, dot, _ink, share = _scale(d)
    return None if share is None else dot


# a "family" row (an ancestor found through its descendants): neutral cells,
# so its ★ check doesn't read as a red "found" lineage (2026-10-14)
S_FAMILY = "background:#f3f4f6;color:#4b5563;box-shadow:inset 0 0 0 1px #d1d5db;"


def _finding_cells(v, f, cities, sel):
    pc = f.get("per_city", {}) or {}
    fam = f.get("status") == "family"
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
            if fam:
                style = S_FAMILY
            fi = ", ".join(city_name(x) for x in d.get("found_in") or [])
            out.append(_cell(label, style, _cell_tip(
                v, city_name(c), chk, f" · checked here because found in {fi}"), c == sel))
            continue
        if chk and n > 0:
            # the same measure as the panel: ★ markers present / measurable
            style, label, _dot, _ink, _share = _scale(chk, found=True)
            if fam:
                style = S_FAMILY
            out.append(_cell(label, style, _cell_tip(
                v, city_name(c), chk,
                f" · deep scan: evidence on {n} {_pl(n, 'sample')}{_ev_quality(d)}"),
                c == sel))
            continue
        if n == 0:
            out.append(_cell("0", S_NAMED, f"{v} · {city_name(c)}: not confirmed — reads point to "
                             f"{v}, but they also fit related lineages, so no sample has evidence "
                             f"specific to {v}", c == sel))
            continue
        tip = (f"{v} · {city_name(c)}: evidence on {n} {_pl(n, 'sample')} · "
               + (f"{len(stars)} ★ {_pl(len(stars), 'marker')} ({', '.join(stars[:5])})"
                  if stars else "no ★ marker")
               + (f" · {len(regs)} genome {_pl(len(regs), 'region')}: "
                  + ", ".join(f"{lo:,}–{hi:,}" for lo, hi in regs[:5]) if regs else ""))
        if n == 1 and d.get("after"):
            tip += " · 1 day " + _after_text(d["after"], html=False)
        out.append(_cell(n, S_FOUND if n >= _num()["ed"] else S_ONEDAY, tip, c == sel))
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
        t = (f"on {day}; {k} later {_pl(k, 'sample')}, none with ≥ {_num()['mc']} reads at its "
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
    "gone": ("background:#9ca3af;", f"covered (≥ {_num()['mc']} reads), not there"),
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


def _ev_quality(d):
    """', all low coverage' / ' (1 well covered)' / '' for the deep scan's
    evidence samples in a city (2026-10-13): evidence only from low-coverage
    samples can come from few virus molecules (PQ.16.1.1 in Basel)."""
    dates = d.get("ev_dates") or []
    q = (d.get("check") or {}).get("quality") or {}
    known = [x for x in dates if x in q]
    if not known:
        return ""
    try:
        from process.sample_qc import well_covered
    except Exception:
        return ""
    k = sum(1 for x in known if well_covered(q[x]))
    if k == 0:
        return ", all low coverage" if len(known) == len(dates) else ""
    return f" ({k} well covered)"


def _finding_evidence(f, sel, gid=""):
    pc = f.get("per_city", {}) or {}
    d = pc.get(sel)
    n_cities = sum(1 for x in pc.values() if int(x.get("days", 0) or 0) >= 1)
    if d is None:
        where = (f" · evidence in {n_cities} other {_pl(n_cities, 'city')}" if n_cities else "")
        return f"<span class='dim'>not seen in this city{where}</span>"
    n = int(d.get("days", 0) or 0)
    if d.get("xcheck"):
        chk = d.get("check") or {}
        fi = ", ".join(city_name(x) for x in d.get("found_in") or [])
        why = (f"not named by this city's deep scan; its ★ markers were checked here because "
               f"it was found in {fi}")
        return (f"<span data-tip='{_e(why)}'>{_e(_summary(chk))}</span>"
                + (_det_link(gid) if gid and chk.get("timeline") else ""))
    if n == 0:
        return ("<span class='dim'>not confirmed — reads point to it but also fit related "
                "lineages; no sample with specific evidence</span>")
    stars, regs = d.get("stars") or [], d.get("regions") or []
    facts = (f"deep scan: {n} {_pl(n, 'sample')} · "
             + (f"{len(stars)} ★ {_pl(len(stars), 'marker')}" if stars else "no ★ marker")
             + f" · {len(regs)} genome {_pl(len(regs), 'region')}"
             + (" · 1 day " + _after_text(d["after"], html=False)
                if n == 1 and d.get("after") else ""))
    tl = (d.get("check") or {}).get("timeline")
    if tl:
        nn = _num()
        scan_tip = (f"The deep scan counts a sample when reads with the lineage's mutations "
                    f"are ≥ {nn['er']} reads and ≥ {nn['es']} % of the sample's reads. The ★ "
                    f"samples are the ★ check's vote on its specific markers (≥ {nn['mp']} "
                    f"present at ≥ {nn['pf']} %): a different, stricter question. " + facts
                    + ". " + _REGION_TIP)
        txt = (f"<span>{_e(_summary(d.get('check') or {}))}</span>"
               f" <span class='dim' data-tip='{_e(scan_tip)}'>· deep scan: evidence on "
               f"{n} {_pl(n, 'sample')}{_e(_ev_quality(d))}</span>")
    else:                                        # results from before 2026-10-02
        txt = (f"{n} {_pl(n, 'sample')} · "
               + (f"{len(stars)} ★ {_pl(len(stars), 'marker')}" if stars else "no ★ marker")
               + f" · <span data-tip='{_e(_REGION_TIP)}'>{len(regs)} genome "
               f"{_pl(len(regs), 'region')}</span>")
        if n == 1 and d.get("after"):
            txt += _after_cal(d["after"]) or (" — " + _after_text(d["after"]))
    # the deep scan's own evidence is thin (2026-10-14: was a bold "weak: no ★
    # marker", which read as if the cells' ★ check had no marker)
    weak = []
    if not stars:
        weak.append("no mutation specific to it")
    if len(regs) <= 1:
        weak.append("one genome region only")
    if weak and n >= _num()["ed"]:
        txt += f"<span class='dim'> — deep scan evidence thin: {', '.join(weak)}</span>"
    return txt + (_det_link(gid) if gid and tl else "")


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
    nn = _num()
    cells = []
    for c in cities:
        n = days.get(c)
        if n is None:
            cells.append(_empty(f"{city_name(c)}: not seen", sel == c))
        elif n == 0:
            cells.append(_cell("0", S_NOMARK, f"{city_name(c)}: seen, but in no sample above "
                               f"the evidence rule ({nn['er']} reads and {nn['es']} % of the sample's reads)",
                               sel == c))
        else:
            cells.append(_cell(n, style_full if n >= nn["ed"] else style_one,
                               f"{city_name(c)}: {tip_word} on {n} {_pl(n, 'sample')} — "
                               f"{' '.join(muts)}", sel == c))
    m = " ".join(muts[:4]) + (" …" if len(muts) > 4 else "")
    return (f"<tr><td style='padding-left:14px'><span class='mono' style='color:{color}' "
            f"data-tip='{_e(' '.join(muts))}'>{_e(m)}</span></td><td class='ov'></td>"
            f"{''.join(cells)}<td class='ev'>{evidence}</td></tr>")


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
                                    + (f"{reads:,} reads, no sample above the evidence rule" if reads
                                       else "not seen"), city == sel))
                continue
            cells.append(_cell(n, S_NEAR if n >= _num()["ed"] else S_NEAR1,
                               f"{c['label']} · {city_name(city)}: {n} {_pl(n, 'sample')} "
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
        ev = (f"{n_sel} {_pl(n_sel, 'sample')} here · " if n_sel else
              "<span class='dim'>not in this city · </span>") + place
        if not gain:
            ev = f"<span class='dim'>{ev}</span>"
        out.append(f"<tr class='near g-{gid} hid'>{tc}<td class='ov'></td>{''.join(cells)}"
                   f"<td class='ev'>{ev}</td></tr>")
    return out


def build(cities, sel, rows, verdicts, findings, current_panel, ot=(),
          novel=(), broad=(), novel_rest=None, near=None, near_min_days=2,
          novel_hotspot=None, data_until=None, novel_info=None,
          broad_total=None, named_hidden=0, quality=None, panel_weekly=None,
          dates=None, clinical=None, pango_parent=None, sample_well=None) -> str:
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
      on >= near_min_days days) that unfolds them.
    broad_total: patterns too broad to name before the list was cut (broad is
      the cut list); named_hidden: found rows "named only" left out by the cap.
      Both show as "showing N of M" (lists.* caps in the config).
    quality: {city: (samples well covered, samples)} in the chosen dates
      (process.sample_qc), shown under each city code."""
    nn = _num()
    near_max = _list_max("near_changes_max", 8)
    near = near or {}
    data_until = data_until or {}
    ot = set(ot)
    n_c = len(cities)
    n_cols = n_c + 3          # tree · overview · cities · details (2026-10-20)
    quality = quality or {}

    def _qtip(c):
        if c not in quality:
            return ""
        k, n = quality[c]
        return (f" · {k} of {n} {_pl(n, 'sample')} in the chosen dates well covered "
                "(median depth and share of the genome at ≥ 100 reads); the others have low "
                "coverage, typical when little virus is in the water")

    # (2026-10-20) no "k/n" under the city code any more: it read like a vote;
    # the samples well covered are in the code's hover and the Evidence header
    # (2026-10-20) tree · overview (all cities, never changes) · cities · details;
    # the per-city words moved into the details, opened by clicking a cell
    head = ("<tr><th>Lineage</th><th style='padding-left:8px'>Overview · all cities</th>"
            + "".join(f"<th class='cc{' sel' if c == sel else ''}'><button data-city='{_e(c)}' "
                      f"data-tip='{_e(city_name(c))}"
                      + (f" · data until {_e(data_until[c])}" if c in data_until else "")
                      + _e(_qtip(c))
                      + f" — click to show this city in the details'>"
                      f"{_e(city_code(c))}</button></th>" for c in cities)
            # the last column takes the spare width, so the tree, the overview
            # and the cells stay together on the left (2026-10-20)
            + "<th style='width:100%'></th></tr>")
    body = []
    folded, folded_rows = [], []
    sname = city_name(sel)
    for r in rows:
        v, kind = r["node"], r["kind"]
        chips = " <span class='chip'>OT</span>" if v in ot and kind not in ("spine", "group") else ""
        if r.get("parents"):
            chips += (f"<span class='par' data-tip='{_e(r.get('parents_tip', ''))}'>"
                      f"{_e(r['parents'])}</span>")
        if kind == "group":
            body.append(f"<tr class='sp'>{_tree_cell(r, GREY, 'spine', f'''<span class='grp-root'>{_e(r['label'])}</span>''', r.get('tip', r['label']))}"
                        f"<td colspan='{n_c + 2}'></td></tr>")
            continue
        if kind == "spine":
            name = f"<span class='spl'>{_e(r['label'])}</span>{chips}"
            body.append(f"<tr class='sp'>{_tree_cell(r, GREY, 'spine', name, r['label'])}"
                        f"<td colspan='{n_c + 2}'></td></tr>")
            continue
        if kind in ("panel", "panel_ot"):
            pc = {c: verdicts[c][v] for c in cities if v in verdicts.get(c, {})}
            # the dot, the name and the overview follow all cities, not the
            # selected one (2026-10-20); amber only for "in panel, not
            # seen in any city"
            ov, ovtip, dcol, ncol, dk = _panel_overview(pc, cities, (panel_weekly or {}).get(v))
            tip = f"{v}: in the panel" + (f" · {ovtip}" if ovtip else "")
            name = f"<span class='nm' style='color:{ncol}'>{_e(v)}</span>{chips}"
            tc = _tree_cell(r, dcol, dk, name, tip)
            tint = None                     # no row tint (2026-10-20)
            dgid = "det-" + re.sub(r"[^A-Za-z0-9]", "_", v)
            cells = _panel_cells(v, pc, cities, sel)
            ov = f"<span data-tip='{_e(ovtip).replace(chr(10), '&#10;')}'>{ov}</span>" if ovtip else ov
            drow = (_detail_row(dgid, v, pc, cities, sel, n_cols, clin=clinical, parent=pango_parent)
                    if any((x or {}).get("timeline") for x in pc.values()) else "")
            ev = _det_link(dgid) if drow else ""
        elif kind == "finding" and v in findings:
            f = findings[v]
            # new / covered / suggested, from all cities (2026-10-20); a
            # suggestion its own mutations don't support where the reads
            # suggested it is folded into one line under the tree
            fkind, ov, ovtip, col, ncol, dk = _finding_overview(
                v, f, cities, current_panel, bool(r.get("has_children")))
            if fkind == "folded":
                # kept, but drawn flat at the bottom in a closed group that the
                # line under the tree opens (2026-10-21: PJ.2.1 had a signal
                # in the wastewater and couldn't be opened while hidden)
                folded.append(v)
                r = dict(r, depth=0, guides=[], last=True, has_children=False)
            fam_desc = f.get("descendants") or []
            tip = f"{v}: suggested by the reads, not in the panel"
            name = f"<span class='nm' style='color:{ncol}'>{_e(v)}</span>{chips}"
            tc = _tree_cell(r, col, dk, name, tip)
            tint = None                     # no row tint (2026-10-20)
            dgid = "det-" + re.sub(r"[^A-Za-z0-9]", "_", v)
            _chk = {c: p_["check"] for c, p_ in (f.get("per_city") or {}).items()
                    if (p_ or {}).get("check")}
            cells = _finding_cells(v, f, cities, sel)
            ev = ""
            family = fkind == "covered"
            its = (f.get("itself") or {}).get(sel)
            if its:
                # Jaccard with the lineages directly below it (found or in
                # the panel; abundance_cooc._family_by_jaccard, 2026-10-14)
                _parts = " · ".join(
                    f"{d} {'can’t tell' if j is None else f'{j:.2f}'}"
                    + (f" ({k})" if k else "")
                    for d, j, k in its[1][:5])
                _calc = (f"{v} itself = its share × (1 − Σ Jaccard) with the nearest lineages below it that "
                         f"are present and in the panel or offered + Add. Jaccard = their share "
                         f"÷ {v}’s share, each read at the ★ markers that every lineage below "
                         f"with reads here still carries: {_parts}")
                _it = ("can’t tell (fewer than 2 well-covered samples)"
                       if its[0] is None and not its[1]
                       else "can’t tell (too few reads, or no marker kept by all lineages below)"
                       if its[0] is None
                       else f"≈ {100 * its[0]:.1f} % of reads")
            # "itself" (how much of it is left after the lineages below it) goes
            # into the overview's hover; the details get it in part B
            if its and fkind == "covered":
                ovtip = (ovtip + " · " if ovtip else "") + (
                    f"{sname}: itself {_it} after the lineages below it — {_calc}")
            ov = f"<span data-tip='{_e(ovtip).replace(chr(10), '&#10;')}'>{ov}</span>" if ovtip else ov
            drow = (_detail_row(dgid, v, _chk, cities, sel, n_cols, found=True, clin=clinical,
                                parent=pango_parent)
                    if any(x.get("timeline") for x in _chk.values()) else "")
            ev = _det_link(dgid) if drow else ""
        else:  # tracked, not selected
            name = f"<span style='color:{TRACKED}'>{_e(v)}</span>{chips}"
            tc = _tree_cell(r, TRACKED, "ring", name, f"{v}: officially tracked, not in your panel")
            ov = _ovb("officially tracked · not in the panel, not in the run", "color:#8b8f97;")
            tint, cells, ev = None, f"<td colspan='{n_c}'></td>", ""
            drow, dgid = "", ""
        bg = f" style='background:{tint}{'22' if kind in ('panel', 'panel_ot') else '14'}'" if tint else ""
        nch = near.get(v) if kind in ("panel", "panel_ot") else None
        if nch:
            gid = "near-" + re.sub(r"[^A-Za-z0-9]", "_", v)
            # the tag follows the chosen city, like the Evidence column
            nch = sorted(nch, key=lambda c: -c["per_city"].get(sel, 0))
            here = [c for c in nch if c["per_city"].get(sel)]
            n_tag = sum(1 for c in here if c["sign"] == "+"
                        and c["per_city"][sel] >= near_min_days)
            if here:
                tip = (f"{v} ± 1 change in {sname} — reads that are {v} with one mutation "
                       f"more (+) or one less (−), in samples with ≥ {nn['er']} such reads and "
                       f"≥ {nn['es']} % of the sample's reads: "
                       + "; ".join(f"{c['label']} ({c['per_city'][sel]} d)"
                                   for c in here[:near_max])
                       + (f"; +{len(here) - near_max} more" if len(here) > near_max else ""))
            else:
                tip = (f"No ± 1 change in any sample in {sname}; in other cities: "
                       + "; ".join(f"{c['label']}" for c in nch[:near_max])
                       + (f"; +{len(nch) - near_max} more" if len(nch) > near_max else ""))
            label = (f"◆ {n_tag} {_pl(n_tag, 'change')}" if n_tag
                     else (f"◆ {len(here)} {_pl(len(here), 'change')} (1 day)" if nn["ed"] == 2
                           else f"◆ {len(here)} {_pl(len(here), 'change')} (< {nn['ed']} days)")
                     if here
                     else "◆ elsewhere")
            cls = "near-tag" + ("" if here else " off")
            tag = (f"<span class='{cls}' data-grp='{gid}' data-tip='{_e(tip)}'>"
                   f"{label} <span class='arr'>▸</span></span>")
            tc = tc.replace("</span></td>", f"</span>{tag}</td>", 1) if tc.endswith("</span></td>") else tc
        # data-gid: clicking a city cell of this row also opens its details
        gid_attr = f" data-gid='{dgid}'" if drow else ""
        if kind == "finding" and v in folded:
            # into the closed "not supported" group at the bottom
            folded_rows.append(f"<tr class='g-fold hid'{gid_attr}>{tc}<td class='ov'>{ov}</td>"
                               f"{cells}<td class='ev'>{ev}</td></tr>")
            if drow:
                folded_rows.append(drow)
            continue
        body.append(f"<tr{bg}{gid_attr}>{tc}<td class='ov'>{ov}</td>{cells}"
                    f"<td class='ev'>{ev}</td></tr>")
        if drow:
            body.append(drow)
        if nch:
            body += _near_rows(r, v, nch[:near_max], cities, sel, gid, n_cols)
            if len(nch) > near_max:
                body.append(f"<tr class='near note g-{gid} hid'><td colspan='{n_cols}'>"
                            f"showing {near_max} of {len(nch)} changes of {_e(v)}</td></tr>")

    if folded:
        body.append(
            f"<tr class='grp' data-grp='fold'><td colspan='{n_cols}'><span data-tip='Names the deep "
            "scan suggested that no city backs with coocs of their own mutations (or only one "
            "sample, with half or more of its mutations absent). Not counted; if one is real, the next "
            "batch moves it into the tree. Click to show their rows.'>"
            f"<span class='arr'>▸</span>{len(folded)} more {_pl(len(folded), 'name')} the deep "
            f"scan suggested, without coocs of {'its' if len(folded) == 1 else 'their'} own: "
            f"{_e(', '.join(folded))}</span></td></tr>")
        body += folded_rows

    if named_hidden:
        n_named = sum(1 for f in findings.values() if f.get("status") == "named")
        body.append(f"<tr class='note'><td colspan='{n_cols}'>"
                    + (f"Not confirmed lineages: showing {n_named} of {n_named + named_hidden}"
                       if n_named else
                       f"{named_hidden} not confirmed {_pl(named_hidden, 'lineage')} hidden — the "
                       "deep scan named them, but their reads also fit related lineages "
                       "(lists.named_only_max)")
                    + "</td></tr>")

    # ── unnamed signal ──
    ed = nn["ed"]
    # recurring (2026-10-27): evidence days on well-covered samples in >= 2
    # different weeks, any city; days on thin samples are shown, not counted
    _w = {tuple(m): _novel_weight(novel_info, sample_well, m) for m, d in novel}
    if sample_well is not None:
        rec = [(m, d) for m, d in novel if _w[tuple(m)][0] >= 2]
        one = [(m, d) for m, d in novel if _w[tuple(m)][0] < 2 and max(d.values(), default=0) > 0]
        rec.sort(key=lambda x: (-_w[tuple(x[0])][1], -_w[tuple(x[0])][0]))
    else:                     # no per-sample coverage: as before
        rec = [(m, d) for m, d in novel if max(d.values(), default=0) >= ed]
        one = [(m, d) for m, d in novel if 0 < max(d.values(), default=0) < ed]
    low = [(m, d) for m, d in novel if max(d.values(), default=0) == 0]
    if novel or broad or novel_rest or novel_hotspot:
        body.append(f"<tr class='sec'><td colspan='{n_cols}'>Unassigned patterns</td></tr>")
        novel_info = novel_info or {}

        def nov_ev_for(m, d):
            info = (novel_info.get(tuple(m)) or {}).get(sel)
            nw, nc = _w.get(tuple(m), (0, 0))
            where = (f"<span class='dim'> · well-covered evidence in {nc} {_pl(nc, 'city')}, "
                     f"{nw} {_pl(nw, 'week')}</span>" if sample_well is not None and nw else "")
            if not info:
                n_other = sum(1 for c, x in d.items() if c != sel)
                return (f"<span class='dim'>not in this city"
                        + (f" · in {n_other} other {_pl(n_other, 'city')}" if n_other else "")
                        + "</span>" + where)
            return _novel_evidence(m, info) + where
        _rlab = ("· on well-covered samples in ≥ 2 weeks · watch the next samples"
                 if sample_well is not None else f"· on ≥ {ed} days in a city · watch the next samples")
        body += _group("nrec", f"<b style='color:{BLUE}'>Unassigned, recurring ({len(rec)})</b> "
                       f"<span class='dim'>{_rlab}</span>", n_cols,
                       [_pattern_row(m, d, cities, sel, S_NOVEL, S_NOVEL1, BLUE, nov_ev_for(m, d), "seen")
                        for m, d in rec], open_=True)
        one_lab = ("seen, not recurring" if sample_well is not None
                   else "1 day" if ed == 2 else f"1–{ed - 1} days")
        body += _group("none", f"<b style='color:{BLUE}'>Unassigned, {one_lab} ({len(one)})</b> "
                       + ("<span class='dim'>· evidence in 1 week only, or only on thin samples</span>"
                          if sample_well is not None else "<span class='dim'>· watch</span>"), n_cols,
                       [_pattern_row(m, d, cities, sel, S_NOVEL, S_NOVEL1, BLUE, nov_ev_for(m, d), "seen")
                        for m, d in one])
        body += _group("nlow", f"<b style='color:{BLUE}'>Unassigned, no evidence day ({len(low)})</b> "
                       f"<span class='dim'>· never {nn['er']} reads and {nn['es']} % of a day</span>",
                       n_cols,
                       [_pattern_row(m, d, cities, sel, S_NOVEL, S_NOVEL1, BLUE, nov_ev_for(m, d), "seen")
                        for m, d in low])
        n_broad = (f"showing {len(broad)} of {broad_total}"
                   if broad_total and broad_total > len(broad) else f"{len(broad)}")
        body += _group("broad", f"<b style='color:{SLATE}'>Fits too many lineages ({n_broad})</b> "
                       "<span class='dim'>· the combination fits many unrelated lineages</span>",
                       n_cols,
                       [_pattern_row(m, d, cities, sel, S_BROAD, _dashed(SLATE), SLATE,
                                     f"<span class='dim'>fits {n} lineages"
                                     + (f" · mostly {_e(a)}" if a else "") + "</span>", "seen")
                        for m, d, n, a in broad])
        if novel_hotspot:
            _tip = ("A position where the city's data shows two or more different new bases "
                    f"(each on ≥ {nn['af']} % of ≥ {nn['mc']} reads) is an error hotspot — a virus has one. "
                    "A mutation there that no pango lineage carries is left out of unassigned "
                    "patterns; a pattern left with one mutation, or one a lineage carries, "
                    "is not unassigned. Lineage mutations are always kept.")
            body.append(f"<tr class='note'><td colspan='{n_cols}'><span data-tip='{_e(_tip)}'>"
                        f"{novel_hotspot[0]:,} patterns ({novel_hotspot[1]} reads) at error-hotspot "
                        "positions — left out ⓘ</span></td></tr>")
        if novel_rest:
            body.append(f"<tr class='note'><td colspan='{n_cols}'>+ more small unassigned patterns "
                        f"({novel_rest} reads, summed over cities and days) — not listed</td></tr>")

    return (CSS + SHORT_LEGEND + _frames_note(findings, panel_weekly, dates) + _add_note(findings)
            + f"<div class='vtw'><table class='vt'><thead>{head}</thead>"
              f"<tbody>{''.join(body)}</tbody></table></div>")


def _novel_weight(novel_info, sample_well, muts):
    """(weeks, cities) with an evidence day on a well-covered sample"""
    from datetime import date as _dd, timedelta as _tt
    wk, cs = set(), set()
    for c, inf in ((novel_info or {}).get(tuple(muts)) or {}).items():
        well = (sample_well or {}).get(c) or {}
        for day, mark in (inf or {}).get("timeline") or []:
            day = str(day)[:10]
            if mark == "day" and well.get(day):
                x = _dd.fromisoformat(day)
                wk.add((x - _tt(days=x.weekday())).isoformat())
                cs.add(c)
    return len(wk), len(cs)


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


def _lc(d, found=False, extra=""):
    st, lb = _scale(d, found=found)[:2]
    return f"<i class='c' style='width:auto;padding:0 6px;{st}{extra}'>{_lab_html(lb)}</i>"


_LEGEND_ROWS = (
    ""
    "<div class='row'><span class='t'>★ check</span>"
    "<span class='i'>each cell is the vote of the lineage's ★ markers, like a score: "
    + f"<b>present : absent</b> (present ≥ {_num()['pf']} % of ≥ {_num()['mc']} reads, absent "
    + f"< {_num()['af']} %), and in brackets the markers that couldn't answer (under "
    + f"{_num()['mc']} reads after weighting by coverage, between {_num()['af']} and "
    + f"{_num()['pf']} %, or reads unlike it)</span>"
    + "<span class='i'>"
    + _lc({"n_markers": 4, "n_measured": 4, "n_present": 3}) + "present: more present than absent, "
    + f"and at least {_min_present()} present</span>"
    + "<span class='i'>" + _lc({"n_markers": 4, "n_measured": 4, "n_present": 0}) + "absent: more absent than present</span>"
    + "<span class='i'>" + _lc({"n_markers": 2, "n_measured": 2, "n_present": 1})
    + _lc({"n_markers": 2, "n_measured": 1, "n_present": 1})
    + f"undecided: a tie, or fewer than {_min_present()} present</span>"
    + "<span class='i'><b style='color:" + AMBER_T + "'>in your panel</b> · present in none of … cities = "
      "a lineage of your panel no city calls present</span>"
    + "<span class='i'>a row with no ★ marker says so once, across the cities</span>"
    + "<span class='i'>details ▸ = samples per city and the markers' reads</span>"
    + "<span class='i'>in hovers: ✓ present (≥ " + f"{_num()['pf']}" + " %) · ✗ absent (< "
    + f"{_num()['af']}" + " %) · ~ in between · ? too few reads (or too few linked reads) · ≠ on reads unlike it · "
      "▧ low-coverage sample</span>"
    + "</div><div class='row'><span class='t'>Found</span>"
    + _li(S_NAMED, "0", "not confirmed (reads also fit related lineages)")
    + "<span class='i'><b style='color:#9ca3af'>covered</b>: its own mutations are on the reads "
      "of lineages below it (itself under the present limit), no + Add — add those instead</span>"
    + _li(S_NOVEL, "3", "unassigned · days")
    + _li(S_BROAD, "2", "fits too many lineages · days")
    + _li(S_NEAR, "3", "◆ panel variant ± 1 change · days (a hint)")
    + "</div><div class='row'><span class='t'>Tree</span>"
    + _ld(f"background:{GREEN}", "in panel, seen somewhere")
    + _ld(f"background:{AMBER}", "in panel, not seen in any city")
    + _ld(f"background:{RED}", "add to panel: confirmed 3 weeks in a row (the 4 weeks up to "
          "each Sunday: all cities pooled seen, ≥ 2 well-covered cities, ≥ 2 weeks, ≤ 10 % of "
          "its mutations missing) · + Add")
    + _ld(f"border:2px solid {RED}", "watch: found, not ready to add; the row says why")
    + _ld("border:2px solid #9ca3af", "suggested: a name from the deep scan without coocs of "
          "its own (most are listed in one closed line under the tree)")
    + _ld("border:2px solid #d1d5db", "covered: its own mutations are on a lineage below it")
    + _ld(f"border:2px solid {TRACKED}", "OT, not selected")
    + "<span class='i'>the dots and the overview column cover all cities and don't change when "
      "you click one · click a city's code or any of its cells to see it in the details · "
      "suggestions their own mutations don't support are listed in one line under the "
      "tree</span></div>")

# The full legend now lives behind "ⓘ How to read this page" above the tabs
# (components.how_to_read, 2026-10-13); the table keeps one short line on top,
# always in view, with the cell colours you need while reading the tree.
SHORT_LEGEND = (
    "<div class='lg lgs'><div class='row'>"
    + "<span class='i'><b>cells</b>: samples with 2 or more separate coocs (its own mutations "
      "together on the same reads) only its family has / samples read</span>"
    + f"<span class='i'><i class='c' style='width:auto;padding:0 6px;{S_YES}'>2/8</i>seen</span>"
    + f"<span class='i'><i class='c' style='width:auto;padding:0 6px;{S_YES}{_LOWQ_CELL}'>1/8</i>"
      "seen, only on thin samples</span>"
    + f"<span class='i'><i class='c' style='width:auto;padding:0 6px;{S_SINGLE}'>0/8</i>"
      "only single mutations seen (weaker)</span>"
    + f"<span class='i'><i class='c' style='width:auto;padding:0 6px;{S_NO}'>0/8</i>not seen</span>"
    + f"<span class='i'><i class='c' style='width:auto;padding:0 6px;{S_OPEN}'>0/4</i>undecided</span>"
    + "<span class='i dim'>results without pair counts: " + _lc({"n_markers": 3, "n_measured": 3, "n_present": 3})
    + "seen : not seen (couldn't answer)</span>"
    + "<span class='i'><span class='cal'><i style='background:" + GREEN + ";" + _LOWQ
    + "'></i></span>low coverage</span>"
    + _ld(f"background:{GREEN}", "in panel")
    + _ld(f"background:{AMBER}", "in panel, not seen in any city")
    + _ld(f"background:{RED}", "add to panel")
    + _ld(f"border:2px solid {RED}", "watch")
    + _ld("border:2px solid #d1d5db", "covered / suggested")
    + "<span class='dim'>click a cell for its details · ⓘ more: “How to read this page” above the tabs</span>"
    + "</div></div>")


def legend_full_html() -> str:
    """Every colour and symbol of the table, for "How to read this page"."""
    return CSS + "<div class='lg'>" + _LEGEND_ROWS + "</div>"


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
       background:#31333f; color:#fff; font-size:12px; pointer-events:none; white-space:pre-line; }
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
            body.append(f"<tr>{_tree_cell(r, TRACKED, 'ring', name, f'{v}: OT, not selected')}</tr>")
            h += 28
    legend = ("<div class='lg'><span class='i'><i class='d' style='background:#1f2430'></i>in your panel</span> · "
              f"<span class='i'><i class='d' style='border:2px solid {TRACKED}'></i>OT, not selected</span> · "
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
  // any cell of a city column picks that city, so there's no need to
  // scroll back to the header codes (2026-10-08)
  // a cell also opens its row's details for that city (2026-10-20): the
  // details stay open across the redraw that the city change triggers
  var td = e.target.closest("td.cc");
  if (td && td.parentNode) {
    var tr = td.parentNode, gid = tr.getAttribute("data-gid");
    var own = Array.prototype.slice.call(tr.querySelectorAll("td.cc"));
    var heads = document.querySelectorAll("th.cc [data-city]");
    var i = own.indexOf(td);
    if (gid) {
      window.vtOpen[gid] = true;
      document.querySelectorAll(".g-" + gid).forEach(function (r) { r.classList.remove("hid"); });
      document.querySelectorAll("[data-grp='" + gid + "'] .arr").forEach(function (a) {
        a.textContent = "▾";
      });
      sendHeight();
    }
    if (own.length === heads.length && i >= 0 && !td.classList.contains("sel")) {
      send("streamlit:setComponentValue",
           {value: {city: heads[i].getAttribute("data-city"), t: Date.now()}, dataType: "json"});
    }
    return;
  }
  var g = e.target.closest("[data-grp]");
  if (g) {
    var id = g.getAttribute("data-grp"), open = false;
    document.querySelectorAll(".g-" + id).forEach(function (r) {
      open = r.classList.toggle("hid") === false;
    });
    document.querySelectorAll("[data-grp='" + id + "'] .arr").forEach(function (a) {
      a.textContent = open ? "▾" : "▸";
    });
    window.vtOpen[id] = open;
    sendHeight();
  }
};

// open details / groups stay open when the table is redrawn (a city click
// sends new html, which replaces the whole table) — 2026-10-13
window.vtOpen = window.vtOpen || {};
window.vtReopen = function () {
  Object.keys(window.vtOpen).forEach(function (id) {
    var open = window.vtOpen[id];
    document.querySelectorAll(".g-" + id).forEach(function (r) {
      r.classList.toggle("hid", !open);
    });
    document.querySelectorAll("[data-grp='" + id + "'] .arr").forEach(function (a) {
      a.textContent = open ? "▾" : "▸";
    });
  });
  if (typeof window.sendHeight === "function") window.sendHeight();
};
if (!window.vtObs) {
  var root = document.getElementById("root") || document.body;
  window.vtObs = new MutationObserver(function (muts) {
    if (muts.some(function (m) { return m.target === root; })) window.vtReopen();
  });
  window.vtObs.observe(root, {childList: true});
}
window.vtReopen();
"""


def render(view_html: str, key: str):
    """Show the view; returns the last click, {"add": node, "t": ms} or
    {"city": name, "t": ms} (it stays the value on later reruns — compare t)."""
    global _COMPONENT
    if _COMPONENT is None:
        import streamlit.components.v1 as components
        _COMPONENT = components.declare_component("variants_table", path=str(_FRONTEND))
    return _COMPONENT(html=view_html, js=CLICK_JS, key=key, default=None)