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
.vt { border-collapse:collapse; font-size:15px; width:100%; }
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
.vt .ovb b { font-weight:600; }
.vt .ovb .add { margin-left:6px; }
.vt td.cc .c { cursor:pointer; }
.vt td.cc:not(.sel):hover { background:rgba(49,51,63,.04); }
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
.vt .dl b { min-width:84px; font-weight:600; }
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
    """'0:24 (15)' with the bracket smaller and lighter: the score stands out."""
    return re.sub(r" \((\d+)\)$", r"<span class='cu'>(\1)</span>", _e(label))


def _cell(label, style, tip, sel) -> str:
    return (f"<td class='cc{' sel' if sel else ''}'><span class='c' style='{style}' "
            f"data-tip='{_e(tip)}'>{_lab_html(label)}</span></td>")


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


def _scale(d, found=False):
    """(cell style, label, dot colour, name colour, share present or None)
    for one city's check result. found: a lineage not in the panel.
    The cell (style, label) shows the answer only (2026-10-20); the dot and
    name colours of the tree are unchanged."""
    n_mk, p, a, u = _vote_counts(d)
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
        if weak or _few_clear(d):
            return (S_YES, lab, ("#fca5a5" if found else "#86efac"),
                    (RED if found else "#15803d"), share)
        return S_YES, lab, col, (RED if found else "#15803d"), share
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


def _cell_tip(v, city, d, extra=""):
    """'XFG · Basel — mixed (1 ✓ 1 ✗ of 2 ★) · 2 well-covered samples · 4214C
    50 % ✓ · 8350C 0 % ✗' (+ why it was checked)."""
    n_mk, p, a, u = _vote_counts(d)
    if n_mk == 0:
        return f"{v} · {city} — no ★ marker"
    word = _WORD.get(_state(d), "too few reads")
    counts = f"{p} ✓ {a} ✗" + (f" {u} ?" if u else "") + f" of {n_mk} ★"
    samples, basis = _basis_samples(d)
    base = (f"{len(samples)} {_pl(len(samples), 'sample')} weighted by coverage (each ≤ "
            f"{_cov_factor():g} × its median depth, ≤ {_cap_reads():,} reads)" if basis == "weighted"
            else f"{len(samples)} well-covered {_pl(len(samples), 'sample')}" if basis == "well"
            else "weak: no coverage data (all counted, ≤ "
                 f"{_cap_reads():,} reads each)")
    if _state(d) == "few":
        base += f" · only {p} ★ measured, {_min_present()} needed"
    ms = _marker_syms(d)
    # what the light shades said until 2026-10-20 is in words: "weak: no
    # coverage data" in base, and on how few samples it was clear
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


def _detail_row(gid, v, per_city, cities, sel, n_cols, found=False):
    """The row's details (hidden until 'details ▸'): one calendar line per
    city, the markers with their real reads in the selected city, and their
    median share across cities."""
    top, low = (RED, None) if found else (GREEN, AMBER)
    lines = []
    for c in cities:
        d = per_city.get(c)
        if not d or not d.get("timeline"):
            continue
        lines.append(f"<div class='dl'><b>{_e(city_name(c))}</b>{_cal_line(d, top, low)}"
                     f"<span class='dim'>{_e(_summary(d))}</span></div>")
    ds = per_city.get(sel) or {}
    mk_rows = []
    import statistics
    for x in ds.get("detail") or []:
        m = x.get("marker")
        shares = [y.get("freq") for dd in per_city.values() for y in (dd.get("detail") or [])
                  if y.get("marker") == m and y.get("freq") is not None]
        med = f"{100 * statistics.median(shares):.0f} %" if shares else "–"
        reads = x.get("cov_raw", x.get("cov", 0)) or 0
        share = f"{100 * x['freq']:.1f} % of {reads:,} reads" if x.get("freq") is not None and reads else "no reads"
        mk_rows.append(f"<tr><td class='mono'>{_ext(_clin_link([m]), m, 'which lineages carry it in clinical sequences (CoV-Spectrum)')}</td><td>{_e(share)}</td>"
                       f"<td>{_e(x.get('why') or '')}</td><td>{med}</td></tr>")
    mk = ""
    if mk_rows:
        basis = {"well": "well-covered samples",
                 "weighted": "all samples, weighted by coverage"}.get(ds.get("basis"), "all samples")
        mk = (f"<table class='mk'><tr><th>★ marker</th><th>{_e(city_name(sel))}, {basis}</th>"
              "<th></th><th>median share, all cities</th></tr>" + "".join(mk_rows) + "</table>")
    body = ("<div class='dim' style='font-size:12px'>Samples per city, newest right · striped = "
            "low coverage (counts less in the vote) · ⋯ = older "
            "low-coverage samples without it, not shown</div>" + "".join(lines) + mk
            + "<div style='font-size:12.5px;margin-top:6px'>"
            + _ext(_ww_link(v, sel), f"{v} in {city_name(sel)} wastewater (GenSpectrum)",
                   "mutation heatmaps of the lineage and its sublineages")
            + ((" · " + _ext(_clin_link([x['marker'] for x in ds.get('detail') or []]),
                             "its ★ markers together in clinical data (CoV-Spectrum)",
                             "sequences carrying all of them: which lineages"))
               if ds.get("detail") else "")
            + "<span class='dim'> · each mutation per week: Signal over time (Look-up)</span></div>")
    return f"<tr class='dt g-{gid} hid'><td colspan='{n_cols}'>{body}</td></tr>"


def _panel_cells(v, per_city, cities, sel):
    res = [d for d in per_city.values() if d]
    if res and all(_vote_counts(d)[0] == 0 for d in res):
        # no ★ marker: the same in every city, so said once (2026-10-20)
        return (f"<td class='nomk' colspan='{len(cities)}'>no mutation of its own — all are "
                "shared with other lineages on your list</td>")
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


def _panel_overview(per_city, cities):
    res = [per_city[c] for c in cities if per_city.get(c)]
    if not res:
        return (_ovb("<b>On your list</b> · no result"), "", GREY, "#6b7280", "dashed")
    if all(_vote_counts(x)[0] == 0 for x in res):
        return (_ovb("<b>On your list</b> · can't be checked", "color:#6b7280;"),
                "no mutation of its own: all its mutations are shared with other lineages on "
                "your list, so its share in Deconvolution can't be checked here",
                GREY, "#6b7280", "filled")
    if _panel_warn(per_city, True):
        return (_ovb("<b>On your list</b> · not seen in any city",
                     f"background:#fffbeb;color:{AMBER_T};"),
                _panel_warn(per_city, True).replace("<b style='color:" + AMBER_T + "'>", "")
                .replace("</b>", "").rstrip(" ·"),
                AMBER, AMBER_T, "filled")
    meas = [c for c in cities if per_city.get(c) and _vote_counts(per_city[c])[0]]
    seen = [c for c in meas if _state(per_city[c]) == "present"]
    where = f" ({', '.join(city_code(c) for c in seen)})" if 0 < len(seen) <= 3 else ""
    return (_ovb(f"<b style='color:#1f2328'>On your list</b> · seen in {len(seen)} of "
                 f"{len(meas)}{where}"),
            "", GREEN if seen else "#d1d5db", "#15803d" if seen else "#6b7280", "filled")


def _finding_status(v, f, cities, has_children=False):
    """'new' / 'possible' / 'covered' / 'named' / 'folded' for a lineage the
    reads suggested (a deep scan find), and the reason for 'possible'.
    folded: its own mutations say not seen in every city where the reads
    suggested it, and no city says seen or unclear (a tie) — a leaf only, so
    no branch loses its children (2026-10-20)."""
    st = f.get("status")
    if st == "family":
        return "covered", ""
    if st == "named":
        return "named", ""
    pc = f.get("per_city") or {}
    chk = {c: (pc.get(c) or {}).get("check") for c in cities if (pc.get(c) or {}).get("check")}
    state = {c: _state(d) for c, d in chk.items()}
    seen = [c for c in cities if state.get(c) == "present"]
    tie = [c for c in cities if state.get(c) == "mixed"]
    few = [c for c in cities if state.get(c) == "few"]
    sugg = {c: int((pc.get(c) or {}).get("days") or 0) for c in cities
            if not (pc.get(c) or {}).get("xcheck") and int((pc.get(c) or {}).get("days") or 0) > 0}
    if f.get("addable"):
        return "new", ""
    if (not seen and not tie and sugg and not has_children
            and all(state.get(c) == "absent" for c in sugg)):
        return "folded", ""
    ed = _num()["ed"]
    one_sample = sugg and max(sugg.values()) < ed and len(sugg) < 2
    code = lambda cs: ", ".join(city_code(c) for c in cs)                    # noqa: E731
    top = max(sugg, key=sugg.get) if sugg else None
    sug_txt = (f"suggested {sugg[top]}× in {city_code(top)}" if top and not one_sample else
               "suggested on 1 sample only" if one_sample else "")
    if seen:
        why = f"seen in {code(seen)}"
    elif tie:
        why = f"unclear in {code(tie)}"
    elif few:
        p = _vote_counts(chk[few[0]])[1]
        why = f"{p} of its own mutations seen in {city_code(few[0])} (needs {_min_present()})"
    elif top and state.get(top) not in ("absent",):
        why = f"couldn't check in {city_code(top)}"
    else:
        why = "its own mutations not seen"
    return "possible", " · ".join(x for x in (why, sug_txt) if x)


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
        seen = [c for c in cities if _state(((pc.get(c) or {}).get("check")) or {}) == "present"]
        btn = ("<span class='added'>added · re-run needed</span>" if v in current_panel
               else f"<button class='add' data-add='{_e(v)}'>+ Add</button>")
        return (kind, _ovb(f"<b>New</b> · seen in {', '.join(city_code(c) for c in seen)}{btn}",
                           "background:#fef2f2;color:#b91c1c;"),
                "suggested by the reads and seen through its own mutations; not on your list",
                RED, RED, "filled")
    if kind == "folded":
        return (kind, "", "", GREY, "#9ca3af", "ring")
    return (kind, _ovb(f"<b style='color:#1f2328'>Possible</b> · {_e(why)}",
                       "color:#4b5563;box-shadow:inset 0 0 0 1px #e5e7eb;"),
            "suggested by the reads, not on your list; its own mutations only partly confirm it",
            "#9ca3af", "#4b5563", "ring")


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
          broad_total=None, named_hidden=0, quality=None) -> str:
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
    folded = []
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
            # selected one (2026-10-20); amber only for "on your list, not
            # seen in any city"
            ov, ovtip, dcol, ncol, dk = _panel_overview(pc, cities)
            tip = f"{v}: on your list" + (f" · {ovtip}" if ovtip else "")
            name = f"<span class='nm' style='color:{ncol}'>{_e(v)}</span>{chips}"
            tc = _tree_cell(r, dcol, dk, name, tip)
            tint = None                     # no row tint (2026-10-20)
            dgid = "det-" + re.sub(r"[^A-Za-z0-9]", "_", v)
            cells = _panel_cells(v, pc, cities, sel)
            ov = f"<span data-tip='{_e(ovtip)}'>{ov}</span>" if ovtip else ov
            drow = (_detail_row(dgid, v, pc, cities, sel, n_cols)
                    if any((x or {}).get("timeline") for x in pc.values()) else "")
            ev = _det_link(dgid) if drow else ""
        elif kind == "finding" and v in findings:
            f = findings[v]
            # new / possible / covered, from all cities (2026-10-20); a
            # suggestion its own mutations don't support where the reads
            # suggested it is folded into one line under the tree
            fkind, ov, ovtip, col, ncol, dk = _finding_overview(
                v, f, cities, current_panel, bool(r.get("has_children")))
            if fkind == "folded":
                folded.append(v)
                continue
            fam_desc = f.get("descendants") or []
            tip = f"{v}: suggested by the reads, not on your list" + (f" · {ovtip}" if ovtip else "")
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
                _calc = (f"{v} itself = 1 − Σ Jaccard with the nearest lineages below it that "
                         f"are present and in the panel or offered + Add. Jaccard = their share "
                         f"÷ {v}’s share, each read at the ★ markers that every lineage below "
                         f"with reads here still carries: {_parts}")
                _it = ("can’t tell (too few reads, or no marker kept by all lineages below)"
                       if its[0] is None
                       else f"≈ {100 * its[0]:.0f} %")
            # "itself" (how much of it is left after the lineages below it) goes
            # into the overview's hover; the details get it in part B
            if its:
                ovtip = (ovtip + " · " if ovtip else "") + (
                    f"{sname}: itself {_it} after the lineages below it — {_calc}")
            ov = f"<span data-tip='{_e(ovtip)}'>{ov}</span>" if ovtip else ov
            drow = (_detail_row(dgid, v, _chk, cities, sel, n_cols, found=True)
                    if any(x.get("timeline") for x in _chk.values()) else "")
            ev = _det_link(dgid) if drow else ""
        else:  # tracked, not selected
            name = f"<span style='color:{TRACKED}'>{_e(v)}</span>{chips}"
            tc = _tree_cell(r, TRACKED, "ring", name, f"{v}: officially tracked, not in your panel")
            ov = _ovb("officially tracked · not on your list, not in the run", "color:#8b8f97;")
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
            f"<tr class='note'><td colspan='{n_cols}'><span data-tip='Suggested by the reads in "
            "some city, but their own mutations say not seen there, and no city sees them or is "
            "unclear. Leaves of the tree only.'>"
            f"▸ {len(folded)} {_pl(len(folded), 'suggestion')} not supported by "
            f"{'its' if len(folded) == 1 else 'their'} own mutations where the reads suggested "
            f"{'it' if len(folded) == 1 else 'them'}: {_e(', '.join(folded))}</span></td></tr>")

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
    rec = [(m, d) for m, d in novel if max(d.values(), default=0) >= ed]
    one = [(m, d) for m, d in novel if 0 < max(d.values(), default=0) < ed]
    low = [(m, d) for m, d in novel if max(d.values(), default=0) == 0]
    if novel or broad or novel_rest or novel_hotspot:
        body.append(f"<tr class='sec'><td colspan='{n_cols}'>Unassigned patterns</td></tr>")
        novel_info = novel_info or {}

        def nov_ev_for(m, d):
            info = (novel_info.get(tuple(m)) or {}).get(sel)
            if not info:
                n_other = sum(1 for c, x in d.items() if c != sel)
                return (f"<span class='dim'>not in this city"
                        + (f" · in {n_other} other {_pl(n_other, 'city')}" if n_other else "")
                        + "</span>")
            return _novel_evidence(m, info)
        body += _group("nrec", f"<b style='color:{BLUE}'>Unassigned, recurring ({len(rec)})</b> "
                       f"<span class='dim'>· on ≥ {ed} days in a city — check the trend: possible new variant</span>", n_cols,
                       [_pattern_row(m, d, cities, sel, S_NOVEL, S_NOVEL1, BLUE, nov_ev_for(m, d), "seen")
                        for m, d in rec], open_=True)
        one_lab = "1 day" if ed == 2 else f"1–{ed - 1} days"
        body += _group("none", f"<b style='color:{BLUE}'>Unassigned, {one_lab} ({len(one)})</b> "
                       "<span class='dim'>· watch</span>", n_cols,
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

    return (CSS + SHORT_LEGEND
            + f"<table class='vt'><thead>{head}</thead><tbody>{''.join(body)}</tbody></table>")


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
    + _ld(f"background:{GREEN}", "on your list, seen somewhere")
    + _ld(f"background:{AMBER}", "on your list, not seen in any city")
    + _ld(f"background:{RED}", "new: suggested by the reads and seen (+ Add)")
    + _ld("border:2px solid #9ca3af", "possible: suggested, only partly seen")
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
    + "<span class='i'>cells: its own mutations seen : not seen (couldn't answer)</span>"
    + "<span class='i'>" + _lc({"n_markers": 3, "n_measured": 3, "n_present": 3}) + "seen</span>"
    + "<span class='i'>" + _lc({"n_markers": 2, "n_measured": 2, "n_present": 0}) + "not seen</span>"
    + "<span class='i'>" + _lc({"n_markers": 2, "n_measured": 1, "n_present": 1})
    + "undecided</span>"
    + "<span class='i'><span class='cal'><i style='background:" + GREEN + ";" + _LOWQ
    + "'></i></span>low coverage</span>"
    + _ld(f"background:{GREEN}", "on your list")
    + _ld(f"background:{AMBER}", "on your list, not seen in any city")
    + _ld(f"background:{RED}", "new: suggested by the reads")
    + _ld("border:2px solid #9ca3af", "possible / covered")
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