"""What the Lineages table decides per lineage, without Streamlit
(2026-10-15, moved from subpages/abundance_cooc.py unchanged), so the page and
the reproducibility study (bench_repro.py) use the same code:

  panel_verdicts   the ★ check of each panel lineage in one city
  plus_add         may a found lineage be offered with + Add
  family_rows      found ancestors carried by lineages below them -> "family"
"""
from __future__ import annotations


def panel_verdicts(cooc_res: dict, variants: list) -> dict:
    """Per panel variant in ONE city: {variant: {"state", "reason"}} from the
    worker's panel_check (reads only, never the abundance).

    state: present / absent / mixed (the marker vote of
    process.cooc.check_verdicts), no_marker (nothing specific to test),
    no_data (no check data, or no marker measurable)."""
    from process.cooc import check_verdicts
    from process.variant_explorer import check_in_data
    pc = (cooc_res or {}).get("panel_check", {}) or {}
    _dates = [str(d)[:10] for d in ((cooc_res or {}).get("dates") or [])]
    sym = {"present": "✓", "absent": "✗", "unmeasured": "?"}
    out = {}
    for v in variants:
        ci = pc.get(v)
        if ci is None:
            out[v] = {"state": "no_data", "reason": "no check data — run the panel check again",
                      "n_present": 0, "n_measured": 0, "n_markers": 0}
            continue
        # one vote per city (process.variant_explorer.check_in_data): every
        # sample in the run's dates votes, each counting at most check.cap_reads
        # reads and at most check.cov_factor x its median depth (2026-10-16;
        # until then the well-covered samples decided alone)
        _cid = check_in_data(ci.get("markers") or [], ci.get("per_date") or {}, _dates,
                             quality=(cooc_res or {}).get("sample_qc"))
        nd, np_, nm = _cid["n_markers"], _cid["n_present"], _cid["n_measured"]
        _sym = {"present": "✓", "absent": "✗"}
        mk = [f"{x['marker']} {(x['freq'] or 0) * 100:.0f}% {_sym.get(x['status'], '?')}"
              if x.get("cov") else f"{x['marker']} no reads ?" for x in (_cid.get("detail") or [])]
        mtxt = " · ".join(mk[:6]) + (f" · +{len(mk) - 6} more" if len(mk) > 6 else "")
        if nd == 0:
            state, reason = "no_marker", ("every mutation is shared with another panel "
                                          "variant or common outside its family")
        elif nm == 0:
            state, reason = "no_data", f"too few reads on its markers — {mtxt}"
        else:
            state = _cid["state"] if _cid["state"] in ("present", "absent", "mixed") else "no_data"
            reason = f"{np_} of {nm} covered markers present — {mtxt}"
        out[v] = {"state": state, "reason": reason,
                  "n_present": np_, "n_measured": nm, "n_markers": nd,
                  # one mark per sampling day, for the Evidence calendar
                  "timeline": _cid["timeline"],
                  # coverage quality of the samples (2026-10-13)
                  "quality": _cid.get("quality") or {}, "n_well": _cid.get("n_well"),
                  "n_samples": _cid.get("n_samples"), "basis": _cid.get("basis"),
                  # each marker's numbers, for the row's details (2026-10-13)
                  "detail": _cid.get("detail") or [],
                  "day_detail": _cid.get("day_detail") or {},
                  # counted reads per sample, for the all-cities week line (2026-10-22)
                  "markers": _cid.get("markers") or [],
                  "samples": _cid.get("samples") or {}}
    return out


def plus_add(per_city: dict, status: str, min_days: int) -> bool:
    """+ Add (2026-10-02): the ★ check confirms it in at least one city, and it
    is more than one sample (>= min_days evidence days in a city, or 1 day in
    >= 2 cities); never for a "named" (not confirmed) or "family" lineage.
    per_city: {city: {"days": int, "check": {"state": …}}}."""
    if status in ("named", "family"):
        return False
    present = any((d.get("check") or {}).get("state") == "present" for d in per_city.values())
    ok = any((d.get("days") or 0) >= min_days for d in per_city.values())
    samples = ok or sum(1 for d in per_city.values() if (d.get("days") or 0) >= 1) >= 2
    return bool(present and samples)


def family_rows(findings: dict, scan_res_all: dict, cr_all: dict,
                       panel: list) -> None:
    """A found lineage whose ★ signal comes from lineages below it
    (2026-10-14). Its ★ markers are carried by all its descendants, so its ★
    check measures the whole family (BA.2.86 "present" in Zürich because JN.1
    and everything after it carry its markers). Per city where X is present:

      1. the nearest lineages below X that will be in the panel: in it
         already, or found and offered + Add (present, not a family row);
         a mixed, absent or family lineage is skipped and the walk goes on
         below it (process.cooc.nearest_below)
      2. each one's share at its markers that every lineage below it with
         reads here (not absent) still carries, and X's share at its markers
         that those carry (process.cooc.kept_markers, lineage_share)
      3. Jaccard = share of D ÷ share of X;  left = 1 − Σ Jaccard;
         X itself = X's share × left: the share of reads that are X and
         not a lineage below it

    on the samples the ★ check uses (the well-covered ones when the city has
    any). "Can't tell" when a share has too few reads or no marker is left.
    X becomes a "family" row (no + Add) when X itself is under the marker
    present limit (check.present_freq) in every city where it is present and
    could be measured, and at least one city could. Lineages lower in the
    tree are decided first, so a family row below X is skipped in X's walk.

    2026-10-23: (a) X itself is a share of reads, compared with the 5 %
    limit; until then the left-over FRACTION (1 − Σ Jaccard) was compared
    with 5 %. (b) A city with fewer than 2 well-covered samples can't tell:
    its shares are jackpots. XFG.14.8.1 (5 Sep – 5 Oct), whose 3 markers TA.1
    all carries, stayed "New" next to TA.1: Basel (0 of 8 well covered) left
    56 % of 34 %, Genève (7 of 8) 16 % of 16 % = 2.6 % of reads."""
    from process.cooc import (_check_cfg, family_jaccard, kept_markers, lineage_share,
                              nearest_below)
    from process.sample_qc import well_covered
    pf = _check_cfg()["present_freq"]
    panel = list(panel or [])
    pv_cache: dict = {}

    def _pstate(loc, v):
        if loc not in pv_cache:
            pv_cache[loc] = panel_verdicts(cr_all.get(loc) or {}, panel) if panel else {}
        return (pv_cache[loc].get(v) or {}).get("state")

    def _seen_markers(loc, v):
        """v's markers its own check says seen in this city (2026-10-23)."""
        if v in panel:
            if loc not in pv_cache:
                pv_cache[loc] = panel_verdicts(cr_all.get(loc) or {}, panel) if panel else {}
            det = (pv_cache[loc].get(v) or {}).get("detail") or []
        else:
            det = (((((findings.get(v) or {}).get("per_city") or {}).get(loc) or {})
                    .get("check") or {}).get("detail") or [])
        return {x.get("marker") for x in det if x.get("status") == "present"}

    def _fstate(loc, v):
        return ((((findings.get(v) or {}).get("per_city") or {}).get(loc) or {})
                .get("check") or {}).get("state")

    # lower lineages first (fewer lineages below them, over all cities)
    nb = {}
    for loc in (scan_res_all or {}):
        for x, ds in (((scan_res_all.get(loc) or {}).get("findings_check") or {})
                      .get("below_all") or {}).items():
            nb[x] = max(nb.get(x, 0), len(ds))
    for x in sorted(findings, key=lambda v: nb.get(v, 0)):
        fx = findings[x]
        if fx.get("status") == "named":
            continue
        itself, measured, below = {}, [], set()
        for loc, px in (fx.get("per_city") or {}).items():
            chk = (px or {}).get("check") or {}
            if chk.get("state") != "present":
                continue
            if chk.get("n_well") is not None and chk["n_well"] < 2:
                itself[loc] = (None, [])          # thin city: can't tell
                continue
            fc = (scan_res_all.get(loc) or {}).get("findings_check") or {}
            ball = fc.get("below_all") or {}
            sigs = {k: set(v) for k, v in (fc.get("sigs") or {}).items()}
            pvar = fc.get("per_variant") or {}
            xd = pvar.get(x) or {}
            if not ball.get(x) or not xd.get("markers"):
                continue
            pc = (cr_all.get(loc) or {}).get("panel_check") or {}

            def take(d, loc=loc):
                if d in panel:
                    return _pstate(loc, d) == "present"
                f = findings.get(d) or {}
                return (_fstate(loc, d) == "present"
                        and f.get("status") not in ("family", "named"))

            def has_reads(e, loc=loc):
                stt = _pstate(loc, e) if e in panel else _fstate(loc, e)
                return stt != "absent"

            ds = nearest_below(x, ball, take)
            if not ds:
                continue
            q = chk.get("quality") or {}
            dates = ({d for d, v in q.items() if well_covered(v)}
                     if chk.get("basis") == "well" else None)
            # the same weights as the ★ check's vote (2026-10-16)
            depth = ({d: float((v or {}).get("median", 0) or 0) for d, v in q.items()}
                     if chk.get("basis") == "weighted" and q else None)
            parts = []
            for d in ds:
                below.add(d)
                dd = pvar.get(d) if d not in panel else pc.get(d)
                dd = dd or {}
                dm = kept_markers(dd.get("markers") or [], d, ball, sigs, has_reads)
                # D's share at the markers its own check sees here: one it
                # doesn't (TA.1's 14808C, 2.5 % where its others are ~20 %)
                # would pull D's share down and leave X a false "itself"
                _seen = _seen_markers(loc, d)
                dm = [m for m in dm if m in _seen] or dm
                xm = kept_markers(xd.get("markers") or [], d, ball, sigs, has_reads)
                sd = lineage_share(dm, dd.get("per_date") or {}, dates, depth=depth)
                sx = lineage_share(xm, xd.get("per_date") or {}, dates, depth=depth)
                j = family_jaccard(sx, sd)
                parts.append((d, j, None if j is None else
                              f"{100 * sd:.0f} % ÷ {100 * sx:.0f} %"))
            sx_all = lineage_share(xd.get("markers") or [], xd.get("per_date") or {},
                                   dates, depth=depth)
            val = (None if sx_all is None or any(j is None for _d, j, _k in parts)
                   else sx_all * max(0.0, 1.0 - sum(j for _d, j, _k in parts)))
            itself[loc] = (val, sorted(parts, key=lambda t: -(t[1] or 0)))
            measured.append(val)
        if itself:
            fx["itself"] = itself
        known = [v for v in measured if v is not None]
        if known and all(v < pf for v in known):
            fx["status"] = "family"
            fx["addable"] = False
            fx["descendants"] = sorted(below)
