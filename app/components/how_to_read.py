"""“ⓘ How to read this page”: one popover button above the result tabs
(2026-10-13). Two tabs:

  Colours               every colour and symbol of the Lineages table, the
                        calendar and the panel check graph
  Rules and thresholds  each rule with its value read from cooc_config.yaml,
                        why it is there and whether it was tested (rules:)

The Lineages table keeps a one-line legend on top (variants_table.SHORT_LEGEND)
so the colours stay in view while reading the tree."""
from __future__ import annotations

import streamlit as st


def _pct(x: float) -> str:
    return f"{x * 100:g} %"


def _notes() -> dict:
    """{rule id: {why, tested}} from the config's rules: section."""
    try:
        from utils.config import get_cooc_setting
        r = get_cooc_setting("rules", default={}) or {}
        return r if isinstance(r, dict) else {}
    except Exception:
        return {}


def rules_rows() -> list:
    """One row per rule: rule · value · why · tested? · config. Values are read
    from the config (process.cooc._check_cfg, process.sample_qc.quality_cfg,
    the deep scan's evidence.*), so the table follows any change there."""
    from process.cooc import _check_cfg
    c = _check_cfg()
    try:
        from process.sample_qc import quality_cfg
        q = quality_cfg()
    except Exception:
        q = {"min_median_depth": 20, "min_share_100x": 0.4}
    try:
        from utils.config import get_cooc_setting
        er = int(get_cooc_setting("evidence.min_reads", 20))
        es = float(get_cooc_setting("evidence.min_share", 0.005))
        ed = int(get_cooc_setting("evidence.min_days", 2))
    except Exception:
        er, es, ed = 20, 0.005, 2
    rows = [
        ("min_cov", "a ★ marker is covered", f"≥ {c['min_cov']} reads", "check.min_cov"),
        ("present_freq", "marker present", f"≥ {_pct(c['present_freq'])} of its reads",
         "check.present_freq"),
        ("link", "… and its reads look like the lineage",
         f"≥ {_pct(c['link_min'])} of ≥ {c['min_link']} reads that reach a neighbouring "
         "mutation; share alone when no read can reach one", "check.link_min, check.min_link"),
        ("absent_freq", "marker absent", f"< {_pct(c['absent_freq'])} of its reads",
         "check.absent_freq"),
        ("confirm_share", "cell present",
         f"≥ {_pct(c['confirm_share'])} of the covered markers, and ≥ {c['min_present']} markers",
         "check.confirm_share, check.min_present"),
        ("notfound_share", "cell absent", f"≤ {_pct(c['notfound_share'])} of the covered markers; "
         "in between: mixed", "check.notfound_share"),
        ("cap_reads", "reads a sample counts per marker", f"at most {int(c['cap_reads']):,}",
         "check.cap_reads"),
        ("quality", "well-covered sample",
         f"median depth ≥ {q['min_median_depth']} and ≥ {_pct(q['min_share_100x'])} of the genome "
         "at ≥ 100 reads; a city with any decides on those alone", "quality.*"),
        ("calendar", "calendar line", f"the newest {c.get('calendar_samples', 12)} samples, plus "
         "older ones where it was present or that are well covered", "check.calendar_samples"),
        ("evidence_day", "deep scan: evidence day",
         f"≥ {er} reads and ≥ {_pct(es)} of the sample", "evidence.min_reads, evidence.min_share"),
        ("min_days", "deep scan: found", f"evidence on ≥ {ed} days", "evidence.min_days"),
    ]
    n = _notes()
    return [{"rule": r, "value": v, "why": (n.get(k) or {}).get("why", ""),
             "tested?": (n.get(k) or {}).get("tested", ""), "config": cfg}
            for k, r, v, cfg in rows]


_GRAPH = ("<div style='font-size:13px;color:#4b5563;margin-top:10px;'>"
          "<b style='font-size:12px;color:#6b7280;text-transform:uppercase;'>Panel check graph</b>"
          " &nbsp;green = explained by your panel · light green = a panel variant ± 1 change · "
          "red = a lineage the deep scan found · grey = unexplained</div>")


def render_how_to_read() -> None:
    """The ⓘ popover button. Call it in the row of the result tabs."""
    with st.popover("ⓘ How to read this page", use_container_width=True):
        t_col, t_rule = st.tabs(["Colours", "Rules and thresholds"])
        with t_col:
            from components.variants_table import legend_full_html
            st.markdown(legend_full_html() + _GRAPH, unsafe_allow_html=True)
        with t_rule:
            st.caption("Read from cooc_config.yaml — change a value there and the page "
                       "follows. “tested?” is updated there too after a check.")
            import pandas as pd
            st.dataframe(pd.DataFrame(rules_rows()), hide_index=True, width="stretch")
