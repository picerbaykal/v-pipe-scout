"""
Primer & Probe Checker
======================
Checks whether SARS-CoV-2 primer/probe binding sites have mutations
in Swiss wastewater data, using real sampling dates only.

Follows the abundance_cooc.py conventions:
- WiseLoculusLapis for LAPIS queries
- asyncio.run() for async calls
- ppc_ prefix for session state keys
- Inline location multiselect (not sidebar)
- Celery for heavy computation (position scanning)
"""

import asyncio
import io
import json
import logging
import os
import re
from datetime import date, datetime, timedelta

import pandas as pd
import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
from celery import Celery
from redis import Redis

from api.wiseloculus import WiseLoculusLapis
from utils.config import get_wiseloculus_url

# ── Module-level Celery / Redis setup (same pattern as abundance_cooc) ─────────
celery_app = Celery("tasks", broker="redis://redis:6379/0")
redis_client = Redis(host="redis", port=6379, db=0)

logger = logging.getLogger(__name__)

# ── Reference genome path ──────────────────────────────────────────────────────
REF_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "SARS-CoV-2.primer.fasta")

# ── Helper functions ───────────────────────────────────────────────────────────

@st.cache_data(ttl=3600)
def load_reference():
    """Load Wuhan reference genome NC_045512.2."""
    if not os.path.exists(REF_PATH):
        return None, f"Reference file not found at: {REF_PATH}"
    record = next(SeqIO.parse(REF_PATH, "fasta"))
    return str(record.seq).upper(), None

@st.cache_data(ttl=3600)
def fetch_locations():
    """Fetch available wastewater locations from LAPIS."""
    try:
        lapis = WiseLoculusLapis(get_wiseloculus_url())
        locations = lapis.fetch_locations()
        return sorted(locations), None
    except Exception as e:
        return None, str(e)

def get_piece_type(name):
    """Detect if a primer/probe piece is forward, reverse or probe."""
    clean = name.split("::")[0].upper()
    if clean.endswith("-P") or clean.endswith("_P"):
        return "probe"
    elif clean.endswith("-F") or clean.endswith("_F") or "_LEFT" in clean:
        return "forward"
    elif clean.endswith("-R") or clean.endswith("_R") or "_RIGHT" in clean:
        return "reverse"
    return "unknown"

def extract_positions_from_header(header):
    """ARTIC format: extract start/end from ::NC_045512.2:47-78"""
    match = re.search(r"NC_045512\.2:(\d+)-(\d+)", header.strip())
    if match:
        return int(match.group(1)), int(match.group(2)), "+"
    return None

def search_reference(primer_seq, ref_seq):
    """CDC format: find primer position by string search."""
    seq_len = len(primer_seq)
    pos = ref_seq.find(primer_seq)
    if pos != -1:
        return pos + 1, pos + seq_len, "+"
    rev_comp = str(Seq(primer_seq).reverse_complement())
    pos = ref_seq.find(rev_comp)
    if pos != -1:
        return pos + 1, pos + seq_len, "-"
    return None

def build_genspectrum_urls(start, end, reference, location=None):
    """Build GenSpectrum Link A and Link B for a primer/probe region."""
    region = reference[start - 1:end]
    positions = list(range(start, end + 1))
    base = "https://genspectrum.org/swiss-wastewater/covid?"
    if location and location != "All sites":
        from urllib.parse import quote_plus
        base += f"locationName={quote_plus(location)}&"
    base += "analysisMode=manual&sequenceType=nucleotide&mutations="
    link_a = base + "%7C".join([f"{pos}{nuc}" for pos, nuc in zip(positions, region)])
    link_b = base + "%7C".join([str(pos) for pos in positions])
    return link_a, link_b

def process_time_series(df, results):
    """
    Process raw mutation DataFrame into monthly proportions per position.
    Returns dict: {position: {primer_name: str, piece_type: str, monthly: {month: {mut: frac}}}}
    """
    if df.empty:
        return {}

    # add month column
    df = df.copy()
    df["month"] = df["date"].dt.to_period("M").astype(str)

    position_data = {}

    for r in results:
        if r["Start"] is None:
            continue
        name = r["Name"]
        start = r["Start"]
        end = r["End"]
        piece_type = get_piece_type(name)

        # filter to positions in this primer/probe region
        region_df = df[
            df["pos"].apply(lambda x: start <= int(''.join(filter(str.isdigit, str(x)))) <= end)
        ]

        if region_df.empty:
            continue

        for pos_mut in region_df["pos"].unique():
            pos_num = int(''.join(filter(str.isdigit, str(pos_mut))))
            if pos_num not in position_data:
                position_data[pos_num] = {
                    "name": name,
                    "piece_type": piece_type,
                    "monthly": {}
                }

            # group by month
            pos_df = region_df[region_df["pos"] == pos_mut]
            for month, group in pos_df.groupby("month"):
                avg_frac = group["frac"].mean()
                if month not in position_data[pos_num]["monthly"]:
                    position_data[pos_num]["monthly"][month] = {}
                position_data[pos_num]["monthly"][month][pos_mut] = avg_frac

    return position_data

def build_html_heatmap(position_data):
    """Build custom HTML heatmap from position data."""
    if not position_data:
        return "<p>No mutations found in primer-probe regions.</p>"

    # collect all months with data
    all_months = sorted(set(
        m for info in position_data.values()
        for m in info["monthly"].keys()
    ))

    if not all_months:
        return "<p>No data available for the selected period.</p>"

    def cell_color(proportion):
        if proportion <= 0:
            return "background:#ffffff; color:#999;"
        elif proportion < 0.1:
            return "background:#fff8f0; color:#bf6000;"
        elif proportion < 0.5:
            return "background:#fac775; color:#412402;"
        elif proportion < 0.8:
            return "background:#d85a30; color:#ffffff;"
        else:
            return "background:#a32d2d; color:#ffffff;"

    def piece_badge(piece_type):
        if piece_type == "probe":
            return '<span style="background:#faeeda;color:#854f0b;font-size:11px;padding:2px 7px;border-radius:4px;font-weight:500;">probe</span>'
        elif piece_type == "forward":
            return '<span style="background:#e6f1fb;color:#185fa5;font-size:11px;padding:2px 7px;border-radius:4px;font-weight:500;">forward</span>'
        elif piece_type == "reverse":
            return '<span style="background:#e1f5ee;color:#0f6e56;font-size:11px;padding:2px 7px;border-radius:4px;font-weight:500;">reverse</span>'
        return ""

    html = """
<style>
.pp-table{border-collapse:collapse;font-size:13px;font-family:sans-serif;}
.pp-table th{text-align:center;padding:8px 10px;color:#666;font-weight:400;border-bottom:1px solid #e5e5e5;white-space:nowrap;}
.pp-table th.pos-col{text-align:left;width:auto;}
.pp-table td{padding:4px 6px;border-bottom:1px solid #f0f0f0;}
.pp-table td.pos-label{padding:6px 10px;white-space:nowrap;padding-right:24px;}
.pp-cell{border-radius:4px;height:32px;display:flex;align-items:center;justify-content:center;font-size:12px;font-weight:500;cursor:default;position:relative;min-width:60px;}
.pp-cell .tooltip{display:none;position:absolute;bottom:calc(100% + 8px);left:50%;transform:translateX(-50%);background:#fff;border:1px solid #ddd;border-radius:6px;padding:10px 14px;min-width:200px;z-index:100;box-shadow:0 4px 12px rgba(0,0,0,0.1);text-align:left;pointer-events:none;}
.pp-cell:hover .tooltip{display:block;}
.pp-cell .tt-title{font-size:12px;font-weight:500;color:#333;margin-bottom:4px;}
.pp-cell .tt-date{font-size:11px;color:#888;margin-bottom:8px;}
.pp-cell .tt-row{display:flex;justify-content:space-between;font-size:12px;color:#333;margin-bottom:3px;}
.pp-cell .tt-div{border-top:1px solid #eee;margin:6px 0;}
.pp-cell .tt-ref{font-size:11px;color:#999;display:flex;justify-content:space-between;}
.pp-legend{display:flex;gap:16px;margin-top:12px;flex-wrap:wrap;align-items:center;}
.pp-legend-item{display:flex;align-items:center;gap:6px;font-size:12px;color:#666;}
</style>
<table class="pp-table">
<thead><tr><th class="pos-col">Position</th>
"""
    for m in all_months:
        html += f"<th>{m}</th>"
    html += "</tr></thead><tbody>"

    for pos in sorted(position_data.keys()):
        info = position_data[pos]
        name = info["name"]
        piece_type = info["piece_type"]
        monthly = info["monthly"]

        html += f'<tr><td class="pos-label">{piece_badge(piece_type)} <span style="color:#333;margin-left:6px;">{pos} ({name})</span></td>'

        for month in all_months:
            mut_values = list(monthly.get(month, {}).items())
            total = min(sum(v for _, v in mut_values), 1.0) if mut_values else 0.0
            style = cell_color(total)
            cell_text = "—" if total <= 0 else f"{total:.0%}"

            if mut_values:
                tooltip_rows = ""
                for mut, val in sorted(mut_values, key=lambda x: -x[1]):
                    ref = str(mut)[0]
                    alt = str(mut)[-1]
                    tooltip_rows += f'<div class="tt-row"><span>{ref}→{alt}</span><span style="font-weight:500;">{val:.1%}</span></div>'
                ref_remaining = max(0, 1.0 - total)
                tooltip_html = f"""<div class="tooltip">
<div class="tt-title">{pos} — {name} ({piece_type})</div>
<div class="tt-date">{month} — total: {total:.1%}</div>
<div class="tt-div"></div>
{tooltip_rows}
<div class="tt-div"></div>
<div class="tt-ref"><span>reference remaining</span><span>{ref_remaining:.1%}</span></div>
</div>"""
            else:
                tooltip_html = ""

            html += f'<td><div class="pp-cell" style="{style}">{cell_text}{tooltip_html}</div></td>'

        html += "</tr>"

    html += """</tbody></table>
<div class="pp-legend">
  <span style="font-size:12px;color:#888;">Mutation prevalence:</span>
  <div class="pp-legend-item"><span style="background:#fff8f0;width:16px;height:16px;display:inline-block;border-radius:3px;border:1px solid #eee;"></span> &lt;10%</div>
  <div class="pp-legend-item"><span style="background:#fac775;width:16px;height:16px;display:inline-block;border-radius:3px;"></span> 10–50%</div>
  <div class="pp-legend-item"><span style="background:#d85a30;width:16px;height:16px;display:inline-block;border-radius:3px;"></span> 50–80%</div>
  <div class="pp-legend-item"><span style="background:#a32d2d;width:16px;height:16px;display:inline-block;border-radius:3px;"></span> &gt;80%</div>
</div>"""
    return html

def build_summary_table(position_data, threshold):
    """Build summary table — one row per primer/probe piece."""
    # group by primer/probe name
    primer_summary = {}
    for pos, info in position_data.items():
        name = info["name"]
        piece_type = info["piece_type"]
        if name not in primer_summary:
            primer_summary[name] = {
                "piece_type": piece_type,
                "worst_mutation": None,
                "worst_proportion": 0.0,
            }
        # find worst recent proportion across all months
        for month, muts in info["monthly"].items():
            for mut, frac in muts.items():
                if frac > primer_summary[name]["worst_proportion"]:
                    primer_summary[name]["worst_proportion"] = frac
                    primer_summary[name]["worst_mutation"] = mut

    rows = []
    for name, summary in primer_summary.items():
        worst_proportion = summary["worst_proportion"]
        worst_mutation = summary["worst_mutation"]
        if worst_proportion >= 0.5:
            status = "🚨 Critical"
        elif worst_proportion >= threshold:
            status = "⚠️ Warning"
        else:
            status = "✅ Good"
        rows.append({
            "Primer/Probe": name,
            "Type": summary["piece_type"],
            "Worst mutation": worst_mutation or "none",
            "Recent proportion": f"{worst_proportion:.1%}" if worst_proportion > 0 else "0.0%",
            "Status": status,
        })

    # add clean primers not in position_data
    return pd.DataFrame(rows)

# ── Main app function ──────────────────────────────────────────────────────────

def app():
    st.title("🔬 Primer & Probe Checker")
    st.markdown(
        """
        Upload a FASTA file containing primer and probe sequences.
        The app checks whether binding sites have mutations in Swiss wastewater data.
        """
    )
    st.divider()

    # ── Load reference ─────────────────────────────────────────────────────────
    reference, ref_error = load_reference()
    if ref_error:
        st.error(f"Could not load reference genome: {ref_error}")
        return

    # ── Fetch locations ────────────────────────────────────────────────────────
    locations, loc_error = fetch_locations()
    if loc_error or not locations:
        st.error(f"Could not load locations from LAPIS: {loc_error}")
        locations = []

    # ── Filters — inline on page (not sidebar) ─────────────────────────────────
    col1, col2, col3 = st.columns([2, 1, 1])

    with col1:
        selected_locations = st.multiselect(
            "Wastewater sites",
            options=locations,
            default=locations[:1] if locations else [],
            key="ppc_locations",
            help="Select one or more wastewater treatment plants"
        )

    with col2:
        date_from = st.date_input(
            "From",
            value=date.today() - timedelta(days=180),
            key="ppc_date_from",
        )

    with col3:
        date_to = st.date_input(
            "To",
            value=date.today(),
            key="ppc_date_to",
        )

    st.divider()

    # ── File upload ────────────────────────────────────────────────────────────
    uploaded_file = st.file_uploader(
        "Upload primer/probe FASTA file",
        type=["fasta", "fa"],
        help="e.g. CDC N1/N2 diagnostic primers, or ARTIC amplicon primers",
        key="ppc_fasta_upload",
    )

    if uploaded_file is None:
        st.info("👆 Upload a FASTA file to get started.")
        return

    if not selected_locations:
        st.warning("Please select at least one wastewater site.")
        return

    # ── Parse FASTA ────────────────────────────────────────────────────────────
    raw_content = uploaded_file.read().decode("utf-8")
    sequences = {}
    for record in SeqIO.parse(io.StringIO(raw_content), "fasta"):
        sequences[record.id] = str(record.seq).upper()

    if not sequences:
        st.error("No sequences found in the FASTA file. Please check the format.")
        return

    with st.expander(f"Show parsed sequences ({len(sequences)} found)"):
        for name, seq in sequences.items():
            clean = name.split("::")[0]
            st.text(f"{clean} — {len(seq)} bp — {get_piece_type(name)}")

    # ── Find positions on reference ────────────────────────────────────────────
    results = []
    for name, seq in sequences.items():
        clean_name = name.split("::")[0]
        positions = extract_positions_from_header(name)
        if positions:
            start, end, strand = positions
            method = "from header"
        else:
            found = search_reference(seq, reference)
            if found:
                start, end, strand = found
                method = "string search"
            else:
                results.append({
                    "Name": clean_name, "Start": None, "End": None,
                    "Strand": None, "Method": "not found",
                    "Link A": None, "Link B": None,
                })
                continue

        link_a, link_b = build_genspectrum_urls(
            start, end, reference,
            location=selected_locations[0] if len(selected_locations) == 1 else None
        )
        results.append({
            "Name": clean_name, "Start": start, "End": end,
            "Strand": strand, "Method": method,
            "Link A": link_a, "Link B": link_b,
        })

    not_found = [r for r in results if r["Start"] is None]
    found = [r for r in results if r["Start"] is not None]

    if not_found:
        st.warning(
            f"{len(not_found)} sequence(s) not found on Wuhan reference "
            "(likely designed for a newer variant genome)."
        )

    if not found:
        st.error("No sequences could be mapped to the reference genome.")
        return

    # ── Collect positions to query ─────────────────────────────────────────────
    all_positions = set()
    for r in found:
        all_positions.update(range(r["Start"], r["End"] + 1))

    # ── Fetch mutation data from LAPIS ─────────────────────────────────────────
    date_range = (
        datetime.combine(date_from, datetime.min.time()),
        datetime.combine(date_to, datetime.max.time()),
    )

    all_dfs = []
    lapis = WiseLoculusLapis(get_wiseloculus_url())

    with st.spinner("Fetching mutation data from LAPIS..."):
        for location in selected_locations:
            try:
                df = asyncio.run(
                    lapis.get_mutations_at_positions(
                        locationName=location,
                        date_range=date_range,
                        positions=all_positions,
                    )
                )
                if not df.empty:
                    df["location"] = location
                    all_dfs.append(df)
            except Exception as e:
                st.error(f"Error fetching data for {location}: {e}")
                logger.error(f"LAPIS error for {location}: {e}")

    if not all_dfs:
        st.warning("No mutation data found for the selected sites and date range.")
        return

    combined_df = pd.concat(all_dfs, ignore_index=True)

    # ── Process time series ────────────────────────────────────────────────────
    position_data = process_time_series(combined_df, found)

    if not position_data:
        st.success("✅ No mutations found in primer-probe regions for the selected period.")
        return

    # ── Heatmap ────────────────────────────────────────────────────────────────
    st.subheader("Mutation prevalence over time")
    html = build_html_heatmap(position_data)
    st.html(html)

    # ── Summary table ──────────────────────────────────────────────────────────
    st.subheader("Summary")
    threshold = st.slider(
        "Warning threshold",
        min_value=0.01, max_value=1.0, value=0.10, step=0.01,
        format="%.2f",
        key="ppc_threshold",
        help="Proportion above which a mutation is flagged as Warning or Critical"
    )
    st.caption("Controls when a mutation is flagged as Warning or Critical in the table below.")

    summary_df = build_summary_table(position_data, threshold)
    st.dataframe(summary_df, use_container_width=True, hide_index=True)

    # ── GenSpectrum links ──────────────────────────────────────────────────────
    st.subheader("GenSpectrum links")

    if len(selected_locations) > 1:
        st.info("Links below use the first selected location for GenSpectrum filtering.")

    for r in found:
        with st.expander(f"🧬 {r['Name']} — positions {r['Start']}–{r['End']} ({r['Strand']} strand)"):
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("**Link A — reference sequence prevalence**")
                st.markdown(f"[Open in GenSpectrum ↗]({r['Link A']})")
                st.caption("Low % = mutation dominant at these positions")
            with col2:
                st.markdown("**Link B — any mutation prevalence**")
                st.markdown(f"[Open in GenSpectrum ↗]({r['Link B']})")
                st.caption("High % = mutation present at these positions")