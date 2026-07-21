"""
Primer & Probe Checker
======================
Checks whether SARS-CoV-2 primer/probe binding sites have mutations
in Swiss wastewater data.

Uses two LAPIS endpoints:
1. /sample/nucleotideMutations — find which mutations exist at primer positions
2. /component/nucleotideMutationsOverTime — fetch monthly time series
"""

import asyncio
import io
import logging
import re
from datetime import date, datetime, timedelta

import aiohttp
import pandas as pd
import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
from dateutil.relativedelta import relativedelta

from utils.config import get_wiseloculus_url

logger = logging.getLogger(__name__)

# ── Helper functions ───────────────────────────────────────────────────────────

def get_piece_type(name):
    clean = name.split("::")[0].upper()
    if "-P-" in clean or "-P_" in clean or clean.endswith("-P") or clean.endswith("_P"):
        return "probe"
    elif clean.endswith("-F") or clean.endswith("_F") or "_LEFT" in clean:
        return "forward"
    elif clean.endswith("-R") or clean.endswith("_R") or "_RIGHT" in clean:
        return "reverse"
    return "unknown"

def extract_positions_from_header(header):
    match = re.search(r"NC_045512\.2:(\d+)-(\d+)", header.strip())
    if match:
        return int(match.group(1)), int(match.group(2)), "+"
    return None

def search_reference(primer_seq, ref_seq):
    seq_len = len(primer_seq)
    pos = ref_seq.find(primer_seq)
    if pos != -1:
        return pos + 1, pos + seq_len, "+"
    rev_comp = str(Seq(primer_seq).reverse_complement())
    pos = ref_seq.find(rev_comp)
    if pos != -1:
        return pos + 1, pos + seq_len, "-"
    return None

def get_primer_letters(primer_seq, start, strand):
    """
    Map each genome position to the letter the primer expects there.
    For reverse strand primers, reverse complement first.
    Returns dict: {genome_position: letter}
    """
    if strand == "-":
        primer_seq = str(Seq(primer_seq).reverse_complement())
    return {
        start + i: primer_seq[i]
        for i in range(len(primer_seq))
    }

from Bio import Align

def align_primer_to_reference(primer_seq, reference, min_identity=0.80):
    """
    Fallback alignment when exact string search fails.
    Uses Biopython PairwiseAligner (local alignment).
    Returns (start, end, strand, identity) or None if below threshold.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    # try forward strand
    alignments = aligner.align(reference, primer_seq)
    if alignments:
        best = alignments[0]
        ref_start = best.aligned[0][0][0]
        ref_end = best.aligned[0][-1][-1]
        matches = sum(
            reference[ref_start + i] == primer_seq[i]
            for i in range(len(primer_seq))
            if ref_start + i < len(reference)
        )
        identity = matches / len(primer_seq)
        if identity >= min_identity:
            return ref_start + 1, ref_start + len(primer_seq), "+", identity

    # try reverse complement
    rev_comp = str(Seq(primer_seq).reverse_complement())
    alignments = aligner.align(reference, rev_comp)
    if alignments:
        best = alignments[0]
        ref_start = best.aligned[0][0][0]
        matches = sum(
            reference[ref_start + i] == rev_comp[i]
            for i in range(len(rev_comp))
            if ref_start + i < len(reference)
        )
        identity = matches / len(rev_comp)
        if identity >= min_identity:
            return ref_start + 1, ref_start + len(primer_seq), "-", identity

    return None

def build_genspectrum_urls(start, end, reference, location=None):
    region = reference[start - 1:end]
    positions = list(range(start, end + 1))
    base = "https://genspectrum.org/swiss-wastewater/covid?"
    if location:
        from urllib.parse import quote_plus
        base += f"locationName={quote_plus(location)}&"
    base += "analysisMode=manual&sequenceType=nucleotide&mutations="
    link_a = base + "%7C".join([f"{pos}{nuc}" for pos, nuc in zip(positions, region)])
    link_b = base + "%7C".join([str(pos) for pos in positions])
    return link_a, link_b

def build_monthly_ranges(date_from, date_to):
    ranges = []
    current = date_from.replace(day=1)
    while current <= date_to:
        last = (current + relativedelta(months=1)) - timedelta(days=1)
        ranges.append({
            "dateFrom": current.strftime("%Y-%m-%d"),
            "dateTo": min(last, date_to).strftime("%Y-%m-%d")
        })
        current = current + relativedelta(months=1)
    return ranges

# ── LAPIS fetch functions ──────────────────────────────────────────────────────

@st.cache_data(ttl=3600)
def fetch_locations_list():
    async def _fetch():
        async with aiohttp.ClientSession() as session:
            async with session.get(
                f"{get_wiseloculus_url()}/sample/aggregated",
                params={"fields": "locationName"},
                headers={"accept": "application/json"}
            ) as response:
                if response.status != 200:
                    return []
                data = await response.json()
                return [
                    d["locationName"]
                    for d in data.get("data", [])
                    if d.get("locationName")
                ]
    try:
        locations = asyncio.run(_fetch())
        return sorted(locations), None
    except Exception as e:
        return [], str(e)

@st.cache_data(ttl=3600)
def fetch_primer_mutations(location, positions_tuple):
    """
    Step 1: fetch which mutations exist at primer/probe positions.
    location=None means all sites.
    """
    async def _fetch():
        params = {"limit": 100000}
        if location:
            params["locationName"] = location
        async with aiohttp.ClientSession() as session:
            async with session.get(
                f"{get_wiseloculus_url()}/sample/nucleotideMutations",
                params=params,
                headers={"accept": "application/json"}
            ) as response:
                if response.status != 200:
                    return []
                data = await response.json()
                positions = set(positions_tuple)
                return [
                    entry["mutation"]
                    for entry in data.get("data", [])
                    if entry.get("position") in positions
                ]
    try:
        return asyncio.run(_fetch())
    except Exception as e:
        logger.error(f"fetch_primer_mutations error: {e}")
        return []

@st.cache_data(ttl=3600)
def fetch_time_series(mutations_tuple, location, date_from_str, date_to_str):
    """
    Step 2: fetch monthly time series for specific mutations.
    location=None means all sites.
    """
    mutations_list = list(mutations_tuple)
    if not mutations_list:
        return {}

    date_from = datetime.strptime(date_from_str, "%Y-%m-%d")
    date_to = datetime.strptime(date_to_str, "%Y-%m-%d")
    date_ranges = build_monthly_ranges(date_from, date_to)

    filters = {}
    if location:
        filters["locationName"] = location

    async def _fetch():
        async with aiohttp.ClientSession(
            timeout=aiohttp.ClientTimeout(total=120)
        ) as session:
            async with session.post(
                f"{get_wiseloculus_url()}/component/nucleotideMutationsOverTime",
                json={
                    "filters": filters,
                    "dateRanges": date_ranges,
                    "includeMutations": mutations_list,
                    "dateField": "samplingDate"
                },
                headers={"content-type": "application/json"}
            ) as response:
                if response.status != 200:
                    logger.error(f"nucleotideMutationsOverTime failed: {response.status}")
                    return {}
                data = (await response.json())["data"]
                mutations = data["mutations"]
                date_labels = [r["dateFrom"][:7] for r in data["dateRanges"]]
                result = {}
                for i, mut in enumerate(mutations):
                    result[mut] = {}
                    for j, label in enumerate(date_labels):
                        cell = data["data"][i][j]
                        if cell["coverage"] > 0:
                            result[mut][label] = {
                                "proportion": cell["count"] / cell["coverage"],
                                "coverage": cell["coverage"],
                                "count": cell["count"],
                            }
                        else:
                            result[mut][label] = {
                                "proportion": 0.0,
                                "coverage": 0,
                                "count": 0,
                            }
                return result

    try:
        return asyncio.run(_fetch())
    except Exception as e:
        logger.error(f"fetch_time_series error: {e}")
        return {}

# ── Data processing ────────────────────────────────────────────────────────────

def process_time_series(time_series, results):
    """Process time series dict into position_data for heatmap."""
    position_data = {}

    for r in results:
        if r["Start"] is None:
            continue
        name = r["Name"]
        start = r["Start"]
        end = r["End"]
        piece_type = get_piece_type(name)

        for mut, proportions in time_series.items():
            try:
                pos_num = int(re.search(r'\d+', mut).group())
            except Exception:
                continue

            if start <= pos_num <= end:
                if pos_num not in position_data:
                    position_data[pos_num] = {
                        "name": name,
                        "piece_type": piece_type,
                        "monthly": {}
                    }
                for month, data in proportions.items():
                    frac = data["proportion"]
                    coverage = data["coverage"]
                    if frac > 0:
                        if month not in position_data[pos_num]["monthly"]:
                            position_data[pos_num]["monthly"][month] = {}
                        position_data[pos_num]["monthly"][month][mut] = {
                            "proportion": frac,
                            "coverage": coverage,
                        }

    return position_data


def get_dominant_letters(time_series, reference, results):
    """
    For each position with a mutation, determine the dominant letter
    in the most recent month.
    Returns dict: {genome_position: dominant_letter}
    """
    dominant = {}

    for mut, proportions in time_series.items():
        # extract position and letters from mutation string e.g. C28311T
        try:
            pos = int(re.search(r'\d+', mut).group())
            ref_letter = mut[0]  # C in C28311T
            alt_letter = mut[-1]  # T in C28311T
        except Exception:
            continue

        # get most recent non-zero proportion
        recent_frac = 0.0
        for month in sorted(proportions.keys(), reverse=True):
            val = proportions[month]
            if isinstance(val, dict):
                frac = val.get("proportion", 0.0)
            else:
                frac = val
            if frac > 0:
                recent_frac = frac
                break

        # dominant letter = whichever is more common
        if recent_frac > 0.5:
            dominant[pos] = alt_letter  # mutation is dominant
        else:
            dominant[pos] = ref_letter  # reference is dominant

    return dominant

# ── Heatmap ────────────────────────────────────────────────────────────────────

def build_html_heatmap(position_data, dominant_letters, results):
    if not position_data:
        return "<p>No mutations found in primer-probe regions.</p>"

    # build primer letters lookup: {pos: {primer_name: letter}}
    primer_letters_by_pos = {}
    for r in results:
        if r["Start"] is None:
            continue
        for pos, letter in r.get("PrimerLetters", {}).items():
            if pos not in primer_letters_by_pos:
                primer_letters_by_pos[pos] = {}
            primer_letters_by_pos[pos][r["Name"]] = letter

    all_months = sorted(set(
        m for info in position_data.values()
        for m in info["monthly"].keys()
    ))

    if not all_months:
        return "<p>No data available for the selected period.</p>"

    def cell_color(total, pos, name, coverage=0, min_coverage=1000):
        if total <= 0:
            return "background:#ffffff; color:#999;"
        if coverage < min_coverage:
            return "background:#e9ecef; color:#6c757d;"
        primer_letter = primer_letters_by_pos.get(pos, {}).get(name)
        dominant = dominant_letters.get(pos)
        if primer_letter and dominant and primer_letter == dominant:
            return "background:#d4edda; color:#155724;"
        elif total < 0.5:
            return "background:#fac775; color:#412402;"
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
            total = min(sum(p for _, p, _ in mut_values), 1.0) if mut_values else 0.0
            max_cov = max((c for _, _, c in mut_values), default=0) if mut_values else 0
            style = cell_color(total, pos, name, max_cov)

            if mut_values:
                tooltip_rows = ""
                for mut, val, cov in sorted(mut_values, key=lambda x: -x[1]):
                    ref = str(mut)[0]
                    alt = str(mut)[-1]
                    tooltip_rows += f'<div class="tt-row"><span>{ref}→{alt}</span><span style="font-weight:500;">{val:.1%}</span></div>'
                ref_remaining = max(0, 1.0 - total)
                tooltip_html = f"""<div class="tooltip">
<div class="tt-title">{pos} — {name} ({piece_type})</div>
<div class="tt-date">{month} — total: {total:.1%}</div>
<div class="tt-div"></div>
{tooltip_rows}
<div class="tt-ref"><span>reference remaining</span><span>{ref_remaining:.1%}</span></div>
<div class="tt-div"></div>
<div class="tt-ref"><span>coverage</span><span>{min_coverage:,} reads</span></div>
</div>"""
            else:
                tooltip_html = ""

            html += f'<td><div class="pp-cell" style="{style}">{cell_text}{tooltip_html}</div></td>'

        html += "</tr>"

    html += """</tbody></table>
<div class="pp-legend">
  <span style="font-size:12px;color:#888;">Cell meaning:</span>
  <div class="pp-legend-item"><span style="background:#ffffff;width:16px;height:16px;display:inline-block;border-radius:3px;border:1px solid #eee;"></span> No mutation</div>
  <div class="pp-legend-item"><span style="background:#d4edda;width:16px;height:16px;display:inline-block;border-radius:3px;"></span> Mutation — primer matches ✅</div>
  <div class="pp-legend-item"><span style="background:#fac775;width:16px;height:16px;display:inline-block;border-radius:3px;"></span> Mismatch &lt;50%</div>
  <div class="pp-legend-item"><span style="background:#a32d2d;width:16px;height:16px;display:inline-block;border-radius:3px;"></span> Mismatch &gt;50%</div>
  <div class="pp-legend-item"><span style="background:#e9ecef;width:16px;height:16px;display:inline-block;border-radius:3px;"></span> Low coverage — unreliable</div>
</div>"""
    return html

# ── Summary table ──────────────────────────────────────────────────────────────
def build_summary_table(position_data, found, threshold, dominant_letters):
    """Build summary table with primer vs circulating virus comparison."""
    primer_summary = {}
    for r in found:
        if r["Start"] is None:
            continue
        name = r["Name"]
        piece_type = get_piece_type(name)
        primer_letters = r.get("PrimerLetters", {})

        primer_summary[name] = {
            "piece_type": piece_type,
            "worst_mutation": None,
            "worst_proportion": 0.0,
            "worst_coverage": 0,
            "primer_match": True,  # assume match until proven otherwise
            "mismatch_positions": [],
        }

        # check each position in this primer/probe against dominant circulating letter
        for pos, primer_letter in primer_letters.items():
            dominant = dominant_letters.get(pos)
            if dominant and dominant != primer_letter:
                primer_summary[name]["primer_match"] = False
                primer_summary[name]["mismatch_positions"].append(
                    f"{pos}({primer_letter}→{dominant})"
                )

    # update with mutation data
    for pos, info in position_data.items():
        name = info["name"]
        if name not in primer_summary:
            continue
        for month, muts in info["monthly"].items():
            for mut, data in muts.items():
                frac = data["proportion"] if isinstance(data, dict) else data
                cov = data.get("coverage", 0) if isinstance(data, dict) else 0
                if frac > primer_summary[name]["worst_proportion"]:
                    primer_summary[name]["worst_proportion"] = frac
                    primer_summary[name]["worst_mutation"] = mut
                    primer_summary[name]["worst_coverage"] = cov

    rows = []
    for name, summary in primer_summary.items():
        worst_proportion = summary["worst_proportion"]
        worst_mutation = summary["worst_mutation"]

        # primer vs virus status
        if not summary["primer_match"]:
            mismatches = ", ".join(summary["mismatch_positions"])
            primer_vs_virus = f"❌ Mismatch at: {mismatches}"
        else:
            primer_vs_virus = "✅ Matches circulating virus"

        # mutation status — based on primer match, not just proportion vs Wuhan
        if not summary["primer_match"]:
            if worst_proportion >= 0.5:
                mut_status = "🚨 Critical"
            else:
                mut_status = "⚠️ Warning"
        else:
            mut_status = "✅ Good"

        # coverage display
        cov = summary['worst_coverage']
        if cov == 0:
            coverage_display = "N/A"
        elif cov < 1000:
            coverage_display = f"⚠️ {cov:,} reads"
        else:
            coverage_display = f"{cov:,} reads"

        rows.append({
            "Primer/Probe": name,
            "Type": summary["piece_type"],
            "Worst mutation": worst_mutation or "none",
            "Recent proportion": f"{worst_proportion:.1%}" if worst_proportion > 0 else "0.0%",
            "Coverage": coverage_display,
            "Primer vs virus": primer_vs_virus,
            "Status": mut_status,
        })

    return pd.DataFrame(rows)


# ── Main app function ──────────────────────────────────────────────────────────

def app():
    st.title("🔬 Primer & Probe Checker")
    st.markdown(
        """
        Upload a reference genome and primer/probe sequences.
        The app checks whether binding sites have mutations in Swiss wastewater data.
        """
    )
    st.divider()

    # ── Fetch locations ────────────────────────────────────────────────────────
    locations, loc_error = fetch_locations_list()
    if loc_error or not locations:
        st.error(f"Could not load locations from LAPIS: {loc_error}")
        locations = []

    # ── Filters — inline on page ───────────────────────────────────────────────
    col1, col2, col3 = st.columns([2, 1, 1])

    with col1:
        selected_location = st.selectbox(
            "Wastewater site",
            options=["All sites"] + locations,
            index=0,
            key="ppc_location",
            help="Select a wastewater treatment plant or all sites combined"
        )

    with col2:
        date_from = st.date_input(
            "From",
            value=date.today() - timedelta(days=210),
            key="ppc_date_from",
        )

    with col3:
        date_to = st.date_input(
            "To",
            value=date.today(),
            key="ppc_date_to",
        )

    # location filter — None means all sites
    location_filter = None if selected_location == "All sites" else selected_location

    st.divider()

    # ── Step 1: Reference genome ───────────────────────────────────────────────
    st.subheader("Step 1 — Upload reference genome")
    ref_file = st.file_uploader(
        "Reference genome (FASTA)",
        type=["fasta", "fa"],
        help="e.g. Wuhan reference NC_045512.2 (29,903 bp)",
        key="ppc_ref_upload",
    )

    if ref_file is None:
        st.info("👆 Upload a reference genome to get started.")
        return

    try:
        ref_content = ref_file.read().decode("utf-8")
        ref_record = next(SeqIO.parse(io.StringIO(ref_content), "fasta"))
        reference = str(ref_record.seq).upper()
    except Exception as e:
        st.error(f"Could not parse file: {e}")
        return

    if len(reference) < 1000:
        st.error(
            f"This file looks like a primer file ({len(reference):,} bp), not a reference genome. "
            "Please upload a full genome reference like NC_045512.2 (29,903 bp)."
        )
        return

    st.success(f"✅ Reference genome accepted: {ref_record.id} ({len(reference):,} bp)")

    # ── Step 2: Primer/probe FASTA ─────────────────────────────────────────────
    st.subheader("Step 2 — Upload primer/probe FASTA")
    primer_file = st.file_uploader(
        "Primer/probe FASTA",
        type=["fasta", "fa"],
        help="e.g. CDC N1/N2 diagnostic primers, or ARTIC amplicon primers",
        key="ppc_fasta_upload",
    )

    if primer_file is None:
        st.info("👆 Now upload your primer/probe FASTA file.")
        return

    # ── Parse primers ──────────────────────────────────────────────────────────
    raw_content = primer_file.read().decode("utf-8")
    sequences = {}
    for record in SeqIO.parse(io.StringIO(raw_content), "fasta"):
        sequences[record.id] = str(record.seq).upper()

    if not sequences:
        st.error("No sequences found in the primer FASTA file. Please check the format.")
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
            found_pos = search_reference(seq, reference)
            if found_pos:
                start, end, strand = found_pos
                method = "string search"
            else:
                # fallback to alignment
                aligned = align_primer_to_reference(seq, reference)
                if aligned:
                    start, end, strand, identity = aligned
                    method = f"alignment ({identity:.0%} identity)"
                else:
                    results.append({
                        "Name": clean_name, "Start": None, "End": None,
                        "Strand": None, "Method": "not found",
                        "PrimerLetters": get_primer_letters(seq, start, strand),
                        "Link A": None, "Link B": None,
                    })
                    continue

        link_a, link_b = build_genspectrum_urls(
            start, end, reference, location=location_filter
        )
        results.append({
            "Name": clean_name, "Start": start, "End": end,
            "Strand": strand, "Method": method,
            "PrimerLetters": get_primer_letters(seq, start, strand),
            "Link A": link_a, "Link B": link_b,
        })

    not_found = [r for r in results if r["Start"] is None]
    found = [r for r in results if r["Start"] is not None]

    if found:
        with st.expander(f"Positions found ({len(found)} sequences)"):
            for r in found:
                st.text(f"{r['Name']} — {r['Start']}–{r['End']} ({r['Strand']}) — {r['Method']}")

    if not_found:
        with st.expander(f"⚠️ Not found on reference ({len(not_found)} sequences)"):
            for r in not_found:
                st.text(f"{r['Name']} — not found (likely designed for a newer variant genome)")

    if not found:
        st.error("No sequences could be mapped to the reference genome.")
        return

    # ── Collect positions ──────────────────────────────────────────────────────
    all_positions = set()
    for r in found:
        all_positions.update(range(r["Start"], r["End"] + 1))

    positions_tuple = tuple(sorted(all_positions))
    date_from_str = date_from.strftime("%Y-%m-%d")
    date_to_str = date_to.strftime("%Y-%m-%d")

    # ── Step 1: find mutations at primer positions ─────────────────────────────
    with st.spinner("Finding mutations at primer-probe positions..."):
        mutations_list = fetch_primer_mutations(location_filter, positions_tuple)

    if not mutations_list:
        st.success("✅ No mutations found at primer-probe positions for the selected site.")
        return

    st.caption(f"Found {len(mutations_list)} mutations at primer-probe positions")

    # ── Step 2: fetch time series ──────────────────────────────────────────────
    with st.spinner("Fetching time series data..."):
        time_series = fetch_time_series(
            tuple(sorted(mutations_list)),
            location_filter,
            date_from_str,
            date_to_str
        )

    if not time_series:
        st.warning("No time series data available for the selected period.")
        return

    # ── Process data ───────────────────────────────────────────────────────────
    position_data = process_time_series(time_series, found)
    dominant_letters = get_dominant_letters(time_series, reference, found)

    if not position_data:
        st.success("✅ No mutations found in primer-probe regions for the selected period.")
        return

    # ── Heatmap ────────────────────────────────────────────────────────────────
    st.subheader("Mutation prevalence over time")
    st.caption(f"Site: {selected_location} | {date_from_str} to {date_to_str}")
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

    summary_df = build_summary_table(position_data, found, threshold, dominant_letters)
    st.dataframe(summary_df, use_container_width=True, hide_index=True)

    # ── GenSpectrum links ──────────────────────────────────────────────────────
    st.subheader("GenSpectrum links")
    if selected_location == "All sites":
        st.info("💡 Links below have no location filter. Select a specific site for location-filtered links.")

    for r in found:
        with st.expander(
            f"🧬 {r['Name']} — positions {r['Start']}–{r['End']} ({r['Strand']} strand)"
        ):
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("**Link A — reference sequence prevalence**")
                st.markdown(f"[Open in GenSpectrum ↗]({r['Link A']})")
                st.caption("Low % = mutation dominant at these positions")
            with col2:
                st.markdown("**Link B — any mutation prevalence**")
                st.markdown(f"[Open in GenSpectrum ↗]({r['Link B']})")
                st.caption("High % = mutation present at these positions")