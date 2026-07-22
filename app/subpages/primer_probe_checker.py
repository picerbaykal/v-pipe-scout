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
import streamlit.components.v1 as components
import os

logger = logging.getLogger(__name__)

# ── Helper functions ───────────────────────────────────────────────────────────

@st.cache_data
def load_reference():
    ref_path = os.path.join(os.path.dirname(__file__), "..", "data", "NC_045512.2.fasta")
    if not os.path.exists(ref_path):
        return None, f"Reference not found at {ref_path}"
    record = next(SeqIO.parse(ref_path, "fasta"))
    return str(record.seq).upper(), None

def get_piece_type(name):
    clean = name.split("::")[0].upper()
    if "-P-" in clean or "-P_" in clean or clean.endswith("-P") or clean.endswith("_P"):
        return "probe"
    elif clean.endswith("-F") or clean.endswith("_F") or "_LEFT" in clean:
        return "forward"
    elif clean.endswith("-R") or clean.endswith("_R") or "_RIGHT" in clean:
        return "reverse"
    return "unknown"

def get_3prime_positions(start, end, strand, n_bases=5):
    """
    Return genome positions corresponding to the 3' end of a primer.
    For forward primers (+): 3' end = last n_bases (highest positions)
    For reverse primers (-): 3' end = first n_bases (lowest positions)
    """
    if strand == "+":
        return set(range(end - n_bases + 1, end + 1))
    else:
        return set(range(start, start + n_bases))

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


def build_primer_heatmap_data(time_series, results, dominant_letters):
    """
    Build heatmap data grouped by primer/probe piece.
    Returns dict: {primer_name: {month: {status, text, proportion}}}

    status: 'clean' | 'match' | 'mismatch'
    text: '—' | 'T=T ✓' | 'T≠C 88%'
    proportion: float 0-1 (for mismatch color intensity)
    """
    heatmap_data = {}

    for r in results:
        if r["Start"] is None:
            continue
        name = r["Name"]
        piece_type = get_piece_type(name)
        primer_letters = r.get("PrimerLetters", {})

        # collect all months from time_series
        all_months = set()
        for mut, proportions in time_series.items():
            all_months.update(proportions.keys())
        all_months = sorted(all_months)

        heatmap_data[name] = {
            "piece_type": piece_type,
            "months": {}
        }

        for month in all_months:
            # find worst mismatch in this primer for this month
            worst_proportion = 0.0
            worst_text = None
            is_match = False
            has_any_mutation = False

            for pos, primer_letter in primer_letters.items():
                # find mutations at this position for this month
                for mut, proportions in time_series.items():
                    try:
                        mut_pos = int(re.search(r'\d+', mut).group())
                    except Exception:
                        continue
                    if mut_pos != pos:
                        continue

                    val = proportions.get(month, {})
                    frac = val.get("proportion", 0.0) if isinstance(val, dict) else val
                    if frac <= 0:
                        continue

                    has_any_mutation = True
                    dominant = dominant_letters.get(pos)
                    ref_letter = mut[0]
                    alt_letter = mut[-1]

                    if dominant and primer_letter == dominant:
                        # primer matches dominant circulating letter
                        is_match = True
                        if frac > worst_proportion:
                            worst_proportion = frac
                            worst_text = f"{primer_letter}={dominant} ✓"
                    else:
                        # primer mismatches dominant circulating letter
                        if frac > worst_proportion:
                            worst_proportion = frac
                            worst_text = f"{primer_letter}≠{alt_letter} {frac:.0%}"

            if not has_any_mutation:
                heatmap_data[name]["months"][month] = {
                    "status": "clean",
                    "text": "—",
                    "proportion": 0.0
                }
            elif is_match and worst_text and "✓" in worst_text:
                heatmap_data[name]["months"][month] = {
                    "status": "match",
                    "text": worst_text,
                    "proportion": worst_proportion
                }
            else:
                heatmap_data[name]["months"][month] = {
                    "status": "mismatch",
                    "text": worst_text or "≠",
                    "proportion": worst_proportion
                }

    return heatmap_data

# ── Heatmap ────────────────────────────────────────────────────────────────────
def build_html_heatmap(heatmap_data):
    """
    Build HTML heatmap — rows = primer/probe pieces, columns = months.
    Color encodes mismatch severity. Green = primer matches virus. White = no mutation.
    """
    if not heatmap_data:
        return "<p>No data available.</p>"

    # collect all months that have any data
    all_months = set()
    for name, info in heatmap_data.items():
        for month, cell in info["months"].items():
            if cell["status"] != "clean":
                all_months.add(month)
    all_months = sorted(all_months)

    if not all_months:
        return "<p>No mutations found in primer-probe regions.</p>"

    def mismatch_color(proportion):
        if proportion < 0.1:
            return "background:#fff8f0;color:#bf6000;"
        elif proportion < 0.5:
            return "background:#fac775;color:#412402;"
        elif proportion < 0.8:
            return "background:#d85a30;color:#ffffff;"
        else:
            return "background:#a32d2d;color:#ffffff;"

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
.pp-table{border-collapse:collapse;font-size:13px;font-family:sans-serif;width:100%;}
.pp-table th{text-align:center;padding:8px 8px;color:#666;font-weight:400;border-bottom:1px solid #e5e5e5;white-space:nowrap;}
.pp-table th.name-col{text-align:left;width:auto;}
.pp-table td{padding:3px 5px;border-bottom:1px solid #f0f0f0;}
.pp-table td.name-col{padding:4px 16px 4px 0;white-space:nowrap;}
.pp-cell{border-radius:4px;height:32px;display:flex;align-items:center;justify-content:center;font-size:11px;font-weight:500;min-width:52px;cursor:default;position:relative;}
.pp-cell .tooltip{display:none;position:absolute;bottom:calc(100% + 8px);left:50%;transform:translateX(-50%);background:#fff;border:1px solid #ddd;border-radius:6px;padding:10px 14px;min-width:180px;z-index:100;box-shadow:0 4px 12px rgba(0,0,0,0.1);text-align:left;pointer-events:none;}
.pp-cell:hover .tooltip{display:block;}
.pp-cell .tt-title{font-size:12px;font-weight:500;color:#333;margin-bottom:4px;}
.pp-cell .tt-row{font-size:12px;color:#555;margin-bottom:2px;}
.pp-legend{display:flex;gap:16px;margin-bottom:12px;flex-wrap:wrap;align-items:center;}
.pp-legend-item{display:flex;align-items:center;gap:6px;font-size:12px;color:#666;}
</style>

<div class="pp-legend">
  <div class="pp-legend-item"><span style="background:#fff;border:0.5px solid #ddd;width:16px;height:16px;display:inline-block;border-radius:3px;"></span> No mutation</div>
  <div class="pp-legend-item"><span style="background:#d4edda;width:16px;height:16px;display:inline-block;border-radius:3px;"></span> Primer matches virus</div>
  <div class="pp-legend-item"><span style="width:48px;height:16px;display:inline-block;border-radius:3px;background:linear-gradient(to right,#fff8f0,#a32d2d);"></span> Mismatch (low → high %)</div>
</div>

<div style="display:flex;gap:16px;align-items:flex-start;">
<div style="flex:1;overflow-x:auto;">
<table class="pp-table">
<thead><tr><th class="name-col">Primer/Probe</th>
"""

    for m in all_months:
        html += f"<th>{m}</th>"
    html += "</tr></thead><tbody>"

    for name, info in heatmap_data.items():
        piece_type = info["piece_type"]
        months = info["months"]

        html += f'<tr><td class="name-col">{piece_badge(piece_type)} <span style="color:#333;margin-left:6px;">{name}</span></td>'

        for month in all_months:
            cell = months.get(month, {"status": "clean", "text": "—", "proportion": 0.0})
            status = cell["status"]
            text = cell["text"]
            proportion = cell["proportion"]

            if status == "clean":
                style = "background:#ffffff;color:#999;border:0.5px solid #eee;"
                tooltip_html = f'<div class="tooltip"><div class="tt-title">{name}</div><div class="tt-row">{month}: no mutations detected</div></div>'
            elif status == "match":
                style = "background:#d4edda;color:#155724;"
                tooltip_html = f'<div class="tooltip"><div class="tt-title">{name}</div><div class="tt-row">{month}: {text}</div><div class="tt-row" style="color:#888;font-size:11px;">Primer matches circulating virus</div></div>'
            else:
                style = mismatch_color(proportion)
                tooltip_html = f'<div class="tooltip"><div class="tt-title">{name}</div><div class="tt-row">{month}: {text}</div><div class="tt-row" style="color:#888;font-size:11px;">Primer mismatches circulating virus</div></div>'

            html += f'<td><div class="pp-cell" style="{style}">{text}{tooltip_html}</div></td>'

        html += "</tr>"

    html += "</tbody></table></div>"

    # color scale bar
    html += """
<div style="display:flex;flex-direction:column;align-items:center;gap:4px;padding-top:44px;">
  <span style="font-size:11px;color:#888;">100%</span>
  <div style="width:14px;height:160px;border-radius:8px;background:linear-gradient(to bottom,#a32d2d,#d85a30,#fac775,#fff8f0,#ffffff);border:0.5px solid #ddd;"></div>
  <span style="font-size:11px;color:#888;">0%</span>
  <span style="font-size:10px;color:#aaa;text-align:center;margin-top:4px;max-width:44px;">mismatch %</span>
</div>
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
            "monthly_proportions": {},
            "three_prime_positions": set(),
            "three_prime_mismatch": False,
        }

        # check each position in this primer/probe against dominant circulating letter
        for pos, primer_letter in primer_letters.items():
            dominant = dominant_letters.get(pos)
            if dominant and dominant != primer_letter:
                primer_summary[name]["primer_match"] = False
                primer_summary[name]["mismatch_positions"].append(
                    f"{pos}({primer_letter}→{dominant})"
                )
                # check if mismatch is at 3' end
                if pos in primer_summary[name]["three_prime_positions"]:
                    primer_summary[name]["three_prime_mismatch"] = True

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

        # status
        if not summary["primer_match"]:
            if worst_proportion >= 0.5:
                status = '<span style="background:#fcebeb;color:#a32d2d;font-size:11px;padding:2px 8px;border-radius:4px;">🚨 Critical</span>'
            else:
                status = '<span style="background:#faeeda;color:#854f0b;font-size:11px;padding:2px 8px;border-radius:4px;">⚠️ Warning</span>'
        else:
            status = '<span style="background:#eaf3de;color:#3b6d11;font-size:11px;padding:2px 8px;border-radius:4px;">✅ Good</span>'

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

def build_summary_html(position_data, found, time_series, dominant_letters):
    """
    Build custom HTML summary table with sparkline trend charts on hover.
    One row per primer/probe piece.
    """

    # ── build primer summary data ─────────────────────────────────────────────
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
            "primer_match": True,
            "mismatch_positions": [],
            "monthly_proportions": {},  # {month: proportion} for sparkline
        }

        # check primer vs dominant circulating letter
        for pos, primer_letter in primer_letters.items():
            dominant = dominant_letters.get(pos)
            if dominant and dominant != primer_letter:
                primer_summary[name]["primer_match"] = False
                primer_summary[name]["mismatch_positions"].append(
                    f"{pos}({primer_letter}→{dominant})"
                )

    # update with mutation data — scan time_series for mutations in each primer region
    for r in found:
        if r["Start"] is None:
            continue
        name = r["Name"]
        strand = r["Strand"]
        start = r["Start"]
        end = r["End"]
        primer_summary[name]["three_prime_positions"] = get_3prime_positions(
            start, end, strand
        )

        for mut, proportions in time_series.items():
            try:
                pos = int(re.search(r'\d+', mut).group())
            except Exception:
                continue

            if not (start <= pos <= end):
                continue

            # find most recent non-zero proportion
            for month in sorted(proportions.keys(), reverse=True):
                val = proportions[month]
                frac = val.get("proportion", 0.0) if isinstance(val, dict) else val
                if frac > 0:
                    if frac > primer_summary[name]["worst_proportion"]:
                        primer_summary[name]["worst_proportion"] = frac
                        primer_summary[name]["worst_mutation"] = mut
                        cov = val.get("coverage", 0) if isinstance(val, dict) else 0
                        primer_summary[name]["worst_coverage"] = cov
                    break

    # build monthly proportions for sparkline
    for name, summary in primer_summary.items():
        worst_mut = summary["worst_mutation"]
        if worst_mut and worst_mut in time_series:
            for month, data in time_series[worst_mut].items():
                frac = data.get("proportion", 0.0) if isinstance(data, dict) else data
                if frac > 0:
                    primer_summary[name]["monthly_proportions"][month] = frac

    # build monthly proportions for sparkline — use worst mutation time series
    for name, summary in primer_summary.items():
        worst_mut = summary["worst_mutation"]
        if worst_mut and worst_mut in time_series:
            for month, data in time_series[worst_mut].items():
                frac = data["proportion"] if isinstance(data, dict) else data
                if frac > 0:
                    primer_summary[name]["monthly_proportions"][month] = frac

    # ── helper functions ───────────────────────────────────────────────────────

    def piece_badge(piece_type):
        if piece_type == "probe":
            return '<span style="background:#faeeda;color:#854f0b;font-size:11px;padding:2px 7px;border-radius:4px;font-weight:500;">probe</span>'
        elif piece_type == "forward":
            return '<span style="background:#e6f1fb;color:#185fa5;font-size:11px;padding:2px 7px;border-radius:4px;font-weight:500;">forward</span>'
        elif piece_type == "reverse":
            return '<span style="background:#e1f5ee;color:#0f6e56;font-size:11px;padding:2px 7px;border-radius:4px;font-weight:500;">reverse</span>'
        return ""

    def build_sparkline(monthly_proportions, width=60, height=24):
        """Build a small inline SVG sparkline."""
        if not monthly_proportions:
            return '<svg width="60" height="24" viewBox="0 0 60 24"><line x1="5" y1="12" x2="55" y2="12" stroke="#ddd" stroke-width="1" stroke-dasharray="3,3"/></svg>'

        months = sorted(monthly_proportions.keys())
        values = [monthly_proportions[m] for m in months]
        n = len(values)

        if n == 1:
            y = height - values[0] * height
            return f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}"><circle cx="{width//2}" cy="{y:.1f}" r="3" fill="#a32d2d"/></svg>'

        # scale x and y
        xs = [5 + (i / (n - 1)) * (width - 10) for i in range(n)]
        ys = [height - 4 - v * (height - 8) for v in values]

        points = " ".join(f"{x:.1f},{y:.1f}" for x, y in zip(xs, ys))
        last_x, last_y = xs[-1], ys[-1]

        # color based on last value
        color = "#a32d2d" if values[-1] > 0.5 else "#d85a30" if values[-1] > 0.1 else "#639922"

        return f'''<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<polyline points="{points}" fill="none" stroke="{color}" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
<circle cx="{last_x:.1f}" cy="{last_y:.1f}" r="2.5" fill="{color}"/>
</svg>'''

    def build_popover(name, monthly_proportions):
        """Build hover popover with full trend chart."""
        if not monthly_proportions:
            return ""

        months = sorted(monthly_proportions.keys())
        values = [monthly_proportions[m] for m in months]
        n = len(values)

        w = 200
        pad_top = 16
        pad_bottom = 12
        pad_x = 12
        chart_h = 80

        if n == 1:
            xs = [pad_x, w - pad_x]
            y = pad_top + (1 - values[0]) * (chart_h - pad_top - pad_bottom)
            ys = [y, y]
        else:
            xs = [pad_x + (i / (n - 1)) * (w - 2 * pad_x) for i in range(n)]
            ys = [pad_top + (1 - v) * (chart_h - pad_top - pad_bottom) for v in values]

        points = " ".join(f"{x:.1f},{y:.1f}" for x, y in zip(xs, ys))
        color = "#a32d2d" if values[-1] > 0.5 else "#d85a30" if values[-1] > 0.1 else "#639922"

        dots = "".join(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3" fill="{color}"/>'
            for x, y in zip(xs, ys)
        )

        mid_y = (pad_top + chart_h - pad_bottom) / 2

        labels = "".join(
            f'<text x="{x:.1f}" y="{chart_h + 12}" text-anchor="middle" font-size="9" fill="#999">{m[5:]}</text>'
            f'<text x="{x:.1f}" y="{chart_h + 22}" text-anchor="middle" font-size="9" fill="#666">{v:.0%}</text>'
            for x, m, v in zip(xs, months, values)
        )

        total_h = chart_h + 28

        return f'''<div class="sparkpop">
    <p style="font-size:11px;font-weight:500;color:#333;margin:0 0 6px;">{name} — trend</p>
    <svg width="{w}" height="{total_h}" viewBox="0 0 {w} {total_h}">
      <line x1="0" y1="{pad_top}" x2="{w}" y2="{pad_top}" stroke="#eee" stroke-width="0.5"/>
      <line x1="0" y1="{mid_y:.1f}" x2="{w}" y2="{mid_y:.1f}" stroke="#eee" stroke-width="0.5" stroke-dasharray="3,3"/>
      <line x1="0" y1="{chart_h - pad_bottom}" x2="{w}" y2="{chart_h - pad_bottom}" stroke="#eee" stroke-width="0.5"/>
      <text x="2" y="{pad_top + 3}" font-size="8" fill="#bbb">100%</text>
      <text x="2" y="{mid_y + 3:.1f}" font-size="8" fill="#bbb">50%</text>
      <text x="2" y="{chart_h - pad_bottom + 3}" font-size="8" fill="#bbb">0%</text>
      <polyline points="{points}" fill="none" stroke="{color}" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
      {dots}
      {labels}
    </svg>
    </div>'''

    # ── build HTML ─────────────────────────────────────────────────────────────

    html = """
<style>
.sum-table{border-collapse:collapse;font-size:13px;font-family:sans-serif;width:100%;}
.sum-table th{text-align:left;padding:8px 12px 8px 0;color:#666;font-weight:400;border-bottom:1px solid #e5e5e5;white-space:nowrap;}
.sum-table td{padding:6px 12px 6px 0;border-bottom:1px solid #f0f0f0;vertical-align:middle;}
.spark-wrap{position:relative;display:inline-block;cursor:default;}
.sparkpop{display:none;position:absolute;top:calc(100% + 8px);left:50%;transform:translateX(-50%);background:#fff;border:1px solid #ddd;border-radius:8px;padding:10px 12px;z-index:100;box-shadow:0 4px 12px rgba(0,0,0,0.1);pointer-events:none;white-space:nowrap;}
.spark-wrap:hover .sparkpop{display:block;}
</style>
<table class="sum-table">
<thead>
<tr>
  <th>Primer/Probe</th>
  <th>Type</th>
  <th>Worst mutation</th>
  <th>Recent %</th>
  <th>Coverage</th>
  <th>Primer vs virus</th>
  <th>Trend</th>
  <th>Status</th>
</tr>
</thead>
<tbody>
"""

    for name, summary in primer_summary.items():
        worst_proportion = summary["worst_proportion"]
        worst_mutation = summary["worst_mutation"]
        monthly = summary["monthly_proportions"]

        if not summary["primer_match"]:
            mismatches = ", ".join(summary["mismatch_positions"])
            three_prime_flag = ""
            if summary.get("three_prime_mismatch"):
                three_prime_flag = ' <span style="background:#fcebeb;color:#a32d2d;font-size:10px;padding:1px 6px;border-radius:3px;font-weight:500;">3′ end ⚠️</span>'
            primer_vs = f'<span style="color:#a32d2d;">❌ Mismatch at: {mismatches}</span>{three_prime_flag}'
        else:
            primer_vs = '<span style="color:#3b6d11;">✅ Matches circulating virus</span>'

        # status
        if not summary["primer_match"]:
            if worst_proportion >= 0.5:
                status = '<span style="background:#fcebeb;color:#a32d2d;font-size:11px;padding:2px 8px;border-radius:4px;">🚨 Critical</span>'
            else:
                status = '<span style="background:#faeeda;color:#854f0b;font-size:11px;padding:2px 8px;border-radius:4px;">⚠️ Warning</span>'
        else:
            status = '<span style="background:#eaf3de;color:#3b6d11;font-size:11px;padding:2px 8px;border-radius:4px;">✅ Good</span>'

        # coverage
        cov = summary["worst_coverage"]
        if cov == 0:
            cov_display = '<span style="color:#999;">N/A</span>'
        elif cov < 1000:
            cov_display = f'<span style="color:#854f0b;">⚠️ {cov:,}</span>'
        else:
            cov_display = f'{cov:,}'

        # sparkline + popover
        sparkline = build_sparkline(monthly)
        popover = build_popover(name, monthly)

        if monthly:
            trend_cell = f'<div class="spark-wrap">{sparkline}{popover}</div>'
        else:
            trend_cell = '<span style="color:#ccc;">—</span>'

        html += f"""<tr>
  <td>{piece_badge(summary["piece_type"])} <span style="margin-left:6px;color:#333;">{name}</span></td>
  <td style="color:#666;">{summary["piece_type"]}</td>
  <td style="color:#333;">{worst_mutation or "none"}</td>
  <td style="color:#333;">{"—" if worst_proportion == 0 else f"{worst_proportion:.1%}"}</td>
  <td style="color:#333;">{cov_display}</td>
  <td>{primer_vs}</td>
  <td>{trend_cell}</td>
  <td>{status}</td>
</tr>
"""

    html += "</tbody></table>"
    return html



# ── Main app function ──────────────────────────────────────────────────────────

def app():
    st.title("🔬 Primer & Probe Checker")
    st.markdown(
        """
        Upload primer/probe sequences.
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

    # ── Load reference ─────────────────────────────────────────────────────────
    reference, ref_error = load_reference()
    if ref_error:
        st.error(ref_error)
        return
    st.caption("Reference: NC_045512.2 (Wuhan, 29,903 bp) — loaded automatically")

    # ── Primer/probe FASTA ─────────────────────────────────────────────
    st.subheader("Upload primer/probe FASTA")
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
    heatmap_data = build_primer_heatmap_data(time_series, found, dominant_letters)
    st.subheader("Mutation prevalence over time")
    st.caption(f"Site: {selected_location} | {date_from_str} to {date_to_str}")
    html = build_html_heatmap(heatmap_data)
    st.html(html)

    # ── Summary table ──────────────────────────────────────────────────────────
    st.subheader("Summary")

    first_mut = list(time_series.keys())[0] if time_series else None

    with st.expander("Summary", expanded=True):
        summary_html = build_summary_html(
            position_data, found, time_series, dominant_letters
        )
        components.html(summary_html, height=400, scrolling=True)

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