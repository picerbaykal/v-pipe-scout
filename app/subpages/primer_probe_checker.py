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

def get_assay_name(name):
    """
    Group pieces into assays (ddPCR) or amplicons (ARTIC) by dropping the role suffix:
    N1-F / N1-P / N1-R -> N1;  SARS-CoV-2_5_LEFT / _RIGHT -> SARS-CoV-2_5.
    Names that don't follow a known pattern form their own group.
    """
    clean = name.split("::")[0]
    clean = re.sub(r"[_-]alt\d*$", "", clean, flags=re.IGNORECASE)
    m = re.match(
        r"^(.*?)[_-](F|R|P|FW|FWD|REV|LEFT|RIGHT|PROBE)$", clean, flags=re.IGNORECASE
    )
    return m.group(1) if m else clean


STATUS_RANK = {"OK": 0, "Watch": 1, "Compromised": 2}
STATUS_STYLE = {
    "OK": ("#e9f5e1", "#3b6d11"),
    "Watch": ("#fdf0dc", "#854f0b"),
    "Compromised": ("#fbe4e4", "#a32d2d"),
}


def classify_piece(summary, mismatch_threshold, watch_threshold, min_coverage, three_prime_n):
    """
    Turn one primer/probe summary into a status + plain-language reason.

    Compromised: probe differs from the circulating virus; or a primer mismatch in the
                 3'-end zone; or 2+ primer mismatches.
    Watch:       1 primer mismatch outside the 3'-end zone; or a minority mutation
                 between the watch and flag thresholds.
    OK:          matches the circulating virus.
    """
    mism = summary["mismatch_positions"]
    is_probe = summary["piece_type"] == "probe"

    if mism:
        listed = ", ".join(mism)
        if is_probe:
            status, reason = "Compromised", f"Probe differs from circulating virus at {listed}"
        elif summary["three_prime_mismatch"]:
            status, reason = "Compromised", f"Mismatch in last {three_prime_n} bases (3′ end): {listed}"
        elif len(mism) >= 2:
            status, reason = "Compromised", f"{len(mism)} mismatches: {listed}"
        else:
            status, reason = "Watch", f"1 mismatch, not at 3′ end: {listed}"
    else:
        frac = summary["worst_proportion"]
        if summary["worst_mutation"] and watch_threshold <= frac < mismatch_threshold:
            status = "Watch"
            reason = (
                f"{summary['worst_mutation']} at {frac:.0%} "
                f"(below the {mismatch_threshold:.0%} flag threshold)"
            )
        elif summary["worst_mutation"]:
            status, reason = "OK", "Matches circulating virus"
        else:
            status, reason = "OK", "No mutations detected"

    cov = summary["worst_coverage"]
    if summary["worst_mutation"] and 0 < cov < min_coverage:
        reason += f" — low coverage ({cov:,} reads), uncertain"
    return status, reason


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

def get_mismatch_position_index(genome_pos, start, end, strand):
    """
    Calculate 1-based position index of a mismatch within the primer sequence.
    For forward primers: index from 5' end (left to right)
    For reverse primers: index from 5' end (right to left on genome)
    Returns (index, length) e.g. (11, 20)
    """
    length = end - start + 1
    if strand == "+":
        index = genome_pos - start + 1
    else:
        index = end - genome_pos + 1
    return index, length

def calculate_gc_content(primer_seq):
    """GC content as a fraction (0-1)."""
    seq = primer_seq.upper()
    if not seq:
        return 0.0
    gc = sum(1 for base in seq if base in "GC")
    return gc / len(seq)

def extract_positions_from_header(header):
    match = re.search(r"NC_045512\.2:(\d+)-(\d+)", header.strip())
    if match:
        start = int(match.group(1))
        end = int(match.group(2))
        # detect strand from primer name
        if "_RIGHT" in header.upper():
            strand = "-"
        else:
            strand = "+"
        return start, end, strand
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
def fetch_primer_mutations(location, positions_tuple, min_proportion=0.01,
                           date_from_str=None, date_to_str=None):
    """
    Step 1: fetch which mutations exist at primer/probe positions.
    location=None means all sites. Only mutations whose share, pooled over the
    selected site and dates, is at least `min_proportion` are returned
    (LAPIS would otherwise apply its own default of 5%).
    Returns a list of mutations, or None if LAPIS could not be queried.
    """
    import logging
    logger = logging.getLogger(__name__)
    async def _fetch():
        params = {"limit": 100000, "minProportion": min_proportion}
        if location:
            params["locationName"] = location
        if date_from_str:
            params["samplingDateFrom"] = date_from_str
        if date_to_str:
            params["samplingDateTo"] = date_to_str
        async with aiohttp.ClientSession() as session:
            async with session.get(
                f"{get_wiseloculus_url()}/sample/nucleotideMutations",
                params=params,
                headers={"accept": "application/json"}
            ) as response:
                if response.status != 200:
                    logger.error(f"fetch_primer_mutations HTTP {response.status}: {await response.text()}")
                    return None
                data = await response.json()
                positions = set(positions_tuple)
                return [
                    f"{entry['mutationFrom']}{entry['position']}{entry['mutationTo']}"
                    for entry in data.get("data", [])
                    if entry.get("position") in positions
                ]
    try:
        return asyncio.run(_fetch())
    except Exception as e:
        logger.error(f"fetch_primer_mutations error: {e}")
        return None

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

    # add main: prefix to mutations
    mutations_with_prefix = [f"main:{m}" for m in mutations_list]
    async def _fetch():
        async with aiohttp.ClientSession(
            timeout=aiohttp.ClientTimeout(total=120)
        ) as session:
            async with session.post(
                f"{get_wiseloculus_url()}/component/nucleotideMutationsOverTime",
                json={
                    "filters": filters,
                    "dateRanges": date_ranges,
                    "includeMutations": mutations_with_prefix,
                    "dateField": "samplingDate"
                },
                headers={"content-type": "application/json"}
            ) as response:
                if response.status != 200:
                    return {}
                data = (await response.json())["data"]
                mutations = data["mutations"]
                date_labels = [r["dateFrom"][:7] for r in data["dateRanges"]]
                result = {}
                for i, mut in enumerate(mutations):
                    # strip main: prefix from mutation key
                    clean_mut = mut.replace("main:", "")
                    result[clean_mut] = {}
                    for j, label in enumerate(date_labels):
                        cell = data["data"][i][j]
                        if cell["coverage"] > 0:
                            result[clean_mut][label] = {
                                "proportion": cell["count"] / cell["coverage"],
                                "coverage": cell["coverage"],
                                "count": cell["count"],
                            }
                        else:
                            result[clean_mut][label] = {
                                "proportion": 0.0,
                                "coverage": 0,
                                "count": 0,
                            }
                return result

    try:
        return asyncio.run(_fetch())
    except Exception as e:
        st.error(f"fetch_time_series exception: {e}")
        import traceback
        st.code(traceback.format_exc())
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


def get_dominant_letters(time_series, reference, results, threshold=0.5):
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
        if recent_frac > threshold:
            dominant[pos] = alt_letter  # mutation is dominant
        else:
            dominant[pos] = ref_letter  # reference is dominant

    return dominant


def build_primer_heatmap_data(time_series, results, dominant_letters, threshold=0.5,
                              min_coverage=1000):
    """
    Build heatmap data grouped by primer/probe piece, judging EACH MONTH on its own.

    For every position in the primer and every month, the circulating letter is the
    mutation if its share that month is above `threshold`, otherwise the reference letter.
    A primer letter that differs from the circulating letter is a mismatch.

    Returns dict: {primer_name: {"piece_type", "months": {month: cell}}}
    cell status: 'clean' | 'match' | 'mismatch'
    cell text:   '—' | 'T ✓ · C 69%' (mutation below threshold) | 'C=C ✓' | 'T≠C 88%'
    proportion:  share of the virus that does NOT match the primer (mismatch colour)
    """
    heatmap_data = {}

    # mutations per genome position, parsed once
    muts_by_pos = {}
    for mut in time_series:
        try:
            mut_pos = int(re.search(r"\d+", mut).group())
        except Exception:
            continue
        muts_by_pos.setdefault(mut_pos, []).append(mut)

    all_months = sorted({m for proportions in time_series.values() for m in proportions})

    for r in results:
        if r["Start"] is None:
            continue
        name = r["Name"]
        primer_letters = r.get("PrimerLetters", {})
        heatmap_data[name] = {"piece_type": get_piece_type(name), "months": {}}

        for month in all_months:
            mismatch_count = 0
            worst_share = 0.0
            worst = None        # (text, pos, circulating letter)
            best_match = None   # (share of mutation, text, note)

            for pos, primer_letter in primer_letters.items():
                # shares of each mutated letter at this position this month
                shares = {}
                counts = {}
                coverage = 0
                ref_letter = None
                for mut in muts_by_pos.get(pos, []):
                    val = time_series[mut].get(month, {})
                    frac = val.get("proportion", 0.0) if isinstance(val, dict) else val
                    if frac > 0:
                        shares[mut[-1]] = shares.get(mut[-1], 0.0) + frac
                        if isinstance(val, dict):
                            counts[mut[-1]] = counts.get(mut[-1], 0) + val.get("count", 0)
                            coverage = max(coverage, val.get("coverage", 0))
                        ref_letter = mut[0]
                if not shares:
                    continue

                ref_share = max(0.0, 1.0 - sum(shares.values()))
                top_alt = max(shares, key=shares.get)
                top_frac = shares[top_alt]
                reads = (top_alt, counts.get(top_alt, 0), coverage, top_frac)
                circulating = top_alt if top_frac > threshold else ref_letter
                primer_share = ref_share if primer_letter == ref_letter else shares.get(primer_letter, 0.0)

                if primer_letter != circulating:
                    mismatch_count += 1
                    share = 1.0 - primer_share
                    if share > worst_share or worst is None:
                        worst_share = share
                        worst = (f"{primer_letter}≠{circulating} {share:.0%}", pos, circulating, reads)
                else:
                    if circulating != ref_letter:
                        # primer already carries the mutation (e.g. redesigned primer)
                        text, note = f"{primer_letter}={circulating} ✓", ""
                    else:
                        # mutation present, but below the threshold
                        # (under 1% it is noise in the cell label; the hover still shows it)
                        pct = f"{top_frac:.0%}" if top_frac >= 0.01 else f"{top_frac:.1%}"
                        text = (
                            f"{primer_letter} ✓ · {top_alt} {pct}"
                            if top_frac >= 0.01 else f"{primer_letter} ✓"
                        )
                        note = f"{pos}: {top_alt} at {pct}, below the {threshold:.0%} threshold"
                    if best_match is None or top_frac > best_match[0]:
                        best_match = (top_frac, text, note, reads)

            if worst is not None:
                heatmap_data[name]["months"][month] = {
                    "status": "mismatch",
                    "text": worst[0],
                    "proportion": worst_share,
                    "pos": worst[1],
                    "alt": worst[2],
                    "mismatch_count": mismatch_count,
                    "reads": worst[3],
                    "low_cov": 0 < worst[3][2] < min_coverage,
                }
            elif best_match is not None:
                heatmap_data[name]["months"][month] = {
                    "status": "match",
                    "text": best_match[1],
                    "proportion": 0.0,
                    "note": best_match[2],
                    "reads": best_match[3],
                    "low_cov": 0 < best_match[3][2] < min_coverage,
                }
            else:
                heatmap_data[name]["months"][month] = {
                    "status": "clean",
                    "text": "—",
                    "proportion": 0.0,
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
.pp-table{border-collapse:collapse;font-size:15px;font-family:sans-serif;width:100%;}
.pp-table th{text-align:center;padding:10px 14px;color:#666;font-weight:400;border-bottom:1px solid #e5e5e5;white-space:nowrap;font-size:14px;}
.pp-table th.name-col{text-align:left;width:auto;}
.pp-table td{padding:3px 5px;border-bottom:1px solid #f0f0f0;}
.pp-table td.name-col{padding:4px 16px 4px 0;white-space:nowrap;}
.pp-cell{border-radius:4px;height:40px;display:flex;align-items:center;justify-content:center;font-size:13px;font-weight:500;min-width:80px;cursor:default;position:relative;}
.pp-cell .tooltip{display:none;position:absolute;top:calc(100% + 8px);left:50%;transform:translateX(-50%);background:#fff;border:1px solid #ddd;border-radius:6px;padding:10px 14px;min-width:180px;z-index:100;box-shadow:0 4px 12px rgba(0,0,0,0.1);text-align:left;pointer-events:none;}
.pp-cell .tooltip.tip-up{top:auto;bottom:calc(100% + 8px);}
.pp-table td:nth-last-child(-n+2) .pp-cell .tooltip{left:auto;right:0;transform:none;}
.pp-cell:hover .tooltip{display:block;}
.pp-cell .tt-title{font-size:12px;font-weight:500;color:#333;margin-bottom:4px;}
.pp-cell .tt-row{font-size:12px;color:#555;margin-bottom:2px;}
.pp-legend{display:flex;gap:16px;margin-bottom:12px;flex-wrap:wrap;align-items:center;}
.pp-legend-item{display:flex;align-items:center;gap:6px;font-size:12px;color:#666;}
</style>

<div class="pp-legend">
  <div class="pp-legend-item"><span style="background:#fff;border:0.5px solid #ddd;width:16px;height:16px;display:inline-block;border-radius:3px;"></span> No mutation</div>
  <div class="pp-legend-item"><span style="background:#d4edda;width:16px;height:16px;display:inline-block;border-radius:3px;"></span> Primer matches virus &nbsp;<span style="color:#999;">(“T ✓ · C 69%” = mutation present but below the threshold)</span></div>
  <div class="pp-legend-item"><span style="width:48px;height:16px;display:inline-block;border-radius:3px;background:linear-gradient(to right,#fff8f0,#a32d2d);"></span> Mismatch (low → high %)</div>
  <div class="pp-legend-item"><span style="width:16px;height:16px;display:inline-block;border-radius:3px;outline:2px dashed #8a8a8a;outline-offset:-3px;"></span> Few reads (below minimum coverage)</div>
</div>

<div style="display:flex;gap:16px;align-items:flex-start;">
<div style="flex:1;overflow-x:auto;">
<table class="pp-table">
<thead><tr><th class="name-col">Primer/Probe</th>
"""

    for m in all_months:
        html += f"<th>{m}</th>"
    html += "</tr></thead><tbody>"

    # count visible rows to decide tooltip direction
    visible_names = [
        nm for nm, inf in heatmap_data.items()
        if not all(c["status"] == "clean" for c in inf["months"].values())
    ]
    total_visible = len(visible_names)

    row_idx = -1
    for name, info in heatmap_data.items():
        piece_type = info["piece_type"]
        months = info["months"]

        # skip rows where all cells are clean
        if all(cell["status"] == "clean" for cell in months.values()):
            continue

        row_idx += 1
        # bottom half rows: tooltip appears ABOVE to stay in view
        flip_up = row_idx >= total_visible / 2
        tip_class = "tooltip tip-up" if flip_up else "tooltip tip-down"

        html += f'<tr><td class="name-col">{piece_badge(piece_type)} <span style="color:#333;margin-left:6px;">{name}</span></td>'

        for month in all_months:
            cell = months.get(month, {"status": "clean", "text": "—", "proportion": 0.0})
            status = cell["status"]
            text = cell["text"]
            proportion = cell["proportion"]
            reads_line = ""
            if cell.get("reads") and cell["reads"][2]:
                r_alt, r_count, r_cov, r_frac = cell["reads"]
                reads_line = (
                    f'<div class="tt-row" style="color:#333;">{r_alt} in {r_count:,} of '
                    f'{r_cov:,} reads ({r_frac:.1%})</div>'
                )
                if cell.get("low_cov"):
                    reads_line += (
                        '<div class="tt-row" style="color:#854f0b;font-size:11px;">'
                        '⚠ fewer reads than the minimum coverage setting</div>'
                    )

            if status == "clean":
                style = "background:#ffffff;color:#999;border:0.5px solid #eee;"
                tooltip_html = f'<div class="{tip_class}"><div class="tt-title">{name}</div><div class="tt-row">{month}: no mutations detected</div></div>'
            elif status == "match":
                style = "background:#d4edda;color:#155724;"
                note_line = (
                    f'<div class="tt-row" style="color:#333;">{cell["note"]}</div>'
                    if cell.get("note") else ""
                )
                tooltip_html = f'<div class="{tip_class}"><div class="tt-title">{name}</div><div class="tt-row">{month}: {text}</div>{note_line}{reads_line}<div class="tt-row" style="color:#888;font-size:11px;">Primer matches circulating virus</div></div>'
            else:
                style = mismatch_color(proportion)
                pos_line = ""
                if cell.get("pos") and cell.get("alt"):
                    ref_at_pos = text.split("≠")[0] if "≠" in text else ""
                    pos_line = f'<div class="tt-row" style="color:#333;font-weight:500;">{cell["pos"]}: {ref_at_pos}→{cell["alt"]}</div>'
                extra = ""
                mc = cell.get("mismatch_count", 0)
                if mc > 1:
                    extra = f'<div class="tt-row" style="color:#a32d2d;font-size:11px;">+{mc - 1} more mismatch(es) in this primer</div>'
                tooltip_html = f'<div class="{tip_class}"><div class="tt-title">{name}</div><div class="tt-row">{month}: {text}</div>{pos_line}{reads_line}{extra}<div class="tt-row" style="color:#888;font-size:11px;">Primer mismatches circulating virus</div></div>'

            if cell.get("low_cov"):
                style += "outline:2px dashed #8a8a8a;outline-offset:-3px;"
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

def build_summary_html(
    position_data,
    found,
    time_series,
    dominant_letters,
    three_prime_n=5,
    min_coverage=1000,
    mismatch_threshold=0.5,
    watch_threshold=0.2,
):
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
        strand = r["Strand"]
        start = r["Start"]
        end = r["End"]

        gc_content = calculate_gc_content("".join(primer_letters.values()))

        primer_summary[name] = {
            "piece_type": piece_type,
            "worst_mutation": None,
            "worst_proportion": 0.0,
            "worst_coverage": 0,
            "primer_match": True,
            "mismatch_positions": [],
            "mismatch_index": None,
            "primer_length": end - start + 1,
            "monthly_proportions": {},
            "three_prime_positions": get_3prime_positions(start, end, strand, n_bases=three_prime_n),
            "three_prime_mismatch": False,
            "gc_content": gc_content,
        }

        # check primer vs dominant circulating letter
        for pos, primer_letter in primer_letters.items():
            dominant = dominant_letters.get(pos)
            if dominant and dominant != primer_letter:
                primer_summary[name]["primer_match"] = False
                primer_summary[name]["mismatch_positions"].append(
                    f"{pos}({primer_letter}→{dominant})"
                )
                if primer_summary[name]["mismatch_index"] is None:
                    idx, length = get_mismatch_position_index(pos, start, end, strand)
                    primer_summary[name]["mismatch_index"] = idx
                if pos in primer_summary[name]["three_prime_positions"]:
                    primer_summary[name]["three_prime_mismatch"] = True


    # update with mutation data — scan time_series for mutations in each primer region
    for r in found:
        if r["Start"] is None:
            continue
        name = r["Name"]
        piece_type = get_piece_type(name)
        primer_letters = r.get("PrimerLetters", {})
        strand = r["Strand"]
        start = r["Start"]
        end = r["End"]
        primer_summary[name]["three_prime_positions"] = get_3prime_positions(
            start, end, strand, n_bases=three_prime_n
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

    # build monthly proportions for sparkline — use worst mutation time series
    for name, summary in primer_summary.items():
        worst_mut = summary["worst_mutation"]
        if worst_mut and worst_mut in time_series:
            for month, data in time_series[worst_mut].items():
                frac = data["proportion"] if isinstance(data, dict) else data
                if frac > 0:
                    primer_summary[name]["monthly_proportions"][month] = frac

    # ── status per piece + verdict per assay ──────────────────────────────────
    assay_verdicts = {}
    for name, summary in primer_summary.items():
        status, reason = classify_piece(
            summary, mismatch_threshold, watch_threshold, min_coverage, three_prime_n
        )
        summary["status"] = status
        summary["reason"] = reason
        assay = get_assay_name(name)
        summary["assay"] = assay
        v = assay_verdicts.setdefault(assay, {"status": "OK", "reason": "", "pieces": []})
        v["pieces"].append(name)
        # keep the first piece with the worst status as the assay's reason
        if not v["reason"] or STATUS_RANK[status] > STATUS_RANK[v["status"]]:
            v["status"] = status
            v["reason"] = f"{name}: {reason}"
    for v in assay_verdicts.values():
        if v["status"] == "OK":
            v["reason"] = "All pieces match the circulating virus"

    # show assay header rows only when grouping actually groups something
    grouped = any(len(v["pieces"]) > 1 for v in assay_verdicts.values())

    def status_pill(status):
        bg, fg = STATUS_STYLE[status]
        return (
            f'<span style="background:{bg};color:{fg};font-size:12px;padding:2px 8px;'
            f'border-radius:10px;font-weight:500;white-space:nowrap;">{status}</span>'
        )

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
        color = "#a32d2d" if values[-1] > mismatch_threshold else "#d85a30" if values[-1] > watch_threshold else "#639922"

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
        color = "#a32d2d" if values[-1] > mismatch_threshold else "#d85a30" if values[-1] > watch_threshold else "#639922"

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


    def build_position_bar(summary, name):
        """Build mini position bar showing mismatch location within primer.
        Uses hardcoded colors (not CSS vars) because this renders inside a
        components.html iframe where Streamlit theme vars are unavailable."""
        length = summary.get("primer_length", 20)
        index = summary.get("mismatch_index")
        mismatch_positions = summary.get("mismatch_positions", [])
        # 3' zone = last `three_prime_n` bases, drawn as a share of the bar
        three_prime_end_pct = min(100, round(100 * three_prime_n / max(length, 1)))

        if index is None:
            dot_color = "#3b6d11"
            dot_left = 50
            zone_html = ""
            label = "\u2713"  # checkmark
            label_color = "#3b6d11"
            tooltip_text = "No mutations detected"
            tooltip_detail = ""
        else:
            dot_left = int((index / length) * 100)
            pct = index / length
            if summary.get("piece_type") == "probe":
                # for a probe, any mismatch weakens the signal regardless of position
                dot_color = "#a32d2d"
                pos_desc = "probe \u2014 any mismatch affects the signal"
                zone_html = ""
            elif summary.get("three_prime_mismatch"):
                dot_color = "#a32d2d"
                pos_desc = (
                    "in last " + str(three_prime_n)
                    + " bases (3\u2032 end) \u2014 severe impact \u26a0\ufe0f"
                )
                zone_html = (
                        '<div style="position:absolute;right:0;top:0;width:'
                        + str(three_prime_end_pct)
                        + '%;height:100%;background:#f7c1c1;border-radius:0 5px 5px 0;"></div>'
                )
            elif pct < 0.25:
                dot_color = "#ba7517"
                pos_desc = "near 5\u2032 end \u2014 minor impact"
                zone_html = ""
            else:
                dot_color = "#d85a30"
                pos_desc = "middle of primer \u2014 moderate impact"
                zone_html = ""
            label = str(index) + "/" + str(length)
            label_color = "#555"
            tooltip_text = ", ".join(mismatch_positions)
            tooltip_detail = "Position " + str(index) + "/" + str(length) + " \u2014 " + pos_desc

        uid = "pb_" + str(abs(hash(name + str(index))))

        dot = (
            '<div style="position:absolute;left:' + str(dot_left)
            + '%;top:50%;transform:translate(-50%,-50%);width:11px;height:11px;background:'
            + dot_color + ';border-radius:50%;border:2px solid #fff;box-shadow:0 0 0 1px '
            + dot_color + ';"></div>'
        )

        parts = [
            '<div style="position:relative;display:inline-block;"',
            ' onmouseenter="document.getElementById(\'', uid, '\').style.display=\'block\'"',
            ' onmouseleave="document.getElementById(\'', uid, '\').style.display=\'none\'">',
            '<div style="display:flex;align-items:center;gap:6px;cursor:default;">',
            '<div style="position:relative;width:72px;height:10px;background:#eeeeee;border-radius:5px;border:0.5px solid #ccc;">',
            zone_html, dot,
            '</div>',
            '<span style="font-size:11px;color:', label_color, ';">', label, '</span>',
            '</div>',
            '<div id="', uid, '" style="display:none;position:absolute;top:calc(100% + 6px);left:0;background:#ffffff;border:1px solid #ccc;border-radius:8px;padding:12px 14px;min-width:200px;z-index:9999;box-shadow:0 6px 20px rgba(0,0,0,0.15);white-space:nowrap;">',
            '<p style="font-size:11px;font-weight:500;color:#333;margin:0 0 8px;">', tooltip_text, '</p>',
            '<div style="position:relative;width:100%;height:10px;background:#eeeeee;border-radius:5px;border:0.5px solid #ccc;margin-bottom:4px;">',
            zone_html, dot,
            '</div>',
            '<div style="display:flex;justify-content:space-between;font-size:10px;color:#999;margin-bottom:8px;"><span>5\u2032</span><span>3\u2032</span></div>',
            '<p style="font-size:11px;color:#666;margin:0;">', tooltip_detail, '</p>',
            '</div>',
            '</div>',
        ]
        return "".join(parts)

    # ── build HTML ─────────────────────────────────────────────────────────────

    html = """
<style>
.sum-table{border-collapse:collapse;font-size:15px;font-family:sans-serif;width:100%;}
.sum-table th{text-align:left;padding:12px 14px 12px 0;color:#666;font-weight:400;border-bottom:1px solid #e5e5e5;white-space:nowrap;font-size:14px;}
.sum-table td{padding:12px 14px 12px 0;border-bottom:1px solid #f0f0f0;vertical-align:middle;}
.sum-table td.name-cell{white-space:nowrap;}
.sum-table tr.grp td{background:#fafafa;padding:10px 12px;border-bottom:1px solid #e5e5e5;}
.reason{font-size:11px;color:#777;margin-top:4px;min-width:140px;max-width:220px;white-space:normal;}
.sum-table .sparkpop{left:auto;right:0;transform:none;}
.sum-table tbody tr:nth-last-child(-n+3) .sparkpop{top:auto;bottom:calc(100% + 8px);}
.sum-table tbody tr:nth-last-child(-n+3) div[id^="pb_"]{top:auto !important;bottom:calc(100% + 6px) !important;}
.spark-wrap{position:relative;display:inline-block;cursor:default;}
.sparkpop{display:none;position:absolute;top:calc(100% + 8px);left:50%;transform:translateX(-50%);background:#fff;border:1px solid #ddd;border-radius:8px;padding:10px 12px;z-index:100;box-shadow:0 4px 12px rgba(0,0,0,0.1);pointer-events:none;white-space:nowrap;}
.spark-wrap:hover .sparkpop{display:block;}
</style>
<table class="sum-table">
<thead>
<tr>
  <th>Primer/Probe</th>
  <th>Status</th>
  <th>Worst mutation</th>
  <th>Recent %</th>
  <th>Coverage</th>
  <th>GC%</th>
  <th>Position</th>
  <th>Trend</th>
</tr>
</thead>
<tbody>
"""

    current_assay = None
    for name, summary in primer_summary.items():
        if grouped and summary["assay"] != current_assay:
            current_assay = summary["assay"]
            v = assay_verdicts[current_assay]
            html += (
                "<tr class='grp'><td colspan='8'>"
                f"<span style='font-weight:600;color:#333;margin-right:10px;'>{current_assay}</span>"
                f"{status_pill(v['status'])}"
                f"<span style='font-size:12px;color:#666;margin-left:10px;'>{v['reason']}</span>"
                "</td></tr>"
            )

        worst_proportion = summary["worst_proportion"]
        worst_mutation = summary["worst_mutation"]
        monthly = summary["monthly_proportions"]

        # coverage
        cov = summary["worst_coverage"]
        if cov == 0:
            cov_display = '<span style="color:#999;">N/A</span>'
        elif cov < min_coverage:
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

        proportion_display = "—" if worst_proportion == 0 else f"{worst_proportion:.1%}"

        gc = summary.get("gc_content", 0)
        gc_pct = f"{gc:.0%}"
        # color hint: 40-60% ideal (green), 30-40/60-70 acceptable (amber), else red
        if 0.40 <= gc <= 0.60:
            gc_color = "#3b6d11"
        elif 0.30 <= gc < 0.40 or 0.60 < gc <= 0.70:
            gc_color = "#ba7517"
        else:
            gc_color = "#a32d2d"
        gc_display = f'<span style="color:{gc_color};">{gc_pct}</span>'

        html += (
            "<tr>"
            f"<td class='name-cell'>{piece_badge(summary['piece_type'])} <span style='margin-left:6px;color:#333;'>{name}</span></td>"
            f"<td>{status_pill(summary['status'])}<div class='reason'>{summary['reason']}</div></td>"
            f"<td style='color:#333;'>{worst_mutation or 'none'}</td>"
            f"<td style='color:#333;'>{proportion_display}</td>"
            f"<td style='color:#333;'>{cov_display}</td>"
            f"<td>{gc_display}</td>"
            f"<td>{build_position_bar(summary, name)}</td>"
            f"<td>{trend_cell}</td>"
            "</tr>"
        )

    html += "</tbody></table>"
    return html, assay_verdicts



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

    # Check URL query params for demo auto-load
    # ?demo=cdc  → load bundled CDC N1/N2 example
    # ?fasta=... → load FASTA content passed directly in the URL
    query_params = st.query_params
    demo_content = None

    if query_params.get("fasta"):
        demo_content = query_params["fasta"]
    elif query_params.get("demo") == "cdc":
        demo_path = os.path.join(
            os.path.dirname(__file__), "..", "data", "cdc_n1_n2.fasta"
        )
        if os.path.exists(demo_path):
            with open(demo_path) as f:
                demo_content = f.read()

    # --- FASTA source state (survives reruns, e.g. when the city filter changes) ---
    if "ppc_source" not in st.session_state:
        st.session_state.ppc_source = None  # "example" | "upload" | None
    if "ppc_uploader_counter" not in st.session_state:
        st.session_state.ppc_uploader_counter = 0

    def _read_cdc_example():
        p = os.path.join(os.path.dirname(__file__), "..", "data", "cdc_n1_n2.fasta")
        if os.path.exists(p):
            with open(p) as f:
                return f.read()
        return None

    primer_file = st.file_uploader(
        "Primer/probe FASTA",
        type=["fasta", "fa"],
        help="e.g. CDC N1/N2 diagnostic primers, or ARTIC amplicon primers",
        key=f"ppc_fasta_upload_{st.session_state.ppc_uploader_counter}",
    )

    # A fresh upload takes over from any active example.
    if primer_file is not None:
        st.session_state.ppc_source = "upload"

    # "Try an example" — load bundled CDC N1/N2 and clear any stale upload.
    load_example = st.button("Try an example (CDC N1/N2)")
    if load_example:
        st.session_state.ppc_source = "example"
        st.session_state.ppc_uploader_counter += 1  # recreate the uploader empty
        st.rerun()

    # Decide the source of the FASTA content.
    if demo_content is not None:  # ?demo=cdc / ?fasta=... (URL)
        raw_content = demo_content
    elif st.session_state.ppc_source == "example":
        raw_content = _read_cdc_example()
    elif st.session_state.ppc_source == "upload" and primer_file is not None:
        raw_content = primer_file.read().decode("utf-8")
    else:
        st.info("👆 Upload a primer/probe FASTA file, or click **Try an example** above.")
        return

    if raw_content is None:
        st.error("Could not load the CDC N1/N2 example file.")
        return

    # ── Parse primers ──────────────────────────────────────────────────────────
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
                        "PrimerLetters": {},  # ← empty dict, no position found
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

    # ── Assessment settings (on the page, pinned while scrolling) ─────────────
    st.markdown(
        """<style>
        /* pin the settings bar: sticky goes on Streamlit's wrapper around the container */
        div[data-testid="stLayoutWrapper"]:has(> .st-key-ppc_settings) {
            position: sticky; top: 3.75rem; z-index: 990;
        }
        .st-key-ppc_settings {
            background: #ffffff; padding: 0.5rem 0.75rem 0.25rem;
            border: 1px solid #e6e6e6; border-radius: 0.5rem;
            box-shadow: 0 2px 6px rgba(0,0,0,0.08);
        }
        </style>""",
        unsafe_allow_html=True,
    )
    with st.container(key="ppc_settings"):
        st.markdown(
            "**Assessment settings** "
            "<span style='color:#888;font-size:0.85em;'>— change these to see how the "
            "results shift; the wastewater data is not re-fetched</span>",
            unsafe_allow_html=True,
        )
        s0, s1, s2, s3, s4 = st.columns(5)
        with s0:
            detect_threshold = st.select_slider(
                "Detect mutations from (%)", options=[0.5, 1, 2, 5, 10], value=1,
                key="ppc_detect_threshold",
                help="Mutations are looked up only if their share, pooled over the selected "
                     "site and dates, reaches this level. Below ~1% sequencing noise increases. "
                     "LAPIS's own default would be 5%.",
            ) / 100
        with s1:
            mismatch_threshold = st.slider(
                "Flag mutations above (%)", 10, 90, 50, 5,
                key="ppc_mismatch_threshold",
                help="A mutation above this share of the virus counts against the primer/probe. "
                     "Heatmap: judged month by month. Summary: most recent month.",
            ) / 100
        with s2:
            watch_threshold = st.slider(
                "Watch threshold (%)", 1, 50, 20, 1,
                key="ppc_watch_threshold",
                help="A minority mutation above this share (but below the flag threshold) "
                     "is flagged early.",
            ) / 100
        with s3:
            three_prime_n = st.slider(
                "3′-end zone (bases)", 1, 10, 5, 1,
                key="ppc_three_prime_n",
                help="A primer mismatch within this many bases of the 3′ end counts as critical.",
            )
        with s4:
            min_coverage = st.number_input(
                "Minimum coverage (reads)", min_value=0, value=1000, step=100,
                key="ppc_min_coverage",
                help="Months or mutations with fewer reads are treated as unreliable.",
            )
        if watch_threshold >= mismatch_threshold:
            st.warning("Watch threshold should be below the flag threshold.")

    # ── Step 1: find mutations at primer positions ─────────────────────────────
    with st.spinner("Finding mutations at primer-probe positions..."):
        mutations_list = fetch_primer_mutations(
            location_filter, positions_tuple, detect_threshold, date_from_str, date_to_str
        )

    if mutations_list is None:
        st.error("Could not query mutations from LAPIS. Please try again later.")
        return
    if not mutations_list:
        st.success(
            f"✅ No mutations at primer-probe positions reach {detect_threshold:.1%} "
            "for the selected site and dates."
        )
        return

    st.caption(
        f"Found {len(mutations_list)} mutations at primer-probe positions "
        f"(share ≥ {detect_threshold:.1%}, pooled over the selected site and dates)"
    )

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
    dominant_letters = get_dominant_letters(
        time_series, reference, found, threshold=mismatch_threshold
    )

    if not position_data:
        st.success("✅ No mutations found in primer-probe regions for the selected period.")
        return

    # ── Heatmap ────────────────────────────────────────────────────────────────
    heatmap_data = build_primer_heatmap_data(
        time_series, found, dominant_letters, threshold=mismatch_threshold,
        min_coverage=min_coverage,
    )
    st.subheader("Mutation prevalence over time")

    # Which requested months actually have data? (zero-coverage months are dropped)
    requested_months = [
        r["dateFrom"][:7]
        for r in build_monthly_ranges(
            datetime.strptime(date_from_str, "%Y-%m-%d"),
            datetime.strptime(date_to_str, "%Y-%m-%d"),
        )
    ]
    months_with_data = sorted({
        m
        for per_month in time_series.values()
        for m, v in per_month.items()
        if v.get("coverage", 0) >= min_coverage
    })
    missing_months = [m for m in requested_months if m not in months_with_data]

    if months_with_data:
        st.caption(
            f"Site: {selected_location} | requested {date_from_str} to {date_to_str} | "
            f"data available: {months_with_data[0]} to {months_with_data[-1]}"
        )
    if missing_months:
        shown = ", ".join(missing_months[:6])
        if len(missing_months) > 6:
            shown += f" … (+{len(missing_months) - 6} more)"
        st.info(
            f"No usable sequencing data (≥{min_coverage:,} reads) for {len(missing_months)} of {len(requested_months)} "
            f"requested month(s): {shown}. Results for these months are missing or unreliable."
        )
    html = build_html_heatmap(heatmap_data)
    n_heatmap_rows = sum(
        1 for info in heatmap_data.values()
        if not all(cell["status"] == "clean" for cell in info["months"].values())
    )
    # legend + header + rows; the colour bar needs ~270px; lower rows' tooltips open
    # upward, so only small tables need extra room below for the top row's tooltip
    heatmap_height = max(110 + n_heatmap_rows * 46, 270) + (30 if n_heatmap_rows >= 4 else 120)
    components.html(html, height=heatmap_height, scrolling=True)
    hidden = [
        name for name, info in heatmap_data.items()
        if all(cell["status"] == "clean" for cell in info["months"].values())
    ]
    if hidden:
        st.caption(
            f"Not shown (no mutation ≥ {detect_threshold:.1%} at any of their positions): "
            + ", ".join(hidden)
        )

    # ── Summary table ──────────────────────────────────────────────────────────
    st.subheader("Summary")

    summary_html, assay_verdicts = build_summary_html(
        position_data, found, time_series, dominant_letters,
        three_prime_n=three_prime_n,
        min_coverage=min_coverage,
        mismatch_threshold=mismatch_threshold,
        watch_threshold=watch_threshold,
    )

    # verdict per assay (ddPCR) / amplicon (ARTIC)
    icons = {"OK": "🟢", "Watch": "🟠", "Compromised": "🔴"}
    if len(assay_verdicts) <= 6:
        cols = st.columns(len(assay_verdicts))
        for col, (assay, v) in zip(cols, assay_verdicts.items()):
            with col:
                with st.container(border=True):
                    st.markdown(f"**{assay}** — {icons[v['status']]} {v['status']}")
                    st.caption(v["reason"])
    else:
        counts = {k: sum(1 for v in assay_verdicts.values() if v["status"] == k) for k in icons}
        st.markdown(
            f"**{len(assay_verdicts)} assays/amplicons:** "
            f"🔴 {counts['Compromised']} compromised · 🟠 {counts['Watch']} watch · 🟢 {counts['OK']} OK"
        )
    st.caption(
        f"Based on the most recent month. **Compromised** = probe differs from the circulating "
        f"virus, or a primer mismatch in the last {three_prime_n} bases (3′ end), or 2+ primer "
        f"mismatches. **Watch** = 1 primer mismatch elsewhere, or a minority mutation "
        f"≥{watch_threshold:.0%}. A mutation counts once it is above {mismatch_threshold:.0%}."
    )

    with st.expander("Details per primer/probe", expanded=True):
        # height follows the table (+ assay header rows); popovers on the last rows open
        # upward, so only short tables need extra room below
        n_rows = len([r for r in found if r["Start"] is not None])
        n_groups = len(assay_verdicts) if any(len(v["pieces"]) > 1 for v in assay_verdicts.values()) else 0
        table_height = 50 + n_rows * 62 + n_groups * 42 + (30 if n_rows >= 5 else 170)
        components.html(summary_html, height=table_height, scrolling=True)

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