"""
components/jaccard_heatmap.py

Jaccard similarity heatmap for selected variants.

Shows pairwise mutation signature overlap between panel variants.
Rendered as HTML/JS via streamlit.components.v1 for precise visual control.

Green = distinct signatures (low overlap).
Red = near-identical signatures (high overlap, harder to discriminate).
"""
from __future__ import annotations

import json
import streamlit as st
import streamlit.components.v1 as components


def compute_jaccard_matrix(
    variants: list[str],
    pango_loader,
) -> list[list[float]]:
    """
    Compute pairwise Jaccard similarity matrix.
    Jaccard(A, B) = |sig(A) ∩ sig(B)| / |sig(A) ∪ sig(B)|
    Diagonal is always 1.0.
    """
    n = len(variants)
    matrix = [[0.0] * n for _ in range(n)]
    for i, v1 in enumerate(variants):
        for j, v2 in enumerate(variants):
            if i == j:
                matrix[i][j] = 1.0
                continue
            s1 = {m for m in pango_loader.get_signature(v1) if not m.endswith('-')}
            s2 = {m for m in pango_loader.get_signature(v2) if not m.endswith('-')}
            if not s1 or not s2:
                matrix[i][j] = 0.0
                continue
            intersection = len(s1 & s2)
            union = len(s1 | s2)
            matrix[i][j] = intersection / union if union > 0 else 0.0
    return matrix


def render_jaccard_heatmap(
    variants: list[str],
    pango_loader,
) -> None:
    """
    Render a Jaccard similarity heatmap for the selected variants.

    Args:
        variants: List of selected pango lineage names
        pango_loader: PangoLoader instance
    """
    if len(variants) < 2:
        st.caption("Select at least 2 variants to see similarity matrix.")
        return

    matrix = compute_jaccard_matrix(variants, pango_loader)

    variants_json = json.dumps(variants)
    matrix_json = json.dumps(matrix)

    cell_px = 52
    label_w = max(len(v) for v in variants) * 7 + 16
    table_w = label_w + len(variants) * (cell_px + 2) + 20
    table_h = len(variants) * (cell_px - 8 + 2) + 80

    html = f"""
<!DOCTYPE html>
<html>
<head>
<style>
*{{box-sizing:border-box;margin:0;padding:0;}}
body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;
     padding:8px 4px;background:transparent;}}
.title{{font-size:13px;font-weight:500;color:#111;margin-bottom:3px;}}
.sub{{font-size:11px;color:#888;margin-bottom:12px;}}
table{{border-collapse:collapse;}}
th{{font-size:11px;font-weight:400;color:#888;padding:3px 6px;white-space:nowrap;}}
th.rh{{text-align:right;padding-right:10px;min-width:{label_w}px;}}
td{{width:{cell_px}px;height:44px;text-align:center;font-size:11px;
    font-weight:500;border:2px solid #fff;}}
.scale-wrap{{display:flex;align-items:center;gap:8px;margin-top:12px;}}
.scale-bar{{height:7px;width:140px;border-radius:3px;}}
.sl{{display:flex;justify-content:space-between;width:140px;
     font-size:10px;color:#888;margin-top:2px;}}
.cap{{font-size:10px;color:#999;margin-top:8px;line-height:1.5;}}
</style>
</head>
<body>
<div class="title">Signature similarity (Jaccard index)</div>
<div class="sub">Pairwise mutation overlap between selected variants</div>
<div id="heatmap"></div>
<div class="scale-wrap">
  <span style="font-size:10px;color:#888;">distinct</span>
  <div>
    <div class="scale-bar" style="background:linear-gradient(to right,#3B6D11,#97C459,#F1EFE8,#EF9F27,#712B13);"></div>
    <div class="sl"><span>0.0</span><span>0.5</span><span>1.0</span></div>
  </div>
  <span style="font-size:10px;color:#888;">identical</span>
</div>
<div class="cap">Green = distinct signatures. Red = near-identical, harder to discriminate in deconvolution. Diagonal = variant vs itself (1.0).</div>

<script>
const variants = {variants_json};
const values = {matrix_json};
const stops = [
  [0.0,  '#3B6D11'],
  [0.25, '#97C459'],
  [0.5,  '#F1EFE8'],
  [0.75, '#EF9F27'],
  [1.0,  '#712B13'],
];

function hexToRgb(h) {{
  return [parseInt(h.slice(1,3),16), parseInt(h.slice(3,5),16), parseInt(h.slice(5,7),16)];
}}
function rgbToHex(r,g,b) {{
  return '#'+[r,g,b].map(v=>Math.round(v).toString(16).padStart(2,'0')).join('');
}}
function lerpColor(c0,c1,f) {{
  const [r0,g0,b0]=hexToRgb(c0), [r1,g1,b1]=hexToRgb(c1);
  return rgbToHex(r0+(r1-r0)*f, g0+(g1-g0)*f, b0+(b1-b0)*f);
}}
function interpolate(t) {{
  for (let i=0;i<stops.length-1;i++) {{
    const [t0,c0]=stops[i], [t1,c1]=stops[i+1];
    if (t<=t1) return lerpColor(c0,c1,(t-t0)/(t1-t0));
  }}
  return stops[stops.length-1][1];
}}
function textColor(hex) {{
  const [r,g,b]=hexToRgb(hex);
  return (0.299*r+0.587*g+0.114*b)>160 ? '#444441' : '#ffffff';
}}

const table = document.createElement('table');
const thead = document.createElement('tr');
const corner = document.createElement('th');
corner.className='rh';
thead.appendChild(corner);
variants.forEach(v=>{{
  const th=document.createElement('th');
  th.textContent=v;
  thead.appendChild(th);
}});
table.appendChild(thead);

values.forEach((row,i)=>{{
  const tr=document.createElement('tr');
  const th=document.createElement('th');
  th.className='rh';
  th.textContent=variants[i];
  tr.appendChild(th);
  row.forEach((val,j)=>{{
    const td=document.createElement('td');
    const bg=interpolate(val);
    td.style.background=bg;
    td.style.color=textColor(bg);
    td.textContent=val.toFixed(2);
    td.title=`${{variants[i]}} ↔ ${{variants[j]}}: ${{val.toFixed(3)}}`;
    tr.appendChild(td);
  }});
  table.appendChild(tr);
}});
document.getElementById('heatmap').appendChild(table);
</script>
</body>
</html>
"""

    height = table_h + 80
    components.html(html, height=height, scrolling=False)