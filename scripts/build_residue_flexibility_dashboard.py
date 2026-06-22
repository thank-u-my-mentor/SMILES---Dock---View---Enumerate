#!/usr/bin/env python3
"""Build a Chinese residue-flexibility dashboard from GROMACS C-alpha RMSF."""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
from pathlib import Path
from statistics import mean, median


AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "ASH": "D",
    "CYS": "C",
    "CYX": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLH": "E",
    "GLY": "G",
    "HIS": "H",
    "HID": "H",
    "HIE": "H",
    "HIP": "H",
    "HD1": "H",
    "HD2": "H",
    "HD3": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "LYN": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def parse_xvg(path: Path) -> list[tuple[int, float]]:
    rows: list[tuple[int, float]] = []
    for line in path.read_text(errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "@")):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        rows.append((int(float(parts[0])), float(parts[1])))
    return rows


def parse_ca_pdb(path: Path) -> list[dict[str, str | int]]:
    residues: list[dict[str, str | int]] = []
    seen: set[tuple[str, int, str]] = set()
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom = line[12:16].strip()
        if atom != "CA":
            continue
        resname = line[17:20].strip()
        if resname not in AA3_TO_1:
            continue
        chain = line[21].strip() or "A"
        try:
            resid = int(line[22:26])
        except ValueError:
            continue
        icode = line[26].strip()
        key = (chain, resid, icode)
        if key in seen:
            continue
        seen.add(key)
        residues.append(
            {
                "chain": chain,
                "resid": resid,
                "icode": icode,
                "resname": resname,
                "aa": AA3_TO_1.get(resname, "X"),
            }
        )
    return residues


def percentile_rank(sorted_values: list[float], value: float) -> float:
    if not sorted_values:
        return 0.0
    below = sum(1 for x in sorted_values if x <= value)
    return 100.0 * below / len(sorted_values)


def classify(percentile: float) -> str:
    if percentile >= 90:
        return "高自由度"
    if percentile >= 75:
        return "偏高"
    if percentile <= 10:
        return "刚性核心"
    if percentile <= 25:
        return "偏低"
    return "中等"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = [
        "rank_flexible",
        "label",
        "chain",
        "resid",
        "icode",
        "resname",
        "aa",
        "md_internal_index",
        "rmsf_nm",
        "rmsf_A",
        "percentile",
        "class",
        "metal_site",
        "mapping_note",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def bar(width: float, cls: str) -> str:
    color = {
        "高自由度": "#c43b3b",
        "偏高": "#e68645",
        "中等": "#5f9ea0",
        "偏低": "#5d7fb9",
        "刚性核心": "#50607a",
    }.get(cls, "#5f9ea0")
    return (
        f'<div class="barbox"><span class="bar" '
        f'style="width:{width:.1f}%; background:{color}"></span></div>'
    )


def write_html(path: Path, rows: list[dict[str, object]], title: str) -> None:
    values = [float(r["rmsf_A"]) for r in rows]
    high_count = sum(1 for r in rows if r["class"] == "高自由度")
    rigid_count = sum(1 for r in rows if r["class"] == "刚性核心")
    metal_count = sum(1 for r in rows if r.get("metal_site"))
    rows_by_resid = sorted(rows, key=lambda r: (str(r["chain"]), int(r["resid"])))
    json_rows = []
    for r in rows:
        json_rows.append(
            {
                "rank": int(r["rank_flexible"]),
                "label": str(r["label"]),
                "chain": str(r["chain"]),
                "resid": int(r["resid"]),
                "icode": str(r["icode"]),
                "resname": str(r["resname"]),
                "aa": str(r["aa"]),
                "md_index": int(r["md_internal_index"]),
                "rmsf": float(r["rmsf_A"]),
                "percentile": float(r["percentile"]),
                "class": str(r["class"]),
                "mapping": str(r["mapping_note"]),
                "metal": bool(r.get("metal_site")),
            }
        )
    json_sequence = [
        {
            "label": str(r["label"]),
            "resid": int(r["resid"]),
            "rmsf": float(r["rmsf_A"]),
            "class": str(r["class"]),
            "metal": bool(r.get("metal_site")),
        }
        for r in rows_by_resid
    ]

    html_text = f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<style>
  :root {{
    --bg: #f6f7f9;
    --panel: #ffffff;
    --ink: #20252d;
    --muted: #667085;
    --line: #d9dee7;
    --blue: #2f5f8f;
    --green: #287a5b;
    --red: #a94040;
    --yellow: #8b6f17;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0;
    background: var(--bg);
    color: var(--ink);
    font-family: Arial, Helvetica, "Microsoft YaHei", sans-serif;
    letter-spacing: 0;
  }}
  header {{
    padding: 18px 22px 12px;
    border-bottom: 1px solid var(--line);
    background: var(--panel);
    position: sticky;
    top: 0;
    z-index: 5;
  }}
  h1 {{ margin: 0 0 4px; font-size: 22px; line-height: 1.2; }}
  .subtitle {{ margin: 0; color: var(--muted); font-size: 13px; }}
  .layout {{
    display: grid;
    grid-template-columns: 390px minmax(0, 1fr);
    gap: 14px;
    padding: 14px;
  }}
  aside, main, .formula {{
    background: var(--panel);
    border: 1px solid var(--line);
    border-radius: 8px;
  }}
  aside {{
    height: calc(100vh - 98px);
    overflow: hidden;
    display: grid;
    grid-template-rows: auto auto 1fr;
  }}
  .controls {{
    padding: 12px;
    border-bottom: 1px solid var(--line);
    display: grid;
    gap: 8px;
  }}
  input, select, button {{
    font: inherit;
    border: 1px solid var(--line);
    border-radius: 6px;
    padding: 8px 10px;
    background: #fff;
    color: var(--ink);
  }}
  .summary {{
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 8px;
    padding: 10px 12px;
    border-bottom: 1px solid var(--line);
    font-size: 12px;
  }}
  .summary strong {{ display: block; font-size: 18px; line-height: 1.1; }}
  .table-wrap {{ overflow: auto; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 12px; }}
  th, td {{
    padding: 8px 9px;
    border-bottom: 1px solid #edf0f4;
    text-align: right;
    white-space: nowrap;
  }}
  th:first-child, td:first-child {{ text-align: left; }}
  th {{ position: sticky; top: 0; background: #f9fafb; z-index: 2; }}
  tr {{ cursor: pointer; }}
  tr:hover {{ background: #f3f6fa; }}
  tr.active {{ background: #e9f1fb; }}
  .tag {{
    display: inline-block;
    border-radius: 999px;
    padding: 1px 6px;
    font-size: 11px;
    border: 1px solid var(--line);
    color: var(--muted);
    margin-left: 5px;
  }}
  .tag.metal {{ color: #8b6f17; border-color: #dccb80; background: #fff8d7; }}
  main {{
    min-height: calc(100vh - 98px);
    padding: 14px;
    display: grid;
    grid-template-rows: auto auto 1fr auto;
    gap: 12px;
  }}
  .formula {{
    padding: 14px 16px;
    line-height: 1.5;
    font-size: 13px;
    background: linear-gradient(180deg, #ffffff 0%, #fbfcfe 100%);
  }}
  .formula h2 {{ margin: 0 0 6px; font-size: 15px; }}
  .math {{
    display: grid;
    grid-template-columns: repeat(3, minmax(0, 1fr));
    gap: 10px;
    margin: 8px 0 6px;
  }}
  .math-card {{
    border: 1px solid #dfe5ee;
    border-radius: 7px;
    padding: 10px 12px;
    background: #ffffff;
    box-shadow: 0 1px 2px rgba(16, 24, 40, 0.04);
    min-width: 0;
  }}
  .math-card strong {{
    display: block;
    font-size: 12px;
    margin-bottom: 7px;
    color: #344054;
  }}
  math {{
    font-family: "Cambria Math", "STIX Two Math", "Times New Roman", serif;
    font-size: 19px;
    overflow-x: auto;
    max-width: 100%;
    color: #1f2937;
  }}
  .note {{ color: var(--muted); font-size: 12px; margin-top: 6px; }}
  .sort-note {{
    margin-top: 8px;
    padding: 8px 10px;
    border-left: 3px solid var(--blue);
    background: #f4f7fb;
    color: #344054;
    font-size: 12px;
  }}
  .detail-grid {{
    display: grid;
    grid-template-columns: repeat(5, minmax(120px, 1fr));
    gap: 8px;
  }}
  .metric {{
    background: #fbfcfd;
    border: 1px solid var(--line);
    border-radius: 7px;
    padding: 9px 10px;
  }}
  .metric span {{
    display: block;
    color: var(--muted);
    font-size: 11px;
    margin-bottom: 4px;
  }}
  .metric strong {{ font-size: 17px; line-height: 1.1; }}
  .charts {{
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 12px;
    min-height: 420px;
  }}
  .chart-card {{
    border: 1px solid var(--line);
    border-radius: 8px;
    background: #fff;
    padding: 10px;
    min-width: 0;
  }}
  .chart-card h2 {{ margin: 0 0 8px; font-size: 15px; }}
  canvas {{ width: 100%; height: 390px; display: block; }}
  .chart-note {{
    margin: 8px 0 0;
    min-height: 34px;
    color: #667085;
    font-size: 12px;
    line-height: 1.42;
  }}
  .files {{ font-size: 12px; color: var(--muted); }}
  @media (max-width: 980px) {{
    .layout {{ grid-template-columns: 1fr; }}
    aside {{ height: 520px; }}
    .charts {{ grid-template-columns: 1fr; }}
    .detail-grid {{ grid-template-columns: repeat(2, 1fr); }}
    .math {{ grid-template-columns: 1fr; }}
  }}
</style>
</head>
<body>
<header>
  <h1>{html.escape(title)}</h1>
  <p class="subtitle">MCPB-Fe(II) bonded model；50 ns production，约 10 ps/帧；PBC molecule centering + Protein_All_CA rot+trans fit 后计算 Cα RMSF。</p>
</header>
<div class="layout">
  <aside>
    <section class="controls">
      <input id="search" placeholder="搜索氨基酸，例如 HID51 或 GLU:A:319">
      <select id="filter">
        <option value="all">全部氨基酸</option>
        <option value="high">高自由度：前 10%</option>
        <option value="rigid">刚性核心：后 10%</option>
        <option value="metal">Fe 配位中心</option>
      </select>
      <select id="sort">
        <option value="rmsf_desc">默认排序：RMSF 从大到小</option>
        <option value="rmsf_asc">RMSF 从小到大</option>
        <option value="percentile_desc">自由度百分位从高到低</option>
        <option value="resid">残基编号从小到大</option>
      </select>
    </section>
    <section class="summary">
      <div><strong id="residueCount">{len(rows)}</strong>个氨基酸</div>
      <div><strong>{high_count}</strong>个高自由度</div>
      <div><strong>{metal_count}</strong>个 Fe 配位</div>
    </section>
    <div class="table-wrap">
      <table>
        <thead><tr><th>氨基酸</th><th>RMSF Å</th><th>百分位</th><th>类别</th></tr></thead>
        <tbody id="residueRows"></tbody>
      </table>
    </div>
  </aside>
  <main>
    <section class="formula">
      <h2>氨基酸自由度 RMSF 的定义</h2>
      <div class="math">
        <div class="math-card">
          <strong>残基 i 的 Cα 位置</strong>
          <math display="block">
            <msub><mi>r</mi><mi>i</mi></msub><mo>(</mo><mi>t</mi><mo>)</mo>
            <mo>=</mo><mtext>Cα coordinate after fit</mtext>
          </math>
        </div>
        <div class="math-card">
          <strong>平均位置</strong>
          <math display="block">
            <mo>&lt;</mo><msub><mi>r</mi><mi>i</mi></msub><mo>&gt;</mo>
            <mo>=</mo>
            <mfrac><mn>1</mn><mi>N</mi></mfrac>
            <munderover><mo>∑</mo><mrow><mi>k</mi><mo>=</mo><mn>1</mn></mrow><mi>N</mi></munderover>
            <msub><mi>r</mi><mi>i</mi></msub><mo>(</mo><msub><mi>t</mi><mi>k</mi></msub><mo>)</mo>
          </math>
        </div>
        <div class="math-card">
          <strong>Cα RMSF / 自由度</strong>
          <math display="block">
            <msub><mtext>RMSF</mtext><mi>i</mi></msub>
            <mo>=</mo>
            <msqrt>
              <mfrac><mn>1</mn><mi>N</mi></mfrac>
              <munderover><mo>∑</mo><mrow><mi>k</mi><mo>=</mo><mn>1</mn></mrow><mi>N</mi></munderover>
              <msup>
                <mrow><mo>|</mo><mo>|</mo><msub><mi>r</mi><mi>i</mi></msub><mo>(</mo><msub><mi>t</mi><mi>k</mi></msub><mo>)</mo><mo>-</mo><mo>&lt;</mo><msub><mi>r</mi><mi>i</mi></msub><mo>&gt;</mo><mo>|</mo><mo>|</mo></mrow>
                <mn>2</mn>
              </msup>
            </msqrt>
          </math>
        </div>
      </div>
      <div class="note">RMSF 越大，表示该残基 Cα 在对齐后的轨迹中摆动越大；RMSF 越小，表示该区域更刚性。这里的自由度是本次 50 ns 经典 MD 内部的相对指标，不等于实验 B-factor，也不直接等价于催化贡献。</div>
      <div class="sort-note">排序逻辑：左侧默认按 <strong>RMSF 从大到小</strong> 排序；类别按体系内部分位数划分，前 10% 记为“高自由度”，后 10% 记为“刚性核心”。Fe 配位中心由 MCPB metal_coordination_candidates.csv 标注。</div>
    </section>
    <section class="detail-grid" id="detailGrid"></section>
    <section class="charts">
      <div class="chart-card">
        <h2 id="sequenceTitle">沿序列的自由度分布</h2>
        <canvas id="sequenceCanvas"></canvas>
        <p class="chart-note" id="sequenceNote">纵轴为 Cα RMSF，单位 Å；横轴为参考 PDB/H++ 文件中的残基顺序。黄色点表示 Fe 配位中心。</p>
      </div>
      <div class="chart-card">
        <h2 id="histTitle">RMSF 分布密度</h2>
        <canvas id="histCanvas"></canvas>
        <p class="chart-note" id="histNote">直方图展示所有氨基酸自由度的分布；虚线标出当前选中的残基。</p>
      </div>
    </section>
    <section class="files">
      数据文件：<code>residue_ca_rmsf.csv</code>, <code>protein_ca_rmsf_residue.xvg</code>, <code>protein_ca_average.pdb</code>。PBC/fit 轨迹：<code>production_50ns_protein_ca_fit.xtc</code>。
    </section>
  </main>
</div>
<script id="flex-data" type="application/json">{json.dumps({"rows": json_rows, "sequence": json_sequence, "summary": {"mean": mean(values), "median": median(values), "max": max(values), "min": min(values), "high": high_count, "rigid": rigid_count}}, ensure_ascii=False)}</script>
<script>
const payload = JSON.parse(document.getElementById('flex-data').textContent);
const rows = payload.rows;
const sequence = payload.sequence;
const summary = payload.summary;
let selected = rows[0];

const searchEl = document.getElementById('search');
const filterEl = document.getElementById('filter');
const sortEl = document.getElementById('sort');
const tbody = document.getElementById('residueRows');
const detailGrid = document.getElementById('detailGrid');
const seqCanvas = document.getElementById('sequenceCanvas');
const histCanvas = document.getElementById('histCanvas');

function clsColor(cls) {{
  if (cls === '高自由度') return '#a94040';
  if (cls === '偏高') return '#d88937';
  if (cls === '刚性核心') return '#50607a';
  if (cls === '偏低') return '#2f5f8f';
  return '#287a5b';
}}

function filteredRows() {{
  const q = searchEl.value.trim().toLowerCase();
  let out = rows.filter(r => {{
    const text = (r.label + ' ' + r.resname + r.aa).toLowerCase();
    if (q && !text.includes(q)) return false;
    if (filterEl.value === 'high' && r.class !== '高自由度') return false;
    if (filterEl.value === 'rigid' && r.class !== '刚性核心') return false;
    if (filterEl.value === 'metal' && !r.metal) return false;
    return true;
  }});
  const s = sortEl.value;
  out.sort((a, b) => {{
    if (s === 'rmsf_asc') return a.rmsf - b.rmsf;
    if (s === 'percentile_desc') return b.percentile - a.percentile;
    if (s === 'resid') return a.resid - b.resid;
    return b.rmsf - a.rmsf;
  }});
  return out;
}}

function renderTable() {{
  const out = filteredRows();
  document.getElementById('residueCount').textContent = out.length;
  tbody.innerHTML = '';
  out.forEach(r => {{
    const tr = document.createElement('tr');
    if (selected && r.label === selected.label) tr.className = 'active';
    const tag = r.metal ? '<span class="tag metal">Fe</span>' : (r.class === '高自由度' ? '<span class="tag">top 10%</span>' : '');
    tr.innerHTML = `<td><b>${{r.label}}</b>${{tag}}</td><td>${{r.rmsf.toFixed(3)}}</td><td>${{r.percentile.toFixed(1)}}</td><td>${{r.class}}</td>`;
    tr.onclick = () => {{ selected = r; renderAll(); }};
    tbody.appendChild(tr);
  }});
}}

function metric(label, value) {{
  return `<div class="metric"><span>${{label}}</span><strong>${{value}}</strong></div>`;
}}

function renderDetails() {{
  detailGrid.innerHTML = [
    metric('当前残基', selected.label + (selected.metal ? ' · Fe' : '')),
    metric('RMSF', selected.rmsf.toFixed(3) + ' Å'),
    metric('自由度百分位', selected.percentile.toFixed(1) + '%'),
    metric('类别', selected.class),
    metric('平均 / 中位数', summary.mean.toFixed(3) + ' / ' + summary.median.toFixed(3) + ' Å')
  ].join('');
}}

function setupCanvas(canvas) {{
  const dpr = window.devicePixelRatio || 1;
  const rect = canvas.getBoundingClientRect();
  canvas.width = Math.max(1, Math.floor(rect.width * dpr));
  canvas.height = Math.max(1, Math.floor(rect.height * dpr));
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return {{ctx, w: rect.width, h: rect.height}};
}}

function drawAxes(ctx, w, h, left, top, right, bottom, ylabel) {{
  ctx.strokeStyle = '#ccd3de';
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(left, top);
  ctx.lineTo(left, h - bottom);
  ctx.lineTo(w - right, h - bottom);
  ctx.stroke();
  ctx.fillStyle = '#667085';
  ctx.font = '12px Arial';
  ctx.save();
  ctx.translate(14, top + 90);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText(ylabel, 0, 0);
  ctx.restore();
}}

function renderSequence() {{
  const {{ctx, w, h}} = setupCanvas(seqCanvas);
  ctx.clearRect(0, 0, w, h);
  const left = 52, right = 18, top = 20, bottom = 42;
  const maxY = Math.max(...sequence.map(d => d.rmsf)) * 1.08;
  drawAxes(ctx, w, h, left, top, right, bottom, 'RMSF (Å)');
  ctx.strokeStyle = '#2f5f8f';
  ctx.lineWidth = 1.8;
  ctx.beginPath();
  sequence.forEach((d, i) => {{
    const x = left + i * (w - left - right) / Math.max(1, sequence.length - 1);
    const y = h - bottom - d.rmsf / maxY * (h - top - bottom);
    if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  }});
  ctx.stroke();
  sequence.forEach((d, i) => {{
    if (!d.metal && d.label !== selected.label) return;
    const x = left + i * (w - left - right) / Math.max(1, sequence.length - 1);
    const y = h - bottom - d.rmsf / maxY * (h - top - bottom);
    ctx.beginPath();
    ctx.fillStyle = d.label === selected.label ? '#a94040' : '#8b6f17';
    ctx.arc(x, y, d.label === selected.label ? 5 : 4, 0, Math.PI * 2);
    ctx.fill();
  }});
  ctx.fillStyle = '#667085';
  ctx.font = '12px Arial';
  ctx.fillText('N-term', left, h - 14);
  ctx.fillText('C-term', w - right - 42, h - 14);
  ctx.fillText(maxY.toFixed(1), 14, top + 4);
  ctx.fillText('0', 32, h - bottom + 4);
}}

function renderHist() {{
  const {{ctx, w, h}} = setupCanvas(histCanvas);
  ctx.clearRect(0, 0, w, h);
  const left = 52, right = 18, top = 20, bottom = 42;
  const vals = rows.map(r => r.rmsf);
  const minV = Math.min(...vals), maxV = Math.max(...vals);
  const bins = 28;
  const counts = new Array(bins).fill(0);
  vals.forEach(v => {{
    const idx = Math.min(bins - 1, Math.floor((v - minV) / Math.max(1e-9, maxV - minV) * bins));
    counts[idx]++;
  }});
  const maxC = Math.max(...counts);
  drawAxes(ctx, w, h, left, top, right, bottom, 'Residues');
  counts.forEach((c, i) => {{
    const x0 = left + i * (w - left - right) / bins;
    const bw = (w - left - right) / bins - 2;
    const bh = c / maxC * (h - top - bottom);
    const hue = 205 - 145 * i / Math.max(1, bins - 1);
    ctx.fillStyle = `hsl(${{hue}}, 48%, 58%)`;
    ctx.fillRect(x0, h - bottom - bh, bw, bh);
  }});
  const sx = left + (selected.rmsf - minV) / Math.max(1e-9, maxV - minV) * (w - left - right);
  ctx.strokeStyle = '#a94040';
  ctx.setLineDash([5, 4]);
  ctx.beginPath();
  ctx.moveTo(sx, top);
  ctx.lineTo(sx, h - bottom);
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = '#667085';
  ctx.font = '12px Arial';
  ctx.fillText(minV.toFixed(2) + ' Å', left, h - 14);
  ctx.fillText(maxV.toFixed(2) + ' Å', w - right - 48, h - 14);
}}

function renderAll() {{
  renderTable();
  renderDetails();
  renderSequence();
  renderHist();
}}

[searchEl, filterEl, sortEl].forEach(el => el.addEventListener('input', renderAll));
window.addEventListener('resize', renderAll);
renderAll();
</script>
</body>
</html>
"""
    path.write_text(html_text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xvg", required=True, type=Path)
    parser.add_argument("--ca-pdb", required=True, type=Path)
    parser.add_argument("--reference-pdb", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--title", default="氨基酸自由度分析")
    parser.add_argument("--metal-csv", type=Path)
    args = parser.parse_args()

    xvg_rows = parse_xvg(args.xvg)
    md_residues = parse_ca_pdb(args.ca_pdb)
    ref_residues = parse_ca_pdb(args.reference_pdb) if args.reference_pdb else []

    if len(xvg_rows) != len(md_residues):
        raise SystemExit(
            f"RMSF row count ({len(xvg_rows)}) != CA PDB count ({len(md_residues)})"
        )
    if ref_residues and len(ref_residues) != len(md_residues):
        raise SystemExit(
            f"Reference CA count ({len(ref_residues)}) != MD CA count ({len(md_residues)})"
        )

    rmsf_values = [v for _, v in xvg_rows]
    sorted_values = sorted(rmsf_values)
    rows: list[dict[str, object]] = []
    metal_ids: set[str] = set()
    if args.metal_csv and args.metal_csv.exists():
        with args.metal_csv.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if row.get("likely_coordination_atom", "").lower() == "true":
                    metal_ids.add(row.get("residue_id", ""))
    for idx, ((internal_index, rmsf_nm), md_res) in enumerate(zip(xvg_rows, md_residues)):
        ref_res = ref_residues[idx] if ref_residues else md_res
        mapping_note = "ok"
        if ref_res["aa"] != md_res["aa"]:
            mapping_note = f"check: MD {md_res['resname']}{md_res['resid']}"
        rmsf_a = rmsf_nm * 10.0
        pct = percentile_rank(sorted_values, rmsf_nm)
        label = f"{ref_res['resname']}:{ref_res['chain']}:{ref_res['resid']}{ref_res['icode']}"
        residue_id = f"{ref_res['resname']}:{ref_res['chain']}:{ref_res['resid']}{ref_res['icode']}"
        rows.append(
            {
                "label": label,
                "chain": ref_res["chain"],
                "resid": ref_res["resid"],
                "icode": ref_res["icode"],
                "resname": ref_res["resname"],
                "aa": ref_res["aa"],
                "md_internal_index": internal_index,
                "rmsf_nm": f"{rmsf_nm:.6f}",
                "rmsf_A": f"{rmsf_a:.4f}",
                "percentile": f"{pct:.2f}",
                "class": classify(pct),
                "mapping_note": mapping_note,
                "metal_site": residue_id in metal_ids,
            }
        )

    rows.sort(key=lambda r: float(r["rmsf_A"]), reverse=True)
    for rank, row in enumerate(rows, 1):
        row["rank_flexible"] = rank

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "residue_ca_rmsf.csv", rows)
    write_html(args.outdir / "index.html", rows, args.title)
    print(f"wrote={args.outdir / 'index.html'}")
    print(f"csv={args.outdir / 'residue_ca_rmsf.csv'}")
    print(f"residues={len(rows)} mean_A={mean([float(r['rmsf_A']) for r in rows]):.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
