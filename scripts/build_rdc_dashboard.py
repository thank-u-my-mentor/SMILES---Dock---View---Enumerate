#!/usr/bin/env python
"""Build an interactive RDC dashboard from residue-ligand MD analysis CSV files."""

from __future__ import annotations

import argparse
import csv
import html
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rdc-dir", required=True, type=Path)
    parser.add_argument("--out", type=Path, help="Output HTML; defaults to <rdc-dir>/index.html")
    parser.add_argument("--title", default="氨基酸-底物 RDC 动态分析 Dashboard")
    parser.add_argument("--subtitle", default="探索性 20 ns MD；每 10 ps 采样一次氨基酸-底物距离。")
    parser.add_argument(
        "--timeseries",
        type=Path,
        help="Distance time-series CSV. Defaults to residue_ligand_distance_timeseries_all.csv, then residue_ligand_distance_timeseries.csv.",
    )
    parser.add_argument("--nearby-csv", type=Path, help="CSV defining residues inside the original contact cutoff")
    parser.add_argument("--velocity-threshold", type=float, default=3.0, help="Reference threshold in V x 10^-2 A/ps units")
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def to_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def read_timeseries(path: Path) -> tuple[list[float], dict[str, list[float]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        residues = header[1:]
        times: list[float] = []
        traces: dict[str, list[float]] = {residue: [] for residue in residues}
        for row in reader:
            if not row:
                continue
            times.append(float(row[0]))
            for residue, value in zip(residues, row[1:]):
                traces[residue].append(float(value) if value else None)  # type: ignore[arg-type]
    return times, traces


def compact_row(row: dict[str, str], nearby: set[str]) -> dict[str, object]:
    residue = row.get("residue", "")
    return {
        "residue": residue,
        "resname": row.get("resname", ""),
        "resid": int(float(row.get("resid") or 0)),
        "segid": row.get("segid", ""),
        "initial": to_float(row.get("initial_distance_A")),
        "final": to_float(row.get("final_distance_A")),
        "mean": to_float(row.get("mean_distance_A")),
        "min": to_float(row.get("min_distance_A")),
        "max": to_float(row.get("max_distance_A")),
        "delta": to_float(row.get("delta_final_minus_initial_A")),
        "slope": to_float(row.get("slope_A_per_ns")),
        "rdc": to_float(row.get("rdc_per_ns")),
        "r2": to_float(row.get("fit_r2")),
        "frames": int(float(row.get("n_frames") or 0)),
        "near10A": residue in nearby,
    }


def metric_summary(rows: list[dict[str, object]], trace_count: int) -> dict[str, object]:
    values = [row for row in rows if isinstance(row.get("rdc"), (float, int))]
    positive = sorted(values, key=lambda row: float(row.get("rdc") or 0.0), reverse=True)
    negative = sorted(values, key=lambda row: float(row.get("rdc") or 0.0))
    return {
        "residueCount": len(rows),
        "traceCount": trace_count,
        "near10ACount": sum(1 for row in rows if row.get("near10A")),
        "topPositive": positive[:5],
        "topNegative": negative[:5],
    }


def read_text_if_exists(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace") if path.exists() else ""


def main() -> None:
    args = parse_args()
    rdc_dir = args.rdc_dir
    out = args.out or (rdc_dir / "index.html")
    all_rows_path = rdc_dir / "residue_ligand_rdc_all.csv"
    nearby_path = args.nearby_csv or (rdc_dir / "residue_ligand_rdc_within_cutoff.csv")
    if args.timeseries:
        timeseries_path = args.timeseries
    else:
        all_ts = rdc_dir / "residue_ligand_distance_timeseries_all.csv"
        timeseries_path = all_ts if all_ts.exists() else rdc_dir / "residue_ligand_distance_timeseries.csv"

    raw_rows = read_csv_rows(all_rows_path)
    nearby = {row.get("residue", "") for row in read_csv_rows(nearby_path)}
    rows = [compact_row(row, nearby) for row in raw_rows]
    rows.sort(key=lambda row: abs(float(row.get("rdc") or 0.0)), reverse=True)
    times, traces = read_timeseries(timeseries_path)
    fe_core_summary = read_csv_rows(rdc_dir / "fe_core_distance_summary.csv")
    metadata = read_text_if_exists(rdc_dir / "rdc_metadata.json")
    summary = metric_summary(rows, len(traces))

    data = {
        "rows": rows,
        "times": times,
        "traces": traces,
        "summary": summary,
        "feCore": fe_core_summary,
        "metadata": metadata,
        "velocityThreshold": args.velocity_threshold,
        "timeseriesFile": timeseries_path.name,
    }
    data_json = json.dumps(data, separators=(",", ":"), ensure_ascii=False)
    title = html.escape(args.title)
    subtitle = html.escape(args.subtitle)

    html_text = f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
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
      font-family: Arial, Helvetica, sans-serif;
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
    h1 {{
      margin: 0 0 4px;
      font-size: 22px;
      line-height: 1.2;
    }}
    .subtitle {{
      margin: 0;
      color: var(--muted);
      font-size: 13px;
    }}
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
    .summary strong {{
      display: block;
      font-size: 18px;
      line-height: 1.1;
    }}
    .table-wrap {{ overflow: auto; }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 12px;
    }}
    th, td {{
      padding: 8px 9px;
      border-bottom: 1px solid #edf0f4;
      text-align: right;
      white-space: nowrap;
    }}
    th:first-child, td:first-child {{ text-align: left; }}
    th {{
      position: sticky;
      top: 0;
      background: #f9fafb;
      z-index: 2;
    }}
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
    main {{
      min-height: calc(100vh - 98px);
      padding: 14px;
      display: grid;
      grid-template-rows: auto auto 1fr;
      gap: 12px;
    }}
    .formula {{
      padding: 14px 16px;
      line-height: 1.5;
      font-size: 13px;
      background:
        linear-gradient(180deg, #ffffff 0%, #fbfcfe 100%);
    }}
    .formula h2 {{
      margin: 0 0 6px;
      font-size: 15px;
    }}
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
    .note {{
      color: var(--muted);
      font-size: 12px;
      margin-top: 6px;
    }}
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
    .metric strong {{
      font-size: 17px;
      line-height: 1.1;
    }}
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
    .chart-card h2 {{
      margin: 0 0 8px;
      font-size: 15px;
    }}
    canvas {{
      width: 100%;
      height: 390px;
      display: block;
    }}
    .chart-note {{
      margin: 8px 0 0;
      min-height: 34px;
      color: #667085;
      font-size: 12px;
      line-height: 1.42;
    }}
    .files {{
      font-size: 12px;
      color: var(--muted);
    }}
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
    <h1>{title}</h1>
    <p class="subtitle">{subtitle}</p>
  </header>
  <div class="layout">
    <aside>
      <section class="controls">
        <input id="search" placeholder="搜索氨基酸，例如 HID160 或 GLU:SYSTEM:319">
        <select id="filter">
          <option value="all">全部氨基酸</option>
          <option value="near">曾进入 10 Å 接触范围</option>
          <option value="pos">RDC 为正：整体远离底物</option>
          <option value="neg">RDC 为负：整体靠近底物</option>
        </select>
        <select id="sort">
          <option value="absrdc">默认排序：|RDC/ns| 从大到小</option>
          <option value="rdc_desc">RDC 从大到小</option>
          <option value="rdc_asc">RDC 从小到大</option>
          <option value="min_distance">最小距离从近到远</option>
          <option value="resid">残基编号从小到大</option>
        </select>
      </section>
      <section class="summary">
        <div><strong id="residueCount">0</strong>个氨基酸</div>
        <div><strong id="nearCount">0</strong>个进入 10 Å</div>
        <div><strong id="traceCount">0</strong>条距离轨迹</div>
      </section>
      <div class="table-wrap">
        <table>
          <thead>
            <tr><th>氨基酸</th><th>RDC/ns</th><th>最近 Å</th><th>Δ Å</th></tr>
          </thead>
          <tbody id="residueRows"></tbody>
        </table>
      </div>
    </aside>
    <main>
      <section class="formula">
        <h2>RDC 与距离波动速度的定义</h2>
        <div class="math">
          <div class="math-card">
            <strong>残基 i 到底物 L 的最近距离</strong>
            <math display="block">
              <msub><mi>d</mi><mi>i</mi></msub><mo>(</mo><mi>t</mi><mo>)</mo>
              <mo>=</mo>
              <munder><mi>min</mi><mrow><mi>a</mi><mo>&#x2208;</mo><mi>i</mi><mo>,</mo><mi>b</mi><mo>&#x2208;</mo><mi>L</mi></mrow></munder>
              <mrow><mo>|</mo><mo>|</mo><msub><mi>r</mi><mi>a</mi></msub><mo>(</mo><mi>t</mi><mo>)</mo><mo>-</mo><msub><mi>r</mi><mi>b</mi></msub><mo>(</mo><mi>t</mi><mo>)</mo><mo>|</mo><mo>|</mo></mrow>
            </math>
          </div>
          <div class="math-card">
            <strong>归一化距离变化率 RDC</strong>
            <math display="block">
              <msub><mtext>RDC</mtext><mi>i</mi></msub>
              <mo>=</mo>
              <mfrac><mn>1</mn><mrow><msub><mi>d</mi><mi>i</mi></msub><mo>(</mo><msub><mi>t</mi><mn>0</mn></msub><mo>)</mo></mrow></mfrac>
              <mo>&#x22C5;</mo>
              <mfrac><mrow><mi>d</mi><msub><mi>d</mi><mi>i</mi></msub><mo>(</mo><mi>t</mi><mo>)</mo></mrow><mrow><mi>d</mi><mi>t</mi></mrow></mfrac>
            </math>
          </div>
          <div class="math-card">
            <strong>相邻帧距离波动速度</strong>
            <math display="block">
              <msub><mi>V</mi><mi>i</mi></msub><mo>(</mo><msub><mi>t</mi><mi>k</mi></msub><mo>)</mo>
              <mo>=</mo>
              <mfrac>
                <mrow><msub><mi>d</mi><mi>i</mi></msub><mo>(</mo><msub><mi>t</mi><mrow><mi>k</mi><mo>+</mo><mn>1</mn></mrow></msub><mo>)</mo><mo>-</mo><msub><mi>d</mi><mi>i</mi></msub><mo>(</mo><msub><mi>t</mi><mi>k</mi></msub><mo>)</mo></mrow>
                <mrow><msub><mi>t</mi><mrow><mi>k</mi><mo>+</mo><mn>1</mn></mrow></msub><mo>-</mo><msub><mi>t</mi><mi>k</mi></msub></mrow>
              </mfrac>
            </math>
          </div>
        </div>
        <div class="note">RDC 在这里表示 20 ns 轨迹中“距离随时间的线性漂移”并按初始距离归一化。RDC 为正通常表示该氨基酸整体远离底物；RDC 为负通常表示整体靠近底物。</div>
        <div class="sort-note">排序逻辑：左侧列表默认按 <strong>|RDC/ns| 从大到小</strong> 排序，也就是优先展示距离变化趋势最强的氨基酸，而不是优先展示最近的氨基酸。进入 10 Å 接触范围的残基会带有 “10 Å” 标签。±{args.velocity_threshold:g} × 10<sup>-2</sup> Å/ps 只是本次 20 ns 短模拟的参考带，不作为催化贡献的硬阈值。</div>
      </section>
      <section class="detail-grid" id="detailGrid"></section>
      <section class="charts">
        <div class="chart-card">
          <h2 id="distanceTitle">距离轨迹</h2>
          <canvas id="distanceCanvas"></canvas>
          <p class="chart-note" id="distanceNote">纵轴距离表示该残基任意原子到底物任意原子的最近原子-原子距离，单位 Å。</p>
        </div>
        <div class="chart-card">
          <h2 id="velocityTitle">距离波动速度分布</h2>
          <canvas id="velocityCanvas"></canvas>
          <p class="chart-note" id="velocityNote">横轴 V 是相邻两帧距离变化除以时间间隔后的速度，并以 ×10⁻² Å/ps 显示。</p>
        </div>
      </section>
      <section class="files">
        数据文件：<code>residue_ligand_rdc_all.csv</code>, <code>{html.escape(timeseries_path.name)}</code>, <code>fe_core_distance_summary.csv</code>。Fe-core 质控数据保留在 CSV 文件中。
      </section>
    </main>
  </div>
  <script id="rdc-data" type="application/json">{data_json}</script>
  <script>
    const DATA = JSON.parse(document.getElementById('rdc-data').textContent);
    const rowsEl = document.getElementById('residueRows');
    const searchEl = document.getElementById('search');
    const filterEl = document.getElementById('filter');
    const sortEl = document.getElementById('sort');
    const detailGrid = document.getElementById('detailGrid');
    const distanceCanvas = document.getElementById('distanceCanvas');
    const velocityCanvas = document.getElementById('velocityCanvas');
    let currentResidue = DATA.rows[0]?.residue || '';

    document.getElementById('residueCount').textContent = DATA.summary.residueCount;
    document.getElementById('nearCount').textContent = DATA.summary.near10ACount;
    document.getElementById('traceCount').textContent = DATA.summary.traceCount;

    function fmt(value, digits = 4) {{
      if (value === null || value === undefined || Number.isNaN(value)) return 'n/a';
      return Number(value).toFixed(digits);
    }}

    function prepareCanvas(canvas) {{
      const ratio = Math.max(1, Math.min(3, window.devicePixelRatio || 1));
      const rect = canvas.getBoundingClientRect();
      const cssWidth = Math.max(320, Math.round(rect.width || canvas.clientWidth || 900));
      const cssHeight = Math.max(300, Math.round(rect.height || canvas.clientHeight || 350));
      const pixelWidth = Math.round(cssWidth * ratio);
      const pixelHeight = Math.round(cssHeight * ratio);
      if (canvas.width !== pixelWidth || canvas.height !== pixelHeight) {{
        canvas.width = pixelWidth;
        canvas.height = pixelHeight;
      }}
      const ctx = canvas.getContext('2d');
      ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
      ctx.imageSmoothingEnabled = true;
      ctx.imageSmoothingQuality = 'high';
      return {{ ctx, w: cssWidth, h: cssHeight }};
    }}

    function filteredRows() {{
      const q = searchEl.value.trim().toLowerCase();
      const f = filterEl.value;
      return DATA.rows.filter(row => {{
        const text = row.residue.toLowerCase() + ' ' + row.resname.toLowerCase() + row.resid;
        if (q && !text.includes(q)) return false;
        if (f === 'near' && !row.near10A) return false;
        if (f === 'pos' && !(row.rdc > 0)) return false;
        if (f === 'neg' && !(row.rdc < 0)) return false;
        return true;
      }});
    }}

    function sortRows(rows) {{
      const mode = sortEl.value;
      const copy = [...rows];
      const value = row => Number(row.rdc || 0);
      if (mode === 'rdc_desc') copy.sort((a, b) => value(b) - value(a));
      else if (mode === 'rdc_asc') copy.sort((a, b) => value(a) - value(b));
      else if (mode === 'min_distance') copy.sort((a, b) => Number(a.min || Infinity) - Number(b.min || Infinity));
      else if (mode === 'resid') copy.sort((a, b) => Number(a.resid || 0) - Number(b.resid || 0));
      else copy.sort((a, b) => Math.abs(value(b)) - Math.abs(value(a)));
      return copy;
    }}

    function renderTable() {{
      const rows = sortRows(filteredRows());
      rowsEl.innerHTML = rows.map(row => `
        <tr data-residue="${{row.residue}}" class="${{row.residue === currentResidue ? 'active' : ''}}">
          <td>${{row.residue}}${{row.near10A ? '<span class="tag">10 Å</span>' : ''}}</td>
          <td>${{fmt(row.rdc, 5)}}</td>
          <td>${{fmt(row.min, 2)}}</td>
          <td>${{fmt(row.delta, 2)}}</td>
        </tr>
      `).join('');
      rowsEl.querySelectorAll('tr').forEach(tr => {{
        tr.addEventListener('click', () => {{
          currentResidue = tr.dataset.residue;
          renderAll();
        }});
      }});
    }}

    function selectedRow() {{
      return DATA.rows.find(row => row.residue === currentResidue) || DATA.rows[0];
    }}

    function renderDetails(row) {{
      const metrics = [
        ['氨基酸', row.residue],
        ['RDC per ns', fmt(row.rdc, 6)],
        ['Slope Å/ns', fmt(row.slope, 5)],
        ['初始距离 Å', fmt(row.initial, 3)],
        ['最终距离 Å', fmt(row.final, 3)],
        ['平均距离 Å', fmt(row.mean, 3)],
        ['最小距离 Å', fmt(row.min, 3)],
        ['最大距离 Å', fmt(row.max, 3)],
        ['Δ 终点-起点 Å', fmt(row.delta, 3)],
        ['线性拟合 R²', fmt(row.r2, 3)],
      ];
      detailGrid.innerHTML = metrics.map(([label, value]) => `
        <div class="metric"><span>${{label}}</span><strong>${{value}}</strong></div>
      `).join('');
    }}

    function niceTicks(min, max, count = 5) {{
      if (!Number.isFinite(min) || !Number.isFinite(max) || min === max) return [min || 0];
      const span = max - min;
      const raw = span / Math.max(1, count - 1);
      const magnitude = Math.pow(10, Math.floor(Math.log10(raw)));
      const normalized = raw / magnitude;
      const step = (normalized <= 1 ? 1 : normalized <= 2 ? 2 : normalized <= 5 ? 5 : 10) * magnitude;
      const start = Math.ceil(min / step) * step;
      const ticks = [];
      for (let v = start; v <= max + step * 0.5; v += step) ticks.push(Number(v.toFixed(8)));
      return ticks;
    }}

    function drawAxes(ctx, w, h, xLabel, yLabel, xMin, xMax, yMin, yMax, xMap, yMap, yTickDigits = 1) {{
      const x0 = 60;
      const y0 = h - 50;
      const x1 = w - 22;
      const y1 = 22;
      ctx.strokeStyle = '#e6eaf0';
      ctx.lineWidth = 1;
      ctx.font = '12px Arial';
      ctx.fillStyle = '#667085';
      const xTicks = niceTicks(xMin, xMax, 5);
      const yTicks = niceTicks(yMin, yMax, 5);
      xTicks.forEach(t => {{
        const x = xMap(t);
        ctx.beginPath();
        ctx.moveTo(x, y1);
        ctx.lineTo(x, y0);
        ctx.stroke();
        ctx.fillText(String(Number(t.toFixed(2))), x - 10, y0 + 18);
      }});
      yTicks.forEach(t => {{
        const y = yMap(t);
        ctx.beginPath();
        ctx.moveTo(x0, y);
        ctx.lineTo(x1, y);
        ctx.stroke();
        ctx.fillText(Number(t).toFixed(yTickDigits), 18, y + 4);
      }});
      ctx.strokeStyle = '#9aa4b2';
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(x0, y1);
      ctx.lineTo(x0, y0);
      ctx.lineTo(x1, y0);
      ctx.stroke();
      ctx.fillStyle = '#667085';
      ctx.font = '13px Arial';
      ctx.fillText(xLabel, w / 2 - 35, h - 12);
      ctx.save();
      ctx.translate(14, h / 2 + 46);
      ctx.rotate(-Math.PI / 2);
      ctx.fillText(yLabel, 0, 0);
      ctx.restore();
    }}

    function drawTrace(row) {{
      const {{ ctx, w, h }} = prepareCanvas(distanceCanvas);
      ctx.clearRect(0, 0, w, h);
      const trace = DATA.traces[row.residue];
      document.getElementById('distanceTitle').textContent = `${{row.residue}} 距离轨迹`;
      document.getElementById('distanceNote').innerHTML = `纵轴距离 d(t) 表示 <strong>${{row.residue}}</strong> 所有原子到底物 LIG 所有原子的最近原子-原子距离，单位 Å；它不是 Cα 距离，也不是残基质心到底物质心距离。`;
      if (!trace) {{
        ctx.fillStyle = '#667085';
        ctx.fillText('这个氨基酸没有可用的距离轨迹。', 80, 80);
        return;
      }}
      const xs = DATA.times;
      const ys = trace.filter(v => Number.isFinite(v));
      const xmin = Math.min(...xs), xmax = Math.max(...xs);
      const ymin = Math.min(...ys), ymax = Math.max(...ys);
      const ypad = Math.max(0.5, (ymax - ymin) * 0.12);
      const yMinPlot = ymin - ypad;
      const yMaxPlot = ymax + ypad;
      const xmap = x => 60 + (x - xmin) / (xmax - xmin) * (w - 82);
      const ymap = y => (h - 50) - (y - yMinPlot) / (yMaxPlot - yMinPlot) * (h - 72);
      drawAxes(ctx, w, h, '时间 (ns)', '最近距离 d(t), Å', xmin, xmax, yMinPlot, yMaxPlot, xmap, ymap, 1);
      ctx.strokeStyle = row.rdc >= 0 ? '#2f5f8f' : '#287a5b';
      ctx.lineWidth = 2.3;
      ctx.beginPath();
      trace.forEach((y, i) => {{
        const x = xmap(xs[i]);
        const yy = ymap(y);
        if (i === 0) ctx.moveTo(x, yy); else ctx.lineTo(x, yy);
      }});
      ctx.stroke();
      ctx.fillStyle = '#20252d';
      ctx.font = '12px Arial';
      ctx.fillText(`最小 ${{fmt(ymin, 2)}} Å   最大 ${{fmt(ymax, 2)}} Å`, 70, 38);
    }}

    function velocities(row) {{
      const trace = DATA.traces[row.residue];
      if (!trace) return [];
      const out = [];
      for (let i = 0; i < trace.length - 1; i++) {{
        const dtPs = (DATA.times[i + 1] - DATA.times[i]) * 1000;
        if (dtPs > 0) out.push(((trace[i + 1] - trace[i]) / dtPs) / 0.01);
      }}
      return out;
    }}

    function drawVelocity(row) {{
      const {{ ctx, w, h }} = prepareCanvas(velocityCanvas);
      ctx.clearRect(0, 0, w, h);
      document.getElementById('velocityTitle').textContent = `${{row.residue}} 距离波动速度分布`;
      document.getElementById('velocityNote').innerHTML = `每个柱子来自相邻两帧的速度 V = Δd/Δt；负值表示该时间间隔内 <strong>${{row.residue}}</strong> 靠近底物，正值表示远离底物。`;
      const vals = velocities(row);
      if (!vals.length) {{
        ctx.fillStyle = '#667085';
        ctx.fillText('这个氨基酸没有可用的速度分布。', 80, 80);
        return;
      }}
      const maxX = Math.max(20, Math.ceil(Math.max(...vals.map(Math.abs))));
      const bins = 54;
      const counts = new Array(bins).fill(0);
      vals.forEach(v => {{
        const idx = Math.max(0, Math.min(bins - 1, Math.floor((v + maxX) / (2 * maxX) * bins)));
        counts[idx]++;
      }});
      const maxCount = Math.max(...counts);
      const x0 = 60, y0 = h - 50, pw = w - 82, ph = h - 72;
      const xmap = x => x0 + (x + maxX) / (2 * maxX) * pw;
      const ymap = y => y0 - y / maxCount * ph;
      drawAxes(ctx, w, h, 'V × 10⁻² (Å/ps)', '计数', -maxX, maxX, 0, maxCount, xmap, ymap, 0);
      counts.forEach((c, i) => {{
        const x = x0 + i / bins * pw;
        const bw = Math.max(2, pw / bins - 1.2);
        const bh = c / maxCount * ph;
        const center = -maxX + (i + 0.5) / bins * (2 * maxX);
        const strength = Math.min(1, Math.abs(center) / Math.max(DATA.velocityThreshold, 1));
        if (Math.abs(center) <= DATA.velocityThreshold) {{
          const light = Math.round(225 - 35 * strength);
          ctx.fillStyle = `rgb(${{light}}, ${{light}}, ${{light - 5}})`;
        }} else if (center < 0) {{
          const g = Math.round(145 - 45 * strength);
          const b = Math.round(185 - 55 * strength);
          ctx.fillStyle = `rgb(50, ${{g}}, ${{b}})`;
        }} else {{
          const r = Math.round(210 + 30 * strength);
          const g = Math.round(135 - 60 * strength);
          ctx.fillStyle = `rgb(${{r}}, ${{g}}, 70)`;
        }}
        ctx.fillRect(x, y0 - bh, bw, bh);
      }});
      const threshold = DATA.velocityThreshold;
      [-threshold, threshold, 0].forEach(v => {{
        const x = xmap(v);
        ctx.strokeStyle = v === 0 ? '#20252d' : '#a94040';
        ctx.setLineDash(v === 0 ? [] : [5, 4]);
        ctx.beginPath();
        ctx.moveTo(x, 20);
        ctx.lineTo(x, y0);
        ctx.stroke();
      }});
      ctx.setLineDash([]);
      ctx.fillStyle = '#20252d';
      ctx.font = '12px Arial';
      const approaching = vals.filter(v => v < 0).length / vals.length;
      ctx.fillText(`n=${{vals.length}} 个速度间隔；靠近底物=${{(approaching * 100).toFixed(1)}}%`, 70, 38);
      ctx.fillStyle = '#667085';
      ctx.fillText('绿色/蓝色：距离缩短；灰色：接近参考带；橙红色：距离增加', 70, 58);
    }}

    function renderAll() {{
      renderTable();
      const row = selectedRow();
      renderDetails(row);
      drawTrace(row);
      drawVelocity(row);
    }}

    searchEl.addEventListener('input', renderTable);
    filterEl.addEventListener('change', renderTable);
    sortEl.addEventListener('change', renderTable);
    renderAll();
  </script>
</body>
</html>
"""
    out.write_text(html_text, encoding="utf-8")
    print(out)
    print(f"rows={len(rows)} traces={len(traces)} times={len(times)}")


if __name__ == "__main__":
    main()
