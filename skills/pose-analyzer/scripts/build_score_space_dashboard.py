#!/usr/bin/env python
"""Build a localhost-friendly molecule space dashboard from score_space_model outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from pathlib import Path


DEFAULT_MODEL_DIR = Path("~/vina_task2/score_space_model")
DEFAULT_DASHBOARD_DIR = Path("~/vina_task2/score_space_dashboard")


def read_csv(path: Path) -> list[dict[str, str]]:
    for encoding in ("utf-8-sig", "utf-8", "gb18030", "cp936", "latin1"):
        try:
            with path.open("r", encoding=encoding, newline="") as handle:
                return list(csv.DictReader(handle))
        except UnicodeDecodeError:
            continue
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle))


def safe_float(value: object) -> float:
    try:
        if value in ("", None):
            return math.nan
        return float(value)
    except Exception:
        return math.nan


def compact_row(row: dict[str, str]) -> dict[str, object]:
    numeric_fields = [
        "official_binding_score",
        "analysis_score",
        "analysis_raw_score",
        "predicted_binding_score",
        "chembl_scaffold_similarity",
        "pure_smiles_cosine1",
        "pure_smiles_cosine2",
        "pure_smiles_pc1",
        "pure_smiles_pc2",
        "structural_interaction_pc1",
        "structural_interaction_pc2",
        "affinity_kcal_mol",
        "inner_rmsd",
        "cnn_pose_score",
        "hbond_count",
        "hydrophobic_count",
        "vdw_contact_count",
        "pi_contact_count",
        "ch_pi_count",
        "qed",
        "mw",
        "logp",
        "hbd",
        "hba",
        "tpsa",
        "rot_bonds",
        "heavy_atoms",
        "aromatic_rings",
        "murcko_heavy_atoms",
        "druglike_refinement_score",
        "prediction_residual",
        "active_learning_rank",
        "frontier_score",
        "novelty_score",
        "uncertainty_proxy",
        "quality_score",
        "druglike_score",
        "contact_score",
        "pocket_atom_count_5a",
        "pocket_residue_count_5a",
        "pocket_rg_5a",
        "pocket_span_5a",
        "pocket_centroid_distance_5a",
        "pocket_hydrophobic_fraction_5a",
        "pocket_polar_fraction_5a",
        "pocket_charged_fraction_5a",
        "pocket_aromatic_fraction_5a",
    ]
    out: dict[str, object] = {
        "seq_id": row.get("seq_id", ""),
        "nickname": row.get("nickname", ""),
        "score_set": row.get("score_set", ""),
        "canonical_smiles": row.get("canonical_smiles", ""),
        "smiles_sanity_status": row.get("smiles_sanity_status", ""),
        "smiles_sanity_reasons": row.get("smiles_sanity_reasons", ""),
        "murcko_scaffold": row.get("murcko_scaffold", ""),
        "analysis_score_source": row.get("analysis_score_source", ""),
    }
    for field in ("selection_reason", "selection_origin"):
        if field in row:
            out[field] = row.get(field, "")
    for field in numeric_fields:
        if field in {"selection_reason", "selection_origin"}:
            continue
        value = safe_float(row.get(field))
        out[field] = None if math.isnan(value) else value
    return out


def build_html() -> str:
    return r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Molecule Space Dashboard</title>
  <link rel="stylesheet" href="static/style.css" />
</head>
<body>
  <header>
    <h1>Molecule Space Dashboard</h1>
    <div class="subtle" id="dashboard-subtitle">Pure-SMILES-space + structural-interaction-space for clustering, dimensionality reduction, and docking-pose triage.</div>
  </header>

  <main>
    <section class="diagnostics">
      <div>
        <h2>Teacher Signal Check</h2>
        <div id="metric-cards" class="metric-grid"></div>
        <div id="model-warning" class="warning">Loading model diagnostics...</div>
      </div>
      <div>
        <h2 id="quality-title">Predicted vs Teacher</h2>
        <svg id="quality-plot" viewBox="0 0 420 260" aria-label="predicted versus teacher score"></svg>
      </div>
      <div>
        <h2>Feature Importance</h2>
        <div id="feature-importance" class="feature-list subtle">Loading...</div>
      </div>
    </section>

    <section class="toolbar">
      <label>Space
        <select id="space-select">
          <option value="pure">Pure-SMILES-space (cosine)</option>
          <option value="pure_pca">Pure-SMILES PCA (legacy)</option>
          <option value="structural">structural-interaction-space</option>
        </select>
      </label>
      <label>Color
        <select id="color-select">
          <option value="affinity_kcal_mol">GNINA/Vina affinity</option>
          <option value="cnn_pose_score">CNN pose score</option>
          <option value="official_binding_score">triage score</option>
          <option value="chembl_scaffold_similarity">ChEMBL scaffold similarity</option>
          <option value="qed">QED</option>
          <option value="mw">MW</option>
          <option value="predicted_binding_score">predicted score (experimental)</option>
        </select>
      </label>
      <label>Search
        <input id="search-box" placeholder="seq_id, nickname, SMILES" />
      </label>
      <label>
        <input id="scored-only" type="checkbox" />
        scored only
      </label>
      <button id="reset-view">Reset View</button>
      <label>Min
        <input id="color-min" type="number" step="0.001" />
      </label>
      <label>Max
        <input id="color-max" type="number" step="0.001" />
      </label>
      <label>
        <input id="color-reverse" type="checkbox" />
        reverse color
      </label>
      <button id="reset-color">Reset Color</button>
      <div id="color-legend" class="color-legend" aria-label="color legend"></div>
    </section>

    <section class="layout">
      <div class="plot-panel">
        <svg id="space-plot" viewBox="0 0 980 680" aria-label="molecule space plot"></svg>
      </div>
      <aside id="detail-panel">
        <h2>Selected Molecule</h2>
        <div id="detail-content" class="subtle">Click a point.</div>
      </aside>
    </section>

    <section class="summary">
      <h2>Human Triage Table</h2>
      <div class="subtle" id="table-note">Sorted by teacher score when present, otherwise by ChEMBL scaffold similarity, QED, and docking context. Predicted score is diagnostic context only.</div>
      <table id="candidate-table">
        <thead>
          <tr>
            <th>seq_id</th>
            <th id="teacher-th">teacher</th>
            <th>pred.</th>
            <th>resid.</th>
            <th>QED</th>
            <th>MW</th>
            <th>ChEMBL scaffold</th>
            <th>affinity</th>
            <th>CNN pose</th>
          </tr>
        </thead>
        <tbody></tbody>
      </table>
    </section>
  </main>

  <script src="static/app.js"></script>
</body>
</html>
"""


def build_css() -> str:
    return r"""body {
  margin: 0;
  font-family: Arial, Helvetica, sans-serif;
  color: #172033;
  background: #f7f8fb;
}
header {
  padding: 18px 24px 12px;
  background: #ffffff;
  border-bottom: 1px solid #d8dee9;
}
h1, h2 { margin: 0 0 8px; }
main { padding: 16px 24px 28px; }
.subtle { color: #607087; font-size: 13px; }
.diagnostics {
  display: grid;
  grid-template-columns: minmax(260px, 1fr) 420px minmax(220px, 0.8fr);
  gap: 16px;
  margin-bottom: 16px;
}
.diagnostics > div {
  background: #ffffff;
  border: 1px solid #d8dee9;
  border-radius: 6px;
  padding: 14px;
}
.metric-grid {
  display: grid;
  grid-template-columns: repeat(2, minmax(90px, 1fr));
  gap: 8px;
  margin: 8px 0 10px;
}
.metric-card {
  border: 1px solid #e3e8f0;
  border-radius: 4px;
  padding: 8px;
  background: #fafbfe;
}
.metric-card .label { color: #607087; font-size: 11px; }
.metric-card .value { font-size: 18px; margin-top: 2px; }
.warning {
  padding: 9px 10px;
  border-radius: 4px;
  border: 1px solid #f0c36d;
  background: #fff8e5;
  color: #604600;
  font-size: 13px;
}
.warning.good {
  border-color: #9bd29b;
  background: #effaf0;
  color: #205522;
}
.feature-row {
  display: grid;
  grid-template-columns: 1fr 58px;
  gap: 8px;
  margin: 5px 0;
}
.feature-bar {
  height: 7px;
  border-radius: 4px;
  background: #d7e3f7;
  margin-top: 3px;
  overflow: hidden;
}
.feature-bar span {
  display: block;
  height: 100%;
  background: #547aa5;
}
.toolbar {
  display: flex;
  flex-wrap: wrap;
  align-items: end;
  gap: 12px;
  padding: 12px;
  background: #ffffff;
  border: 1px solid #d8dee9;
  border-radius: 6px;
}
label { display: grid; gap: 4px; font-size: 12px; color: #3f4f65; }
select, input, button {
  height: 32px;
  border: 1px solid #bcc6d4;
  border-radius: 4px;
  background: #ffffff;
  color: #172033;
  padding: 0 9px;
}
button { cursor: pointer; }
.color-legend {
  min-width: 260px;
  align-self: stretch;
  display: grid;
  align-content: center;
  gap: 5px;
  padding-left: 4px;
}
.legend-title { font-size: 12px; color: #3f4f65; }
.legend-bar {
  height: 12px;
  border-radius: 3px;
  border: 1px solid #cbd5e1;
  overflow: hidden;
  display: flex;
}
.legend-stop { flex: 1; }
.legend-labels {
  display: flex;
  justify-content: space-between;
  color: #607087;
  font-size: 11px;
}
.layout {
  display: grid;
  grid-template-columns: minmax(420px, 1fr) 360px;
  gap: 16px;
  margin-top: 16px;
}
.plot-panel, #detail-panel, .summary {
  background: #ffffff;
  border: 1px solid #d8dee9;
  border-radius: 6px;
}
.plot-panel { min-height: 690px; overflow: hidden; }
#space-plot { width: 100%; height: 690px; display: block; }
#detail-panel { padding: 16px; overflow-wrap: anywhere; }
.detail-grid {
  display: grid;
  grid-template-columns: 132px 1fr;
  gap: 7px 10px;
  font-size: 13px;
}
.key { color: #607087; }
.summary { margin-top: 16px; padding: 16px; }
table { width: 100%; border-collapse: collapse; font-size: 13px; }
th, td { border-bottom: 1px solid #e3e8f0; padding: 7px 8px; text-align: left; }
th { color: #3f4f65; background: #f8fafc; }
.point { cursor: pointer; }
.recommendation-ring {
  fill: none;
  stroke: #f5b700;
  stroke-width: 3.2;
  opacity: 0.95;
  filter: drop-shadow(0 0 5px rgba(245, 183, 0, 0.85));
  pointer-events: none;
}
.frontier-ring {
  fill: none;
  stroke: #2dd4bf;
  stroke-width: 2;
  opacity: 0.72;
  pointer-events: none;
}
.axis-label { fill: #607087; font-size: 12px; }
.quality-axis { stroke: #c6d0df; stroke-width: 1; }
.quality-point { opacity: 0.82; }
@media (max-width: 900px) {
  .diagnostics { grid-template-columns: 1fr; }
  .layout { grid-template-columns: 1fr; }
  #detail-panel { min-height: 240px; }
}
"""


def build_js() -> str:
    return r"""let DATA = [];
let METRICS = {};
let ACTIVE = {frontier_seq_ids: [], official_score_seq_ids: [], frontier_seeds: [], official_score_recommendations: []};
let userColorRange = {};
let selectedId = null;
let view = {scale: 1, dx: 0, dy: 0};
const svg = document.getElementById("space-plot");
const qualitySvg = document.getElementById("quality-plot");
const detail = document.getElementById("detail-content");
const tableBody = document.querySelector("#candidate-table tbody");

function num(v) { return typeof v === "number" && Number.isFinite(v); }
function fmt(v, n=3) { return num(v) ? v.toFixed(n) : ""; }
function teacherLabel() { return (METRICS && METRICS.score_label) || "official binding score"; }
function unsupervisedMode() { return (METRICS && METRICS.analysis_mode) === "unsupervised"; }
function teacherShortLabel() {
  const label = teacherLabel();
  return label.length > 22 ? "teacher" : label;
}
function officialRecommendationIds() {
  return (ACTIVE.official_score_seq_ids || []).slice(0, 3);
}
function officialRecommendationRows() {
  const allowed = new Set(officialRecommendationIds());
  return (ACTIVE.official_score_recommendations || []).filter(r => allowed.has(r.seq_id)).slice(0, 3);
}
function activeInfo(seqId) {
  const official = officialRecommendationRows().find(r => r.seq_id === seqId);
  const frontier = (ACTIVE.frontier_seeds || []).find(r => r.seq_id === seqId);
  return official || frontier || null;
}
function isOfficialRecommendation(seqId) { return officialRecommendationIds().includes(seqId); }
function isFrontierSeed(seqId) { return (ACTIVE.frontier_seq_ids || []).includes(seqId); }

function colorScale(value, field) {
  if (!num(value)) return "#9aa0a6";
  if (field === "mw") {
    if (value < 250 || value > 600) return "#e15759";
    if (value >= 300 && value <= 500) return "#59a14f";
    return "#edc948";
  }
  const range = colorRange(field);
  const u = userColorRange[field];
  const minV = (u && num(u.min)) ? u.min : range.min;
  const maxV = (u && num(u.max)) ? u.max : range.max;
  const reverse = (u && u.reverse !== undefined) ? u.reverse : (range.reverse || false);
  return continuousColor(value, minV, maxV, reverse);
}

function colorRange(field) {
  const vals = DATA.map(r => r[field]).filter(num).sort((a,b) => a-b);
  if (!vals.length) return {min: 0, max: 1};
  if (field === "official_binding_score") {
    return {min: Math.min(0, vals[0]), max: Math.max(0.4, vals[vals.length - 1])};
  }
  if (field === "affinity_kcal_mol") {
    return {min: vals[0], max: vals[vals.length - 1], reverse: true};
  }
  return {min: vals[0], max: vals[vals.length - 1]};
}

function continuousColor(value, minV, maxV, reverse=false) {
  let t = Math.max(0, Math.min(1, (value - minV) / ((maxV - minV) || 1)));
  if (reverse) t = 1 - t;
  const stops = [
    [231, 76, 60],
    [246, 189, 96],
    [91, 155, 213],
    [75, 192, 125],
  ];
  const scaled = t * (stops.length - 1);
  const i = Math.min(stops.length - 2, Math.floor(scaled));
  const f = scaled - i;
  const c = stops[i].map((v, idx) => Math.round(v + (stops[i + 1][idx] - v) * f));
  return `rgb(${c[0]},${c[1]},${c[2]})`;
}

function legendFor(field) {
  if (field === "mw") {
    return {
      title: "Molecular Weight",
      colors: ["#e15759", "#edc948", "#59a14f", "#edc948", "#e15759"],
      labels: ["<<250", "250-300", "300-500", "500-600", ">600"],
      note: "Optimal oral-drug window (300-500) is green."
    };
  }
  const title = {
    chembl_scaffold_similarity: "ChEMBL scaffold similarity",
    qed: "QED",
    official_binding_score: teacherLabel(),
    cnn_pose_score: "CNN pose score",
    predicted_binding_score: "predicted score (experimental)",
    affinity_kcal_mol: "GNINA/Vina affinity (kcal/mol)"
  }[field] || field;
  const range = colorRange(field);
  const u = userColorRange[field];
  const minV = (u && num(u.min)) ? u.min : range.min;
  const maxV = (u && num(u.max)) ? u.max : range.max;
  const reverse = (u && u.reverse !== undefined) ? u.reverse : (range.reverse || false);
  const colors = reverse
    ? ["rgb(75,192,125)", "rgb(91,155,213)", "rgb(246,189,96)", "rgb(231,76,60)"]
    : ["rgb(231,76,60)", "rgb(246,189,96)", "rgb(91,155,213)", "rgb(75,192,125)"];
  return {
    title,
    colors,
    labels: [fmt(minV), "", "", fmt(maxV)],
    note: field === "predicted_binding_score"
      ? "Continuous scale from current predicted-score distribution; low R^2 means diagnostic only."
      : "Continuous scale from current data range; gray means missing value."
  };
}

function drawLegend() {
  const field = document.getElementById("color-select").value;
  const legend = legendFor(field);
  const target = document.getElementById("color-legend");
  target.innerHTML = `<div class="legend-title">${legend.title}</div>
    <div class="legend-bar">${legend.colors.map(color => `<span class="legend-stop" style="background:${color}"></span>`).join("")}</div>
    <div class="legend-labels">${legend.labels.map(label => `<span>${label}</span>`).join("")}</div>
    <div class="subtle">${legend.note}</div>`;
}

function syncColorControls() {
  const field = document.getElementById("color-select").value;
  const range = colorRange(field);
  const u = userColorRange[field] || {};
  const minInput = document.getElementById("color-min");
  const maxInput = document.getElementById("color-max");
  const reverseInput = document.getElementById("color-reverse");
  minInput.value = num(u.min) ? u.min : range.min;
  maxInput.value = num(u.max) ? u.max : range.max;
  reverseInput.checked = u.reverse !== undefined ? !!u.reverse : !!range.reverse;
}

function updateColorOverride() {
  const field = document.getElementById("color-select").value;
  const minV = Number(document.getElementById("color-min").value);
  const maxV = Number(document.getElementById("color-max").value);
  userColorRange[field] = {
    min: Number.isFinite(minV) ? minV : undefined,
    max: Number.isFinite(maxV) ? maxV : undefined,
    reverse: document.getElementById("color-reverse").checked,
  };
  draw();
}

function currentRows() {
  const query = document.getElementById("search-box").value.toLowerCase().trim();
  const scoredOnly = document.getElementById("scored-only").checked;
  return DATA.filter(row => {
    if (scoredOnly && !num(row.official_binding_score)) return false;
    if (!query) return true;
    return [row.seq_id, row.nickname, row.canonical_smiles, row.murcko_scaffold]
      .some(v => String(v || "").toLowerCase().includes(query));
  });
}

function coords(row) {
  const space = document.getElementById("space-select").value;
  if (space === "structural") return [row.structural_interaction_pc1, row.structural_interaction_pc2];
  if (space === "pure_pca") return [row.pure_smiles_pc1, row.pure_smiles_pc2];
  if (!num(row.pure_smiles_cosine1) || !num(row.pure_smiles_cosine2)) return [row.pure_smiles_pc1, row.pure_smiles_pc2];
  return [row.pure_smiles_cosine1, row.pure_smiles_cosine2];
}

function draw() {
  drawLegend();
  const rows = currentRows();
  const colorField = document.getElementById("color-select").value;
  const points = rows.map(r => [r, ...coords(r)]).filter(([,x,y]) => num(x) && num(y));
  const xs = points.map(p => p[1]);
  const ys = points.map(p => p[2]);
  if (!points.length) {
    svg.innerHTML = `<rect x="0" y="0" width="980" height="680" fill="#fff"/><text x="44" y="44" class="axis-label">No matching points.</text>`;
    tableBody.innerHTML = "";
    return;
  }
  const minX = Math.min(...xs), maxX = Math.max(...xs);
  const minY = Math.min(...ys), maxY = Math.max(...ys);
  const W = 980, H = 680, P = 44;
  function sx(x) { return P + ((x - minX) / ((maxX - minX) || 1)) * (W - 2*P) * view.scale + view.dx; }
  function sy(y) { return H - P - ((y - minY) / ((maxY - minY) || 1)) * (H - 2*P) * view.scale + view.dy; }
  svg.innerHTML = `<rect x="0" y="0" width="${W}" height="${H}" fill="#fff"/>
    <text x="${P}" y="24" class="axis-label">${document.getElementById("space-select").selectedOptions[0].text}</text>
    <text x="${P + 250}" y="24" class="axis-label">gold outer ring = top teacher-score picks (default <=3); teal ring = frontier seed</text>`;
  for (const [row, x, y] of points) {
    const scored = num(row.official_binding_score);
    const r = scored ? 6 : 3.2;
    const stroke = row.seq_id === selectedId ? "#111827" : (scored ? "#253145" : "none");
    const sw = row.seq_id === selectedId ? 2.4 : 1.1;
    const circle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
    if (isOfficialRecommendation(row.seq_id) || isFrontierSeed(row.seq_id)) {
      const ring = document.createElementNS("http://www.w3.org/2000/svg", "circle");
      ring.setAttribute("cx", sx(x));
      ring.setAttribute("cy", sy(y));
      ring.setAttribute("r", isOfficialRecommendation(row.seq_id) ? r + 7 : r + 4.5);
      ring.classList.add(isOfficialRecommendation(row.seq_id) ? "recommendation-ring" : "frontier-ring");
      svg.appendChild(ring);
    }
    circle.setAttribute("cx", sx(x));
    circle.setAttribute("cy", sy(y));
    circle.setAttribute("r", r);
    circle.setAttribute("fill", colorScale(row[colorField], colorField));
    circle.setAttribute("stroke", stroke);
    circle.setAttribute("stroke-width", sw);
    circle.setAttribute("opacity", scored ? "0.9" : "0.56");
    circle.classList.add("point");
    circle.addEventListener("click", () => selectRow(row.seq_id));
    const title = document.createElementNS("http://www.w3.org/2000/svg", "title");
    title.textContent = `${row.seq_id} teacher=${fmt(row.official_binding_score)} pred=${fmt(row.predicted_binding_score)} QED=${fmt(row.qed)} MW=${fmt(row.mw)} ${row.canonical_smiles}`;
    circle.appendChild(title);
    svg.appendChild(circle);
  }
  drawTable(rows);
}

function renderDiagnostics() {
  const label = teacherLabel();
  document.getElementById("dashboard-subtitle").textContent = unsupervisedMode()
    ? "Pure-SMILES-space + structural-interaction-space for clustering, dimensionality reduction, affinity, CNN pose score, and contact triage."
    : "Pure-SMILES-space + structural-interaction-space. The selected teacher signal drives model diagnostics; prediction is diagnostic only when model quality is weak.";
  document.getElementById("quality-title").textContent = unsupervisedMode() ? "Affinity vs CNN Pose" : `Predicted vs ${teacherShortLabel()}`;
  document.getElementById("teacher-th").textContent = teacherShortLabel();
  document.getElementById("table-note").textContent = `Sorted by ${label} when present, otherwise by ChEMBL scaffold similarity, QED, and docking context. Predicted score is diagnostic context only.`;
  const cards = document.getElementById("metric-cards");
  const warning = document.getElementById("model-warning");
  const r2 = METRICS.test_r2;
  const mae = METRICS.test_mae;
  const cardData = [
    [unsupervisedMode() ? "rows" : "scored", unsupervisedMode() ? METRICS.space_rows : METRICS.scored_rows],
    ["train", unsupervisedMode() ? "" : METRICS.train_rows],
    ["test", unsupervisedMode() ? "" : METRICS.test_rows],
    ["test R^2", num(r2) ? r2.toFixed(3) : ""],
    ["test MAE", num(mae) ? mae.toFixed(4) : ""],
    ["CV R^2", num(METRICS.cv_r2_mean) ? `${METRICS.cv_r2_mean.toFixed(3)} +/- ${(METRICS.cv_r2_std || 0).toFixed(3)}` : ""],
    ["CV MAE", num(METRICS.cv_mae_mean) ? `${METRICS.cv_mae_mean.toFixed(4)} +/- ${(METRICS.cv_mae_std || 0).toFixed(4)}` : ""],
    ["status", METRICS.model_interpretation || ""],
  ];
  cards.innerHTML = cardData.map(([label, value]) => `<div class="metric-card"><div class="label">${label}</div><div class="value">${value ?? ""}</div></div>`).join("");
  const poor = !num(r2) || r2 < 0.3 || !num(METRICS.cv_r2_mean) || METRICS.cv_r2_mean < 0.3;
  warning.className = poor ? "warning" : "warning good";
  warning.textContent = unsupervisedMode()
    ? "No supervised binding-score model was trained. Use the point clouds, affinity, CNN pose score, inner RMSD, contacts, and descriptors as clustering/triage signals."
    : poor
    ? `Prediction is not reliable here. R^2=${num(r2) ? r2.toFixed(3) : "NA"} means this model should be treated as a failure diagnostic, not a binding-score predictor.`
    : `Prediction has a usable screening signal, but still needs external validation. R^2=${r2.toFixed(3)}.`;
  drawQualityPlot();
  drawFeatureImportance();
}

function drawQualityPlot() {
  const rows = unsupervisedMode()
    ? DATA.filter(r => num(r.affinity_kcal_mol) && num(r.cnn_pose_score))
    : DATA.filter(r => num(r.official_binding_score) && num(r.predicted_binding_score));
  const W = 420, H = 260, P = 38;
  qualitySvg.innerHTML = `<rect x="0" y="0" width="${W}" height="${H}" fill="#fff"/>`;
  if (!rows.length) {
    qualitySvg.innerHTML += `<text x="${P}" y="${P}" class="axis-label">No paired affinity/CNN pose data.</text>`;
    return;
  }
  const values = unsupervisedMode()
    ? rows.flatMap(r => [r.cnn_pose_score, -r.affinity_kcal_mol])
    : rows.flatMap(r => [r.official_binding_score, r.predicted_binding_score]);
  const minV = Math.min(...values), maxV = Math.max(...values);
  function s(v) { return P + ((v - minV) / ((maxV - minV) || 1)) * (W - 2 * P); }
  function y(v) { return H - P - ((v - minV) / ((maxV - minV) || 1)) * (H - 2 * P); }
  qualitySvg.innerHTML += `<line x1="${P}" y1="${H-P}" x2="${W-P}" y2="${P}" class="quality-axis"/>
    <text x="${P}" y="${H-8}" class="axis-label">${unsupervisedMode() ? "CNN pose" : "predicted"}</text>
    <text x="6" y="${P}" class="axis-label">${unsupervisedMode() ? "-affinity" : teacherShortLabel()}</text>`;
  for (const row of rows) {
    const c = document.createElementNS("http://www.w3.org/2000/svg", "circle");
    c.setAttribute("cx", s(unsupervisedMode() ? row.cnn_pose_score : row.predicted_binding_score));
    c.setAttribute("cy", y(unsupervisedMode() ? -row.affinity_kcal_mol : row.official_binding_score));
    c.setAttribute("r", row.score_set === "test" ? 5.5 : 3.8);
    c.setAttribute("fill", row.score_set === "test" ? "#d95f02" : "#547aa5");
    c.setAttribute("stroke", "#263445");
    c.setAttribute("stroke-width", "0.8");
    c.classList.add("quality-point");
    c.addEventListener("click", () => selectRow(row.seq_id));
    const title = document.createElementNS("http://www.w3.org/2000/svg", "title");
    title.textContent = unsupervisedMode()
      ? `${row.seq_id} affinity=${fmt(row.affinity_kcal_mol,3)} cnn_pose=${fmt(row.cnn_pose_score,4)}`
      : `${row.seq_id} set=${row.score_set} teacher=${fmt(row.official_binding_score,4)} predicted=${fmt(row.predicted_binding_score,4)}`;
    c.appendChild(title);
    qualitySvg.appendChild(c);
  }
}

function drawFeatureImportance() {
  const target = document.getElementById("feature-importance");
  const items = Array.isArray(METRICS.feature_importance) ? METRICS.feature_importance.slice(0, 12) : [];
  if (!items.length) {
    target.textContent = "No feature importance available.";
    return;
  }
  const maxVal = Math.max(...items.map(i => i.importance || 0)) || 1;
  target.innerHTML = items.map(item => `<div class="feature-row">
    <div><div>${item.feature}</div><div class="feature-bar"><span style="width:${100 * (item.importance || 0) / maxVal}%"></span></div></div>
    <div>${fmt(item.importance, 3)}</div>
  </div>`).join("");
}

function selectRow(seqId) {
  selectedId = seqId;
  const row = DATA.find(r => r.seq_id === seqId);
  if (!row) return;
  const active = activeInfo(seqId) || {};
  const keys = [
    "seq_id", "nickname", "score_set", "official_binding_score", "analysis_raw_score", "analysis_score_source", "predicted_binding_score",
    "prediction_residual", "chembl_scaffold_similarity", "qed", "mw", "logp", "tpsa", "rot_bonds",
    "affinity_kcal_mol", "cnn_pose_score", "hbond_count", "hydrophobic_count", "pi_contact_count", "ch_pi_count",
    "inner_rmsd", "druglike_refinement_score",
    "pocket_atom_count_5a", "pocket_residue_count_5a", "pocket_rg_5a", "pocket_span_5a",
    "pocket_centroid_distance_5a", "pocket_hydrophobic_fraction_5a", "pocket_polar_fraction_5a",
    "pocket_charged_fraction_5a", "pocket_aromatic_fraction_5a",
    "active_learning_rank", "frontier_score", "novelty_score", "uncertainty_proxy",
    "selection_origin", "selection_reason",
    "smiles_sanity_status", "smiles_sanity_reasons",
    "pure_smiles_cosine1", "pure_smiles_cosine2", "pure_smiles_pc1", "pure_smiles_pc2", "structural_interaction_pc1", "structural_interaction_pc2",
    "murcko_scaffold", "canonical_smiles"
  ];
  const merged = {...row, ...active};
  const badge = isOfficialRecommendation(seqId) ? `<div class="warning good">Top teacher-score pick for inspection.</div>` : (isFrontierSeed(seqId) ? `<div class="warning good">Frontier seed for exploration/refinement.</div>` : "");
  detail.innerHTML = badge + `<div class="detail-grid">` + keys.map(k => {
    const v = merged[k];
    return `<div class="key">${k}</div><div>${num(v) ? fmt(v, 4) : (v ?? "")}</div>`;
  }).join("") + `</div>`;
  draw();
}

function drawTable(rows) {
  const sorted = [...rows].sort((a,b) => {
    const ar = isOfficialRecommendation(a.seq_id) ? 1 : 0;
    const br = isOfficialRecommendation(b.seq_id) ? 1 : 0;
    if (br !== ar) return br - ar;
    const ao = num(a.official_binding_score) ? a.official_binding_score : -Infinity;
    const bo = num(b.official_binding_score) ? b.official_binding_score : -Infinity;
    if (bo !== ao) return bo - ao;
    const aa = num(a.affinity_kcal_mol) ? a.affinity_kcal_mol : Infinity;
    const ba = num(b.affinity_kcal_mol) ? b.affinity_kcal_mol : Infinity;
    if (aa !== ba) return aa - ba;
    const ac = num(a.cnn_pose_score) ? a.cnn_pose_score : -Infinity;
    const bc = num(b.cnn_pose_score) ? b.cnn_pose_score : -Infinity;
    if (bc !== ac) return bc - ac;
    return (b.chembl_scaffold_similarity || 0) - (a.chembl_scaffold_similarity || 0);
  }).slice(0, 80);
  tableBody.innerHTML = sorted.map(row => `<tr data-id="${row.seq_id}">
    <td>${isOfficialRecommendation(row.seq_id) ? "[teacher] " : (isFrontierSeed(row.seq_id) ? "[frontier] " : "")}${row.seq_id}</td>
    <td>${fmt(row.official_binding_score)}</td>
    <td>${fmt(row.predicted_binding_score)}</td>
    <td>${num(row.official_binding_score) && num(row.predicted_binding_score) ? fmt(row.official_binding_score - row.predicted_binding_score) : ""}</td>
    <td>${fmt(row.qed)}</td>
    <td>${fmt(row.mw,1)}</td>
    <td>${fmt(row.chembl_scaffold_similarity)}</td>
    <td>${fmt(row.affinity_kcal_mol)}</td>
    <td>${fmt(row.cnn_pose_score)}</td>
  </tr>`).join("");
  tableBody.querySelectorAll("tr").forEach(tr => tr.addEventListener("click", () => selectRow(tr.dataset.id)));
}

svg.addEventListener("wheel", ev => {
  ev.preventDefault();
  const factor = ev.deltaY < 0 ? 1.12 : 0.9;
  view.scale = Math.max(0.4, Math.min(8, view.scale * factor));
  draw();
});
let dragging = false, last = null;
svg.addEventListener("mousedown", ev => { dragging = true; last = [ev.clientX, ev.clientY]; });
window.addEventListener("mouseup", () => dragging = false);
window.addEventListener("mousemove", ev => {
  if (!dragging || !last) return;
  view.dx += ev.clientX - last[0];
  view.dy += ev.clientY - last[1];
  last = [ev.clientX, ev.clientY];
  draw();
});

for (const id of ["space-select", "search-box", "scored-only"]) {
  document.getElementById(id).addEventListener("input", draw);
}
document.getElementById("color-select").addEventListener("input", () => { syncColorControls(); draw(); });
for (const id of ["color-min", "color-max", "color-reverse"]) {
  document.getElementById(id).addEventListener("input", updateColorOverride);
}
document.getElementById("reset-color").addEventListener("click", () => {
  delete userColorRange[document.getElementById("color-select").value];
  syncColorControls();
  draw();
});
document.getElementById("reset-view").addEventListener("click", () => { view = {scale:1, dx:0, dy:0}; draw(); });

Promise.all([
  fetch("data/molecules.json").then(r => r.json()),
  fetch("data/model_metrics.json").then(r => r.ok ? r.json() : {}).catch(() => ({})),
  fetch("data/active_learning_recommendations.json").then(r => r.ok ? r.json() : {frontier_seq_ids: [], official_score_seq_ids: []}).catch(() => ({frontier_seq_ids: [], official_score_seq_ids: []})),
]).then(([data, metrics, active]) => {
  DATA = data;
  METRICS = metrics || {};
  ACTIVE = active || {frontier_seq_ids: [], official_score_seq_ids: [], frontier_seeds: [], official_score_recommendations: []};
  ACTIVE.official_score_seq_ids = (ACTIVE.official_score_seq_ids && ACTIVE.official_score_seq_ids.length)
    ? ACTIVE.official_score_seq_ids.slice(0, 3)
    : ["__none__"];
  
  // 限制 official recommendations 最多 3 个（按 official binding score 降序，优先保留有 score 的）
  if (ACTIVE.official_score_seq_ids && ACTIVE.official_score_seq_ids.length) {
    const mapped = ACTIVE.official_score_seq_ids.map(id => DATA.find(r => r.seq_id === id)).filter(Boolean);
    const scored = mapped.filter(r => num(r.official_binding_score)).sort((a, b) => b.official_binding_score - a.official_binding_score);
    const unscored = mapped.filter(r => !num(r.official_binding_score));
    ACTIVE.official_score_seq_ids = scored.slice(0, 3).concat(unscored).slice(0, 3).map(r => r.seq_id);
  } else {
    const scored = DATA.filter(r => num(r.official_binding_score)).sort((a, b) => b.official_binding_score - a.official_binding_score);
    // If no external recommendations exist, highlight the top teacher-score rows.
    ACTIVE.official_score_seq_ids = scored.slice(0, 3).map(r => r.seq_id);
  }
  // 同步 recommendations 对象数组
  if (ACTIVE.official_score_recommendations) {
    ACTIVE.official_score_recommendations = ACTIVE.official_score_recommendations.filter(r => ACTIVE.official_score_seq_ids.includes(r.seq_id));
  }
  
  ACTIVE.official_score_seq_ids = (ACTIVE.official_score_seq_ids || []).filter(id => id !== "__none__").slice(0, 3);
  ACTIVE.official_score_recommendations = (ACTIVE.official_score_recommendations || [])
    .filter(r => ACTIVE.official_score_seq_ids.includes(r.seq_id))
    .slice(0, 3);
  syncColorControls();
  renderDiagnostics();
  draw();
});
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-matrix", type=Path, default=DEFAULT_MODEL_DIR / "binding_score_feature_matrix.csv")
    parser.add_argument("--metrics-json", type=Path, help="Optional binding_score_model_metrics.json; defaults to the feature matrix directory")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_DASHBOARD_DIR)
    args = parser.parse_args()

    matrix = args.feature_matrix.expanduser().resolve()
    metrics_path = (args.metrics_json.expanduser().resolve() if args.metrics_json else matrix.parent / "binding_score_model_metrics.json")
    rows = [compact_row(row) for row in read_csv(matrix)]
    outdir = args.outdir.expanduser().resolve()
    data_dir = outdir / "data"
    static_dir = outdir / "static"
    data_dir.mkdir(parents=True, exist_ok=True)
    static_dir.mkdir(parents=True, exist_ok=True)
    (outdir / "index.html").write_text(build_html(), encoding="utf-8")
    (static_dir / "style.css").write_text(build_css(), encoding="utf-8")
    (static_dir / "app.js").write_text(build_js(), encoding="utf-8")
    with (data_dir / "molecules.json").open("w", encoding="utf-8") as handle:
        json.dump(rows, handle, ensure_ascii=False)
    if metrics_path.exists():
        shutil.copyfile(metrics_path, data_dir / "model_metrics.json")
    else:
        (data_dir / "model_metrics.json").write_text("{}", encoding="utf-8")
    shutil.copyfile(matrix, data_dir / "binding_score_feature_matrix.csv")
    print(f"wrote={outdir} rows={len(rows)}")
    print(f"serve: cd {outdir} && python -m http.server 8765")
    print("open: http://127.0.0.1:8765")


if __name__ == "__main__":
    main()
