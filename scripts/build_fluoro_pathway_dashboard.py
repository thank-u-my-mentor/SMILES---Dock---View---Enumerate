#!/usr/bin/env python
"""Build a Metal-F-style dashboard and identity matrix from Table2.csv."""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
from pathlib import Path
from statistics import mean


AA_RE = re.compile(r"[^A-Z*]")


CLASS_INFO = {
    "PoxB-like": {
        "color": "#5DA5DA",
        "stage": "A",
        "module": "Quinone-linked oxidative decarboxylation",
        "reaction": "Pyruvate/2-oxo acid -> acetate-like product + CO2, electrons to quinone",
        "interpretation": "Closest to E. coli PoxB logic. Strong for asking whether fluoro-pyruvate-like substrates can be oxidatively decarboxylated into small fluorinated acid/acetyl equivalents.",
        "caution": "Many rows are annotation-only homologs; high identity to PoxB helps, but substrate scope still needs assay.",
    },
    "POX/SpxB-like": {
        "color": "#60BD68",
        "stage": "B",
        "module": "O2/phosphate pyruvate oxidase branch",
        "reaction": "Pyruvate + phosphate + O2 -> acetyl phosphate + CO2 + H2O2",
        "interpretation": "KEGG-like acetyl-phosphate branch. Potentially useful if the fluorinated pathway needs activated acyl phosphate chemistry or oxidative flux.",
        "caution": "H2O2 generation and oxygen coupling can complicate downstream enzyme compatibility.",
    },
    "PDC-like": {
        "color": "#F17CB0",
        "stage": "C",
        "module": "TPP decarboxylase aldehyde branch",
        "reaction": "Pyruvate/2-oxo acid -> aldehyde + CO2",
        "interpretation": "Classical pyruvate decarboxylase logic; useful as a decarboxylation branch from 2-oxo acids toward aldehydes.",
        "caution": "PDC family can be promiscuous, but fluorinated substrate acceptance is not guaranteed from EC alone.",
    },
    "Direct_RHEA54368": {
        "color": "#B276B2",
        "stage": "C+",
        "module": "Direct Rhea-supported 2-oxo-acid decarboxylase",
        "reaction": "Rhea 54368-associated 2-oxo-acid decarboxylation annotations",
        "interpretation": "These are still PDC/2ODC-like, but the table marks them as directly connected to the target Rhea reaction set.",
        "caution": "Very close isozyme pairs may be redundant for ordering.",
    },
    "AceE/PDH-E1-like": {
        "color": "#FAA43A",
        "stage": "D",
        "module": "C-C forming / PDH-E1-like TPP branch",
        "reaction": "Pyruvate-derived hydroxyethyl/acetolactate chemistry",
        "interpretation": "A pathway branch rather than a single clean pyruvate oxidase replacement; useful for thinking about condensation/C-C formation around fluorinated 2-oxo acids.",
        "caution": "Flagged as not a single-enzyme target reaction in the table; lower priority unless the desired route explicitly needs ALS/PDH-E1 chemistry.",
    },
}


PRIORITY_COLORS = {
    "positive_control": "#2F5F7F",
    "high": "#3E8F65",
    "medium": "#C17C21",
    "low": "#9B4A45",
}


IDENTITY_STOPS = [
    (0, (201, 76, 76)),
    (25, (236, 167, 106)),
    (50, (242, 231, 166)),
    (75, (155, 203, 156)),
    (90, (99, 182, 194)),
    (100, (63, 119, 181)),
]


def clean_seq(value: str) -> str:
    return AA_RE.sub("", (value or "").upper())


def short_gene(row: dict[str, str]) -> str:
    genes = (row.get("gene_names") or "").strip()
    return genes.split()[0] if genes else row["accession"]


def short_name(row: dict[str, str]) -> str:
    gene = short_gene(row)
    protein = row.get("protein_name", "")
    if gene and gene != row["accession"]:
        return gene
    return protein.split("(")[0].strip()[:36] or row["accession"]


def num(value: str, default: float = 0.0) -> float:
    try:
        if value in ("", "NA", None):
            return default
        return float(value)
    except Exception:
        return default


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = [row for row in csv.DictReader(handle)]
    out = []
    for idx, row in enumerate(rows, start=1):
        seq = clean_seq(row.get("amino_acid_sequence", ""))
        if not row.get("accession") or not seq:
            continue
        row = dict(row)
        row["rank"] = idx
        row["sequence"] = seq
        row["gene_label"] = short_gene(row)
        row["display_label"] = f"{row['gene_label']} | {row['accession']}"
        row["short_name"] = short_name(row)
        row["score"] = num(row.get("total_balanced_score"))
        row["poxb_score"] = num(row.get("total_poxB_similarity_score"))
        row["functional_score"] = num(row.get("total_functional_acetoin_score"))
        row["class_color"] = CLASS_INFO.get(row.get("class_final", ""), {}).get("color", "#999999")
        out.append(row)
    return out


def write_fasta(rows: list[dict[str, str]], path: Path) -> None:
    lines = []
    for row in rows:
        header = (
            f">{row['accession']}|gene={row['gene_label']}|class={row.get('class_final','')}"
            f"|ec={row.get('ec_numbers','')}|organism={row.get('organism_name','')}"
        )
        lines.append(header)
        seq = row["sequence"]
        lines.extend(seq[i : i + 80] for i in range(0, len(seq), 80))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    records[current] = "".join(chunks)
                current = line[1:].split("|", 1)[0].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current:
        records[current] = "".join(chunks)
    return records


def pair_identity(a: str, b: str) -> float:
    same = 0
    denom = 0
    for x, y in zip(a, b):
        if x == "-" and y == "-":
            continue
        denom += 1
        if x == y and x != "-":
            same += 1
    return (same / denom * 100.0) if denom else 0.0


def compute_identities(rows: list[dict[str, str]], alignment: Path | None) -> tuple[list[dict], dict[str, dict]]:
    if alignment and alignment.exists():
        seqs = read_fasta(alignment)
        source = "MSA"
    else:
        seqs = {row["accession"]: row["sequence"] for row in rows}
        source = "raw-position"

    pairs: list[dict] = []
    nearest: dict[str, dict] = {}
    for i, row_a in enumerate(rows):
        acc_a = row_a["accession"]
        best = {"partner": "", "identity": -1.0}
        for j, row_b in enumerate(rows):
            if i == j:
                continue
            acc_b = row_b["accession"]
            value = pair_identity(seqs.get(acc_a, ""), seqs.get(acc_b, ""))
            if value > best["identity"]:
                best = {"partner": acc_b, "identity": value}
            if i < j:
                pairs.append(
                    {
                        "a": acc_a,
                        "b": acc_b,
                        "identity": round(value, 2),
                        "class_a": row_a.get("class_final", ""),
                        "class_b": row_b.get("class_final", ""),
                    }
                )
        nearest[acc_a] = best
    for row in rows:
        info = nearest.get(row["accession"], {})
        row["nearest_identity_partner"] = info.get("partner", "")
        row["nearest_identity_percent"] = round(info.get("identity", 0.0), 2)
    return pairs, {"identity_source": source}


def class_summary(rows: list[dict[str, str]]) -> list[dict]:
    order = ["PoxB-like", "POX/SpxB-like", "PDC-like", "Direct_RHEA54368", "AceE/PDH-E1-like"]
    summary = []
    for cls in order:
        members = [r for r in rows if r.get("class_final") == cls]
        if not members:
            continue
        scores = [num(r.get("total_balanced_score")) for r in members]
        info = CLASS_INFO[cls]
        summary.append(
            {
                "class": cls,
                "count": len(members),
                "mean_score": round(mean(scores), 2),
                "top": max(members, key=lambda r: num(r.get("total_balanced_score")))["accession"],
                **info,
            }
        )
    return summary


def records_for_json(rows: list[dict[str, str]]) -> list[dict]:
    keep = [
        "rank",
        "accession",
        "entry_name",
        "reviewed",
        "protein_name",
        "gene_names",
        "gene_label",
        "display_label",
        "short_name",
        "organism_name",
        "taxonomy_id",
        "length",
        "class_final",
        "class_initial",
        "class_color",
        "ec_numbers",
        "rhea_ids",
        "source_query_category",
        "is_seed",
        "cluster70",
        "cluster90",
        "poxB_pident",
        "poxB_qcov",
        "poxB_similarity_group",
        "has_complete_TPP_architecture",
        "has_pox_or_pdc_related_domain",
        "pfam_ids",
        "interpro_ids",
        "best_hmm_family",
        "hmm_support_level",
        "structure_similarity_group",
        "priority_group",
        "risk_flags",
        "rationale",
        "score",
        "poxb_score",
        "functional_score",
        "domain_score",
        "sequence_score",
        "hmm_score",
        "structure_score",
        "practicality_diversity_score",
        "nearest_identity_partner",
        "nearest_identity_percent",
        "sequence",
    ]
    out = []
    for row in rows:
        item = {key: row.get(key, "") for key in keep}
        item["protein_fasta"] = f">{row['accession']} {row.get('protein_name','')} OS={row.get('organism_name','')}\n" + "\n".join(
            row["sequence"][i : i + 80] for i in range(0, len(row["sequence"]), 80)
        )
        out.append(item)
    return out


def write_identity_csv(rows: list[dict[str, str]], pairs: list[dict], path: Path) -> None:
    ids = [row["accession"] for row in rows]
    matrix = {acc: {acc2: "" for acc2 in ids} for acc in ids}
    for acc in ids:
        matrix[acc][acc] = "100.00"
    for pair in pairs:
        matrix[pair["a"]][pair["b"]] = f"{pair['identity']:.2f}"
        matrix[pair["b"]][pair["a"]] = f"{pair['identity']:.2f}"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["accession"] + ids)
        for acc in ids:
            writer.writerow([acc] + [matrix[acc][b] for b in ids])


def build_data(rows: list[dict[str, str]], pairs: list[dict], meta: dict) -> dict:
    values = [p["identity"] for p in pairs]
    duplicate_pairs = [p for p in pairs if p["identity"] >= 90]
    data = {
        "records": records_for_json(rows),
        "identity_pairs": pairs,
        "class_summary": class_summary(rows),
        "pathway_modules": class_summary(rows),
        "meta": {
            "row_count": len(rows),
            "pair_count": len(pairs),
            "identity_source": meta.get("identity_source", ""),
            "mean_identity": round(mean(values), 2) if values else 0,
            "max_identity": round(max(values), 2) if values else 0,
            "min_identity": round(min(values), 2) if values else 0,
            "duplicate_pair_count_90": len(duplicate_pairs),
            "source_csv": "Table2.csv",
        },
        "priority_colors": PRIORITY_COLORS,
        "class_info": CLASS_INFO,
    }
    return data


COMMON_CSS = r"""
:root {
  --bg:#f6f6f2; --panel:#ffffff; --ink:#202326; --muted:#687076; --line:#d8ddd6;
  --soft:#eef1ed; --blue:#2f5f7f; --green:#3e8f65; --amber:#c17c21; --red:#9b4a45;
  --shadow:0 12px 30px rgba(33,37,41,.08);
}
* { box-sizing:border-box; }
body { margin:0; font-family:Arial,"Microsoft YaHei",sans-serif; background:var(--bg); color:var(--ink); letter-spacing:0; }
header { position:sticky; top:0; z-index:5; background:#fff; border-bottom:1px solid var(--line); padding:16px 22px 14px; }
h1 { margin:0 0 10px; font-size:24px; letter-spacing:0; }
h2 { margin:0 0 8px; font-size:22px; letter-spacing:0; }
h3 { margin:14px 0 8px; font-size:16px; letter-spacing:0; }
.controls { display:grid; grid-template-columns:minmax(220px,1.8fr) repeat(4,minmax(130px,1fr)); gap:8px; }
input,select,button { border:1px solid var(--line); border-radius:6px; background:#fff; padding:8px 10px; font-size:14px; }
.stats { display:flex; flex-wrap:wrap; gap:8px; margin-top:10px; color:var(--muted); font-size:13px; }
.stat { background:var(--soft); border-radius:6px; padding:5px 8px; }
.file-links { display:flex; flex-wrap:wrap; gap:8px; margin-top:10px; font-size:13px; }
.file-links a { border:1px solid var(--line); border-radius:6px; background:#fff; padding:6px 9px; text-decoration:none; }
main { display:grid; grid-template-columns:minmax(360px,.92fr) minmax(520px,1.08fr); gap:16px; padding:16px 22px 24px; }
.list { display:flex; flex-direction:column; gap:10px; }
.card { background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:12px; cursor:pointer; box-shadow:var(--shadow); }
.card.active { border-color:var(--blue); box-shadow:0 0 0 2px rgba(47,95,127,.14); }
.row { display:flex; align-items:flex-start; justify-content:space-between; gap:10px; }
.title { font-weight:700; font-size:17px; }
.meta { color:var(--muted); font-size:12px; line-height:1.45; }
.chips { display:flex; flex-wrap:wrap; gap:6px; margin-top:9px; }
.chip { background:var(--soft); border:1px solid transparent; border-radius:999px; padding:3px 8px; font-size:12px; }
.chip.high,.chip.positive_control { color:#fff; background:var(--green); }
.chip.medium { color:#fff; background:var(--amber); }
.detail { background:var(--panel); border:1px solid var(--line); border-radius:8px; min-height:360px; padding:16px; box-shadow:var(--shadow); }
.toolbar { display:flex; flex-wrap:wrap; gap:8px; margin:12px 0; }
.toolbar button.active { background:var(--blue); color:#fff; border-color:var(--blue); }
.kv { display:grid; grid-template-columns:180px 1fr; gap:8px 12px; font-size:14px; margin-top:12px; }
.kv div:nth-child(odd) { color:var(--muted); }
.module-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(220px,1fr)); gap:10px; margin-top:12px; }
.module { border:1px solid var(--line); border-radius:8px; padding:12px; background:#fbfcfb; }
.module b { display:block; margin-bottom:4px; }
.pathway-flow { display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:8px; margin:12px 0; }
.flow-node { border:1px solid var(--line); border-radius:8px; background:#fff; padding:10px; min-height:82px; }
.flow-node small { display:block; color:var(--muted); margin-top:4px; }
pre { white-space:pre-wrap; overflow:auto; background:#f8f8f5; border:1px solid var(--line); border-radius:6px; padding:10px; font:12px/1.45 Consolas,monospace; max-height:420px; }
a { color:var(--blue); }
.identity-summary { display:grid; grid-template-columns:repeat(5,minmax(92px,1fr)); gap:8px; margin:8px 0 10px; }
.identity-summary .stat { text-align:center; }
.identity-heatmap-wrap { overflow:auto; max-height:520px; border:1px solid var(--line); border-radius:8px; background:#fff; margin-top:8px; }
table.identity-heatmap { border-collapse:separate; border-spacing:0; font-size:11px; min-width:max-content; }
.identity-heatmap th,.identity-heatmap td { border-right:1px solid rgba(255,255,255,.42); border-bottom:1px solid rgba(255,255,255,.42); min-width:28px; height:26px; text-align:center; }
.identity-corner { position:sticky; top:0; left:0; z-index:4; min-width:132px; background:#fff; color:var(--muted); border-right:1px solid var(--line); border-bottom:1px solid var(--line); }
.identity-col-head { position:sticky; top:0; z-index:3; height:132px; min-width:28px; vertical-align:bottom; background:#fff; border-bottom:1px solid var(--line); }
.identity-col-head > div { writing-mode:vertical-rl; transform:rotate(180deg); white-space:nowrap; max-height:126px; overflow:hidden; padding:4px 2px; color:var(--muted); font-weight:400; }
.identity-row-head { position:sticky; left:0; z-index:2; min-width:132px; max-width:132px; background:#fff; text-align:right; padding:0 6px; border-right:1px solid var(--line); white-space:nowrap; overflow:hidden; text-overflow:ellipsis; color:var(--muted); font-weight:400; }
.identity-cell { cursor:pointer; color:transparent; }
.identity-cell.show-number { color:#202326; font-size:10px; }
.identity-cell.high-number { color:#fff; }
.identity-cell:hover { outline:2px solid #202326; outline-offset:-2px; }
.identity-cell.diagonal { box-shadow:inset 0 0 0 1px rgba(0,0,0,.22); }
.identity-cell.focus { box-shadow:inset 0 0 0 2px #202326; }
.identity-pair-detail { display:grid; grid-template-columns:136px 1fr; gap:6px 10px; margin:10px 0 12px; font-size:13px; line-height:1.35; border:1px solid var(--line); border-radius:8px; padding:10px; background:#fff; }
.identity-pair-detail div:nth-child(odd) { color:var(--muted); }
.identity-legend { display:grid; grid-template-columns:36px 1fr 42px; gap:8px; align-items:start; font-size:12px; color:var(--muted); margin-top:8px; }
.identity-ramp-wrap { position:relative; padding-bottom:22px; }
.identity-ramp { height:12px; border-radius:999px; border:1px solid var(--line); background:linear-gradient(90deg,#C94C4C 0%,#ECA76A 25%,#F2E7A6 50%,#9BCB9C 75%,#63B6C2 90%,#3F77B5 100%); }
.identity-tick { position:absolute; top:15px; transform:translateX(-50%); white-space:nowrap; font-size:11px; }
.identity-tick::before { content:""; display:block; width:1px; height:6px; margin:0 auto 2px; background:var(--line); }
@media (max-width:980px) {
  main { grid-template-columns:1fr; }
  .controls { grid-template-columns:1fr 1fr; }
  .identity-summary { grid-template-columns:repeat(2,1fr); }
}
"""


COMMON_JS = r"""
const records = DATA.records;
const identityPairs = DATA.identity_pairs;
let filtered = [...records];
let selectedId = records[0]?.accession || "";
let activeTab = "overview";
const $ = id => document.getElementById(id);
function esc(s){return String(s ?? "").replace(/[&<>"']/g, m=>({"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;","'":"&#39;"}[m]));}
function pct(v){const n=Number(v); return Number.isFinite(n)?`${n.toFixed(1)}%`:"";}
function score(v){const n=Number(v); return Number.isFinite(n)?n.toFixed(1):"";}
function rec(id){return records.find(r=>r.accession===id);}
function label(c){return c ? `${c.gene_label || c.short_name} | ${c.accession}` : "";}
function num(v){const n=Number(v); return Number.isFinite(n)?n:null;}
const IDENTITY_STOPS=[[0,[201,76,76]],[25,[236,167,106]],[50,[242,231,166]],[75,[155,203,156]],[90,[99,182,194]],[100,[63,119,181]]];
function mix(a,b,t){return Math.round(a+(b-a)*t);}
function identityColor(v){
  const x=Math.max(0,Math.min(100,Number(v)||0));
  for(let i=1;i<IDENTITY_STOPS.length;i++){
    const [rv,rc]=IDENTITY_STOPS[i], [lv,lc]=IDENTITY_STOPS[i-1];
    if(x<=rv){
      const t=(x-lv)/Math.max(1,rv-lv);
      return `rgb(${mix(lc[0],rc[0],t)}, ${mix(lc[1],rc[1],t)}, ${mix(lc[2],rc[2],t)})`;
    }
  }
  return "rgb(63,119,181)";
}
function pairKey(a,b){return [a,b].sort().join("__");}
const pairMap = new Map();
identityPairs.forEach(p=>pairMap.set(pairKey(p.a,p.b),p.identity));
function identityFor(a,b){if(a===b) return 100; const v=pairMap.get(pairKey(a,b)); return v===undefined?null:Number(v);}
function identityStats(ids){
  const vals=[];
  for(let i=0;i<ids.length;i++) for(let j=i+1;j<ids.length;j++){const v=identityFor(ids[i],ids[j]); if(v!==null) vals.push(v);}
  return {
    pairs: vals.length,
    mean: vals.length?vals.reduce((a,b)=>a+b,0)/vals.length:0,
    max: vals.length?Math.max(...vals):0,
    min: vals.length?Math.min(...vals):0,
    ge90: vals.filter(v=>v>=90).length
  };
}
function identityPairDetailHtml(a,b){
  const ca=rec(a), cb=rec(b), v=identityFor(a,b);
  const value=v===null?"n/a":`${v.toFixed(2)}%`;
  return `<div>Pair</div><div><b>${esc(label(ca))} vs ${esc(label(cb))}</b></div>
    <div>Identity</div><div><span class="chip">${esc(value)}</span></div>
    <div>A class</div><div>${esc(ca?.class_final || "")}</div>
    <div>B class</div><div>${esc(cb?.class_final || "")}</div>
    <div>A organism</div><div>${esc(ca?.organism_name || "")}</div>
    <div>B organism</div><div>${esc(cb?.organism_name || "")}</div>
    <div>A EC/Rhea</div><div>${esc(ca?.ec_numbers || "")} ${esc(ca?.rhea_ids || "")}</div>
    <div>B EC/Rhea</div><div>${esc(cb?.ec_numbers || "")} ${esc(cb?.rhea_ids || "")}</div>`;
}
function renderIdentityPairDetail(a,b){
  const el=$("identityPairDetail");
  if(el) el.innerHTML=identityPairDetailHtml(a,b);
}
function identityHeatmapHtml(anchorId, inputIds){
  const ids=(inputIds && inputIds.length ? inputIds : records.map(r=>r.accession));
  const stats=identityStats(ids);
  const head=`<tr><th class="identity-corner">${ids.length} seq</th>`+ids.map(id=>{
    const c=rec(id);
    return `<th class="identity-col-head" title="${esc(label(c))}"><div>${esc(c?.gene_label || c?.accession || id)}</div></th>`;
  }).join("")+`</tr>`;
  const rows=ids.map(a=>{
    const ca=rec(a);
    const cells=ids.map(b=>{
      const v=identityFor(a,b);
      const text=v===null?"":(v>=70?v.toFixed(0):"");
      const bg=v===null?"#f1f1ec":identityColor(v);
      const cls=["identity-cell",a===b?"diagonal":"",a===anchorId||b===anchorId?"focus":"",text?"show-number":"",v>=94?"high-number":""].filter(Boolean).join(" ");
      return `<td class="${cls}" style="background:${bg}" data-a="${esc(a)}" data-b="${esc(b)}" title="${esc(label(rec(a)))} vs ${esc(label(rec(b)))}: ${v===null?"n/a":v.toFixed(2)+"%"}">${esc(text)}</td>`;
    }).join("");
    return `<tr><th class="identity-row-head" title="${esc(label(ca))}">${esc(ca?.gene_label || a)}</th>${cells}</tr>`;
  }).join("");
  setTimeout(()=>{
    document.querySelectorAll(".identity-cell").forEach(td=>td.addEventListener("click",()=>renderIdentityPairDetail(td.dataset.a,td.dataset.b)));
    if(anchorId) renderIdentityPairDetail(anchorId,anchorId); else if(ids[0]) renderIdentityPairDetail(ids[0],ids[0]);
  },0);
  return `<div class="identity-summary">
      <span class="stat">Sequences ${ids.length}</span>
      <span class="stat">Pairs ${stats.pairs}</span>
      <span class="stat">mean ${stats.mean.toFixed(1)}%</span>
      <span class="stat">max ${stats.max.toFixed(1)}%</span>
      <span class="stat">>=90% ${stats.ge90}</span>
    </div>
    <div id="identityPairDetail" class="identity-pair-detail"></div>
    <div class="identity-heatmap-wrap"><table class="identity-heatmap">${head}${rows}</table></div>
    <div class="identity-legend"><span>low</span><div class="identity-ramp-wrap"><div class="identity-ramp"></div>
      <span class="identity-tick" style="left:0%">0%</span><span class="identity-tick" style="left:25%">25%</span>
      <span class="identity-tick" style="left:50%">50%</span><span class="identity-tick" style="left:75%">75%</span>
      <span class="identity-tick" style="left:90%">90%</span><span class="identity-tick" style="left:100%">100%</span>
      </div><span>high</span></div>`;
}
function unique(field){return [...new Set(records.map(r=>r[field]).filter(Boolean))].sort();}
function fillSelect(id, values){
  const el=$(id); if(!el) return;
  values.forEach(v=>{const opt=document.createElement("option"); opt.value=v; opt.textContent=v; el.appendChild(opt);});
}
function classChip(c){
  const color=c.class_color || "#999";
  return `<span class="chip" style="background:${esc(color)};color:#fff">${esc(c.class_final || "class?")}</span>`;
}
function applyFilters(){
  const q=($("q")?.value||"").toLowerCase();
  const cls=$("classFilter")?.value||"";
  const pri=$("priorityFilter")?.value||"";
  const hmm=$("hmmFilter")?.value||"";
  const minScore=Number($("minScore")?.value||0);
  filtered=records.filter(c=>{
    const blob=[c.accession,c.entry_name,c.protein_name,c.gene_names,c.organism_name,c.class_final,c.ec_numbers,c.rhea_ids,c.rationale,c.risk_flags].join(" ").toLowerCase();
    return (!q||blob.includes(q)) && (!cls||c.class_final===cls) && (!pri||c.priority_group===pri) && (!hmm||c.best_hmm_family===hmm) && (Number(c.score||0)>=minScore);
  });
  if(!filtered.find(c=>c.accession===selectedId)) selectedId=filtered[0]?.accession || records[0]?.accession || "";
}
function renderStats(){
  const total=filtered.length;
  const classes=unique("class_final").map(cls=>`${cls} ${filtered.filter(c=>c.class_final===cls).length}`).join(" | ");
  const high=filtered.filter(c=>c.priority_group==="high"||c.priority_group==="positive_control").length;
  const med=filtered.filter(c=>c.priority_group==="medium").length;
  $("stats").innerHTML=[
    `Showing ${total} / ${records.length}`,
    `high/control ${high}`,
    `medium ${med}`,
    `identity>=90 pairs ${DATA.meta.duplicate_pair_count_90}`,
    `MSA source ${DATA.meta.identity_source}`,
    classes
  ].filter(Boolean).map(x=>`<span class="stat">${esc(x)}</span>`).join("");
}
function renderList(){
  $("list").innerHTML=filtered.map(c=>`<article class="card ${c.accession===selectedId?"active":""}" onclick="selectedId='${esc(c.accession)}';renderDetail();renderList();">
    <div class="row"><div><div class="title">${esc(c.gene_label || c.accession)} <span class="meta">${esc(c.accession)}</span></div>
    <div class="meta">${esc(c.organism_name || "")}</div></div><span class="chip ${esc(c.priority_group)}">${esc(c.priority_group || "")}</span></div>
    <div class="chips">${classChip(c)}<span class="chip">score ${esc(score(c.score))}</span><span class="chip">EC ${esc(c.ec_numbers || "n/a")}</span><span class="chip">nearest ${esc(c.nearest_identity_partner || "")} ${esc(pct(c.nearest_identity_percent))}</span></div>
  </article>`).join("");
}
function current(){return records.find(c=>c.accession===selectedId)||filtered[0]||records[0]||null;}
function tabButton(tab,labelText){return `<button onclick="activeTab='${tab}';renderDetail();" class="${activeTab===tab?"active":""}">${labelText}</button>`;}
function overviewTab(c){
  return `<div class="kv">
    <div>Protein</div><div>${esc(c.protein_name)}</div>
    <div>Organism</div><div>${esc(c.organism_name)}</div>
    <div>Class / EC</div><div>${classChip(c)} ${esc(c.ec_numbers || "")}</div>
    <div>Rhea IDs</div><div>${esc(c.rhea_ids || "")}</div>
    <div>Priority / score</div><div>${esc(c.priority_group)} - balanced ${esc(score(c.score))} - functional ${esc(score(c.functional_score))} - PoxB similarity ${esc(score(c.poxb_score))}</div>
    <div>HMM / domain</div><div>${esc(c.best_hmm_family)} - ${esc(c.hmm_support_level)} - complete TPP ${esc(c.has_complete_TPP_architecture)}</div>
    <div>Clusters</div><div>70% ${esc(c.cluster70)} - 90% ${esc(c.cluster90)}</div>
    <div>Nearest identity</div><div>${esc(c.nearest_identity_partner)} - ${esc(pct(c.nearest_identity_percent))}</div>
    <div>Risk flags</div><div>${esc(c.risk_flags || "none")}</div>
  </div><h3>Rationale</h3><p>${esc(c.rationale || "")}</p>`;
}
function pathwayTab(c){
  const info=DATA.class_info[c.class_final] || {};
  const modules=DATA.pathway_modules.map(m=>`<div class="module" style="border-top:5px solid ${esc(m.color)}">
    <b>${esc(m.stage)} - ${esc(m.class)}</b>
    <div class="meta">${esc(m.module)} - ${esc(m.count)} enzymes - mean score ${esc(m.mean_score)}</div>
    <p>${esc(m.reaction)}</p><p>${esc(m.interpretation)}</p><p class="meta">${esc(m.caution)}</p>
  </div>`).join("");
  return `<h3>Selected pathway logic</h3>
    <div class="module" style="border-top:5px solid ${esc(info.color||"#999")}"><b>${esc(c.class_final)}</b><p>${esc(info.reaction||"")}</p><p>${esc(info.interpretation||"")}</p><p class="meta">${esc(info.caution||"")}</p></div>
    <h3>KEGG-like pathway abstraction</h3>
    <div class="pathway-flow">
      <div class="flow-node"><b>fluoro-/nonfluoro 2-oxo acid pool</b><small>pyruvate, fluoropyruvate-like substrates, branched 2-oxo acids</small></div>
      <div class="flow-node"><b>oxidative decarboxylation</b><small>PoxB-like and POX/SpxB-like branches</small></div>
      <div class="flow-node"><b>aldehyde branch</b><small>PDC/direct Rhea 2ODC branch</small></div>
      <div class="flow-node"><b>C-C / activated acyl branch</b><small>ALS/PDH-E1-like branch, acetolactate/acetyl phosphate logic</small></div>
    </div>
    <h3>Modules in this CSV</h3><div class="module-grid">${modules}</div>`;
}
function fastaTab(c){return `<h3>Protein FASTA</h3><pre>${esc(c.protein_fasta || "")}</pre>`;}
function identityTab(c){return `<h3>Identity matrix</h3>${identityHeatmapHtml(c.accession, records.map(r=>r.accession))}`;}
function notesTab(c){
  return `<h3>How to analyze as a multi-gene/pathway problem</h3>
  <div class="module-grid">
    <div class="module"><b>1. EC/Rhea pathway graph</b><p>Use EC and Rhea IDs as reaction nodes, then place each enzyme class on pyruvate/2-oxo-acid, aldehyde, acetyl-phosphate, and C-C forming branches.</p></div>
    <div class="module"><b>2. SSN + GNN logic</b><p>Keep the sequence similarity network separate from genome-neighborhood context. Similar enzymes can sit in different metabolic neighborhoods.</p></div>
    <div class="module"><b>3. Gene-neighborhood window</b><p>For bacterial rows, inspect +/-10 genes for transporters, dehydrogenases, thioesterases, halogenation/fluorination genes, TPP enzyme partners, and regulators.</p></div>
    <div class="module"><b>4. Synteny and redundancy</b><p>Use cluster70/cluster90 and identity>=90 pairs to avoid ordering multiple near-duplicates unless they represent different pathway neighborhoods.</p></div>
  </div>
  <h3>Current CSV interpretation</h3>
  <p>The panel is a TPP-enzyme-centered pathway panel. It does not need more homolog discovery now; the useful question is which reaction module each enzyme supports and whether the module connects to a plausible fluorinated 2-oxo-acid route.</p>`;
}
function renderDetail(){
  const c=current();
  if(!c){$("detail").innerHTML="<p>No records.</p>"; return;}
  const body=activeTab==="pathway"?pathwayTab(c):activeTab==="fasta"?fastaTab(c):activeTab==="identity"?identityTab(c):activeTab==="notes"?notesTab(c):overviewTab(c);
  $("detail").innerHTML=`<h2>${esc(label(c))}</h2>
    <div class="meta">${esc(c.class_final)} - ${esc(c.organism_name || "")} - score ${esc(score(c.score))}</div>
    <div class="toolbar">${tabButton("overview","Overview")}${tabButton("pathway","Pathway")}${tabButton("fasta","FASTA")}${tabButton("identity","Identity")}${tabButton("notes","Analysis notes")}</div>
    ${body}`;
}
function render(){applyFilters(); renderStats(); renderList(); renderDetail();}
"""


def html_shell(title: str, body: str, data: dict, extra_js: str = "") -> str:
    return (
        "<!doctype html>\n<html lang=\"zh-CN\"><head><meta charset=\"utf-8\">"
        f"<title>{html.escape(title)}</title><style>{COMMON_CSS}</style></head><body>"
        f"{body}<script>const DATA = {json.dumps(data, ensure_ascii=False)};\n{COMMON_JS}\n{extra_js}</script>"
        "</body></html>\n"
    )


def write_dashboard(data: dict, path: Path) -> None:
    body = r"""
<header>
  <h1>Fluorinated Pathway Enzyme Dashboard</h1>
  <div class="controls">
    <input id="q" placeholder="Search accession, gene, organism, EC, Rhea, rationale">
    <select id="classFilter"><option value="">All classes</option></select>
    <select id="priorityFilter"><option value="">All priorities</option></select>
    <select id="hmmFilter"><option value="">All HMM families</option></select>
    <select id="minScore"><option value="0">score >= 0</option><option value="70">score >= 70</option><option value="80">score >= 80</option><option value="85">score >= 85</option></select>
  </div>
  <div class="stats" id="stats"></div>
  <div class="file-links">
    <a href="identity_matrix.html">Standalone identity matrix</a>
    <a href="table2_enzymes.fasta">Protein FASTA</a>
    <a href="identity_matrix.csv">Identity CSV</a>
    <a href="pathway_relationship_summary.md">Pathway summary</a>
  </div>
</header>
<main>
  <section class="list" id="list"></section>
  <section class="detail" id="detail"></section>
</main>
"""
    extra_js = r"""
fillSelect("classFilter", unique("class_final"));
fillSelect("priorityFilter", unique("priority_group"));
fillSelect("hmmFilter", unique("best_hmm_family"));
["q","classFilter","priorityFilter","hmmFilter","minScore"].forEach(id=>$(id)?.addEventListener("input",render));
["classFilter","priorityFilter","hmmFilter","minScore"].forEach(id=>$(id)?.addEventListener("change",render));
render();
"""
    path.write_text(html_shell("Fluorinated Pathway Enzyme Dashboard", body, data, extra_js), encoding="utf-8")


def write_identity_html(data: dict, path: Path) -> None:
    body = r"""
<header>
  <h1>Table2 Identity Matrix</h1>
  <div class="controls">
    <input id="q" placeholder="Search accession, gene, organism, EC, Rhea">
    <select id="classFilter"><option value="">All classes</option></select>
    <select id="priorityFilter"><option value="">All priorities</option></select>
    <select id="hmmFilter"><option value="">All HMM families</option></select>
    <select id="minScore"><option value="0">score >= 0</option><option value="70">score >= 70</option><option value="80">score >= 80</option><option value="85">score >= 85</option></select>
  </div>
  <div class="stats" id="stats"></div>
  <div class="file-links"><a href="fluorinated_pathway_enzyme_dashboard.html">Back to dashboard</a><a href="identity_matrix.csv">Identity CSV</a><a href="pathway_relationship_summary.md">Pathway summary</a></div>
</header>
<main style="display:block">
  <section class="detail" id="matrix"></section>
</main>
"""
    extra_js = r"""
fillSelect("classFilter", unique("class_final"));
fillSelect("priorityFilter", unique("priority_group"));
fillSelect("hmmFilter", unique("best_hmm_family"));
function renderMatrixOnly(){
  applyFilters(); renderStats();
  $("matrix").innerHTML=identityHeatmapHtml(filtered[0]?.accession || records[0]?.accession, filtered.map(r=>r.accession));
}
["q","classFilter","priorityFilter","hmmFilter","minScore"].forEach(id=>$(id)?.addEventListener("input",renderMatrixOnly));
["classFilter","priorityFilter","hmmFilter","minScore"].forEach(id=>$(id)?.addEventListener("change",renderMatrixOnly));
renderMatrixOnly();
"""
    path.write_text(html_shell("Table2 Identity Matrix", body, data, extra_js), encoding="utf-8")


def write_summary(data: dict, path: Path) -> None:
    lines = [
        "# Table2 fluorinated-pathway enzyme relationship summary",
        "",
        f"- Records: {data['meta']['row_count']}",
        f"- Pairwise identities: {data['meta']['pair_count']} pairs from {data['meta']['identity_source']}",
        f"- Mean identity: {data['meta']['mean_identity']}%",
        f"- Max identity: {data['meta']['max_identity']}%",
        f"- Identity >=90% pairs: {data['meta']['duplicate_pair_count_90']}",
        "",
        "## Module view",
        "",
    ]
    for module in data["pathway_modules"]:
        lines.extend(
            [
                f"### {module['stage']} - {module['class']} ({module['count']} enzymes)",
                f"- Reaction logic: {module['reaction']}",
                f"- Interpretation: {module['interpretation']}",
                f"- Caution: {module['caution']}",
                "",
            ]
        )
    lines.extend(
        [
            "## Practical interpretation",
            "",
            "- Treat this CSV as a TPP-enzyme-centered pathway panel, not as a set requiring more homolog discovery.",
            "- Use EC/Rhea to place enzymes on a KEGG-like reaction graph, then use identity/cluster70/cluster90 to avoid redundant ordering.",
            "- For bacterial candidates, a next useful layer is genome-neighborhood analysis around each accession: inspect nearby transporters, dehydrogenases, thioesterases, regulators, fluorination/halogenation genes, and TPP-partner enzymes.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--alignment")
    args = parser.parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = read_rows(input_path)
    rows.sort(key=lambda r: (["PoxB-like", "POX/SpxB-like", "PDC-like", "Direct_RHEA54368", "AceE/PDH-E1-like"].index(r["class_final"]) if r["class_final"] in CLASS_INFO else 99, -num(r.get("total_balanced_score")), r["accession"]))

    fasta_path = outdir / "table2_enzymes.fasta"
    write_fasta(rows, fasta_path)
    alignment_path = Path(args.alignment) if args.alignment else None
    pairs, meta = compute_identities(rows, alignment_path)
    data = build_data(rows, pairs, meta)

    write_identity_csv(rows, pairs, outdir / "identity_matrix.csv")
    write_dashboard(data, outdir / "fluorinated_pathway_enzyme_dashboard.html")
    write_identity_html(data, outdir / "identity_matrix.html")
    write_summary(data, outdir / "pathway_relationship_summary.md")
    (outdir / "dashboard_data.json").write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"records={len(rows)}")
    print(f"fasta={fasta_path}")
    print(f"dashboard={outdir / 'fluorinated_pathway_enzyme_dashboard.html'}")
    print(f"identity={outdir / 'identity_matrix.html'}")


if __name__ == "__main__":
    main()
