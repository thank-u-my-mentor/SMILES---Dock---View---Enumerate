#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


sys.path.insert(0, str(Path(__file__).resolve().parent))
from mf_plasmid_dashboard_tools import read_xls  # noqa: E402
from lipase_order_dashboard_builder import clean_dna, needleman_identity, translate_dna, wrap_fasta  # noqa: E402


def clean_text(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "").replace("\xa0", " ")).strip()


def norm_acc(value: str) -> str:
    return clean_text(value).split()[0].strip().strip(".")


def row_value(row: list[Any], idx: int) -> str:
    return clean_text(row[idx]) if idx < len(row) else ""


def extract_data_block(html: str) -> tuple[dict[str, Any], int, int]:
    marker = "const DATA = "
    start = html.find(marker)
    if start < 0:
        raise ValueError("Could not find const DATA block")
    json_start = start + len(marker)
    end = html.find(";\nconst candidates", json_start)
    if end < 0:
        end = html.find(";\r\nconst candidates", json_start)
    if end < 0:
        raise ValueError("Could not find end of const DATA block")
    return json.loads(html[json_start:end]), json_start, end


def read_order_range(path: Path, start_gene: str, end_gene: str) -> list[dict[str, Any]]:
    sheets = read_xls(path)
    rows = next(iter(sheets.values()))
    out: list[dict[str, Any]] = []
    in_range = False
    for excel_row, row in enumerate(rows, start=1):
        gene = row_value(row, 1)
        if gene == start_gene:
            in_range = True
        if not in_range:
            continue
        accession = norm_acc(row_value(row, 22))
        dna = clean_dna(row_value(row, 5))
        if gene and accession and dna:
            protein, terminal_stop, internal_stop_count, trailing = translate_dna(dna)
            out.append(
                {
                    "excel_row": excel_row,
                    "serial": row_value(row, 0),
                    "gene_name": gene,
                    "accession": accession,
                    "five_prime_site": row_value(row, 2),
                    "three_prime_site": row_value(row, 8),
                    "optimization_host": row_value(row, 9),
                    "standard_vector": row_value(row, 12),
                    "vector_name": row_value(row, 13),
                    "antibiotic": row_value(row, 14),
                    "plasmid_prep": row_value(row, 19),
                    "endotoxin": row_value(row, 20),
                    "order_note": row_value(row, 24),
                    "dna_sequence": dna,
                    "protein_sequence_from_order": protein,
                    "dna_length": len(dna),
                    "protein_length_from_order": len(protein),
                    "terminal_stop": terminal_stop,
                    "internal_stop_count": internal_stop_count,
                    "incomplete_codon_bases": trailing,
                }
            )
        if gene == end_gene:
            break
    return out


def merge_candidates(data: dict[str, Any], order_rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[str], list[str]]:
    order_by_acc = {row["accession"]: row for row in order_rows}
    candidates_by_acc = {norm_acc(row.get("accession", "")): row for row in data.get("candidates", [])}
    kept: list[dict[str, Any]] = []
    missing_from_dashboard: list[str] = []
    for order in order_rows:
        acc = order["accession"]
        candidate = dict(candidates_by_acc.get(acc) or {})
        if not candidate:
            missing_from_dashboard.append(acc)
            candidate = {
                "accession": acc,
                "rank": order.get("serial", ""),
                "protein_name": "4-hydroxy-2-oxovalerate aldolase family protein",
                "organism": "",
                "target_ec": "4.1.3.39",
                "all_ec_numbers": "4.1.3.39",
                "source_bucket": "Unknown_source",
                "kingdom": "Bacteria",
                "phylum_class": "",
                "reviewed": "",
                "length": str(order.get("protein_length_from_order", "")),
                "sequence_status": "translated_from_order_sheet",
                "protein_sequence": order.get("protein_sequence_from_order", ""),
                "fasta_header": f"{acc} {order.get('gene_name', '')}",
                "uniprot_url": f"https://www.uniprot.org/uniprotkb/{acc}/entry",
                "construct_decision": "REVIEW",
                "construct_score": "",
                "construct_strategy": "codon_optimize_for_E_coli_and_order_synthetic_gene",
                "construct_reasons": "from plasmid order sheet",
                "construct_cautions": "",
            }
        candidate["original_rank"] = clean_text(candidate.get("rank", ""))
        candidate["rank"] = clean_text(order.get("serial", "")) or candidate["original_rank"]
        candidate["gene_name"] = order["gene_name"]
        candidate["gene_label"] = order["gene_name"]
        candidate["display_name"] = order["gene_name"]
        candidate["order_sheet"] = order
        candidate["order_serial"] = order.get("serial", "")
        candidate["order_excel_row"] = str(order.get("excel_row", ""))
        candidate["order_dna_length"] = str(order.get("dna_length", ""))
        candidate["order_protein_length"] = str(order.get("protein_length_from_order", ""))
        candidate["aldolase_class"] = "Type II aldolase"
        candidate["primary_activity"] = "HOA-family activity: 4-hydroxy-2-oxovalerate aldolase"
        candidate["metal_dependency"] = "Metal-dependent class II aldolase assignment; HOA-family activity is treated here as divalent-metal dependent."
        candidate["classification_note"] = "Classified as Type II aldolase for this library because the target HOA activity is metal-dependent."
        candidate["order_protein_fasta"] = (
            f">{order['gene_name']}|{acc}|translated_from_order_sheet|row={order['excel_row']}\n"
            f"{wrap_fasta(order.get('protein_sequence_from_order', ''))}"
        )
        candidate["order_dna_fasta"] = (
            f">{order['gene_name']}|{acc}|coding_sequence_from_order_sheet|row={order['excel_row']}\n"
            f"{wrap_fasta(order.get('dna_sequence', ''))}"
        )
        warnings: list[str] = []
        old_seq = re.sub(r"[^A-Za-z]", "", str(candidate.get("protein_sequence", ""))).upper()
        new_seq = order.get("protein_sequence_from_order", "")
        if old_seq and new_seq and old_seq != new_seq:
            if len(old_seq) != len(new_seq):
                warnings.append(f"Dashboard protein length {len(old_seq)}; order-sheet translation length {len(new_seq)}")
            else:
                mismatches = sum(1 for a, b in zip(old_seq, new_seq) if a != b)
                if mismatches:
                    warnings.append(f"Order-sheet translated protein differs at {mismatches} residue(s)")
        if order.get("internal_stop_count"):
            warnings.append(f"{order['internal_stop_count']} internal stop codon(s) in order sequence")
        if order.get("incomplete_codon_bases"):
            warnings.append(f"{order['incomplete_codon_bases']} trailing base(s) ignored in translation")
        candidate["order_warnings"] = warnings
        kept.append(candidate)
    used = {row["accession"] for row in order_rows}
    removed = [acc for acc in candidates_by_acc if acc and acc not in used]
    return kept, missing_from_dashboard, removed


def build_identity(candidates: list[dict[str, Any]]) -> tuple[list[str], list[list[float]], list[dict[str, Any]]]:
    ids = [row["accession"] for row in candidates]
    matrix = [[100.0 if i == j else 0.0 for j in range(len(ids))] for i in range(len(ids))]
    pairs: list[dict[str, Any]] = []
    for i, a in enumerate(candidates):
        for j in range(i + 1, len(candidates)):
            b = candidates[j]
            frac, identical, compared = needleman_identity(str(a.get("protein_sequence", "")), str(b.get("protein_sequence", "")))
            pct = round(frac * 100, 4)
            matrix[i][j] = matrix[j][i] = pct
            pairs.append(
                {
                    "query": a["accession"],
                    "target": b["accession"],
                    "query_gene": a.get("gene_name", ""),
                    "target_gene": b.get("gene_name", ""),
                    "identity_percent": pct,
                    "identity_fraction": round(frac, 6),
                    "identical_sites": identical,
                    "compared_sites": compared,
                    "query_header": a.get("fasta_header", ""),
                    "target_header": b.get("fasta_header", ""),
                }
            )
    pairs.sort(key=lambda row: (-float(row["identity_percent"]), row["query_gene"], row["target_gene"]))
    return ids, matrix, pairs


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(outdir: Path, source_html: Path, source_xls: Path, data: dict[str, Any]) -> dict[str, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    html_path = outdir / "HOA_2型醛缩酶.html"
    json_path = outdir / "hoa_type_ii_aldolase_dashboard_data.json"
    protein_fasta = outdir / "hoa_type_ii_aldolase_order_proteins.fasta"
    dna_fasta = outdir / "hoa_type_ii_aldolase_order_coding_sequences.fasta"
    matrix_csv = outdir / "hoa_type_ii_aldolase_identity_matrix.csv"
    records_csv = outdir / "hoa_type_ii_aldolase_order_records.csv"

    html_path.write_text(HTML_TEMPLATE.replace("__DATA_JSON__", json.dumps(data, ensure_ascii=False).replace("</", "<\\/")), encoding="utf-8")
    json_path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    protein_fasta.write_text(
        "\n\n".join(row.get("order_protein_fasta", "") for row in data["candidates"]) + "\n",
        encoding="utf-8",
    )
    dna_fasta.write_text(
        "\n\n".join(row.get("order_dna_fasta", "") for row in data["candidates"]) + "\n",
        encoding="utf-8",
    )
    with matrix_csv.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        labels = [row.get("gene_name", row.get("accession", "")) for row in data["candidates"]]
        writer.writerow(["gene_name", *labels])
        for label, row in zip(labels, data["identity_matrix"]):
            writer.writerow([label, *row])
    record_rows = []
    for row in data["candidates"]:
        order = row.get("order_sheet", {})
        record_rows.append(
            {
                "serial": row.get("order_serial", ""),
                "excel_row": row.get("order_excel_row", ""),
                "gene_name": row.get("gene_name", ""),
                "accession": row.get("accession", ""),
                "aldolase_class": row.get("aldolase_class", ""),
                "primary_activity": row.get("primary_activity", ""),
                "target_ec": row.get("target_ec", ""),
                "organism": row.get("organism", ""),
                "protein_name": row.get("protein_name", ""),
                "pdb_ids": row.get("pdb_ids", ""),
                "vector": order.get("standard_vector", ""),
                "antibiotic": order.get("antibiotic", ""),
                "dna_length": order.get("dna_length", ""),
                "translated_aa_length": order.get("protein_length_from_order", ""),
                "warnings": "; ".join(row.get("order_warnings", [])),
            }
        )
    write_csv(records_csv, record_rows)
    return {
        "html": html_path,
        "json": json_path,
        "protein_fasta": protein_fasta,
        "dna_fasta": dna_fasta,
        "matrix_csv": matrix_csv,
        "records_csv": records_csv,
    }


def build(args: argparse.Namespace) -> dict[str, Any]:
    source_html = args.html
    html = source_html.read_text(encoding="utf-8")
    old_data, _, _ = extract_data_block(html)
    order_rows = read_order_range(args.xls, args.start_gene, args.end_gene)
    candidates, missing_from_dashboard, removed_from_order_range = merge_candidates(old_data, order_rows)
    identity_ids, identity_matrix, identity_pairs = build_identity(candidates)
    data = {
        "title": args.title,
        "source_html": str(source_html),
        "source_xls": str(args.xls),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "range": {
            "start_gene": args.start_gene,
            "end_gene": args.end_gene,
            "order_rows_read": len(order_rows),
        },
        "summary": {
            "candidate_count": len(candidates),
            "identity_pair_count": len(identity_pairs),
            "with_pdb": sum(1 for row in candidates if clean_text(row.get("pdb_ids", ""))),
            "with_order_warnings": sum(1 for row in candidates if row.get("order_warnings")),
            "missing_from_dashboard": missing_from_dashboard,
            "removed_not_in_order_range": removed_from_order_range,
        },
        "candidates": candidates,
        "identity_ids": identity_ids,
        "identity_matrix": identity_matrix,
        "pairs": identity_pairs,
        "pdbFiles": old_data.get("pdbFiles", []),
    }
    outputs = write_outputs(args.outdir, source_html, args.xls, data)
    return {
        "ok": True,
        "candidate_count": len(candidates),
        "identity_pairs": len(identity_pairs),
        "missing_from_dashboard": missing_from_dashboard,
        "removed_not_in_order_range": removed_from_order_range,
        "outputs": {key: str(value) for key, value in outputs.items()},
    }


def check(args: argparse.Namespace) -> dict[str, Any]:
    html_path = args.outdir / "HOA_2型醛缩酶.html"
    json_path = args.outdir / "hoa_type_ii_aldolase_dashboard_data.json"
    html = html_path.read_text(encoding="utf-8")
    data = json.loads(json_path.read_text(encoding="utf-8"))
    problems: list[str] = []
    candidates = data.get("candidates", [])
    if len(candidates) != args.expected_candidates:
        problems.append(f"Expected {args.expected_candidates} candidates, found {len(candidates)}")
    if any(not row.get("gene_name") for row in candidates):
        problems.append("At least one candidate is missing gene_name")
    if any(row.get("aldolase_class") != "Type II aldolase" for row in candidates):
        problems.append("At least one candidate is not marked Type II aldolase")
    matrix = data.get("identity_matrix", [])
    if len(matrix) != len(candidates) or any(len(row) != len(candidates) for row in matrix):
        problems.append("Identity matrix dimensions do not match candidate count")
    for marker in ["identity-heatmap", "Type II aldolase", "Metal-dependent", "Order Sheet", "FASTA"]:
        if marker not in html:
            problems.append(f"HTML missing marker: {marker}")
    return {"ok": not problems, "problems": problems, "candidates": len(candidates), "html": str(html_path)}


HTML_TEMPLATE = r"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>HOA Type II Aldolase Candidate Dashboard</title>
<style>
:root {
  --bg:#f6f7f2; --panel:#fff; --ink:#202326; --muted:#65706d; --line:#d8ded6;
  --blue:#315f9c; --teal:#257b80; --green:#4f8d60; --amber:#a8742c; --red:#b65348; --soft:#eef1ed;
  --shadow:0 1px 2px rgba(20,30,24,.06),0 8px 24px rgba(20,30,24,.07);
}
* { box-sizing:border-box; }
body { margin:0; font-family:Arial,"Microsoft YaHei",sans-serif; background:var(--bg); color:var(--ink); letter-spacing:0; }
header { position:sticky; top:0; z-index:10; background:#fff; border-bottom:1px solid var(--line); padding:14px 18px 12px; }
h1 { margin:0 0 8px; font-size:23px; line-height:1.18; }
.sub { display:flex; gap:9px; flex-wrap:wrap; color:var(--muted); font-size:13px; }
.controls { margin-top:12px; display:grid; grid-template-columns:minmax(240px,1.5fr) repeat(4,minmax(120px,1fr)); gap:8px; }
input,select,button { border:1px solid var(--line); border-radius:6px; background:#fff; padding:8px 10px; font:14px Arial,"Microsoft YaHei",sans-serif; color:var(--ink); }
button { cursor:pointer; }
button.active { background:var(--blue); color:#fff; border-color:var(--blue); }
main { display:grid; grid-template-columns:minmax(320px,430px) minmax(0,1fr); gap:14px; padding:14px 18px 22px; align-items:start; }
.list { display:flex; flex-direction:column; gap:8px; max-height:calc(100vh - 150px); overflow:auto; padding-right:4px; }
.card { background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:10px; cursor:pointer; }
.card.active { border-color:var(--teal); box-shadow:0 0 0 2px rgba(37,123,128,.15); }
.row { display:flex; justify-content:space-between; gap:8px; align-items:flex-start; }
.gene { font-weight:700; font-size:17px; overflow-wrap:anywhere; }
.acc { display:block; color:var(--muted); font-weight:400; font-size:12px; margin-top:2px; }
.meta { color:var(--muted); font-size:12px; line-height:1.42; margin-top:5px; overflow-wrap:anywhere; }
.chips { display:flex; flex-wrap:wrap; gap:5px; margin-top:8px; }
.chip { display:inline-flex; align-items:center; min-height:20px; border-radius:999px; padding:2px 7px; background:var(--soft); color:#43504b; font-size:12px; white-space:nowrap; }
.BUY,.green { background:#dceee1; color:#245d35; }
.REVIEW,.amber { background:#f5ead1; color:#805819; }
.LOW_PRIORITY,.red { background:#f5dfda; color:#8b382e; }
.blue { background:#e1eafa; color:#244f88; }
.detail { min-width:0; background:var(--panel); border:1px solid var(--line); border-radius:8px; box-shadow:var(--shadow); }
.detail-head { padding:14px 16px 12px; border-bottom:1px solid var(--line); }
.detail-head h2 { margin:0; font-size:22px; line-height:1.18; overflow-wrap:anywhere; }
.tabs { display:flex; gap:6px; flex-wrap:wrap; padding:10px 12px; border-bottom:1px solid var(--line); background:#fbfcf8; }
.content { padding:14px 16px 18px; }
.kv { display:grid; grid-template-columns:minmax(140px,190px) minmax(0,1fr); gap:8px 12px; font-size:14px; line-height:1.4; }
.kv div:nth-child(odd) { color:var(--muted); }
a { color:var(--blue); text-decoration:none; }
a:hover { text-decoration:underline; }
pre { margin:0; white-space:pre-wrap; overflow:auto; overflow-wrap:anywhere; max-height:520px; border:1px solid var(--line); border-radius:8px; background:#f6f7f2; padding:11px; font:12px/1.45 Consolas,"Courier New",monospace; }
textarea { width:100%; min-height:130px; resize:vertical; border:1px solid var(--line); border-radius:8px; padding:10px; font:13px/1.45 Arial,"Microsoft YaHei",sans-serif; }
.seq-grid { display:grid; gap:12px; }
.identity-wrap { overflow:auto; max-height:600px; border:1px solid var(--line); border-radius:8px; background:#fff; }
table.identity-heatmap { border-collapse:separate; border-spacing:0; font-size:11px; min-width:max-content; }
.identity-heatmap th,.identity-heatmap td { border-right:1px solid rgba(255,255,255,.5); border-bottom:1px solid rgba(255,255,255,.5); width:34px; min-width:34px; height:28px; text-align:center; }
.corner { position:sticky; top:0; left:0; z-index:4; min-width:132px; background:#fff; color:var(--muted); border-right:1px solid var(--line); border-bottom:1px solid var(--line); }
.col-head { position:sticky; top:0; z-index:3; height:124px; vertical-align:bottom; background:#fff; border-bottom:1px solid var(--line); }
.col-head > div { writing-mode:vertical-rl; transform:rotate(180deg); white-space:nowrap; max-height:118px; overflow:hidden; padding:4px 2px; color:var(--muted); font-weight:400; }
.row-head { position:sticky; left:0; z-index:2; min-width:132px; max-width:132px; background:#fff; text-align:right; padding:0 7px; border-right:1px solid var(--line); color:var(--muted); font-weight:400; white-space:nowrap; overflow:hidden; text-overflow:ellipsis; }
.cell { cursor:pointer; color:transparent; }
.cell.show { color:#202326; font-size:10px; }
.cell.dark { color:#fff; }
.cell.focus { box-shadow:inset 0 0 0 2px #202326; }
.cell:hover { outline:2px solid #202326; outline-offset:-2px; }
.pair-detail { display:grid; grid-template-columns:130px 1fr; gap:6px 10px; border:1px solid var(--line); border-radius:8px; padding:10px; margin:0 0 10px; font-size:13px; }
.pair-detail div:nth-child(odd) { color:var(--muted); }
.legend { display:grid; grid-template-columns:36px 1fr 42px; gap:8px; align-items:start; font-size:12px; color:var(--muted); margin-top:9px; }
.ramp-wrap { position:relative; padding-bottom:22px; }
.ramp { height:12px; border-radius:999px; border:1px solid var(--line); background:linear-gradient(90deg,#b65348 0%,#d9a048 35%,#e8dc91 55%,#81b57f 78%,#257b80 100%); }
.tick { position:absolute; top:15px; transform:translateX(-50%); white-space:nowrap; font-size:11px; }
.tick::before { content:""; display:block; width:1px; height:6px; margin:0 auto 2px; background:var(--line); }
.warning { border-left:3px solid var(--red); background:#fbf1ef; padding:7px 9px; border-radius:5px; font-size:13px; margin-top:7px; }
@media (max-width: 980px) {
  main { grid-template-columns:1fr; }
  .list { max-height:none; }
  .controls { grid-template-columns:1fr 1fr; }
  .kv { grid-template-columns:1fr; }
}
</style>
</head>
<body>
<header>
  <h1 id="title"></h1>
  <div class="sub" id="subtitle"></div>
  <div class="controls">
    <input id="q" placeholder="Search gene, accession, organism, EC, PDB, DOI, expression">
    <select id="decision"></select>
    <select id="source"></select>
    <select id="structure"></select>
    <select id="sort">
      <option value="order">Order sheet</option>
      <option value="gene">Gene name</option>
      <option value="identity">Nearest identity</option>
      <option value="length">Length</option>
    </select>
  </div>
</header>
<main>
  <section class="list" id="list"></section>
  <section class="detail" id="detail"></section>
</main>
<script>
const DATA = __DATA_JSON__;
const candidates = DATA.candidates || [];
const matrix = DATA.identity_matrix || [];
const pairs = DATA.pairs || [];
let selected = candidates[0]?.accession || "";
let tab = "overview";
const $ = id => document.getElementById(id);
function esc(value){return String(value ?? "").replace(/[&<>"']/g,ch=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));}
function uniq(field){return [...new Set(candidates.map(c=>c[field]).filter(Boolean))].sort();}
function fillSelect(id,label,values){$(id).innerHTML=`<option value="">${esc(label)}</option>`+values.map(v=>`<option value="${esc(v)}">${esc(v)}</option>`).join("");}
function splitSemi(v){return String(v||"").split(/[;,]/).map(x=>x.trim()).filter(Boolean).filter((x,i,a)=>a.indexOf(x)===i);}
function accUrl(acc){return `https://www.uniprot.org/uniprotkb/${encodeURIComponent(String(acc||"").replace(/\..*$/,""))}/entry`;}
function pdbLinks(v){const xs=splitSemi(v); return xs.length?xs.map(x=>`<a target="_blank" href="https://www.rcsb.org/structure/${esc(x)}">${esc(x)}</a>`).join("; "):"";}
function doiLinks(v){const xs=splitSemi(v); return xs.length?xs.map(x=>`<a target="_blank" href="https://doi.org/${esc(x)}">${esc(x)}</a>`).join("; "):"";}
function gene(c){return c?.gene_name || c?.uniprot_gene_names || c?.accession || "";}
function byAcc(acc){return candidates.find(c=>c.accession===acc) || candidates[0];}
function idx(acc){return candidates.findIndex(c=>c.accession===acc);}
function passes(c){
  const q=$("q").value.toLowerCase().trim();
  const d=$("decision").value, s=$("source").value, st=$("structure").value;
  const blob=[gene(c),c.accession,c.organism,c.target_ec,c.protein_name,c.pdb_ids,c.uniprot_dois,c.pdb_expression_hosts,c.construct_reasons,c.construct_cautions,c.aldolase_class,c.primary_activity,c.metal_dependency].join(" ").toLowerCase();
  return (!q||blob.includes(q))&&(!d||c.construct_decision===d)&&(!s||c.source_bucket===s)&&(!st||c.structure_status===st);
}
function sortedRows(){
  const rows=candidates.filter(passes);
  const mode=$("sort").value;
  if(mode==="gene") rows.sort((a,b)=>gene(a).localeCompare(gene(b)));
  if(mode==="identity") rows.sort((a,b)=>Number(b.nearest_panel_identity_percent||0)-Number(a.nearest_panel_identity_percent||0));
  if(mode==="length") rows.sort((a,b)=>Number(b.length||0)-Number(a.length||0));
  if(mode==="order") rows.sort((a,b)=>Number(a.order_serial||999)-Number(b.order_serial||999));
  return rows;
}
function renderList(){
  const rows=sortedRows();
  if(!rows.some(c=>c.accession===selected) && rows[0]) selected=rows[0].accession;
  $("list").innerHTML=rows.map(c=>`<article class="card ${c.accession===selected?"active":""}" data-acc="${esc(c.accession)}">
    <div class="row"><div class="gene">${esc(gene(c))}<span class="acc">${esc(c.accession)} · row ${esc(c.order_excel_row)} · serial ${esc(c.order_serial)}</span></div><span class="chip ${esc(c.construct_decision)}">${esc(c.construct_decision||"")}</span></div>
    <div class="meta">${esc(c.organism || "organism not resolved")} | ${esc(c.protein_name || "")}</div>
    <div class="chips"><span class="chip blue">Type II aldolase</span><span class="chip green">metal-dependent</span><span class="chip">${esc(c.target_ec||"EC ?")}</span><span class="chip">${esc(c.structure_status==="rcsb_structure_found"?"PDB "+c.pdb_ids:"no PDB")}</span></div>
  </article>`).join("") || "<p class='meta'>No matching candidates.</p>";
  document.querySelectorAll(".card").forEach(el=>el.addEventListener("click",()=>{selected=el.dataset.acc;render();}));
}
function field(k,v){return `<div>${esc(k)}</div><div>${v || ""}</div>`;}
function overview(c){
  const warn=(c.order_warnings||[]).map(w=>`<div class="warning">${esc(w)}</div>`).join("");
  return `<div class="kv">
    ${field("Gene name", `<b>${esc(gene(c))}</b>`)}
    ${field("UniProt", `<a target="_blank" href="${esc(c.uniprot_url || accUrl(c.accession))}">${esc(c.accession)}</a>`)}
    ${field("Aldolase class", `<span class="chip blue">${esc(c.aldolase_class)}</span> <span class="chip green">Metal-dependent</span>`)}
    ${field("Primary activity", esc(c.primary_activity || ""))}
    ${field("Classification note", esc(c.classification_note || ""))}
    ${field("Protein/enzyme", esc(c.protein_name || ""))}
    ${field("EC / selection", `${esc(c.target_ec || "")} | ${esc(c.selection_type || "")} rank ${esc(c.rank || "")}`)}
    ${field("Organism", `${esc(c.organism || "")} | ${esc(c.source_bucket || "")} | ${esc(c.phylum_class || "")}`)}
    ${field("Length / sequence", `${esc(c.length || "")} aa | ${esc(c.sequence_status || "")}`)}
    ${field("Order sheet", `serial ${esc(c.order_serial || "")} | Excel row ${esc(c.order_excel_row || "")}`)}
    ${field("Nearest panel identity", `${esc(c.nearest_panel_identity_percent || "")}% vs ${esc(c.nearest_panel_identity_partner || "")}`)}
    ${field("Construct decision", `<span class="chip ${esc(c.construct_decision)}">${esc(c.construct_decision || "")}</span> score ${esc(c.construct_score || "")} | ${esc(c.construct_strategy || "")}`)}
    ${field("Construct evidence", esc(c.construct_reasons || ""))}
    ${field("Cautions", esc(c.construct_cautions || ""))}
  </div>${warn}`;
}
function structure(c){
  return `<div class="kv">
    ${field("Structure status", esc(c.structure_status || ""))}
    ${field("PDB IDs", pdbLinks(c.pdb_ids))}
    ${field("PDB DOI", doiLinks(c.pdb_dois))}
    ${field("PDB PubMed", esc(c.pdb_pubmeds || ""))}
    ${field("PDB title", esc(c.pdb_citation_titles || ""))}
    ${field("Method / resolution", esc([c.pdb_methods,c.pdb_resolutions].filter(Boolean).join(" | ")))}
    ${field("PDB source organism", esc(c.pdb_source_organisms || ""))}
    ${field("PDB expression host", esc(c.pdb_expression_hosts || ""))}
    ${field("Native/heterologous mode", esc(c.native_vs_heterologous_evidence || ""))}
    ${field("UniProt DOI", doiLinks(c.uniprot_dois))}
    ${field("UniProt titles", esc(c.uniprot_reference_titles || ""))}
    ${field("Expression hints", esc(c.uniprot_expression_reference_hints || ""))}
    ${field("Metal note", esc(c.metal_dependency || ""))}
  </div>`;
}
function fasta(c){
  return `<div class="seq-grid">
    <div><h3>Dashboard Protein FASTA</h3><pre>${esc(c.protein_fasta || "No sequence")}</pre></div>
    <div><h3>Order Sheet Translated Protein FASTA</h3><pre>${esc(c.order_protein_fasta || "")}</pre></div>
    <div><h3>Order Sheet Coding DNA FASTA</h3><pre>${esc(c.order_dna_fasta || "")}</pre></div>
  </div>`;
}
function colorFor(v){
  const x=Math.max(0,Math.min(100,+v||0));
  const stops=[[0,[182,83,72]],[35,[217,160,72]],[55,[232,220,145]],[78,[129,181,127]],[100,[37,123,128]]];
  for(let i=1;i<stops.length;i++){
    const [rv,rc]=stops[i], [lv,lc]=stops[i-1];
    if(x<=rv){const t=(x-lv)/Math.max(1,rv-lv); const mix=k=>Math.round(lc[k]+(rc[k]-lc[k])*t); return `rgb(${mix(0)}, ${mix(1)}, ${mix(2)})`;}
  }
  return "rgb(37,123,128)";
}
function pairFor(a,b){
  if(a===b){const c=byAcc(a); return {query:a,target:b,query_gene:gene(c),target_gene:gene(c),identity_percent:100,identical_sites:c.length,compared_sites:c.length};}
  return pairs.find(p=>(p.query===a&&p.target===b)||(p.query===b&&p.target===a));
}
function pairDetail(a,b){
  const p=pairFor(a,b)||{}; const ca=byAcc(a), cb=byAcc(b);
  return `<div>Pair</div><div><b>${esc(gene(ca))} vs ${esc(gene(cb))}</b></div>
    <div>Accessions</div><div>${esc(a)} vs ${esc(b)}</div>
    <div>Identity</div><div><span class="chip ${(+p.identity_percent||0)>=70?"green":"amber"}">${(+p.identity_percent||0).toFixed(2)}%</span></div>
    <div>Compared sites</div><div>${esc(p.identical_sites ?? "")} identical / ${esc(p.compared_sites ?? "")}</div>
    <div>A organism</div><div>${esc(ca.organism || "")}</div>
    <div>B organism</div><div>${esc(cb.organism || "")}</div>`;
}
function identity(c){
  const selectedIndex=idx(c.accession);
  const labels=candidates.map(gene);
  const ids=candidates.map(x=>x.accession);
  const head=`<tr><th class="corner">${ids.length} seq</th>${labels.map(x=>`<th class="col-head" title="${esc(x)}"><div>${esc(x)}</div></th>`).join("")}</tr>`;
  const body=ids.map((a,i)=>`<tr><th class="row-head" title="${esc(labels[i])}">${esc(labels[i])}</th>${ids.map((b,j)=>{
    const v=matrix[i]?.[j] ?? 0; const text=i===j||v>=70?Math.round(v):""; const cls=["cell",text!==""?"show":"",v>=90?"dark":"",(i===selectedIndex||j===selectedIndex)?"focus":""].filter(Boolean).join(" ");
    return `<td class="${cls}" style="background:${colorFor(v)}" data-a="${esc(a)}" data-b="${esc(b)}" title="${esc(labels[i])} vs ${esc(labels[j])}: ${(+v).toFixed(2)}%">${esc(text)}</td>`;
  }).join("")}</tr>`).join("");
  const rows=pairs.filter(p=>p.query===c.accession||p.target===c.accession).slice(0,36).map(p=>{
    const other=p.query===c.accession?p.target:p.query; const oc=byAcc(other);
    return `<div class="card"><b>${esc(gene(c))} vs ${esc(gene(oc))}</b><br>${(+p.identity_percent).toFixed(2)}% | ${esc(c.accession)} vs ${esc(other)}</div>`;
  }).join("");
  setTimeout(()=>{document.querySelectorAll(".cell").forEach(td=>td.addEventListener("click",()=>{$("pairDetail").innerHTML=pairDetail(td.dataset.a,td.dataset.b);})); $("pairDetail").innerHTML=pairDetail(c.accession,c.accession);},0);
  return `<div id="pairDetail" class="pair-detail"></div><div class="identity-wrap"><table class="identity-heatmap">${head}${body}</table></div>
    <div class="legend"><span>low</span><div class="ramp-wrap"><div class="ramp"></div><span class="tick" style="left:0%">0%</span><span class="tick" style="left:50%">50%</span><span class="tick" style="left:75%">75%</span><span class="tick" style="left:90%">90%</span><span class="tick" style="left:100%">100%</span></div><span>high</span></div>
    <h3>Pairs for ${esc(gene(c))}</h3>${rows || "<p class='meta'>No identity pairs.</p>"}`;
}
function orderSheet(c){
  const o=c.order_sheet||{};
  return `<div class="kv">
    ${field("Serial / Excel row", `${esc(o.serial || "")} / ${esc(o.excel_row || "")}`)}
    ${field("Gene name", esc(o.gene_name || gene(c)))}
    ${field("Order accession", esc(o.accession || c.accession))}
    ${field("5 prime site", esc(o.five_prime_site || ""))}
    ${field("3 prime site", esc(o.three_prime_site || ""))}
    ${field("Optimization host", esc(o.optimization_host || ""))}
    ${field("Vector", esc(o.standard_vector || ""))}
    ${field("Vector name", esc(o.vector_name || ""))}
    ${field("Antibiotic", esc(o.antibiotic || ""))}
    ${field("DNA / translated protein", `${esc(o.dna_length || "")} nt / ${esc(o.protein_length_from_order || "")} aa`)}
    ${field("Plasmid prep", esc(o.plasmid_prep || ""))}
    ${field("Endotoxin", esc(o.endotoxin || ""))}
    ${field("Order note", esc(o.order_note || ""))}
  </div>`;
}
function pdbFiles(c){
  const ids=splitSemi(c.pdb_ids);
  const rows=(DATA.pdbFiles||[]).filter(p=>ids.includes(p.pdb_id));
  if(!rows.length) return "<p class='meta'>No downloaded PDB file for this candidate.</p>";
  return rows.map(r=>`<div class="kv">${field("PDB",esc(r.pdb_id))}${field("Status",esc(r.status))}${field("File",esc(r.file))}${field("RCSB",`<a target="_blank" href="${esc(r.rcsb_url)}">${esc(r.rcsb_url)}</a>`)}</div>`).join("");
}
function notes(c){
  const key="hoa_type_ii_aldolase_note_"+c.accession;
  const val=localStorage.getItem(key)||"";
  setTimeout(()=>{const box=document.getElementById("note"); if(box) box.addEventListener("input",()=>localStorage.setItem(key,box.value));},0);
  return `<textarea id="note" placeholder="Browser-local note for construct design or ordering">${esc(val)}</textarea>`;
}
function tabButton(id,label){return `<button class="${tab===id?"active":""}" data-tab="${id}">${esc(label)}</button>`;}
function renderDetail(){
  const c=byAcc(selected);
  if(!c){$("detail").innerHTML="<p>No candidates.</p>";return;}
  const body=tab==="structure"?structure(c):tab==="fasta"?fasta(c):tab==="identity"?identity(c):tab==="order"?orderSheet(c):tab==="pdbfiles"?pdbFiles(c):tab==="notes"?notes(c):overview(c);
  $("detail").innerHTML=`<div class="detail-head"><h2>${esc(gene(c))} <span class="chip ${esc(c.construct_decision)}">${esc(c.construct_decision||"")}</span></h2>
    <div class="meta">${esc(c.accession)} | ${esc(c.aldolase_class)} | ${esc(c.primary_activity)}</div>
    <div class="chips"><span class="chip green">metal-dependent</span><span class="chip">${esc(c.target_ec||"")}</span><span class="chip">row ${esc(c.order_excel_row||"")}</span></div></div>
    <nav class="tabs">${tabButton("overview","Overview")}${tabButton("structure","PDB/DOI/Expression")}${tabButton("fasta","FASTA")}${tabButton("identity","Identity Matrix")}${tabButton("order","Order Sheet")}${tabButton("pdbfiles","PDB files")}${tabButton("notes","Notes")}</nav>
    <div class="content">${body}</div>`;
  document.querySelectorAll(".tabs button").forEach(btn=>btn.addEventListener("click",()=>{tab=btn.dataset.tab;renderDetail();}));
}
function render(){
  $("title").textContent=DATA.title || "HOA Type II Aldolase Candidate Dashboard";
  $("subtitle").innerHTML=[`${DATA.summary.candidate_count} order-sheet candidates`,`Type II aldolase / metal-dependent HOA-family`,`range ${esc(DATA.range.start_gene)} to ${esc(DATA.range.end_gene)}`].map(x=>`<span>${x}</span>`).join("");
  renderList();
  renderDetail();
}
function setup(){
  fillSelect("decision","All decisions",uniq("construct_decision"));
  fillSelect("source","All sources",uniq("source_bucket"));
  fillSelect("structure","All structures",uniq("structure_status"));
  ["q","decision","source","structure","sort"].forEach(id=>$(id).addEventListener("input",render));
  render();
}
setup();
</script>
</body>
</html>
"""


def main() -> None:
    parser = argparse.ArgumentParser(description="Update the HOA dashboard from the Aldolase/Fluorinase plasmid-order XLS.")
    sub = parser.add_subparsers(dest="command", required=True)
    build_parser = sub.add_parser("build")
    build_parser.add_argument("--html", type=Path, required=True)
    build_parser.add_argument("--xls", type=Path, required=True)
    build_parser.add_argument("--outdir", type=Path, required=True)
    build_parser.add_argument("--title", default="HOA Type II Aldolase Candidate Dashboard")
    build_parser.add_argument("--start-gene", default="dmpG")
    build_parser.add_argument("--end-gene", default="BPHF_RHOJR")
    check_parser = sub.add_parser("check")
    check_parser.add_argument("--outdir", type=Path, required=True)
    check_parser.add_argument("--expected-candidates", type=int, default=18)
    args = parser.parse_args()
    result = build(args) if args.command == "build" else check(args)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
