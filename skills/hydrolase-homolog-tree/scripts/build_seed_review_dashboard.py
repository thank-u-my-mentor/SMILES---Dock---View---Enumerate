#!/usr/bin/env python3
"""
Build a small dashboard for reviewing seed enzyme provenance, expression setup,
and optional solubility annotations.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
from pathlib import Path
from typing import Dict, List


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def key_value(row: Dict[str, str], *names: str) -> str:
    lower = {k.lower(): k for k in row}
    for name in names:
        key = lower.get(name.lower())
        if key and str(row.get(key, "")).strip():
            return str(row[key]).strip()
    return ""


FIELD_ALIASES = {
    "family": ["family", "Family"],
    "enzyme": ["enzyme_name", "酶名称", "enzyme", "id"],
    "organism": ["source_organism", "来源生物", "organism"],
    "promiscuity": ["known_promiscuity", "已知混杂性"],
    "uniprot": ["uniprot_id", "Uniprot id", "UniProt ID"],
    "doi": ["key_doi", "关键DOI", "doi", "DOI"],
}


def display_value(row: Dict[str, str], logical_key: str) -> str:
    return key_value(row, *FIELD_ALIASES.get(logical_key, [logical_key]))


def merge_solubility(records: List[Dict[str, str]], sol_rows: List[Dict[str, str]]) -> None:
    if not sol_rows:
        return
    by_id = {}
    for row in sol_rows:
        for key in [
            key_value(row, "id"),
            key_value(row, "source_id"),
            key_value(row, "uniprot_id", "Uniprot id"),
        ]:
            if key:
                by_id[key.strip()] = row
    for rec in records:
        rec_id = key_value(rec, "Uniprot id", "uniprot_id", "id")
        sol = by_id.get(rec_id)
        if not sol:
            continue
        rec["solubility_score"] = key_value(sol, "solubility_score", "score", "probability")
        rec["solubility_call"] = key_value(sol, "solubility_call", "prediction", "call")


def write_dashboard(records: List[Dict[str, str]], output: Path) -> None:
    data = json.dumps(records, ensure_ascii=False)
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Seed Enzyme Review</title>
<style>
body {{ margin:0; font-family: Arial, sans-serif; background:#f7f8fa; color:#1f2937; }}
header {{ padding:18px 24px; background:#ffffff; border-bottom:1px solid #d8dee9; position:sticky; top:0; z-index:2; }}
h1 {{ margin:0 0 10px; font-size:20px; }}
.controls {{ display:flex; gap:10px; flex-wrap:wrap; }}
input, select {{ padding:8px 10px; border:1px solid #b8c0cc; border-radius:4px; background:white; }}
main {{ display:grid; grid-template-columns: 440px 1fr; gap:16px; padding:16px; }}
.list {{ display:flex; flex-direction:column; gap:8px; }}
.item {{ background:white; border:1px solid #d8dee9; border-radius:6px; padding:10px; cursor:pointer; }}
.item.active {{ outline:2px solid #2563eb; }}
.title {{ font-weight:700; font-size:14px; }}
.meta {{ color:#5b6472; font-size:12px; margin-top:4px; }}
.tags {{ display:flex; gap:6px; flex-wrap:wrap; margin-top:8px; }}
.tag {{ border-radius:999px; padding:2px 7px; font-size:12px; background:#e5e7eb; }}
.good {{ background:#dcfce7; color:#166534; }}
.warn {{ background:#fef3c7; color:#92400e; }}
.bad {{ background:#fee2e2; color:#991b1b; }}
.panel {{ background:white; border:1px solid #d8dee9; border-radius:6px; padding:16px; min-height:420px; }}
.grid {{ display:grid; grid-template-columns: 180px 1fr; gap:8px 14px; align-items:start; }}
.k {{ color:#5b6472; font-size:12px; text-transform:uppercase; }}
.v {{ white-space:pre-wrap; overflow-wrap:anywhere; }}
a {{ color:#1d4ed8; }}
@media (max-width: 900px) {{ main {{ grid-template-columns:1fr; }} }}
</style>
</head>
<body>
<header>
<h1>Seed Enzyme Review</h1>
<div class="controls">
<input id="q" placeholder="Search enzyme, organism, DOI, vector, tag">
<select id="confidence">
<option value="">All confidence</option>
<option value="high">High</option>
<option value="medium">Medium</option>
<option value="low">Low</option>
<option value="missing">Missing evidence</option>
</select>
<select id="vector">
<option value="">All vectors</option>
<option value="pet28">pET28-like</option>
<option value="pet23">pET23-like</option>
<option value="other">Other/unknown</option>
</select>
</div>
</header>
<main>
<section class="list" id="list"></section>
<section class="panel" id="detail"></section>
</main>
<script>
const DATA = {data};
const fields = [
  ['family','Family'], ['enzyme','Enzyme'], ['organism','Organism'], ['uniprot','UniProt'],
  ['doi','DOI'], ['paper_title','Paper title'], ['expression_host','Host'],
  ['expression_strain','Strain'], ['plasmid_vector','Vector'], ['promoter','Promoter'],
  ['antibiotic','Antibiotic'], ['tag_type','Tag'], ['tag_position','Tag position'],
  ['tag_cleaved','Tag cleaved'], ['expression_temperature_c','Temp C'], ['inducer','Inducer'],
  ['soluble_expression_reported','Soluble reported'], ['solubility_call','Predicted soluble'],
  ['solubility_score','Solubility score'], ['purification_method','Purification'],
  ['evidence_quote','Evidence quote'], ['evidence_source','Evidence source'],
  ['provenance_confidence','Confidence'], ['notes_for_cloning','Notes']
];
let selected = 0;
function esc(x) {{ return String(x || '').replace(/[&<>"']/g, c => ({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}}[c])); }}
const aliases = {{
  family: ['family','Family'],
  enzyme: ['enzyme_name','酶名称','enzyme','id'],
  organism: ['source_organism','来源生物','organism'],
  promiscuity: ['known_promiscuity','已知混杂性'],
  uniprot: ['uniprot_id','Uniprot id','UniProt ID'],
  doi: ['key_doi','关键DOI','doi','DOI']
}};
function value(row, key) {{
  const list = aliases[key] || [key];
  for (const k of list) if (row[k]) return row[k];
  return '';
}}
function doiLink(doi) {{ return doi ? `<a target="_blank" href="https://doi.org/${{esc(doi)}}">${{esc(doi)}}</a>` : ''; }}
function cls(row) {{
  const c = (row.provenance_confidence || '').toLowerCase();
  if (c === 'high') return 'good';
  if (c === 'medium') return 'warn';
  if (c === 'low') return 'bad';
  return 'warn';
}}
function passes(row) {{
  const q = document.getElementById('q').value.toLowerCase();
  const conf = document.getElementById('confidence').value;
  const vf = document.getElementById('vector').value;
  const blob = JSON.stringify(row).toLowerCase();
  if (q && !blob.includes(q)) return false;
  const c = (row.provenance_confidence || '').toLowerCase();
  if (conf === 'missing' && c) return false;
  if (conf && conf !== 'missing' && c !== conf) return false;
  const v = (row.plasmid_vector || '').toLowerCase();
  if (vf === 'pet28' && !v.includes('pet-28') && !v.includes('pet28')) return false;
  if (vf === 'pet23' && !v.includes('pet-23') && !v.includes('pet23')) return false;
  if (vf === 'other' && (v.includes('pet-28') || v.includes('pet28') || v.includes('pet-23') || v.includes('pet23'))) return false;
  return true;
}}
function render() {{
  const rows = DATA.filter(passes);
  const list = document.getElementById('list');
  list.innerHTML = rows.map((r, i) => `
    <div class="item ${{i===selected?'active':''}}" onclick="selected=${{i}}; renderDetail(DATA.filter(passes)[selected]); render();">
      <div class="title">${{esc(value(r,'enzyme') || 'Untitled')}}</div>
      <div class="meta">${{esc(value(r,'organism'))}} · ${{esc(value(r,'uniprot'))}}</div>
      <div class="tags">
        <span class="tag ${{cls(r)}}">${{esc(r.provenance_confidence || 'missing evidence')}}</span>
        <span class="tag">${{esc(r.plasmid_vector || 'vector ?')}}</span>
        <span class="tag">${{esc(r.tag_position || 'tag ?')}}</span>
        <span class="tag">${{esc(r.solubility_call || 'solubility ?')}}</span>
      </div>
    </div>`).join('');
  if (selected >= rows.length) selected = 0;
  renderDetail(rows[selected]);
}}
function renderDetail(row) {{
  const d = document.getElementById('detail');
  if (!row) {{ d.innerHTML = 'No records match.'; return; }}
  d.innerHTML = '<div class="grid">' + fields.map(([key,label]) => {{
    const raw = value(row, key);
    const val = key === 'doi' ? doiLink(raw) : esc(raw);
    return `<div class="k">${{esc(label)}}</div><div class="v">${{val}}</div>`;
  }}).join('') + '</div>';
}}
for (const id of ['q','confidence','vector']) document.getElementById(id).addEventListener('input', () => {{ selected=0; render(); }});
render();
</script>
</body>
</html>
"""
    output.write_text(page, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build seed provenance review dashboard.")
    parser.add_argument("seed_provenance_csv", help="CSV produced by extract_seed_metadata_from_xlsx.py and manually/LLM enriched.")
    parser.add_argument("--solubility-csv", default=None, help="Optional normalized solubility CSV.")
    parser.add_argument("--output", default="seed_review_dashboard.html")
    args = parser.parse_args()

    records = read_csv(Path(args.seed_provenance_csv))
    if args.solubility_csv:
        merge_solubility(records, read_csv(Path(args.solubility_csv)))
    write_dashboard(records, Path(args.output))
    print(f"Wrote {args.output} with {len(records)} records.")


if __name__ == "__main__":
    main()
