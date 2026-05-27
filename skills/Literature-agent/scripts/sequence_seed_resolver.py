#!/usr/bin/env python3
"""Build a human-readable sequence seed review table from literature-agent CSVs.

Default mode is offline: infer organism/source clues and concrete enzyme IDs from
existing CSV evidence, then write an HTML review dashboard plus an enhanced CSV.
Network sequence fetching is intentionally not enabled here yet; use the emitted
database query links for manual review or add explicit fetch modes later.
"""

from __future__ import annotations

import argparse
import csv
import html
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.parse import quote_plus


GENERIC_NAME_PATTERNS = [
    r"\benzymes?\b",
    r"\bfamily members?\b",
    r"\bpanel\b",
    r"\bspecific variants? not named\b",
    r"\bspecific enzyme unclear\b",
    r"\bspecific not named\b",
    r"\brounds? of protein engineering\b",
    r"\bprotein engineering\b",
    r"\bFMNsq\b",
    r"\bflavin semiquinone\b",
    r"\bTyr\s*[->→]\s*Phe\b",
    r"\bCys-incorporated\b",
    r"\bvariants?\b",
    r"\bmutants?\b",
    r"\bwild type\b",
    r"\bapo-",
    r"\bholo-",
    r"\bbiosynthesis enzyme\b",
    r"\bflavoenzymes?\b",
    r"\boxidases? and dehydrogenases?\b",
    r"\bdehydrogenases?\b",
    r"\breductases?\b",
]

REVIEW_HINTS = [
    "review", "strategies", "advances", "versatile catalysts", "biochemical implications",
    "accounts", "perspective", "concept", "on their way",
]

ORGANISM_ALIASES = {
    "CvFAP": "Chlorella variabilis",
    "Chlorella variabilis FAP": "Chlorella variabilis",
    "GluER": "Gluconobacter",
    "GsOYE": "Galdieria sulphuraria",
    "OYE1": "Saccharomyces pastorianus",
    "PaDADH": "Pseudomonas aeruginosa",
    "Pseudomonas aeruginosa D-arginine dehydrogenase": "Pseudomonas aeruginosa",
    "avenolide biosynthetic flavoenzyme": "Streptomyces avermitilis",
}

KNOWN_ORGANISM_NAMES = [
    "Bacillus amyloliquefaciens",
    "Enterobacter cloacae",
    "Galdieria sulphuraria",
    "Chlorella variabilis",
    "Gluconobacter",
    "Saccharomyces pastorianus",
    "Nicotiana tabacum",
    "Escherichia coli",
    "Pseudomonas aeruginosa",
    "Streptomyces avermitilis",
]

NAME_PATTERNS = [
    r"\b[A-Z][a-z]{1,4}ER(?:-[A-Z0-9]+)*\b",
    r"\b[A-Z][a-z]OYE\d*(?:-[A-Z0-9]+)*\b",
    r"\bOYE\d+(?:-[A-Z0-9]+)*\b",
    r"\bCvFAP(?:-[A-Z][0-9]+[A-Z])*\b",
    r"\bFAP(?:-[A-Z][0-9]+[A-Z])*\b",
    r"\b[A-Z][a-z]{1,5}FDH(?:-[A-Z0-9]+)*\b",
    r"\b[A-Z][a-z]{1,5}NTR(?:-[A-Z0-9]+)*\b",
    r"\bPaDADH(?:-[A-Z0-9]+)*\b",
    r"\bPseudomonas aeruginosa D-arginine dehydrogenase\b",
    r"\bovenolide biosynthetic flavoenzyme\b",
    r"\b[A-Z][A-Za-z0-9]{1,8}(?:ase|ER|OYE|FAP|FDH|NTR)(?:-[A-Z0-9]+)*\b",
]

MUTATION_RE = re.compile(r"\b[A-Z][0-9]{1,4}[A-Z]\b")


def read_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def write_csv(path: Path, rows: Sequence[Dict[str, str]], fields: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def norm_space(value: str) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip()


def split_multi(value: str) -> List[str]:
    parts = re.split(r"\s*\|\s*|[;,]\s*", str(value or ""))
    return [norm_space(p) for p in parts if norm_space(p)]


def clean_candidate_name(value: str) -> str:
    value = norm_space(value)
    value = re.sub(r"^\(?\s*(engineered|evolved|directed evolution|mutant|variant)\s+", "", value, flags=re.I)
    value = re.sub(r"\s+from\s+[A-Z][a-z]+(?:\s+[a-z][a-z-]+){0,3}\b", "", value)
    value = re.sub(r"\bspecific\s+variants?\s+not\s+named(?:\s+in\s+abstract)?\b", "specific not named", value, flags=re.I)
    value = re.sub(r"\b(?:mutants?|variants?)\b", "", value, flags=re.I)
    value = re.sub(r"\b[A-Z][0-9]{1,4}[A-Z]\s+mutant\b", "", value)
    value = re.sub(r"\bTyr\s*[->→]\s*Phe/Trp\s+mutants?\b", "", value, flags=re.I)
    value = value.strip(" .,:;()[]")
    value = norm_space(value)
    return value


def is_generic_name(value: str, title: str = "") -> bool:
    low = f"{value} {title}".lower()
    if not value or value.lower() in {"unclear", "nan", "none", "n/a", "specific not named", "specific not named in abstract"}:
        return True
    if any(phrase in value.lower() for phrase in [
        "rounds of protein engineering", "specific not named", "specific enzyme unclear",
        "flavin semiquinone", "fmnsq", "cys-incorporated", "tyr→phe", "tyr->phe",
    ]):
        return True
    if re.fullmatch(r"[A-Z][0-9]{1,4}[A-Z](?:\s+mutant)?", value):
        return True
    if re.fullmatch(r"(?:specific\s+)?variants?(?:\s+not\s+named(?:\s+in\s+abstract)?)?", value, flags=re.I):
        return True
    if re.fullmatch(r"(?:mutants?|wild type|apo-[A-Za-z0-9]+|holo-[A-Za-z0-9]+)", value, flags=re.I):
        return True
    if any(re.search(pattern, value.lower()) for pattern in GENERIC_NAME_PATTERNS):
        return True
    if len(value) > 120 or value.count("|") >= 2:
        return True
    if any(hint in low for hint in REVIEW_HINTS):
        return True
    return False


def extract_specific_names(row: Dict[str, str]) -> List[str]:
    text = " ".join(
        norm_space(row.get(k, ""))
        for k in ["enzyme_name", "enzyme_name_or_target", "title", "abstract", "evidence_summary", "metadata_keywords"]
    )
    names: List[str] = []
    for field in ["enzyme_name", "enzyme_name_or_target"]:
        for part in split_multi(row.get(field, "")):
            cleaned = clean_candidate_name(part)
            if cleaned and not is_generic_name(cleaned):
                names.append(cleaned)
    for pattern in NAME_PATTERNS:
        for match in re.findall(pattern, text):
            cleaned = clean_candidate_name(match)
            if cleaned and not is_generic_name(cleaned):
                names.append(cleaned)
    return prefer_specific_names(dedupe(names))[:10]


def prefer_specific_names(names: Sequence[str]) -> List[str]:
    scored = []
    for name in names:
        low = name.lower()
        score = 0
        if re.search(r"\b(?:[A-Z][a-z]{1,4}ER|[A-Z][a-z]OYE\d*|OYE\d+|CvFAP|GsOYE|OaER)\b", name):
            score -= 4
        if any(alias.lower() in low for alias in ORGANISM_ALIASES):
            score -= 2
        if " from " in low:
            score -= 1
        if any(generic in low for generic in ["fatty acid photodecarboxylase", "photodecarboxylase", "ered/oye"]):
            score += 2
        if len(name) > 70:
            score += 3
        scored.append((score, name))
    return [name for _, name in sorted(scored, key=lambda item: (item[0], item[1].lower()))]


def parent_enzyme_name(name: str) -> str:
    value = norm_space(name)
    value = re.sub(r"\s*\([^)]*\)", "", value)
    value = re.sub(r"\s+from\s+[A-Z][a-z]+(?:\s+[a-z][a-z-]+){0,3}\b", "", value)
    value = re.sub(r"(?:-[A-Z][0-9]{1,4}[A-Z])+$", "", value)
    return value.strip(" -")


def extract_mutations(names: Sequence[str], row: Dict[str, str]) -> List[str]:
    text = " ".join(list(names) + [row.get("enzyme_name", ""), row.get("enzyme_name_or_target", ""), row.get("abstract", ""), row.get("evidence_summary", "")])
    return dedupe(MUTATION_RE.findall(text))


def extract_organisms(row: Dict[str, str], names: Sequence[str]) -> Tuple[List[str], List[str]]:
    text = " ".join(
        norm_space(row.get(k, ""))
        for k in ["enzyme_name", "enzyme_name_or_target", "title", "abstract", "evidence_summary", "metadata_keywords"]
    )
    organisms: List[str] = []
    evidence: List[str] = []
    for name in names:
        for key, organism in ORGANISM_ALIASES.items():
            if key.lower() in name.lower() or key.lower() in text.lower():
                organisms.append(organism)
                evidence.append(f"alias:{key}")
    for organism in KNOWN_ORGANISM_NAMES:
        if organism.lower() in text.lower():
            organisms.append(organism)
            evidence.append(f"text:{organism}")
    for match in re.findall(r"\bfrom\s+([A-Z][a-z]+(?:\s+[a-z][a-z-]+){0,2})", text):
        if not looks_like_false_organism(match):
            organisms.append(match)
            evidence.append(f"from:{match}")
    for match in re.findall(r"\b([A-Z][a-z]+(?:\s+[a-z][a-z-]+){0,2})\s+(?:FAP|OYE|ERED|ene-reductase|nitroreductase|monooxygenase|halogenase)\b", text):
        if not looks_like_false_organism(match):
            organisms.append(match)
            evidence.append(f"prefix:{match}")
    return dedupe(organisms), dedupe(evidence)


def looks_like_false_organism(value: str) -> bool:
    low = value.lower().strip()
    bad = {
        "the", "this", "these", "photoenzymatic", "flavin", "fatty acid", "old yellow",
        "single electron", "directed evolution", "natural redox", "functional group",
        "journal of", "nature communications",
        "enzyme", "mass spectra of", "photodecarboxylase", "fatty acid photodecarboxylase",
        "th", "abstract photoenzymatic", "light footprinting", "collision cross",
        "nitronate", "enterobacter cloacae enables", "saccharomyces cerevisiae and",
    }
    if len(low) < 5:
        return True
    if low in bad:
        return True
    if any(low.startswith(prefix) for prefix in ["the ", "this ", "these ", "mass ", "abstract "]):
        return True
    if any(word in low for word in ["spectra", "reaction", "substrate", "protein engineering"]):
        return True
    return False


def dedupe(values: Iterable[str]) -> List[str]:
    seen = set()
    out = []
    for value in values:
        value = norm_space(value)
        key = value.lower()
        if not value or key in seen:
            continue
        seen.add(key)
        out.append(value)
    return out


def classify_row(row: Dict[str, str], names: Sequence[str], organisms: Sequence[str]) -> Tuple[str, str]:
    title = row.get("title", "")
    base_name = row.get("enzyme_name") or row.get("enzyme_name_or_target") or ""
    if any(hint in title.lower() for hint in REVIEW_HINTS):
        return "review_or_family_context", "Review/background paper; use for leads, not direct FASTA seed."
    if not names and is_generic_name(base_name, title):
        return "family_only", "No concrete enzyme symbol/name extracted."
    if names and organisms:
        return "specific_with_organism", "Concrete enzyme plus organism/source evidence."
    if names:
        return "specific_missing_organism", "Concrete enzyme candidate but organism/source needs database or paper check."
    return "needs_manual_review", "Insufficient sequence-resolution evidence."


def query_urls(name: str, organism: str) -> Dict[str, str]:
    terms = " ".join(t for t in [parent_enzyme_name(name), organism] if t)
    protein_terms = terms or name or organism
    dna_terms = " ".join(t for t in [parent_enzyme_name(name), organism, "gene"] if t)
    return {
        "uniprot_url": "https://www.uniprot.org/uniprotkb?query=" + quote_plus(protein_terms),
        "ncbi_protein_url": "https://www.ncbi.nlm.nih.gov/protein/?term=" + quote_plus(protein_terms),
        "ncbi_nucleotide_url": "https://www.ncbi.nlm.nih.gov/nuccore/?term=" + quote_plus(dna_terms),
        "ena_url": "https://www.ebi.ac.uk/ena/browser/text-search?query=" + quote_plus(dna_terms),
    }


def combine_rows(enzyme_rows: List[Dict[str, str]], homolog_rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    rows = []
    by_doi = defaultdict(dict)
    for row in homolog_rows:
        doi = norm_space(row.get("doi", "")).lower()
        if doi:
            by_doi[doi].update(row)
    for row in enzyme_rows:
        merged = dict(by_doi.get(norm_space(row.get("doi", "")).lower(), {}))
        merged.update(row)
        rows.append(merged)
    for row in homolog_rows:
        doi = norm_space(row.get("doi", "")).lower()
        if doi and not any(norm_space(r.get("doi", "")).lower() == doi for r in enzyme_rows):
            rows.append(dict(row))
    return rows


def build_candidates(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    out = []
    for row in rows:
        names = extract_specific_names(row)
        mutations = extract_mutations(names, row)
        organisms, organism_evidence = extract_organisms(row, names)
        status, reason = classify_row(row, names, organisms)
        name_for_query = parent_enzyme_name(names[0]) if names else clean_candidate_name(row.get("enzyme_name") or row.get("enzyme_name_or_target", ""))
        organism_for_query = organisms[0] if organisms else ""
        urls = query_urls(name_for_query, organism_for_query)
        out.append({
            "resolution_status": status,
            "resolution_reason": reason,
            "specific_enzyme_names": " | ".join(names),
            "parent_enzyme_names": " | ".join(dedupe(parent_enzyme_name(n) for n in names)),
            "mutations_or_variants": " | ".join(mutations),
            "organism_candidates": " | ".join(organisms),
            "organism_evidence": " | ".join(organism_evidence),
            "protein_query": " ".join(t for t in [name_for_query, organism_for_query] if t),
            "dna_query": " ".join(t for t in [name_for_query, organism_for_query, "gene"] if t),
            **urls,
            "enzyme_family": row.get("enzyme_family", ""),
            "enzyme_name_original": row.get("enzyme_name") or row.get("enzyme_name_or_target", ""),
            "doi": row.get("doi", ""),
            "title": row.get("title", ""),
            "year": row.get("year", ""),
            "journal": row.get("journal", ""),
            "search_track": row.get("search_track", ""),
            "flavin_cofactor": row.get("flavin_cofactor", ""),
            "reaction_type": row.get("reaction_type", ""),
            "new_to_nature_reaction": row.get("new_to_nature_reaction", ""),
            "characterized_enzyme_evidence": row.get("characterized_enzyme_evidence", ""),
            "structure_similarity_hint": row.get("structure_similarity_hint", ""),
            "criteria_status": row.get("criteria_status", ""),
            "manual_review_priority": row.get("manual_review_priority", row.get("seed_priority", "")),
            "pdb_ids": row.get("pdb_ids", ""),
            "uniprot_id": row.get("uniprot_id", ""),
            "evidence_summary": row.get("evidence_summary") or row.get("abstract", ""),
        })
    return sorted(dedupe_candidates(out), key=sort_key)


def canonical_enzyme_name(row: Dict[str, str]) -> str:
    blob = " | ".join([
        row.get("specific_enzyme_names", ""),
        row.get("parent_enzyme_names", ""),
        row.get("enzyme_name_original", ""),
    ])
    priority = [
        r"\bCvFAP\b", r"\bGluER\b", r"\bGsOYE\b", r"\bOaER\b", r"\bOYE1\b", r"\bPaDADH\b",
        r"\bPqsL\b", r"\b[A-Z][a-z]{1,4}ER\b", r"\b[A-Z][a-z]OYE\d*\b",
        r"\b[A-Z][a-z]{1,5}FDH\b", r"\b[A-Z][a-z]{1,5}NTR\b",
        r"\bPseudomonas aeruginosa D-arginine dehydrogenase\b",
        r"\bovenolide biosynthetic flavoenzyme\b",
    ]
    for pattern in priority:
        match = re.search(pattern, blob)
        if match:
            return match.group(0)
    parents = split_multi(row.get("parent_enzyme_names", ""))
    if parents:
        return parents[0]
    names = split_multi(row.get("specific_enzyme_names", ""))
    if names:
        return parent_enzyme_name(names[0])
    return parent_enzyme_name(row.get("enzyme_name_original", ""))


def dedupe_candidates(rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    grouped: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    status_rank = {
        "specific_with_organism": 0,
        "specific_missing_organism": 1,
        "needs_manual_review": 2,
        "family_only": 3,
        "review_or_family_context": 4,
    }
    for row in rows:
        canonical = canonical_enzyme_name(row)
        organism = split_multi(row.get("organism_candidates", ""))
        organism_key = organism[0] if organism else ""
        key = (norm_space(row.get("doi", "")).lower(), canonical.lower(), organism_key.lower())
        if key not in grouped:
            merged = dict(row)
            merged["canonical_enzyme"] = canonical
            grouped[key] = merged
            continue
        current = grouped[key]
        if status_rank.get(row["resolution_status"], 9) < status_rank.get(current["resolution_status"], 9):
            current["resolution_status"] = row["resolution_status"]
            current["resolution_reason"] = row["resolution_reason"]
        for field in [
            "specific_enzyme_names", "parent_enzyme_names", "mutations_or_variants",
            "organism_candidates", "organism_evidence", "pdb_ids", "uniprot_id",
        ]:
            current[field] = " | ".join(dedupe(split_multi(current.get(field, "")) + split_multi(row.get(field, ""))))
        if len(norm_space(row.get("evidence_summary", ""))) > len(norm_space(current.get("evidence_summary", ""))):
            current["evidence_summary"] = row.get("evidence_summary", "")
        if not current.get("title") and row.get("title"):
            current["title"] = row["title"]
    for row in grouped.values():
        canonical = row.get("canonical_enzyme") or canonical_enzyme_name(row)
        organisms = split_multi(row.get("organism_candidates", ""))
        organism = organisms[0] if organisms else ""
        urls = query_urls(canonical, organism)
        row.update(urls)
        row["protein_query"] = " ".join(t for t in [canonical, organism] if t)
        row["dna_query"] = " ".join(t for t in [canonical, organism, "gene"] if t)
    return list(grouped.values())


def sort_key(row: Dict[str, str]) -> Tuple[int, int, int, str]:
    status_score = {
        "specific_with_organism": 0,
        "specific_missing_organism": 1,
        "needs_manual_review": 2,
        "family_only": 3,
        "review_or_family_context": 4,
    }.get(row["resolution_status"], 5)
    priority_score = 0 if row.get("manual_review_priority") == "high" else 1
    new_score = 0 if str(row.get("new_to_nature_reaction", "")).lower() == "true" else 1
    return (status_score, priority_score, new_score, row.get("title", ""))


def render_html(rows: Sequence[Dict[str, str]], output_path: Path) -> None:
    counts = defaultdict(int)
    for row in rows:
        counts[row["resolution_status"]] += 1
    cards = "\n".join(
        f"<div class='stat'><b>{html.escape(k)}</b><span>{v}</span></div>"
        for k, v in sorted(counts.items())
    )
    table_rows = "\n".join(render_row(row) for row in rows)
    doc = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Sequence Seed Review</title>
<style>
:root {{ color-scheme: light; --ink:#1f2933; --muted:#64748b; --line:#d7dee8; --bg:#f7f9fb; --panel:#ffffff; --accent:#006d77; --warn:#a15c00; --bad:#a11d33; }}
body {{ margin:0; font:14px/1.45 system-ui, -apple-system, Segoe UI, sans-serif; color:var(--ink); background:var(--bg); }}
header {{ padding:22px 28px 14px; background:var(--panel); border-bottom:1px solid var(--line); position:sticky; top:0; z-index:4; }}
h1 {{ margin:0 0 8px; font-size:22px; letter-spacing:0; }}
.sub {{ color:var(--muted); max-width:1100px; }}
.stats {{ display:flex; gap:10px; flex-wrap:wrap; margin-top:14px; }}
.stat {{ border:1px solid var(--line); background:#fbfdff; padding:8px 10px; border-radius:6px; min-width:160px; display:flex; justify-content:space-between; gap:16px; }}
.toolbar {{ display:flex; gap:10px; flex-wrap:wrap; padding:12px 28px; background:#eef4f6; border-bottom:1px solid var(--line); }}
input, select {{ height:34px; border:1px solid var(--line); border-radius:6px; padding:0 10px; background:white; }}
input {{ min-width:320px; flex:1; }}
main {{ padding:18px 28px 32px; }}
table {{ border-collapse:collapse; width:100%; background:var(--panel); border:1px solid var(--line); }}
th, td {{ border-bottom:1px solid var(--line); vertical-align:top; padding:9px 10px; }}
th {{ text-align:left; background:#f1f5f8; position:sticky; top:109px; z-index:3; font-size:12px; color:#334155; }}
tr:hover td {{ background:#fbfdff; }}
.status {{ display:inline-block; padding:2px 7px; border-radius:999px; font-size:12px; border:1px solid var(--line); background:#f8fafc; white-space:nowrap; }}
.specific_with_organism {{ color:#075e54; border-color:#86c5b9; background:#eaf8f5; }}
.specific_missing_organism {{ color:#7a4d00; border-color:#e4bd75; background:#fff7e6; }}
.family_only, .review_or_family_context {{ color:#7c2432; border-color:#e5a3af; background:#fff0f3; }}
.small {{ color:var(--muted); font-size:12px; }}
.title {{ max-width:340px; }}
.evidence {{ max-width:420px; }}
a {{ color:var(--accent); text-decoration:none; }}
a:hover {{ text-decoration:underline; }}
.links a {{ display:block; margin-bottom:3px; }}
</style>
</head>
<body>
<header>
  <h1>Sequence Seed Review</h1>
  <div class="sub">Offline review of enzyme-to-sequence readiness. Use specific enzyme names and organism/source evidence before trusting any protein or DNA FASTA hit.</div>
  <div class="stats">{cards}</div>
</header>
<div class="toolbar">
  <input id="q" placeholder="Search title, enzyme, organism, DOI, reaction">
  <select id="status"><option value="">All statuses</option>{''.join(f'<option>{html.escape(k)}</option>' for k in sorted(counts))}</select>
</div>
<main>
<table id="tbl">
<thead><tr>
<th>Status</th><th>Specific Enzyme</th><th>Organism Evidence</th><th>Queries</th><th>Paper / Evidence</th><th>Reaction</th>
</tr></thead>
<tbody>
{table_rows}
</tbody>
</table>
</main>
<script>
const q = document.getElementById('q');
const status = document.getElementById('status');
const rows = [...document.querySelectorAll('tbody tr')];
function apply() {{
  const needle = q.value.toLowerCase();
  const st = status.value;
  rows.forEach(row => {{
    const okText = !needle || row.innerText.toLowerCase().includes(needle);
    const okStatus = !st || row.dataset.status === st;
    row.style.display = okText && okStatus ? '' : 'none';
  }});
}}
q.addEventListener('input', apply);
status.addEventListener('input', apply);
</script>
</body></html>"""
    output_path.write_text(doc, encoding="utf-8")


def render_row(row: Dict[str, str]) -> str:
    doi = row.get("doi", "")
    doi_link = f"<a href='https://doi.org/{html.escape(doi)}'>{html.escape(doi)}</a>" if doi else ""
    evidence = norm_space(row.get("evidence_summary", ""))[:520]
    return f"""<tr data-status="{html.escape(row['resolution_status'])}">
<td><span class="status {html.escape(row['resolution_status'])}">{html.escape(row['resolution_status'])}</span><div class="small">{html.escape(row['resolution_reason'])}</div></td>
<td><b>{html.escape(row.get('canonical_enzyme') or row.get('specific_enzyme_names') or row.get('enzyme_name_original',''))}</b><div class="small">Names: {html.escape(row.get('specific_enzyme_names',''))}</div><div class="small">Parent: {html.escape(row.get('parent_enzyme_names',''))}</div><div class="small">Variants: {html.escape(row.get('mutations_or_variants',''))}</div><div class="small">Family: {html.escape(row.get('enzyme_family',''))}</div></td>
<td>{html.escape(row.get('organism_candidates',''))}<div class="small">{html.escape(row.get('organism_evidence',''))}</div></td>
<td class="links"><a href="{html.escape(row['uniprot_url'])}">UniProt protein</a><a href="{html.escape(row['ncbi_protein_url'])}">NCBI Protein</a><a href="{html.escape(row['ncbi_nucleotide_url'])}">NCBI Nucleotide/gene</a><a href="{html.escape(row['ena_url'])}">ENA text search</a><div class="small">Protein: {html.escape(row.get('protein_query',''))}</div><div class="small">DNA: {html.escape(row.get('dna_query',''))}</div></td>
<td class="title"><b>{html.escape(row.get('title',''))}</b><div>{doi_link}</div><div class="small">{html.escape(str(row.get('year','')))} {html.escape(row.get('journal',''))}</div><div class="evidence small">{html.escape(evidence)}</div></td>
<td>{html.escape(row.get('reaction_type',''))}<div class="small">new-to-nature={html.escape(str(row.get('new_to_nature_reaction','')))}</div><div class="small">{html.escape(row.get('characterized_enzyme_evidence',''))}</div></td>
</tr>"""


def main() -> None:
    parser = argparse.ArgumentParser(description="Build sequence seed review HTML from literature-agent CSV outputs.")
    parser.add_argument("--enzyme-csv", required=True, help="Path to enzyme_seed_candidates.csv")
    parser.add_argument("--homolog-csv", help="Optional path to homolog_candidate_seeds.csv")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--prefix", default="sequence_seed", help="Output file prefix")
    args = parser.parse_args()

    enzyme_rows = read_csv(Path(args.enzyme_csv))
    homolog_rows = read_csv(Path(args.homolog_csv)) if args.homolog_csv else []
    rows = build_candidates(combine_rows(enzyme_rows, homolog_rows))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fields = [
        "resolution_status", "resolution_reason", "specific_enzyme_names",
        "canonical_enzyme", "parent_enzyme_names", "mutations_or_variants", "organism_candidates",
        "organism_evidence", "protein_query", "dna_query", "uniprot_url",
        "ncbi_protein_url", "ncbi_nucleotide_url", "ena_url", "enzyme_family",
        "enzyme_name_original", "doi", "title", "year", "journal", "search_track",
        "flavin_cofactor", "reaction_type", "new_to_nature_reaction",
        "characterized_enzyme_evidence", "structure_similarity_hint", "criteria_status",
        "manual_review_priority", "pdb_ids", "uniprot_id", "evidence_summary",
    ]
    csv_path = outdir / f"{args.prefix}_candidates.csv"
    html_path = outdir / f"{args.prefix}_review.html"
    write_csv(csv_path, rows, fields)
    render_html(rows, html_path)
    print(f"[SequenceSeed] Wrote {csv_path} ({len(rows)} rows)")
    print(f"[SequenceSeed] Wrote {html_path}")


if __name__ == "__main__":
    main()
