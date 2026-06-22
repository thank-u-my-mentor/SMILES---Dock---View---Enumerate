#!/usr/bin/env python3
"""Build a plasmid-order candidate review for Metal-F enzyme hits."""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq


TARGET_DOCX = Path("/mnt/e/QZHAO-LAB/M-F酶.docx")
PROJECT = Path("/home/qin/Metal-F_project")
INTERACTIONS = PROJECT / "Metal-F_interactions"
DEFAULT_XLSX = INTERACTIONS / "Metal-F_element_ligand_rescan_full_merged.xlsx"
DEFAULT_FASTA = INTERACTIONS / "Metal-F_element_ligand_rescan_full_merged_protein.fasta"
DEFAULT_OUTDIR = Path("/mnt/e/QZHAO-LAB/Metal-F_plasmid_candidate_review")

sys.path.insert(0, str(PROJECT))
import enrich_metal_f_hits_with_sequences as seqmeta  # noqa: E402


MANUAL_EXPRESSION_EVIDENCE = {
    "4H4C": {
        "manual_decision_override": "AVOID_RISK",
        "literature_expression_evidence": (
            "Recombinant ispH genes were expressed in E. coli XL-1 Blue with pACYC184iscSfdx, "
            "a helper plasmid specifying iron-sulfur cluster assembly enzymes."
        ),
        "literature_expression_strain": "E. coli XL-1 Blue + pACYC184iscSfdx",
        "literature_vector_evidence": "Fe-S helper plasmid required; not a standard BL21(DE3)-only expression setup.",
        "literature_evidence_source": (
            "Wiley SI for 10.1002/anie.201208469: "
            "https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fanie.201208469&file=anie_201208469_sm_miscellaneous_information.pdf"
        ),
        "expression_curation_status": "user_curated_si_special_system",
    },
    "7U1Y": {
        "manual_decision_override": "BUY",
        "literature_expression_evidence": "Reported expression in E. coli BL21(DE3).",
        "literature_expression_strain": "E. coli BL21(DE3)",
        "literature_vector_evidence": "",
        "literature_evidence_source": "user-curated literature evidence; RCSB host also reports BL21(DE3).",
        "expression_curation_status": "user_curated_bl21",
    },
    "6RWZ": {
        "manual_decision_override": "AVOID_RISK",
        "literature_expression_evidence": "ACS SI indicates Rosetta cells and a non-simple N-terminal plasmid/vector setup.",
        "literature_expression_strain": "E. coli Rosetta",
        "literature_vector_evidence": "Likely not a simple pET22-style construct; may require pET28/pET32b-like N-terminal design.",
        "literature_evidence_source": "ACS SI for 10.1021/acscatal.0c04500; user-curated expression note.",
        "expression_curation_status": "user_curated_special_vector_or_strain",
    },
    "1T5G": {
        "manual_decision_override": "BUY",
        "literature_expression_evidence": "User-curated evidence indicates BL21(DE3)-compatible expression.",
        "literature_expression_strain": "E. coli BL21(DE3)",
        "literature_vector_evidence": "",
        "literature_evidence_source": "user-curated literature evidence for 10.1021/bi0491705.",
        "expression_curation_status": "user_curated_bl21",
    },
    "4JB4": {
        "manual_decision_override": "BUY",
        "literature_expression_evidence": "User-curated evidence indicates formal BL21(DE3) expression.",
        "literature_expression_strain": "E. coli BL21(DE3)",
        "literature_vector_evidence": "",
        "literature_evidence_source": "user-curated literature evidence for 10.1021/bi400220w.",
        "expression_curation_status": "user_curated_bl21",
    },
    "5U8Z": {
        "manual_decision_override": "BUY",
        "literature_expression_evidence": (
            "CAO1 coding sequence was synthesized, cloned into pET3a without fusion tags, "
            "and transformed into T7 Express BL21 E. coli for protein expression."
        ),
        "literature_expression_strain": "T7 Express BL21 E. coli",
        "literature_vector_evidence": "pET3a, no fusion tag.",
        "literature_evidence_source": "ACS SI for 10.1021/acs.biochem.7b00251; user-provided SI excerpt.",
        "expression_curation_status": "user_curated_bl21",
    },
    "1BS3": {
        "manual_decision_override": "AVOID_RISK",
        "literature_expression_evidence": "Native protein/purification evidence; not a recombinant BL21(DE3) expression target.",
        "literature_expression_strain": "native protein",
        "literature_vector_evidence": "No simple recombinant plasmid expression evidence for this candidate.",
        "literature_evidence_source": "user-curated native-protein triage.",
        "expression_curation_status": "user_curated_native_avoid",
    },
}


DEFAULT_EVIDENCE_FIELDS = {
    "manual_decision_override": "",
    "literature_expression_evidence": "",
    "literature_expression_strain": "",
    "literature_vector_evidence": "",
    "literature_evidence_source": "",
    "expression_curation_status": "rcsb_metadata_only",
}


def read_url_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "metal-f-plasmid-review/1.0"})
    last_error = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                return json.loads(response.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, ConnectionResetError) as exc:
            last_error = exc
            time.sleep(1.5 * (attempt + 1))
    raise last_error


def extract_docx_text(path: Path) -> str:
    with zipfile.ZipFile(path) as zf:
        xml = zf.read("word/document.xml")
    root = ET.fromstring(xml)
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    parts = []
    for text_node in root.findall(".//w:t", ns):
        if text_node.text:
            parts.append(text_node.text)
    return "\n".join(parts)


def extract_pdb_ids_from_docx(path: Path) -> list[str]:
    text = extract_docx_text(path)
    ids = re.findall(r"\b[0-9][A-Za-z0-9]{3}\b", text.upper())
    return sorted(dict.fromkeys(ids))


def parse_fasta(path: Path) -> dict[str, dict]:
    records = {}
    header = ""
    seq_parts = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    add_fasta_record(records, header, "".join(seq_parts))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header:
        add_fasta_record(records, header, "".join(seq_parts))
    return records


def parse_pdb_prefixed_fasta(path: Path) -> dict[str, list[dict]]:
    records: dict[str, list[dict]] = {}
    if not path.exists():
        return records
    header = ""
    seq_parts = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    add_pdb_prefixed_record(records, header, "".join(seq_parts))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header:
        add_pdb_prefixed_record(records, header, "".join(seq_parts))
    return records


def add_pdb_prefixed_record(records: dict[str, list[dict]], header: str, sequence: str) -> None:
    pdb_id = header.split("|", 1)[0].upper()
    if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", pdb_id):
        records.setdefault(pdb_id, []).append({"header": header, "sequence": sequence, "length": len(sequence)})


def add_fasta_record(records: dict[str, dict], header: str, sequence: str) -> None:
    pdb_entity = header.split("|", 1)[0]
    pdb_id = pdb_entity.split("_", 1)[0].upper()
    meta = {"header": header, "sequence": sequence, "pdb_id": pdb_id, "length": len(sequence)}
    for item in header.split("|")[1:]:
        if "=" in item:
            key, value = item.split("=", 1)
            meta[key] = value
    records[pdb_id] = meta


def lineage_kingdom(rows: list[dict]) -> str:
    values = []
    for row in rows or []:
        for lin in row.get("taxonomy_lineage") or []:
            if lin.get("name"):
                values.append(lin["name"].lower())
        for key in ("scientific_name", "ncbi_scientific_name", "pdbx_organism_scientific"):
            if row.get(key):
                values.append(str(row[key]).lower())
    exact = set(values)
    if "bacteria" in exact:
        return "Bacteria"
    if "archaea" in exact:
        return "Archaea"
    if "viridiplantae" in exact or "chloroplastida" in exact:
        return "Plant"
    if "fungi" in exact:
        return "Fungi"
    if "metazoa" in exact or "animalia" in exact or any(x in exact for x in ["homo sapiens", "rattus norvegicus"]):
        return "Animal"
    if "eukaryota" in exact or "eukarya" in exact:
        return "Eukaryota"
    return "Unknown"


def source_names(rows: list[dict]) -> str:
    names = []
    for row in rows or []:
        name = row.get("scientific_name") or row.get("ncbi_scientific_name") or row.get("pdbx_organism_scientific")
        if name and name not in names:
            names.append(name)
    return ";".join(names) if names else "unknown"


def host_names(entity: dict) -> str:
    names = []
    for row in entity.get("rcsb_entity_host_organism") or []:
        name = row.get("scientific_name") or row.get("ncbi_scientific_name") or row.get("pdbx_organism_scientific")
        if name and name not in names:
            names.append(name)
    for row in entity.get("entity_src_gen") or []:
        for key in ("pdbx_host_org_scientific_name", "host_org_scientific_name"):
            name = row.get(key)
            if name and name not in names:
                names.append(name)
    return ";".join(names) if names else "unknown"


def polymer_metadata(pdb_id: str) -> list[dict]:
    entry = read_url_json(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
    entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids") or []
    out = []
    for entity_id in entity_ids:
        data = read_url_json(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}")
        poly = data.get("entity_poly") or {}
        if poly.get("rcsb_entity_polymer_type") != "Protein":
            continue
        ids = data.get("rcsb_polymer_entity_container_identifiers") or {}
        src_rows = data.get("rcsb_entity_source_organism") or data.get("entity_src_nat") or []
        seq = poly.get("pdbx_seq_one_letter_code_can") or poly.get("pdbx_seq_one_letter_code") or ""
        seq = "".join(str(seq).split())
        out.append(
            {
                "entity_id": str(entity_id),
                "chains": ",".join(ids.get("auth_asym_ids") or ids.get("asym_ids") or []),
                "uniprot_ids": ";".join(ids.get("uniprot_ids") or []),
                "source_organism_corrected": source_names(src_rows),
                "source_kingdom_corrected": lineage_kingdom(src_rows),
                "expression_system_corrected": host_names(data),
                "sequence": seq,
                "length": len(seq),
            }
        )
    return out


def identity(seq_a: str, seq_b: str, aligner: PairwiseAligner) -> float:
    if not seq_a or not seq_b:
        return 0.0
    aln = aligner.align(seq_a, seq_b)[0]
    matches = 0
    aligned_len = 0
    for (a0, a1), (b0, b1) in zip(aln.aligned[0], aln.aligned[1]):
        chunk_a = seq_a[a0:a1]
        chunk_b = seq_b[b0:b1]
        for aa, bb in zip(chunk_a, chunk_b):
            aligned_len += 1
            if aa == bb:
                matches += 1
    denom = max(len(seq_a), len(seq_b))
    return matches / denom if denom else 0.0


def ecoli_expression(value: str) -> bool:
    return "escherichia coli" in str(value or "").lower()


def source_priority(value: str) -> int:
    value = str(value or "Unknown")
    if value == "Bacteria":
        return 0
    if value == "Unknown":
        return 2
    return 1


def expression_priority(value: str) -> int:
    value = str(value or "unknown")
    if ecoli_expression(value):
        return 0
    if value.lower() == "unknown":
        return 1
    return 2


def classify_risk(row: dict) -> tuple[str, str, int]:
    protein_entities = int(row.get("protein_entity_count_corrected") or 0)
    source_kingdom = str(row.get("source_kingdom_corrected") or "Unknown")
    expr = str(row.get("expression_system_corrected") or "unknown").lower()
    source = str(row.get("source_organism_corrected") or "unknown").lower()
    representative = str(row.get("cluster_representative_pdb") or row.get("pdb_id"))
    manual_override = str(row.get("manual_decision_override") or "").strip()
    curation_status = str(row.get("expression_curation_status") or "")
    lit_evidence = str(row.get("literature_expression_evidence") or "")

    reasons = []
    score = 0
    if protein_entities > 1:
        reasons.append("multi-subunit/complex protein entity")
        score += 4
    if curation_status and curation_status != "rcsb_metadata_only":
        reasons.append(f"expression evidence: {curation_status}")
    if lit_evidence:
        reasons.append(lit_evidence)
    if "escherichia coli" in expr:
        reasons.append("reported E. coli expression")
        score -= 3
    elif expr == "unknown":
        reasons.append("unknown expression system")
        score += 1
    else:
        reasons.append(f"non-E.coli expression host: {row.get('expression_system_corrected')}")
        score += 3
    if source_kingdom != "Bacteria":
        reasons.append(f"non-bacterial source: {source_kingdom}")
        score += 2
    if "homo sapiens" in source:
        reasons.append("human source, do not treat as E. coli")
        score += 3
    if str(row.get("pdb_id")) != representative:
        reasons.append(f"redundant with selected representative {representative}")
        score += 4

    if str(row.get("pdb_id")) != representative:
        decision = "AVOID_DUPLICATE"
    elif manual_override:
        decision = manual_override
        if manual_override.startswith("AVOID"):
            score += 4
        elif manual_override == "BUY":
            score -= 4
    elif protein_entities > 1:
        decision = "AVOID_COMPLEX"
    elif "homo sapiens" in source:
        decision = "AVOID_RISK"
    elif source_kingdom == "Bacteria" and "escherichia coli" in expr:
        decision = "BUY"
    elif source_kingdom == "Bacteria":
        decision = "REVIEW"
    elif source_kingdom != "Bacteria" and "escherichia coli" in expr:
        decision = "REVIEW"
    else:
        decision = "AVOID_RISK"
    return decision, "; ".join(reasons), score


def final_cluster_representatives(frame: pd.DataFrame, pairs: list[dict]) -> dict[str, str]:
    parent = {pdb_id: pdb_id for pdb_id in frame["pdb_id"]}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for _, sub in frame.groupby("uniprot_key"):
        ids = [p for p in sub["pdb_id"] if str(sub["uniprot_key"].iloc[0]) not in {"", "nan", "None"}]
        for pdb_id in ids[1:]:
            union(ids[0], pdb_id)
    for pair in pairs:
        if pair["identity"] >= 0.9:
            union(pair["pdb_id_a"], pair["pdb_id_b"])

    clusters: dict[str, list[str]] = {}
    for pdb_id in frame["pdb_id"]:
        clusters.setdefault(find(pdb_id), []).append(pdb_id)

    records = frame.set_index("pdb_id").to_dict("index")
    representatives = {}
    for members in clusters.values():
        def key(pdb_id: str):
            row = records[pdb_id]
            return (
                int(row.get("protein_entity_count_corrected") or 99),
                expression_priority(row.get("expression_system_corrected", "")),
                source_priority(row.get("source_kingdom_corrected", "")),
                float(row.get("nearest_F_metal_A") or 999),
                pdb_id,
            )

        rep = sorted(members, key=key)[0]
        for pdb_id in members:
            representatives[pdb_id] = rep
    return representatives


def json_scalar(value):
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, (int, float, bool)):
        return value
    return str(value)


def wrap_sequence(sequence: str, width: int = 80) -> str:
    sequence = "".join(str(sequence or "").split())
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def normalize_aa(sequence: str) -> str:
    return "".join(str(sequence or "").replace("*", "").split()).upper()


def translate_dna(sequence: str) -> str:
    dna = re.sub(r"[^ACGTUacgtu]", "", str(sequence or "")).upper().replace("U", "T")
    if len(dna) < 3:
        return ""
    trimmed = dna[: len(dna) - (len(dna) % 3)]
    return str(Seq(trimmed).translate(to_stop=False))


def best_protein_window_identity(translated: str, protein: str) -> tuple[float, int, int, str]:
    translated = normalize_aa(translated)
    protein = normalize_aa(protein)
    if not translated or not protein:
        return 0.0, -1, -1, ""
    if protein in translated:
        start = translated.index(protein)
        return 1.0, start, start + len(protein), "protein_exact_substring_of_dna_translation"
    if translated in protein:
        start = protein.index(translated)
        return len(translated) / len(protein), 0, len(translated), "dna_translation_substring_of_protein"
    best_identity = 0.0
    best_start = 0
    best_end = 0
    if len(translated) >= len(protein):
        window = len(protein)
        for start in range(0, len(translated) - window + 1):
            chunk = translated[start : start + window]
            ident = sum(a == b for a, b in zip(chunk, protein)) / len(protein)
            if ident > best_identity:
                best_identity = ident
                best_start = start
                best_end = start + window
        return best_identity, best_start, best_end, "best_translation_window"
    window = len(translated)
    for start in range(0, len(protein) - window + 1):
        chunk = protein[start : start + window]
        ident = sum(a == b for a, b in zip(translated, chunk)) / len(protein)
        if ident > best_identity:
            best_identity = ident
            best_start = 0
            best_end = len(translated)
    return best_identity, best_start, best_end, "best_protein_window"


def dna_qc_for_record(dna_sequence: str, protein_sequence: str) -> dict:
    dna = re.sub(r"[^ACGTUacgtu]", "", str(dna_sequence or "")).upper().replace("U", "T")
    translated = translate_dna(dna).rstrip("*")
    identity, aa_start, aa_end, match_mode = best_protein_window_identity(translated, protein_sequence)
    internal_stops = max(0, translate_dna(dna).count("*") - (1 if translate_dna(dna).endswith("*") else 0))
    clean_translation = normalize_aa(translated)
    clean_protein = normalize_aa(protein_sequence)
    exact_or_clean_partial = (
        (match_mode == "protein_exact_substring_of_dna_translation" and identity >= 0.999)
        or (
            match_mode == "dna_translation_substring_of_protein"
            and clean_translation
            and clean_translation in clean_protein
            and identity >= 0.98
        )
    )
    if len(dna) % 3 != 0 or internal_stops:
        status = "review_frame_or_stop"
    elif exact_or_clean_partial:
        status = "pass_exact_or_pdb_extra_residue"
    elif identity >= 0.95:
        status = "review_near_match"
    else:
        status = "review_mismatch"
    trimmed_dna = ""
    if status.startswith("pass") and aa_start >= 0 and aa_end > aa_start and len(dna) >= aa_end * 3:
        trimmed_dna = dna[aa_start * 3 : aa_end * 3]
    return {
        "dna_length": len(dna),
        "dna_mod3": len(dna) % 3,
        "dna_start_codon": dna[:3],
        "dna_terminal_codon": dna[-3:] if len(dna) >= 3 else "",
        "dna_internal_stop_count": internal_stops,
        "translated_aa_length": len(translated),
        "protein_length": len(normalize_aa(protein_sequence)),
        "dna_translation_identity_to_pdb_protein": identity,
        "dna_match_aa_start": aa_start,
        "dna_match_aa_end": aa_end,
        "dna_match_mode": match_mode,
        "dna_qc_status": status,
        "pdb_matched_dna": trimmed_dna if status.startswith("pass") else "",
    }


def html_report(
    frame: pd.DataFrame,
    identity_pairs: list[dict],
    doc_ids: list[str],
    dna_index_rows: list[dict] | None = None,
    dna_fasta_records: dict[str, list[dict]] | None = None,
    scene_image_rows: list[dict] | None = None,
) -> str:
    records = []
    for raw in frame.fillna("").to_dict("records"):
        record = {key: json_scalar(value) for key, value in raw.items()}
        pdb_id = str(record.get("pdb_id", "")).upper()
        sequence = str(record.get("sequence", ""))
        header = str(record.get("fasta_header", "")) or pdb_id
        record["pdb_id"] = pdb_id
        record["pdb_url"] = f"https://www.rcsb.org/structure/{pdb_id}" if pdb_id else ""
        record["protein_length"] = len("".join(sequence.split()))
        record["protein_fasta"] = f">{header}\n{wrap_sequence(sequence)}" if sequence else ""
        records.append(record)

    pairs = []
    for raw in identity_pairs or []:
        pair = {key: json_scalar(value) for key, value in raw.items()}
        try:
            pair["identity"] = float(pair.get("identity") or 0)
        except Exception:
            pair["identity"] = 0.0
        pair["identity_percent"] = f"{pair['identity'] * 100:.1f}%"
        pairs.append(pair)

    dna_by_pdb: dict[str, list[dict]] = {}
    for raw in dna_index_rows or []:
        row = {key: json_scalar(value) for key, value in raw.items()}
        pdb_id = str(row.get("pdb_id", "")).upper()
        if pdb_id:
            dna_by_pdb.setdefault(pdb_id, []).append(row)
    dna_fasta_by_pdb: dict[str, list[dict]] = {}
    for pdb_id, fasta_rows in (dna_fasta_records or {}).items():
        dna_fasta_by_pdb[pdb_id] = [
            {**row, "fasta": f">{row.get('header','')}\n{wrap_sequence(row.get('sequence',''))}"}
            for row in fasta_rows
        ]
    scene_images: dict[str, list[dict]] = {}
    for raw in scene_image_rows or []:
        row = {key: json_scalar(value) for key, value in raw.items()}
        pdb_id = str(row.get("pdb_id", "")).upper()
        if pdb_id:
            scene_images.setdefault(pdb_id, []).append(row)

    payload = {
        "doc_ids": doc_ids,
        "candidates": records,
        "identity_pairs": pairs,
        "dna_index": dna_by_pdb,
        "dna_fasta": dna_fasta_by_pdb,
        "scene_images": scene_images,
    }
    data_json = json.dumps(payload, ensure_ascii=False).replace("</", "<\\/")
    template = """<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Metal-F Plasmid Candidate Dashboard</title>
  <style>
    :root {
      --bg:#f6f6f2; --panel:#ffffff; --ink:#202326; --muted:#687076; --line:#d8ddd6;
      --blue:#2f5f7f; --green:#3f7155; --amber:#9a6b2f; --red:#9a4a43;
      --pink:#d88aa1; --orange:#c77a35; --soft:#eef1ed; --shadow:0 12px 30px rgba(33,37,41,.08);
    }
    * { box-sizing:border-box; }
    body { margin:0; font-family:Arial,"Microsoft YaHei",sans-serif; background:var(--bg); color:var(--ink); }
    header { position:sticky; top:0; z-index:5; background:#fff; border-bottom:1px solid var(--line); padding:16px 22px 14px; }
    h1 { margin:0 0 10px; font-size:24px; letter-spacing:0; }
    .controls { display:grid; grid-template-columns:1.45fr repeat(5,minmax(128px,1fr)); gap:10px; }
    input,select { width:100%; border:1px solid var(--line); border-radius:6px; background:#fff; padding:9px 10px; font-size:14px; }
    .stats { display:flex; flex-wrap:wrap; gap:8px; margin-top:10px; color:var(--muted); font-size:13px; }
    .stat { background:var(--soft); border-radius:6px; padding:5px 8px; }
    main { display:grid; grid-template-columns:minmax(360px,.92fr) minmax(460px,1.08fr); gap:16px; padding:16px 22px 24px; }
    .list { display:flex; flex-direction:column; gap:10px; max-height:calc(100vh - 150px); overflow:auto; padding-right:4px; }
    .card { background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:12px; cursor:pointer; }
    .card.active { border-color:var(--blue); box-shadow:0 0 0 2px rgba(47,95,127,.14); }
    .row { display:flex; align-items:center; justify-content:space-between; gap:10px; }
    .pdb { font-weight:700; font-size:18px; }
    .meta { color:var(--muted); font-size:12px; line-height:1.45; }
    .chips { display:flex; flex-wrap:wrap; gap:6px; margin-top:9px; }
    .chip { background:var(--soft); border-radius:999px; padding:3px 8px; font-size:12px; }
    .BUY { color:#fff; background:var(--green); }
    .REVIEW { color:#fff; background:var(--amber); }
    .AVOID_DUPLICATE,.AVOID_COMPLEX,.AVOID_RISK { color:#fff; background:var(--red); }
    .Bacteria { color:#fff; background:var(--pink); }
    .Archaea { color:#fff; background:var(--orange); }
    .Plant,.Fungi { color:#143f2b; background:#b8d9c4; }
    .detail { background:var(--panel); border:1px solid var(--line); border-radius:8px; min-height:360px; padding:16px; box-shadow:var(--shadow); }
    .detail h2 { margin:0 0 4px; font-size:22px; }
    .toolbar { display:flex; gap:8px; flex-wrap:wrap; margin:12px 0; }
    button { border:1px solid var(--line); border-radius:6px; background:#fff; padding:8px 10px; cursor:pointer; font-size:14px; }
    button.active { background:var(--blue); color:#fff; border-color:var(--blue); }
    .kv { display:grid; grid-template-columns:168px 1fr; gap:8px 12px; font-size:14px; margin-top:12px; }
    .kv div:nth-child(odd) { color:var(--muted); }
    .notes { margin-top:12px; line-height:1.55; }
    textarea { width:100%; min-height:120px; resize:vertical; border:1px solid var(--line); border-radius:6px; padding:10px; font:13px/1.45 Arial,"Microsoft YaHei",sans-serif; }
    .hint { color:var(--muted); font-size:12px; margin-top:8px; line-height:1.45; }
    pre { white-space:pre-wrap; overflow:auto; background:#f8f8f5; border:1px solid var(--line); border-radius:6px; padding:10px; font:12px/1.45 Consolas,monospace; }
    .pair { border:1px solid var(--line); border-radius:6px; padding:8px; margin:8px 0; }
    .pair.high { border-color:var(--red); background:#fff1ef; }
    .scene-img { width:100%; max-height:720px; object-fit:contain; border:1px solid var(--line); border-radius:8px; background:#fff; }
    a { color:var(--blue); }
    @media (max-width:980px) {
      main { grid-template-columns:1fr; }
      .controls { grid-template-columns:1fr 1fr; }
      .list { max-height:520px; }
      .kv { grid-template-columns:1fr; }
    }
  </style>
</head>
<body>
  <header>
    <h1>Metal-F Plasmid Candidate Dashboard</h1>
    <div class="controls">
      <input id="q" placeholder="Search PDB, UniProt, organism, ligand, metal, note">
      <select id="decision"></select>
      <select id="kingdom"></select>
      <select id="expression"></select>
      <select id="metal"></select>
      <select id="entity"></select>
    </div>
    <div class="stats" id="stats"></div>
  </header>
  <main>
    <section class="list" id="list"></section>
    <section class="detail" id="detail"></section>
  </main>
<script>
const DATA = __DATA__;
const candidates = DATA.candidates || [];
const identityPairs = DATA.identity_pairs || [];
let filtered = candidates.slice();
let selectedId = "";
let activeTab = "overview";
const $ = id => document.getElementById(id);
function esc(s){
  return String(s ?? "").replace(/[&<>"']/g,c=>{
    if(c==="&") return "&amp;";
    if(c==="<") return "&lt;";
    if(c===">") return "&gt;";
    if(c==='"') return "&quot;";
    return "&#39;";
  });
}
function num(v){const n=parseFloat(v); return Number.isFinite(n)?n:null;}
function pct(v){const n=num(v); return n===null?"":(n*100).toFixed(1)+"%";}
function dist(v){const n=num(v); return n===null?"":n.toFixed(2)+" A";}
function exprBucket(c){
  const v=String(c.expression_system_corrected||"").toLowerCase();
  if(v.includes("escherichia coli")) return "E. coli";
  if(!v || v==="unknown") return "unknown";
  return "other host";
}
function entityBucket(c){
  const n=parseInt(c.protein_entity_count_corrected||0,10);
  if(n>1) return "complex";
  if(n===1) return "single protein";
  return "unknown";
}
function uniq(field){
  return [...new Set(candidates.map(c=>c[field]).filter(Boolean))].sort();
}
function fillSelect(id,label,values){
  $(id).innerHTML = `<option value="">${esc(label)}</option>` + values.map(v=>`<option value="${esc(v)}">${esc(v)}</option>`).join("");
}
function setupFilters(){
  fillSelect("decision","All decisions",uniq("final_decision"));
  fillSelect("kingdom","All kingdoms",uniq("source_kingdom_corrected"));
  fillSelect("expression","All expression",["E. coli","unknown","other host"]);
  fillSelect("metal","All metals",uniq("metal_comp_ids"));
  fillSelect("entity","All entity types",["single protein","complex","unknown"]);
  ["q","decision","kingdom","expression","metal","entity"].forEach(id=>$(id).addEventListener("input",applyFilters));
}
function applyFilters(){
  const q=$("q").value.toLowerCase().trim();
  const decision=$("decision").value, kingdom=$("kingdom").value, expression=$("expression").value, metal=$("metal").value, entity=$("entity").value;
  filtered=candidates.filter(c=>{
    const blob=[c.pdb_id,c.uniprot_ids_primary,c.source_organism_corrected,c.expression_system_corrected,c.fluorine_comp_ids,c.metal_comp_ids,c.risk_notes,c.doi,c.nearest_F_atom,c.nearest_metal_atom].join(" ").toLowerCase();
    return (!q||blob.includes(q)) && (!decision||c.final_decision===decision) && (!kingdom||c.source_kingdom_corrected===kingdom) && (!expression||exprBucket(c)===expression) && (!metal||c.metal_comp_ids===metal) && (!entity||entityBucket(c)===entity);
  });
  if(!filtered.some(c=>c.pdb_id===selectedId)) selectedId = filtered[0]?.pdb_id || "";
  render();
}
function renderStats(){
  const total=filtered.length;
  const buy=filtered.filter(c=>c.final_decision==="BUY").length;
  const review=filtered.filter(c=>c.final_decision==="REVIEW").length;
  const duplicate=filtered.filter(c=>c.final_decision==="AVOID_DUPLICATE").length;
  const complex=filtered.filter(c=>c.final_decision==="AVOID_COMPLEX").length;
  const avoid=filtered.filter(c=>String(c.final_decision||"").startsWith("AVOID")).length;
  const nondup=filtered.filter(c=>c.final_decision!=="AVOID_DUPLICATE").length;
  $("stats").innerHTML=[`Showing ${total} / ${candidates.length}`,`BUY ${buy}`,`REVIEW ${review}`,`AVOID ${avoid}`,`complex ${complex}`,`identity>90 duplicate ${duplicate}`,`PSE nonduplicate ${nondup}`].map(x=>`<span class="stat">${esc(x)}</span>`).join("");
}
function renderList(){
  $("list").innerHTML=filtered.map(c=>{
    const active=c.pdb_id===selectedId?"active":"";
    return `<article class="card ${active}" data-pdb="${esc(c.pdb_id)}">
      <div class="row"><div class="pdb">${esc(c.pdb_id)}</div><span class="chip ${esc(c.final_decision)}">${esc(c.final_decision)}</span></div>
      <div class="meta">${esc(c.source_organism_corrected||"unknown")} · ${esc(c.uniprot_ids_primary||"no UniProt")} · ${dist(c.nearest_F_metal_A)}</div>
      <div class="chips">
        <span class="chip ${esc(c.source_kingdom_corrected)}">${esc(c.source_kingdom_corrected||"Unknown")}</span>
        <span class="chip">${esc(exprBucket(c))}</span>
        <span class="chip">${esc(c.fluorine_comp_ids||"F?")}</span>
        <span class="chip">${esc(c.metal_comp_ids||"metal?")}</span>
        <span class="chip">max id ${esc(pct(c.max_identity_to_other_selected))}</span>
      </div>
    </article>`;
  }).join("");
  document.querySelectorAll(".card").forEach(el=>el.addEventListener("click",()=>{selectedId=el.dataset.pdb; activeTab="overview"; render();}));
}
function current(){return candidates.find(c=>c.pdb_id===selectedId)||filtered[0]||candidates[0]||null;}
function pairsFor(pdb){
  return identityPairs.filter(p=>p.pdb_id_a===pdb||p.pdb_id_b===pdb).sort((a,b)=>(b.identity||0)-(a.identity||0));
}
function setTab(tab){activeTab=tab; renderDetail();}
function tabButton(id,label){return `<button class="${activeTab===id?"active":""}" onclick="setTab('${id}')">${label}</button>`;}
function overview(c){
  const doi=c.doi?`<a href="https://doi.org/${esc(c.doi)}" target="_blank">${esc(c.doi)}</a>`:"";
  return `<div class="kv">
    <div>PDB</div><div><a href="${esc(c.pdb_url)}" target="_blank">${esc(c.pdb_id)}</a></div>
    <div>Decision</div><div><span class="chip ${esc(c.final_decision)}">${esc(c.final_decision)}</span></div>
    <div>Nearest F-metal</div><div>${dist(c.nearest_F_metal_A)} · ${esc(c.nearest_F_atom)} -- ${esc(c.nearest_metal_atom)}</div>
    <div>Ligand / metal</div><div>${esc(c.fluorine_comp_ids)} / ${esc(c.metal_comp_ids)}</div>
    <div>Source organism</div><div>${esc(c.source_organism_corrected)} · ${esc(c.source_kingdom_corrected)}</div>
    <div>Expression host</div><div>${esc(c.expression_system_corrected)}</div>
    <div>Literature strain</div><div>${esc(c.literature_expression_strain || "")}</div>
    <div>Vector/helper</div><div>${esc(c.literature_vector_evidence || "")}</div>
    <div>Evidence source</div><div>${esc(c.literature_evidence_source || "")}</div>
    <div>Protein entities</div><div>${esc(c.protein_entity_count_corrected)} · ${esc(entityBucket(c))}</div>
    <div>UniProt</div><div>${esc(c.uniprot_ids_primary || c.uniprot_ids || "")}</div>
    <div>DOI</div><div>${doi}</div>
    <div>Representative</div><div>${esc(c.cluster_representative_pdb)} · partner ${esc(c.nearest_identity_partner || "")} · max identity ${esc(pct(c.max_identity_to_other_selected))}</div>
  </div><div class="notes">${esc([c.literature_expression_evidence,c.risk_notes].filter(Boolean).join("\\n"))}</div>`;
}
function notesTab(c){
  const key="metal_f_candidate_note_"+c.pdb_id;
  const value=localStorage.getItem(key)||"";
  setTimeout(()=>{
    const box=document.getElementById("userNote");
    if(box) box.addEventListener("input",()=>localStorage.setItem(key,box.value));
  },0);
  return `<h3>Browser note for ${esc(c.pdb_id)}</h3><textarea id="userNote" placeholder="Add your plasmid-ordering note here. Saved only in this browser.">${esc(value)}</textarea><div class="hint">This note is stored in browser localStorage. It does not modify the CSV/XLSX on disk.</div>`;
}
function fastaTab(c){
  const dnaRows=(DATA.dna_index||{})[c.pdb_id]||[];
  const dnaFasta=(DATA.dna_fasta||{})[c.pdb_id]||[];
  const dna=dnaRows.length?dnaRows.map(r=>`${r.dna_status||""} | qc ${r.dna_qc_status||""} | id ${r.dna_translation_identity_to_pdb_protein||""} | UniProt ${r.uniprot_id||""} | EMBL ${r.embl_accession||""} | coding ${r.coding_accession||""} | length ${r.dna_length||""}`).join("\\n"):"No CDS/DNA record was fetched for this candidate in the BUY-only DNA sidecar.";
  const dnaSeq=dnaFasta.length?dnaFasta.map(r=>r.fasta||"").join("\\n\\n"):"No CDS/DNA FASTA sequence embedded for this candidate.";
  return `<h3>Protein FASTA</h3><pre>${esc(c.protein_fasta||"No protein sequence")}</pre><h3>CDS/DNA index</h3><pre>${esc(dna)}</pre><h3>CDS/DNA FASTA</h3><pre>${esc(dnaSeq)}</pre>`;
}
function sceneTab(c){
  const rows=(DATA.scene_images||{})[c.pdb_id]||[];
  if(!rows.length) return "<p class='meta'>No rendered PyMOL scene image for this candidate.</p>";
  return rows.map(r=>`<h3>${esc(r.scene||c.pdb_id)}</h3><a href="${esc(r.png_relpath)}" target="_blank"><img class="scene-img" src="${esc(r.png_relpath)}" loading="lazy" alt="${esc(c.pdb_id)} scene"></a><div class="meta">${esc(r.fluorine_comp_ids||"")} / ${esc(r.metal_comp_ids||"")} · ${esc(r.nearest_F_metal_A||"")} A</div>`).join("");
}
function identityTab(c){
  const rows=pairsFor(c.pdb_id);
  if(!rows.length) return "<p class='meta'>No pairwise identity rows.</p>";
  return rows.map(p=>{
    const other=p.pdb_id_a===c.pdb_id?p.pdb_id_b:p.pdb_id_a;
    const high=(p.identity||0)>=0.9?"high":"";
    return `<div class="pair ${high}"><b>${esc(c.pdb_id)} vs ${esc(other)}</b><br>${esc(p.identity_percent||pct(p.identity))}</div>`;
  }).join("");
}
function renderDetail(){
  const c=current();
  if(!c){$("detail").innerHTML="<p>No candidates.</p>"; return;}
  const body=activeTab==="fasta"?fastaTab(c):activeTab==="identity"?identityTab(c):activeTab==="scene"?sceneTab(c):activeTab==="notes"?notesTab(c):overview(c);
  $("detail").innerHTML=`<h2>${esc(c.pdb_id)} <span class="chip ${esc(c.final_decision)}">${esc(c.final_decision)}</span></h2>
    <div class="meta">${esc(c.source_organism_corrected||"unknown")} · ${esc(c.expression_system_corrected||"unknown")} · ${dist(c.nearest_F_metal_A)}</div>
    <div class="toolbar">${tabButton("overview","Overview")}${tabButton("scene","Scene")}${tabButton("fasta","FASTA")}${tabButton("identity","Identity")}${tabButton("notes","Notes")}</div>
    ${body}`;
}
function render(){renderStats(); renderList(); renderDetail();}
setupFilters();
selectedId=candidates[0]?.pdb_id||"";
render();
</script>
</body>
</html>"""
    return template.replace("__DATA__", data_json)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--docx", default=str(TARGET_DOCX))
    parser.add_argument("--xlsx", default=str(DEFAULT_XLSX))
    parser.add_argument("--protein-fasta", default=str(DEFAULT_FASTA))
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR))
    args = parser.parse_args()

    doc_ids = extract_pdb_ids_from_docx(Path(args.docx))
    xlsx = Path(args.xlsx)
    geometry = pd.read_excel(xlsx, sheet_name="geometry_hits").fillna("")
    geometry["pdb_id"] = geometry["pdb_id"].astype(str).str.upper()
    selected = geometry[geometry["pdb_id"].isin(doc_ids)].copy()
    fasta = parse_fasta(Path(args.protein_fasta))

    meta_cache = {}
    rows = []
    for _, row in selected.iterrows():
        pdb_id = row["pdb_id"]
        entities = meta_cache.get(pdb_id)
        if entities is None:
            entities = polymer_metadata(pdb_id)
            meta_cache[pdb_id] = entities
        proteins = entities
        primary = proteins[0] if proteins else {}
        out = row.to_dict()
        out["protein_entity_count_corrected"] = len(proteins)
        out["source_organism_corrected"] = ";".join(dict.fromkeys(e.get("source_organism_corrected", "") for e in proteins if e.get("source_organism_corrected"))) or "unknown"
        kingdoms = [e.get("source_kingdom_corrected", "Unknown") for e in proteins if e.get("source_kingdom_corrected")]
        out["source_kingdom_corrected"] = ";".join(dict.fromkeys(kingdoms)) if kingdoms else "Unknown"
        out["expression_system_corrected"] = ";".join(dict.fromkeys(e.get("expression_system_corrected", "") for e in proteins if e.get("expression_system_corrected"))) or "unknown"
        out["uniprot_ids_primary"] = primary.get("uniprot_ids", out.get("uniprot_ids", ""))
        if pdb_id in fasta:
            out["sequence"] = fasta[pdb_id]["sequence"]
            out["fasta_header"] = fasta[pdb_id]["header"]
        elif primary.get("sequence"):
            out["sequence"] = primary["sequence"]
            out["fasta_header"] = f"{pdb_id}_{primary.get('entity_id','1')}|chains={primary.get('chains','')}|uniprot={primary.get('uniprot_ids','')}"
        else:
            out["sequence"] = ""
            out["fasta_header"] = ""
        evidence = dict(DEFAULT_EVIDENCE_FIELDS)
        evidence.update(MANUAL_EXPRESSION_EVIDENCE.get(pdb_id, {}))
        out.update(evidence)
        rows.append(out)

    review = pd.DataFrame(rows)
    if review.empty:
        raise SystemExit("No selected PDB IDs from docx were found in geometry_hits.")

    review["_dist"] = pd.to_numeric(review["nearest_F_metal_A"], errors="coerce")
    review["uniprot_key"] = review["uniprot_ids_primary"].astype(str).str.split(";").str[0]
    review = review.sort_values(["uniprot_key", "_dist", "pdb_id"], na_position="last")
    review["uniprot_duplicate_rank"] = review.groupby("uniprot_key").cumcount() + 1

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0
    ids = list(review["pdb_id"])
    seqs = dict(zip(review["pdb_id"], review["sequence"]))
    max_identity = {pdb_id: 0.0 for pdb_id in ids}
    nearest_identity_partner = {pdb_id: "" for pdb_id in ids}
    pairs = []
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            ident = identity(seqs.get(a, ""), seqs.get(b, ""), aligner)
            pairs.append({"pdb_id_a": a, "pdb_id_b": b, "identity": ident})
            if ident > max_identity[a]:
                max_identity[a] = ident
                nearest_identity_partner[a] = b
            if ident > max_identity[b]:
                max_identity[b] = ident
                nearest_identity_partner[b] = a
    review["max_identity_to_other_selected"] = review["pdb_id"].map(max_identity)
    review["nearest_identity_partner"] = review["pdb_id"].map(nearest_identity_partner)
    reps = final_cluster_representatives(review, pairs)
    review["cluster_representative_pdb"] = review["pdb_id"].map(reps)

    classified = review.apply(lambda row: classify_risk(row.to_dict()), axis=1)
    review["final_decision"] = [item[0] for item in classified]
    review["risk_notes"] = [item[1] for item in classified]
    review["risk_score"] = [item[2] for item in classified]
    decision_order = {"BUY": 0, "REVIEW": 1, "AVOID_DUPLICATE": 2, "AVOID_COMPLEX": 3, "AVOID_RISK": 4}
    review["_decision_order"] = review["final_decision"].map(decision_order).fillna(9)
    review = review.sort_values(["_decision_order", "risk_score", "_dist", "pdb_id"], ascending=[True, True, True, True])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "Metal-F_plasmid_candidate_review.csv"
    html_path = outdir / "Metal-F_plasmid_candidate_review.html"
    fasta_path = outdir / "Metal-F_plasmid_candidate_review_selected.fasta"
    dna_fasta_path = outdir / "Metal-F_plasmid_candidate_review_BUY_cds_dna.fasta"
    dna_matched_fasta_path = outdir / "Metal-F_plasmid_candidate_review_BUY_pdb_matched_cds_dna.fasta"
    dna_index_path = outdir / "Metal-F_plasmid_candidate_review_BUY_dna_index.csv"
    pairs_path = outdir / "Metal-F_plasmid_candidate_identity_pairs.csv"
    review.drop(columns=["_dist"], errors="ignore").to_csv(csv_path, index=False)
    pd.DataFrame(pairs).sort_values("identity", ascending=False).to_csv(pairs_path, index=False)
    fasta_records_for_html = {}
    with open(fasta_path, "w", encoding="utf-8") as handle:
        for _, row in review.iterrows():
            if row.get("final_decision") != "BUY":
                continue
            seq = str(row.get("sequence", ""))
            header = str(row.get("fasta_header", "")) or str(row.get("pdb_id", ""))
            if not seq:
                continue
            handle.write(f">{header}\n")
            handle.write("\n".join(seq[i : i + 80] for i in range(0, len(seq), 80)) + "\n")
            fasta_records_for_html[str(row["pdb_id"]).upper()] = {"header": header, "sequence": seq}

    dna_rows = []
    with open(dna_fasta_path, "w", encoding="utf-8") as handle:
        for _, row in review.iterrows():
            if row.get("final_decision") != "BUY":
                continue
            for uniprot_id in [x for x in str(row.get("uniprot_ids_primary", "")).split(";") if x]:
                fetched = False
                try:
                    refs = seqmeta.uniprot_embl_refs(uniprot_id)
                except Exception as exc:
                    dna_rows.append({"pdb_id": row["pdb_id"], "uniprot_id": uniprot_id, "dna_status": "uniprot_error", "error": str(exc)})
                    continue
                for ref in refs:
                    if fetched:
                        break
                    coding_rows = seqmeta.ena_coding_rows_by_protein_id(ref.get("protein_id", ""))
                    coding = coding_rows[0].get("accession", "") if coding_rows else ""
                    records = seqmeta.fetch_ena_fasta(coding or ref.get("embl_accession", ""))
                    if not records:
                        continue
                    ena_header, seq = records[0]
                    header = (
                        f"{row['pdb_id']}|uniprot={uniprot_id}|embl={ref.get('embl_accession','')}"
                        f"|coding={coding}|protein_id={ref.get('protein_id','')}"
                    )
                    handle.write(f">{header}\n")
                    handle.write("\n".join(seq[i : i + 80] for i in range(0, len(seq), 80)) + "\n")
                    qc = dna_qc_for_record(seq, str(row.get("sequence", "")))
                    dna_rows.append(
                        {
                            "pdb_id": row["pdb_id"],
                            "uniprot_id": uniprot_id,
                            "embl_accession": ref.get("embl_accession", ""),
                            "coding_accession": coding,
                            "dna_status": "fetched",
                            "ena_header": ena_header,
                            **{key: value for key, value in qc.items() if key != "pdb_matched_dna"},
                        }
                    )
                    if qc.get("pdb_matched_dna"):
                        dna_rows[-1]["pdb_matched_dna_length"] = len(qc["pdb_matched_dna"])
                    fetched = True
                if not fetched:
                    dna_rows.append({"pdb_id": row["pdb_id"], "uniprot_id": uniprot_id, "dna_status": "not_found"})
    pd.DataFrame(dna_rows).to_csv(dna_index_path, index=False)
    dna_fasta_records = parse_pdb_prefixed_fasta(dna_fasta_path)
    protein_by_pdb = {str(row.get("pdb_id", "")).upper(): str(row.get("sequence", "")) for _, row in review.iterrows()}
    with open(dna_matched_fasta_path, "w", encoding="utf-8") as handle:
        for pdb_id, dna_records in sorted(dna_fasta_records.items()):
            protein_sequence = protein_by_pdb.get(pdb_id, "")
            for record in dna_records:
                qc = dna_qc_for_record(record["sequence"], protein_sequence)
                if not qc.get("pdb_matched_dna"):
                    continue
                header = (
                    f"{record['header']}|pdb_matched=YES|aa_start={qc['dna_match_aa_start']}"
                    f"|aa_end={qc['dna_match_aa_end']}|identity={qc['dna_translation_identity_to_pdb_protein']:.3f}"
                )
                seq = qc["pdb_matched_dna"]
                handle.write(f">{header}\n")
                handle.write("\n".join(seq[i : i + 80] for i in range(0, len(seq), 80)) + "\n")
    scene_index_path = outdir / "scene_image_index.csv"
    scene_rows = []
    if scene_index_path.exists():
        with open(scene_index_path, newline="", encoding="utf-8", errors="replace") as handle:
            scene_rows = list(csv.DictReader(handle))
    html_path.write_text(html_report(review, pairs, doc_ids, dna_rows, dna_fasta_records, scene_rows), encoding="utf-8")
    print(f"docx_pdb_ids: {len(doc_ids)} {','.join(doc_ids)}")
    print(f"review_rows: {len(review)}")
    print(f"csv: {csv_path}")
    print(f"html: {html_path}")
    print(f"selected_fasta: {fasta_path}")
    print(f"buy_dna_fasta: {dna_fasta_path}")
    print(f"buy_pdb_matched_dna_fasta: {dna_matched_fasta_path}")
    print(f"buy_dna_index: {dna_index_path}")
    print(f"identity_pairs: {pairs_path}")


if __name__ == "__main__":
    main()
