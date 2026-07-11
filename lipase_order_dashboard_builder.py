#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


sys.path.insert(0, str(Path(__file__).resolve().parent))
from mf_plasmid_dashboard_tools import read_xls  # noqa: E402


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass
class LiteratureMeta:
    title: str = ""
    doi: str = ""
    url: str = ""
    enzyme: str = ""
    organism: str = ""
    uniprot: str = ""
    promiscuity: str = ""
    seed_notes: str = ""
    tags: str = ""


@dataclass
class FastaMeta:
    header: str = ""
    accession: str = ""
    protein_name: str = ""
    organism: str = ""
    gene_names: str = ""
    sequence: str = ""
    pdb_ids: list[str] = field(default_factory=list)


@dataclass
class OrderRecord:
    excel_row: int
    serial: str
    gene_name: str
    accession_raw: str
    accession_tokens: list[str]
    uniprot_ids: list[str]
    ncbi_ids: list[str]
    pdb_ids: list[str]
    five_prime_site: str
    three_prime_site: str
    optimization_host: str
    standard_vector: str
    vector_name: str
    antibiotic: str
    plasmid_prep: str
    endotoxin: str
    order_note: str
    dna_sequence: str
    protein_sequence: str
    dna_length: int
    protein_length: int
    terminal_stop: bool
    internal_stop_count: int
    incomplete_codon_bases: int
    protein_name: str = ""
    organism: str = ""
    gene_names_from_db: str = ""
    enzyme_class: str = ""
    ec: str = ""
    literature_title: str = ""
    doi: str = ""
    source_url: str = ""
    promiscuity: str = ""
    seed_notes: str = ""
    metadata_sources: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def clean_text(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "").replace("\xa0", " ")).strip()


def clean_dna(value: Any) -> str:
    return re.sub(r"[^ACGTUacgtu]", "", str(value or "")).upper().replace("U", "T")


def clean_protein(value: Any) -> str:
    return re.sub(r"[^A-Za-z*]", "", str(value or "")).upper()


def normalize_accession(value: str) -> str:
    value = clean_text(value).strip().strip(".")
    if not value:
        return ""
    value = value.replace("WP-", "WP_").replace("XP-", "XP_")
    if "_" in value and not value.startswith(("WP_", "XP_", "NP_", "YP_", "AP_", "KAF")):
        first = value.split("_", 1)[0]
        if re.match(r"^[A-Z0-9]{6,12}$", first):
            return first
    return value


def accession_keys(value: str) -> list[str]:
    base = normalize_accession(value)
    keys = []
    for key in [value, base, base.split(".", 1)[0] if "." in base else ""]:
        key = clean_text(key)
        if key and key not in keys:
            keys.append(key)
    return keys


def split_accession_tokens(text: str) -> list[str]:
    parts = [clean_text(x) for x in re.split(r"[/;,\s]+", text or "") if clean_text(x)]
    out: list[str] = []
    seen: set[str] = set()
    for part in parts:
        part = part.strip()
        if part and part not in seen:
            seen.add(part)
            out.append(part)
    return out


def is_pdb_id(token: str) -> bool:
    return bool(re.fullmatch(r"[0-9][A-Za-z0-9]{3}", token or ""))


def is_ncbi_id(token: str) -> bool:
    return bool(re.fullmatch(r"(?:WP|XP|NP|YP|AP)_\d+(?:\.\d+)?", token or "")) or bool(
        re.fullmatch(r"[A-Z]{3}\d{7}(?:\.\d+)?", token or "")
    )


def translate_dna(dna: str) -> tuple[str, bool, int, int]:
    dna = clean_dna(dna)
    protein: list[str] = []
    terminal_stop = False
    internal_stop_count = 0
    complete_len = len(dna) - (len(dna) % 3)
    for i in range(0, complete_len, 3):
        codon = dna[i : i + 3]
        aa = CODON_TABLE.get(codon, "X")
        if aa == "*":
            if i == complete_len - 3:
                terminal_stop = True
            else:
                internal_stop_count += 1
            continue
        protein.append(aa)
    return "".join(protein), terminal_stop, internal_stop_count, len(dna) % 3


def wrap_fasta(seq: str, width: int = 70) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def record_protein_fasta(record: OrderRecord) -> str:
    meta = [record.gene_name, record.accession_raw, f"order_row={record.excel_row}", f"aa={record.protein_length}"]
    if record.organism:
        meta.append(f"organism={record.organism}")
    return f">{' | '.join(meta)}\n{wrap_fasta(record.protein_sequence)}"


def record_dna_fasta(record: OrderRecord) -> str:
    meta = [record.gene_name, record.accession_raw, f"order_row={record.excel_row}", f"nt={record.dna_length}"]
    return f">{' | '.join(meta)}\n{wrap_fasta(record.dna_sequence)}"


def find_header_row(rows: list[list[Any]]) -> int:
    for idx, row in enumerate(rows):
        values = [clean_text(v).lower() for v in row]
        if "uniprot id" in values:
            return idx
    raise ValueError("Could not find order sheet header row containing 'uniprot ID'")


def row_value(row: list[Any], col: int) -> str:
    return clean_text(row[col]) if col < len(row) else ""


def parse_order_xls(path: Path) -> list[OrderRecord]:
    sheets = read_xls(path)
    if not sheets:
        raise ValueError(f"No sheets found in {path}")
    rows = next(iter(sheets.values()))
    header_idx = find_header_row(rows)
    records: list[OrderRecord] = []
    for excel_row, row in enumerate(rows[header_idx + 1 :], start=header_idx + 2):
        serial = row_value(row, 0)
        gene_name = row_value(row, 1)
        dna = clean_dna(row_value(row, 5))
        accession_raw = row_value(row, 22)
        if not gene_name or not dna or not re.fullmatch(r"[ACGT]+", dna):
            continue
        tokens = split_accession_tokens(accession_raw)
        pdb_ids = [t.upper() for t in tokens if is_pdb_id(t)]
        ncbi_ids = [normalize_accession(t) for t in tokens if is_ncbi_id(normalize_accession(t))]
        uniprot_ids = [normalize_accession(t) for t in tokens if t.upper() not in pdb_ids and normalize_accession(t) not in ncbi_ids]
        protein, terminal_stop, internal_stops, trailing = translate_dna(dna)
        warnings: list[str] = []
        if not protein:
            warnings.append("No protein translation from order DNA")
        if dna and not dna.startswith("ATG"):
            warnings.append("Coding sequence does not start with ATG")
        if internal_stops:
            warnings.append(f"{internal_stops} internal stop codon(s) after translation")
        if trailing:
            warnings.append(f"{trailing} trailing base(s) ignored during translation")
        records.append(
            OrderRecord(
                excel_row=excel_row,
                serial=serial,
                gene_name=gene_name,
                accession_raw=accession_raw,
                accession_tokens=tokens,
                uniprot_ids=uniprot_ids,
                ncbi_ids=ncbi_ids,
                pdb_ids=pdb_ids,
                five_prime_site=row_value(row, 2),
                three_prime_site=row_value(row, 8),
                optimization_host=row_value(row, 9),
                standard_vector=row_value(row, 12),
                vector_name=row_value(row, 13),
                antibiotic=row_value(row, 14),
                plasmid_prep=row_value(row, 19),
                endotoxin=row_value(row, 20),
                order_note=row_value(row, 24),
                dna_sequence=dna,
                protein_sequence=protein,
                dna_length=len(dna),
                protein_length=len(protein),
                terminal_stop=terminal_stop,
                internal_stop_count=internal_stops,
                incomplete_codon_bases=trailing,
                warnings=warnings,
            )
        )
    return records


def parse_literature_csv(path: Path | None) -> dict[str, LiteratureMeta]:
    if not path or not path.exists():
        return {}
    out: dict[str, LiteratureMeta] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            note = row.get("Abstract Note", "") or ""
            values: dict[str, str] = {}
            for line in note.splitlines():
                if "=" in line:
                    key, value = line.split("=", 1)
                    values[key.strip().lower()] = value.strip()
            acc = clean_text(values.get("uniprot", ""))
            if not acc:
                continue
            meta = LiteratureMeta(
                title=clean_text(row.get("Title", "")),
                doi=clean_text(row.get("DOI", "")),
                url=clean_text(row.get("Url", "")),
                enzyme=clean_text(values.get("enzyme", "")),
                organism=clean_text(values.get("organism", "")),
                uniprot=acc,
                promiscuity=clean_text(values.get("promiscuity", "")),
                seed_notes=clean_text(values.get("seed_notes", "")),
                tags=clean_text(row.get("Tags", "")),
            )
            for key in accession_keys(acc):
                out[key] = meta
    return out


def parse_fasta(path: Path | None) -> dict[str, FastaMeta]:
    if not path or not path.exists():
        return {}
    text = path.read_text(encoding="utf-8", errors="replace")
    metas: dict[str, FastaMeta] = {}
    header = ""
    chunks: list[str] = []

    def emit() -> None:
        nonlocal header, chunks
        if not header:
            return
        seq = clean_protein("".join(chunks))
        accession = fasta_accession(header)
        protein_name = fasta_protein_name(header)
        organism = fasta_organism(header)
        gene_names = fasta_gene_names(header)
        pdb_ids = sorted({x.upper() for x in re.findall(r"PDB:([0-9][A-Za-z0-9]{3})", header)})
        meta = FastaMeta(
            header=header,
            accession=accession,
            protein_name=protein_name,
            organism=organism,
            gene_names=gene_names,
            sequence=seq,
            pdb_ids=pdb_ids,
        )
        keys = set(accession_keys(accession))
        first = header.split()[0]
        for token in first.split("|"):
            keys.update(accession_keys(token))
        if "_" in first:
            keys.update(accession_keys(first.split("_", 1)[0]))
        for key in keys:
            if key:
                metas.setdefault(key, meta)
        header = ""
        chunks = []

    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            emit()
            header = line[1:].strip()
        else:
            chunks.append(line)
    emit()
    return metas


def fasta_accession(header: str) -> str:
    first = header.split()[0] if header else ""
    parts = first.split("|")
    if len(parts) >= 3 and parts[0] in {"sp", "tr"}:
        return parts[1]
    if parts:
        return parts[0]
    return first


def fasta_protein_name(header: str) -> str:
    first = header.split()[0] if header else ""
    rest = header[len(first) :].strip()
    rest = re.sub(r"\s+OS=.*$", "", rest).strip()
    return rest


def fasta_organism(header: str) -> str:
    match = re.search(r"\bOS=(.*?)(?:\s+OX=|\s+GN=|\s+PE=|\s+SV=|$)", header)
    return clean_text(match.group(1)) if match else ""


def fasta_gene_names(header: str) -> str:
    match = re.search(r"\bGN=(.*?)(?:\s+PE=|\s+SV=|$)", header)
    return clean_text(match.group(1)) if match else ""


def metadata_for_record(
    record: OrderRecord,
    literature: dict[str, LiteratureMeta],
    fasta_meta: dict[str, FastaMeta],
) -> None:
    keys: list[str] = []
    for token in record.accession_tokens:
        keys.extend(accession_keys(token))
    seen: set[str] = set()
    keys = [k for k in keys if not (k in seen or seen.add(k))]
    lit = next((literature[k] for k in keys if k in literature), None)
    fst = next((fasta_meta[k] for k in keys if k in fasta_meta), None)
    if lit:
        record.protein_name = lit.enzyme or record.protein_name
        record.organism = lit.organism or record.organism
        record.literature_title = lit.title
        record.doi = lit.doi
        record.source_url = lit.url
        record.promiscuity = lit.promiscuity
        record.seed_notes = lit.seed_notes
        record.metadata_sources.append("literature_seed_csv")
    if fst:
        record.protein_name = record.protein_name or fst.protein_name
        record.organism = record.organism or fst.organism
        record.gene_names_from_db = fst.gene_names
        for pdb in fst.pdb_ids:
            if pdb not in record.pdb_ids:
                record.pdb_ids.append(pdb)
        record.metadata_sources.append("project_fasta")
        if fst.sequence and record.protein_sequence and fst.sequence != record.protein_sequence:
            if len(fst.sequence) != len(record.protein_sequence):
                record.warnings.append(
                    f"Translated order protein length {record.protein_length}; project FASTA length {len(fst.sequence)}"
                )
    if not record.protein_name:
        record.protein_name = record.gene_name
    if not record.organism:
        record.organism = ""
    record.pdb_ids = sorted(set(record.pdb_ids))


def needleman_identity(a: str, b: str) -> tuple[float, int, int]:
    a = clean_protein(a).replace("*", "")
    b = clean_protein(b).replace("*", "")
    if not a or not b:
        return 0.0, 0, 0
    n, m = len(a), len(b)
    match_score = 2
    mismatch_score = -1
    gap_score = -2
    prev = [j * gap_score for j in range(m + 1)]
    trace: list[bytearray] = [bytearray([2] * (m + 1))]
    trace[0][0] = 0
    for i in range(1, n + 1):
        curr = [i * gap_score] + [0] * m
        dirs = bytearray(m + 1)
        dirs[0] = 1
        ai = a[i - 1]
        for j in range(1, m + 1):
            diag = prev[j - 1] + (match_score if ai == b[j - 1] else mismatch_score)
            up = prev[j] + gap_score
            left = curr[j - 1] + gap_score
            if diag >= up and diag >= left:
                curr[j] = diag
                dirs[j] = 0
            elif up >= left:
                curr[j] = up
                dirs[j] = 1
            else:
                curr[j] = left
                dirs[j] = 2
        trace.append(dirs)
        prev = curr
    i, j = n, m
    matches = 0
    compared = 0
    while i > 0 or j > 0:
        direction = trace[i][j] if i >= 0 and j >= 0 else 0
        if i > 0 and j > 0 and direction == 0:
            compared += 1
            if a[i - 1] == b[j - 1]:
                matches += 1
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or direction == 1):
            i -= 1
        else:
            j -= 1
    fraction = matches / compared if compared else 0.0
    return fraction, matches, compared


def build_identity(records: list[OrderRecord]) -> tuple[list[list[float]], list[dict[str, Any]]]:
    n = len(records)
    matrix = [[100.0 if i == j and records[i].protein_sequence else 0.0 for j in range(n)] for i in range(n)]
    pairs: list[dict[str, Any]] = []
    for i in range(n):
        for j in range(i + 1, n):
            frac, identical, compared = needleman_identity(records[i].protein_sequence, records[j].protein_sequence)
            pct = round(frac * 100, 4)
            matrix[i][j] = matrix[j][i] = pct
            pairs.append(
                {
                    "query_gene": records[i].gene_name,
                    "target_gene": records[j].gene_name,
                    "query_accession": records[i].accession_raw,
                    "target_accession": records[j].accession_raw,
                    "identity_percent": pct,
                    "identity_fraction": round(frac, 6),
                    "identical_sites": identical,
                    "compared_sites": compared,
                }
            )
    pairs.sort(key=lambda row: (-float(row["identity_percent"]), row["query_gene"], row["target_gene"]))
    return matrix, pairs


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(
    outdir: Path,
    title: str,
    source_xls: Path,
    records: list[OrderRecord],
    matrix: list[list[float]],
    pairs: list[dict[str, Any]],
) -> dict[str, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    protein_fasta = outdir / "lipase_order_proteins.fasta"
    dna_fasta = outdir / "lipase_order_coding_sequences.fasta"
    records_csv = outdir / "lipase_order_records.csv"
    pairwise_csv = outdir / "lipase_order_pairwise_identity.csv"
    matrix_csv = outdir / "lipase_order_identity_matrix.csv"
    data_json = outdir / "lipase_order_dashboard_data.json"
    html_path = outdir / "Lipase_plasmid_order_dashboard.html"

    protein_fasta.write_text("\n\n".join(record_protein_fasta(r) for r in records) + "\n", encoding="utf-8")
    dna_fasta.write_text("\n\n".join(record_dna_fasta(r) for r in records) + "\n", encoding="utf-8")

    record_rows = []
    for rec in records:
        row = asdict(rec)
        for key in ("accession_tokens", "uniprot_ids", "ncbi_ids", "pdb_ids", "metadata_sources", "warnings"):
            row[key] = "; ".join(row[key])
        row.pop("dna_sequence", None)
        row.pop("protein_sequence", None)
        record_rows.append(row)
    write_csv(records_csv, record_rows, list(record_rows[0].keys()) if record_rows else [])
    write_csv(pairwise_csv, pairs, list(pairs[0].keys()) if pairs else [])
    with matrix_csv.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["gene_name", *[r.gene_name for r in records]])
        for rec, row in zip(records, matrix):
            writer.writerow([rec.gene_name, *row])

    data = {
        "title": title,
        "source_xls": str(source_xls),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "records": [asdict(r) for r in records],
        "identity_matrix": matrix,
        "identity_pairs": pairs,
        "protein_fasta": protein_fasta.read_text(encoding="utf-8"),
        "dna_fasta": dna_fasta.read_text(encoding="utf-8"),
        "summary": {
            "record_count": len(records),
            "pair_count": len(pairs),
            "with_literature": sum(1 for r in records if "literature_seed_csv" in r.metadata_sources),
            "with_project_fasta": sum(1 for r in records if "project_fasta" in r.metadata_sources),
            "with_pdb_id": sum(1 for r in records if r.pdb_ids),
            "warnings": sum(len(r.warnings) for r in records),
        },
    }
    data_json.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    html_text = HTML_TEMPLATE.replace("__DATA_JSON__", json.dumps(data, ensure_ascii=False).replace("</", "<\\/"))
    html_path.write_text(html_text, encoding="utf-8")

    return {
        "html": html_path,
        "records_csv": records_csv,
        "pairwise_csv": pairwise_csv,
        "matrix_csv": matrix_csv,
        "protein_fasta": protein_fasta,
        "dna_fasta": dna_fasta,
        "data_json": data_json,
    }


def build(args: argparse.Namespace) -> dict[str, Any]:
    records = parse_order_xls(args.xls)
    literature = parse_literature_csv(args.literature_csv)
    fasta_meta = parse_fasta(args.project_fasta)
    for record in records:
        metadata_for_record(record, literature, fasta_meta)
    matrix, pairs = build_identity(records)
    outputs = write_outputs(args.outdir, args.title, args.xls, records, matrix, pairs)
    return {
        "records": len(records),
        "pairs": len(pairs),
        "with_literature": sum(1 for r in records if "literature_seed_csv" in r.metadata_sources),
        "with_project_fasta": sum(1 for r in records if "project_fasta" in r.metadata_sources),
        "with_pdb_id": sum(1 for r in records if r.pdb_ids),
        "warnings": sum(len(r.warnings) for r in records),
        "outputs": {k: str(v) for k, v in outputs.items()},
    }


def check(args: argparse.Namespace) -> dict[str, Any]:
    data_path = args.outdir / "lipase_order_dashboard_data.json"
    html_path = args.outdir / "Lipase_plasmid_order_dashboard.html"
    data = json.loads(data_path.read_text(encoding="utf-8"))
    records = data.get("records", [])
    matrix = data.get("identity_matrix", [])
    pairs = data.get("identity_pairs", [])
    problems: list[str] = []
    if len(records) != args.expected_records:
        problems.append(f"Expected {args.expected_records} records, found {len(records)}")
    if len(matrix) != len(records) or any(len(row) != len(records) for row in matrix):
        problems.append("Identity matrix dimensions do not match record count")
    expected_pairs = len(records) * (len(records) - 1) // 2
    if len(pairs) != expected_pairs:
        problems.append(f"Expected {expected_pairs} pairwise rows, found {len(pairs)}")
    for rec in records:
        if not rec.get("gene_name"):
            problems.append("A record is missing gene_name")
        if not rec.get("protein_sequence"):
            problems.append(f"{rec.get('gene_name', '?')} missing translated protein")
    html_text = html_path.read_text(encoding="utf-8")
    for marker in ["identity-heatmap", "Protein FASTA", "Coding DNA FASTA", "Order Sheet"]:
        if marker not in html_text:
            problems.append(f"HTML missing marker: {marker}")
    return {
        "ok": not problems,
        "problems": problems,
        "records": len(records),
        "pairs": len(pairs),
        "html": str(html_path),
    }


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Lipase Plasmid Order Dashboard</title>
<style>
:root {
  --bg:#f6f7f2;
  --panel:#ffffff;
  --ink:#202323;
  --muted:#65706d;
  --line:#d8ded6;
  --teal:#257b80;
  --blue:#315f9c;
  --amber:#d49b3d;
  --rust:#b85c4f;
  --green:#5f966f;
  --shadow:0 1px 2px rgba(20,30,24,.06),0 8px 24px rgba(20,30,24,.07);
}
* { box-sizing:border-box; }
body {
  margin:0;
  font-family:Arial, Helvetica, sans-serif;
  background:var(--bg);
  color:var(--ink);
  letter-spacing:0;
}
header {
  position:sticky;
  top:0;
  z-index:10;
  background:rgba(246,247,242,.97);
  border-bottom:1px solid var(--line);
  padding:14px 18px 12px;
}
h1 { margin:0 0 8px; font-size:22px; line-height:1.15; }
.sub { color:var(--muted); font-size:13px; display:flex; gap:10px; flex-wrap:wrap; }
.controls { margin-top:12px; display:flex; gap:8px; flex-wrap:wrap; align-items:center; }
input, select, button {
  font:inherit;
  border:1px solid var(--line);
  border-radius:6px;
  background:#fff;
  color:var(--ink);
}
input, select { padding:8px 10px; min-height:36px; }
input { min-width:260px; flex:1 1 320px; }
button { padding:7px 10px; min-height:34px; cursor:pointer; }
button.active { border-color:var(--teal); background:#e7f1ef; color:#174f53; }
main {
  display:grid;
  grid-template-columns:minmax(300px, 430px) minmax(0, 1fr);
  gap:14px;
  padding:14px;
  align-items:start;
}
.list { display:flex; flex-direction:column; gap:8px; }
.item {
  background:var(--panel);
  border:1px solid var(--line);
  border-radius:8px;
  padding:10px;
  cursor:pointer;
  box-shadow:0 1px 0 rgba(0,0,0,.02);
}
.item:hover { border-color:#a9b7af; }
.item.active { outline:2px solid var(--teal); }
.item-top { display:flex; gap:8px; justify-content:space-between; align-items:start; }
.gene { font-weight:700; font-size:15px; overflow-wrap:anywhere; }
.meta { color:var(--muted); font-size:12px; line-height:1.4; margin-top:4px; overflow-wrap:anywhere; }
.chips { display:flex; gap:5px; flex-wrap:wrap; margin-top:8px; }
.chip {
  display:inline-flex;
  align-items:center;
  border-radius:999px;
  padding:2px 7px;
  min-height:20px;
  background:#eef0ea;
  color:#4a534f;
  font-size:12px;
  line-height:1.2;
  white-space:nowrap;
}
.chip.green { background:#dceee1; color:#245d35; }
.chip.amber { background:#f5ead1; color:#805819; }
.chip.blue { background:#e1eafa; color:#244f88; }
.chip.rust { background:#f5dfda; color:#8b382e; }
.detail {
  min-width:0;
  background:var(--panel);
  border:1px solid var(--line);
  border-radius:8px;
  box-shadow:var(--shadow);
}
.detail-head { padding:14px 16px 12px; border-bottom:1px solid var(--line); }
.detail-head h2 { margin:0; font-size:22px; line-height:1.2; overflow-wrap:anywhere; }
.tabs {
  display:flex;
  gap:6px;
  flex-wrap:wrap;
  padding:10px 12px;
  border-bottom:1px solid var(--line);
  background:#fbfcf8;
}
.content { padding:14px 16px 18px; }
.grid {
  display:grid;
  grid-template-columns:minmax(130px, 180px) minmax(0, 1fr);
  gap:8px 12px;
  align-items:start;
}
.k { color:var(--muted); font-size:12px; text-transform:uppercase; }
.v { overflow-wrap:anywhere; white-space:pre-wrap; }
a { color:var(--blue); text-decoration:none; }
a:hover { text-decoration:underline; }
pre {
  margin:0;
  white-space:pre-wrap;
  overflow-wrap:anywhere;
  font:12px/1.45 Consolas, "Courier New", monospace;
  background:#f5f6f1;
  border:1px solid var(--line);
  border-radius:8px;
  padding:12px;
  max-height:520px;
  overflow:auto;
}
.seq-grid { display:grid; grid-template-columns:1fr; gap:12px; }
.summary-grid {
  display:grid;
  grid-template-columns:repeat(6, minmax(90px, 1fr));
  gap:8px;
  margin-bottom:12px;
}
.stat {
  border:1px solid var(--line);
  border-radius:8px;
  padding:9px;
  background:#fbfcf8;
}
.stat b { display:block; font-size:18px; margin-bottom:2px; }
.stat span { color:var(--muted); font-size:12px; }
.identity-wrap {
  overflow:auto;
  max-height:620px;
  border:1px solid var(--line);
  border-radius:8px;
  background:#fff;
}
table.identity-heatmap {
  border-collapse:separate;
  border-spacing:0;
  font-size:11px;
  min-width:max-content;
}
.identity-heatmap th, .identity-heatmap td {
  border-right:1px solid rgba(255,255,255,.55);
  border-bottom:1px solid rgba(255,255,255,.55);
  width:34px;
  min-width:34px;
  height:28px;
  text-align:center;
}
.corner {
  position:sticky;
  top:0;
  left:0;
  z-index:4;
  min-width:132px;
  background:#fff;
  color:var(--muted);
  border-right:1px solid var(--line);
  border-bottom:1px solid var(--line);
}
.col-head {
  position:sticky;
  top:0;
  z-index:3;
  height:124px;
  vertical-align:bottom;
  background:#fff;
  border-bottom:1px solid var(--line);
}
.col-head > div {
  writing-mode:vertical-rl;
  transform:rotate(180deg);
  white-space:nowrap;
  max-height:118px;
  overflow:hidden;
  padding:4px 2px;
  color:var(--muted);
  font-weight:400;
}
.row-head {
  position:sticky;
  left:0;
  z-index:2;
  min-width:132px;
  max-width:132px;
  background:#fff;
  text-align:right;
  padding:0 7px;
  border-right:1px solid var(--line);
  color:var(--muted);
  font-weight:400;
  white-space:nowrap;
  overflow:hidden;
  text-overflow:ellipsis;
}
.cell { cursor:pointer; color:transparent; }
.cell.show { color:#202323; font-size:10px; }
.cell.dark { color:#fff; }
.cell.focus { box-shadow:inset 0 0 0 2px #202323; }
.cell:hover { outline:2px solid #202323; outline-offset:-2px; }
.pair-detail {
  display:grid;
  grid-template-columns:130px 1fr;
  gap:6px 10px;
  border:1px solid var(--line);
  border-radius:8px;
  padding:10px;
  margin:0 0 10px;
  font-size:13px;
}
.pair-detail div:nth-child(odd) { color:var(--muted); }
.legend {
  display:grid;
  grid-template-columns:36px 1fr 42px;
  gap:8px;
  align-items:start;
  font-size:12px;
  color:var(--muted);
  margin-top:9px;
}
.ramp-wrap { position:relative; padding-bottom:22px; }
.ramp {
  height:12px;
  border-radius:999px;
  border:1px solid var(--line);
  background:linear-gradient(90deg,#b85c4f 0%,#d49b3d 35%,#e9dc91 55%,#75a976 78%,#257b80 100%);
}
.tick { position:absolute; top:15px; transform:translateX(-50%); white-space:nowrap; font-size:11px; }
.tick::before { content:""; display:block; width:1px; height:6px; margin:0 auto 2px; background:var(--line); }
.warn-list { display:flex; flex-direction:column; gap:6px; }
.warning { border-left:3px solid var(--rust); background:#fbf1ef; padding:7px 9px; border-radius:5px; font-size:13px; }
.empty { color:var(--muted); }
@media (max-width: 940px) {
  main { grid-template-columns:1fr; }
  .summary-grid { grid-template-columns:repeat(2, minmax(0, 1fr)); }
  .grid { grid-template-columns:1fr; }
  .k { margin-top:8px; }
}
</style>
</head>
<body>
<header>
  <h1 id="title">Lipase Plasmid Order Dashboard</h1>
  <div class="sub" id="subtitle"></div>
  <div class="controls">
    <input id="search" placeholder="Search gene, accession, organism, enzyme, note">
    <select id="metaFilter">
      <option value="">All metadata</option>
      <option value="literature_seed_csv">Has literature seed</option>
      <option value="project_fasta">Has project FASTA</option>
      <option value="pdb">Has PDB ID</option>
      <option value="warning">Has warning</option>
    </select>
    <select id="sort">
      <option value="serial">Order sheet</option>
      <option value="gene">Gene name</option>
      <option value="length">Protein length</option>
    </select>
  </div>
</header>
<main>
  <section class="list" id="list"></section>
  <section class="detail" id="detail"></section>
</main>
<script>
const DATA = __DATA_JSON__;
const records = DATA.records || [];
const matrix = DATA.identity_matrix || [];
const pairs = DATA.identity_pairs || [];
let selectedGene = records[0]?.gene_name || "";
let tab = "overview";

function esc(value) {
  return String(value ?? "").replace(/[&<>"']/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));
}
function byId(id) { return document.getElementById(id); }
function recIndex(gene) { return records.findIndex(r => r.gene_name === gene); }
function recByGene(gene) { return records[recIndex(gene)] || records[0]; }
function unique(arr) { return [...new Set((arr || []).filter(Boolean))]; }
function plainList(arr) { return unique(arr).join(", "); }
function fmtPct(v) { return Number.isFinite(+v) ? `${(+v).toFixed(2)}%` : "n/a"; }
function accessionUrl(token) {
  if (!token) return "";
  if (/^[0-9][A-Za-z0-9]{3}$/.test(token)) return `https://www.rcsb.org/structure/${token.toUpperCase()}`;
  if (/^(WP|XP|NP|YP|AP)_/.test(token) || /^[A-Z]{3}\d{7}/.test(token)) return `https://www.ncbi.nlm.nih.gov/protein/${encodeURIComponent(token)}`;
  if (/^UPI/.test(token)) return `https://www.uniprot.org/uniparc/${encodeURIComponent(token)}`;
  return `https://www.uniprot.org/uniprotkb/${encodeURIComponent(token.replace(/\..*$/, ""))}/entry`;
}
function accessionLinks(tokens) {
  const values = unique(tokens || []);
  if (!values.length) return '<span class="empty">none</span>';
  return values.map(t => `<a target="_blank" href="${esc(accessionUrl(t))}">${esc(t)}</a>`).join(", ");
}
function metadataChip(r) {
  const chips = [];
  if ((r.metadata_sources || []).includes("literature_seed_csv")) chips.push('<span class="chip green">literature</span>');
  if ((r.metadata_sources || []).includes("project_fasta")) chips.push('<span class="chip blue">FASTA</span>');
  if ((r.pdb_ids || []).length) chips.push('<span class="chip amber">PDB</span>');
  if ((r.warnings || []).length) chips.push('<span class="chip rust">warning</span>');
  return chips.join("");
}
function passes(r) {
  const q = byId("search").value.trim().toLowerCase();
  const f = byId("metaFilter").value;
  const blob = JSON.stringify(r).toLowerCase();
  if (q && !blob.includes(q)) return false;
  if (f === "pdb" && !(r.pdb_ids || []).length) return false;
  if (f === "warning" && !(r.warnings || []).length) return false;
  if (f && !["pdb","warning"].includes(f) && !(r.metadata_sources || []).includes(f)) return false;
  return true;
}
function sortedRows() {
  const s = byId("sort").value;
  const rows = records.filter(passes);
  if (s === "gene") rows.sort((a,b) => a.gene_name.localeCompare(b.gene_name));
  if (s === "length") rows.sort((a,b) => b.protein_length - a.protein_length || a.gene_name.localeCompare(b.gene_name));
  if (s === "serial") rows.sort((a,b) => Number(a.serial || 0) - Number(b.serial || 0));
  return rows;
}
function colorForIdentity(v) {
  const x = Math.max(0, Math.min(100, +v || 0));
  const stops = [
    [0,[184,92,79]],[35,[212,155,61]],[55,[233,220,145]],[78,[117,169,118]],[100,[37,123,128]]
  ];
  for (let i=1; i<stops.length; i++) {
    const [rv, rc] = stops[i], [lv, lc] = stops[i-1];
    if (x <= rv) {
      const t = (x - lv) / Math.max(1, rv - lv);
      const mix = k => Math.round(lc[k] + (rc[k] - lc[k]) * t);
      return `rgb(${mix(0)}, ${mix(1)}, ${mix(2)})`;
    }
  }
  return "rgb(37,123,128)";
}
function renderList() {
  const rows = sortedRows();
  if (!rows.some(r => r.gene_name === selectedGene) && rows[0]) selectedGene = rows[0].gene_name;
  byId("list").innerHTML = rows.map(r => `
    <article class="item ${r.gene_name === selectedGene ? "active" : ""}" data-gene="${esc(r.gene_name)}">
      <div class="item-top">
        <div>
          <div class="gene">${esc(r.gene_name)}</div>
          <div class="meta">${esc(r.protein_name || r.gene_name)}</div>
        </div>
        <span class="chip">${esc(r.serial || "")}</span>
      </div>
      <div class="meta">${esc(r.organism || "organism not locally resolved")} | ${esc(r.accession_raw || "no accession")} | ${r.protein_length} aa</div>
      <div class="chips">${metadataChip(r) || '<span class="chip">order sheet</span>'}</div>
    </article>`).join("") || '<div class="empty">No matching records.</div>';
  document.querySelectorAll(".item").forEach(el => el.addEventListener("click", () => {
    selectedGene = el.dataset.gene;
    render();
  }));
}
function field(label, value) {
  const rendered = value || '<span class="empty">blank</span>';
  return `<div class="k">${esc(label)}</div><div class="v">${rendered}</div>`;
}
function overview(r) {
  return `<div class="summary-grid">
    <div class="stat"><b>${esc(DATA.summary.record_count)}</b><span>records</span></div>
    <div class="stat"><b>${esc(DATA.summary.pair_count)}</b><span>identity pairs</span></div>
    <div class="stat"><b>${esc(DATA.summary.with_literature)}</b><span>literature seed</span></div>
    <div class="stat"><b>${esc(DATA.summary.with_project_fasta)}</b><span>project FASTA</span></div>
    <div class="stat"><b>${esc(DATA.summary.with_pdb_id)}</b><span>with PDB</span></div>
    <div class="stat"><b>${esc(DATA.summary.warnings)}</b><span>warnings</span></div>
  </div>
  <div class="grid">
    ${field("Gene name", `<b>${esc(r.gene_name)}</b>`)}
    ${field("Protein/enzyme", esc(r.protein_name || r.gene_name))}
    ${field("Organism", esc(r.organism || ""))}
    ${field("Order IDs", accessionLinks(r.accession_tokens))}
    ${field("UniProt-like IDs", accessionLinks(r.uniprot_ids))}
    ${field("NCBI IDs", accessionLinks(r.ncbi_ids))}
    ${field("PDB IDs", accessionLinks(r.pdb_ids))}
    ${field("Length", `${esc(r.protein_length)} aa / ${esc(r.dna_length)} nt`)}
    ${field("Metadata sources", esc(plainList(r.metadata_sources)))}
    ${field("Promiscuity", esc(r.promiscuity || ""))}
    ${field("Seed notes", esc(r.seed_notes || ""))}
    ${field("Literature", r.doi ? `<a target="_blank" href="https://doi.org/${esc(r.doi)}">${esc(r.doi)}</a>` : esc(r.literature_title || ""))}
  </div>
  ${(r.warnings || []).length ? `<h3>Warnings</h3><div class="warn-list">${r.warnings.map(w => `<div class="warning">${esc(w)}</div>`).join("")}</div>` : ""}`;
}
function fastaTab(r) {
  const protein = `>${r.gene_name} | ${r.accession_raw} | order_row=${r.excel_row} | aa=${r.protein_length}\n${wrap(r.protein_sequence)}`;
  const dna = `>${r.gene_name} | ${r.accession_raw} | order_row=${r.excel_row} | nt=${r.dna_length}\n${wrap(r.dna_sequence)}`;
  return `<div class="seq-grid">
    <div><h3>Protein FASTA</h3><pre>${esc(protein)}</pre></div>
    <div><h3>Coding DNA FASTA</h3><pre>${esc(dna)}</pre></div>
    <div><h3>All Protein FASTA</h3><pre>${esc(DATA.protein_fasta || "")}</pre></div>
  </div>`;
}
function wrap(seq) {
  return String(seq || "").replace(/(.{70})/g, "$1\n").trim();
}
function pairFor(a, b) {
  if (a === b) return {identity_percent:100, identical_sites:records[recIndex(a)]?.protein_length || 0, compared_sites:records[recIndex(a)]?.protein_length || 0};
  return pairs.find(p => (p.query_gene === a && p.target_gene === b) || (p.query_gene === b && p.target_gene === a));
}
function pairDetail(a, b) {
  const p = pairFor(a, b) || {};
  const ra = recByGene(a), rb = recByGene(b);
  return `<div>Pair</div><div><b>${esc(a)} vs ${esc(b)}</b></div>
    <div>Identity</div><div><span class="chip ${(+p.identity_percent || 0) >= 80 ? "green" : "amber"}">${fmtPct(p.identity_percent)}</span></div>
    <div>Compared sites</div><div>${esc(p.identical_sites ?? "")} identical / ${esc(p.compared_sites ?? "")} aligned residue pairs</div>
    <div>A accession</div><div>${esc(ra.accession_raw || "")}</div>
    <div>B accession</div><div>${esc(rb.accession_raw || "")}</div>
    <div>A organism</div><div>${esc(ra.organism || "")}</div>
    <div>B organism</div><div>${esc(rb.organism || "")}</div>`;
}
function identityTab(r) {
  const labels = records.map(x => x.gene_name);
  const selectedIndex = recIndex(r.gene_name);
  const header = `<tr><th class="corner">${labels.length} seq</th>${labels.map(g => `<th class="col-head" title="${esc(g)}"><div>${esc(g)}</div></th>`).join("")}</tr>`;
  const body = labels.map((g, i) => {
    const cells = labels.map((h, j) => {
      const v = matrix[i]?.[j] ?? 0;
      const text = i === j || v >= 70 ? Math.round(v) : "";
      const cls = ["cell", text !== "" ? "show" : "", v >= 90 ? "dark" : "", (i === selectedIndex || j === selectedIndex) ? "focus" : ""].filter(Boolean).join(" ");
      return `<td class="${cls}" style="background:${colorForIdentity(v)}" data-a="${esc(g)}" data-b="${esc(h)}" title="${esc(g)} vs ${esc(h)}: ${fmtPct(v)}">${esc(text)}</td>`;
    }).join("");
    return `<tr><th class="row-head" title="${esc(g)}">${esc(g)}</th>${cells}</tr>`;
  }).join("");
  setTimeout(() => {
    document.querySelectorAll(".cell").forEach(td => td.addEventListener("click", () => {
      byId("pairDetail").innerHTML = pairDetail(td.dataset.a, td.dataset.b);
    }));
    byId("pairDetail").innerHTML = pairDetail(r.gene_name, r.gene_name);
  }, 0);
  return `<div id="pairDetail" class="pair-detail"></div>
    <div class="identity-wrap"><table class="identity-heatmap">${header}${body}</table></div>
    <div class="legend"><span>low</span><div class="ramp-wrap"><div class="ramp"></div>
      <span class="tick" style="left:0%">0%</span><span class="tick" style="left:50%">50%</span>
      <span class="tick" style="left:75%">75%</span><span class="tick" style="left:90%">90%</span>
      <span class="tick" style="left:100%">100%</span></div><span>high</span></div>`;
}
function orderTab(r) {
  return `<div class="grid">
    ${field("Excel row", esc(r.excel_row))}
    ${field("Serial", esc(r.serial))}
    ${field("5 prime site", esc(r.five_prime_site))}
    ${field("3 prime site", esc(r.three_prime_site))}
    ${field("Optimization host", esc(r.optimization_host))}
    ${field("Standard vector", esc(r.standard_vector))}
    ${field("Vector name", esc(r.vector_name))}
    ${field("Antibiotic", esc(r.antibiotic))}
    ${field("Plasmid prep", esc(r.plasmid_prep))}
    ${field("Endotoxin", esc(r.endotoxin))}
    ${field("Order sheet note", esc(r.order_note))}
    ${field("Terminal stop", esc(r.terminal_stop ? "yes" : "no"))}
    ${field("Internal stops", esc(r.internal_stop_count))}
    ${field("Incomplete codon bases", esc(r.incomplete_codon_bases))}
  </div>
  <h3>Raw coding sequence</h3><pre>${esc(wrap(r.dna_sequence))}</pre>`;
}
function detailContent(r) {
  if (tab === "fasta") return fastaTab(r);
  if (tab === "identity") return identityTab(r);
  if (tab === "order") return orderTab(r);
  return overview(r);
}
function renderDetail() {
  const r = recByGene(selectedGene);
  byId("detail").innerHTML = `<div class="detail-head">
    <h2>${esc(r.gene_name)}</h2>
    <div class="meta">${esc(r.protein_name || r.gene_name)} | ${esc(r.accession_raw || "no accession")} | ${esc(r.protein_length)} aa</div>
    <div class="chips">${metadataChip(r) || '<span class="chip">order sheet</span>'}</div>
  </div>
  <nav class="tabs">
    ${["overview","fasta","identity","order"].map(t => `<button class="${tab === t ? "active" : ""}" data-tab="${t}">${esc(t === "fasta" ? "FASTA" : t === "identity" ? "Identity" : t === "order" ? "Order Sheet" : "Overview")}</button>`).join("")}
  </nav>
  <div class="content">${detailContent(r)}</div>`;
  document.querySelectorAll(".tabs button").forEach(btn => btn.addEventListener("click", () => {
    tab = btn.dataset.tab;
    renderDetail();
  }));
}
function render() {
  byId("title").textContent = DATA.title || "Lipase Plasmid Order Dashboard";
  byId("subtitle").innerHTML = `<span>${esc(DATA.summary.record_count)} genes from order sheet</span><span>${esc(DATA.summary.pair_count)} identity pairs</span><span>source: ${esc(DATA.source_xls || "")}</span>`;
  renderList();
  renderDetail();
}
["search","metaFilter","sort"].forEach(id => byId(id).addEventListener("input", render));
render();
</script>
</body>
</html>
"""


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a plasmid-order enzyme dashboard from a GenScript-style .xls file.")
    sub = parser.add_subparsers(dest="command", required=True)
    build_parser = sub.add_parser("build")
    build_parser.add_argument("--xls", type=Path, required=True)
    build_parser.add_argument("--outdir", type=Path, required=True)
    build_parser.add_argument("--title", default="Lipase Plasmid Order Dashboard")
    build_parser.add_argument("--literature-csv", type=Path)
    build_parser.add_argument("--project-fasta", type=Path)
    check_parser = sub.add_parser("check")
    check_parser.add_argument("--outdir", type=Path, required=True)
    check_parser.add_argument("--expected-records", type=int, default=29)
    args = parser.parse_args()
    if args.command == "build":
        result = build(args)
    else:
        result = check(args)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
