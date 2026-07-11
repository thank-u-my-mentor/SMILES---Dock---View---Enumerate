#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


DEFAULT_LOCAL_CSV = Path("E:/QZHAO-LAB/Metal-F_plasmid_candidate_review/Natrue-Fluorine-Enzyme.csv")
DEFAULT_LOCAL_FASTA = Path("E:/Codex/nature_fluorine_enzyme_hmmer/corrected_aligned_identity/verified_26_all.fasta")
USER_AGENT = "Codex enzyme-library-dashboard/1.0"


@dataclass
class Record:
    input_id: str
    accession: str = ""
    accession_type: str = ""
    source: str = ""
    source_url: str = ""
    protein_name: str = ""
    organism: str = ""
    gene_name: str = ""
    length: int = 0
    sequence: str = ""
    fasta_header: str = ""
    level: str = ""
    nearest_seed: str = ""
    nearest_seed_identity: str = ""
    motif_score: float | None = None
    motif_positions: dict[str, int] = field(default_factory=dict)
    pdb_ids: list[str] = field(default_factory=list)
    nearest_pdb_accession: str = ""
    nearest_pdb_identity: float | None = None
    nearest_pdb_ids: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    extra: dict[str, str] = field(default_factory=dict)


def clean_text(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip()


def clean_seq(value: Any) -> str:
    return re.sub(r"[^A-Za-z]", "", str(value or "")).upper()


def split_ids(text: str) -> list[str]:
    values = [x.strip() for x in re.split(r"[;,\s]+", text or "") if x.strip()]
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out


def classify_accession(acc: str) -> str:
    if acc.startswith(("WP_", "XP_", "NP_", "YP_", "AP_")) or re.match(r"^[A-Z]{2,4}_?\d+(?:\.\d+)?$", acc):
        return "ncbi_protein"
    return "uniprot"


def request_text(url: str, timeout: int = 30, attempts: int = 3) -> str:
    last_error: Exception | None = None
    for attempt in range(attempts):
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError:
            raise
        except Exception as exc:
            last_error = exc
            if attempt + 1 >= attempts:
                break
            time.sleep(0.8 * (attempt + 1))
    assert last_error is not None
    raise last_error


def request_json(url: str, timeout: int = 30) -> dict[str, Any]:
    return json.loads(request_text(url, timeout=timeout))


def parse_bracket_organism(text: str) -> tuple[str, str]:
    value = clean_text(text)
    match = re.search(r"\s+\[([^\[\]]+)\]\s*$", value)
    if not match:
        return value, ""
    return value[: match.start()].strip(), match.group(1).strip()


def parse_ncbi_fasta(text: str, requested_id: str) -> Record | None:
    header = ""
    parts: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">") and not header:
            header = line[1:]
        elif line.startswith(">"):
            break
        else:
            parts.append(line)
    seq = clean_seq("".join(parts))
    if not header or not seq:
        return None
    tokens = header.split(maxsplit=1)
    acc = tokens[0] if tokens else requested_id
    desc = tokens[1] if len(tokens) > 1 else ""
    protein, organism = parse_bracket_organism(desc)
    return Record(
        input_id=requested_id,
        accession=acc,
        accession_type="ncbi_protein",
        source="ncbi_efetch",
        source_url=ncbi_fasta_url(requested_id),
        protein_name=protein,
        organism=organism,
        length=len(seq),
        sequence=seq,
        fasta_header=header,
    )


def ncbi_fasta_url(acc: str) -> str:
    query = urllib.parse.urlencode({"db": "protein", "id": acc, "rettype": "fasta", "retmode": "text"})
    return f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{query}"


def uniprot_json_url(acc: str) -> str:
    return f"https://rest.uniprot.org/uniprotkb/{urllib.parse.quote(acc)}.json"


def uniprot_full_name(desc: dict[str, Any]) -> str:
    for key in ("recommendedName", "submissionNames"):
        value = desc.get(key)
        if isinstance(value, dict):
            full = value.get("fullName", {})
            if isinstance(full, dict) and full.get("value"):
                return clean_text(full.get("value"))
        if isinstance(value, list):
            for item in value:
                full = item.get("fullName", {}) if isinstance(item, dict) else {}
                if isinstance(full, dict) and full.get("value"):
                    return clean_text(full.get("value"))
    for item in desc.get("alternativeNames", []) or []:
        full = item.get("fullName", {}) if isinstance(item, dict) else {}
        if isinstance(full, dict) and full.get("value"):
            return clean_text(full.get("value"))
    return ""


def parse_uniprot_record(data: dict[str, Any], requested_id: str) -> Record:
    acc = clean_text(data.get("primaryAccession")) or requested_id
    seq = clean_seq((data.get("sequence") or {}).get("value"))
    genes = data.get("genes") or []
    gene_name = ""
    for gene in genes:
        if not isinstance(gene, dict):
            continue
        for key in ("geneName", "orderedLocusNames", "orfNames", "synonyms"):
            value = gene.get(key)
            if isinstance(value, dict) and value.get("value"):
                gene_name = clean_text(value.get("value"))
                break
            if isinstance(value, list) and value:
                first = value[0]
                if isinstance(first, dict) and first.get("value"):
                    gene_name = clean_text(first.get("value"))
                    break
        if gene_name:
            break
    pdb_ids = sorted(
        {
            clean_text(x.get("id")).upper()
            for x in data.get("uniProtKBCrossReferences", []) or []
            if isinstance(x, dict) and x.get("database") == "PDB" and x.get("id")
        }
    )
    organism = clean_text((data.get("organism") or {}).get("scientificName"))
    protein = uniprot_full_name(data.get("proteinDescription") or {})
    return Record(
        input_id=requested_id,
        accession=acc,
        accession_type="uniprot",
        source="uniprot_rest",
        source_url=uniprot_json_url(requested_id),
        protein_name=protein,
        organism=organism,
        gene_name=gene_name,
        length=len(seq),
        sequence=seq,
        fasta_header=clean_text(data.get("uniProtkbId")) or acc,
        pdb_ids=pdb_ids,
        extra={"entry_type": clean_text(data.get("entryType"))},
    )


def fetch_uniprot(acc: str) -> Record | None:
    try:
        return parse_uniprot_record(request_json(uniprot_json_url(acc)), acc)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return None
        raise


def fetch_ncbi(acc: str) -> Record | None:
    text = request_text(ncbi_fasta_url(acc))
    if "Error" in text and not text.lstrip().startswith(">"):
        return None
    return parse_ncbi_fasta(text, acc)


def read_local_csv(path: Path) -> dict[str, Record]:
    records: dict[str, Record] = {}
    if not path.exists():
        return records
    with path.open("r", newline="", encoding="utf-8-sig", errors="replace") as handle:
        for row in csv.DictReader(handle):
            acc = clean_text(row.get("accession"))
            seq = clean_seq(row.get("sequence"))
            if not acc or not seq:
                continue
            protein = clean_text(row.get("protein_name"))
            organism = clean_text(row.get("organism"))
            if not organism:
                protein, organism = parse_bracket_organism(protein)
            motif_positions = {}
            raw_positions = clean_text(row.get("motif_positions"))
            if raw_positions:
                try:
                    parsed = json.loads(raw_positions)
                    motif_positions = {clean_text(k): int(v) for k, v in parsed.items()}
                except Exception:
                    motif_positions = {}
            records[acc] = Record(
                input_id=acc,
                accession=acc,
                accession_type=classify_accession(acc),
                source="local_csv",
                source_url=str(path),
                protein_name=protein,
                organism=organism,
                length=int(float(row.get("length") or len(seq))),
                sequence=seq,
                fasta_header=acc,
                level=clean_text(row.get("level")),
                nearest_seed=clean_text(row.get("nearest_seed")),
                nearest_seed_identity=clean_text(row.get("nearest_seed_identity")),
                motif_score=parse_float(row.get("motif_active_site_score")),
                motif_positions=motif_positions,
                extra={
                    "rank": clean_text(row.get("rank")),
                    "final_score": clean_text(row.get("final_score")),
                    "domain_class": clean_text(row.get("domain_class")),
                    "hmmsearch_evalue": clean_text(row.get("hmmsearch_evalue")),
                    "hmmsearch_bitscore": clean_text(row.get("hmmsearch_bitscore")),
                },
            )
    return records


def parse_float(value: Any) -> float | None:
    try:
        if value == "" or value is None:
            return None
        return float(value)
    except Exception:
        return None


def clean_accession_from_header(header: str) -> str:
    first = str(header or "").split()[0]
    if "|" in first:
        for part in first.split("|"):
            if re.match(r"^[A-Za-z0-9_]+(?:\.\d+)?(?:-\d+)?$", part):
                return part
    return first


def parse_header_kv(header: str) -> dict[str, str]:
    meta: dict[str, str] = {}
    for part in str(header or "").split("|")[1:]:
        if "=" in part:
            key, value = part.split("=", 1)
            meta[clean_text(key)] = clean_text(value)
    return meta


def read_local_fasta(path: Path) -> dict[str, Record]:
    records: dict[str, Record] = {}
    if not path.exists():
        return records
    header = ""
    parts: list[str] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                add_fasta_record(records, header, "".join(parts), path)
            header = line[1:]
            parts = []
        else:
            parts.append(line)
    if header:
        add_fasta_record(records, header, "".join(parts), path)
    return records


def add_fasta_record(records: dict[str, Record], header: str, sequence: str, path: Path) -> None:
    acc = clean_accession_from_header(header)
    seq = clean_seq(sequence)
    if not acc or not seq:
        return
    meta = parse_header_kv(header)
    records[acc] = Record(
        input_id=acc,
        accession=acc,
        accession_type=classify_accession(acc),
        source="local_fasta",
        source_url=str(path),
        protein_name=meta.get("name", ""),
        organism=meta.get("organism", ""),
        length=len(seq),
        sequence=seq,
        fasta_header=header,
        level=meta.get("level", ""),
        extra={"rank": meta.get("rank", "")},
    )


def merge_record(base: Record, incoming: Record) -> Record:
    for attr in (
        "accession",
        "accession_type",
        "source",
        "source_url",
        "protein_name",
        "organism",
        "gene_name",
        "length",
        "sequence",
        "fasta_header",
        "level",
        "nearest_seed",
        "nearest_seed_identity",
    ):
        value = getattr(incoming, attr)
        if value and not getattr(base, attr):
            setattr(base, attr, value)
    if incoming.motif_score is not None and base.motif_score is None:
        base.motif_score = incoming.motif_score
    if incoming.motif_positions and not base.motif_positions:
        base.motif_positions = incoming.motif_positions
    if incoming.pdb_ids:
        base.pdb_ids = sorted(set(base.pdb_ids) | set(incoming.pdb_ids))
    base.extra.update({k: v for k, v in incoming.extra.items() if v and not base.extra.get(k)})
    return base


def resolve_records(ids: list[str], local_csv: Path, local_fasta: Path, fetch: bool, sleep_s: float) -> list[Record]:
    csv_records = read_local_csv(local_csv)
    fasta_records = read_local_fasta(local_fasta)
    out: list[Record] = []
    for acc in ids:
        record = Record(input_id=acc, accession=acc, accession_type=classify_accession(acc))
        if acc in csv_records:
            record = merge_record(record, csv_records[acc])
        if acc in fasta_records:
            record = merge_record(record, fasta_records[acc])
        if fetch and classify_accession(acc) == "uniprot":
            try:
                fetched = fetch_uniprot(acc)
            except Exception as exc:
                fetched = None
                record.warnings.append(f"UniProt fetch failed: {type(exc).__name__}")
            if fetched:
                fetched.input_id = acc
                record = merge_record(record, fetched)
                # Prefer current UniProt metadata/PDB over older local headers.
                record.protein_name = fetched.protein_name or record.protein_name
                record.organism = fetched.organism or record.organism
                record.gene_name = fetched.gene_name or record.gene_name
                record.pdb_ids = fetched.pdb_ids or record.pdb_ids
                record.source = "local+uniprot" if record.sequence and record.source != "uniprot_rest" else fetched.source
                record.source_url = fetched.source_url
            elif not record.sequence:
                record.warnings.append("UniProt lookup returned no record")
            time.sleep(sleep_s)
        elif fetch and classify_accession(acc) == "ncbi_protein" and not record.sequence:
            try:
                fetched = fetch_ncbi(acc)
            except Exception as exc:
                fetched = None
                record.warnings.append(f"NCBI fetch failed: {type(exc).__name__}")
            if fetched:
                record = merge_record(record, fetched)
            else:
                record.warnings.append("NCBI efetch returned no FASTA")
            time.sleep(sleep_s)
        if record.sequence:
            record.length = len(record.sequence)
        else:
            record.warnings.append("No sequence found")
        if not record.protein_name:
            record.protein_name = record.accession
        out.append(record)
    return out


def motif_positions(sequence: str, motifs: list[str]) -> dict[str, int]:
    seq = clean_seq(sequence)
    return {motif: seq.find(motif) for motif in motifs}


def motif_score(positions: dict[str, int]) -> float | None:
    if not positions:
        return None
    return round(100.0 * sum(1 for value in positions.values() if value >= 0) / len(positions), 2)


def apply_motif_scores(records: list[Record], motifs: list[str]) -> None:
    if not motifs:
        return
    for rec in records:
        if not rec.sequence:
            continue
        rec.motif_positions = motif_positions(rec.sequence, motifs)
        rec.motif_score = motif_score(rec.motif_positions)
        if rec.motif_score == 0 and "No configured motifs found" not in rec.warnings:
            rec.warnings.append("No configured motifs found")


def needleman_identity(a: str, b: str) -> tuple[float, int, int]:
    a = clean_seq(a)
    b = clean_seq(b)
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


def build_identity(records: list[Record]) -> tuple[list[list[float]], list[dict[str, Any]]]:
    n = len(records)
    matrix = [[100.0 if i == j and records[i].sequence else 0.0 for j in range(n)] for i in range(n)]
    pairs: list[dict[str, Any]] = []
    for i in range(n):
        for j in range(i + 1, n):
            frac, identical, compared = needleman_identity(records[i].sequence, records[j].sequence)
            pct = round(frac * 100, 4)
            matrix[i][j] = matrix[j][i] = pct
            pairs.append(
                {
                    "query": records[i].accession,
                    "target": records[j].accession,
                    "identity_percent": pct,
                    "identity_fraction": round(frac, 6),
                    "identical_sites": identical,
                    "compared_sites": compared,
                    "query_header": records[i].fasta_header,
                    "target_header": records[j].fasta_header,
                }
            )
    pairs.sort(key=lambda row: (-float(row["identity_percent"]), row["query"], row["target"]))
    return matrix, pairs


def assign_nearest_pdb(records: list[Record], matrix: list[list[float]]) -> None:
    structural = [idx for idx, rec in enumerate(records) if rec.pdb_ids]
    for i, rec in enumerate(records):
        if rec.pdb_ids:
            rec.nearest_pdb_accession = rec.accession
            rec.nearest_pdb_identity = 100.0
            rec.nearest_pdb_ids = rec.pdb_ids
            continue
        if not structural:
            continue
        best = max(structural, key=lambda idx: matrix[i][idx])
        rec.nearest_pdb_accession = records[best].accession
        rec.nearest_pdb_identity = matrix[i][best]
        rec.nearest_pdb_ids = records[best].pdb_ids


def safe_filename_label(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "record"


def write_fasta(records: list[Record], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for rec in records:
            if not rec.sequence:
                continue
            label = rec.gene_name or rec.accession
            org = rec.organism or "unknown"
            name = rec.protein_name.replace("|", ";")
            handle.write(f">{rec.accession}|label={label}|organism={org}|name={name}\n")
            for i in range(0, len(rec.sequence), 80):
                handle.write(rec.sequence[i : i + 80] + "\n")


def record_to_row(rec: Record) -> dict[str, Any]:
    row = {
        "input_id": rec.input_id,
        "accession": rec.accession,
        "accession_type": rec.accession_type,
        "label": rec.gene_name or rec.accession,
        "gene_name": rec.gene_name,
        "protein_name": rec.protein_name,
        "organism": rec.organism,
        "length": rec.length,
        "source": rec.source,
        "source_url": rec.source_url,
        "level": rec.level,
        "nearest_seed": rec.nearest_seed,
        "nearest_seed_identity": rec.nearest_seed_identity,
        "motif_score": "" if rec.motif_score is None else rec.motif_score,
        "motif_positions": json.dumps(rec.motif_positions, sort_keys=True),
        "pdb_ids": ";".join(rec.pdb_ids),
        "nearest_pdb_accession": rec.nearest_pdb_accession,
        "nearest_pdb_identity": "" if rec.nearest_pdb_identity is None else rec.nearest_pdb_identity,
        "nearest_pdb_ids": ";".join(rec.nearest_pdb_ids),
        "warnings": "; ".join(rec.warnings),
        "md5": hashlib.md5(rec.sequence.encode()).hexdigest() if rec.sequence else "",
        "sequence": rec.sequence,
    }
    row.update({f"extra_{k}": v for k, v in sorted(rec.extra.items())})
    return row


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def clustered_order(matrix: list[list[float]]) -> list[int]:
    n = len(matrix)
    if n <= 2:
        return list(range(n))
    averages = [
        (sum(matrix[i][j] for j in range(n) if j != i) / max(1, n - 1), i)
        for i in range(n)
    ]
    start = max(averages)[1]
    order = [start]
    remaining = set(range(n))
    remaining.remove(start)
    while remaining:
        last = order[-1]
        nxt = max(remaining, key=lambda idx: (matrix[last][idx], -idx))
        remaining.remove(nxt)
        order.append(nxt)
    return order


def reorder_matrix(matrix: list[list[float]], order: list[int]) -> list[list[float]]:
    return [[matrix[i][j] for j in order] for i in order]


def record_to_meta(rec: Record) -> dict[str, Any]:
    return {
        "input_id": rec.input_id,
        "accession": rec.accession,
        "label": rec.gene_name or rec.accession,
        "accession_type": rec.accession_type,
        "protein_name": rec.protein_name,
        "organism": rec.organism,
        "gene_name": rec.gene_name,
        "length": rec.length,
        "source": rec.source,
        "source_url": rec.source_url,
        "level": rec.level,
        "nearest_seed": rec.nearest_seed,
        "nearest_seed_identity": rec.nearest_seed_identity,
        "motif_score": rec.motif_score,
        "motif_positions": rec.motif_positions,
        "pdb_ids": rec.pdb_ids,
        "nearest_pdb_accession": rec.nearest_pdb_accession,
        "nearest_pdb_identity": rec.nearest_pdb_identity,
        "nearest_pdb_ids": rec.nearest_pdb_ids,
        "warnings": rec.warnings,
        "fasta_header": rec.fasta_header,
        "sequence": rec.sequence,
    }


def build_html_data(records: list[Record], matrix: list[list[float]], pairs: list[dict[str, Any]], title: str) -> dict[str, Any]:
    c_order = clustered_order(matrix)
    original_meta = [record_to_meta(rec) for rec in records]
    clustered_records = [records[i] for i in c_order]
    top_similar = pairs[: min(30, len(pairs))]
    top_distant = sorted(pairs, key=lambda row: (float(row["identity_percent"]), row["query"], row["target"]))[: min(30, len(pairs))]
    return {
        "title": title,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "original": {
            "ids": [rec.accession for rec in records],
            "meta": original_meta,
            "matrix": matrix,
        },
        "clustered": {
            "ids": [rec.accession for rec in clustered_records],
            "meta": [record_to_meta(rec) for rec in clustered_records],
            "matrix": reorder_matrix(matrix, c_order),
        },
        "pairs": pairs,
        "topSimilar": top_similar,
        "topDistant": top_distant,
        "summary": {
            "records": len(records),
            "resolved_sequences": sum(1 for rec in records if rec.sequence),
            "direct_pdb_records": sum(1 for rec in records if rec.pdb_ids),
            "records_with_pdb_template": sum(1 for rec in records if rec.nearest_pdb_ids),
            "warnings": sum(len(rec.warnings) for rec in records),
        },
    }


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>__TITLE__</title>
  <style>
    :root {
      --bg:#f7f8f4; --panel:#fff; --ink:#202326; --muted:#687076; --line:#d9dfd7;
      --blue:#2f5f7f; --green:#3f7155; --amber:#9a6b2f; --red:#9a4a43; --soft:#eef2ec;
      --shadow:0 10px 24px rgba(33,37,41,.08);
    }
    * { box-sizing:border-box; }
    body { margin:0; background:var(--bg); color:var(--ink); font-family:Arial,"Microsoft YaHei",sans-serif; }
    header { position:sticky; top:0; z-index:20; background:#fff; border-bottom:1px solid var(--line); padding:12px 16px; }
    h1 { margin:0 0 8px; font-size:23px; letter-spacing:0; }
    button, input, select { font:inherit; }
    .topbar { display:grid; grid-template-columns:minmax(220px,1fr) auto; gap:10px; align-items:center; }
    .tabs { display:flex; gap:6px; flex-wrap:wrap; }
    .tab { border:1px solid var(--line); background:#fff; color:var(--ink); border-radius:6px; padding:7px 10px; cursor:pointer; }
    .tab.active { background:var(--ink); color:#fff; border-color:var(--ink); }
    .controls { display:grid; grid-template-columns:minmax(180px,1fr) repeat(4,minmax(112px,150px)); gap:8px; align-items:center; margin-top:8px; }
    input,select { width:100%; border:1px solid var(--line); border-radius:6px; background:#fff; padding:7px 9px; font-size:14px; }
    .stats { display:flex; flex-wrap:wrap; gap:7px; margin-top:8px; color:var(--muted); font-size:13px; }
    .stat { background:var(--soft); border-radius:6px; padding:4px 8px; }
    main { padding:14px 16px 22px; }
    .view { display:none; }
    .view.active { display:block; }
    .overview { display:grid; grid-template-columns:1fr; gap:12px; }
    .table-wrap { border:1px solid var(--line); background:#fff; border-radius:8px; overflow:auto; box-shadow:var(--shadow); }
    table { border-collapse:separate; border-spacing:0; width:100%; }
    th, td { border-bottom:1px solid var(--line); padding:8px 9px; text-align:left; vertical-align:top; font-size:13px; }
    th { position:sticky; top:0; background:#fff; z-index:2; color:var(--muted); font-weight:600; }
    tr:last-child td { border-bottom:0; }
    .muted { color:var(--muted); }
    .mono { font-family:Consolas,Menlo,monospace; }
    .badge { display:inline-block; border-radius:999px; background:var(--soft); padding:2px 7px; margin:1px 2px 1px 0; font-size:12px; white-space:nowrap; }
    .badge.warn { background:#faeadf; color:#7a3a22; }
    .badge.pdb { background:#e8f0f5; color:#214b65; }
    .matrix-layout { display:grid; grid-template-columns:minmax(0,1fr) 390px; gap:12px; }
    .heatmap-wrap { border:1px solid var(--line); background:#fff; border-radius:8px; overflow:auto; max-height:calc(100vh - 170px); box-shadow:var(--shadow); }
    table.heatmap { border-collapse:separate; border-spacing:0; width:auto; font-size:11px; }
    .heatmap th,.heatmap td { border-right:1px solid rgba(255,255,255,.38); border-bottom:1px solid rgba(255,255,255,.38); min-width:28px; height:25px; text-align:center; padding:0; }
    .heatmap th { background:#fff; color:var(--muted); font-weight:400; }
    .corner { position:sticky; top:0; left:0; z-index:5; min-width:128px; border-right:1px solid var(--line); border-bottom:1px solid var(--line); }
    .col-head { position:sticky; top:0; z-index:4; height:128px; min-width:28px; vertical-align:bottom; border-bottom:1px solid var(--line); }
    .col-head > div { writing-mode:vertical-rl; transform:rotate(180deg); white-space:nowrap; max-height:120px; overflow:hidden; padding:4px 2px; }
    .row-head { position:sticky; left:0; z-index:3; min-width:128px; max-width:128px; background:#fff; text-align:right; padding:0 7px; border-right:1px solid var(--line); white-space:nowrap; overflow:hidden; text-overflow:ellipsis; }
    .cell { cursor:pointer; color:transparent; }
    .cell:hover { outline:2px solid #202326; outline-offset:-2px; }
    .diagonal { box-shadow:inset 0 0 0 1px rgba(0,0,0,.2); }
    .side { border:1px solid var(--line); background:#fff; border-radius:8px; padding:12px; display:flex; flex-direction:column; gap:12px; max-height:calc(100vh - 170px); overflow:auto; box-shadow:var(--shadow); }
    h2 { margin:0 0 8px; font-size:16px; }
    .kv { display:grid; grid-template-columns:120px 1fr; gap:6px 9px; font-size:13px; line-height:1.35; }
    .kv div:nth-child(odd) { color:var(--muted); }
    .pair-list { display:flex; flex-direction:column; gap:6px; }
    .pair { border:1px solid var(--line); border-radius:6px; padding:7px; cursor:pointer; background:#fff; }
    .pair:hover { border-color:var(--blue); }
    .bar { height:7px; border-radius:999px; margin-top:5px; background:var(--soft); overflow:hidden; }
    .bar span { display:block; height:100%; }
    .legend { display:grid; grid-template-columns:42px 1fr 42px; gap:8px; align-items:start; font-size:12px; color:var(--muted); }
    .ramp-wrap { position:relative; padding-bottom:24px; }
    .ramp { height:12px; width:100%; border-radius:999px; background:linear-gradient(90deg,#C94C4C 0%,#ECA76A 25%,#F2E7A6 50%,#9BCB9C 75%,#63B6C2 90%,#3F77B5 100%); border:1px solid var(--line); }
    .tick { position:absolute; top:16px; transform:translateX(-50%); white-space:nowrap; font-size:11px; color:var(--muted); }
    .tick::before { content:""; display:block; width:1px; height:6px; margin:0 auto 2px; background:var(--line); }
    a { color:var(--blue); text-decoration:none; }
    a:hover { text-decoration:underline; }
    pre { margin:0; white-space:pre-wrap; word-break:break-word; background:#fff; border:1px solid var(--line); border-radius:8px; padding:12px; font-size:12px; line-height:1.45; box-shadow:var(--shadow); }
    @media (max-width:1100px) {
      .matrix-layout { grid-template-columns:1fr; }
      .side, .heatmap-wrap { max-height:none; }
      .controls { grid-template-columns:1fr 1fr; }
      .topbar { grid-template-columns:1fr; }
    }
  </style>
</head>
<body>
<header>
  <div class="topbar">
    <h1>__TITLE__</h1>
    <nav class="tabs">
      <button class="tab active" data-view="overview">Overview</button>
      <button class="tab" data-view="identity">Identity Matrix</button>
      <button class="tab" data-view="fasta">FASTA</button>
    </nav>
  </div>
  <div class="controls">
    <input id="q" placeholder="Search accession, organism, protein, PDB">
    <select id="order"><option value="clustered">Clustered order</option><option value="original">Input order</option></select>
    <select id="label"><option value="label">Gene/accession labels</option><option value="organism">Organism labels</option><option value="protein_name">Protein labels</option></select>
    <select id="threshold"><option value="0">Show all cells</option><option value="70">Highlight >=70%</option><option value="80">Highlight >=80%</option><option value="90">Highlight >=90%</option><option value="95">Highlight >=95%</option></select>
    <select id="pairMode"><option value="similar">Top similar pairs</option><option value="distant">Most distant pairs</option></select>
  </div>
  <div class="stats" id="stats"></div>
</header>
<main>
  <section id="overview" class="view active overview">
    <div class="table-wrap"><table id="recordTable"></table></div>
  </section>
  <section id="identity" class="view">
    <div class="matrix-layout">
      <div class="heatmap-wrap"><table class="heatmap" id="heatmap"></table></div>
      <aside class="side">
        <section><h2>Selected Pair</h2><div id="detail" class="kv"></div></section>
        <section>
          <h2>Color Scale</h2>
          <div class="legend">
            <span>low</span><div class="ramp-wrap"><div class="ramp"></div>
              <span class="tick" style="left:0%">0%</span><span class="tick" style="left:25%">25%</span>
              <span class="tick" style="left:50%">50%</span><span class="tick" style="left:75%">75%</span>
              <span class="tick" style="left:90%">90%</span><span class="tick" style="left:100%">100%</span>
            </div><span>high</span>
          </div>
        </section>
        <section><h2 id="pairsTitle">Top similar pairs</h2><div id="pairs" class="pair-list"></div></section>
      </aside>
    </div>
  </section>
  <section id="fasta" class="view"><pre id="fastaText"></pre></section>
</main>
<script id="enzyme-data" type="application/json">__DATA_JSON__</script>
<script>
const DATA = JSON.parse(document.getElementById("enzyme-data").textContent);
const $ = id => document.getElementById(id);
let view = DATA.clustered;
let selected = null;
function esc(s) {
  return String(s ?? "").replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
}
const COLOR_STOPS = [[0,[201,76,76]],[25,[236,167,106]],[50,[242,231,166]],[75,[155,203,156]],[90,[99,182,194]],[100,[63,119,181]]];
function mix(a,b,t){ return Math.round(a+(b-a)*t); }
function color(v){
  const x=Math.max(0,Math.min(100,Number(v)||0));
  for(let i=1;i<COLOR_STOPS.length;i++){
    const [rv,rc]=COLOR_STOPS[i], [lv,lc]=COLOR_STOPS[i-1];
    if(x<=rv){ const t=(x-lv)/Math.max(1,rv-lv); return `rgb(${mix(lc[0],rc[0],t)}, ${mix(lc[1],rc[1],t)}, ${mix(lc[2],rc[2],t)})`; }
  }
  return "rgb(63,119,181)";
}
function textColor(v){ return Number(v)>=94 ? "#fff" : "transparent"; }
function labelFor(m){
  const mode=$("label").value;
  if(mode==="organism") return m.organism || m.accession;
  if(mode==="protein_name") return m.protein_name || m.accession;
  return m.label || m.accession;
}
function searchable(m){ return [m.input_id,m.accession,m.label,m.gene_name,m.organism,m.protein_name,(m.pdb_ids||[]).join(" "),m.nearest_pdb_accession,(m.nearest_pdb_ids||[]).join(" "),m.source].join(" ").toLowerCase(); }
function filteredIndices(){
  const q=$("q").value.trim().toLowerCase();
  const idx=[]; view.meta.forEach((m,i)=>{ if(!q || searchable(m).includes(q)) idx.push(i); });
  return idx;
}
function pdbLinks(ids){
  return (ids||[]).map(p=>`<a href="https://www.rcsb.org/structure/${esc(p)}" target="_blank">${esc(p)}</a>`).join(" ");
}
function accessionLink(m){
  if(m.accession_type==="uniprot") return `<a href="https://www.uniprot.org/uniprotkb/${esc(m.accession)}/entry" target="_blank">${esc(m.accession)}</a>`;
  return `<a href="https://www.ncbi.nlm.nih.gov/protein/${esc(m.accession)}" target="_blank">${esc(m.accession)}</a>`;
}
function renderStats(indices){
  const vals=[];
  for(let a=0;a<indices.length;a++){ for(let b=a+1;b<indices.length;b++) vals.push(view.matrix[indices[a]][indices[b]]); }
  const avg=vals.length ? vals.reduce((a,b)=>a+b,0)/vals.length : 0;
  const direct=view.meta.filter(m=>(m.pdb_ids||[]).length).length;
  const templ=view.meta.filter(m=>(m.nearest_pdb_ids||[]).length).length;
  $("stats").innerHTML=[`Sequences ${indices.length} / ${view.ids.length}`,`Pairs ${vals.length}`,`>=90% ${vals.filter(v=>v>=90).length}`,`mean ${avg.toFixed(1)}%`,`direct PDB ${direct}`,`PDB template ${templ}`].map(x=>`<span class="stat">${esc(x)}</span>`).join("");
}
function renderOverview(){
  const q=$("q").value.trim().toLowerCase();
  const rows=DATA.original.meta.filter(m=>!q || searchable(m).includes(q));
  const head="<tr><th>#</th><th>Accession</th><th>Protein</th><th>Organism</th><th>Length</th><th>Motif</th><th>PDB</th><th>Template</th><th>Source</th><th>Notes</th></tr>";
  const body=rows.map((m,i)=>{
    const motif=m.motif_score==null ? "" : `${Number(m.motif_score).toFixed(0)}%`;
    const warn=(m.warnings||[]).map(w=>`<span class="badge warn">${esc(w)}</span>`).join("");
    const pdb=(m.pdb_ids||[]).length ? pdbLinks(m.pdb_ids) : '<span class="muted">none</span>';
    const template=(m.nearest_pdb_ids||[]).length ? `${esc(m.nearest_pdb_accession)} <span class="muted">${Number(m.nearest_pdb_identity||0).toFixed(1)}%</span><br>${pdbLinks(m.nearest_pdb_ids.slice(0,8))}` : '<span class="muted">none</span>';
    return `<tr><td>${i+1}</td><td class="mono">${accessionLink(m)}</td><td><b>${esc(m.label||m.accession)}</b><br><span class="muted">${esc(m.protein_name||"")}</span></td><td>${esc(m.organism||"")}</td><td>${esc(m.length)}</td><td>${esc(motif)}</td><td>${pdb}</td><td>${template}</td><td>${esc(m.source||"")}</td><td>${warn}</td></tr>`;
  }).join("");
  $("recordTable").innerHTML=head+body;
}
function renderHeatmap(){
  view=DATA[$("order").value];
  const indices=filteredIndices();
  const threshold=Number($("threshold").value||0);
  renderStats(indices);
  const head=`<tr><th class="corner">${indices.length} seq</th>`+indices.map(i=>`<th class="col-head" title="${esc(view.meta[i].protein_name||"")}"><div>${esc(labelFor(view.meta[i]))}</div></th>`).join("")+"</tr>";
  const rows=indices.map(i=>{
    const cells=indices.map(j=>{
      const v=view.matrix[i][j];
      const visible=threshold===0 || v>=threshold || i===j;
      const bg=visible ? color(v) : "#f1f1ec";
      return `<td class="cell ${i===j?"diagonal":""}" style="background:${bg};color:${textColor(v)}" data-i="${i}" data-j="${j}" title="${esc(view.ids[i])} vs ${esc(view.ids[j])}: ${v.toFixed(2)}%">${v>=94 ? v.toFixed(0) : ""}</td>`;
    }).join("");
    return `<tr><th class="row-head" title="${esc(view.meta[i].protein_name||"")}">${esc(labelFor(view.meta[i]))}</th>${cells}</tr>`;
  }).join("");
  $("heatmap").innerHTML=head+rows;
  document.querySelectorAll(".cell").forEach(td=>td.addEventListener("click",()=>selectPair(Number(td.dataset.i),Number(td.dataset.j))));
  if(!selected && indices.length) selectPair(indices[0],indices[0],false);
}
function selectPair(i,j,rerender=true){
  selected=[i,j];
  const a=view.meta[i], b=view.meta[j], v=view.matrix[i][j];
  $("detail").innerHTML=[
    ["Query", accessionLink(a)],["Target", accessionLink(b)],["Identity", `${v.toFixed(2)}%`],
    ["Query protein", esc(a.protein_name||"")],["Target protein", esc(b.protein_name||"")],
    ["Query PDB", pdbLinks(a.pdb_ids)||'<span class="muted">none</span>'],["Target PDB", pdbLinks(b.pdb_ids)||'<span class="muted">none</span>']
  ].map(([k,val])=>`<div>${esc(k)}</div><div>${val}</div>`).join("");
  if(rerender) renderPairs();
}
function renderPairs(){
  const mode=$("pairMode").value;
  const list=mode==="similar" ? DATA.topSimilar : DATA.topDistant;
  $("pairsTitle").textContent=mode==="similar" ? "Top similar pairs" : "Most distant pairs";
  $("pairs").innerHTML=list.map(p=>{
    const pct=Number(p.identity_percent)||0;
    return `<div class="pair" data-a="${esc(p.query)}" data-b="${esc(p.target)}"><b>${esc(p.query)} vs ${esc(p.target)}</b><br><span class="muted">${pct.toFixed(2)}%</span><div class="bar"><span style="width:${pct}%;background:${color(pct)}"></span></div></div>`;
  }).join("");
  document.querySelectorAll(".pair").forEach(el=>el.addEventListener("click",()=>{
    const a=view.ids.indexOf(el.dataset.a), b=view.ids.indexOf(el.dataset.b);
    if(a>=0 && b>=0) selectPair(a,b);
  }));
}
function renderFasta(){
  $("fastaText").textContent=DATA.original.meta.filter(m=>m.sequence).map(m=>{
    const label=m.label||m.accession, org=m.organism||"unknown", name=(m.protein_name||"").replace(/\|/g,";");
    const seq=String(m.sequence||"").replace(/(.{1,80})/g,"$1\n").trim();
    return `>${m.accession}|label=${label}|organism=${org}|name=${name}\n${seq}`;
  }).join("\n");
}
function renderAll(){ renderOverview(); renderHeatmap(); renderPairs(); renderFasta(); }
document.querySelectorAll(".tab").forEach(btn=>btn.addEventListener("click",()=>{
  document.querySelectorAll(".tab").forEach(x=>x.classList.remove("active")); btn.classList.add("active");
  document.querySelectorAll(".view").forEach(x=>x.classList.remove("active")); $(btn.dataset.view).classList.add("active");
}));
["q","order","label","threshold","pairMode"].forEach(id=>$(id).addEventListener("input",renderAll));
renderAll();
</script>
</body>
</html>
"""


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>__TITLE__</title>
  <style>
    :root {
      --bg:#f6f6f2; --panel:#fff; --ink:#202326; --muted:#687076; --line:#d8ddd6;
      --blue:#2f5f7f; --green:#3f7155; --amber:#9a6b2f; --red:#9a4a43; --soft:#eef1ed;
      --shadow:0 12px 30px rgba(33,37,41,.08);
    }
    * { box-sizing:border-box; }
    body { margin:0; background:var(--bg); color:var(--ink); font-family:Arial,"Microsoft YaHei",sans-serif; }
    header { position:sticky; top:0; z-index:20; background:#fff; border-bottom:1px solid var(--line); padding:14px 18px 12px; }
    h1 { margin:0 0 10px; font-size:24px; letter-spacing:0; }
    h2 { margin:0 0 4px; font-size:22px; letter-spacing:0; }
    h3 { margin:14px 0 8px; font-size:15px; }
    button,input,select { font:inherit; }
    .controls { display:grid; grid-template-columns:minmax(220px,1.4fr) repeat(3,minmax(130px,180px)); gap:10px; align-items:center; }
    input,select { width:100%; border:1px solid var(--line); border-radius:6px; background:#fff; padding:8px 10px; font-size:14px; }
    .stats { display:flex; flex-wrap:wrap; gap:8px; margin-top:10px; color:var(--muted); font-size:13px; }
    .stat { background:var(--soft); border-radius:6px; padding:5px 8px; }
    main { display:grid; grid-template-columns:360px minmax(0,1fr); gap:14px; padding:14px 18px 22px; }
    .list-panel,.detail { background:var(--panel); border:1px solid var(--line); border-radius:8px; box-shadow:var(--shadow); }
    .list-panel { max-height:calc(100vh - 142px); overflow:auto; padding:10px; }
    .detail { min-height:520px; padding:16px; overflow:hidden; }
    .card { border:1px solid var(--line); border-radius:8px; padding:10px; margin-bottom:8px; cursor:pointer; background:#fff; }
    .card:hover { border-color:var(--blue); }
    .card.active { border-color:var(--blue); box-shadow:inset 4px 0 0 var(--blue); }
    .card.warn { background:#fff8f4; }
    .card-title { display:flex; justify-content:space-between; gap:8px; align-items:flex-start; font-weight:700; }
    .card-title span:first-child { overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }
    .meta { color:var(--muted); font-size:12px; line-height:1.35; }
    .chips { display:flex; flex-wrap:wrap; gap:4px; margin-top:7px; }
    .chip { display:inline-block; border-radius:999px; padding:2px 7px; background:var(--soft); font-size:12px; white-space:nowrap; }
    .chip.pdb { background:#e7f0f5; color:#214b65; }
    .chip.warn { background:#faeadf; color:#7a3a22; }
    .chip.good { background:#e6f1e8; color:#2d5d3f; }
    .toolbar { display:flex; flex-wrap:wrap; gap:6px; margin:12px 0 14px; }
    .tab { border:1px solid var(--line); background:#fff; border-radius:6px; padding:7px 10px; cursor:pointer; }
    .tab.active { background:var(--ink); color:#fff; border-color:var(--ink); }
    .kv { display:grid; grid-template-columns:150px minmax(0,1fr); gap:7px 10px; font-size:13px; line-height:1.4; }
    .kv div:nth-child(odd) { color:var(--muted); }
    .two-col { display:grid; grid-template-columns:minmax(0,1fr) minmax(280px,.55fr); gap:14px; align-items:start; }
    .table-wrap { border:1px solid var(--line); border-radius:8px; overflow:auto; background:#fff; }
    table { border-collapse:separate; border-spacing:0; width:100%; }
    th,td { border-bottom:1px solid var(--line); padding:7px 8px; text-align:left; vertical-align:top; font-size:12px; }
    th { position:sticky; top:0; z-index:2; background:#fff; color:var(--muted); font-weight:600; }
    tr:last-child td { border-bottom:0; }
    a { color:var(--blue); text-decoration:none; }
    a:hover { text-decoration:underline; }
    .mono { font-family:Consolas,Menlo,monospace; }
    pre { margin:0; white-space:pre-wrap; word-break:break-word; background:#fff; border:1px solid var(--line); border-radius:8px; padding:12px; font-size:12px; line-height:1.45; max-height:420px; overflow:auto; }
    .identity-summary { display:grid; grid-template-columns:repeat(5,minmax(92px,1fr)); gap:8px; margin:8px 0 10px; }
    .identity-summary .stat { text-align:center; }
    .identity-layout { display:grid; grid-template-columns:minmax(0,1fr) 340px; gap:12px; align-items:start; }
    .identity-heatmap-wrap { overflow:auto; max-height:560px; border:1px solid var(--line); border-radius:8px; background:#fff; }
    table.identity-heatmap { border-collapse:separate; border-spacing:0; font-size:10px; width:auto; min-width:max-content; }
    .identity-heatmap th,.identity-heatmap td { border-right:1px solid rgba(255,255,255,.42); border-bottom:1px solid rgba(255,255,255,.42); min-width:34px; height:28px; text-align:center; padding:0 2px; }
    .identity-corner { position:sticky; top:0; left:0; z-index:5; min-width:142px; background:#fff; color:var(--muted); border-right:1px solid var(--line); border-bottom:1px solid var(--line); }
    .identity-col-head { position:sticky; top:0; z-index:4; height:142px; min-width:34px; vertical-align:bottom; background:#fff; border-bottom:1px solid var(--line); }
    .identity-col-head > div { writing-mode:vertical-rl; transform:rotate(180deg); white-space:nowrap; max-height:134px; overflow:hidden; padding:4px 2px; color:var(--muted); font-weight:400; }
    .identity-row-head { position:sticky; left:0; z-index:3; min-width:142px; max-width:142px; background:#fff; text-align:right; padding:0 7px; border-right:1px solid var(--line); white-space:nowrap; overflow:hidden; text-overflow:ellipsis; color:var(--muted); font-weight:400; }
    .identity-cell { cursor:pointer; font-size:10px; }
    .identity-cell:hover { outline:2px solid #202326; outline-offset:-2px; }
    .identity-cell.diagonal { box-shadow:inset 0 0 0 1px rgba(0,0,0,.24); font-weight:700; }
    .identity-cell.focus { box-shadow:inset 0 0 0 2px #202326; }
    .identity-pair-detail { display:grid; grid-template-columns:116px 1fr; gap:6px 9px; font-size:13px; line-height:1.35; border:1px solid var(--line); border-radius:8px; padding:10px; background:#fff; }
    .identity-pair-detail div:nth-child(odd) { color:var(--muted); }
    .pair-list { display:flex; flex-direction:column; gap:6px; max-height:280px; overflow:auto; }
    .pair { border:1px solid var(--line); border-radius:6px; padding:7px; cursor:pointer; background:#fff; font-size:12px; }
    .pair:hover { border-color:var(--blue); }
    .bar { height:7px; border-radius:999px; margin-top:5px; background:var(--soft); overflow:hidden; }
    .bar span { display:block; height:100%; }
    .legend { display:grid; grid-template-columns:36px 1fr 42px; gap:8px; align-items:start; font-size:12px; color:var(--muted); margin-top:8px; }
    .ramp-wrap { position:relative; padding-bottom:22px; }
    .ramp { height:12px; border-radius:999px; border:1px solid var(--line); background:linear-gradient(90deg,#C94C4C 0%,#ECA76A 25%,#F2E7A6 50%,#9BCB9C 75%,#63B6C2 90%,#3F77B5 100%); }
    .tick { position:absolute; top:15px; transform:translateX(-50%); white-space:nowrap; font-size:11px; }
    .tick::before { content:""; display:block; width:1px; height:6px; margin:0 auto 2px; background:var(--line); }
    @media (max-width:1050px) {
      main,.two-col,.identity-layout { grid-template-columns:1fr; }
      .list-panel { max-height:none; }
      .controls { grid-template-columns:1fr 1fr; }
      .detail { min-height:0; }
    }
  </style>
</head>
<body>
<header>
  <h1>__TITLE__</h1>
  <div class="controls">
    <input id="q" placeholder="Search accession, gene, organism, protein, PDB">
    <select id="order"><option value="original">Input order</option><option value="clustered">Clustered identity order</option></select>
    <select id="label"><option value="accession">Accession labels</option><option value="label">Gene/accession labels</option><option value="organism">Organism labels</option></select>
    <select id="threshold"><option value="0">All matrix cells</option><option value="70">Dim below 70%</option><option value="80">Dim below 80%</option><option value="90">Dim below 90%</option><option value="95">Dim below 95%</option></select>
  </div>
  <div class="stats" id="stats"></div>
</header>
<main>
  <section class="list-panel" id="list"></section>
  <section class="detail" id="detail"></section>
</main>
<script id="enzyme-data" type="application/json">__DATA_JSON__</script>
<script>
const DATA = JSON.parse(document.getElementById("enzyme-data").textContent);
const $ = id => document.getElementById(id);
let selectedAcc = DATA.original.ids[0] || "";
let activeTab = "overview";
let selectedPair = null;
function esc(s) {
  return String(s ?? "").replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
}
function metaByAcc(acc) { return DATA.original.meta.find(m=>m.accession===acc) || DATA.clustered.meta.find(m=>m.accession===acc) || DATA.original.meta[0]; }
function currentView() { return DATA[$("order").value] || DATA.original; }
function displayLabel(m) {
  const mode = $("label").value;
  if (mode === "organism") return m.organism || m.accession;
  if (mode === "label") return m.label || m.accession;
  return m.accession;
}
function searchable(m) {
  return [m.input_id,m.accession,m.label,m.gene_name,m.organism,m.protein_name,(m.pdb_ids||[]).join(" "),m.nearest_pdb_accession,(m.nearest_pdb_ids||[]).join(" "),m.source,(m.warnings||[]).join(" ")].join(" ").toLowerCase();
}
function filteredOriginal() {
  const q = $("q").value.trim().toLowerCase();
  return DATA.original.meta.filter(m => !q || searchable(m).includes(q));
}
function filteredViewIndices() {
  const q = $("q").value.trim().toLowerCase();
  const view = currentView();
  const out = [];
  view.meta.forEach((m,i)=>{ if(!q || searchable(m).includes(q)) out.push(i); });
  return out;
}
function accessionLink(m) {
  if (!m) return "";
  if (m.accession_type === "uniprot") return `<a href="https://www.uniprot.org/uniprotkb/${esc(m.accession)}/entry" target="_blank">${esc(m.accession)}</a>`;
  return `<a href="https://www.ncbi.nlm.nih.gov/protein/${esc(m.accession)}" target="_blank">${esc(m.accession)}</a>`;
}
function pdbLinks(ids) {
  return (ids || []).map(p=>`<a href="https://www.rcsb.org/structure/${esc(p)}" target="_blank">${esc(p)}</a>`).join(" ");
}
function motifText(m) {
  if (m.motif_score == null) return "not scored";
  const pos = Object.entries(m.motif_positions || {}).map(([k,v])=>`${k}:${v >= 0 ? v + 1 : "absent"}`).join(" · ");
  return `${Number(m.motif_score).toFixed(0)}%${pos ? " · " + pos : ""}`;
}
function wrapSeq(seq) {
  return String(seq || "").match(/.{1,80}/g)?.join("\n") || "";
}
function fastaFor(m) {
  const header = m.fasta_header || `${m.accession}|label=${m.label||m.accession}|organism=${m.organism||"unknown"}|name=${(m.protein_name||"").replace(/\|/g,";")}`;
  return `>${header}\n${wrapSeq(m.sequence)}`;
}
function allFasta() {
  return DATA.original.meta.filter(m=>m.sequence).map(fastaFor).join("\n\n");
}
const COLOR_STOPS = [[0,[201,76,76]],[25,[236,167,106]],[50,[242,231,166]],[75,[155,203,156]],[90,[99,182,194]],[100,[63,119,181]]];
function mix(a,b,t){ return Math.round(a+(b-a)*t); }
function rgbFor(v) {
  const x=Math.max(0,Math.min(100,Number(v)||0));
  for(let i=1;i<COLOR_STOPS.length;i++){
    const [rv,rc]=COLOR_STOPS[i], [lv,lc]=COLOR_STOPS[i-1];
    if(x<=rv){ const t=(x-lv)/Math.max(1,rv-lv); return [mix(lc[0],rc[0],t),mix(lc[1],rc[1],t),mix(lc[2],rc[2],t)]; }
  }
  return [63,119,181];
}
function color(v) { const c = rgbFor(v); return `rgb(${c[0]}, ${c[1]}, ${c[2]})`; }
function textColor(v) {
  const [r,g,b] = rgbFor(v);
  const luminance = (0.299*r + 0.587*g + 0.114*b) / 255;
  return luminance < 0.52 ? "#fff" : "#202326";
}
function pairIdentity(a,b) {
  const p = DATA.pairs.find(x => (x.query===a && x.target===b) || (x.query===b && x.target===a));
  return p ? Number(p.identity_percent) : (a===b ? 100 : null);
}
function pairRowsFor(acc) {
  return DATA.pairs.filter(p=>p.query===acc || p.target===acc).sort((a,b)=>Number(b.identity_percent)-Number(a.identity_percent));
}
function renderStats() {
  const rows = filteredOriginal();
  const vals = DATA.pairs.map(p=>Number(p.identity_percent)).filter(v=>Number.isFinite(v));
  const mean = vals.length ? vals.reduce((a,b)=>a+b,0)/vals.length : 0;
  const direct = DATA.original.meta.filter(m=>(m.pdb_ids||[]).length).length;
  const templ = DATA.original.meta.filter(m=>(m.nearest_pdb_ids||[]).length).length;
  const warnings = DATA.original.meta.filter(m=>(m.warnings||[]).length).length;
  $("stats").innerHTML = [
    `Showing ${rows.length} / ${DATA.original.meta.length}`,
    `FASTA ${DATA.summary.resolved_sequences}`,
    `pairs ${DATA.pairs.length}`,
    `mean ${mean.toFixed(1)}%`,
    `direct PDB ${direct}`,
    `PDB template ${templ}`,
    `warnings ${warnings}`
  ].map(x=>`<span class="stat">${esc(x)}</span>`).join("");
}
function renderList() {
  const rows = filteredOriginal();
  $("list").innerHTML = rows.map(m=>{
    const warn = (m.warnings||[]).length ? " warn" : "";
    const active = m.accession === selectedAcc ? " active" : "";
    const pdb = (m.pdb_ids||[]).length ? `<span class="chip pdb">PDB ${(m.pdb_ids||[]).length}</span>` : "";
    const motif = m.motif_score == null ? "" : `<span class="chip ${m.motif_score ? "good" : "warn"}">motif ${Number(m.motif_score).toFixed(0)}%</span>`;
    const template = m.nearest_pdb_accession ? `<span class="chip">template ${esc(m.nearest_pdb_accession)} ${Number(m.nearest_pdb_identity||0).toFixed(1)}%</span>` : "";
    return `<article class="card${active}${warn}" data-acc="${esc(m.accession)}">
      <div class="card-title"><span>${esc(displayLabel(m))}</span><span class="meta">${esc(m.length)} aa</span></div>
      <div class="meta mono">${esc(m.accession)} · ${esc(m.accession_type||"")}</div>
      <div class="meta">${esc(m.protein_name||"")}</div>
      <div class="meta">${esc(m.organism||"")}</div>
      <div class="chips">${motif}${pdb}${template}${(m.warnings||[]).map(w=>`<span class="chip warn">${esc(w)}</span>`).join("")}</div>
    </article>`;
  }).join("") || "<p class='meta'>No records match the search.</p>";
  document.querySelectorAll(".card").forEach(el=>el.addEventListener("click",()=>{
    selectedAcc = el.dataset.acc;
    activeTab = "overview";
    selectedPair = null;
    render();
  }));
}
function tabButton(name,label) {
  return `<button class="tab ${activeTab===name ? "active" : ""}" data-tab="${name}">${esc(label)}</button>`;
}
function overviewHtml(m) {
  const warnings = (m.warnings||[]).map(w=>`<span class="chip warn">${esc(w)}</span>`).join("") || '<span class="meta">none</span>';
  const directPdb = (m.pdb_ids||[]).length ? pdbLinks(m.pdb_ids) : '<span class="meta">none</span>';
  const template = (m.nearest_pdb_ids||[]).length ? `${accessionLink(metaByAcc(m.nearest_pdb_accession))} <span class="chip">${Number(m.nearest_pdb_identity||0).toFixed(1)}%</span><br>${pdbLinks(m.nearest_pdb_ids)}` : '<span class="meta">none</span>';
  return `<div class="two-col">
    <div class="kv">
      <div>Accession</div><div class="mono">${accessionLink(m)}</div>
      <div>Label</div><div>${esc(m.label||m.accession)}</div>
      <div>Protein</div><div>${esc(m.protein_name||"")}</div>
      <div>Organism</div><div>${esc(m.organism||"")}</div>
      <div>Length</div><div>${esc(m.length)} aa</div>
      <div>Source</div><div>${esc(m.source||"")}</div>
      <div>Motif</div><div>${esc(motifText(m))}</div>
      <div>Direct PDB</div><div>${directPdb}</div>
      <div>Nearest PDB template</div><div>${template}</div>
      <div>Warnings</div><div>${warnings}</div>
    </div>
    <pre>${esc(fastaFor(m))}</pre>
  </div>`;
}
function fastaHtml(m) {
  const rows = DATA.original.meta.map((x,i)=>`<tr><td>${i+1}</td><td class="mono">${accessionLink(x)}</td><td>${esc(x.label||x.accession)}</td><td>${esc(x.length)}</td><td>${esc(x.motif_score == null ? "" : Number(x.motif_score).toFixed(0)+"%")}</td><td>${esc(x.organism||"")}</td></tr>`).join("");
  return `<div class="two-col">
    <section>
      <h3>Selected protein FASTA</h3>
      <pre>${esc(fastaFor(m))}</pre>
      <h3>All protein FASTA</h3>
      <pre>${esc(allFasta())}</pre>
    </section>
    <section>
      <h3>FASTA index</h3>
      <div class="table-wrap"><table><tr><th>#</th><th>Accession</th><th>Label</th><th>aa</th><th>Motif</th><th>Organism</th></tr>${rows}</table></div>
      <h3>Selected metadata</h3>
      <div class="kv">
        <div>Header</div><div class="mono">${esc(m.fasta_header||"")}</div>
        <div>MD5 source</div><div>${esc(m.source||"")}</div>
        <div>PDB template</div><div>${esc(m.nearest_pdb_accession||"")} ${m.nearest_pdb_identity == null ? "" : Number(m.nearest_pdb_identity).toFixed(1)+"%"}</div>
        <div>Motif positions</div><div>${esc(motifText(m))}</div>
      </div>
    </section>
  </div>`;
}
function pairDetailHtml(a,b) {
  const ma = metaByAcc(a), mb = metaByAcc(b), v = pairIdentity(a,b);
  return `<div>Query</div><div>${accessionLink(ma)} · ${esc(ma.label||"")}</div>
    <div>Target</div><div>${accessionLink(mb)} · ${esc(mb.label||"")}</div>
    <div>Identity</div><div><span class="chip">${v == null ? "NA" : v.toFixed(2)+"%"}</span></div>
    <div>Query organism</div><div>${esc(ma.organism||"")}</div>
    <div>Target organism</div><div>${esc(mb.organism||"")}</div>
    <div>Query PDB</div><div>${pdbLinks(ma.pdb_ids)||'<span class="meta">none</span>'}</div>
    <div>Target PDB</div><div>${pdbLinks(mb.pdb_ids)||'<span class="meta">none</span>'}</div>`;
}
function renderPairDetail(a,b) {
  const el = $("pairDetail");
  if (el) el.innerHTML = pairDetailHtml(a,b);
}
function heatmapHtml(m) {
  const view = currentView();
  const indices = filteredViewIndices();
  const threshold = Number($("threshold").value || 0);
  const ids = indices.map(i=>view.ids[i]);
  const selected = ids.includes(m.accession) ? m.accession : ids[0] || m.accession;
  const vals = [];
  for (let a=0;a<indices.length;a++) for (let b=a+1;b<indices.length;b++) vals.push(view.matrix[indices[a]][indices[b]]);
  const avg = vals.length ? vals.reduce((a,b)=>a+b,0)/vals.length : 0;
  const summary = `<div class="identity-summary">${[
    `seq ${indices.length}`,
    `pairs ${vals.length}`,
    `>=90% ${vals.filter(v=>v>=90).length}`,
    `<80% ${vals.filter(v=>v<80).length}`,
    `mean ${avg.toFixed(1)}%`
  ].map(x=>`<span class="stat">${esc(x)}</span>`).join("")}</div>`;
  const head = `<tr><th class="identity-corner">${indices.length} seq</th>` + indices.map(i=>`<th class="identity-col-head" title="${esc(view.meta[i].protein_name||"")}"><div>${esc(displayLabel(view.meta[i]))}</div></th>`).join("") + `</tr>`;
  const rows = indices.map(i=>{
    const a = view.ids[i];
    const cells = indices.map(j=>{
      const b = view.ids[j];
      const v = Number(view.matrix[i][j] || 0);
      const dim = threshold && v < threshold && i !== j;
      const bg = dim ? "#f1f1ec" : color(v);
      const cls = ["identity-cell", i===j ? "diagonal" : "", (a===selected || b===selected) ? "focus" : ""].filter(Boolean).join(" ");
      return `<td class="${cls}" style="background:${bg};color:${dim ? "#9aa19a" : textColor(v)}" data-a="${esc(a)}" data-b="${esc(b)}" title="${esc(a)} vs ${esc(b)}: ${v.toFixed(2)}%">${v.toFixed(1)}</td>`;
    }).join("");
    return `<tr><th class="identity-row-head" title="${esc(view.meta[i].protein_name||"")}">${esc(displayLabel(view.meta[i]))}</th>${cells}</tr>`;
  }).join("");
  const pairs = pairRowsFor(m.accession).slice(0,18).map(p=>{
    const other = p.query === m.accession ? p.target : p.query;
    const pct = Number(p.identity_percent);
    return `<div class="pair" data-a="${esc(m.accession)}" data-b="${esc(other)}"><b>${esc(m.accession)} vs ${esc(other)}</b><br><span class="meta">${pct.toFixed(2)}%</span><div class="bar"><span style="width:${pct}%;background:${color(pct)}"></span></div></div>`;
  }).join("");
  const selectedB = selectedPair ? selectedPair[1] : selected;
  return `${summary}<div class="identity-layout">
    <section>
      <div class="identity-heatmap-wrap"><table class="identity-heatmap">${head}${rows}</table></div>
      <div class="legend"><span>low</span><div class="ramp-wrap"><div class="ramp"></div>
        <span class="tick" style="left:0%">0%</span><span class="tick" style="left:50%">50%</span><span class="tick" style="left:75%">75%</span><span class="tick" style="left:90%">90%</span><span class="tick" style="left:100%">100%</span>
      </div><span>high</span></div>
    </section>
    <aside>
      <h3>Selected pair</h3>
      <div id="pairDetail" class="identity-pair-detail">${pairDetailHtml(selected, selectedB)}</div>
      <h3>Pairs for ${esc(displayLabel(m))}</h3>
      <div class="pair-list">${pairs || "<p class='meta'>No pairwise identity rows.</p>"}</div>
    </aside>
  </div>`;
}
function renderDetail() {
  const m = metaByAcc(selectedAcc);
  if (!m) { $("detail").innerHTML = "<p>No records.</p>"; return; }
  const titleWarn = (m.warnings||[]).length ? `<span class="chip warn">${esc((m.warnings||[])[0])}</span>` : "";
  const body = activeTab === "fasta" ? fastaHtml(m) : activeTab === "identity" ? heatmapHtml(m) : overviewHtml(m);
  $("detail").innerHTML = `<h2>${esc(displayLabel(m))} ${titleWarn}</h2>
    <div class="meta">${esc(m.protein_name||"")} · ${esc(m.organism||"")}</div>
    <div class="toolbar">${tabButton("overview","Overview")}${tabButton("fasta","FASTA")}${tabButton("identity","Identity")}</div>
    ${body}`;
  document.querySelectorAll(".tab").forEach(btn=>btn.addEventListener("click",()=>{
    activeTab = btn.dataset.tab;
    selectedPair = null;
    renderDetail();
  }));
  document.querySelectorAll(".identity-cell").forEach(td=>td.addEventListener("click",()=>{
    selectedPair = [td.dataset.a, td.dataset.b];
    renderPairDetail(td.dataset.a, td.dataset.b);
  }));
  document.querySelectorAll(".pair").forEach(el=>el.addEventListener("click",()=>{
    selectedPair = [el.dataset.a, el.dataset.b];
    renderPairDetail(el.dataset.a, el.dataset.b);
  }));
}
function render() {
  if (!metaByAcc(selectedAcc)) selectedAcc = DATA.original.ids[0] || "";
  renderStats();
  renderList();
  renderDetail();
}
["q","order","label","threshold"].forEach(id=>$(id).addEventListener("input",render));
render();
</script>
</body>
</html>
"""


def write_html(data: dict[str, Any], path: Path, title: str) -> None:
    json_data = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    rendered = HTML_TEMPLATE.replace("__TITLE__", html.escape(title)).replace("__DATA_JSON__", json_data)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(rendered)


def write_json(data: dict[str, Any], path: Path) -> None:
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def command_build(args: argparse.Namespace) -> None:
    ids = split_ids(args.ids or "")
    if args.ids_file:
        ids.extend(split_ids(Path(args.ids_file).read_text(encoding="utf-8", errors="replace")))
        ids = split_ids(" ".join(ids))
    if not ids:
        raise SystemExit("No IDs supplied")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    motifs = split_ids(args.motifs)
    records = resolve_records(ids, Path(args.local_csv), Path(args.local_fasta), args.fetch, args.sleep)
    apply_motif_scores(records, motifs)
    matrix, pairs = build_identity(records)
    assign_nearest_pdb(records, matrix)

    metadata_rows = [record_to_row(rec) for rec in records]
    pair_rows = pairs
    write_csv(metadata_rows, outdir / "enzyme_library_metadata.csv")
    write_csv(pair_rows, outdir / "pairwise_identity.csv")
    write_fasta(records, outdir / "enzyme_library.fasta")
    data = build_html_data(records, matrix, pairs, args.title)
    write_json(data, outdir / "enzyme_library.json")
    write_html(data, outdir / "enzyme_library_dashboard.html", args.title)

    unresolved = [rec.input_id for rec in records if not rec.sequence]
    warnings = sum(len(rec.warnings) for rec in records)
    print(f"records={len(records)}")
    print(f"resolved_sequences={len(records) - len(unresolved)}")
    print(f"direct_pdb_records={sum(1 for rec in records if rec.pdb_ids)}")
    print(f"records_with_pdb_template={sum(1 for rec in records if rec.nearest_pdb_ids)}")
    print(f"pairwise_comparisons={len(pairs)}")
    print(f"warnings={warnings}")
    if unresolved:
        print("unresolved=" + ",".join(unresolved))
    print(f"html={outdir / 'enzyme_library_dashboard.html'}")
    print(f"metadata={outdir / 'enzyme_library_metadata.csv'}")
    print(f"fasta={outdir / 'enzyme_library.fasta'}")


def command_check(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    data_path = Path(args.json) if args.json else outdir / "enzyme_library.json"
    html_path = Path(args.html) if args.html else outdir / "enzyme_library_dashboard.html"
    fasta_path = Path(args.fasta) if args.fasta else outdir / "enzyme_library.fasta"
    pairs_path = Path(args.pairs) if args.pairs else outdir / "pairwise_identity.csv"
    failures: list[str] = []
    for path in (data_path, html_path, fasta_path, pairs_path):
        if not path.exists():
            failures.append(f"missing {path}")
    if failures:
        print("\n".join(failures))
        raise SystemExit(2)
    data = json.loads(data_path.read_text(encoding="utf-8"))
    meta = data.get("original", {}).get("meta", [])
    matrix = data.get("original", {}).get("matrix", [])
    n = len(meta)
    if len(matrix) != n or any(len(row) != n for row in matrix):
        failures.append("identity matrix is not square")
    unresolved = [m.get("accession") for m in meta if not m.get("sequence")]
    if unresolved:
        failures.append("unresolved sequences: " + ",".join(unresolved))
    expected_pairs = n * (n - 1) // 2
    actual_pairs = sum(1 for _ in pairs_path.open("r", encoding="utf-8-sig", errors="replace")) - 1
    if actual_pairs != expected_pairs:
        failures.append(f"pairwise CSV rows {actual_pairs} != expected {expected_pairs}")
    html_text = html_path.read_text(encoding="utf-8", errors="replace")
    has_render_code = "renderHeatmap" in html_text or ("identity-heatmap" in html_text and "renderDetail" in html_text)
    if 'id="enzyme-data"' not in html_text or not has_render_code:
        failures.append("HTML does not contain embedded data/render code")
    fasta_records = sum(1 for line in fasta_path.read_text(encoding="utf-8", errors="replace").splitlines() if line.startswith(">"))
    if fasta_records != n - len(unresolved):
        failures.append(f"FASTA records {fasta_records} != resolved sequence count {n - len(unresolved)}")
    pdb_templates = sum(1 for m in meta if m.get("nearest_pdb_ids"))
    if failures:
        for failure in failures:
            print(f"FAIL {failure}")
        raise SystemExit(2)
    print(f"OK records={n}")
    print(f"OK pairwise_rows={actual_pairs}")
    print(f"OK fasta_records={fasta_records}")
    print(f"OK pdb_template_records={pdb_templates}")
    print(f"OK html={html_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build an ID-list enzyme library dashboard with FASTA, PDB links, and identity matrix.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    build = sub.add_parser("build")
    build.add_argument("--ids", default="", help="Semicolon/comma/whitespace separated UniProt or NCBI protein IDs.")
    build.add_argument("--ids-file", default="", help="Optional text file with IDs.")
    build.add_argument("--outdir", required=True)
    build.add_argument("--title", default="Enzyme Library Dashboard")
    build.add_argument("--local-csv", default=str(DEFAULT_LOCAL_CSV))
    build.add_argument("--local-fasta", default=str(DEFAULT_LOCAL_FASTA))
    build.add_argument("--motifs", default="", help="Optional motif list, e.g. 'GTTDDS,APNNGLL,FADAG'.")
    build.add_argument("--fetch", action="store_true", help="Fetch missing FASTA and UniProt PDB metadata from public APIs.")
    build.add_argument("--sleep", type=float, default=0.15, help="Delay between public API requests.")
    build.set_defaults(func=command_build)

    check = sub.add_parser("check")
    check.add_argument("--outdir", required=True)
    check.add_argument("--json", default="")
    check.add_argument("--html", default="")
    check.add_argument("--fasta", default="")
    check.add_argument("--pairs", default="")
    check.set_defaults(func=command_check)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
