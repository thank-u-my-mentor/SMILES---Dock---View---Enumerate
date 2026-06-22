#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


KINGDOM_KEYWORDS = [
    ("Bacteria", ["bacteria", "bacterium", "escherichia", "pseudomonas", "staphylococcus", "xanthomonas", "bacillus", "streptomyces", "clostridium", "salmonella", "klebsiella", "mycobacterium", "lactobacillus", "propionibacterium", "sporosarcina", "rhizobium"]),
    ("Archaea", ["archaea", "archaeon", "methano", "halobacter"]),
    ("Fungi", ["fungi", "fungus", "saccharomyces", "aspergillus", "candida", "neurospora", "fusarium", "penicillium", "steccherinum"]),
    ("Plant", ["viridiplantae", "plantae", "arabidopsis", "oryza sativa", "zea mays", "canavalia", "phaseolus"]),
    ("Animal", ["metazoa", "animalia", "homo sapiens", "rattus", "mus musculus", "drosophila", "danio", "xenopus", "bos taurus"]),
    ("Virus", ["virus", "viruses", "phage"]),
    ("Metagenome", ["metagenome", "environmental sample"]),
]


def infer_kingdom(text: str) -> str:
    lower = (text or "").lower()
    for kingdom, keys in KINGDOM_KEYWORDS:
        if any(key in lower for key in keys):
            return kingdom
    return "Unknown"


def read_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: List[Dict[str, str]], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean_seq(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq or "").upper()


def wrap(seq: str, width: int = 80) -> Iterable[str]:
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


def primary_uniprot(text: str) -> str:
    value = (text or "").strip()
    if not value:
        return ""
    value = value.replace(",", ";")
    return value.split(";")[0].strip()


def safe_id(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.:-]", "_", text or "")
    return text or "seq"


def parse_fasta(path: Path) -> List[Tuple[str, str]]:
    rows: List[Tuple[str, str]] = []
    if not path.exists():
        return rows
    header = ""
    seq: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header and seq:
                    rows.append((header, clean_seq("".join(seq))))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line)
    if header and seq:
        rows.append((header, clean_seq("".join(seq))))
    return rows


def parse_header_meta(header: str) -> Dict[str, str]:
    meta = {"raw_header": header}
    first = header.split()[0] if header else ""
    meta["header_id"] = first
    if "|" in first:
        for part in first.split("|"):
            if part.startswith("uniprot="):
                meta["uniprot_ids"] = part.split("=", 1)[1]
            elif part.startswith("organism="):
                meta["organism"] = part.split("=", 1)[1]
            elif part.startswith("length="):
                meta["sequence_length"] = part.split("=", 1)[1]
            elif part.startswith("chains="):
                meta["chains"] = part.split("=", 1)[1]
        m = re.match(r"^([0-9A-Za-z]{4})_([^|]+)", first)
        if m:
            meta["pdb_id"] = m.group(1).upper()
            meta["entity_id"] = m.group(2)
    if " OS=" in header:
        meta["organism"] = header.split(" OS=", 1)[1].split(" OX=", 1)[0].strip()
    return meta


def load_literature(project: Path) -> Dict[str, Dict[str, str]]:
    by_pdb: Dict[str, Dict[str, str]] = {}
    for path in [project / "Metal-F_literature_seed.csv", *project.glob("Metal-*/*_literature_seed.csv")]:
        set_name = "Metal-F" if path.parent == project else path.parent.name.replace("Metal-", "")
        for row in read_csv(path):
            pdbs = re.split(r"[;,]\s*", row.get("PDB ID", "") or "")
            for pdb in pdbs:
                pdb = pdb.strip().upper()
                if not pdb:
                    continue
                out = dict(row)
                out["literature_set"] = set_name
                by_pdb[pdb] = out
    return by_pdb


def collect_csv_rows(project: Path, min_length: int, exclude_pdb: set[str]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    files = [project / "Metal-F_protein_entities.csv", *project.glob("Metal-*/*_protein_entities_kept.csv")]
    literature = load_literature(project)
    for path in files:
        set_name = "Metal-F" if path.parent == project else path.parent.name.replace("Metal-", "")
        for row in read_csv(path):
            seq = clean_seq(row.get("sequence", ""))
            if len(seq) < min_length:
                continue
            pdb_id = (row.get("pdb_id", "") or "").upper()
            if pdb_id in exclude_pdb:
                continue
            lit = literature.get(pdb_id, {})
            uid = primary_uniprot(row.get("uniprot_ids", ""))
            seq_id = uid or f"{pdb_id}_{row.get('entity_id', '1')}"
            rows.append({
                "seed_id": safe_id(seq_id),
                "set_name": set_name,
                "source_file": str(path),
                "pdb_id": pdb_id,
                "entity_id": row.get("entity_id", ""),
                "chains": row.get("chains", ""),
                "sequence_length": str(len(seq)),
                "uniprot_ids": row.get("uniprot_ids", ""),
                "organism": row.get("organism", ""),
                "taxonomy_id": row.get("taxonomy_id", ""),
                "kingdom": infer_kingdom(row.get("organism", "")),
                "doi": lit.get("DOI", ""),
                "title": lit.get("Title", ""),
                "pdb_url": f"https://www.rcsb.org/structure/{pdb_id}" if pdb_id else "",
                "sequence": seq,
            })
    return rows


def collect_fasta_rows(project: Path, min_length: int, exclude_pdb: set[str]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    files = [project / "Metal-F_protein_entities.fasta", *project.glob("Metal-*/*_hmmer_full_length_seed.fasta")]
    for path in files:
        set_name = "Metal-F" if path.parent == project else path.parent.name.replace("Metal-", "")
        for header, seq in parse_fasta(path):
            if len(seq) < min_length:
                continue
            meta = parse_header_meta(header)
            pdb_id = meta.get("pdb_id", "").upper()
            if pdb_id in exclude_pdb:
                continue
            uid = primary_uniprot(meta.get("uniprot_ids", ""))
            seq_id = uid or meta.get("header_id") or hashlib.sha1(seq.encode()).hexdigest()[:12]
            rows.append({
                "seed_id": safe_id(seq_id),
                "set_name": set_name,
                "source_file": str(path),
                "pdb_id": pdb_id,
                "entity_id": meta.get("entity_id", ""),
                "chains": meta.get("chains", ""),
                "sequence_length": str(len(seq)),
                "uniprot_ids": meta.get("uniprot_ids", uid),
                "organism": meta.get("organism", ""),
                "taxonomy_id": "",
                "kingdom": infer_kingdom(meta.get("organism", "")),
                "doi": "",
                "title": "",
                "pdb_url": f"https://www.rcsb.org/structure/{meta.get('pdb_id', '')}" if meta.get("pdb_id") else "",
                "sequence": seq,
            })
    return rows


def dedupe_rows(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    best: Dict[str, Dict[str, str]] = {}
    for row in rows:
        uid = primary_uniprot(row.get("uniprot_ids", ""))
        key = f"uid:{uid}" if uid else f"seq:{hashlib.sha1(row['sequence'].encode()).hexdigest()}"
        old = best.get(key)
        if old is None or int(row["sequence_length"]) > int(old["sequence_length"]):
            best[key] = row
        elif old is not None:
            old_sets = set(filter(None, old.get("set_name", "").split(";")))
            old_sets.add(row.get("set_name", ""))
            old["set_name"] = ";".join(sorted(old_sets))
    out = sorted(best.values(), key=lambda r: (r.get("kingdom", ""), r.get("seed_id", "")))
    used = set()
    for row in out:
        base = row["seed_id"]
        seq_id = base
        suffix = 2
        while seq_id in used:
            seq_id = f"{base}_dup{suffix}"
            suffix += 1
        used.add(seq_id)
        row["seed_id"] = seq_id
    return out


def write_fasta(path: Path, rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            desc = (
                f"set={row.get('set_name','')} pdb={row.get('pdb_id','')} "
                f"entity={row.get('entity_id','')} chains={row.get('chains','')} "
                f"uniprot={row.get('uniprot_ids','')} kingdom={row.get('kingdom','Unknown')} "
                f"organism={row.get('organism','')} length={row.get('sequence_length','')}"
            )
            handle.write(f">{row['seed_id']} {desc}\n")
            for line in wrap(row["sequence"]):
                handle.write(line + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a broad Metal-F tree seed FASTA from local project records.")
    parser.add_argument("--metal-project", default="/home/qin/Metal-F_project")
    parser.add_argument("--outdir", default="/mnt/e/Tree-Metal-F")
    parser.add_argument("--min-length", type=int, default=120)
    parser.add_argument("--exclude-pdb", nargs="*", default=[])
    args = parser.parse_args()

    project = Path(args.metal_project)
    outdir = Path(args.outdir)
    exclude_pdb = {p.upper() for p in args.exclude_pdb}
    rows = collect_csv_rows(project, args.min_length, exclude_pdb) + collect_fasta_rows(project, args.min_length, exclude_pdb)
    deduped = dedupe_rows(rows)
    fields = [
        "seed_id", "set_name", "source_file", "pdb_id", "entity_id", "chains",
        "sequence_length", "uniprot_ids", "organism", "taxonomy_id", "kingdom",
        "doi", "title", "pdb_url", "sequence",
    ]
    write_fasta(outdir / "tree_core.fasta", deduped)
    write_csv(outdir / "tree_core_metadata.csv", deduped, fields)
    write_csv(outdir / "tree_core_all_source_rows.csv", rows, fields)
    print(f"Collected source rows: {len(rows)}")
    print(f"Wrote {outdir / 'tree_core.fasta'} with {len(deduped)} deduplicated records.")
    print(f"Wrote {outdir / 'tree_core_metadata.csv'}")


if __name__ == "__main__":
    main()
