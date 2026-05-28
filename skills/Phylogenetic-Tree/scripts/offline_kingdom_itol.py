#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List


KINGDOM_COLORS = {
    "Bacteria": "#4E79A7",
    "Archaea": "#B07AA1",
    "Fungi": "#59A14F",
    "Plant": "#8CD17D",
    "Animal": "#F28E2B",
    "Virus": "#E15759",
    "Metagenome": "#9C755F",
    "Eukaryota": "#EDC948",
    "Unknown": "#BDBDBD",
}

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


def read_fasta_headers(path: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    if not path.exists():
        return rows
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.startswith(">"):
                continue
            header = raw[1:].strip()
            seq_id = header.split()[0]
            rows.append({"id": seq_id, "header": header})
    return rows


def load_seed_metadata(path: Path) -> Dict[str, Dict[str, str]]:
    metadata: Dict[str, Dict[str, str]] = {}
    for row in read_csv(path):
        ids = {
            row.get("seed_id", ""),
            row.get("uniprot_ids", "").split(";")[0].strip(),
            row.get("pdb_id", ""),
        }
        for key in ids:
            if key:
                metadata[key] = row
    return metadata


def row_id(row: Dict[str, str]) -> str:
    return row.get("id") or row.get("seed_id") or row.get("uniprot_id") or row.get("accession") or row.get("Entry") or ""


def annotate_rows(
    nodes: Iterable[Dict[str, str]],
    seed_metadata: Dict[str, Dict[str, str]],
) -> List[Dict[str, str]]:
    out: List[Dict[str, str]] = []
    for row in nodes:
        node_id = row_id(row)
        organism = row.get("organism") or row.get("Organism") or row.get("header") or ""
        kingdom = row.get("kingdom") or row.get("Kingdom") or ""
        if not kingdom or kingdom == "Unknown":
            md = seed_metadata.get(node_id) or seed_metadata.get(row.get("uniprot_id", ""))
            if md:
                organism = organism or md.get("organism", "")
                kingdom = md.get("kingdom", "") or infer_kingdom(organism)
        if not kingdom or kingdom == "Unknown":
            kingdom = infer_kingdom(organism)
        out_row = dict(row)
        out_row["id"] = node_id
        out_row["organism"] = organism or row.get("organism", "")
        out_row["kingdom"] = kingdom or "Unknown"
        out.append(out_row)
    return out


def write_csv(path: Path, rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row.keys()})
    preferred = ["id", "uniprot_id", "kingdom", "organism", "is_core", "length"]
    fields = [f for f in preferred if f in fields] + [f for f in fields if f not in preferred]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_itol(path: Path, rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    colors = KINGDOM_COLORS
    with path.open("w", encoding="utf-8") as handle:
        handle.write("DATASET_COLORSTRIP\nSEPARATOR TAB\n")
        handle.write("DATASET_LABEL\tKingdom_offline\nCOLOR\t#000000\nSTRIP_WIDTH\t25\nMARGIN\t5\n")
        handle.write("LEGEND_TITLE\tKingdom\n")
        handle.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(colors)) + "\n")
        handle.write("LEGEND_COLORS\t" + "\t".join(colors.values()) + "\n")
        handle.write("LEGEND_LABELS\t" + "\t".join(colors.keys()) + "\n")
        handle.write("DATA\n")
        for row in rows:
            node_id = row.get("id", "")
            kingdom = row.get("kingdom") or "Unknown"
            if not node_id:
                continue
            handle.write(f"{node_id}\t{colors.get(kingdom, colors['Unknown'])}\t{kingdom}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Offline kingdom annotation for tree representatives.")
    parser.add_argument("--nodes", default="/mnt/e/Tree-Metal-F/ssn/ssn_nodes.csv")
    parser.add_argument("--fasta", default="/mnt/e/Tree-Metal-F/ssn/representatives.fasta")
    parser.add_argument("--seed-metadata", default="/mnt/e/Tree-Metal-F/seed_metadata.csv")
    parser.add_argument("--outdir", default="/mnt/e/Tree-Metal-F/tree_analysis")
    args = parser.parse_args()

    node_rows = read_csv(Path(args.nodes))
    if not node_rows:
        node_rows = read_fasta_headers(Path(args.fasta))
    seed_metadata = load_seed_metadata(Path(args.seed_metadata))
    rows = annotate_rows(node_rows, seed_metadata)
    outdir = Path(args.outdir)
    write_csv(outdir / "nodes_with_offline_kingdom.csv", rows)
    write_itol(outdir / "itol_kingdom_color_strip_offline.txt", rows)
    unknown = sum(1 for row in rows if row.get("kingdom") == "Unknown")
    print(f"Wrote {outdir / 'nodes_with_offline_kingdom.csv'}")
    print(f"Wrote {outdir / 'itol_kingdom_color_strip_offline.txt'}")
    print(f"Rows: {len(rows)}; Unknown kingdom: {unknown}")


if __name__ == "__main__":
    main()
