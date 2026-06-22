#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List


KINGDOM_COLORS = {
    "Plant": "#A8D5BA",
    "Animal": "#8FB8E6",
    "Fungi": "#59A14F",
    "Bacteria": "#E8A0B0",
    "Archaea": "#F4C2A1",
    "Protist": "#F9E79F",
    "Metagenome": "#D5D5D5",
    "Eukaryota": "#EDC948",
    "Unknown": "#AAAAAA",
}

DEFAULT_KINGDOM_STYLE = {
    "strip_width": "25",
    "margin": "5",
    "border_width": "1",
    "border_color": "#000000",
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


def normalize_kingdom(kingdom: str, organism: str = "", phylum: str = "") -> str:
    text = " ".join([kingdom or "", organism or "", phylum or ""]).lower()
    if any(key in text for key in ["streptophyta", "viridiplantae", "plantae", "embryophyta", "tracheophyta"]):
        return "Plant"
    if any(key in text for key in ["fungi", "ascomycota", "basidiomycota", "saccharomyces", "fusarium", "steccherinum"]):
        return "Fungi"
    if any(key in text for key in ["metazoa", "animalia", "chordata", "arthropoda", "mammalia", "homo sapiens", "rattus"]):
        return "Animal"
    if "bacteria" in text:
        return "Bacteria"
    if "archaea" in text:
        return "Archaea"
    if "virus" in text or "viruses" in text:
        return "Virus"
    if kingdom and kingdom not in {"Unknown", "Eukaryota"}:
        return kingdom
    inferred = infer_kingdom(" ".join([organism or "", phylum or ""]))
    return inferred


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
        phylum = row.get("phylum") or row.get("phylum_class") or row.get("Phylum") or ""
        if not kingdom or kingdom == "Unknown":
            md = seed_metadata.get(node_id) or seed_metadata.get(row.get("uniprot_id", ""))
            if md:
                organism = organism or md.get("organism", "")
                kingdom = md.get("kingdom", "") or infer_kingdom(organism)
                phylum = phylum or md.get("phylum", "") or md.get("phylum_class", "")
        kingdom = normalize_kingdom(kingdom, organism, phylum)
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


def write_itol(
    path: Path,
    rows: List[Dict[str, str]],
    *,
    strip_width: str = DEFAULT_KINGDOM_STYLE["strip_width"],
    margin: str = DEFAULT_KINGDOM_STYLE["margin"],
    border_width: str = DEFAULT_KINGDOM_STYLE["border_width"],
    border_color: str = DEFAULT_KINGDOM_STYLE["border_color"],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    colors = KINGDOM_COLORS
    with path.open("w", encoding="utf-8") as handle:
        handle.write("DATASET_COLORSTRIP\nSEPARATOR TAB\n")
        handle.write("DATASET_LABEL\tKingdom_offline\nCOLOR\t#000000\n")
        handle.write(f"STRIP_WIDTH\t{strip_width}\n")
        handle.write(f"MARGIN\t{margin}\n")
        handle.write(f"BORDER_WIDTH\t{border_width}\n")
        handle.write(f"BORDER_COLOR\t{border_color}\n")
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
    parser.add_argument("--kingdom-strip-width", default=DEFAULT_KINGDOM_STYLE["strip_width"])
    parser.add_argument("--kingdom-margin", default=DEFAULT_KINGDOM_STYLE["margin"])
    parser.add_argument("--kingdom-border-width", default=DEFAULT_KINGDOM_STYLE["border_width"])
    parser.add_argument("--kingdom-border-color", default=DEFAULT_KINGDOM_STYLE["border_color"])
    args = parser.parse_args()

    node_rows = read_csv(Path(args.nodes))
    if not node_rows:
        node_rows = read_fasta_headers(Path(args.fasta))
    seed_metadata = load_seed_metadata(Path(args.seed_metadata))
    rows = annotate_rows(node_rows, seed_metadata)
    outdir = Path(args.outdir)
    write_csv(outdir / "nodes_with_offline_kingdom.csv", rows)
    write_itol(
        outdir / "itol_kingdom_color_strip_offline.txt",
        rows,
        strip_width=args.kingdom_strip_width,
        margin=args.kingdom_margin,
        border_width=args.kingdom_border_width,
        border_color=args.kingdom_border_color,
    )
    unknown = sum(1 for row in rows if row.get("kingdom") == "Unknown")
    print(f"Wrote {outdir / 'nodes_with_offline_kingdom.csv'}")
    print(f"Wrote {outdir / 'itol_kingdom_color_strip_offline.txt'}")
    print(f"Rows: {len(rows)}; Unknown kingdom: {unknown}")


if __name__ == "__main__":
    main()
