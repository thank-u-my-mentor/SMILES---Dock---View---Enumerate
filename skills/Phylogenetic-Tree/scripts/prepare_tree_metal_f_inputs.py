#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List


KINGDOM_KEYWORDS = [
    (
        "Bacteria",
        [
            "bacter",
            "escherichia",
            "pseudomonas",
            "staphylococcus",
            "xanthomonas",
            "bacillus",
            "streptomyces",
            "clostridium",
            "salmonella",
            "klebsiella",
            "mycobacterium",
            "lactobacillus",
        ],
    ),
    ("Archaea", ["archaea", "archaeon", "methano", "halobacter"]),
    ("Fungi", ["fung", "saccharomyces", "aspergillus", "candida", "neurospora", "fusarium", "penicillium"]),
    ("Plant", ["viridiplantae", "plantae", "arabidopsis", "oryza sativa", "zea mays"]),
    ("Animal", ["homo sapiens", "rattus", "mus musculus", "drosophila", "danio", "xenopus", "bos taurus"]),
    ("Virus", ["virus", "phage"]),
]


def infer_kingdom(organism: str) -> str:
    text = (organism or "").lower()
    for kingdom, keys in KINGDOM_KEYWORDS:
        if any(key in text for key in keys):
            return kingdom
    return "Unknown"


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: List[Dict[str, str]], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def wrap(seq: str, width: int = 80) -> Iterable[str]:
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


def parse_int(value: str) -> int:
    try:
        return int(float(str(value)))
    except Exception:
        return 0


def safe_id(row: Dict[str, str]) -> str:
    uid = (row.get("uniprot_ids") or row.get("uniprot_id") or "").split(";")[0].strip()
    if uid:
        return uid
    pdb = row.get("pdb_id", "PDB")
    entity = row.get("entity_id", "1")
    return f"{pdb}_{entity}"


def collect_seed_rows(metal_project: Path, sets: List[str], min_length: int) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for name in sets:
        d = metal_project / f"Metal-{name}"
        kept = d / f"Metal-{name}_protein_entities_kept.csv"
        geom = d / f"Metal-{name}_geometry_summary.csv"
        if not kept.exists():
            continue
        geom_by_pdb = {row.get("pdb_id", ""): row for row in read_csv(geom)} if geom.exists() else {}
        for row in read_csv(kept):
            length = parse_int(row.get("sequence_length", "0"))
            if length < min_length:
                continue
            pdb_id = row.get("pdb_id", "")
            geom_row = geom_by_pdb.get(pdb_id, {})
            organism = row.get("organism", "")
            out = dict(row)
            out["set_name"] = name
            out["seed_id"] = safe_id(row)
            out["kingdom"] = infer_kingdom(organism)
            out["nearest_ligandF_metal_A"] = geom_row.get("nearest_ligandF_metal_A", "")
            out["doi"] = geom_row.get("doi", "")
            out["pdb_url"] = geom_row.get("pdb_url", f"https://www.rcsb.org/structure/{pdb_id}" if pdb_id else "")
            rows.append(out)
    return rows


def dedupe_rows(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    seen = set()
    out = []
    for row in rows:
        key = (row.get("uniprot_ids", ""), row.get("sequence", ""))
        if key in seen:
            continue
        seen.add(key)
        out.append(row)
    return out


def write_fasta(path: Path, rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    used = set()
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            seq_id = row["seed_id"]
            base = seq_id
            suffix = 2
            while seq_id in used:
                seq_id = f"{base}_dup{suffix}"
                suffix += 1
            used.add(seq_id)
            desc = (
                f"pdb={row.get('pdb_id','')} entity={row.get('entity_id','')} "
                f"set={row.get('set_name','')} kingdom={row.get('kingdom','Unknown')} "
                f"organism={row.get('organism','')}"
            )
            handle.write(f">{seq_id} {desc}\n")
            for line in wrap(re.sub(r"[^A-Za-z]", "", row.get("sequence", "")).upper()):
                handle.write(line + "\n")


def write_kingdom_itol(path: Path, rows: List[Dict[str, str]]) -> None:
    colors = {
        "Plant": "#A8D5BA",
        "Animal": "#8FB8E6",
        "Fungi": "#59A14F",
        "Bacteria": "#E8A0B0",
        "Archaea": "#F4C2A1",
        "Protist": "#F9E79F",
        "Virus": "#E15759",
        "Metagenome": "#D5D5D5",
        "Eukaryota": "#EDC948",
        "Unknown": "#AAAAAA",
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("DATASET_COLORSTRIP\nSEPARATOR TAB\n")
        handle.write("DATASET_LABEL\tKingdom\nCOLOR\t#000000\nSTRIP_WIDTH\t25\nMARGIN\t5\n")
        handle.write("LEGEND_TITLE\tKingdom\n")
        handle.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(colors)) + "\n")
        handle.write("LEGEND_COLORS\t" + "\t".join(colors.values()) + "\n")
        handle.write("LEGEND_LABELS\t" + "\t".join(colors.keys()) + "\n")
        handle.write("DATA\n")
        for row in rows:
            kingdom = row.get("kingdom") or "Unknown"
            handle.write(f"{row['seed_id']}\t{colors.get(kingdom, colors['Unknown'])}\t{kingdom}\n")


def write_soluprot_inputs(outdir: Path, rows: List[Dict[str, str]]) -> None:
    sol_dir = outdir / "soluprot_inputs"
    sol_dir.mkdir(parents=True, exist_ok=True)
    write_fasta(sol_dir / "representatives_for_soluprot.fasta", rows)
    with (sol_dir / "soluprot_id_mapping.tsv").open("w", encoding="utf-8") as handle:
        handle.write("old_id\tclean_id\n")
        for row in rows:
            handle.write(f"{row['seed_id']}\t{row['seed_id']}\n")
    (sol_dir / "README_soluprot_next_steps.txt").write_text(
        "Upload representatives_for_soluprot.fasta to SoluProt/NetSolP manually.\n"
        "Save the downloaded prediction CSV as input.csv in this directory.\n"
        "Then run:\n\n"
        "python /mnt/e/Codex/skills/Phylogenetic-Tree/scripts/make_soluprot_annotation.py "
        "--input /mnt/e/Tree-Metal-F/soluprot_inputs/input.csv "
        "--outdir /mnt/e/Tree-Metal-F/tree_analysis\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metal-project", default="/home/qin/Metal-F_project")
    parser.add_argument("--outdir", default="/mnt/e/Tree-Metal-F")
    parser.add_argument("--sets", nargs="+", default=["CF3", "C-F", "fluoroaryl"])
    parser.add_argument("--min-length", type=int, default=150)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = dedupe_rows(collect_seed_rows(Path(args.metal_project), args.sets, args.min_length))
    fields = [
        "seed_id", "set_name", "pdb_id", "entity_id", "chains", "sequence_length",
        "uniprot_ids", "organism", "kingdom", "nearest_ligandF_metal_A", "doi",
        "pdb_url", "sequence",
    ]
    write_csv(outdir / "seed_metadata.csv", rows, fields)
    write_fasta(outdir / "core.fasta", rows)
    write_kingdom_itol(outdir / "itol_seed_kingdom_colorstrip.txt", rows)
    write_soluprot_inputs(outdir, rows)
    (outdir / "README_next_steps.txt").write_text(
        "Core inputs prepared for a general phylogenetic-tree workflow.\n"
        "Use core.fasta as curated seeds. Add HMMER homologs to all.fasta before SSN/tree.\n"
        "Kingdom labels are lightweight, inferred from existing organism metadata only.\n",
        encoding="utf-8",
    )
    print(f"Wrote {outdir / 'core.fasta'} with {len(rows)} deduplicated seed records.")
    print(f"Wrote {outdir / 'seed_metadata.csv'}")
    print(f"Wrote {outdir / 'soluprot_inputs/representatives_for_soluprot.fasta'}")


if __name__ == "__main__":
    main()
