#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Tuple


def simplify_id(seq_id: str) -> str:
    text = seq_id.strip()
    if "|" in text:
        left, right = text.split("|", 1)
        if right == left or right.startswith(left + "_"):
            return right
    return text


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header = ""
    seq: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header and seq:
                    records.append((header.split()[0], re.sub(r"[^A-Za-z]", "", "".join(seq)).upper()))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header and seq:
        records.append((header.split()[0], re.sub(r"[^A-Za-z]", "", "".join(seq)).upper()))
    return records


def write_fasta(records: List[Tuple[str, str]], path: Path, keep_original_ids: bool) -> List[Tuple[str, str]]:
    mapping: List[Tuple[str, str]] = []
    used = set()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for old_id, seq in records:
            clean_id = old_id if keep_original_ids else simplify_id(old_id)
            base = clean_id
            suffix = 2
            while clean_id in used:
                clean_id = f"{base}_dup{suffix}"
                suffix += 1
            used.add(clean_id)
            mapping.append((old_id, clean_id))
            handle.write(f">{clean_id}\n")
            for idx in range(0, len(seq), 80):
                handle.write(seq[idx:idx + 80] + "\n")
    return mapping


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare tree representative FASTA for manual SoluProt/NetSolP.")
    parser.add_argument(
        "representatives_fasta",
        nargs="?",
        default="/mnt/e/Tree-Metal-F/ssn/representatives.fasta",
    )
    parser.add_argument("--outdir", default="/mnt/e/Tree-Metal-F/soluprot_inputs")
    parser.add_argument("--keep-original-ids", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    records = read_fasta(Path(args.representatives_fasta))
    mapping = write_fasta(records, outdir / "representatives_for_soluprot.fasta", args.keep_original_ids)
    with (outdir / "soluprot_id_mapping.tsv").open("w", encoding="utf-8") as handle:
        handle.write("old_id\tclean_id\n")
        for old_id, clean_id in mapping:
            handle.write(f"{old_id}\t{clean_id}\n")
    (outdir / "README_soluprot_next_steps.txt").write_text(
        "This FASTA should be generated after SSN representative selection and tree construction.\n"
        "Upload representatives_for_soluprot.fasta to SoluProt/NetSolP manually.\n"
        "Save the downloaded prediction CSV as input.csv in this directory.\n"
        "Then run:\n\n"
        "python /mnt/e/Codex/skills/Phylogenetic-Tree/scripts/make_soluprot_annotation.py "
        "--input /mnt/e/Tree-Metal-F/soluprot_inputs/input.csv "
        "--outdir /mnt/e/Tree-Metal-F/tree_analysis "
        "--nodes /mnt/e/Tree-Metal-F/ssn/ssn_nodes.csv "
        "--id-mapping /mnt/e/Tree-Metal-F/soluprot_inputs/soluprot_id_mapping.tsv\n",
        encoding="utf-8",
    )
    print(f"Wrote {outdir / 'representatives_for_soluprot.fasta'} with {len(records)} records.")
    print(f"Wrote {outdir / 'soluprot_id_mapping.tsv'}")


if __name__ == "__main__":
    main()
