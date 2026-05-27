#!/usr/bin/env python3
"""
Prepare clean tree representative FASTA for SoluProt and document next steps.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple


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
                    records.append((header, re.sub(r"[^A-Za-z]", "", "".join(seq)).upper()))
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if header and seq:
        records.append((header, re.sub(r"[^A-Za-z]", "", "".join(seq)).upper()))
    return records


def write_fasta(records: List[Tuple[str, str]], path: Path, simplify: bool) -> List[Tuple[str, str]]:
    mapping: List[Tuple[str, str]] = []
    used = set()
    with path.open("w", encoding="utf-8") as handle:
        for old, seq in records:
            new = simplify_id(old) if simplify else old
            base = new
            suffix = 2
            while new in used:
                new = f"{base}_dup{suffix}"
                suffix += 1
            used.add(new)
            mapping.append((old, new))
            handle.write(f">{new}\n")
            for idx in range(0, len(seq), 80):
                handle.write(seq[idx:idx + 80] + "\n")
    return mapping


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare representative FASTA for SoluProt.")
    parser.add_argument("representatives_fasta")
    parser.add_argument("--outdir", default="soluprot_inputs")
    parser.add_argument("--keep-original-ids", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    records = read_fasta(Path(args.representatives_fasta))
    mapping = write_fasta(records, outdir / "representatives_for_soluprot.fasta", simplify=not args.keep_original_ids)
    with (outdir / "soluprot_id_mapping.tsv").open("w", encoding="utf-8") as handle:
        handle.write("old_id\tclean_id\n")
        for old, new in mapping:
            handle.write(f"{old}\t{new}\n")
    (outdir / "README_soluprot_next_steps.txt").write_text(
        "Run SoluProt on representatives_for_soluprot.fasta, then feed its CSV output to:\n\n"
        "python /mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/annotate_solubility_itol.py \\\n"
        "  soluprot_predictions.csv \\\n"
        "  --nodes /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/ssn_nodes.csv \\\n"
        "  --outdir /mnt/e/alpha_beta_hydrolase/hydrolase_repro_ssn/tree_analysis \\\n"
        "  --tool-label SoluProt\n\n"
        "Upload itol_solubility_simplebar.txt and itol_solubility_gradient_strip.txt to iTOL.\n",
        encoding="utf-8",
    )
    print(f"Wrote {outdir / 'representatives_for_soluprot.fasta'} with {len(records)} records.")
    print(f"Wrote {outdir / 'soluprot_id_mapping.tsv'}")


if __name__ == "__main__":
    main()
