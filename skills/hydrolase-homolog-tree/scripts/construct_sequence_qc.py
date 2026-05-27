#!/usr/bin/env python3
"""QC FASTA sequences before plasmid design or expression screening."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


PROTEIN_STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
PROTEIN_AMBIGUOUS = set("XBZJUO*")
DNA_CANONICAL = set("ACGT")
DNA_IUPAC = set("ACGTRYSWKMBDHVN")
HYDROPHOBIC = set("AILMFWYVC")
SMALL = set("AVGSTC")


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    current_id = ""
    chunks: List[str] = []
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    records.append((current_id, "".join(chunks)))
                current_id = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(re.sub(r"\s+", "", line))
    if current_id:
        records.append((current_id, "".join(chunks)))
    return records


def guess_sequence_type(seq: str) -> str:
    letters = [ch for ch in seq.upper() if ch.isalpha()]
    if not letters:
        return "protein"
    dna_like = sum(1 for ch in letters if ch in DNA_IUPAC)
    return "dna" if dna_like / len(letters) >= 0.9 else "protein"


def run_lengths(seq: str, allowed: Iterable[str]) -> Tuple[int, int, int]:
    allowed_set = set(allowed)
    best_len = best_start = 0
    current_start = current_len = 0
    for idx, ch in enumerate(seq, start=1):
        if ch in allowed_set:
            if current_len == 0:
                current_start = idx
            current_len += 1
            if current_len > best_len:
                best_len = current_len
                best_start = current_start
        else:
            current_len = 0
    return best_len, best_start, best_start + best_len - 1 if best_len else 0


def cleavage_candidates(seq: str) -> List[int]:
    out: List[int] = []
    # Signal peptide cleavage often follows the (-3, -1) small-residue rule.
    for cut in range(16, min(len(seq), 35) + 1):
        if seq[cut - 3 - 1] in SMALL and seq[cut - 1 - 1] in SMALL and seq[cut - 2 - 1] != "P":
            out.append(cut)
    return out


def protein_qc(seq_id: str, seq: str) -> Dict[str, str]:
    clean = re.sub(r"[^A-Za-z*]", "", seq).upper()
    invalid = sorted(set(clean) - PROTEIN_STANDARD - PROTEIN_AMBIGUOUS)
    ambiguous = sorted(set(clean) & PROTEIN_AMBIGUOUS)
    nterm = clean[:40]
    cterm = clean[-40:]
    n_run, n_start, n_end = run_lengths(nterm, HYDROPHOBIC)
    c_run, c_start, c_end = run_lengths(cterm, HYDROPHOBIC)
    charges = sum(1 for ch in clean[:8] if ch in {"K", "R"})
    cleavages = cleavage_candidates(clean)

    signal_like = "no"
    if len(clean) >= 45 and n_run >= 8 and (charges >= 1 or n_run >= 12):
        signal_like = "yes"
    terminal_anchor = "yes" if c_run >= 15 else "no"

    issues: List[str] = []
    if "X" in ambiguous:
        issues.append("contains_X_unknown_residue")
    other_amb = [ch for ch in ambiguous if ch != "X"]
    if other_amb:
        issues.append("contains_nonstandard_residue:" + "".join(other_amb))
    if invalid:
        issues.append("contains_invalid_residue:" + "".join(invalid))
    if "*" in clean[:-1]:
        issues.append("internal_stop_symbol")
    if signal_like == "yes":
        issues.append("n_terminal_signal_peptide_like")
    if terminal_anchor == "yes":
        issues.append("c_terminal_hydrophobic_anchor_like")

    return {
        "id": seq_id,
        "sequence_type": "protein",
        "length": str(len(clean)),
        "bad_characters": "".join(invalid),
        "ambiguous_characters": "".join(ambiguous),
        "x_count": str(clean.count("X")),
        "internal_stop_count": str(clean[:-1].count("*")),
        "nterm_positive_charges_1_8": str(charges),
        "nterm_hydrophobic_run_len_1_40": str(n_run),
        "nterm_hydrophobic_run_pos_1_40": f"{n_start}-{n_end}" if n_run else "",
        "signal_peptide_like": signal_like,
        "possible_signalp_cleavage_after": ",".join(map(str, cleavages[:5])),
        "cterm_hydrophobic_run_len_last40": str(c_run),
        "cterm_hydrophobic_run_pos_last40": f"{c_start}-{c_end}" if c_run else "",
        "cterm_anchor_like": terminal_anchor,
        "issue_summary": ";".join(issues) if issues else "ok",
    }


def dna_qc(seq_id: str, seq: str, coding: bool) -> Dict[str, str]:
    clean = re.sub(r"\s+", "", seq).upper()
    letters = re.sub(r"[^A-Z]", "", clean)
    invalid = sorted(set(letters) - DNA_IUPAC)
    ambiguous = sorted(set(letters) & (DNA_IUPAC - DNA_CANONICAL))
    issues: List[str] = []
    if invalid:
        issues.append("contains_invalid_dna_character:" + "".join(invalid))
    if ambiguous:
        issues.append("contains_ambiguous_iupac_base:" + "".join(ambiguous))
    if coding:
        if len(letters) % 3:
            issues.append("coding_length_not_multiple_of_3")
        codons = [letters[i : i + 3] for i in range(0, len(letters) - 2, 3)]
        stops = [i + 1 for i, codon in enumerate(codons[:-1]) if codon in {"TAA", "TAG", "TGA"}]
        if stops:
            issues.append("internal_stop_codon")
    return {
        "id": seq_id,
        "sequence_type": "dna",
        "length": str(len(letters)),
        "bad_characters": "".join(invalid),
        "ambiguous_characters": "".join(ambiguous),
        "x_count": str(letters.count("X")),
        "internal_stop_count": "" if not coding else str(len([c for c in re.findall(r"...", letters[:-3]) if c in {"TAA", "TAG", "TGA"}])),
        "nterm_positive_charges_1_8": "",
        "nterm_hydrophobic_run_len_1_40": "",
        "nterm_hydrophobic_run_pos_1_40": "",
        "signal_peptide_like": "",
        "possible_signalp_cleavage_after": "",
        "cterm_hydrophobic_run_len_last40": "",
        "cterm_hydrophobic_run_pos_last40": "",
        "cterm_anchor_like": "",
        "issue_summary": ";".join(issues) if issues else "ok",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Check FASTA records for construct-design risk flags.")
    parser.add_argument("fasta")
    parser.add_argument("--output", default=None, help="Output CSV. Defaults to <fasta>.construct_qc.csv.")
    parser.add_argument("--sequence-type", choices=["auto", "protein", "dna"], default="auto")
    parser.add_argument("--dna-coding", action="store_true", help="Treat DNA records as coding sequences and check frame/stops.")
    args = parser.parse_args()

    fasta = Path(args.fasta)
    outpath = Path(args.output) if args.output else fasta.with_suffix(fasta.suffix + ".construct_qc.csv")
    records = read_fasta(fasta)
    rows: List[Dict[str, str]] = []
    for seq_id, seq in records:
        seq_type = guess_sequence_type(seq) if args.sequence_type == "auto" else args.sequence_type
        if seq_type == "dna":
            rows.append(dna_qc(seq_id, seq, args.dna_coding))
        else:
            rows.append(protein_qc(seq_id, seq))

    headers = [
        "id",
        "sequence_type",
        "length",
        "bad_characters",
        "ambiguous_characters",
        "x_count",
        "internal_stop_count",
        "nterm_positive_charges_1_8",
        "nterm_hydrophobic_run_len_1_40",
        "nterm_hydrophobic_run_pos_1_40",
        "signal_peptide_like",
        "possible_signalp_cleavage_after",
        "cterm_hydrophobic_run_len_last40",
        "cterm_hydrophobic_run_pos_last40",
        "cterm_anchor_like",
        "issue_summary",
    ]
    with outpath.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)

    flagged = [row for row in rows if row["issue_summary"] != "ok"]
    print(f"Wrote {outpath}")
    print(f"Records: {len(rows)}")
    print(f"Flagged: {len(flagged)}")
    if flagged:
        for row in flagged[:20]:
            print(f"{row['id']}\t{row['issue_summary']}")


if __name__ == "__main__":
    main()
