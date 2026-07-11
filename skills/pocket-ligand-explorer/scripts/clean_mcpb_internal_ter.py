#!/usr/bin/env python3
"""Remove spurious internal protein TER records from MCPB-generated PDB files.

MCPB.py may write TER records between many protein fragments after residue
truncation/reassembly.  tleap interprets each TER as a real peptide terminus and
adds OXT/H1/H2/H3, which can overlap the next residue and explode GROMACS energy
minimization.  This script removes TER records between consecutive protein ATOM
records in the same chain when the residue numbers are close.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_atom(line: str) -> tuple[str, str, int] | None:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    chain = line[21].strip()
    resname = line[17:20].strip()
    try:
        resseq = int(line[22:26])
    except ValueError:
        return None
    return chain, resname, resseq


def is_protein_atom(line: str) -> bool:
    return line.startswith("ATOM")


def clean_lines(lines: list[str], max_gap: int) -> tuple[list[str], int, int]:
    kept: list[str] = []
    removed = 0
    kept_ter = 0
    n = len(lines)
    for i, line in enumerate(lines):
        if not line.startswith("TER"):
            kept.append(line)
            continue

        prev_line = None
        for j in range(len(kept) - 1, -1, -1):
            if kept[j].startswith(("ATOM", "HETATM")):
                prev_line = kept[j]
                break
        next_line = None
        for j in range(i + 1, n):
            if lines[j].startswith(("ATOM", "HETATM")):
                next_line = lines[j]
                break
            if lines[j].startswith("END"):
                break

        prev_info = parse_atom(prev_line) if prev_line else None
        next_info = parse_atom(next_line) if next_line else None
        if (
            prev_line
            and next_line
            and is_protein_atom(prev_line)
            and is_protein_atom(next_line)
            and prev_info
            and next_info
            and prev_info[0] == next_info[0]
            and abs(next_info[2] - prev_info[2]) <= max_gap
        ):
            removed += 1
            continue

        kept.append(line)
        kept_ter += 1
    return kept, removed, kept_ter


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-pdb", required=True)
    parser.add_argument("--out-pdb", required=True)
    parser.add_argument(
        "--max-gap",
        type=int,
        default=2,
        help="Remove TER between same-chain protein residues whose residue ids differ by at most this value.",
    )
    args = parser.parse_args()

    in_pdb = Path(args.in_pdb)
    out_pdb = Path(args.out_pdb)
    lines = in_pdb.read_text().splitlines(keepends=True)
    cleaned, removed, kept_ter = clean_lines(lines, args.max_gap)
    out_pdb.write_text("".join(cleaned))
    print(f"input={in_pdb}")
    print(f"output={out_pdb}")
    print(f"removed_internal_protein_TER={removed}")
    print(f"kept_TER={kept_ter}")


if __name__ == "__main__":
    main()
