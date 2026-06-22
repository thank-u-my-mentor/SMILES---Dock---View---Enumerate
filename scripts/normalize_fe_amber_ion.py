#!/usr/bin/env python3
"""Normalize single-atom iron HETATM records for Amber/tleap ion templates."""

from __future__ import annotations

import argparse
from pathlib import Path


def normalize_line(line: str, residue_name: str, atom_name: str) -> str:
    if not line.startswith(("ATOM  ", "HETATM")):
        return line
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    current_atom = line[12:16].strip().upper()
    current_residue = line[17:20].strip().upper()
    if element != "FE" and current_atom not in {"FE", "FE2", "FE3"} and current_residue not in {"FE", "FE2", "FE3"}:
        return line

    padded = line.rstrip("\n").ljust(80)
    # PDB fixed columns: atom name 13-16, residue name 18-20, element 77-78.
    padded = padded[:12] + f"{atom_name:>4}" + padded[16:17] + f"{residue_name:>3}" + padded[20:76] + "FE" + padded[78:]
    return padded.rstrip() + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--residue-name", default="FE", help="Amber ion residue name; FE is Fe3+, FE2 is Fe2+.")
    parser.add_argument("--atom-name", default="FE", help="Amber ion atom name; FE is Fe3+, FE2 is Fe2+.")
    args = parser.parse_args()

    out = args.out or args.pdb
    changed = 0
    output_lines = []
    for line in args.pdb.read_text().splitlines(True):
        new_line = normalize_line(line, args.residue_name, args.atom_name)
        if new_line != line:
            changed += 1
        output_lines.append(new_line)
    out.write_text("".join(output_lines))
    print(f"pdb={out}")
    print(f"normalized_fe_records={changed}")
    print(f"amber_ion={args.residue_name}/{args.atom_name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
