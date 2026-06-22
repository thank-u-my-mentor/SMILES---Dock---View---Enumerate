#!/usr/bin/env python3
"""Apply a PLE protonation_map.csv residue renaming table to a PDB."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--map", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    mapping: dict[tuple[str, str], str] = {}
    with args.map.open(newline="", encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle):
            chain = (row.get("chain") or "").strip()
            resseq = (row.get("resseq") or row.get("resid") or "").strip()
            resname = (row.get("resname") or row.get("amber_resname") or "").strip().upper()
            if chain and resseq and resname:
                mapping[(chain, resseq)] = resname

    out_lines = []
    renamed: dict[tuple[str, str, str, str], int] = {}
    for line in args.pdb.read_text(errors="replace").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            chain = line[21:22].strip()
            resseq = line[22:26].strip()
            new = mapping.get((chain, resseq))
            if new:
                old = line[17:20].strip()
                if old != new:
                    renamed[(chain, resseq, old, new)] = renamed.get((chain, resseq, old, new), 0) + 1
                    line = f"{line[:17]}{new:>3}{line[20:]}"
        out_lines.append(line)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    for (chain, resseq, old, new), count in sorted(renamed.items()):
        print(f"renamed {old} {chain}{resseq} -> {new} atoms={count}")
    print(f"out={args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
