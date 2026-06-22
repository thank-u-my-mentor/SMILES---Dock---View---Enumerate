#!/usr/bin/env python
"""Write a small GROMACS index for protein/ligand/metal visualization."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gro", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    return parser.parse_args()


def chunks(values: list[int], size: int = 15):
    for i in range(0, len(values), size):
        yield values[i : i + size]


def write_group(lines: list[str], name: str, atoms: list[int]) -> None:
    lines.append(f"[ {name} ]")
    for chunk in chunks(atoms):
        lines.append(" ".join(f"{value:6d}" for value in chunk))
    lines.append("")


def main() -> None:
    args = parse_args()
    gro_lines = args.gro.read_text(errors="replace").splitlines()
    atom_lines = gro_lines[2:-1]
    protein = []
    ligand = []
    fe = []
    water_ions = []
    for line in atom_lines:
        atomnr = int(line[15:20])
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        if resname == "LIG":
            ligand.append(atomnr)
        elif resname == "FE" or atomname.upper() == "FE":
            fe.append(atomnr)
        elif resname in {"WAT", "SOL", "NA", "CL", "Na+", "Cl-"}:
            water_ions.append(atomnr)
        else:
            protein.append(atomnr)

    lines: list[str] = []
    write_group(lines, "Protein", protein)
    write_group(lines, "LIG", ligand)
    write_group(lines, "FE", fe)
    write_group(lines, "Protein_LIG_FE", sorted(protein + ligand + fe))
    write_group(lines, "Water_Ions", water_ions)
    args.out.write_text("\n".join(lines), encoding="utf-8")
    print(args.out)
    print(f"protein={len(protein)} ligand={len(ligand)} fe={len(fe)} water_ions={len(water_ions)}")


if __name__ == "__main__":
    main()
