#!/usr/bin/env python3
"""Audit a ligand PDB geometry against the ligand's mol2 bond table."""

from __future__ import annotations

import argparse
import math
from pathlib import Path


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def parse_mol2(path: Path) -> tuple[dict[int, str], list[tuple[str, str, str]]]:
    atoms: dict[int, str] = {}
    bonds: list[tuple[str, str, str]] = []
    section = ""
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("@<TRIPOS>"):
            section = line.strip()
            continue
        parts = line.split()
        if not parts:
            continue
        if section == "@<TRIPOS>ATOM" and len(parts) >= 6:
            atoms[int(parts[0])] = parts[1]
        elif section == "@<TRIPOS>BOND" and len(parts) >= 4:
            a = atoms.get(int(parts[1]))
            b = atoms.get(int(parts[2]))
            if a is None or b is None:
                raise SystemExit(f"Bond references unknown atom ids: {line}")
            bonds.append((a, b, parts[3]))
    if not atoms or not bonds:
        raise SystemExit(f"No mol2 atoms/bonds parsed from {path}")
    return atoms, bonds


def parse_pdb_residue(path: Path, resname: str) -> dict[str, tuple[float, float, float]]:
    coords: dict[str, tuple[float, float, float]] = {}
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip() != resname:
            continue
        atom_name = line[12:16].strip()
        coords[atom_name] = tuple(float(line[index : index + 8]) for index in (30, 38, 46))
    if not coords:
        raise SystemExit(f"No residue {resname!r} found in {path}")
    return coords


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mol2", required=True)
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--resname", required=True)
    parser.add_argument("--reference-pdb")
    parser.add_argument("--out")
    parser.add_argument("--long-bond-cutoff", type=float, default=2.2)
    parser.add_argument("--short-bond-cutoff", type=float, default=0.85)
    parser.add_argument("--large-change-cutoff", type=float, default=0.35)
    args = parser.parse_args()

    _, bonds = parse_mol2(Path(args.mol2))
    coords = parse_pdb_residue(Path(args.pdb), args.resname)
    ref_coords = parse_pdb_residue(Path(args.reference_pdb), args.resname) if args.reference_pdb else None

    rows = ["bond\torder\tpdb_A\treference_A\tdelta_A\tflag"]
    problem_count = 0
    for atom_a, atom_b, order in bonds:
        if atom_a not in coords or atom_b not in coords:
            rows.append(f"{atom_a}-{atom_b}\t{order}\tNA\tNA\tNA\tMISSING_PDB_ATOM")
            problem_count += 1
            continue
        observed = distance(coords[atom_a], coords[atom_b])
        ref_value = None
        delta = None
        if ref_coords is not None and atom_a in ref_coords and atom_b in ref_coords:
            ref_value = distance(ref_coords[atom_a], ref_coords[atom_b])
            delta = observed - ref_value
        flags: list[str] = []
        if observed > args.long_bond_cutoff:
            flags.append("BROKEN_LONG_BOND")
        if observed < args.short_bond_cutoff:
            flags.append("TOO_SHORT")
        if delta is not None and abs(delta) > args.large_change_cutoff:
            flags.append("LARGE_CHANGE")
        if flags:
            problem_count += 1
        ref_text = f"{ref_value:.4f}" if ref_value is not None else "NA"
        delta_text = f"{delta:+.4f}" if delta is not None else "NA"
        rows.append(
            f"{atom_a}-{atom_b}\t{order}\t{observed:.4f}\t"
            f"{ref_text}\t{delta_text}\t{','.join(flags)}"
        )

    text = "\n".join(rows) + "\n"
    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    print(f"problem_count={problem_count}")
    return 2 if problem_count else 0


if __name__ == "__main__":
    raise SystemExit(main())
