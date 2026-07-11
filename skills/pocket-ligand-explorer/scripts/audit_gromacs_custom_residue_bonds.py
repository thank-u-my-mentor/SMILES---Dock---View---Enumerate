#!/usr/bin/env python3
"""Audit custom residue covalent bonds in GROMACS gro/top files."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


DEFAULT_BONDS = {
    "HD1": [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB"), ("CB", "CG")],
    "HD2": [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB"), ("CB", "CG")],
    "GU1": [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD")],
}


def parse_gro(path: Path) -> list[dict[str, object]]:
    atoms: list[dict[str, object]] = []
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for line in lines[2:-1]:
        if len(line) < 44:
            continue
        try:
            atoms.append(
                {
                    "resid": int(line[:5]),
                    "resname": line[5:10].strip(),
                    "atom": line[10:15].strip(),
                    "atomnr": int(line[15:20]),
                    "x": float(line[20:28]) * 10.0,
                    "y": float(line[28:36]) * 10.0,
                    "z": float(line[36:44]) * 10.0,
                }
            )
        except ValueError:
            continue
    return atoms


def dist(a: dict[str, object], b: dict[str, object]) -> float:
    return math.sqrt(
        (float(a["x"]) - float(b["x"])) ** 2
        + (float(a["y"]) - float(b["y"])) ** 2
        + (float(a["z"]) - float(b["z"])) ** 2
    )


def parse_top_bonds(path: Path) -> tuple[dict[tuple[int, str, str], int], set[frozenset[int]]]:
    atom_map: dict[tuple[int, str, str], int] = {}
    bonds: set[frozenset[int]] = set()
    section = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(";"):
            continue
        match = re.match(r"\[\s*(\w+)\s*\]", stripped)
        if match:
            section = match.group(1)
            continue
        if section == "atoms":
            parts = stripped.split()
            if len(parts) >= 5 and parts[0].isdigit():
                atom_map[(int(parts[2]), parts[3], parts[4])] = int(parts[0])
        elif section == "bonds":
            parts = stripped.split()
            if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
                bonds.add(frozenset((int(parts[0]), int(parts[1]))))
    return atom_map, bonds


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gro", type=Path, required=True)
    parser.add_argument("--top", type=Path, required=True)
    parser.add_argument("--resname", action="append", help="Custom residue name to audit. Repeatable.")
    args = parser.parse_args()

    wanted = args.resname or list(DEFAULT_BONDS)
    gro_atoms = parse_gro(args.gro)
    atom_map, top_bonds = parse_top_bonds(args.top)
    status = 0
    for resname in wanted:
        residues = sorted({int(a["resid"]) for a in gro_atoms if a["resname"] == resname})
        if not residues:
            print(f"residue={resname}:not_found")
            status = 2
            continue
        for resid in residues:
            residue = [a for a in gro_atoms if int(a["resid"]) == resid and a["resname"] == resname]
            by_atom = {str(a["atom"]): a for a in residue}
            print(f"residue={resname}:{resid} atoms={len(residue)}")
            for a, b in DEFAULT_BONDS.get(resname, [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB")]):
                if a not in by_atom or b not in by_atom:
                    continue
                d = dist(by_atom[a], by_atom[b])
                atom_a = atom_map.get((resid, resname, a))
                atom_b = atom_map.get((resid, resname, b))
                has_top_bond = atom_a is not None and atom_b is not None and frozenset((atom_a, atom_b)) in top_bonds
                ok = has_top_bond
                if (a, b) == ("CA", "CB") and not (1.42 <= d <= 1.70):
                    ok = False
                print(
                    f"  {a}-{b} distance_A={d:.3f} "
                    f"topology_bond={has_top_bond} status={'ok' if ok else 'warn'}"
                )
                if not ok:
                    status = 2
    return status


if __name__ == "__main__":
    raise SystemExit(main())
