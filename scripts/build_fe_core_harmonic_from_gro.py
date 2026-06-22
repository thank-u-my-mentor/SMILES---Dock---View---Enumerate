#!/usr/bin/env python3
"""Find Fe coordination atoms in a GRO file and emit harmonic-topology bonds."""

from __future__ import annotations

import argparse
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Atom:
    resid: int
    resname: str
    atom: str
    index: int
    x: float
    y: float
    z: float


def parse_gro(path: Path) -> list[Atom]:
    lines = path.read_text().splitlines()
    atoms: list[Atom] = []
    for line in lines[2:-1]:
        if len(line) < 44:
            continue
        atoms.append(
            Atom(
                resid=int(line[0:5]),
                resname=line[5:10].strip(),
                atom=line[10:15].strip(),
                index=int(line[15:20]),
                x=float(line[20:28]),
                y=float(line[28:36]),
                z=float(line[36:44]),
            )
        )
    return atoms


def dist_nm(a: Atom, b: Atom) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gro", type=Path, required=True)
    parser.add_argument("--input-top", type=Path, required=True)
    parser.add_argument("--output-top", type=Path, required=True)
    parser.add_argument("--helper", type=Path, default=Path("/mnt/e/Codex/scripts/make_fe_core_harmonic_top.py"))
    parser.add_argument("--cutoff-nm", type=float, default=0.30)
    parser.add_argument("--force", type=float, default=100000.0)
    args = parser.parse_args()

    atoms = parse_gro(args.gro)
    fe_atoms = [a for a in atoms if a.atom.upper() in {"FE", "FE2", "FE3"} and a.resname.upper() in {"FE", "FE2", "FE3"}]
    if not fe_atoms:
        fe_atoms = [a for a in atoms if a.atom.upper() in {"FE", "FE2", "FE3"}]
    if len(fe_atoms) != 1:
        print(f"Expected exactly one Fe atom; found {len(fe_atoms)}", file=sys.stderr)
        for atom in fe_atoms[:20]:
            print(atom, file=sys.stderr)
        return 2
    fe = fe_atoms[0]

    candidates = []
    for atom in atoms:
        if atom.index == fe.index:
            continue
        if atom.atom.upper() not in {"NE2", "ND1", "OE1", "OE2", "OD1", "OD2"}:
            continue
        d = dist_nm(fe, atom)
        if d <= args.cutoff_nm:
            candidates.append((d, atom))
    candidates.sort(key=lambda item: item[0])

    selected = []
    seen = set()
    for d, atom in candidates:
        key = (atom.resid, atom.resname)
        if key in seen:
            continue
        if atom.resname.upper() in {"HID", "HIE", "HIP", "HIS", "GLU", "ASP"}:
            selected.append((d, atom))
            seen.add(key)
        if len(selected) >= 3:
            break
    if len(selected) < 3:
        print(f"Found only {len(selected)} Fe contacts within {args.cutoff_nm:.3f} nm", file=sys.stderr)
        for d, atom in candidates:
            print(f"{atom.index} {atom.resname}{atom.resid}:{atom.atom} {d:.4f} nm", file=sys.stderr)
        return 3

    cmd = [
        sys.executable,
        str(args.helper),
        "--input-top",
        str(args.input_top),
        "--output-top",
        str(args.output_top),
    ]
    print(f"fe_atom={fe.index} {fe.resname}{fe.resid}:{fe.atom}")
    for d, atom in selected:
        comment = f"Fe3core_{atom.resname}{atom.resid}_{atom.atom}"
        print(f"contact={atom.index} {atom.resname}{atom.resid}:{atom.atom} distance_nm={d:.4f}")
        cmd.extend(["--bond", f"{fe.index},{atom.index},{d:.4f},{args.force:.1f},{comment}"])

    subprocess.run(cmd, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
