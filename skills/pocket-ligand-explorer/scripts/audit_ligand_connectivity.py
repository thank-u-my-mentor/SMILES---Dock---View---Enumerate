#!/usr/bin/env python3
"""Audit PDB CONECT/PyMOL-style close contacts versus mol2 bonds for one ligand."""

from __future__ import annotations

import argparse
import math
from pathlib import Path


def parse_pdb(path: Path, resname: str, resseq: str, chain: str | None):
    atoms = []
    conect: dict[int, list[int]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        rec = line[:6].strip().upper()
        if rec in {"ATOM", "HETATM"}:
            name = line[12:16].strip()
            res = line[17:20].strip()
            ch = line[21].strip()
            rs = line[22:26].strip()
            if res == resname and rs == resseq and (chain is None or ch == chain):
                atoms.append(
                    {
                        "serial": int(line[6:11]),
                        "name": name,
                        "resname": res,
                        "chain": ch,
                        "resseq": rs,
                        "x": float(line[30:38]),
                        "y": float(line[38:46]),
                        "z": float(line[46:54]),
                    }
                )
        elif rec == "CONECT":
            parts = line.split()
            if len(parts) >= 2:
                try:
                    first = int(parts[1])
                    conect.setdefault(first, [])
                    for token in parts[2:]:
                        conect[first].append(int(token))
                except ValueError:
                    pass
    return atoms, conect


def parse_mol2(path: Path):
    atoms = {}
    bonds = []
    section = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("@<TRIPOS>"):
            section = line.strip().upper()
            continue
        if not line.strip():
            continue
        if section == "@<TRIPOS>ATOM":
            parts = line.split()
            if len(parts) >= 6:
                atoms[int(parts[0])] = {"name": parts[1], "type": parts[5]}
        elif section == "@<TRIPOS>BOND":
            parts = line.split()
            if len(parts) >= 4:
                bonds.append((int(parts[1]), int(parts[2]), parts[3]))
    return atoms, bonds


def dist(a, b) -> float:
    return math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--mol2", type=Path, required=True)
    parser.add_argument("--resname", default="UNK")
    parser.add_argument("--resseq", default="501")
    parser.add_argument("--chain", default="A")
    parser.add_argument("--atom-a", default="C1")
    parser.add_argument("--atom-b", default="C3")
    parser.add_argument("--close-cutoff", type=float, default=1.9)
    args = parser.parse_args()

    atoms, conect = parse_pdb(args.pdb, args.resname, args.resseq, args.chain)
    by_serial = {a["serial"]: a for a in atoms}
    by_name = {a["name"]: a for a in atoms}
    print(f"pdb={args.pdb}")
    print(f"mol2={args.mol2}")
    print(f"ligand={args.resname} {args.chain}:{args.resseq}")
    print(f"pdb_atoms={len(atoms)}")
    print("pdb_atom_names=" + ",".join(a["name"] for a in atoms))

    atom_a = by_name.get(args.atom_a)
    atom_b = by_name.get(args.atom_b)
    if atom_a and atom_b:
        print(f"pdb_distance_{args.atom_a}_{args.atom_b}_A={dist(atom_a, atom_b):.3f}")
    else:
        print(f"pdb_distance_{args.atom_a}_{args.atom_b}_A=missing_atom")

    pdb_bonds = set()
    for src, neighs in conect.items():
        if src not in by_serial:
            continue
        for dst in neighs:
            if dst in by_serial:
                pdb_bonds.add(tuple(sorted((by_serial[src]["name"], by_serial[dst]["name"]))))
    target_pair = tuple(sorted((args.atom_a, args.atom_b)))
    print("pdb_conect_internal=" + (";".join(f"{x}-{y}" for x, y in sorted(pdb_bonds)) or "none"))
    print(f"pdb_conect_has_{args.atom_a}_{args.atom_b}={str(target_pair in pdb_bonds).lower()}")

    close_pairs = []
    for i, atom_i in enumerate(atoms):
        for atom_j in atoms[i + 1 :]:
            d = dist(atom_i, atom_j)
            if d <= args.close_cutoff:
                close_pairs.append((d, atom_i["name"], atom_j["name"]))
    print(
        "pdb_close_pairs_le_%.2fA=%s"
        % (
            args.close_cutoff,
            ";".join(f"{x}-{y}:{d:.3f}" for d, x, y in sorted(close_pairs)) or "none",
        )
    )

    mol2_atoms, mol2_bonds = parse_mol2(args.mol2)
    mol2_pairs = []
    for i, j, typ in mol2_bonds:
        if i in mol2_atoms and j in mol2_atoms:
            n1 = mol2_atoms[i]["name"]
            n2 = mol2_atoms[j]["name"]
            mol2_pairs.append((tuple(sorted((n1, n2))), typ))
    print(f"mol2_bonds={len(mol2_pairs)}")
    print(f"mol2_has_{args.atom_a}_{args.atom_b}={str(any(pair == target_pair for pair, _ in mol2_pairs)).lower()}")
    print(
        "mol2_bonds_touching_%s_or_%s=%s"
        % (
            args.atom_a,
            args.atom_b,
            ";".join(
                f"{x}-{y}:{typ}"
                for (x, y), typ in sorted(mol2_pairs)
                if args.atom_a in (x, y) or args.atom_b in (x, y)
            )
            or "none",
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
