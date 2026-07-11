#!/usr/bin/env python3
"""Patch a ligand mol2 atom name/coordinate and one bond using PDB coordinates.

This is intentionally small and explicit.  It is useful for hand-built radical
or nitrene-like intermediates where antechamber/reduce would add the wrong
hydrogen, but MCPB.py still needs a mol2 topology whose atom names and bond
graph match the curated PDB.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def pdb_coords(path: Path, resname: str) -> dict[str, tuple[float, float, float]]:
    coords: dict[str, tuple[float, float, float]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip() != resname:
            continue
        coords[line[12:16].strip()] = (
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54]),
        )
    return coords


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mol2", type=Path, required=True)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--resname", required=True)
    parser.add_argument("--rename", nargs=2, metavar=("OLD", "NEW"), required=True)
    parser.add_argument("--replace-bond", nargs=2, metavar=("ATOM_A", "ATOM_B"), required=True)
    parser.add_argument("--new-type", default=None, help="Optional atom type for the renamed atom.")
    parser.add_argument("--set-type", action="append", default=[], metavar="ATOM=TYPE")
    args = parser.parse_args()

    coords = pdb_coords(args.pdb, args.resname)
    old_name, new_name = args.rename
    type_updates: dict[str, str] = {}
    for item in args.set_type:
        if "=" not in item:
            raise SystemExit(f"--set-type expects ATOM=TYPE, got {item!r}")
        name, atom_type = item.split("=", 1)
        type_updates[name] = atom_type
    if args.new_type:
        type_updates[new_name] = args.new_type
    required = {new_name, args.replace_bond[0], args.replace_bond[1]}
    missing = sorted(required - set(coords))
    if missing:
        raise SystemExit(f"{args.pdb}: missing {args.resname} coordinates for {missing}")

    lines = args.mol2.read_text(encoding="utf-8", errors="replace").splitlines()
    out: list[str] = []
    section = None
    id_by_name: dict[str, int] = {}
    atom_count = 0
    charge_sum = 0.0

    for line in lines:
        if line.startswith("@<TRIPOS>"):
            section = line.strip()
            out.append(line)
            continue
        if section == "@<TRIPOS>ATOM" and line.strip():
            parts = line.split()
            if len(parts) < 9:
                out.append(line)
                continue
            atom_id = int(parts[0])
            name = new_name if parts[1] == old_name else parts[1]
            atom_type = type_updates.get(name, parts[5])
            if name in coords:
                x, y, z = coords[name]
            else:
                x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            charge = float(parts[8])
            id_by_name[name] = atom_id
            charge_sum += charge
            atom_count += 1
            out.append(
                f"{atom_id:7d} {name:<8s} {x:10.4f} {y:10.4f} {z:10.4f} "
                f"{atom_type:<8s} {parts[6]:>3s} {args.resname:<8s} {charge:10.6f}"
            )
            continue
        if section == "@<TRIPOS>BOND" and line.strip():
            parts = line.split()
            if len(parts) >= 4:
                a = int(parts[1])
                b = int(parts[2])
                renamed_id = id_by_name.get(new_name)
                if renamed_id is not None and renamed_id in {a, b}:
                    left = id_by_name[args.replace_bond[0]]
                    right = id_by_name[args.replace_bond[1]]
                    out.append(f"{int(parts[0]):6d} {left:5d} {right:5d} 1   ")
                    continue
            out.append(line)
            continue
        out.append(line)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(out).rstrip() + "\n", encoding="utf-8", newline="\n")
    print(f"wrote={args.out}")
    print(f"atoms={atom_count}")
    print(f"charge_sum={charge_sum:.6f}")
    print(f"renamed={old_name}->{new_name}")
    print(f"replace_bond={args.replace_bond[0]}-{args.replace_bond[1]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
