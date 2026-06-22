#!/usr/bin/env python3
"""Build a heavy-atom SDF by transplanting docked PDB coordinates.

The reference SDF supplies atom order and bond orders. The pose PDB supplies
the protein-frame coordinates for the heavy atoms with names like C1, O6, N8.
Hydrogens are intentionally omitted so OpenBabel can rebuild them near the
docked heavy atoms.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_pdb_coords(path: Path) -> dict[int, tuple[str, float, float, float]]:
    coords: dict[int, tuple[str, float, float, float]] = {}
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        digits = "".join(ch for ch in name if ch.isdigit())
        if not digits:
            continue
        idx = int(digits)
        elem = (line[76:78].strip() or "".join(ch for ch in name if ch.isalpha())).title()
        coords[idx] = (elem, float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return coords


def parse_sdf(path: Path):
    lines = path.read_text().splitlines()
    counts = lines[3]
    natoms = int(counts[0:3])
    nbonds = int(counts[3:6])
    atoms = []
    for i in range(natoms):
        line = lines[4 + i]
        atoms.append(
            {
                "x": float(line[0:10]),
                "y": float(line[10:20]),
                "z": float(line[20:30]),
                "elem": line[31:34].strip(),
            }
        )
    bonds = []
    for i in range(nbonds):
        line = lines[4 + natoms + i]
        bonds.append((int(line[0:3]), int(line[3:6]), int(line[6:9])))
    return atoms, bonds


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--template-sdf", required=True, type=Path)
    ap.add_argument("--pose-pdb", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--title", default="LIG_pose_heavy")
    args = ap.parse_args()

    atoms, bonds = parse_sdf(args.template_sdf)
    pdb_coords = parse_pdb_coords(args.pose_pdb)
    heavy_indices = [i + 1 for i, atom in enumerate(atoms) if atom["elem"].upper() != "H"]

    missing = [idx for idx in heavy_indices if idx not in pdb_coords]
    if missing:
        raise SystemExit(f"Missing docked coordinates for heavy atom indices: {missing}")

    index_map = {old: new for new, old in enumerate(heavy_indices, start=1)}
    heavy_bonds = [
        (index_map[a], index_map[b], order)
        for a, b, order in bonds
        if a in index_map and b in index_map
    ]

    out_lines = [
        args.title,
        "  Codex pose-coordinate transplant",
        "",
        f"{len(heavy_indices):>3}{len(heavy_bonds):>3}  0  0  0  0            999 V2000",
    ]
    for old_idx in heavy_indices:
        elem, x, y, z = pdb_coords[old_idx]
        tmpl_elem = atoms[old_idx - 1]["elem"]
        if elem.upper()[0] != tmpl_elem.upper()[0]:
            raise SystemExit(f"Atom element mismatch at {old_idx}: PDB {elem}, SDF {tmpl_elem}")
        out_lines.append(
            f"{x:10.4f}{y:10.4f}{z:10.4f} {tmpl_elem:<3} 0  0  0  0  0  0  0  0  0  0  0  0"
        )
    for a, b, order in heavy_bonds:
        out_lines.append(f"{a:>3}{b:>3}{order:>3}  0")
    out_lines.extend(["M  END", "$$$$"])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(out_lines) + "\n")
    print(f"wrote={args.out}")
    print(f"heavy_atoms={len(heavy_indices)} heavy_bonds={len(heavy_bonds)}")


if __name__ == "__main__":
    main()
