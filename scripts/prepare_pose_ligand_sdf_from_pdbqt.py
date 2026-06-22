#!/usr/bin/env python3
"""Write a pose SDF in template-SDF atom order using PDBQT atom mapping.

Docking engines often reorder atoms in PDBQT. This script maps template SDF
atoms to the input PDBQT by matching the original 3D coordinates, then applies
coordinates from a selected output PDBQT MODEL. The output keeps SDF bond
orders and writes heavy atoms only; add hydrogens afterwards with OpenBabel.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


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
                "idx": i + 1,
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


def parse_pdbqt_atoms(path: Path, model: int | None = None):
    atoms = []
    current_model = None
    for line in path.read_text().splitlines():
        if line.startswith("MODEL"):
            parts = line.split()
            current_model = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else None
            continue
        if line.startswith("ENDMDL") and model is not None and current_model == model:
            break
        if model is not None and current_model != model:
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        ad_type = line[77:].strip().split()[-1] if len(line) > 77 and line[77:].strip() else name[:1]
        elem = {"A": "C", "C": "C", "N": "N", "NA": "N", "OA": "O", "O": "O", "S": "S", "SA": "S", "HD": "H", "H": "H"}.get(ad_type.upper(), ad_type[:1].upper())
        atoms.append(
            {
                "serial": int(line[6:11]),
                "name": name,
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "elem": elem,
            }
        )
    if not atoms:
        raise SystemExit(f"No atoms parsed from {path}")
    return atoms


def dist(a, b) -> float:
    return math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--template-sdf", required=True, type=Path)
    ap.add_argument("--input-pdbqt", required=True, type=Path)
    ap.add_argument("--pose-pdbqt", required=True, type=Path)
    ap.add_argument("--model", type=int, default=1)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--title", default="LIG_pose_heavy")
    ap.add_argument("--tolerance", type=float, default=0.01)
    args = ap.parse_args()

    sdf_atoms, sdf_bonds = parse_sdf(args.template_sdf)
    input_atoms = parse_pdbqt_atoms(args.input_pdbqt, model=None)
    pose_atoms = parse_pdbqt_atoms(args.pose_pdbqt, model=args.model)
    if len(input_atoms) != len(pose_atoms):
        raise SystemExit(f"Input/pose PDBQT atom counts differ: {len(input_atoms)} vs {len(pose_atoms)}")

    heavy_sdf = [a for a in sdf_atoms if a["elem"].upper() != "H"]
    available = set(range(len(input_atoms)))
    sdf_to_pdbqt: dict[int, int] = {}
    max_error = 0.0
    for sdf_atom in heavy_sdf:
        candidates = [
            (dist(sdf_atom, input_atoms[i]), i)
            for i in available
            if input_atoms[i]["elem"].upper()[0] == sdf_atom["elem"].upper()[0]
        ]
        if not candidates:
            raise SystemExit(f"No PDBQT candidate for SDF atom {sdf_atom['idx']} {sdf_atom['elem']}")
        err, best = min(candidates)
        if err > args.tolerance:
            raise SystemExit(f"Mapping error too large for SDF atom {sdf_atom['idx']}: {err:.4f} A")
        sdf_to_pdbqt[sdf_atom["idx"]] = best
        available.remove(best)
        max_error = max(max_error, err)

    heavy_indices = [a["idx"] for a in heavy_sdf]
    index_map = {old: new for new, old in enumerate(heavy_indices, start=1)}
    heavy_bonds = [
        (index_map[a], index_map[b], order)
        for a, b, order in sdf_bonds
        if a in index_map and b in index_map
    ]

    out_lines = [
        args.title,
        "  Codex pdbqt-mapped pose coordinates",
        "",
        f"{len(heavy_indices):>3}{len(heavy_bonds):>3}  0  0  0  0            999 V2000",
    ]
    for old_idx in heavy_indices:
        sdf_atom = sdf_atoms[old_idx - 1]
        pose_atom = pose_atoms[sdf_to_pdbqt[old_idx]]
        out_lines.append(
            f"{pose_atom['x']:10.4f}{pose_atom['y']:10.4f}{pose_atom['z']:10.4f} {sdf_atom['elem']:<3} 0  0  0  0  0  0  0  0  0  0  0  0"
        )
    for a, b, order in heavy_bonds:
        out_lines.append(f"{a:>3}{b:>3}{order:>3}  0")
    out_lines.extend(["M  END", "$$$$"])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(out_lines) + "\n")
    print(f"wrote={args.out}")
    print(f"heavy_atoms={len(heavy_indices)} heavy_bonds={len(heavy_bonds)} max_mapping_error={max_error:.5f}")


if __name__ == "__main__":
    main()
