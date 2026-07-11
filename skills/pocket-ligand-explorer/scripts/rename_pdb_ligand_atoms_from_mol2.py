#!/usr/bin/env python3
"""Rename ligand PDB atom names to match a mol2 template without moving atoms.

The PDB file is treated as a coordinate container. The mol2 file is treated as
the chemistry/name template. By default, hydrogen parent atoms are inferred from
distances so broken or misleading PDB CONECT records do not drive the mapping.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
import math


def atom_name(line: str) -> str:
    return line[12:16].strip()


def resname(line: str) -> str:
    return line[17:20].strip()


def serial(line: str) -> int:
    return int(line[6:11])


def xyz(line: str) -> tuple[float, float, float]:
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def with_atom_name(line: str, name: str) -> str:
    if len(name) > 4:
        raise ValueError(f"PDB atom names cannot exceed 4 characters: {name}")
    return line[:12] + f"{name:>4}" + line[16:]


def parse_mol2_atom_bonds(path: Path) -> tuple[dict[int, str], list[tuple[int, int]]]:
    atoms: dict[int, str] = {}
    bonds: list[tuple[int, int]] = []
    section: str | None = None
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("@<TRIPOS>"):
            section = line
            continue
        if section == "@<TRIPOS>ATOM":
            parts = line.split()
            atoms[int(parts[0])] = parts[1]
        elif section == "@<TRIPOS>BOND":
            parts = line.split()
            bonds.append((int(parts[1]), int(parts[2])))
    return atoms, bonds


def hydrogen_parent_map(names: dict[int, str], bonds: list[tuple[int, int]]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = defaultdict(list)
    for a, b in bonds:
        an = names.get(a, "")
        bn = names.get(b, "")
        if an.upper().startswith("H") and not bn.upper().startswith("H"):
            out[bn].append(an)
        elif bn.upper().startswith("H") and not an.upper().startswith("H"):
            out[an].append(bn)
    return {parent: sorted(children) for parent, children in out.items()}


def pdb_hydrogen_parent_map(lines: list[str], ligand_resname: str) -> dict[str, list[tuple[int, str]]]:
    serial_to_name: dict[int, str] = {}
    ligand_serials: set[int] = set()
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            try:
                ser = serial(line)
            except ValueError:
                continue
            serial_to_name[ser] = atom_name(line)
            if resname(line) == ligand_resname:
                ligand_serials.add(ser)

    out: dict[str, dict[int, str]] = defaultdict(dict)
    for line in lines:
        if not line.startswith("CONECT"):
            continue
        fields = line.split()
        if len(fields) < 3:
            continue
        try:
            nums = [int(item) for item in fields[1:]]
        except ValueError:
            continue
        src = nums[0]
        if src not in ligand_serials:
            continue
        src_name = serial_to_name.get(src, "")
        for dst in nums[1:]:
            if dst not in ligand_serials:
                continue
            dst_name = serial_to_name.get(dst, "")
            if src_name.upper().startswith("H") and not dst_name.upper().startswith("H"):
                out[dst_name][src] = src_name
            elif dst_name.upper().startswith("H") and not src_name.upper().startswith("H"):
                out[src_name][dst] = dst_name
    return {
        parent: sorted(children.items(), key=lambda item: (item[1], item[0]))
        for parent, children in out.items()
    }


def pdb_hydrogen_parent_map_by_distance(
    lines: list[str],
    ligand_resname: str,
    max_distance: float,
) -> dict[str, list[tuple[int, str]]]:
    heavy: list[tuple[int, str, tuple[float, float, float]]] = []
    hydrogens: list[tuple[int, str, tuple[float, float, float]]] = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")) or resname(line) != ligand_resname:
            continue
        name = atom_name(line)
        item = (serial(line), name, xyz(line))
        if name.upper().startswith("H"):
            hydrogens.append(item)
        else:
            heavy.append(item)
    out: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for h_serial, h_name, h_xyz in hydrogens:
        best: tuple[float, int, str] | None = None
        for heavy_serial, heavy_name, heavy_xyz in heavy:
            d = distance(h_xyz, heavy_xyz)
            candidate = (d, heavy_serial, heavy_name)
            if best is None or candidate < best:
                best = candidate
        if best is None:
            continue
        d, _, parent = best
        if d <= max_distance:
            out[parent].append((h_serial, h_name))
        else:
            out[f"UNASSIGNED>{max_distance:.2f}A"].append((h_serial, h_name))
    return {parent: sorted(children, key=lambda item: item[0]) for parent, children in out.items()}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Rename PDB ligand atom names to match a mol2 template using heavy-atom hydrogen parent bonds."
    )
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--mol2", type=Path, required=True)
    parser.add_argument("--resname", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--parent-mode", choices=["distance", "conect"], default="distance")
    parser.add_argument("--max-h-distance", type=float, default=1.35)
    args = parser.parse_args()

    lines = args.pdb.read_text(encoding="utf-8", errors="replace").splitlines()
    mol2_atoms, mol2_bonds = parse_mol2_atom_bonds(args.mol2)
    mol2_h = hydrogen_parent_map(mol2_atoms, mol2_bonds)
    if args.parent_mode == "distance":
        pdb_h = pdb_hydrogen_parent_map_by_distance(lines, args.resname, args.max_h_distance)
    else:
        pdb_h = pdb_hydrogen_parent_map(lines, args.resname)

    rename_by_serial: dict[int, str] = {}
    report_rows: list[dict[str, object]] = []
    for parent, pdb_children in sorted(pdb_h.items()):
        mol2_children = mol2_h.get(parent, [])
        if not mol2_children:
            report_rows.append({
                "parent": parent,
                "status": "missing_parent_in_mol2",
                "pdb_h": [name for _, name in pdb_children],
                "mol2_h": mol2_children,
            })
            continue
        if len(pdb_children) != len(mol2_children):
            report_rows.append({
                "parent": parent,
                "status": "hydrogen_count_mismatch",
                "pdb_h": [name for _, name in pdb_children],
                "mol2_h": mol2_children,
            })
        for (ser, old_name), new_name in zip(pdb_children, mol2_children):
            rename_by_serial[ser] = new_name
            report_rows.append({
                "parent": parent,
                "serial": ser,
                "old": old_name,
                "new": new_name,
                "status": "renamed" if old_name != new_name else "unchanged",
            })

    out_lines: list[str] = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            try:
                ser = serial(line)
            except ValueError:
                out_lines.append(line)
                continue
            if ser in rename_by_serial:
                line = with_atom_name(line, rename_by_serial[ser])
        out_lines.append(line)
    args.out.write_text("\n".join(out_lines) + "\n", encoding="utf-8", newline="\n")

    report = args.report or args.out.with_suffix(args.out.suffix + ".atom_rename_report.json")
    report.write_text(json.dumps(report_rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"out={args.out}")
    print(f"report={report}")
    for row in report_rows:
        if row.get("status") == "renamed":
            print(f"{row['parent']}: {row['old']} -> {row['new']} serial={row['serial']}")
        elif row.get("status") not in {"unchanged"}:
            print(f"warning={row}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
