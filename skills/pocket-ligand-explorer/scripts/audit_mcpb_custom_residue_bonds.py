#!/usr/bin/env python3
"""Audit and optionally patch MCPB custom amino-acid mol2 connectivity.

MCPB.py can generate custom residue templates such as HD1/HD2/GU1 for
metal-coordinating residues.  tleap will accept a syntactically valid mol2 even
if an internal covalent bond such as CA-CB is missing, which can destroy a
production MD trajectory without an obvious tleap error.  This script checks the
minimum backbone/side-chain bonds and, when requested, patches missing single
bonds in the mol2 files.
"""

from __future__ import annotations

import argparse
import math
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


HIS_LIKE_PREFIXES = ("HD", "HE", "HP", "HI")
GLU_LIKE_PREFIXES = ("GU", "GLU", "GLH")
ASP_LIKE_PREFIXES = ("AD", "ASD", "ASP", "ASH")


@dataclass
class Mol2Data:
    path: Path
    lines: list[str]
    atom_start: int | None
    bond_start: int | None
    sub_start: int | None
    atoms: dict[str, int]
    atom_names_by_id: dict[int, str]
    bonds: set[frozenset[int]]
    counts_line_index: int | None
    bond_count: int


def parse_mol2(path: Path) -> Mol2Data:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    atom_start = bond_start = sub_start = None
    counts_line_index = None
    bond_count = 0
    for idx, line in enumerate(lines):
        if line.startswith("@<TRIPOS>ATOM"):
            atom_start = idx
        elif line.startswith("@<TRIPOS>BOND"):
            bond_start = idx
        elif line.startswith("@<TRIPOS>SUBSTRUCTURE"):
            sub_start = idx

    for idx, line in enumerate(lines):
        if idx <= 1:
            continue
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
            counts_line_index = idx
            bond_count = int(parts[1])
            break

    atoms: dict[str, int] = {}
    atom_names_by_id: dict[int, str] = {}
    if atom_start is not None:
        end = bond_start if bond_start is not None else len(lines)
        for line in lines[atom_start + 1 : end]:
            if not line.strip() or line.startswith("@<TRIPOS>"):
                break
            parts = line.split()
            if len(parts) < 2 or not parts[0].isdigit():
                continue
            atom_id = int(parts[0])
            atom_name = parts[1]
            atoms[atom_name] = atom_id
            atom_names_by_id[atom_id] = atom_name

    bonds: set[frozenset[int]] = set()
    if bond_start is not None:
        end = sub_start if sub_start is not None else len(lines)
        for line in lines[bond_start + 1 : end]:
            if not line.strip() or line.startswith("@<TRIPOS>"):
                break
            parts = line.split()
            if len(parts) < 4 or not parts[0].isdigit():
                continue
            try:
                bonds.add(frozenset((int(parts[1]), int(parts[2]))))
            except ValueError:
                continue

    return Mol2Data(
        path=path,
        lines=lines,
        atom_start=atom_start,
        bond_start=bond_start,
        sub_start=sub_start,
        atoms=atoms,
        atom_names_by_id=atom_names_by_id,
        bonds=bonds,
        counts_line_index=counts_line_index,
        bond_count=bond_count,
    )


def required_bonds(mol_name: str, atoms: dict[str, int]) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    if {"N", "CA", "C", "O"}.issubset(atoms):
        pairs.extend([("N", "CA"), ("CA", "C"), ("C", "O")])
    if {"CA", "CB"}.issubset(atoms):
        pairs.append(("CA", "CB"))
    if {"CB", "CG"}.issubset(atoms) and mol_name.upper().startswith(HIS_LIKE_PREFIXES + GLU_LIKE_PREFIXES):
        pairs.append(("CB", "CG"))
    if {"CG", "CD"}.issubset(atoms) and mol_name.upper().startswith(GLU_LIKE_PREFIXES):
        pairs.append(("CG", "CD"))
    if {"CB", "CG"}.issubset(atoms) and mol_name.upper().startswith(ASP_LIKE_PREFIXES):
        pairs.append(("CB", "CG"))
    return pairs


def has_bond(data: Mol2Data, a: str, b: str) -> bool:
    return frozenset((data.atoms[a], data.atoms[b])) in data.bonds


def add_single_bonds(data: Mol2Data, missing: list[tuple[str, str]], backup: bool) -> None:
    if not missing:
        return
    if data.bond_start is None:
        raise SystemExit(f"{data.path}: no @<TRIPOS>BOND section")
    if data.counts_line_index is None:
        raise SystemExit(f"{data.path}: could not locate mol2 counts line")
    if data.sub_start is None:
        insert_at = len(data.lines)
    else:
        insert_at = data.sub_start

    if backup:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_path = data.path.with_suffix(data.path.suffix + f".bak_missing_bonds_{stamp}")
        shutil.copy2(data.path, backup_path)
        print(f"backup={backup_path}")

    new_lines = list(data.lines)
    next_id = 1
    for line in data.lines[data.bond_start + 1 : insert_at]:
        parts = line.split()
        if parts and parts[0].isdigit():
            next_id = max(next_id, int(parts[0]) + 1)

    additions: list[str] = []
    for a, b in missing:
        additions.append(f"{next_id:6d} {data.atoms[a]:5d} {data.atoms[b]:5d} 1")
        print(f"patched={data.path.name} bond={a}-{b} atom_ids={data.atoms[a]}-{data.atoms[b]}")
        next_id += 1

    new_lines[insert_at:insert_at] = additions
    counts = new_lines[data.counts_line_index].split()
    counts[1] = str(int(counts[1]) + len(additions))
    new_lines[data.counts_line_index] = " ".join(counts)
    data.path.write_text("\n".join(new_lines) + "\n", encoding="utf-8")


def parse_pdb_atoms(path: Path) -> dict[tuple[str, str, int], dict[str, tuple[float, float, float]]]:
    residues: dict[tuple[str, str, int], dict[str, tuple[float, float, float]]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            resseq = int(line[22:26])
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        atom = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21:22].strip()
        residues.setdefault((resname, chain, resseq), {})[atom] = xyz
    return residues


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def audit_pdb(path: Path, custom_names: set[str]) -> int:
    warnings = 0
    residues = parse_pdb_atoms(path)
    for (resname, chain, resseq), atoms in sorted(residues.items(), key=lambda x: (x[0][1], x[0][2], x[0][0])):
        if resname.upper() not in custom_names:
            continue
        for a, b in required_bonds(resname, {name: idx for idx, name in enumerate(atoms, start=1)}):
            if a not in atoms or b not in atoms:
                continue
            d = distance(atoms[a], atoms[b])
            status = "ok"
            if (a, b) == ("CA", "CB") and not (1.42 <= d <= 1.70):
                status = "warn"
                warnings += 1
            print(f"pdb={path.name}\tres={resname}:{chain}:{resseq}\tbond={a}-{b}\tdistance_A={d:.3f}\tstatus={status}")
    return warnings


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mcpb-dir", type=Path, required=True)
    parser.add_argument("--files", nargs="*", help="Specific mol2 filenames to audit or patch.")
    parser.add_argument("--pdb", type=Path, help="Optional MCPB PDB to audit for covalent distances.")
    parser.add_argument("--fix", action="store_true", help="Patch missing required mol2 single bonds.")
    parser.add_argument("--no-backup", action="store_true", help="Do not write .bak files before --fix.")
    args = parser.parse_args()

    files = [args.mcpb_dir / name for name in args.files] if args.files else sorted(args.mcpb_dir.glob("*.mol2"))
    total_missing = 0
    for path in files:
        if not path.exists() or not path.name.lower().endswith(".mol2"):
            continue
        data = parse_mol2(path)
        mol_name = path.stem
        missing: list[tuple[str, str]] = []
        for a, b in required_bonds(mol_name, data.atoms):
            present = has_bond(data, a, b)
            print(f"mol2={path.name}\tbond={a}-{b}\tstatus={'present' if present else 'missing'}")
            if not present:
                missing.append((a, b))
        total_missing += len(missing)
        if args.fix:
            add_single_bonds(data, missing, backup=not args.no_backup)

    custom_names = {path.stem.upper() for path in files}
    pdb_warnings = audit_pdb(args.pdb, custom_names) if args.pdb else 0
    print(f"summary_missing_mol2_bonds={total_missing}")
    if args.pdb:
        print(f"summary_pdb_distance_warnings={pdb_warnings}")
    mol2_failed = total_missing and not args.fix
    return 2 if mol2_failed or pdb_warnings else 0


if __name__ == "__main__":
    raise SystemExit(main())
