#!/usr/bin/env python3
"""Build a human visual-check PDB for a MCPB small model.

MCPB.py already writes ``*_small.pdb``.  This script keeps that atom collection
but adds conservative visualization CONECT records:

* mol2-defined bonds for non-standard residues such as UNL/UNK/ACT/HOH
* Fe-to-nearby N/O/S donor links for visual inspection only

It does not change the chemistry used by MCPB or Amber.  Ligand/radical
hydrogen states remain exactly as supplied in the small PDB and mol2 files.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


DONOR_ELEMENTS = {"N", "O", "S"}
METAL_ELEMENTS = {"FE", "ZN", "CU", "MN", "MG", "CO", "NI"}
SKIP_MOL2 = {"FE"}


@dataclass(frozen=True)
class Atom:
    serial: int
    name: str
    resname: str
    chain: str
    resseq: str
    icode: str
    element: str
    xyz: tuple[float, float, float]
    line: str

    @property
    def resid(self) -> tuple[str, str, str, str]:
        return (self.resname, self.chain, self.resseq, self.icode)


def guess_element(atom_name: str, element_field: str = "") -> str:
    raw = "".join(ch for ch in element_field.strip() if ch.isalpha()).upper()
    if raw:
        return raw
    letters = "".join(ch for ch in atom_name.strip() if ch.isalpha()).upper()
    if letters.startswith("FE"):
        return "FE"
    if len(letters) >= 2 and letters[:2] in METAL_ELEMENTS | {"CL", "BR"}:
        return letters[:2]
    return letters[:1]


def parse_pdb(path: Path) -> tuple[list[str], list[Atom]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    atoms: list[Atom] = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            atoms.append(
                Atom(
                    serial=int(line[6:11]),
                    name=line[12:16].strip(),
                    resname=line[17:20].strip(),
                    chain=line[21:22].strip(),
                    resseq=line[22:26].strip(),
                    icode=line[26:27].strip(),
                    element=guess_element(line[12:16], line[76:78] if len(line) >= 78 else ""),
                    xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
                    line=(line + " " * 80)[:80],
                )
            )
        except ValueError:
            continue
    if not atoms:
        raise SystemExit(f"No atoms parsed from {path}")
    return lines, atoms


def distance(a: Atom, b: Atom) -> float:
    return math.sqrt(sum((a.xyz[i] - b.xyz[i]) ** 2 for i in range(3)))


def grouped_residues(atoms: list[Atom]) -> dict[tuple[str, str, str, str], list[Atom]]:
    groups: dict[tuple[str, str, str, str], list[Atom]] = {}
    for atom in atoms:
        groups.setdefault(atom.resid, []).append(atom)
    return groups


def parse_mol2_bonds(path: Path) -> list[tuple[str, str, str]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    section = ""
    atom_names: dict[int, str] = {}
    bonds: list[tuple[str, str, str]] = []
    for line in lines:
        if line.startswith("@<TRIPOS>"):
            section = line.strip()
            continue
        if section == "@<TRIPOS>ATOM":
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit():
                atom_names[int(parts[0])] = parts[1]
        elif section == "@<TRIPOS>BOND":
            parts = line.split()
            if len(parts) >= 4 and parts[0].isdigit():
                a = atom_names.get(int(parts[1]))
                b = atom_names.get(int(parts[2]))
                if a and b:
                    bonds.append((a, b, parts[3]))
    return bonds


def add_edge(edges: dict[int, set[int]], a: int, b: int) -> None:
    if a == b:
        return
    edges.setdefault(a, set()).add(b)
    edges.setdefault(b, set()).add(a)


def build_edges(
    atoms: list[Atom],
    *,
    mol2_dir: Path | None,
    donor_cutoff: float,
) -> tuple[dict[int, set[int]], list[dict[str, str]]]:
    edges: dict[int, set[int]] = {}
    report: list[dict[str, str]] = []
    residues = grouped_residues(atoms)

    if mol2_dir:
        for resid, residue_atoms in residues.items():
            resname = resid[0]
            if resname.upper() in SKIP_MOL2:
                continue
            mol2 = mol2_dir / f"{resname}.mol2"
            if not mol2.exists():
                continue
            atom_by_name = {atom.name: atom for atom in residue_atoms}
            for a_name, b_name, bond_type in parse_mol2_bonds(mol2):
                a = atom_by_name.get(a_name)
                b = atom_by_name.get(b_name)
                if not a or not b:
                    report.append(
                        {
                            "kind": "mol2_bond",
                            "residue": f"{resname}:{resid[1]}:{resid[2]}",
                            "atom_a": a_name,
                            "atom_b": b_name,
                            "distance_A": "",
                            "status": "missing_atom_in_small_pdb",
                            "note": mol2.name,
                        }
                    )
                    continue
                add_edge(edges, a.serial, b.serial)
                report.append(
                    {
                        "kind": "mol2_bond",
                        "residue": f"{resname}:{resid[1]}:{resid[2]}",
                        "atom_a": a.name,
                        "atom_b": b.name,
                        "distance_A": f"{distance(a, b):.3f}",
                        "status": "visual_conect_added",
                        "note": f"{mol2.name}:{bond_type}",
                    }
                )

    metals = [atom for atom in atoms if atom.element.upper() in METAL_ELEMENTS or atom.name.upper() in METAL_ELEMENTS]
    for metal in metals:
        for atom in atoms:
            if atom.serial == metal.serial or atom.element.upper() not in DONOR_ELEMENTS:
                continue
            d = distance(metal, atom)
            if d <= donor_cutoff:
                add_edge(edges, metal.serial, atom.serial)
                report.append(
                    {
                        "kind": "metal_donor",
                        "residue": f"{atom.resname}:{atom.chain}:{atom.resseq}",
                        "atom_a": metal.name,
                        "atom_b": atom.name,
                        "distance_A": f"{d:.3f}",
                        "status": "visual_conect_added",
                        "note": "Fe-donor visual link only; not an Amber/GAMESS parameter",
                    }
                )
    return edges, report


def conect_lines(edges: dict[int, set[int]]) -> list[str]:
    lines: list[str] = []
    for serial in sorted(edges):
        partners = sorted(partner for partner in edges[serial] if partner > serial)
        if not partners:
            continue
        for idx in range(0, len(partners), 4):
            chunk = partners[idx : idx + 4]
            lines.append("CONECT" + f"{serial:5d}" + "".join(f"{partner:5d}" for partner in chunk))
    return lines


def write_report(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["kind", "residue", "atom_a", "atom_b", "distance_A", "status", "note"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--small-pdb", type=Path, required=True)
    parser.add_argument("--mol2-dir", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--donor-cutoff", type=float, default=2.8)
    parser.add_argument("--model-charge", type=int)
    parser.add_argument("--mult", type=int)
    args = parser.parse_args()

    lines, atoms = parse_pdb(args.small_pdb)
    edges, report = build_edges(atoms, mol2_dir=args.mol2_dir, donor_cutoff=args.donor_cutoff)

    remarks = [
        "REMARK PLE MCPB SMALL MODEL VISUAL CHECK",
        f"REMARK SOURCE_SMALL_PDB {args.small_pdb}",
        "REMARK CONECT FROM LIGAND MOL2 AND METAL-DONOR DISTANCE ONLY",
        "REMARK DO NOT USE THIS FILE AS A FORCE-FIELD PARAMETER FILE",
        "REMARK LIGAND/RADICAL HYDROGENS ARE READ-ONLY AND NOT AUTO-REPAIRED",
    ]
    if args.model_charge is not None and args.mult is not None:
        remarks.append(f"REMARK MODEL_CHARGE {args.model_charge} MULT {args.mult}")
    atom_lines = [line for line in lines if line.startswith(("ATOM", "HETATM", "TER"))]
    output = remarks + atom_lines + conect_lines(edges) + ["END"]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(output).rstrip() + "\n")
    report_path = args.report or args.out.with_suffix(args.out.suffix + ".report.tsv")
    write_report(report_path, report)
    print(f"visual_pdb={args.out}")
    print(f"report={report_path}")
    print(f"visual_conect_edges={sum(len(v) for v in edges.values()) // 2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
