#!/usr/bin/env python3
"""Audit Fe-coordinating histidine protonation in MCPB/Amber handoffs.

Amber standard histidine naming is tautomer-specific:

* HID: ND1-H present, NE2 unprotonated.
* HIE: NE2-H present, ND1 unprotonated.
* HIP: both nitrogens protonated and formally positive.

For a neutral Fe-bound histidine, the Fe-donor nitrogen should usually be the
unprotonated imidazole N.  This script reports the geometric Fe donor and the
observed N-H placement in PDB and/or mol2 files so a MCPB custom residue such
as HD1/HD2 is not silently accepted with impossible or unintended chemistry.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


HIS_NAMES = {"HIS", "HID", "HIE", "HIP", "HD1", "HD2", "HE1", "HE2", "HP1", "HP2"}
N_ATOMS = {"ND1", "NE2"}
HIS_RING_C = {"CD2", "CE1"}


@dataclass
class Atom:
    serial: int
    name: str
    resname: str
    chain: str
    resseq: int
    icode: str
    element: str
    xyz: tuple[float, float, float]

    @property
    def resid(self) -> tuple[str, str, int, str]:
        return (self.resname, self.chain, self.resseq, self.icode)


def dist(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def guess_element(atom_name: str, element_field: str = "") -> str:
    if element_field.strip():
        return element_field.strip().upper()
    name = atom_name.strip()
    if not name:
        return ""
    if len(name) >= 2 and name[:2].upper() in {"FE", "ZN", "CU", "MN", "MG", "CA", "CO", "NI"}:
        return name[:2].upper()
    return name[0].upper()


def read_pdb(path: Path) -> list[Atom]:
    atoms: list[Atom] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            serial = int(line[6:11])
            name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21:22].strip()
            resseq = int(line[22:26])
            icode = line[26:27].strip()
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            element = guess_element(name, line[76:78] if len(line) >= 78 else "")
        except ValueError:
            continue
        atoms.append(Atom(serial, name, resname, chain, resseq, icode, element, xyz))
    return atoms


def grouped_residues(atoms: list[Atom]) -> dict[tuple[str, str, int, str], list[Atom]]:
    residues: dict[tuple[str, str, int, str], list[Atom]] = {}
    for atom in atoms:
        residues.setdefault(atom.resid, []).append(atom)
    return residues


def nearest_h_to_n(res_atoms: list[Atom], n_atom: Atom, cutoff: float = 1.25) -> tuple[str, float] | tuple[str, None]:
    best: tuple[str, float] | None = None
    for atom in res_atoms:
        if atom.element != "H":
            continue
        d = dist(n_atom.xyz, atom.xyz)
        if d <= cutoff and (best is None or d < best[1]):
            best = (atom.name, d)
    return best if best else ("", None)


def expected_tautomer(donor: str) -> str:
    if donor == "NE2":
        return "HID"
    if donor == "ND1":
        return "HIE"
    return "unknown"


def donor_status(donor: str, nd1_has_h: bool, ne2_has_h: bool, resname: str) -> str:
    if donor == "NE2":
        if nd1_has_h and not ne2_has_h:
            return "ok_ne2_bound_neutral_hid"
        if ne2_has_h:
            return "warn_fe_donor_ne2_is_protonated"
        if not nd1_has_h and not ne2_has_h:
            return "warn_no_imidazole_nh_possible_imidazolate"
    if donor == "ND1":
        if ne2_has_h and not nd1_has_h:
            return "ok_nd1_bound_neutral_hie"
        if nd1_has_h:
            return "warn_fe_donor_nd1_is_protonated"
        if not nd1_has_h and not ne2_has_h:
            return "warn_no_imidazole_nh_possible_imidazolate"
    if resname.upper() == "HIP" and nd1_has_h and ne2_has_h:
        return "warn_hip_both_n_protonated_for_metal_site"
    return "review"


def audit_pdb(path: Path, metal_element: str, cutoff: float) -> list[dict[str, str]]:
    atoms = read_pdb(path)
    metals = [a for a in atoms if a.element.upper() == metal_element.upper() or a.name.upper() == metal_element.upper()]
    residues = grouped_residues(atoms)
    rows: list[dict[str, str]] = []
    for metal in metals:
        for resid, res_atoms in residues.items():
            resname = resid[0].upper()
            if resname not in HIS_NAMES:
                continue
            by_name = {a.name: a for a in res_atoms}
            if "ND1" not in by_name or "NE2" not in by_name:
                continue
            d_nd1 = dist(metal.xyz, by_name["ND1"].xyz)
            d_ne2 = dist(metal.xyz, by_name["NE2"].xyz)
            donor = "ND1" if d_nd1 <= d_ne2 else "NE2"
            donor_dist = min(d_nd1, d_ne2)
            if donor_dist > cutoff:
                continue
            nd1_h, nd1_h_d = nearest_h_to_n(res_atoms, by_name["ND1"])
            ne2_h, ne2_h_d = nearest_h_to_n(res_atoms, by_name["NE2"])
            ring_h = sum(1 for atom in res_atoms if atom.element == "H" and atom.name not in {nd1_h, ne2_h})
            rows.append(
                {
                    "source": str(path),
                    "kind": "pdb",
                    "metal": f"{metal.name}:{metal.chain}:{metal.resseq}",
                    "residue": f"{resid[0]}:{resid[1]}:{resid[2]}{resid[3]}",
                    "donor": donor,
                    "donor_distance_A": f"{donor_dist:.3f}",
                    "Fe_ND1_A": f"{d_nd1:.3f}",
                    "Fe_NE2_A": f"{d_ne2:.3f}",
                    "expected_neutral_amber_name": expected_tautomer(donor),
                    "ND1_H": nd1_h,
                    "ND1_H_distance_A": "" if nd1_h_d is None else f"{nd1_h_d:.3f}",
                    "NE2_H": ne2_h,
                    "NE2_H_distance_A": "" if ne2_h_d is None else f"{ne2_h_d:.3f}",
                    "ring_or_other_H_count": str(ring_h),
                    "status": donor_status(donor, bool(nd1_h), bool(ne2_h), resid[0]),
                }
            )
    return rows


@dataclass
class Mol2:
    atoms: dict[int, tuple[str, str]]
    atom_id_by_name: dict[str, int]
    bonds: set[frozenset[int]]


def read_mol2(path: Path) -> Mol2:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    section = ""
    atoms: dict[int, tuple[str, str]] = {}
    atom_id_by_name: dict[str, int] = {}
    bonds: set[frozenset[int]] = set()
    for line in lines:
        if line.startswith("@<TRIPOS>"):
            section = line.strip()
            continue
        parts = line.split()
        if section == "@<TRIPOS>ATOM" and len(parts) >= 6 and parts[0].isdigit():
            atom_id = int(parts[0])
            name = parts[1]
            atom_type = parts[5]
            atoms[atom_id] = (name, atom_type)
            atom_id_by_name[name] = atom_id
        elif section == "@<TRIPOS>BOND" and len(parts) >= 4 and parts[0].isdigit():
            try:
                bonds.add(frozenset((int(parts[1]), int(parts[2]))))
            except ValueError:
                pass
    return Mol2(atoms=atoms, atom_id_by_name=atom_id_by_name, bonds=bonds)


def atom_element_from_mol2_name_type(name: str, atom_type: str) -> str:
    if atom_type:
        token = atom_type.split(".", 1)[0]
        if token:
            return token.upper()
    return guess_element(name)


def bonded_hydrogens(mol2: Mol2, atom_name: str) -> list[str]:
    if atom_name not in mol2.atom_id_by_name:
        return []
    atom_id = mol2.atom_id_by_name[atom_name]
    result: list[str] = []
    for bond in mol2.bonds:
        if atom_id not in bond:
            continue
        other = next(iter(bond - {atom_id}))
        name, atom_type = mol2.atoms.get(other, ("", ""))
        if atom_element_from_mol2_name_type(name, atom_type) == "H":
            result.append(name)
    return sorted(result)


def audit_mol2(path: Path, donor_hint: str | None) -> list[dict[str, str]]:
    mol2 = read_mol2(path)
    if "ND1" not in mol2.atom_id_by_name or "NE2" not in mol2.atom_id_by_name:
        return []
    nd1_hs = bonded_hydrogens(mol2, "ND1")
    ne2_hs = bonded_hydrogens(mol2, "NE2")
    donor = donor_hint or "unknown"
    ring_h = 0
    for atom_id, (name, atom_type) in mol2.atoms.items():
        if atom_element_from_mol2_name_type(name, atom_type) == "H" and name not in set(nd1_hs + ne2_hs):
            ring_h += 1
    return [
        {
            "source": str(path),
            "kind": "mol2",
            "metal": "",
            "residue": path.stem,
            "donor": donor,
            "donor_distance_A": "",
            "Fe_ND1_A": "",
            "Fe_NE2_A": "",
            "expected_neutral_amber_name": expected_tautomer(donor),
            "ND1_H": ",".join(nd1_hs),
            "ND1_H_distance_A": "",
            "NE2_H": ",".join(ne2_hs),
            "NE2_H_distance_A": "",
            "ring_or_other_H_count": str(ring_h),
            "status": donor_status(donor, bool(nd1_hs), bool(ne2_hs), path.stem),
        }
    ]


def print_rows(rows: list[dict[str, str]]) -> None:
    columns = [
        "kind",
        "source",
        "metal",
        "residue",
        "donor",
        "donor_distance_A",
        "Fe_ND1_A",
        "Fe_NE2_A",
        "expected_neutral_amber_name",
        "ND1_H",
        "NE2_H",
        "ring_or_other_H_count",
        "status",
    ]
    for row in rows:
        print("\t".join(row.get(col, "") for col in columns))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, help="PDB to audit for Fe-bound histidines.")
    parser.add_argument("--mol2", type=Path, nargs="*", help="HD*/HID-like mol2 files to audit.")
    parser.add_argument("--metal-element", default="FE")
    parser.add_argument("--cutoff", type=float, default=2.8)
    parser.add_argument(
        "--mol2-donor",
        choices=["ND1", "NE2"],
        help="Donor atom expected for all mol2 files when no PDB geometry is available.",
    )
    parser.add_argument("--out-csv", type=Path)
    args = parser.parse_args()

    rows: list[dict[str, str]] = []
    if args.pdb:
        rows.extend(audit_pdb(args.pdb, args.metal_element, args.cutoff))
    for mol2_path in args.mol2 or []:
        rows.extend(audit_mol2(mol2_path, args.mol2_donor))

    print_rows(rows)
    if args.out_csv:
        args.out_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.out_csv.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()) if rows else ["status"])
            writer.writeheader()
            writer.writerows(rows)
    bad = [row for row in rows if row.get("status", "").startswith("warn")]
    return 2 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
