#!/usr/bin/env python3
"""Audit MCPB small-model hydrogens before GAMESS OPT/Hessian.

This audit is deliberately conservative: it does not require every atom to have
hydrogen.  It only warns when atoms that should normally carry hydrogens in the
MCPB small model are missing them, or when non-standard ligand atom sets disagree
with their mol2 templates.  If warning_count > --max-warnings, the script exits
with code 2 so downstream OPT/Hessian scripts can hard-stop.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


HIS_NAMES = {"HIS", "HID", "HIE", "HIP", "HD1", "HD2", "HE1", "HE2", "HP1", "HP2"}
GLU_LIKE = {"GLU", "GLH", "GU1", "GU2", "GUX"}
ASP_LIKE = {"ASP", "ASH", "AD1", "AD2", "ASD"}
IGNORE_MOL2 = {"FE"}
ATOMIC_NUMBERS = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "S": 16,
    "P": 15,
    "F": 9,
    "CL": 17,
    "BR": 35,
    "I": 53,
    "FE": 26,
    "ZN": 30,
    "CU": 29,
    "MN": 25,
    "MG": 12,
    "CA": 20,
    "CO": 27,
    "NI": 28,
}


@dataclass(frozen=True)
class Atom:
    serial: int
    name: str
    resname: str
    chain: str
    resseq: str
    element: str
    xyz: tuple[float, float, float]

    @property
    def resid(self) -> tuple[str, str, str]:
        return (self.resname, self.chain, self.resseq)


def guess_element(atom_name: str, element_field: str = "") -> str:
    raw = "".join(ch for ch in element_field.strip() if ch.isalpha()).upper()
    if raw:
        return raw
    letters = "".join(ch for ch in atom_name.strip() if ch.isalpha()).upper()
    if letters.startswith("FE"):
        return "FE"
    if letters.startswith(("CL", "BR", "ZN", "MG", "MN", "CU", "CO", "NI", "CA")):
        return letters[:2]
    return letters[:1]


def read_pdb(path: Path) -> list[Atom]:
    atoms: list[Atom] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
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
                    element=guess_element(line[12:16].strip(), line[76:78] if len(line) >= 78 else ""),
                    xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
                )
            )
        except ValueError:
            continue
    if not atoms:
        raise SystemExit(f"No PDB atoms parsed from {path}")
    return atoms


def grouped_residues(atoms: list[Atom]) -> dict[tuple[str, str, str], list[Atom]]:
    residues: dict[tuple[str, str, str], list[Atom]] = {}
    for atom in atoms:
        residues.setdefault(atom.resid, []).append(atom)
    return residues


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def hydrogens_near(atom: Atom, residue_atoms: list[Atom], cutoff: float) -> list[Atom]:
    return [
        other
        for other in residue_atoms
        if other.element == "H" and distance(atom.xyz, other.xyz) <= cutoff
    ]


def add_count_warning(
    rows: list[dict[str, str]],
    *,
    source: Path,
    residue: tuple[str, str, str],
    atom: Atom,
    expected: int,
    observed: int,
    reason: str,
) -> None:
    if observed >= expected:
        rows.append(
            {
                "level": "ok",
                "source": str(source),
                "residue": f"{residue[0]}:{residue[1]}:{residue[2]}",
                "atom": atom.name,
                "expected_H": str(expected),
                "observed_H": str(observed),
                "status": "ok",
                "reason": reason,
            }
        )
        return
    rows.append(
        {
            "level": "warning",
            "source": str(source),
            "residue": f"{residue[0]}:{residue[1]}:{residue[2]}",
            "atom": atom.name,
            "expected_H": str(expected),
            "observed_H": str(observed),
            "status": "missing_expected_hydrogen",
            "reason": reason,
        }
    )


def nearest_metal_donor(atom_by_name: dict[str, Atom], metals: list[Atom]) -> str:
    if not metals or "ND1" not in atom_by_name or "NE2" not in atom_by_name:
        return ""
    nd1 = min(distance(atom_by_name["ND1"].xyz, metal.xyz) for metal in metals)
    ne2 = min(distance(atom_by_name["NE2"].xyz, metal.xyz) for metal in metals)
    return "ND1" if nd1 < ne2 else "NE2"


def audit_histidine(
    rows: list[dict[str, str]],
    *,
    source: Path,
    resid: tuple[str, str, str],
    residue_atoms: list[Atom],
    metals: list[Atom],
    h_cutoff: float,
) -> None:
    atom_by_name = {atom.name: atom for atom in residue_atoms}
    donor = nearest_metal_donor(atom_by_name, metals)
    nd1_h = hydrogens_near(atom_by_name["ND1"], residue_atoms, h_cutoff) if "ND1" in atom_by_name else []
    ne2_h = hydrogens_near(atom_by_name["NE2"], residue_atoms, h_cutoff) if "NE2" in atom_by_name else []
    if donor == "NE2":
        if not nd1_h:
            rows.append(
                {
                    "level": "warning",
                    "source": str(source),
                    "residue": f"{resid[0]}:{resid[1]}:{resid[2]}",
                    "atom": "ND1",
                    "expected_H": "1",
                    "observed_H": "0",
                    "status": "missing_opposite_imidazole_h",
                    "reason": "Fe-NE2 donor neutral HID should keep ND1-H; NE2-H is not expected",
                }
            )
        if ne2_h:
            rows.append(
                {
                    "level": "warning",
                    "source": str(source),
                    "residue": f"{resid[0]}:{resid[1]}:{resid[2]}",
                    "atom": "NE2",
                    "expected_H": "0",
                    "observed_H": str(len(ne2_h)),
                    "status": "fe_donor_n_is_protonated",
                    "reason": "Fe donor NE2 should normally be unprotonated",
                }
            )
    elif donor == "ND1":
        if not ne2_h:
            rows.append(
                {
                    "level": "warning",
                    "source": str(source),
                    "residue": f"{resid[0]}:{resid[1]}:{resid[2]}",
                    "atom": "NE2",
                    "expected_H": "1",
                    "observed_H": "0",
                    "status": "missing_opposite_imidazole_h",
                    "reason": "Fe-ND1 donor neutral HIE should keep NE2-H; ND1-H is not expected",
                }
            )
        if nd1_h:
            rows.append(
                {
                    "level": "warning",
                    "source": str(source),
                    "residue": f"{resid[0]}:{resid[1]}:{resid[2]}",
                    "atom": "ND1",
                    "expected_H": "0",
                    "observed_H": str(len(nd1_h)),
                    "status": "fe_donor_n_is_protonated",
                    "reason": "Fe donor ND1 should normally be unprotonated",
                }
            )

    for atom_name in ("CH3", "CB"):
        if atom_name in atom_by_name:
            expected = 3 if atom_name == "CH3" else 2
            observed = len(hydrogens_near(atom_by_name[atom_name], residue_atoms, h_cutoff))
            add_count_warning(
                rows,
                source=source,
                residue=resid,
                atom=atom_by_name[atom_name],
                expected=expected,
                observed=observed,
                reason=f"histidine small-model {atom_name} should retain normal aliphatic hydrogens",
            )

    for atom_name in ("CD2", "CE1"):
        if atom_name in atom_by_name:
            observed = len(hydrogens_near(atom_by_name[atom_name], residue_atoms, h_cutoff))
            add_count_warning(
                rows,
                source=source,
                residue=resid,
                atom=atom_by_name[atom_name],
                expected=1,
                observed=observed,
                reason=f"histidine imidazole carbon {atom_name} should normally carry one H",
            )


def audit_acidic_side_chain(
    rows: list[dict[str, str]],
    *,
    source: Path,
    resid: tuple[str, str, str],
    residue_atoms: list[Atom],
    h_cutoff: float,
) -> None:
    atom_by_name = {atom.name: atom for atom in residue_atoms}
    if "CH3" in atom_by_name:
        add_count_warning(
            rows,
            source=source,
            residue=resid,
            atom=atom_by_name["CH3"],
            expected=3,
            observed=len(hydrogens_near(atom_by_name["CH3"], residue_atoms, h_cutoff)),
            reason="MCPB methyl cap should have 3 hydrogens",
        )
    if "CB" in atom_by_name:
        add_count_warning(
            rows,
            source=source,
            residue=resid,
            atom=atom_by_name["CB"],
            expected=2,
            observed=len(hydrogens_near(atom_by_name["CB"], residue_atoms, h_cutoff)),
            reason="acidic side-chain CB should normally be CH2",
        )
    if resid[0].upper() in GLU_LIKE and "CG" in atom_by_name:
        add_count_warning(
            rows,
            source=source,
            residue=resid,
            atom=atom_by_name["CG"],
            expected=2,
            observed=len(hydrogens_near(atom_by_name["CG"], residue_atoms, h_cutoff)),
            reason="glutamate side-chain CG should normally be CH2",
        )


def parse_mol2_atom_names(path: Path) -> dict[str, str]:
    section = ""
    atoms: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("@<TRIPOS>"):
            section = line.strip()
            continue
        if section != "@<TRIPOS>ATOM":
            continue
        parts = line.split()
        if len(parts) >= 6 and parts[0].isdigit():
            atoms[parts[1]] = guess_element(parts[1], parts[5].split(".", 1)[0])
    return atoms


def audit_ligand_templates(
    rows: list[dict[str, str]],
    *,
    source: Path,
    residues: dict[tuple[str, str, str], list[Atom]],
    mol2_dir: Path,
) -> None:
    for resid, residue_atoms in residues.items():
        resname = resid[0]
        upper = resname.upper()
        if upper in HIS_NAMES or upper in GLU_LIKE or upper in ASP_LIKE or upper in {"HOH", "WAT"}:
            continue
        if upper in IGNORE_MOL2:
            continue
        mol2_path = mol2_dir / f"{resname}.mol2"
        if not mol2_path.exists():
            continue
        mol2_atoms = parse_mol2_atom_names(mol2_path)
        pdb_names = {atom.name for atom in residue_atoms}
        mol2_names = set(mol2_atoms)
        missing = sorted(mol2_names - pdb_names)
        extra = sorted(pdb_names - mol2_names)
        pdb_h = sorted(atom.name for atom in residue_atoms if atom.element == "H")
        mol2_h = sorted(name for name, element in mol2_atoms.items() if element == "H")
        if missing or extra:
            rows.append(
                {
                    "level": "warning",
                    "source": str(source),
                    "residue": f"{resid[0]}:{resid[1]}:{resid[2]}",
                    "atom": "*",
                    "expected_H": ",".join(mol2_h),
                    "observed_H": ",".join(pdb_h),
                    "status": "ligand_atom_set_differs_from_mol2",
                    "reason": f"missing_from_pdb={','.join(missing)}; extra_in_pdb={','.join(extra)}",
                }
            )
        else:
            rows.append(
                {
                    "level": "ok",
                    "source": str(source),
                    "residue": f"{resid[0]}:{resid[1]}:{resid[2]}",
                    "atom": "*",
                    "expected_H": ",".join(mol2_h),
                    "observed_H": ",".join(pdb_h),
                    "status": "ok",
                    "reason": f"ligand atom names match {mol2_path.name}",
                }
            )


def print_rows(rows: list[dict[str, str]]) -> None:
    columns = ["level", "source", "residue", "atom", "expected_H", "observed_H", "status", "reason"]
    print("\t".join(columns))
    for row in rows:
        print("\t".join(row.get(col, "") for col in columns))


def electron_count(atoms: list[Atom], charge: int) -> int:
    total = 0
    for atom in atoms:
        elem = atom.element.upper()
        if elem not in ATOMIC_NUMBERS:
            raise SystemExit(f"Unknown element {elem!r} for atom {atom.name} in {atom.resname}:{atom.chain}:{atom.resseq}")
        total += ATOMIC_NUMBERS[elem]
    return total - charge


def multiplicity_compatible(electrons: int, mult: int) -> bool:
    return (electrons % 2 == 0 and mult % 2 == 1) or (electrons % 2 == 1 and mult % 2 == 0)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--small-pdb", type=Path, required=True)
    parser.add_argument("--mol2-dir", type=Path, help="Directory containing ACT/UNL/UNK mol2 templates for ligand atom-name checks.")
    parser.add_argument("--out", type=Path, help="Optional TSV output path.")
    parser.add_argument("--max-warnings", type=int, default=2, help="Exit with error if warning_count is greater than this value.")
    parser.add_argument("--h-cutoff", type=float, default=1.25, help="Distance cutoff for assigning H to a heavy atom.")
    parser.add_argument("--model-charge", type=int, help="Optional MCPB small-model total charge for electron/multiplicity parity audit.")
    parser.add_argument("--mult", type=int, help="Optional MCPB small-model spin multiplicity for electron/multiplicity parity audit.")
    args = parser.parse_args()

    atoms = read_pdb(args.small_pdb)
    residues = grouped_residues(atoms)
    metals = [atom for atom in atoms if atom.element == "FE" or atom.name.upper() == "FE"]

    rows: list[dict[str, str]] = []
    for resid, residue_atoms in residues.items():
        resname = resid[0].upper()
        if resname in HIS_NAMES:
            audit_histidine(rows, source=args.small_pdb, resid=resid, residue_atoms=residue_atoms, metals=metals, h_cutoff=args.h_cutoff)
        elif resname in GLU_LIKE or resname in ASP_LIKE:
            audit_acidic_side_chain(rows, source=args.small_pdb, resid=resid, residue_atoms=residue_atoms, h_cutoff=args.h_cutoff)

    if args.mol2_dir:
        audit_ligand_templates(rows, source=args.small_pdb, residues=residues, mol2_dir=args.mol2_dir)

    print_rows(rows)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["level", "source", "residue", "atom", "expected_H", "observed_H", "status", "reason"], delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    warnings = [row for row in rows if row.get("level") == "warning"]
    print(f"warning_count={len(warnings)} max_warnings={args.max_warnings}")
    parity_ok = True
    if args.model_charge is not None or args.mult is not None:
        if args.model_charge is None or args.mult is None:
            raise SystemExit("--model-charge and --mult must be provided together")
        electrons = electron_count(atoms, args.model_charge)
        parity_ok = multiplicity_compatible(electrons, args.mult)
        print(f"electron_count={electrons} model_charge={args.model_charge} mult={args.mult} parity_ok={str(parity_ok).lower()}")
        if not parity_ok:
            print("verdict=error_electron_multiplicity_parity_mismatch")
            return 3
    if len(warnings) > args.max_warnings:
        print("verdict=error_too_many_mcpb_small_model_hydrogen_warnings")
        return 2
    print("verdict=ok_or_review")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
