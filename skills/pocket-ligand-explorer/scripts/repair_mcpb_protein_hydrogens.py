#!/usr/bin/env python3
"""Protein-only hydrogen repair before MCPB.py step 1.

This helper is intentionally narrow.  It repairs standard protein residues in
the Fe first shell before MCPB truncates/caps them, and it never edits ligand or
radical residues such as UNL/UNK/ACT.  Non-standard residues should be curated
manually and audited against their mol2 templates instead.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import replace
from pathlib import Path

from fix_fe_histidine_protonation import Atom
from fix_fe_histidine_protonation import fmt_atom
from fix_fe_histidine_protonation import infer_element
from fix_fe_histidine_protonation import repair_histidines
from fix_fe_histidine_protonation import single_h
from fix_fe_histidine_protonation import two_h


PROTEIN_RECORD = "ATOM"
HIS_NAMES = {"HIS", "HID", "HIE", "HIP", "HD1", "HD2", "HE1", "HE2"}
ACIDIC_NAMES = {"ASP", "ASH", "GLU", "GLH"}
DONOR_ELEMENTS = {"N", "O", "S"}


def parse_pdb_lines(lines: list[str]) -> list[Atom]:
    atoms: list[Atom] = []
    for idx, line in enumerate(lines):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            atoms.append(
                Atom(
                    index=idx,
                    record=line[:6].strip() or "ATOM",
                    serial=int(line[6:11]),
                    name=line[12:16].strip(),
                    resname=line[17:20].strip(),
                    chain=line[21:22].strip(),
                    resseq=int(line[22:26]),
                    icode=line[26:27].strip(),
                    element=infer_element(line[12:16], line[76:78] if len(line) >= 78 else ""),
                    xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
                    line=(line + " " * 80)[:80],
                )
            )
        except ValueError:
            continue
    return atoms


def v_sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def norm(v: tuple[float, float, float]) -> float:
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def distance(a: Atom, b: Atom) -> float:
    return norm(v_sub(a.xyz, b.xyz))


def group_residues(atoms: list[Atom]) -> dict[tuple[str, int, str], list[Atom]]:
    residues: dict[tuple[str, int, str], list[Atom]] = {}
    for atom in atoms:
        residues.setdefault(atom.resid, []).append(atom)
    return residues


def h_near(atom: Atom, atoms: list[Atom], cutoff: float = 1.25) -> list[Atom]:
    return [other for other in atoms if other.element == "H" and norm(v_sub(other.xyz, atom.xyz)) <= cutoff]


def parse_residue_spec(spec: str) -> tuple[str, int, str]:
    parts = spec.split(":")
    if len(parts) == 2:
        chain, resseq = parts
        icode = ""
    elif len(parts) == 3:
        chain, resseq, icode = parts
    else:
        raise argparse.ArgumentTypeError(f"Residue must look like A:349 or A:349:A, got {spec!r}")
    try:
        return (chain, int(resseq), icode)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Bad residue number in {spec!r}") from exc


def select_metal_site_residues(
    residues: dict[tuple[str, int, str], list[Atom]],
    atoms: list[Atom],
    *,
    metal_element: str,
    cutoff: float,
) -> set[tuple[str, int, str]]:
    metals = [
        atom
        for atom in atoms
        if atom.element.upper() == metal_element.upper() or atom.name.upper() == metal_element.upper()
    ]
    selected: set[tuple[str, int, str]] = set()
    for resid, residue_atoms in residues.items():
        if not residue_atoms or residue_atoms[0].record != PROTEIN_RECORD:
            continue
        if residue_atoms[0].resname.upper() not in HIS_NAMES | ACIDIC_NAMES:
            continue
        for atom in residue_atoms:
            if atom.element.upper() not in DONOR_ELEMENTS:
                continue
            if any(distance(atom, metal) <= cutoff for metal in metals):
                selected.add(resid)
                break
    return selected


def add_after(insert_after: dict[int, list[Atom]], parent: Atom, atom: Atom) -> None:
    insert_after.setdefault(parent.index, []).append(atom)


def repair_acidic_residue(
    residue_atoms: list[Atom],
    *,
    serial_start: int,
    insert_after: dict[int, list[Atom]],
) -> tuple[int, list[str]]:
    by = {atom.name: atom for atom in residue_atoms}
    added: list[str] = []
    serial = serial_start

    def add(parent_name: str, atom: Atom) -> None:
        nonlocal serial
        add_after(insert_after, by[parent_name], atom)
        added.append(atom.name)
        serial = max(serial, atom.serial)

    if "CA" in by and "HA" not in by and {"N", "C", "CB"}.issubset(by):
        serial += 1
        add("CA", single_h(by["CA"], [by["N"], by["C"], by["CB"]], name="HA", length=1.09, serial=serial))

    if "CB" in by and {"CA", "CG"}.issubset(by):
        existing = [atom for atom in residue_atoms if atom.element == "H" and norm(v_sub(atom.xyz, by["CB"].xyz)) <= 1.20]
        if len(existing) < 2:
            new_h = two_h(by["CB"], [by["CA"], by["CG"]], names=("HB2", "HB3"), length=1.09, serial=serial + 1)
            serial += 2
            add("CB", new_h[0])
            add("CB", new_h[1])

    if residue_atoms[0].resname.upper() in {"GLU", "GLH"} and "CG" in by and {"CB", "CD"}.issubset(by):
        existing = [atom for atom in residue_atoms if atom.element == "H" and norm(v_sub(atom.xyz, by["CG"].xyz)) <= 1.20]
        if len(existing) < 2:
            new_h = two_h(by["CG"], [by["CB"], by["CD"]], names=("HG2", "HG3"), length=1.09, serial=serial + 1)
            serial += 2
            add("CG", new_h[0])
            add("CG", new_h[1])

    return serial, added


def repair_protein_hydrogens(
    lines: list[str],
    *,
    metal_element: str = "FE",
    cutoff: float = 2.8,
    residues_to_repair: set[tuple[str, int, str]] | None = None,
    all_protein: bool = False,
    fix_fe_his: bool = True,
) -> tuple[list[str], list[dict[str, str]]]:
    report: list[dict[str, str]] = []

    if fix_fe_his:
        atoms_for_his = parse_pdb_lines(lines)
        lines, his_report = repair_histidines(
            lines,
            atoms_for_his,
            metal_element=metal_element,
            cutoff=cutoff,
            add_backbone_h=True,
        )
        for row in his_report:
            row = dict(row)
            row["repair_scope"] = "protein_fe_histidine"
            report.append(row)

    atoms = parse_pdb_lines(lines)
    residues = group_residues(atoms)
    selected = set(residues_to_repair or set())
    if all_protein:
        selected.update(
            resid
            for resid, residue_atoms in residues.items()
            if residue_atoms and residue_atoms[0].record == PROTEIN_RECORD and residue_atoms[0].resname.upper() in ACIDIC_NAMES
        )
    else:
        selected.update(select_metal_site_residues(residues, atoms, metal_element=metal_element, cutoff=cutoff))

    max_serial = max((atom.serial for atom in atoms), default=0)
    insert_after: dict[int, list[Atom]] = {}
    acid_rows: list[dict[str, str]] = []
    for resid in sorted(selected, key=lambda item: (item[0], item[1], item[2])):
        residue_atoms = residues.get(resid, [])
        if not residue_atoms:
            continue
        first = residue_atoms[0]
        if first.record != PROTEIN_RECORD:
            acid_rows.append(
                {
                    "repair_scope": "skip_non_protein",
                    "residue": f"{first.resname}:{first.chain}:{first.resseq}{first.icode}",
                    "status": "skipped",
                    "added_atoms": "",
                    "note": "HETATM/nonstandard residue is never auto-repaired",
                }
            )
            continue
        if first.resname.upper() not in ACIDIC_NAMES:
            continue
        max_serial, added = repair_acidic_residue(residue_atoms, serial_start=max_serial, insert_after=insert_after)
        acid_rows.append(
            {
                "repair_scope": "protein_acidic_sidechain",
                "residue": f"{first.resname}:{first.chain}:{first.resseq}{first.icode}",
                "status": "repaired" if added else "already_ok_or_no_action",
                "added_atoms": ",".join(added),
                "note": "Only CA/CB/CG aliphatic H atoms are repaired; carboxylate O atoms are not protonated automatically.",
            }
        )

    output: list[str] = []
    for idx, line in enumerate(lines):
        output.append(line)
        if idx in insert_after:
            output.extend(fmt_atom(replace(atom, record=PROTEIN_RECORD)) for atom in insert_after[idx])

    report.extend(acid_rows)
    return output, report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--metal-element", default="FE")
    parser.add_argument("--cutoff", type=float, default=2.8)
    parser.add_argument("--residue", action="append", type=parse_residue_spec, default=[])
    parser.add_argument("--all-protein", action="store_true", help="Repair all standard protein GLU/ASP residues instead of only the metal site.")
    parser.add_argument("--no-fix-fe-his", action="store_true", help="Skip Fe-bound histidine tautomer repair.")
    args = parser.parse_args()

    lines = args.pdb.read_text(encoding="utf-8", errors="replace").splitlines()
    output, report = repair_protein_hydrogens(
        lines,
        metal_element=args.metal_element,
        cutoff=args.cutoff,
        residues_to_repair=set(args.residue),
        all_protein=args.all_protein,
        fix_fe_his=not args.no_fix_fe_his,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(output).rstrip() + "\n")

    report_path = args.report or args.out.with_suffix(args.out.suffix + ".protein_h_repair.tsv")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["repair_scope", "residue", "expected_resname", "metal", "donor", "donor_distance_A", "added_atoms", "removed_atoms", "status", "note"]
    with report_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(report)

    print(f"out={args.out}")
    print(f"report={report_path}")
    print(f"repair_rows={len(report)}")
    print("ligand_policy=read_only_no_auto_hydrogen_repair")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
