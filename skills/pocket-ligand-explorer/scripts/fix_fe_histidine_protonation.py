#!/usr/bin/env python3
"""Repair Fe-bound histidine tautomers before MCPB.py step 1.

Rule used by default:

* If Fe binds NE2, rename the residue to HID, keep NE2 unprotonated, and add
  ND1-H if missing.
* If Fe binds ND1, rename the residue to HIE, keep ND1 unprotonated, and add
  NE2-H if missing.

The script also adds missing imidazole C-H atoms and common side-chain/backbone
hydrogens for the repaired histidines.  This is a geometry preflight helper for
MCPB/Amber handoffs; always inspect the output PDB and report before long QM.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, replace
from pathlib import Path


HIS_NAMES = {"HIS", "HID", "HIE", "HIP", "HD1", "HD2", "HE1", "HE2"}
METAL_ELEMENTS = {"FE", "ZN", "CU", "MN", "MG", "CO", "NI"}


@dataclass
class Atom:
    index: int
    record: str
    serial: int
    name: str
    resname: str
    chain: str
    resseq: int
    icode: str
    element: str
    xyz: tuple[float, float, float]
    line: str

    @property
    def resid(self) -> tuple[str, int, str]:
        return (self.chain, self.resseq, self.icode)


def infer_element(name: str, element_field: str = "") -> str:
    raw = element_field.strip().upper().replace("+", "").replace("-", "")
    if raw:
        return raw
    letters = "".join(ch for ch in name.strip() if ch.isalpha()).upper()
    if letters[:2] in METAL_ELEMENTS:
        return letters[:2]
    return letters[:1]


def parse_pdb(path: Path) -> tuple[list[str], list[Atom]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
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
    return lines, atoms


def fmt_atom(atom: Atom) -> str:
    altloc = " "
    chain = atom.chain[:1] if atom.chain else " "
    icode = atom.icode[:1] if atom.icode else " "
    element = (atom.element or infer_element(atom.name)).upper()
    return (
        f"{atom.record:<6}{atom.serial:5d} {atom.name:>4s}{altloc}"
        f"{atom.resname:>3s} {chain}{atom.resseq:4d}{icode}   "
        f"{atom.xyz[0]:8.3f}{atom.xyz[1]:8.3f}{atom.xyz[2]:8.3f}"
        f"  1.00  0.00          {element:>2s}  "
    )


def v_add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def v_sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def v_scale(a, s: float):
    return (a[0] * s, a[1] * s, a[2] * s)


def dot(a, b) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(a) -> float:
    return math.sqrt(dot(a, a))


def unit(a):
    n = norm(a)
    if n < 1e-8:
        return (1.0, 0.0, 0.0)
    return (a[0] / n, a[1] / n, a[2] / n)


def dist(a: Atom, b: Atom) -> float:
    return norm(v_sub(a.xyz, b.xyz))


def perpendicular_basis(axis):
    axis = unit(axis)
    trial = (1.0, 0.0, 0.0) if abs(axis[0]) < 0.85 else (0.0, 1.0, 0.0)
    p = unit(cross(axis, trial))
    q = unit(cross(axis, p))
    return p, q


def single_h(parent: Atom, neighbors: list[Atom], *, name: str, length: float = 1.01, serial: int) -> Atom:
    direction = (0.0, 0.0, 0.0)
    for neighbor in neighbors:
        direction = v_add(direction, unit(v_sub(parent.xyz, neighbor.xyz)))
    direction = unit(direction)
    return replace(
        parent,
        serial=serial,
        name=name,
        element="H",
        xyz=v_add(parent.xyz, v_scale(direction, length)),
        record="ATOM",
    )


def two_h(parent: Atom, neighbors: list[Atom], *, names: tuple[str, str], length: float, serial: int) -> list[Atom]:
    away = (0.0, 0.0, 0.0)
    axis_refs = []
    for neighbor in neighbors:
        vec = unit(v_sub(parent.xyz, neighbor.xyz))
        away = v_add(away, vec)
        axis_refs.append(vec)
    away = unit(away)
    if len(axis_refs) >= 2:
        perp = unit(cross(axis_refs[0], axis_refs[1]))
        if norm(perp) < 1e-8:
            perp, _ = perpendicular_basis(away)
    else:
        perp, _ = perpendicular_basis(away)
    dirs = [unit(v_add(v_scale(away, 0.72), v_scale(perp, 0.69))), unit(v_add(v_scale(away, 0.72), v_scale(perp, -0.69)))]
    return [
        replace(parent, serial=serial, name=names[0], element="H", xyz=v_add(parent.xyz, v_scale(dirs[0], length)), record="ATOM"),
        replace(parent, serial=serial + 1, name=names[1], element="H", xyz=v_add(parent.xyz, v_scale(dirs[1], length)), record="ATOM"),
    ]


def group_residues(atoms: list[Atom]) -> dict[tuple[str, int, str], list[Atom]]:
    residues: dict[tuple[str, int, str], list[Atom]] = {}
    for atom in atoms:
        residues.setdefault(atom.resid, []).append(atom)
    return residues


def h_near(atom: Atom, atoms: list[Atom], cutoff: float = 1.25) -> list[Atom]:
    return [other for other in atoms if other.element == "H" and norm(v_sub(atom.xyz, other.xyz)) <= cutoff]


def has_atom(atoms_by_name: dict[str, Atom], name: str) -> bool:
    return name in atoms_by_name


def repair_histidines(lines: list[str], atoms: list[Atom], *, metal_element: str, cutoff: float, add_backbone_h: bool) -> tuple[list[str], list[dict[str, str]]]:
    residues = group_residues(atoms)
    metals = [atom for atom in atoms if atom.element.upper() == metal_element.upper() or atom.name.upper() == metal_element.upper()]
    max_serial = max((atom.serial for atom in atoms), default=0)
    replacements: dict[int, str | None] = {}
    insert_after: dict[int, list[Atom]] = {}
    report: list[dict[str, str]] = []

    for resid, res_atoms in sorted(residues.items(), key=lambda item: item[0]):
        resname = res_atoms[0].resname.upper()
        if resname not in HIS_NAMES:
            continue
        by = {atom.name: atom for atom in res_atoms}
        if "ND1" not in by or "NE2" not in by:
            continue
        nearest: tuple[float, Atom, str] | None = None
        for metal in metals:
            d_nd1 = dist(metal, by["ND1"])
            d_ne2 = dist(metal, by["NE2"])
            donor = "ND1" if d_nd1 <= d_ne2 else "NE2"
            d = min(d_nd1, d_ne2)
            if d <= cutoff and (nearest is None or d < nearest[0]):
                nearest = (d, metal, donor)
        if nearest is None:
            continue

        donor_dist, metal, donor = nearest
        expected = "HID" if donor == "NE2" else "HIE"
        acceptor_h_name = "HD1" if expected == "HID" else "HE2"
        donor_h_name = "HE2" if donor == "NE2" else "HD1"
        acceptor_n = by["ND1"] if expected == "HID" else by["NE2"]
        donor_n = by[donor]
        added: list[str] = []
        removed: list[str] = []

        # Rename the whole residue to the expected neutral Amber tautomer.
        for atom in res_atoms:
            new_atom = replace(atom, resname=expected)
            replacements[atom.index] = fmt_atom(new_atom)

        # Remove hydrogens attached to the Fe-donor nitrogen.
        for h in h_near(donor_n, res_atoms):
            if h.name == donor_h_name or h.name.startswith("H"):
                replacements[h.index] = None
                removed.append(h.name)

        def add_after(parent_name: str, new_atom: Atom) -> None:
            parent = by[parent_name]
            # Hydrogens are constructed from the original atom records, whose
            # residue name may still be HIS.  Keep inserted atoms in the same
            # Amber tautomer residue as the renamed heavy atoms.
            insert_after.setdefault(parent.index, []).append(replace(new_atom, resname=expected))

        if not h_near(acceptor_n, res_atoms) and acceptor_h_name not in by:
            max_serial += 1
            neighbors = [by[name] for name in (("CG", "CE1") if acceptor_n.name == "ND1" else ("CD2", "CE1")) if name in by]
            add_after(acceptor_n.name, single_h(acceptor_n, neighbors, name=acceptor_h_name, serial=max_serial))
            added.append(acceptor_h_name)

        for parent_name, h_name, neighbor_names in [
            ("CD2", "HD2", ("CG", "NE2")),
            ("CE1", "HE1", ("ND1", "NE2")),
        ]:
            if parent_name in by and h_name not in by and not h_near(by[parent_name], res_atoms, cutoff=1.20):
                max_serial += 1
                add_after(parent_name, single_h(by[parent_name], [by[n] for n in neighbor_names if n in by], name=h_name, length=1.08, serial=max_serial))
                added.append(h_name)

        if "CB" in by and "CA" in by and "CG" in by:
            existing_hb = [atom for atom in res_atoms if atom.name in {"HB2", "HB3", "1HB", "2HB"} or (atom.element == "H" and norm(v_sub(atom.xyz, by["CB"].xyz)) <= 1.20)]
            if len(existing_hb) < 2:
                h_names = ("HB2", "HB3")
                new_hs = two_h(by["CB"], [by["CA"], by["CG"]], names=h_names, length=1.09, serial=max_serial + 1)
                max_serial += 2
                add_after("CB", new_hs[0])
                add_after("CB", new_hs[1])
                added.extend(h_names)

        if "CA" in by and "HA" not in by and {"N", "C", "CB"}.issubset(by):
            max_serial += 1
            add_after("CA", single_h(by["CA"], [by["N"], by["C"], by["CB"]], name="HA", length=1.09, serial=max_serial))
            added.append("HA")

        if add_backbone_h and "N" in by and "H" not in by:
            prev_c = None
            prev_res = residues.get((resid[0], resid[1] - 1, ""))
            if prev_res:
                prev_by = {atom.name: atom for atom in prev_res}
                prev_c = prev_by.get("C")
            neighbors = [by["CA"]]
            if prev_c is not None:
                neighbors.append(prev_c)
            max_serial += 1
            add_after("N", single_h(by["N"], neighbors, name="H", length=1.01, serial=max_serial))
            added.append("H")

        report.append(
            {
                "residue": f"{res_atoms[0].resname}:{resid[0]}:{resid[1]}{resid[2]}",
                "expected_resname": expected,
                "metal": f"{metal.name}:{metal.chain}:{metal.resseq}",
                "donor": donor,
                "donor_distance_A": f"{donor_dist:.3f}",
                "added_atoms": ",".join(added),
                "removed_atoms": ",".join(removed),
                "status": "repaired",
            }
        )

    output: list[str] = []
    for idx, line in enumerate(lines):
        if idx in replacements:
            if replacements[idx] is not None:
                output.append(replacements[idx] or line)
            # None means the atom was removed.
        else:
            output.append(line)
        if idx in insert_after:
            output.extend(fmt_atom(atom) for atom in insert_after[idx])
    return output, report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--metal-element", default="FE")
    parser.add_argument("--cutoff", type=float, default=2.8)
    parser.add_argument("--no-backbone-h", action="store_true")
    args = parser.parse_args()

    lines, atoms = parse_pdb(args.pdb)
    output, report = repair_histidines(
        lines,
        atoms,
        metal_element=args.metal_element,
        cutoff=args.cutoff,
        add_backbone_h=not args.no_backbone_h,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(output).rstrip() + "\n", encoding="utf-8", newline="\n")
    report_path = args.report or args.out.with_suffix(args.out.suffix + ".fe_his_repair.tsv")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = ["residue", "expected_resname", "metal", "donor", "donor_distance_A", "added_atoms", "removed_atoms", "status"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(report)
    print(f"out={args.out}")
    print(f"report={report_path}")
    print(f"repaired_histidines={len(report)}")
    return 0 if report else 1


if __name__ == "__main__":
    raise SystemExit(main())
