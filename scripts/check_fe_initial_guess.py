from __future__ import annotations

import math
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Atom:
    serial: int
    record: str
    name: str
    resn: str
    chain: str
    resi: str
    x: float
    y: float
    z: float
    elem: str


def dist(a: Atom, b: Atom) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def parse_pdb(path: Path) -> tuple[list[Atom], list[tuple[int, int]]]:
    atoms: list[Atom] = []
    bonds: set[tuple[int, int]] = set()
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            raw_elem = line[76:80].strip() if len(line) >= 78 else ""
            elem = "".join(ch for ch in raw_elem if ch.isalpha())
            if not elem:
                name = line[12:16].strip()
                elem = "".join(ch for ch in name if ch.isalpha())[:1]
            atoms.append(
                Atom(
                    serial=int(line[6:11]),
                    record=line[0:6].strip(),
                    name=line[12:16].strip(),
                    resn=line[17:20].strip(),
                    chain=line[21].strip(),
                    resi=line[22:26].strip(),
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    elem=elem.capitalize(),
                )
            )
        elif line.startswith("CONECT"):
            parts = line.split()
            if len(parts) >= 3:
                src = int(parts[1])
                for token in parts[2:]:
                    try:
                        dst = int(token)
                    except ValueError:
                        continue
                    bonds.add(tuple(sorted((src, dst))))
    return atoms, sorted(bonds)


def atom_label(atom: Atom) -> str:
    chain = atom.chain or "-"
    return f"{atom.serial}:{atom.resn}:{chain}:{atom.resi}:{atom.name}:{atom.elem}"


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: check_fe_initial_guess.py input.pdb")
        return 2

    atoms, bonds = parse_pdb(Path(sys.argv[1]))
    by_serial = {atom.serial: atom for atom in atoms}
    metals = [a for a in atoms if a.elem.upper() == "FE" or a.resn.upper() in {"FE", "FE2", "FE3"}]

    print(f"atoms={len(atoms)} bonds={len(bonds)} fe_atoms={len(metals)}")
    print()

    donor_elems = {"N", "O", "S"}
    for metal in metals:
        print(f"FE_SITE {atom_label(metal)}")
        neighbors = [
            (dist(metal, atom), atom)
            for atom in atoms
            if atom.serial != metal.serial and atom.elem in donor_elems and dist(metal, atom) <= 3.2
        ]
        neighbors.sort(key=lambda row: row[0])
        for d, atom in neighbors:
            print(f"  {d:6.3f} A  {atom_label(atom)}")
        if not neighbors:
            print("  no N/O/S donors within 3.2 A")
        print()

    residues = defaultdict(list)
    for atom in atoms:
        residues[(atom.resn, atom.chain or "-", atom.resi)].append(atom)

    print("LIGAND_INTERNAL_BONDS")
    for res_key in sorted(k for k in residues if k[0] in {"ACE", "ACT", "UNL", "LIG", "CAR", "NHB"}):
        res_atoms = {a.serial for a in residues[res_key]}
        print(f"RES {res_key[0]} chain={res_key[1]} resi={res_key[2]} atoms={len(res_atoms)}")
        local_bonds = [(a, b) for a, b in bonds if a in res_atoms and b in res_atoms]
        if not local_bonds:
            print("  no CONECT bonds")
            continue
        for a_serial, b_serial in local_bonds:
            a = by_serial[a_serial]
            b = by_serial[b_serial]
            d = dist(a, b)
            flag = "  !!long" if d > 2.1 and not {a.elem, b.elem} == {"Fe"} else ""
            print(f"  {d:6.3f} A  {atom_label(a)} -- {atom_label(b)}{flag}")
        print()

    print("FORMAT_WARNINGS")
    for atom in atoms:
        if atom.resn == "ACE":
            print(
                "  ACE residue name is an Amber N-terminal acetyl cap name; "
                "rename ligand acetate/carboxylate before Amber/tleap."
            )
            break
    for atom in atoms:
        if len(atom.elem) > 2:
            print(f"  suspicious element field: {atom_label(atom)}")
    if len(metals) > 1:
        print("  multiple Fe atoms found; later QM/MCPB/MD setup must choose chain A or B explicitly.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
