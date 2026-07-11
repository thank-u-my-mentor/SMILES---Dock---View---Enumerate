#!/usr/bin/env python3
"""Audit MCPB metal-site residue bonds from an Amber topology."""

from __future__ import annotations

import argparse
from collections import defaultdict

import parmed as pmd


EXPECTED_METAL_BONDS = [
    ("UL1", "N1", "FE1", "FE"),
    ("HD1", "NE2", "FE1", "FE"),
    ("HD2", "NE2", "FE1", "FE"),
    ("GU1", "OE1", "FE1", "FE"),
    ("AT1", "O2", "FE1", "FE"),
    ("HH1", "O", "FE1", "FE"),
]

EXPECTED_INTERNAL = {
    "HD1": [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB"), ("CB", "CG"), ("CG", "ND1"), ("CG", "CD2"), ("ND1", "CE1"), ("CE1", "NE2"), ("CD2", "NE2")],
    "HD2": [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB"), ("CB", "CG"), ("CG", "ND1"), ("CG", "CD2"), ("ND1", "CE1"), ("CE1", "NE2"), ("CD2", "NE2")],
    "GU1": [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "OE2")],
}


def atom_key(atom):
    return f"{atom.residue.idx + 1}:{atom.residue.name}:{atom.name}"


def distance(a, b):
    dx = a.xx - b.xx
    dy = a.xy - b.xy
    dz = a.xz - b.xz
    return (dx * dx + dy * dy + dz * dz) ** 0.5


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("prmtop")
    parser.add_argument("inpcrd")
    args = parser.parse_args()

    struct = pmd.load_file(args.prmtop, args.inpcrd)
    by_res = defaultdict(list)
    for atom in struct.atoms:
        by_res[atom.residue.name].append(atom)

    print(f"atoms={len(struct.atoms)} residues={len(struct.residues)} bonds={len(struct.bonds)}")

    res_summary = defaultdict(int)
    for res in struct.residues:
        res_summary[res.name] += 1
    for name in ["UL1", "HD1", "HD2", "GU1", "FE1", "AT1", "HH1"]:
        print(f"residue_count {name} {res_summary[name]}")

    bond_pairs = set()
    for bond in struct.bonds:
        a1, a2 = bond.atom1, bond.atom2
        bond_pairs.add(tuple(sorted((a1.idx, a2.idx))))

    def has_bond(a, b):
        return tuple(sorted((a.idx, b.idx))) in bond_pairs

    def find_atom(resname, atomname):
        hits = [a for a in by_res[resname] if a.name == atomname]
        if len(hits) != 1:
            raise RuntimeError(f"Expected one atom {resname}:{atomname}, found {len(hits)}")
        return hits[0]

    print("\n[metal bonds]")
    ok = True
    for r1, a1, r2, a2 in EXPECTED_METAL_BONDS:
        atom1 = find_atom(r1, a1)
        atom2 = find_atom(r2, a2)
        hb = has_bond(atom1, atom2)
        ok = ok and hb
        print(f"{atom_key(atom1)} -- {atom_key(atom2)} bond={hb} dist={distance(atom1, atom2):.3f} A")

    print("\n[internal residue bonds]")
    for resname, pairs in EXPECTED_INTERNAL.items():
        atoms = {a.name: a for a in by_res[resname]}
        for a1, a2 in pairs:
            atom1, atom2 = atoms[a1], atoms[a2]
            hb = has_bond(atom1, atom2)
            ok = ok and hb
            print(f"{resname}:{a1}-{a2} bond={hb} dist={distance(atom1, atom2):.3f} A")

    print("\n[mainchain external bonds]")
    for resname in ["HD1", "HD2", "GU1"]:
        atoms = {a.name: a for a in by_res[resname]}
        n_atom = atoms["N"]
        c_atom = atoms["C"]
        n_ext = [b for b in n_atom.bond_partners if b.residue is not n_atom.residue]
        c_ext = [b for b in c_atom.bond_partners if b.residue is not c_atom.residue]
        ok = ok and bool(n_ext) and bool(c_ext)
        print(f"{resname}: N external partners = {', '.join(atom_key(x) for x in n_ext) or 'NONE'}")
        print(f"{resname}: C external partners = {', '.join(atom_key(x) for x in c_ext) or 'NONE'}")

    print(f"\nOVERALL={'OK' if ok else 'FAIL'}")
    return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
