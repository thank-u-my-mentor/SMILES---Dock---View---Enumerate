#!/usr/bin/env python3
"""Merge MCPB protein/metal PDB with a GAFF mol2 ligand for tleap.

The ligand PDB is written from mol2 atom names and coordinates so tleap can
match the loaded LIG template exactly. Duplicate protein atom records are
removed by (chain, residue number, insertion code, residue name, atom name).
"""

from __future__ import annotations

import argparse
from pathlib import Path


TYPE_TO_ELEMENT = {
    "c": "C",
    "c1": "C",
    "c2": "C",
    "c3": "C",
    "ca": "C",
    "cc": "C",
    "cd": "C",
    "ce": "C",
    "cf": "C",
    "n": "N",
    "n1": "N",
    "n2": "N",
    "n3": "N",
    "n4": "N",
    "na": "N",
    "nb": "N",
    "nc": "N",
    "nd": "N",
    "ne": "N",
    "nf": "N",
    "nh": "N",
    "no": "N",
    "ns": "N",
    "o": "O",
    "oh": "O",
    "os": "O",
    "s": "S",
    "ss": "S",
    "sh": "S",
    "p": "P",
    "h": "H",
    "h1": "H",
    "h2": "H",
    "h3": "H",
    "h4": "H",
    "h5": "H",
    "ha": "H",
    "hc": "H",
    "hn": "H",
    "ho": "H",
    "hs": "H",
}


def element_from_type(atom_type: str, atom_name: str) -> str:
    key = atom_type.lower()
    if key in TYPE_TO_ELEMENT:
        return TYPE_TO_ELEMENT[key]
    letters = "".join(ch for ch in atom_name if ch.isalpha())
    return (letters[:1] or atom_type[:1]).upper()


def parse_mol2_atoms(path: Path):
    atoms = []
    in_atoms = False
    for line in path.read_text().splitlines():
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if line.startswith("@<TRIPOS>") and in_atoms:
            break
        if not in_atoms or not line.strip():
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        atoms.append(
            {
                "idx": int(parts[0]),
                "name": parts[1],
                "x": float(parts[2]),
                "y": float(parts[3]),
                "z": float(parts[4]),
                "type": parts[5],
                "elem": element_from_type(parts[5], parts[1]),
            }
        )
    if not atoms:
        raise SystemExit(f"No atoms parsed from {path}")
    return atoms


def pdb_atom_line(serial: int, name: str, resname: str, chain: str, resid: int, x: float, y: float, z: float, elem: str) -> str:
    # PDB atom name alignment: one-letter elements usually start in column 14.
    atom_name = f" {name:<3}" if len(elem) == 1 and len(name) < 4 else f"{name:<4}"
    return (
        f"HETATM{serial:5d} {atom_name}{resname:>4s} {chain:1s}{resid:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {elem:>2s}"
    )


def clean_pdb_lines(path: Path):
    seen = set()
    kept = []
    dropped = []
    max_serial = 0
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            atom_name = line[12:16]
            resname = line[17:20]
            atom_name_clean = atom_name.strip()
            resname_clean = resname.strip()
            if (resname_clean == "HIE" and atom_name_clean == "HD1") or (
                resname_clean == "HID" and atom_name_clean == "HE2"
            ):
                dropped.append(line)
                continue
            chain = line[21:22]
            resid = line[22:26]
            icode = line[26:27]
            key = (chain, resid, icode, resname, atom_name)
            if key in seen:
                dropped.append(line)
                continue
            seen.add(key)
            max_serial = max(max_serial, int(line[6:11]))
            kept.append(line)
        elif line.startswith(("TER", "END")):
            continue
        else:
            kept.append(line)
    return kept, dropped, max_serial


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mcpb-pdb", required=True, type=Path)
    ap.add_argument("--ligand-mol2", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--lig-resname", default="LIG")
    ap.add_argument("--lig-chain", default="A")
    ap.add_argument("--lig-resid", default=510, type=int)
    args = ap.parse_args()

    protein_lines, dropped, max_serial = clean_pdb_lines(args.mcpb_pdb)
    lig_atoms = parse_mol2_atoms(args.ligand_mol2)
    out_lines = list(protein_lines)
    if out_lines and not out_lines[-1].startswith("TER"):
        out_lines.append("TER")
    for offset, atom in enumerate(lig_atoms, start=1):
        out_lines.append(
            pdb_atom_line(
                max_serial + offset,
                atom["name"],
                args.lig_resname,
                args.lig_chain,
                args.lig_resid,
                atom["x"],
                atom["y"],
                atom["z"],
                atom["elem"],
            )
        )
    out_lines.extend(["TER", "END"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(out_lines) + "\n")
    print(f"wrote={args.out}")
    print(f"ligand_atoms={len(lig_atoms)}")
    print(f"dropped_duplicate_atoms={len(dropped)}")
    if dropped:
        dup_path = args.out.with_suffix(".dropped_duplicates.txt")
        dup_path.write_text("\n".join(dropped) + "\n")
        print(f"dropped_duplicates={dup_path}")


if __name__ == "__main__":
    main()
