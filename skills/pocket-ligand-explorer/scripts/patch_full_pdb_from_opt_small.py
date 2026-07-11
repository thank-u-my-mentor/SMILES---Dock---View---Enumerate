#!/usr/bin/env python3
"""Patch a full MCPB PDB preview from an optimized small-model PDB.

This is a pre-MD audit helper. It does not replace MCPB.py step 4 output.
It helps catch bad coordinate grafts before tleap/GROMACS work starts.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path


Vec = tuple[float, float, float]
Key = tuple[str, str, str, str]


@dataclass
class PdbAtom:
    line_index: int
    record: str
    atom: str
    resn: str
    chain: str
    resi: str
    xyz: Vec
    element: str

    @property
    def key(self) -> Key:
        return (self.resn, self.chain, self.resi, self.atom)

    @property
    def label(self) -> str:
        return f"{self.resn}:{self.chain}:{self.resi}:{self.atom}"


def distance(a: Vec, b: Vec) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def parse_pdb(path: Path) -> tuple[list[str], dict[Key, PdbAtom]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    atoms: dict[Key, PdbAtom] = {}
    for i, line in enumerate(lines):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        atom = PdbAtom(
            line_index=i,
            record=line[:6].strip(),
            atom=line[12:16].strip(),
            resn=line[17:20].strip(),
            chain=line[21:22].strip(),
            resi=line[22:26].strip(),
            xyz=xyz,
            element=(line[76:78].strip() or line[12:16].strip()[0]).upper(),
        )
        atoms[atom.key] = atom
    return lines, atoms


def replace_xyz(line: str, xyz: Vec) -> str:
    padded = (line + " " * 80)[:80]
    return f"{padded[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{padded[54:]}"


def write_patched_full_pdb(full_lines: list[str], full_atoms: dict[Key, PdbAtom], opt_atoms: dict[Key, PdbAtom], out: Path) -> list[tuple[Key, Vec, Vec, float]]:
    patch_by_index: dict[int, Vec] = {}
    rows: list[tuple[Key, Vec, Vec, float]] = []
    for key, opt_atom in opt_atoms.items():
        full_atom = full_atoms.get(key)
        if full_atom is None:
            continue
        patch_by_index[full_atom.line_index] = opt_atom.xyz
        rows.append((key, full_atom.xyz, opt_atom.xyz, distance(full_atom.xyz, opt_atom.xyz)))

    out_lines: list[str] = []
    for i, line in enumerate(full_lines):
        xyz = patch_by_index.get(i)
        out_lines.append(replace_xyz(line, xyz) if xyz else line)
    out.write_text("\n".join(out_lines).rstrip() + "\n", encoding="utf-8")
    return sorted(rows)


def atom_by(parts: tuple[str, str, str, str], atoms: dict[Key, PdbAtom]) -> PdbAtom | None:
    return atoms.get(parts)


def bond_rows(full_atoms: dict[Key, PdbAtom], patched_atoms: dict[Key, PdbAtom]) -> list[str]:
    checks: list[tuple[str, Key, Key, float, float]] = []

    def add(label: str, resn: str, chain: str, resi: str, a: str, b: str, lo: float, hi: float) -> None:
        checks.append((label, (resn, chain, resi, a), (resn, chain, resi, b), lo, hi))

    for resi in ("187", "270"):
        for a, b, lo, hi in [
            ("N", "CA", 1.25, 1.65),
            ("CA", "C", 1.25, 1.70),
            ("C", "O", 1.05, 1.40),
            ("CA", "CB", 1.35, 1.75),
            ("CB", "CG", 1.35, 1.75),
            ("CG", "ND1", 1.20, 1.55),
            ("CG", "CD2", 1.20, 1.55),
            ("ND1", "CE1", 1.20, 1.55),
            ("CE1", "NE2", 1.20, 1.55),
            ("CD2", "NE2", 1.20, 1.55),
        ]:
            add(f"HID{resi}:{a}-{b}", "HID", "A", resi, a, b, lo, hi)

    for a, b, lo, hi in [
        ("N", "CA", 1.25, 1.65),
        ("CA", "C", 1.25, 1.70),
        ("C", "O", 1.05, 1.40),
        ("CA", "CB", 1.35, 1.75),
        ("CB", "CG", 1.35, 1.75),
        ("CG", "CD", 1.35, 1.75),
        ("CD", "OE1", 1.15, 1.40),
        ("CD", "OE2", 1.15, 1.40),
    ]:
        add(f"GLU349:{a}-{b}", "GLU", "A", "349", a, b, lo, hi)

    rows = ["label\tinitial_A\tpreview_A\tdelta_A\tallowed_A\tflag"]
    for label, a_key, b_key, lo, hi in checks:
        a0 = atom_by(a_key, full_atoms)
        b0 = atom_by(b_key, full_atoms)
        a1 = atom_by(a_key, patched_atoms)
        b1 = atom_by(b_key, patched_atoms)
        if not a0 or not b0 or not a1 or not b1:
            rows.append(f"{label}\tNA\tNA\tNA\t{lo:.2f}-{hi:.2f}\tmissing_atom")
            continue
        d0 = distance(a0.xyz, b0.xyz)
        d1 = distance(a1.xyz, b1.xyz)
        flag = "" if lo <= d1 <= hi else "out_of_range"
        rows.append(f"{label}\t{d0:.4f}\t{d1:.4f}\t{(d1-d0):+.4f}\t{lo:.2f}-{hi:.2f}\t{flag}")
    return rows


def fe_donor_rows(full_atoms: dict[Key, PdbAtom], patched_atoms: dict[Key, PdbAtom], cutoff: float) -> list[str]:
    fe_keys = [key for key in full_atoms if key[0].upper() in {"FE", "FE2", "FE3"} or key[3].upper() == "FE"]
    rows = ["atom\tinitial_Fe_distance_A\tpreview_Fe_distance_A\tdelta_A"]
    if not fe_keys:
        return rows + ["no_fe\tNA\tNA\tNA"]
    fe_key = fe_keys[0]
    fe0 = full_atoms[fe_key]
    fe1 = patched_atoms.get(fe_key, fe0)
    candidate_keys: set[Key] = set()
    for key, atom in full_atoms.items():
        if key == fe_key or atom.element not in {"N", "O", "S"}:
            continue
        d0 = distance(fe0.xyz, atom.xyz)
        d1 = distance(fe1.xyz, patched_atoms.get(key, atom).xyz)
        if min(d0, d1) <= cutoff:
            candidate_keys.add(key)
    for key in sorted(candidate_keys):
        atom0 = full_atoms[key]
        atom1 = patched_atoms.get(key, atom0)
        d0 = distance(fe0.xyz, atom0.xyz)
        d1 = distance(fe1.xyz, atom1.xyz)
        rows.append(f"{atom0.label}\t{d0:.4f}\t{d1:.4f}\t{(d1-d0):+.4f}")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full-pdb", type=Path, required=True)
    parser.add_argument("--initial-small-pdb", type=Path, required=True)
    parser.add_argument("--optimized-small-pdb", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--out-name", default="mcpb_original_opt_nserch008_preview.pdb")
    parser.add_argument("--fe-cutoff", type=float, default=4.0)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    full_lines, full_atoms = parse_pdb(args.full_pdb)
    _, initial_atoms = parse_pdb(args.initial_small_pdb)
    _, opt_atoms = parse_pdb(args.optimized_small_pdb)

    preview_pdb = args.outdir / args.out_name
    mapping = write_patched_full_pdb(full_lines, full_atoms, opt_atoms, preview_pdb)
    _, patched_atoms = parse_pdb(preview_pdb)

    mapping_rows = ["atom\tinitial_x\tinitial_y\tinitial_z\tpreview_x\tpreview_y\tpreview_z\tdisplacement_A"]
    for key, xyz0, xyz1, delta in mapping:
        mapping_rows.append(
            f"{':'.join(key)}\t{xyz0[0]:.6f}\t{xyz0[1]:.6f}\t{xyz0[2]:.6f}\t"
            f"{xyz1[0]:.6f}\t{xyz1[1]:.6f}\t{xyz1[2]:.6f}\t{delta:.6f}"
        )
    (args.outdir / "full_pdb_patch_mapping.tsv").write_text("\n".join(mapping_rows) + "\n", encoding="utf-8")
    (args.outdir / "protein_key_bond_audit.tsv").write_text(
        "\n".join(bond_rows(full_atoms, patched_atoms)) + "\n", encoding="utf-8"
    )
    (args.outdir / "fe_donor_preview_compare.tsv").write_text(
        "\n".join(fe_donor_rows(full_atoms, patched_atoms, args.fe_cutoff)) + "\n", encoding="utf-8"
    )
    print(f"preview_pdb={preview_pdb}")
    print(f"patched_atoms={len(mapping)}")
    print(f"outdir={args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
