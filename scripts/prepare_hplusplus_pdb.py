#!/usr/bin/env python
"""Prepare conservative PDB variants for H++ submission."""

from __future__ import annotations

import argparse
from pathlib import Path


STANDARD_AA = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}
WATER_NAMES = {"HOH", "WAT", "DOD", "SOL"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--chain", default="A")
    return parser.parse_args()


def atom_serialized(line: str, serial: int) -> str:
    return f"{line[:6]}{serial:5d}{line[11:80].rstrip():<69}"


def rename_fe_line(line: str) -> str:
    # H++ compatibility is not guaranteed for metals, but FE/FE is less exotic
    # than the CCD residue name FE2 if the user wants a quick manual test.
    chars = list(f"{line[:80]:<80}")
    chars[12:16] = f"{'FE':>4}"
    chars[17:20] = f"{'FE':>3}"
    chars[76:78] = f"{'FE':>2}"
    return "".join(chars).rstrip()


def collect_chain_atoms(lines: list[str], chain: str) -> list[str]:
    atoms: list[str] = []
    for line in lines:
        if not line.startswith("ATOM  "):
            continue
        if line[21:22] != chain:
            continue
        resname = line[17:20].strip()
        altloc = line[16:17]
        if resname not in STANDARD_AA:
            continue
        if altloc not in (" ", "A"):
            continue
        if altloc == "A":
            line = line[:16] + " " + line[17:]
        atoms.append(line[:80].rstrip())
    return atoms


def collect_hetatm(lines: list[str], chain: str, *, keep_resnames: set[str], keep_elements: set[str]) -> list[str]:
    atoms: list[str] = []
    for line in lines:
        if not line.startswith("HETATM"):
            continue
        if line[21:22] != chain:
            continue
        resname = line[17:20].strip()
        element = line[76:78].strip().upper()
        if resname in WATER_NAMES:
            continue
        if resname in keep_resnames or element in keep_elements:
            atoms.append(line[:80].rstrip())
    return atoms


def residue_gaps(atoms: list[str]) -> tuple[list[int], list[tuple[int, int]]]:
    residues: list[int] = []
    seen: set[int] = set()
    for line in atoms:
        if line[12:16].strip() != "CA":
            continue
        try:
            resid = int(line[22:26])
        except ValueError:
            continue
        if resid in seen:
            continue
        seen.add(resid)
        residues.append(resid)
    gaps = [(a, b) for a, b in zip(residues, residues[1:]) if b != a + 1]
    return residues, gaps


def write_pdb(path: Path, atoms: list[str], hetatm: list[str] | None = None, *, rename_fe: bool = False) -> None:
    serial = 1
    out_lines: list[str] = []
    for line in atoms:
        out_lines.append(atom_serialized(line, serial))
        serial += 1
    if atoms:
        last = atoms[-1]
        out_lines.append(f"TER   {serial:5d}      {last[17:20]} {last[21:22]}{last[22:26]}")
        serial += 1
    for line in hetatm or []:
        if rename_fe and line[76:78].strip().upper() == "FE":
            line = rename_fe_line(line)
        out_lines.append(atom_serialized(line, serial))
        serial += 1
    out_lines.append("END")
    path.write_text("\n".join(out_lines) + "\n", encoding="ascii")


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    lines = args.pdb.read_text(errors="replace").splitlines()
    atoms = collect_chain_atoms(lines, args.chain)
    residues, gaps = residue_gaps(atoms)

    stem = args.pdb.stem
    protein_only = args.outdir / f"{stem}_chain{args.chain}_Hpp_ready_protein_only.pdb"
    keep_fe = args.outdir / f"{stem}_chain{args.chain}_Hpp_try_keepFe_renamedFE.pdb"
    keep_fe_ntd = args.outdir / f"{stem}_chain{args.chain}_Hpp_try_keepFe_NTD_original_names.pdb"

    fe_atoms = collect_hetatm(lines, args.chain, keep_resnames=set(), keep_elements={"FE"})
    fe_ntd_atoms = collect_hetatm(lines, args.chain, keep_resnames={"NTD", "FE2"}, keep_elements={"FE"})

    write_pdb(protein_only, atoms)
    write_pdb(keep_fe, atoms, fe_atoms, rename_fe=True)
    write_pdb(keep_fe_ntd, atoms, fe_ntd_atoms, rename_fe=False)

    report = args.outdir / f"{stem}_chain{args.chain}_Hpp_prep_report.txt"
    report.write_text(
        "\n".join(
            [
                f"source={args.pdb}",
                f"chain={args.chain}",
                f"protein_residue_count={len(residues)}",
                f"first_residue={residues[0] if residues else 'NA'}",
                f"last_residue={residues[-1] if residues else 'NA'}",
                f"internal_ca_gaps={gaps if gaps else 'none'}",
                f"safe_hplusplus_file={protein_only}",
                f"try_keep_fe_file={keep_fe}",
                f"try_keep_fe_ntd_file={keep_fe_ntd}",
                "notes=The safest H++ input is protein-only. H++ may reject Fe or NTD because they are non-standard HETATM records; if so, run H++ on the protein-only file and merge protonation choices back onto the Fe/ligand receptor later.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(report)
    print(protein_only)
    print(keep_fe)
    print(keep_fe_ntd)


if __name__ == "__main__":
    main()
