#!/usr/bin/env python3
"""Prepare GYM 4X8B chain-A + pose1 files for H++/MCPB/MD handoff."""

from __future__ import annotations

import argparse
from pathlib import Path

from rdkit import Chem


def pdb_atom_line_with_serial(line: str, serial: int) -> str:
    return f"{line[:6]}{serial:5d}{line[11:]}"


def keep_chain_a_receptor(line: str) -> bool:
    if line.startswith("ATOM  "):
        return line[21:22] == "A"
    if not line.startswith("HETATM"):
        return False
    resname = line[17:20].strip().upper()
    chain = line[21:22]
    if chain != "A":
        return False
    if resname in {"GOL", "HOH", "WAT"}:
        return False
    return resname in {"FE", "FE2", "FE3", "MG", "CA", "CL"}


def write_chain_a_receptor(source: Path, out: Path) -> None:
    serial = 1
    out_lines: list[str] = []
    last_was_atom = False
    for raw in source.read_text(errors="replace").splitlines():
        if keep_chain_a_receptor(raw):
            if raw.startswith("HETATM") and last_was_atom and (not out_lines or not out_lines[-1].startswith("TER")):
                out_lines.append(f"TER   {serial:5d}")
                serial += 1
            out_lines.append(pdb_atom_line_with_serial(raw.rstrip(), serial))
            serial += 1
            last_was_atom = raw.startswith("ATOM  ")
            continue
        if raw.startswith("TER") and last_was_atom:
            out_lines.append(f"TER   {serial:5d}")
            serial += 1
            last_was_atom = False
    if out_lines and not out_lines[-1].startswith("TER"):
        out_lines.append(f"TER   {serial:5d}")
    out_lines.append("END")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(out_lines) + "\n", encoding="utf-8")


def pdbqt_mode_atoms(path: Path, mode: int) -> list[tuple[str, float, float, float]]:
    atoms: list[tuple[str, float, float, float]] = []
    current = 1
    in_mode = mode == 1
    saw_model = False
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("MODEL"):
            saw_model = True
            try:
                current = int(line.split()[1])
            except Exception:
                current += 1
            in_mode = current == mode
            continue
        if line.startswith("ENDMDL"):
            if in_mode:
                break
            in_mode = False
            continue
        if saw_model and not in_mode:
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        if name.upper().startswith("H"):
            continue
        atom_type = (line[70:].split() or [""])[-1]
        element = atom_type
        if element == "A":
            element = "C"
        if len(element) > 1:
            element = element[:2].capitalize()
        else:
            element = element.upper()
        atoms.append((element, float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not atoms:
        raise RuntimeError(f"no atoms found for mode {mode} in {path}")
    return atoms


def write_ligand_pose_pdb(pdbqt: Path, sdf: Path, out: Path, mode: int, resname: str = "LIG") -> None:
    atoms = pdbqt_mode_atoms(pdbqt, mode)
    mol = next((m for m in Chem.SDMolSupplier(str(sdf), removeHs=False) if m is not None), None)
    if mol is None:
        raise RuntimeError(f"could not read ligand SDF: {sdf}")
    heavy = [a for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    if len(heavy) != len(atoms):
        raise RuntimeError(f"heavy atom count mismatch: pdbqt={len(atoms)} sdf={len(heavy)}")
    lines = []
    for idx, (atom, (element, x, y, z)) in enumerate(zip(heavy, atoms), start=1):
        symbol = atom.GetSymbol()
        name = f"{symbol}{idx}"
        lines.append(
            f"HETATM{idx:5d} {name:<4s} {resname:>3s} Z   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {symbol:>2s}"
        )
    lines.append(f"TER   {len(lines) + 1:5d}")
    lines.append("END")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--pose-pdbqt", type=Path, required=True)
    parser.add_argument("--pose-sdf", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--mode", type=int, default=1)
    args = parser.parse_args()

    receptor = args.outdir / "receptor_chainA_noGOL_noHOH_keepFe.pdb"
    ligand = args.outdir / f"ligand_S000001_mode{args.mode}.pdb"
    complex_pdb = args.outdir / f"complex_chainA_Fe_LIG_mode{args.mode}_for_mcpb.pdb"
    write_chain_a_receptor(args.pdb, receptor)
    write_ligand_pose_pdb(args.pose_pdbqt, args.pose_sdf, ligand, args.mode)
    receptor_text = receptor.read_text()
    ligand_lines = [line for line in ligand.read_text().splitlines() if line != "END"]
    complex_pdb.write_text(receptor_text.replace("END\n", "") + "\n".join(ligand_lines) + "\nEND\n", encoding="utf-8")
    print(f"receptor={receptor}")
    print(f"ligand={ligand}")
    print(f"complex={complex_pdb}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
