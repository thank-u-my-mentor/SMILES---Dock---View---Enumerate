#!/usr/bin/env python3
"""Create small Fe-binding ligand building blocks with OpenBabel."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


LIGANDS = [
    {
        "name": "acetate_raw_acetic_acid",
        "resname": "ACE",
        "smiles": "CC(O)=O",
        "note": "neutral acetic acid from the user SMILES; not ideal for Fe coordination until deprotonated",
    },
    {
        "name": "acetate_fe_bound_carboxylate",
        "resname": "ACO",
        "smiles": "CC(=O)[O-]",
        "note": "deprotonated carboxylate; use one O atom as Fe donor, or both O atoms for bidentate acetate-like binding",
    },
    {
        "name": "nhbutanamide_raw_radical_hint",
        "resname": "NHB",
        "smiles": "[NH]C(CCCC)=O",
        "note": "user radical-like amide SMILES; OpenBabel geometry only, radical electronic state still needs QM/topology",
    },
    {
        "name": "nhbutanamide_fe_bound_n_deprot",
        "resname": "NBR",
        "smiles": "[N-]C(CCCC)=O",
        "note": "N-deprotonated Fe-donor proxy; use for PyMOL placement, not final radical force field",
    },
]


def run(cmd: list[str]) -> None:
    print("$ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def convert_smiles(name: str, smiles: str, outdir: Path) -> None:
    smi = outdir / f"{name}.smi"
    smi.write_text(smiles + "\n", encoding="utf-8")
    sdf = outdir / f"{name}.sdf"
    mol2 = outdir / f"{name}.mol2"
    pdb = outdir / f"{name}.pdb"
    run(["obabel", "-ismi", str(smi), "-osdf", "-O", str(sdf), "--gen3d", "--best", "-h"])
    run(["obabel", str(sdf), "-O", str(mol2)])
    run(["obabel", str(sdf), "-O", str(pdb)])


def display_path(path: Path) -> str:
    text = path.as_posix()
    if text.startswith("/mnt/") and len(text) > 6 and text[6] == "/":
        drive = text[5].upper()
        return f"{drive}:{text[6:]}"
    return text


def write_pymol_helper(outdir: Path, complex_structure: Path | None) -> None:
    load_complex = f"load {display_path(complex_structure)}, complex\n" if complex_structure else ""
    loads = "\n".join(f"load {item['name']}.pdb, {item['name']}" for item in LIGANDS)
    text = f"""# Load small Fe-binding ligand building blocks for manual placement.
cd {display_path(outdir)}
{load_complex}{loads}

hide everything
show cartoon, polymer.protein
color gray70, polymer.protein
show spheres, elem Fe
set sphere_scale, 0.45, elem Fe
show sticks, organic
util.cbag organic

# Suggested manual workflow in PyMOL:
# 1. Use Mouse > 3 Button Editing, or enable editing mode.
# 2. Translate/rotate acetate_fe_bound_carboxylate and nhbutanamide_fe_bound_n_deprot near Fe.
# 3. For visual-only coordination lines after placement, use:
#      dist fe_acetate_O, elem Fe, acetate_fe_bound_carboxylate and elem O
#      dist fe_nitrene_N, elem Fe, nhbutanamide_fe_bound_n_deprot and elem N
# 4. Save the placed complex as PDB for MCPB/QM preparation:
#      save placed_fe_ligand_complex.pdb
zoom elem Fe or organic, 8
"""
    (outdir / "load_fe_ligand_building_blocks.pml").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("fe_ligand_building_blocks"))
    parser.add_argument("--complex-structure", type=Path, help="Optional structure to load with the ligand blocks")
    args = parser.parse_args()
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for item in LIGANDS:
        convert_smiles(item["name"], item["smiles"], outdir)
        rows.append(item)
    with (outdir / "ligand_building_blocks_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "resname", "smiles", "note"])
        writer.writeheader()
        writer.writerows(rows)
    write_pymol_helper(outdir, args.complex_structure)
    print(f"outdir={outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
