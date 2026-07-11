#!/usr/bin/env python3
"""Export a GAMESS-optimized MCPB small model and audit ligand identity."""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from export_gamess_optimized_model import (  # noqa: E402
    ModelAtom,
    apply_transform,
    best_rigid_transform,
    parse_coords_for_nserch,
)


Vec = tuple[float, float, float]


@dataclass(frozen=True)
class PdbAtom:
    serial: int
    name: str
    resname: str
    chain: str
    resseq: int
    element: str
    xyz: Vec
    record: str

    @property
    def label(self) -> str:
        return f"{self.resname}:{self.chain}:{self.resseq}:{self.name}"


def distance(a: Vec, b: Vec) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def parse_pdb(path: Path) -> list[PdbAtom]:
    atoms: list[PdbAtom] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        padded = (line + " " * 80)[:80]
        element = padded[76:78].strip()
        if not element:
            element = re.sub(r"[^A-Za-z]", "", padded[12:16]).strip()[:2]
        element = "Fe" if element.upper() == "FE" else element[:1].upper()
        atoms.append(
            PdbAtom(
                serial=int(padded[6:11]),
                name=padded[12:16].strip(),
                resname=padded[17:20].strip(),
                chain=padded[21].strip() or "A",
                resseq=int(padded[22:26]),
                element=element,
                xyz=(
                    float(padded[30:38]),
                    float(padded[38:46]),
                    float(padded[46:54]),
                ),
                record=padded[:6].strip() or "HETATM",
            )
        )
    if not atoms:
        raise SystemExit(f"No atoms parsed from {path}")
    return atoms


def to_model_atoms(pdb_atoms: list[PdbAtom]) -> list[ModelAtom]:
    return [
        ModelAtom(atom.element.upper(), atom.xyz, atom.label)
        for atom in pdb_atoms
    ]


def write_pdb(path: Path, template: list[PdbAtom], coords: list[Vec], title: str) -> None:
    lines = [f"REMARK {title}"]
    for idx, (atom, xyz) in enumerate(zip(template, coords), start=1):
        record = "ATOM" if atom.record == "ATOM" else "HETATM"
        elem = "Fe" if atom.element.upper() == "FE" else atom.element[:1].upper()
        lines.append(
            f"{record:<6}{idx:5d} {atom.name:>4s} {atom.resname:>3s} {atom.chain:1s}"
            f"{atom.resseq:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            f"  1.00  0.00          {elem:>2s}"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_mol2(path: Path) -> list[tuple[str, str, str]]:
    atom_names: dict[int, str] = {}
    bonds: list[tuple[str, str, str]] = []
    section = ""
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("@<TRIPOS>"):
            section = line.strip()
            continue
        parts = line.split()
        if not parts:
            continue
        if section == "@<TRIPOS>ATOM" and len(parts) >= 6:
            atom_names[int(parts[0])] = parts[1]
        elif section == "@<TRIPOS>BOND" and len(parts) >= 4:
            a = atom_names.get(int(parts[1]))
            b = atom_names.get(int(parts[2]))
            if a and b:
                bonds.append((a, b, parts[3]))
    if not bonds:
        raise SystemExit(f"No bonds parsed from {path}")
    return bonds


def coords_by_atom_name(atoms: list[PdbAtom], coords: list[Vec], resname: str) -> dict[str, Vec]:
    found: dict[str, Vec] = {}
    for atom, xyz in zip(atoms, coords):
        if atom.resname == resname:
            found[atom.name] = xyz
    if not found:
        raise SystemExit(f"No residue {resname} found in model")
    return found


def audit_mol2_bonds(
    *,
    mol2: Path,
    resname: str,
    atoms: list[PdbAtom],
    initial_coords: list[Vec],
    optimized_coords: list[Vec],
    long_cutoff: float,
    short_cutoff: float,
    large_change_cutoff: float,
) -> tuple[list[str], int]:
    bonds = parse_mol2(mol2)
    initial = coords_by_atom_name(atoms, initial_coords, resname)
    optimized = coords_by_atom_name(atoms, optimized_coords, resname)
    rows = ["resname\tbond\torder\tinitial_A\toptimized_A\tdelta_A\tflag"]
    problems = 0
    for atom_a, atom_b, order in bonds:
        flags: list[str] = []
        if atom_a not in optimized or atom_b not in optimized:
            rows.append(f"{resname}\t{atom_a}-{atom_b}\t{order}\tNA\tNA\tNA\tMISSING_OPT_ATOM")
            problems += 1
            continue
        if atom_a not in initial or atom_b not in initial:
            rows.append(f"{resname}\t{atom_a}-{atom_b}\t{order}\tNA\tNA\tNA\tMISSING_INITIAL_ATOM")
            problems += 1
            continue
        d0 = distance(initial[atom_a], initial[atom_b])
        d1 = distance(optimized[atom_a], optimized[atom_b])
        delta = d1 - d0
        if d1 > long_cutoff:
            flags.append("BROKEN_LONG_BOND")
        if d1 < short_cutoff:
            flags.append("TOO_SHORT")
        if abs(delta) > large_change_cutoff:
            flags.append("LARGE_CHANGE")
        if flags:
            problems += 1
        rows.append(
            f"{resname}\t{atom_a}-{atom_b}\t{order}\t{d0:.4f}\t{d1:.4f}\t{delta:+.4f}\t{','.join(flags)}"
        )
    return rows, problems


def write_fe_distances(path: Path, atoms: list[PdbAtom], initial: list[Vec], optimized: list[Vec], cutoff: float) -> None:
    fe_indices = [idx for idx, atom in enumerate(atoms) if atom.element.upper() == "FE" or atom.name.upper() == "FE"]
    if not fe_indices:
        path.write_text("No Fe atom found.\n", encoding="utf-8")
        return
    fe_idx = fe_indices[0]
    rows = ["atom\tinitial_Fe_distance_A\toptimized_Fe_distance_A\tdelta_A"]
    for idx, atom in enumerate(atoms):
        if idx == fe_idx:
            continue
        if atom.element.upper() not in {"N", "O", "S"}:
            continue
        d0 = distance(initial[fe_idx], initial[idx])
        d1 = distance(optimized[fe_idx], optimized[idx])
        if min(d0, d1) <= cutoff:
            rows.append(f"{atom.label}\t{d0:.4f}\t{d1:.4f}\t{(d1 - d0):+.4f}")
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def write_displacements(path: Path, atoms: list[PdbAtom], initial: list[Vec], optimized: list[Vec]) -> None:
    rows = ["atom\tinitial_x\tinitial_y\tinitial_z\toptimized_x\toptimized_y\toptimized_z\tdisplacement_A"]
    for atom, xyz0, xyz1 in zip(atoms, initial, optimized):
        rows.append(
            f"{atom.label}\t{xyz0[0]:.6f}\t{xyz0[1]:.6f}\t{xyz0[2]:.6f}\t"
            f"{xyz1[0]:.6f}\t{xyz1[1]:.6f}\t{xyz1[2]:.6f}\t{distance(xyz0, xyz1):.6f}"
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template-pdb", type=Path, required=True)
    parser.add_argument("--opt-log", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--nserch", type=int)
    parser.add_argument("--ligand", nargs=3, action="append", metavar=("RESNAME", "MOL2", "LABEL"), default=[])
    parser.add_argument("--long-bond-cutoff", type=float, default=2.2)
    parser.add_argument("--short-bond-cutoff", type=float, default=0.85)
    parser.add_argument("--large-change-cutoff", type=float, default=0.35)
    parser.add_argument("--fe-cutoff", type=float, default=3.2)
    args = parser.parse_args()

    atoms = parse_pdb(args.template_pdb)
    initial_model = to_model_atoms(atoms)
    nserch, final_model, metadata = parse_coords_for_nserch(args.opt_log, len(atoms), args.nserch)
    rot, trans, rmsd = best_rigid_transform([atom.xyz for atom in final_model], [atom.xyz for atom in initial_model])
    initial_coords = [atom.xyz for atom in initial_model]
    optimized_coords = [apply_transform(atom.xyz, rot, trans) for atom in final_model]

    args.outdir.mkdir(parents=True, exist_ok=True)
    initial_pdb = args.outdir / "initial_small_model.pdb"
    optimized_pdb = args.outdir / "optimized_small_model_aligned.pdb"
    write_pdb(initial_pdb, atoms, initial_coords, "Initial MCPB small model coordinates")
    write_pdb(
        optimized_pdb,
        atoms,
        optimized_coords,
        f"GAMESS optimized MCPB small model aligned to initial frame, NSERCH={nserch}",
    )
    write_displacements(args.outdir / "atom_displacement.tsv", atoms, initial_coords, optimized_coords)
    write_fe_distances(args.outdir / "fe_donor_distance_compare.tsv", atoms, initial_coords, optimized_coords, args.fe_cutoff)

    total_problems = 0
    audit_files: list[str] = []
    for resname, mol2_text, label in args.ligand:
        rows, problems = audit_mol2_bonds(
            mol2=Path(mol2_text),
            resname=resname,
            atoms=atoms,
            initial_coords=initial_coords,
            optimized_coords=optimized_coords,
            long_cutoff=args.long_bond_cutoff,
            short_cutoff=args.short_bond_cutoff,
            large_change_cutoff=args.large_change_cutoff,
        )
        total_problems += problems
        audit_path = args.outdir / f"{resname}_mol2_bond_identity_audit.tsv"
        audit_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
        audit_files.append(f"- `{audit_path.name}`: {label}, problem_count={problems}")

    summary = [
        "# MCPB OPT Chemical Identity Export",
        "",
        f"- Optimized geometry: `NSERCH={nserch}`",
        f"- Energy: `{metadata.get('energy_hartree', 'NA')}` Hartree",
        f"- GRAD. MAX / RMS: `{metadata.get('grad_max', 'NA')}` / `{metadata.get('grad_rms', 'NA')}`",
        f"- S-squared: `{metadata.get('s2', 'NA')}`",
        f"- Alignment RMSD back to initial frame: `{rmsd:.6f} A`",
        f"- Ligand identity problem count: `{total_problems}`",
        "",
        "## Files",
        "",
        f"- `{initial_pdb.name}`: optimization 前的小模型 PDB。",
        f"- `{optimized_pdb.name}`: optimization 后的小模型 PDB，已刚体对齐回原始坐标系。",
        "- `atom_displacement.tsv`: 每个原子的位移。",
        "- `fe_donor_distance_compare.tsv`: Fe 到近邻 N/O/S donor 的距离变化。",
        *audit_files,
        "",
        "## Interpretation",
        "",
        "如果 ligand identity problem count 为 0，说明 mol2 中定义的 UNK/ACT 共价键在优化后没有出现过长、过短或大幅改变；",
        "也就是说这次 OPT 没有像先前 bidentate 模型那样把 ligand 优化成断键/反应产物状态。",
    ]
    (args.outdir / "README_chemical_identity.md").write_text(
        "\n".join(summary) + "\n",
        encoding="utf-8",
    )

    print(f"outdir={args.outdir}")
    print(f"nserch={nserch}")
    print(f"energy_hartree={metadata.get('energy_hartree', 'NA')}")
    print(f"grad_max={metadata.get('grad_max', 'NA')}")
    print(f"grad_rms={metadata.get('grad_rms', 'NA')}")
    print(f"s2={metadata.get('s2', 'NA')}")
    print(f"alignment_rmsd_A={rmsd:.6f}")
    print(f"optimized_pdb={optimized_pdb}")
    print(f"ligand_identity_problem_count={total_problems}")
    return 2 if total_problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
