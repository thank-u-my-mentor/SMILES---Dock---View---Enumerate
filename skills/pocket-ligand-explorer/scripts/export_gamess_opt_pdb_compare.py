#!/usr/bin/env python3
"""Export initial/final GAMESS small-model PDBs and a geometry comparison table."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys
from dataclasses import dataclass

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from export_gamess_optimized_model import (  # noqa: E402
    ModelAtom,
    apply_transform,
    best_rigid_transform,
    parse_coords_for_nserch,
    parse_xyz,
)


@dataclass
class TemplateAtom:
    name: str
    resn: str
    chain: str
    resi: str
    element: str


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def element_from_name(name: str) -> str:
    letters = "".join(ch for ch in name if ch.isalpha()).upper()
    if letters.startswith("FE"):
        return "Fe"
    if letters:
        return letters[0].capitalize()
    return "C"


def parse_pdb_template(path: Path | None) -> list[TemplateAtom] | None:
    if path is None:
        return None
    atoms: list[TemplateAtom] = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        padded = (line + " " * 80)[:80]
        name = padded[12:16].strip()
        resn = padded[17:20].strip()
        chain = padded[21].strip() or "A"
        resi = padded[22:26].strip() or "1"
        elem = padded[76:78].strip()
        if not elem:
            elem = element_from_name(name)
        if elem.upper() == "FE":
            elem = "Fe"
        else:
            elem = elem.capitalize()
        atoms.append(TemplateAtom(name=name, resn=resn, chain=chain, resi=resi, element=elem))
    return atoms


def atom_metadata(atom, index: int, template: list[TemplateAtom] | None = None) -> tuple[str, str, str, str, str]:
    if template is not None and index - 1 < len(template):
        item = template[index - 1]
        return item.resn[:3], item.resi, item.name[:4], item.chain, item.element
    token = atom.note.split()[0] if atom.note.split() else ""
    parts = token.split(":")
    if len(parts) == 3:
        resn, resi, name = parts
        chain = "A"
    else:
        resn, resi, name, chain = "SM", "1", atom.element, "A"
    elem = "Fe" if atom.element.upper() == "FE" else atom.element.capitalize()
    return resn[:3], resi, name[:4], chain, elem or element_from_name(name)


def write_pdb(path: Path, atoms, coords, *, title: str, template: list[TemplateAtom] | None = None) -> None:
    lines = [f"REMARK {title}"]
    for serial, (atom, coord_atom) in enumerate(zip(atoms, coords), start=1):
        resn, resi, name, chain, elem = atom_metadata(atom, serial, template)
        try:
            resseq = int("".join(ch for ch in resi if ch.isdigit()) or "1")
        except ValueError:
            resseq = serial
        x, y, z = coord_atom.xyz if hasattr(coord_atom, "xyz") else coord_atom
        lines.append(
            f"HETATM{serial:5d} {name:>4s} {resn:>3s} {chain:1s}"
            f"{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"  1.00  0.00          {elem:>2s}"
        )
    lines.append("END")
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")


def write_tsv(
    path: Path,
    atoms,
    initial_coords,
    final_coords,
    *,
    template: list[TemplateAtom] | None = None,
) -> list[tuple[str, float]]:
    rows: list[str] = [
        "index\tresidue\tatom\telement\tinitial_x\tinitial_y\tinitial_z\t"
        "optimized_x\toptimized_y\toptimized_z\tdisplacement_A"
    ]
    displacements: list[tuple[str, float]] = []
    for index, (atom, init_atom, final_atom) in enumerate(
        zip(atoms, initial_coords, final_coords),
        start=1,
    ):
        resn, resi, name, chain, elem = atom_metadata(atom, index, template)
        label = f"{resn}:{chain}{resi}:{name}"
        d = distance(init_atom.xyz, final_atom.xyz)
        displacements.append((label, d))
        rows.append(
            f"{index}\t{resn}:{chain}{resi}\t{name}\t{elem}\t"
            f"{init_atom.xyz[0]:.6f}\t{init_atom.xyz[1]:.6f}\t{init_atom.xyz[2]:.6f}\t"
            f"{final_atom.xyz[0]:.6f}\t{final_atom.xyz[1]:.6f}\t{final_atom.xyz[2]:.6f}\t{d:.6f}"
        )
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(rows) + "\n")
    return displacements


def write_distance_compare(
    path: Path,
    atoms,
    initial_coords,
    final_coords,
    *,
    template: list[TemplateAtom] | None = None,
) -> None:
    labels = []
    fe_index = None
    for index, atom in enumerate(atoms):
        resn, resi, name, chain, _ = atom_metadata(atom, index + 1, template)
        label = f"{resn}:{chain}{resi}:{name}"
        labels.append(label)
        if atom.element.upper() == "FE" and fe_index is None:
            fe_index = index
    if fe_index is None:
        with path.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("No Fe atom found.\n")
        return
    rows = ["atom\tinitial_Fe_distance_A\toptimized_Fe_distance_A\tdelta_A"]
    fe0 = initial_coords[fe_index].xyz
    fe1 = final_coords[fe_index].xyz
    donor_tokens = ("NE2", "ND1", "OE1", "OE2", "O1", "O2", "N1", "N2", ":O")
    for index, label in enumerate(labels):
        if index == fe_index:
            continue
        if not any(token in label for token in donor_tokens):
            continue
        d0 = distance(fe0, initial_coords[index].xyz)
        d1 = distance(fe1, final_coords[index].xyz)
        if min(d0, d1) <= 3.2:
            rows.append(f"{label}\t{d0:.4f}\t{d1:.4f}\t{(d1 - d0):+.4f}")
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(rows) + "\n")


def write_markdown(path: Path, *, metadata: dict, displacements: list[tuple[str, float]]) -> None:
    top = sorted(displacements, key=lambda item: item[1], reverse=True)[:12]
    lines = [
        "# GAMESS OPT Before/After Export",
        "",
        f"- Evaluated NSERCH: `{int(metadata.get('nserch', -1))}`",
        f"- Energy: `{metadata.get('energy_hartree', 'NA')}` Hartree",
        f"- GRAD.MAX: `{metadata.get('grad_max', 'NA')}`",
        f"- RMS gradient: `{metadata.get('grad_rms', 'NA')}`",
        f"- S^2: `{metadata.get('s2', 'NA')}`",
        "",
        "## Files",
        "",
        "- `initial_small_model.pdb`: OPT 前的小模型坐标。",
        "- `optimized_small_model.pdb`: GAMESS 最终已评估收敛坐标。",
        "- `atom_displacement.tsv`: 每个原子的位移。",
        "- `fe_donor_distance_compare.tsv`: Fe 到近邻 donor 的距离变化。",
        "",
        "## Largest Atom Displacements",
        "",
        "| atom | displacement (A) |",
        "|---|---:|",
    ]
    for label, disp in top:
        lines.append(f"| `{label}` | {disp:.4f} |")
    lines += [
        "",
        "PyMOL quick view:",
        "",
        "```pml",
        "load initial_small_model.pdb, before_opt",
        "load optimized_small_model.pdb, after_opt",
        "hide everything",
        "show sticks, before_opt or after_opt",
        "color gray70, before_opt and elem C",
        "color palegreen, after_opt and elem C",
        "color blue, elem N",
        "color red, elem O",
        "show spheres, elem Fe",
        "set sphere_scale, 0.25, elem Fe",
        "orient",
        "```",
    ]
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--initial-xyz", required=True)
    parser.add_argument("--opt-log", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--nserch", type=int, default=None)
    parser.add_argument("--template-pdb", default=None)
    args = parser.parse_args()

    initial_xyz = Path(args.initial_xyz)
    opt_log = Path(args.opt_log)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    atoms = parse_xyz(initial_xyz)
    template = parse_pdb_template(Path(args.template_pdb)) if args.template_pdb else None
    if template is not None and len(template) != len(atoms):
        raise SystemExit(
            f"Template atom count {len(template)} does not match initial xyz atom count {len(atoms)}"
        )
    final_nserch, final_coords, metadata = parse_coords_for_nserch(
        opt_log,
        len(atoms),
        args.nserch,
    )
    initial_points = [atom.xyz for atom in atoms]
    final_points = [atom.xyz for atom in final_coords]
    rot, trans, align_rmsd = best_rigid_transform(final_points, initial_points)
    aligned_final_coords = [
        ModelAtom(atom.element, apply_transform(coord_atom.xyz, rot, trans), atom.note)
        for atom, coord_atom in zip(atoms, final_coords)
    ]
    metadata["nserch"] = final_nserch
    metadata["alignment_rmsd"] = align_rmsd
    write_pdb(
        outdir / "initial_small_model.pdb",
        atoms,
        atoms,
        title=f"Initial coordinates from {initial_xyz.name}",
        template=template,
    )
    write_pdb(
        outdir / "optimized_small_model.pdb",
        atoms,
        final_coords,
        title=f"Raw optimized evaluated NSERCH={final_nserch} from {opt_log.name}",
        template=template,
    )
    write_pdb(
        outdir / "optimized_small_model_aligned.pdb",
        atoms,
        aligned_final_coords,
        title=f"Aligned optimized evaluated NSERCH={final_nserch} from {opt_log.name}",
        template=template,
    )
    displacements = write_tsv(
        outdir / "atom_displacement.tsv",
        atoms,
        atoms,
        aligned_final_coords,
        template=template,
    )
    write_distance_compare(
        outdir / "fe_donor_distance_compare.tsv",
        atoms,
        atoms,
        aligned_final_coords,
        template=template,
    )
    write_markdown(outdir / "README_opt_compare.md", metadata=metadata, displacements=displacements)
    print(f"outdir={outdir}")
    print(f"nserch={final_nserch}")
    print(f"energy_hartree={metadata.get('energy_hartree', 'NA')}")
    print(f"grad_max={metadata.get('grad_max', 'NA')}")
    print(f"grad_rms={metadata.get('grad_rms', 'NA')}")
    print(f"s2={metadata.get('s2', 'NA')}")
    print(f"alignment_rmsd_A={align_rmsd}")
    print(f"initial_pdb={outdir / 'initial_small_model.pdb'}")
    print(f"optimized_pdb={outdir / 'optimized_small_model.pdb'}")
    print(f"optimized_aligned_pdb={outdir / 'optimized_small_model_aligned.pdb'}")
    print(f"comparison={outdir / 'atom_displacement.tsv'}")
    print(f"fe_distances={outdir / 'fe_donor_distance_compare.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
