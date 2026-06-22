from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Atom:
    serial: int
    name: str
    resn: str
    resi: str
    element: str
    x: float
    y: float
    z: float


def v_add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def v_sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def v_mul(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def v_dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def v_cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def v_norm(a):
    return math.sqrt(v_dot(a, a))


def v_unit(a):
    norm = v_norm(a)
    if norm < 1.0e-8:
        raise ValueError("zero vector")
    return v_mul(a, 1.0 / norm)


def parse_pdb_atoms(path: Path) -> list[Atom]:
    atoms: list[Atom] = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        element = line[76:78].strip() or line[12:16].strip()[0]
        atoms.append(
            Atom(
                serial=int(line[6:11]),
                name=line[12:16].strip(),
                resn=line[17:20].strip(),
                resi=line[22:26].strip() or "1",
                element=element,
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
            )
        )
    return atoms


def format_pdb_atom(i: int, atom: Atom) -> str:
    return (
        f"HETATM{i:5d} {atom.name:<4s} {atom.resn:>3s} A{int(atom.resi):4d}    "
        f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
        f"  1.00  0.00          {atom.element:>2s}"
    )


def atom_with(atom: Atom, serial: int, name: str, element: str, coord=None) -> Atom:
    x, y, z = coord if coord is not None else (atom.x, atom.y, atom.z)
    return Atom(serial, name, "UNL", "1", element, x, y, z)


def build_phenyl_unl(unl_atoms: list[Atom]) -> tuple[list[Atom], list[tuple[int, int, int]]]:
    heavy = [atom for atom in unl_atoms if atom.element.upper() != "H"]
    hydrogens = [atom for atom in unl_atoms if atom.element.upper() == "H"]
    if len(heavy) < 7:
        raise ValueError(f"expected at least 7 heavy atoms in UNL, found {len(heavy)}")

    # Original exported UNL heavy order:
    # N, carbonyl C, alpha C, beta C, gamma C, terminal methyl C, carbonyl O.
    n_atom, c1, c2, c3, c4, ipso_old, o_atom = heavy[:7]
    p_c4 = (c4.x, c4.y, c4.z)
    p_ipso = (ipso_old.x, ipso_old.y, ipso_old.z)
    p_c3 = (c3.x, c3.y, c3.z)

    chain_axis = v_unit(v_sub(p_ipso, p_c4))
    prev_axis = v_unit(v_sub(p_c4, p_c3))
    in_plane = v_sub(prev_axis, v_mul(chain_axis, v_dot(prev_axis, chain_axis)))
    if v_norm(in_plane) < 1.0e-4:
        trial = (0.0, 0.0, 1.0)
        if abs(v_dot(trial, chain_axis)) > 0.9:
            trial = (0.0, 1.0, 0.0)
        in_plane = v_cross(chain_axis, trial)
    in_plane = v_unit(in_plane)

    ring_radius = 1.397
    ch_bond = 1.083
    center = v_add(p_ipso, v_mul(chain_axis, ring_radius))
    e1 = v_mul(chain_axis, -1.0)  # center -> ipso
    e2 = in_plane

    ring_positions = []
    for angle_deg in (0, 60, 120, 180, 240, 300):
        theta = math.radians(angle_deg)
        direction = v_add(v_mul(e1, math.cos(theta)), v_mul(e2, math.sin(theta)))
        ring_positions.append(v_add(center, v_mul(direction, ring_radius)))

    # Keep all non-terminal-methyl atoms from the original UNL. The three H atoms
    # bonded to the old terminal methyl are the last three hydrogens in this PDB.
    kept_atoms = [
        atom_with(n_atom, 1, "N1", "N"),
        atom_with(c1, 2, "C1", "C"),
        atom_with(c2, 3, "C2", "C"),
        atom_with(c3, 4, "C3", "C"),
        atom_with(c4, 5, "C4", "C"),
        atom_with(ipso_old, 6, "C5", "C", ring_positions[0]),
        atom_with(o_atom, 7, "O1", "O"),
    ]

    for name, old in zip(("HN1", "H2A", "H2B", "H3A", "H3B", "H4A", "H4B"), hydrogens[:7]):
        kept_atoms.append(atom_with(old, len(kept_atoms) + 1, name, "H"))

    for i, coord in enumerate(ring_positions[1:], start=6):
        kept_atoms.append(Atom(len(kept_atoms) + 1, f"C{i}", "UNL", "1", "C", *coord))

    for i, coord in enumerate(ring_positions[1:], start=6):
        direction = v_unit(v_sub(coord, center))
        h_coord = v_add(coord, v_mul(direction, ch_bond))
        kept_atoms.append(Atom(len(kept_atoms) + 1, f"H{i}", "UNL", "1", "H", *h_coord))

    # MOL-style bond types: 1 single, 2 double, 4 aromatic.
    bonds = [
        (1, 2, 1),
        (2, 3, 1),
        (2, 7, 2),
        (3, 4, 1),
        (4, 5, 1),
        (5, 6, 1),
        (1, 8, 1),
        (3, 9, 1),
        (3, 10, 1),
        (4, 11, 1),
        (4, 12, 1),
        (5, 13, 1),
        (5, 14, 1),
        (6, 15, 4),
        (15, 16, 4),
        (16, 17, 4),
        (17, 18, 4),
        (18, 19, 4),
        (19, 6, 4),
        (15, 20, 1),
        (16, 21, 1),
        (17, 22, 1),
        (18, 23, 1),
        (19, 24, 1),
    ]
    return kept_atoms, bonds


def write_pdb(path: Path, atoms: list[Atom], bonds: list[tuple[int, int, int]]) -> None:
    lines = [format_pdb_atom(i, atom) for i, atom in enumerate(atoms, start=1)]
    by_atom: dict[int, list[int]] = {}
    for a, b, _order in bonds:
        by_atom.setdefault(a, []).append(b)
        by_atom.setdefault(b, []).append(a)
    for atom_index in sorted(by_atom):
        bonded = "".join(f"{x:5d}" for x in sorted(by_atom[atom_index]))
        lines.append(f"CONECT{atom_index:5d}{bonded}")
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def write_mol(path: Path, atoms: list[Atom], bonds: list[tuple[int, int, int]]) -> None:
    lines = [
        "UNL_terminal_phenyl",
        "  PLE",
        "",
        f"{len(atoms):3d}{len(bonds):3d}  0  0  0  0            999 V2000",
    ]
    for atom in atoms:
        lines.append(
            f"{atom.x:10.4f}{atom.y:10.4f}{atom.z:10.4f} "
            f"{atom.element:<3s} 0  0  0  0  0  0  0  0  0  0  0  0"
        )
    for a, b, order in bonds:
        lines.append(f"{a:3d}{b:3d}{order:3d}  0  0  0  0")
    lines.append("M  END")
    path.write_text("\n".join(lines) + "\n")


def write_pml(path: Path, fixed_pdb: Path, fixed_pse: Path) -> None:
    fixed_pdb_win = str(fixed_pdb).replace("/mnt/e/", "E:/")
    fixed_pse_win = str(fixed_pse).replace("/mnt/e/", "E:/")
    path.write_text(
        "\n".join(
            [
                "# Run after loading E:/TJ/616.pse if you want to patch manually.",
                "remove resn UNL",
                f"load {fixed_pdb_win}, UNL_phenyl",
                "hide everything, UNL_phenyl",
                "show sticks, UNL_phenyl",
                "color palegreen, UNL_phenyl and elem C",
                "show spheres, UNL_phenyl and elem N+O",
                "set sphere_scale, 0.22, UNL_phenyl and elem N+O",
                "show_metal",
                f"save {fixed_pse_win}",
                "",
            ]
        )
    )


def update_pse(original_pse: Path, fixed_pdb: Path, fixed_pse: Path) -> None:
    import pymol

    pymol.finish_launching(["pymol", "-cq"])
    from pymol import cmd

    cmd.load(str(original_pse))
    cmd.remove("resn UNL")
    cmd.load(str(fixed_pdb), "UNL_phenyl")
    cmd.hide("everything", "UNL_phenyl")
    cmd.show("sticks", "UNL_phenyl")
    cmd.color("palegreen", "UNL_phenyl and elem C")
    cmd.show("spheres", "UNL_phenyl and elem N+O")
    cmd.set("sphere_scale", 0.22, "UNL_phenyl and elem N+O")
    if "show_metal" in cmd.keyword:
        cmd.do("show_metal")
    cmd.save(str(fixed_pse))
    cmd.quit()


def main() -> int:
    if len(sys.argv) != 4:
        print("usage: fix_unl_terminal_phenyl.py original.pse exported_UNL.pdb outdir")
        return 2

    original_pse = Path(sys.argv[1])
    unl_pdb = Path(sys.argv[2])
    outdir = Path(sys.argv[3])
    outdir.mkdir(parents=True, exist_ok=True)

    atoms = parse_pdb_atoms(unl_pdb)
    fixed_atoms, bonds = build_phenyl_unl(atoms)

    fixed_pdb = outdir / "UNL_phenyl.pdb"
    fixed_mol = outdir / "UNL_phenyl.mol"
    fixed_sdf = outdir / "UNL_phenyl.sdf"
    fixed_pse = outdir / "616_UNL_phenyl_fixed.pse"
    patch_pml = outdir / "patch_616_UNL_to_phenyl.pml"

    write_pdb(fixed_pdb, fixed_atoms, bonds)
    write_mol(fixed_mol, fixed_atoms, bonds)
    fixed_sdf.write_text(fixed_mol.read_text() + "$$$$\n")
    write_pml(patch_pml, fixed_pdb, fixed_pse)
    update_pse(original_pse, fixed_pdb, fixed_pse)

    print(f"fixed_pdb={fixed_pdb}")
    print(f"fixed_sdf={fixed_sdf}")
    print(f"fixed_mol={fixed_mol}")
    print(f"patch_pml={patch_pml}")
    print(f"fixed_pse={fixed_pse}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
