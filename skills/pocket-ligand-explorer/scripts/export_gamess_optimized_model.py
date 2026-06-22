#!/usr/bin/env python3
"""Export an evaluated GAMESS optimized geometry and graft it onto a full PDB.

This is intended for Fe-radical MCPB preparation. GAMESS may rotate/translate
the input coordinates internally, so the evaluated optimized geometry is aligned
back to the original truncated model before it is written as a PDB or patched
into the full protein structure.
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path


Vec = tuple[float, float, float]


@dataclass
class ModelAtom:
    element: str
    xyz: Vec
    note: str = ""


@dataclass
class PdbAtom:
    line_index: int
    line: str
    serial: int
    name: str
    resn: str
    chain: str
    resi: str
    element: str
    xyz: Vec

    @property
    def key(self) -> tuple[str, str, str, str]:
        return (self.resn, self.chain, self.resi, self.name)


DONOR_KEYS = [
    ("HIS", "A", "187", "NE2"),
    ("HIS", "A", "270", "NE2"),
    ("GLU", "A", "349", "OE1"),
    ("ACT", "A", "501", "O2"),
    ("UNL", "A", "1", "N1"),
    ("HOH", "B", "875", "O"),
]
FE_KEY = ("FE2", "A", "431", "FE")


def v_add(a: Vec, b: Vec) -> Vec:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def v_sub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def v_scale(a: Vec, scalar: float) -> Vec:
    return (a[0] * scalar, a[1] * scalar, a[2] * scalar)


def v_dot(a: Vec, b: Vec) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def distance(a: Vec, b: Vec) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def centroid(points: list[Vec]) -> Vec:
    n = float(len(points))
    return (
        sum(point[0] for point in points) / n,
        sum(point[1] for point in points) / n,
        sum(point[2] for point in points) / n,
    )


def mat_vec_mul(mat: list[list[float]], vec: Vec) -> Vec:
    return (
        mat[0][0] * vec[0] + mat[0][1] * vec[1] + mat[0][2] * vec[2],
        mat[1][0] * vec[0] + mat[1][1] * vec[1] + mat[1][2] * vec[2],
        mat[2][0] * vec[0] + mat[2][1] * vec[1] + mat[2][2] * vec[2],
    )


def normalize4(vec: list[float]) -> list[float]:
    norm = math.sqrt(sum(value * value for value in vec))
    if norm < 1e-14:
        return [1.0, 0.0, 0.0, 0.0]
    return [value / norm for value in vec]


def largest_eigenvector_4x4(mat: list[list[float]]) -> list[float]:
    vec = [1.0, 0.0, 0.0, 0.0]
    for _ in range(200):
        nxt = [
            sum(mat[row][col] * vec[col] for col in range(4))
            for row in range(4)
        ]
        nxt = normalize4(nxt)
        if sum(abs(nxt[i] - vec[i]) for i in range(4)) < 1e-12:
            return nxt
        vec = nxt
    return vec


def quaternion_rotation(source: list[Vec], target: list[Vec]) -> tuple[list[list[float]], Vec, float]:
    """Return R,t,rmsd such that R*source + t best matches target."""
    if len(source) != len(target) or len(source) < 3:
        raise ValueError("Need equal source/target point lists with at least three atoms")
    cs = centroid(source)
    ct = centroid(target)
    cov = [[0.0, 0.0, 0.0] for _ in range(3)]
    for src, dst in zip(source, target):
        p = v_sub(src, cs)
        q = v_sub(dst, ct)
        for i in range(3):
            for j in range(3):
                cov[i][j] += p[i] * q[j]
    sxx, sxy, sxz = cov[0]
    syx, syy, syz = cov[1]
    szx, szy, szz = cov[2]
    trace = sxx + syy + szz
    nmat = [
        [trace, syz - szy, szx - sxz, sxy - syx],
        [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ]
    w, x, y, z = largest_eigenvector_4x4(nmat)
    rot = [
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
    ]
    trans = v_sub(ct, mat_vec_mul(rot, cs))
    fitted = [v_add(mat_vec_mul(rot, point), trans) for point in source]
    rmsd = math.sqrt(sum(distance(a, b) ** 2 for a, b in zip(fitted, target)) / len(source))
    return rot, trans, rmsd


def best_rigid_transform(source: list[Vec], target: list[Vec]) -> tuple[list[list[float]], Vec, float]:
    rot, trans, rmsd = quaternion_rotation(source, target)
    # Try reversed covariance convention too; pure-Python quaternion conventions
    # are easy to transpose accidentally, so choose the transform with lower RMSD.
    rot2, trans2, rmsd2 = quaternion_rotation(target, source)
    rot2_t = [[rot2[j][i] for j in range(3)] for i in range(3)]
    trans2_inv = mat_vec_mul(rot2_t, v_scale(trans2, -1.0))
    fitted = [v_add(mat_vec_mul(rot2_t, point), trans2_inv) for point in source]
    rmsd2_inv = math.sqrt(sum(distance(a, b) ** 2 for a, b in zip(fitted, target)) / len(source))
    if rmsd2_inv < rmsd:
        return rot2_t, trans2_inv, rmsd2_inv
    return rot, trans, rmsd


def apply_transform(point: Vec, rot: list[list[float]], trans: Vec) -> Vec:
    return v_add(mat_vec_mul(rot, point), trans)


def parse_xyz(path: Path) -> list[ModelAtom]:
    atoms: list[ModelAtom] = []
    lines = path.read_text(errors="replace").splitlines()
    for line in lines[2:]:
        if not line.strip():
            continue
        body, _, note = line.partition("#")
        parts = body.split()
        if len(parts) < 4:
            continue
        atoms.append(
            ModelAtom(
                parts[0].upper(),
                (float(parts[1]), float(parts[2]), float(parts[3])),
                note.strip(),
            )
        )
    return atoms


def parse_coordinate_block(lines: list[str], start: int, atom_count: int) -> list[ModelAtom]:
    for i in range(start, len(lines)):
        if "COORDINATES OF ALL ATOMS ARE (ANGS)" not in lines[i]:
            continue
        atoms: list[ModelAtom] = []
        for line in lines[i + 1 :]:
            parts = line.split()
            if len(parts) >= 5 and re.fullmatch(r"[A-Za-z][A-Za-z0-9]*", parts[0]):
                try:
                    atoms.append(
                        ModelAtom(
                            parts[0].upper(),
                            (float(parts[2]), float(parts[3]), float(parts[4])),
                        )
                    )
                except ValueError:
                    pass
                if len(atoms) == atom_count:
                    return atoms
            elif atoms:
                if len(atoms) == atom_count:
                    return atoms
                atoms = []
        break
    raise SystemExit(f"Could not parse a {atom_count}-atom GAMESS coordinate block")


def find_evaluated_nserch(lines: list[str], requested: int | None = None) -> int:
    searches = [
        int(match.group(1))
        for line in lines
        for match in [re.search(r"NSERCH:\s*(\d+)\s+E=", line)]
        if match
    ]
    if requested is not None:
        if requested not in searches:
            raise SystemExit(f"Requested NSERCH={requested} was not evaluated in the log")
        return requested
    if not searches:
        raise SystemExit("No evaluated NSERCH lines were found in the GAMESS log")
    return searches[-1]


def parse_coords_for_nserch(log: Path, atom_count: int, nserch: int | None = None) -> tuple[int, list[ModelAtom], dict[str, float | str]]:
    lines = log.read_text(errors="replace").splitlines()
    chosen = find_evaluated_nserch(lines, nserch)
    start = None
    pattern = re.compile(rf"BEGINNING GEOMETRY SEARCH POINT NSERCH=\s*{chosen}\b")
    for i, line in enumerate(lines):
        if pattern.search(line):
            start = i
    if start is None:
        raise SystemExit(f"Could not find coordinate block for evaluated NSERCH={chosen}")
    coords = parse_coordinate_block(lines, start, atom_count)
    metadata: dict[str, float | str] = {"nserch": float(chosen)}
    for line in lines:
        match = re.search(
            rf"NSERCH:\s*{chosen}\s+E=\s*(-?\d+\.\d+)\s+GRAD\. MAX=\s*([0-9.Ee+-]+)\s+R\.M\.S\.=\s*([0-9.Ee+-]+)",
            line,
        )
        if match:
            metadata["energy_hartree"] = float(match.group(1))
            metadata["grad_max"] = float(match.group(2))
            metadata["grad_rms"] = float(match.group(3))
    s2_values = [
        float(match.group(1))
        for line in lines
        for match in [re.search(r"S-SQUARED\s*=\s*([0-9.]+)", line)]
        if match
    ]
    if s2_values:
        metadata["s2"] = s2_values[-1]
    if "FAILURE TO LOCATE STATIONARY POINT" in "\n".join(lines):
        metadata["optimizer"] = "not_converged_nstep_limit"
    else:
        metadata["optimizer"] = "converged_or_no_failure_marker"
    return chosen, coords, metadata


def parse_pdb(path: Path) -> tuple[list[str], dict[tuple[str, str, str, str], PdbAtom]]:
    lines = path.read_text(errors="replace").splitlines()
    atoms: dict[tuple[str, str, str, str], PdbAtom] = {}
    for index, line in enumerate(lines):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        padded = (line + " " * 80)[:80]
        raw_elem = padded[76:80].strip()
        elem = "".join(ch for ch in raw_elem if ch.isalpha())
        if not elem:
            elem = "".join(ch for ch in padded[12:16].strip() if ch.isalpha())[:1]
        atom = PdbAtom(
            line_index=index,
            line=padded,
            serial=int(padded[6:11]),
            name=padded[12:16].strip(),
            resn=padded[17:20].strip(),
            chain=padded[21].strip(),
            resi=padded[22:26].strip(),
            element="FE" if elem.upper() == "FE" else elem.upper(),
            xyz=(float(padded[30:38]), float(padded[38:46]), float(padded[46:54])),
        )
        atoms[atom.key] = atom
    return lines, atoms


def key_from_note(note: str) -> tuple[str, str, str, str] | None:
    token = note.split()[0] if note.split() else ""
    parts = token.split(":")
    if len(parts) != 3:
        return None
    resn, resi, atom = parts
    if resn == "FE":
        return ("FE2", "A", resi, atom)
    if resn in {"H1", "H2"} and atom in {"CB", "CG", "ND1", "CD2", "CE1", "NE2"}:
        return ("HIS", "A", resi, atom)
    if resn == "EAC" and atom in {"CG", "CD", "OE1", "OE2"}:
        return ("GLU", "A", resi, atom)
    if resn == "ACT":
        return ("ACT", "A", resi, atom)
    if resn == "UNL":
        return ("UNL", "A", resi, atom)
    if resn == "HOH" and atom == "O":
        return ("HOH", "B", resi, atom)
    return None


def atom_label(key: tuple[str, str, str, str]) -> str:
    return f"{key[0]} {key[1]}:{key[2]} {key[3]}"


def replace_xyz(line: str, xyz: Vec) -> str:
    padded = (line + " " * 80)[:80]
    return f"{padded[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{padded[54:]}"


def write_model_pdb(path: Path, atoms: list[ModelAtom], coords: list[Vec]) -> None:
    lines: list[str] = []
    for serial, (atom, xyz) in enumerate(zip(atoms, coords), start=1):
        token = atom.note.split()[0] if atom.note.split() else f"UNK:1:{atom.element}"
        parts = token.split(":")
        resn, resi, name = (parts + ["UNK", "1", atom.element])[:3]
        chain = "A"
        if resn == "HOH":
            chain = "B"
        if resn in {"H1", "H2"}:
            resn = resn
        element = "Fe" if atom.element.upper() == "FE" else atom.element.capitalize()
        lines.append(
            f"HETATM{serial:5d} {name[:4]:>4s} {resn[:3]:>3s} {chain:1s}"
            f"{int(resi):4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            f"  1.00  0.00          {element:>2s}"
        )
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")


def write_model_xyz(path: Path, atoms: list[ModelAtom], coords: list[Vec], title: str) -> None:
    lines = [str(len(atoms)), title]
    for atom, xyz in zip(atoms, coords):
        element = "Fe" if atom.element.upper() == "FE" else atom.element.capitalize()
        suffix = f"  # {atom.note}" if atom.note else ""
        lines.append(f"{element:<2s} {xyz[0]:14.6f} {xyz[1]:14.6f} {xyz[2]:14.6f}{suffix}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_patched_full_pdb(
    path: Path,
    lines: list[str],
    full_atoms: dict[tuple[str, str, str, str], PdbAtom],
    patch_coords: dict[tuple[str, str, str, str], Vec],
    *,
    remove_keys: set[tuple[str, str, str, str]],
) -> None:
    patch_by_index = {full_atoms[key].line_index: xyz for key, xyz in patch_coords.items() if key in full_atoms}
    remove_indices = {full_atoms[key].line_index for key in remove_keys if key in full_atoms}
    out: list[str] = []
    for index, line in enumerate(lines):
        if index in remove_indices:
            continue
        if index in patch_by_index:
            out.append(replace_xyz(line, patch_by_index[index]))
        else:
            out.append(line)
    path.write_text("\n".join(out).rstrip() + "\n", encoding="utf-8")


def dist_for(keys: dict[tuple[str, str, str, str], Vec], donor: tuple[str, str, str, str]) -> float | None:
    fe = keys.get(FE_KEY)
    point = keys.get(donor)
    if fe is None or point is None:
        return None
    return distance(fe, point)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", required=True, help="GAMESS spin-screen directory")
    parser.add_argument("--log", default="logs/FEIII_radical_M7_fastopt.log")
    parser.add_argument("--initial-xyz", default="feiii_radical_truncated_model.xyz")
    parser.add_argument(
        "--full-pdb",
        default=r"E:\TJ\feiii_radical_preflight\complex_AFeIII_UNL_C4radical_dropH4A.pdb",
        help="Full protein PDB used for local coordinate grafting",
    )
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--nserch", type=int, default=None)
    args = parser.parse_args()

    workdir = Path(args.workdir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    initial_atoms = parse_xyz(workdir / args.initial_xyz)
    nserch0, coords0, _ = parse_coords_for_nserch(workdir / args.log, len(initial_atoms), 0)
    chosen_nserch, opt_atoms, metadata = parse_coords_for_nserch(workdir / args.log, len(initial_atoms), args.nserch)
    if len(opt_atoms) != len(initial_atoms):
        raise SystemExit("Optimized atom count does not match initial model atom count")

    rot, trans, align_rmsd = best_rigid_transform(
        [atom.xyz for atom in coords0],
        [atom.xyz for atom in initial_atoms],
    )
    aligned_coords = [apply_transform(atom.xyz, rot, trans) for atom in opt_atoms]

    write_model_pdb(outdir / "M7_optimized_small_model_gamess_frame.pdb", initial_atoms, [atom.xyz for atom in opt_atoms])
    write_model_xyz(
        outdir / "M7_optimized_small_model_gamess_frame.xyz",
        initial_atoms,
        [atom.xyz for atom in opt_atoms],
        f"M7 optimized GAMESS frame NSERCH={chosen_nserch}",
    )
    write_model_pdb(outdir / "M7_optimized_small_model_aligned_to_initial.pdb", initial_atoms, aligned_coords)
    write_model_xyz(
        outdir / "M7_optimized_small_model_aligned_to_initial.xyz",
        initial_atoms,
        aligned_coords,
        f"M7 optimized aligned to initial PDB frame NSERCH={chosen_nserch}",
    )

    full_lines, full_atoms = parse_pdb(Path(args.full_pdb))
    patch_coords: dict[tuple[str, str, str, str], Vec] = {}
    initial_key_coords: dict[tuple[str, str, str, str], Vec] = {}
    optimized_key_coords: dict[tuple[str, str, str, str], Vec] = {}
    mapping_rows: list[tuple[tuple[str, str, str, str], Vec, Vec, float]] = []
    for model_atom, aligned in zip(initial_atoms, aligned_coords):
        key = key_from_note(model_atom.note)
        if key is None:
            continue
        optimized_key_coords[key] = aligned
        if key in full_atoms:
            initial_key_coords[key] = full_atoms[key].xyz
            patch_coords[key] = aligned
            mapping_rows.append((key, full_atoms[key].xyz, aligned, distance(full_atoms[key].xyz, aligned)))

    remove_keys = {("UNL", "A", "1", "H4A")}
    patched_pdb = outdir / "complex_AFeIII_UNL_C4radical_M7_metalcenter_patch.pdb"
    write_patched_full_pdb(patched_pdb, full_lines, full_atoms, patch_coords, remove_keys=remove_keys)
    alias_pdb = outdir / "0616_initial_guess_ACT_M7_metalcenter_patch.pdb"
    alias_pdb.write_text(patched_pdb.read_text(encoding="utf-8"), encoding="utf-8")

    compare_tsv = outdir / "M7_initial_vs_optimized_compare.tsv"
    with compare_tsv.open("w", encoding="utf-8", newline="") as handle:
        handle.write("section\tatom\tinitial_x\tinitial_y\tinitial_z\tm7_x\tm7_y\tm7_z\tdisplacement_A\n")
        for key, initial, optimized, delta in sorted(mapping_rows):
            handle.write(
                f"mapped_atom\t{atom_label(key)}\t"
                f"{initial[0]:.6f}\t{initial[1]:.6f}\t{initial[2]:.6f}\t"
                f"{optimized[0]:.6f}\t{optimized[1]:.6f}\t{optimized[2]:.6f}\t{delta:.6f}\n"
            )
        handle.write("section\tatom\tinitial_distance_A\tm7_distance_A\tdelta_A\n")
        for donor in DONOR_KEYS:
            initial_dist = dist_for(initial_key_coords, donor)
            optimized_dist = dist_for(optimized_key_coords, donor)
            if initial_dist is None or optimized_dist is None:
                continue
            handle.write(
                f"fe_donor_distance\t{atom_label(donor)}\t"
                f"{initial_dist:.6f}\t{optimized_dist:.6f}\t{optimized_dist - initial_dist:.6f}\n"
            )

    donor_lines: list[str] = []
    for donor in DONOR_KEYS:
        initial_dist = dist_for(initial_key_coords, donor)
        optimized_dist = dist_for(optimized_key_coords, donor)
        if initial_dist is None or optimized_dist is None:
            continue
        donor_lines.append(
            f"| {atom_label(donor)} | {initial_dist:.3f} | {optimized_dist:.3f} | {optimized_dist - initial_dist:+.3f} |"
        )

    top_moves = sorted(mapping_rows, key=lambda row: row[3], reverse=True)[:15]
    move_lines = [
        f"| {atom_label(key)} | {delta:.3f} | {initial[0]:.3f},{initial[1]:.3f},{initial[2]:.3f} | {optimized[0]:.3f},{optimized[1]:.3f},{optimized[2]:.3f} |"
        for key, initial, optimized, delta in top_moves
    ]
    report = outdir / "M7_export_report.md"
    if metadata.get("optimizer") == "not_converged_nstep_limit":
        geometry_note = (
            "Important: the final failure block in the GAMESS log contains a next "
            f"predicted geometry with unknown energy/gradient. This export uses the "
            f"evaluated `NSERCH={chosen_nserch}` geometry instead."
        )
    else:
        geometry_note = (
            f"Important: this export uses the converged/evaluated `NSERCH={chosen_nserch}` "
            "geometry from the GAMESS log."
        )
    report.write_text(
        "\n".join(
            [
                "# M7 Optimized Geometry Export",
                "",
                f"- Source full PDB: `{Path(args.full_pdb)}`",
                f"- GAMESS log: `{workdir / args.log}`",
                f"- Evaluated geometry exported: `NSERCH={chosen_nserch}`",
                f"- Energy: `{metadata.get('energy_hartree', 'NA')}` Hartree",
                f"- S-squared: `{metadata.get('s2', 'NA')}`",
                f"- Gradient max/RMS: `{metadata.get('grad_max', 'NA')}` / `{metadata.get('grad_rms', 'NA')}`",
                f"- Optimizer status: `{metadata.get('optimizer', 'NA')}`",
                f"- GAMESS NSERCH=0 to initial-frame alignment RMSD: `{align_rmsd:.6f} A`",
                "",
                geometry_note,
                "",
                "## Output Files",
                "",
                "- `M7_optimized_small_model_gamess_frame.pdb`: evaluated QM model in GAMESS internal frame.",
                "- `M7_optimized_small_model_aligned_to_initial.pdb`: evaluated QM model aligned back to the initial PDB frame.",
                "- `complex_AFeIII_UNL_C4radical_M7_metalcenter_patch.pdb`: full protein with safely mapped metal-center atoms replaced by M7 aligned coordinates.",
                "- `0616_initial_guess_ACT_M7_metalcenter_patch.pdb`: same full-protein patch with a convenient name.",
                "- `M7_initial_vs_optimized_compare.tsv`: atom displacement and Fe-donor distance table.",
                "",
                "## Fe-Donor Distance Changes",
                "",
                "| donor | initial A | M7 aligned A | delta A |",
                "|---|---:|---:|---:|",
                *donor_lines,
                "",
                "## Largest Mapped Atom Movements",
                "",
                "| atom | displacement A | initial xyz | M7 aligned xyz |",
                "|---|---:|---|---|",
                *move_lines,
                "",
                "## Caveats",
                "",
                "- This is a local QM geometry graft, not a full-protein quantum optimization.",
                "- His/Glu protein atoms are only patched for atoms that exist in both the full PDB and the truncated QM model; methyl-cap hydrogens are not inserted into the protein.",
                "- `UNL A 1 H4A` is removed in the radical preflight model and remains absent in the patched full PDB.",
                "- Before production MD, use the M7 small-model Hessian/MCPB.py workflow to generate bonded Fe parameters, then run classical minimization in Amber/GROMACS.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"outdir={outdir}")
    print(f"nserch={chosen_nserch}")
    print(f"alignment_rmsd_A={align_rmsd:.6f}")
    print(f"small_model_aligned={outdir / 'M7_optimized_small_model_aligned_to_initial.pdb'}")
    print(f"patched_full_pdb={patched_pdb}")
    print(f"compare_tsv={compare_tsv}")
    print(f"report={report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
