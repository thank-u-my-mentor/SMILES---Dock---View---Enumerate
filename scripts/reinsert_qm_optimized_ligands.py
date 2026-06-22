from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class AtomLine:
    line: str
    serial: int
    name: str
    resn: str
    chain: str
    resi: str
    elem: str
    xyz: tuple[float, float, float]


def parse_pdb_atoms(path: Path) -> list[AtomLine]:
    atoms: list[AtomLine] = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        elem = (line[76:78].strip() or "".join(ch for ch in line[12:16].strip() if ch.isalpha())[:1]).capitalize()
        atoms.append(
            AtomLine(
                line=line,
                serial=int(line[6:11]),
                name=line[12:16].strip(),
                resn=line[17:20].strip(),
                chain=line[21].strip(),
                resi=line[22:26].strip(),
                elem=elem,
                xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
            )
        )
    return atoms


def v_add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def v_sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def v_mul(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def centroid(points: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    n = float(len(points))
    return (sum(p[0] for p in points) / n, sum(p[1] for p in points) / n, sum(p[2] for p in points) / n)


def mat_vec(m, v):
    return (
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    )


def fit_transform(mobile: list[tuple[float, float, float]], target: list[tuple[float, float, float]]):
    # Use numpy if available; the MD/base environments have it, and it keeps the
    # alignment code small and auditable.
    import numpy as np

    mob = np.array(mobile, dtype=float)
    ref = np.array(target, dtype=float)
    mob_c = mob.mean(axis=0)
    ref_c = ref.mean(axis=0)
    mob0 = mob - mob_c
    ref0 = ref - ref_c
    h = mob0.T @ ref0
    u, _s, vt = np.linalg.svd(h)
    r = vt.T @ u.T
    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = vt.T @ u.T
    aligned = (mob0 @ r.T) + ref_c
    rmsd = float(np.sqrt(((aligned - ref) ** 2).sum() / len(mobile)))
    return r.tolist(), tuple(mob_c), tuple(ref_c), rmsd


def apply_transform(xyz, rotation, mobile_centroid, target_centroid):
    shifted = v_sub(xyz, mobile_centroid)
    rotated = mat_vec(rotation, shifted)
    return v_add(rotated, target_centroid)


def replace_xyz(line: str, xyz: tuple[float, float, float]) -> str:
    padded = (line + " " * 80)[:80]
    return f"{padded[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{padded[54:]}"


def select_residue(atoms: list[AtomLine], resn: str, chain: str | None = None, resi: str | None = None) -> list[AtomLine]:
    out = [a for a in atoms if a.resn == resn]
    if chain is not None:
        out = [a for a in out if a.chain == chain]
    if resi is not None:
        out = [a for a in out if a.resi.strip() == str(resi).strip()]
    return out


def build_name_map(atoms: list[AtomLine]) -> dict[str, AtomLine]:
    return {a.name: a for a in atoms}


def main() -> int:
    if len(sys.argv) != 6:
        print(
            "usage: reinsert_qm_optimized_ligands.py complex.pdb act_opt.pdb unl_opt.pdb out_complex.pdb report.txt"
        )
        return 2

    complex_pdb = Path(sys.argv[1])
    act_opt_pdb = Path(sys.argv[2])
    unl_opt_pdb = Path(sys.argv[3])
    out_complex = Path(sys.argv[4])
    report = Path(sys.argv[5])
    out_complex.parent.mkdir(parents=True, exist_ok=True)
    report.parent.mkdir(parents=True, exist_ok=True)

    complex_atoms = parse_pdb_atoms(complex_pdb)
    act_ref = select_residue(complex_atoms, "ACT", "A", "501")
    unl_ref = select_residue(complex_atoms, "UNL", "A", "1")
    act_opt = parse_pdb_atoms(act_opt_pdb)
    unl_opt = parse_pdb_atoms(unl_opt_pdb)

    replacements: dict[int, tuple[float, float, float]] = {}
    report_lines = []
    for label, ref_atoms, opt_atoms in [
        ("ACT", act_ref, act_opt),
        ("UNL", unl_ref, unl_opt),
    ]:
        ref_by_name = build_name_map(ref_atoms)
        opt_by_name = build_name_map(opt_atoms)
        common_heavy = [
            name
            for name in opt_by_name
            if name in ref_by_name and opt_by_name[name].elem.upper() != "H"
        ]
        if len(common_heavy) < 3:
            raise ValueError(f"{label}: not enough common heavy atoms for alignment: {common_heavy}")
        rotation, mob_c, ref_c, rmsd = fit_transform(
            [opt_by_name[name].xyz for name in common_heavy],
            [ref_by_name[name].xyz for name in common_heavy],
        )
        report_lines.append(f"{label}\tcommon_heavy={','.join(common_heavy)}\talignment_rmsd={rmsd:.4f}")
        for name, opt_atom in opt_by_name.items():
            if name not in ref_by_name:
                continue
            replacements[ref_by_name[name].serial] = apply_transform(opt_atom.xyz, rotation, mob_c, ref_c)

    out_lines = []
    for line in complex_pdb.read_text(errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            serial = int(line[6:11])
            if serial in replacements:
                line = replace_xyz(line, replacements[serial])
        out_lines.append(line)
    out_complex.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    report.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(f"out_complex={out_complex}")
    print(f"report={report}")
    print("\n".join(report_lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
