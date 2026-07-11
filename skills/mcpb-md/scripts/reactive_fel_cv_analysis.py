#!/usr/bin/env python3
"""Compute mechanistic CVs and 2D free-energy landscapes for reactive MD."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import calc_angles, calc_dihedrals


K_B_KCAL = 0.00198720425864083
PHE_RING_NAMES = ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--top", required=True)
    parser.add_argument("--traj", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--phe-resid", type=int, default=322)
    return parser.parse_args()


def one_atom(u: mda.Universe, selection: str):
    ag = u.select_atoms(selection)
    if len(ag) != 1:
        raise SystemExit(f"Selection {selection!r} matched {len(ag)} atoms, expected 1")
    return ag[0]


def atoms_by_names(u: mda.Universe, selection_prefix: str, names: list[str]):
    atoms = []
    for name in names:
        atoms.append(one_atom(u, f"{selection_prefix} and name {name}"))
    return atoms


def dist(a, b) -> float:
    return float(np.linalg.norm(a.position - b.position))


def angle_deg(a, b, c) -> float:
    return float(np.degrees(calc_angles(a.position[None, :], b.position[None, :], c.position[None, :])[0]))


def dihedral_deg(a, b, c, d) -> float:
    return float(np.degrees(calc_dihedrals(a.position[None, :], b.position[None, :], c.position[None, :], d.position[None, :])[0]))


def normalize(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    if n == 0:
        return v
    return v / n


def ring_geometry(ring_atoms, point: np.ndarray) -> tuple[float, float, float, float, np.ndarray, np.ndarray]:
    coords = np.array([atom.position for atom in ring_atoms], dtype=float)
    centroid = coords.mean(axis=0)
    # Plane normal from SVD is more stable than relying on a single cross product.
    _, _, vh = np.linalg.svd(coords - centroid)
    normal = normalize(vh[-1])
    vec = point - centroid
    centroid_distance = float(np.linalg.norm(vec))
    signed_height = float(np.dot(vec, normal))
    height_abs = abs(signed_height)
    lateral_offset = math.sqrt(max(centroid_distance * centroid_distance - height_abs * height_abs, 0.0))
    normal_angle = float(np.degrees(math.acos(max(-1.0, min(1.0, abs(signed_height) / centroid_distance))))) if centroid_distance else 0.0
    return centroid_distance, signed_height, lateral_offset, normal_angle, centroid, normal


def c4_face_height(c3, c4, c5, o1) -> float:
    normal = normalize(np.cross(c3.position - c4.position, c5.position - c4.position))
    return float(np.dot(o1.position - c4.position, normal))


def c4_plane_signed_height(c3, c4, c5, atom) -> float:
    normal = normalize(np.cross(c3.position - c4.position, c5.position - c4.position))
    return float(np.dot(atom.position - c4.position, normal))


def plane_signed_height(origin, atom_a, atom_b, point_atom) -> float:
    normal = normalize(np.cross(atom_a.position - origin.position, atom_b.position - origin.position))
    return float(np.dot(point_atom.position - origin.position, normal))


def point_plane_signed_height(origin, atom_a, atom_b, point: np.ndarray) -> float:
    normal = normalize(np.cross(atom_a.position - origin.position, atom_b.position - origin.position))
    return float(np.dot(point - origin.position, normal))


def signed_tetrahedral_volume(a, b, c, d) -> float:
    """Signed volume for four ligand atoms in priority order a,b,c,d.

    Here used as a pseudo-CIP descriptor for product-like C4 chirality:
    O1 > C5(phenyl) > C3(chain) > H02. The sign must be calibrated to
    R/S using an explicit product model.
    """
    mat = np.vstack(
        [
            a.position - d.position,
            b.position - d.position,
            c.position - d.position,
        ]
    )
    return float(np.linalg.det(mat) / 6.0)


def normal_angle_deg(vec_a: np.ndarray, vec_b: np.ndarray, absolute: bool = True) -> float:
    va = normalize(vec_a)
    vb = normalize(vec_b)
    dot = float(np.dot(va, vb))
    if absolute:
        dot = abs(dot)
    dot = max(-1.0, min(1.0, dot))
    return float(np.degrees(math.acos(dot)))


def compute_cvs(args: argparse.Namespace) -> pd.DataFrame:
    u = mda.Universe(args.top, args.traj)
    atoms = {
        "C1": one_atom(u, "resname UL1 and name C1"),
        "N1": one_atom(u, "resname UL1 and name N1"),
        "O1": one_atom(u, "resname UL1 and name O1"),
        "C2": one_atom(u, "resname UL1 and name C2"),
        "C3": one_atom(u, "resname UL1 and name C3"),
        "C4": one_atom(u, "resname UL1 and name C4"),
        "C5": one_atom(u, "resname UL1 and name C5"),
        "C6": one_atom(u, "resname UL1 and name C6"),
        "C10": one_atom(u, "resname UL1 and name C10"),
        "H02": one_atom(u, "resname UL1 and name H02"),
        "FE": one_atom(u, "resname FE1 and name FE"),
    }
    ring_atoms = atoms_by_names(u, f"resid {args.phe_resid} and resname PHE", PHE_RING_NAMES)

    rows = []
    for ts in u.trajectory[:: args.stride]:
        c4 = atoms["C4"]
        o1 = atoms["O1"]
        ring_dist, ring_height, ring_lateral, ring_angle, phe_centroid, phe_normal = ring_geometry(ring_atoms, c4.position)
        face_h = c4_face_height(atoms["C3"], c4, atoms["C5"], o1)
        h02_face_h = c4_plane_signed_height(atoms["C3"], c4, atoms["C5"], atoms["H02"])
        lactone_plane_h02 = plane_signed_height(c4, o1, atoms["C1"], atoms["H02"])
        lactone_plane_c5 = plane_signed_height(c4, o1, atoms["C1"], atoms["C5"])
        lactone_plane_c3 = plane_signed_height(c4, o1, atoms["C1"], atoms["C3"])
        lactone_plane_phe_centroid = point_plane_signed_height(c4, o1, atoms["C1"], phe_centroid)
        c4_radical_normal = normalize(np.cross(atoms["C3"].position - c4.position, atoms["C5"].position - c4.position))
        rows.append(
            {
                "frame": int(ts.frame),
                "time_ps": float(ts.time),
                "time_ns": float(ts.time / 1000.0),
                "C4_O1_A": dist(c4, o1),
                "C4_face_signed_height_A": face_h,
                "C4_face_abs_height_A": abs(face_h),
                "H02_face_signed_height_A": h02_face_h,
                "H02_face_abs_height_A": abs(h02_face_h),
                "O1_H02_same_face_product_A2": face_h * h02_face_h,
                "O1_H02_same_face_binary": int(np.sign(face_h) == np.sign(h02_face_h)),
                "pseudo_CIP_O1_C5_C3_H02_chiral_volume_A3": signed_tetrahedral_volume(o1, atoms["C5"], atoms["C3"], atoms["H02"]),
                "lactone_plane_C4_O1_C1_H02_signed_height_A": lactone_plane_h02,
                "lactone_plane_C4_O1_C1_C5_signed_height_A": lactone_plane_c5,
                "lactone_plane_C4_O1_C1_C3_signed_height_A": lactone_plane_c3,
                "lactone_plane_C4_O1_C1_Phe_centroid_signed_height_A": lactone_plane_phe_centroid,
                "lactone_plane_H02_C5_same_side_product_A2": lactone_plane_h02 * lactone_plane_c5,
                "dihedral_C3_C4_C5_O1_deg": dihedral_deg(atoms["C3"], c4, atoms["C5"], o1),
                "dihedral_C3_C4_O1_C1_deg": dihedral_deg(atoms["C3"], c4, o1, atoms["C1"]),
                "dihedral_C1_O1_C4_C5_deg": dihedral_deg(atoms["C1"], o1, c4, atoms["C5"]),
                "dihedral_C1_O1_C4_H02_deg": dihedral_deg(atoms["C1"], o1, c4, atoms["H02"]),
                "dihedral_C1_O1_C4_C3_deg": dihedral_deg(atoms["C1"], o1, c4, atoms["C3"]),
                "dihedral_O1_C4_C5_C6_deg": dihedral_deg(o1, c4, atoms["C5"], atoms["C6"]),
                "dihedral_C3_C4_C5_C6_deg": dihedral_deg(atoms["C3"], c4, atoms["C5"], atoms["C6"]),
                "dihedral_O1_C1_C2_C3_deg": dihedral_deg(o1, atoms["C1"], atoms["C2"], atoms["C3"]),
                "dihedral_C1_C2_C3_C4_deg": dihedral_deg(atoms["C1"], atoms["C2"], atoms["C3"], c4),
                "dihedral_N1_C1_C2_C3_deg": dihedral_deg(atoms["N1"], atoms["C1"], atoms["C2"], atoms["C3"]),
                "dihedral_FE_N1_C1_O1_deg": dihedral_deg(atoms["FE"], atoms["N1"], atoms["C1"], o1),
                "dihedral_FE_N1_C1_C2_deg": dihedral_deg(atoms["FE"], atoms["N1"], atoms["C1"], atoms["C2"]),
                "angle_C3_C4_O1_deg": angle_deg(atoms["C3"], c4, o1),
                "angle_C5_C4_O1_deg": angle_deg(atoms["C5"], c4, o1),
                "angle_C4_O1_C1_deg": angle_deg(c4, o1, atoms["C1"]),
                "PHE322_ring_centroid_to_C4_A": ring_dist,
                "PHE322_C4_signed_height_to_ring_A": ring_height,
                "PHE322_C4_abs_height_to_ring_A": abs(ring_height),
                "PHE322_C4_lateral_offset_A": ring_lateral,
                "PHE322_C4_ring_normal_angle_deg": ring_angle,
                "C4_radical_plane_normal_vs_PHE322_ring_normal_angle_deg": normal_angle_deg(c4_radical_normal, phe_normal, absolute=True),
                "FE_C4_A": dist(atoms["FE"], c4),
                "FE_O1_A": dist(atoms["FE"], o1),
            }
        )
    return pd.DataFrame(rows)


def free_energy_grid(x, y, bins=80, temperature=300.0):
    hist, xedges, yedges = np.histogram2d(x, y, bins=bins)
    prob = hist / hist.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        fe = -K_B_KCAL * temperature * np.log(prob)
    finite = np.isfinite(fe)
    fe[finite] -= np.nanmin(fe[finite])
    fe[~finite] = np.nan
    xc = 0.5 * (xedges[:-1] + xedges[1:])
    yc = 0.5 * (yedges[:-1] + yedges[1:])
    return xc, yc, fe.T


def plot_fel(df: pd.DataFrame, xcol: str, ycol: str, out_png: Path, title: str, temperature: float) -> None:
    x = df[xcol].to_numpy(dtype=float)
    y = df[ycol].to_numpy(dtype=float)
    xc, yc, fe = free_energy_grid(x, y, bins=85, temperature=temperature)
    fig, ax = plt.subplots(figsize=(6.5, 5.4))
    im = ax.contourf(xc, yc, fe, levels=np.linspace(0, np.nanpercentile(fe, 95), 24), cmap="viridis_r")
    ax.scatter(x, y, s=1.0, c="k", alpha=0.05, linewidths=0)
    ax.set_xlabel(xcol)
    ax.set_ylabel(ycol)
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Free energy (kcal/mol, relative)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=240)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    df = compute_cvs(args)
    df.to_csv(outdir / "reactive_geometry_CVs.csv", index=False)

    plot_specs = [
        ("C4_O1_A", "C4_face_signed_height_A", "fel_C4O1_vs_signed_face_height.png", "Reaction distance vs O1 approach face"),
        ("C4_O1_A", "H02_face_signed_height_A", "fel_C4O1_vs_H02_face_height.png", "Reaction distance vs benzylic H face"),
        ("C4_O1_A", "O1_H02_same_face_product_A2", "fel_C4O1_vs_O1_H02_same_face_product.png", "Reaction distance vs O1/H02 face relation"),
        ("C4_O1_A", "pseudo_CIP_O1_C5_C3_H02_chiral_volume_A3", "fel_C4O1_vs_pseudo_CIP_chiral_volume.png", "Reaction distance vs pseudo-CIP chiral volume"),
        ("C4_O1_A", "lactone_plane_H02_C5_same_side_product_A2", "fel_C4O1_vs_lactone_plane_H02_C5_product.png", "Reaction distance vs H/Ph side of C4-O1-C1 plane"),
        ("dihedral_C1_O1_C4_C5_deg", "dihedral_C1_O1_C4_H02_deg", "fel_C1O1C4C5_vs_C1O1C4H02.png", "Product-like O1-C4 dihedrals"),
        ("C4_O1_A", "dihedral_C3_C4_C5_O1_deg", "fel_C4O1_vs_face_dihedral.png", "Reaction distance vs face dihedral"),
        ("C4_O1_A", "PHE322_ring_centroid_to_C4_A", "fel_C4O1_vs_PHE322_ring_centroid_C4.png", "Reaction distance vs PHE322-C4 centroid distance"),
        ("PHE322_C4_abs_height_to_ring_A", "PHE322_C4_lateral_offset_A", "fel_PHE322_height_vs_lateral_offset.png", "C4 position relative to PHE322 ring plane"),
        ("dihedral_FE_N1_C1_O1_deg", "dihedral_O1_C1_C2_C3_deg", "fel_FeN1C1O1_vs_O1C1C2C3.png", "Fe-bound carbonyl orientation vs chain fold"),
    ]
    for xcol, ycol, png, title in plot_specs:
        plot_fel(df, xcol, ycol, outdir / png, title, args.temperature)

    summary = [
        "# Reactive CV/FEL Analysis",
        "",
        "Recommended CVs:",
        "- `C4_O1_A`: reaction distance for C4 to O1 lactonization.",
        "- `C4_face_signed_height_A`: signed O1 approach side relative to the C3-C4-C5 local plane. Its sign can be calibrated to pro-R/pro-S using a product/QM structure.",
        "- `pseudo_CIP_O1_C5_C3_H02_chiral_volume_A3`: product-like C4 signed volume with tentative CIP order O1 > C5(Ph) > C3(chain) > H02. Calibrate the sign to R/S before naming.",
        "- `lactone_plane_C4_O1_C1_*`: side of H02, C5(Ph), C3, and phenyl centroid relative to the product-like C4-O1-C1 plane.",
        "- `PHE322_ring_centroid_to_C4_A`, `PHE322_C4_abs_height_to_ring_A`, and `PHE322_C4_lateral_offset_A`: pi-radical contact diagnostics.",
        "- `dihedral_FE_N1_C1_O1_deg` and `dihedral_FE_N1_C1_C2_deg`: Fe-bound ligand orientation controls, not direct R/S labels.",
        "",
        "Important: the sign of the face/chirality CV is a convention until calibrated against an actual product stereochemical assignment.",
        "",
        "PHE322 ring atoms used: CG, CD1, CE1, CZ, CE2, CD2.",
    ]
    (outdir / "README_reactive_CV_FEL.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(outdir)


if __name__ == "__main__":
    main()
