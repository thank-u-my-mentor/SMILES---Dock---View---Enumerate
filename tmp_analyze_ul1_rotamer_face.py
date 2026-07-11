#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import calc_dihedrals


TOP = Path("/mnt/e/TJ/260706_mcpb_m2_200ns/md_200ns/pbc/md_200ns_prot_mcpb_center_start.gro")
TRAJ = Path("/mnt/e/TJ/260706_mcpb_m2_200ns/md_200ns/pbc/md_200ns_prot_mcpb_fit.xtc")
OUTDIR = Path("/mnt/e/TJ/260706_mcpb_m2_200ns/md_200ns/analysis_reactive_C4O1_20260708/ul1_rotamer_face_analysis")


def one(u: mda.Universe, selection: str):
    ag = u.select_atoms(selection)
    if len(ag) != 1:
        raise SystemExit(f"Selection {selection!r} matched {len(ag)} atoms, expected 1")
    return ag[0]


def dihedral_deg(a, b, c, d) -> float:
    val = calc_dihedrals(
        a.position[None, :],
        b.position[None, :],
        c.position[None, :],
        d.position[None, :],
    )[0]
    return float(np.degrees(val))


def norm(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n == 0:
        return v
    return v / n


def signed_height_to_c3_c4_c5_plane(c3, c4, c5, point) -> float:
    # Plane normal with origin at C4; sign depends on atom order, but same/opposite face does not.
    normal = norm(np.cross(c3.position - c4.position, c5.position - c4.position))
    return float(np.dot(point.position - c4.position, normal))


def circular_r(deg: np.ndarray) -> float:
    rad = np.deg2rad(deg)
    return float(np.hypot(np.mean(np.cos(rad)), np.mean(np.sin(rad))))


def circular_mean(deg: np.ndarray) -> float:
    rad = np.deg2rad(deg)
    return float(np.rad2deg(np.arctan2(np.mean(np.sin(rad)), np.mean(np.cos(rad)))))


def rotamer_bin(deg: float) -> str:
    if -60.0 <= deg < 60.0:
        return "around_0"
    if 60.0 <= deg <= 180.0:
        return "positive_60_to_180"
    return "negative_minus180_to_minus60"


def write_selected_pdbs(u: mda.Universe, selected: pd.DataFrame, outdir: Path) -> None:
    ul1_dir = outdir / "selected_ul1_pdbs"
    site_dir = outdir / "selected_active_site_pdbs"
    ul1_dir.mkdir(parents=True, exist_ok=True)
    site_dir.mkdir(parents=True, exist_ok=True)
    ul1 = u.select_atoms("resname UL1")
    site = u.select_atoms("resname UL1 FE1 HD1 HD2 GU1 AT1 HH1 or resid 322 333 255 256")
    for row in selected.itertuples(index=False):
        u.trajectory[int(row.frame)]
        tag = f"frame{int(row.frame):05d}_time{float(row.time_ns):07.3f}ns_C4O1_{float(row.C4_O1_A):.2f}_tauC3C4_{float(row.tau_C2_C3_C4_C5_deg):.1f}_{row.O1_H02_face_relation}"
        ul1.write(str(ul1_dir / f"{tag}_UL1.pdb"))
        site.write(str(site_dir / f"{tag}_active_site.pdb"))


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name, sub in [
        ("all_frames", df),
        ("C4_O1_le_3A", df[df["C4_O1_A"] <= 3.0]),
        ("C4_O1_le_3p5A", df[df["C4_O1_A"] <= 3.5]),
    ]:
        if len(sub) == 0:
            continue
        tau = sub["tau_C2_C3_C4_C5_deg"].to_numpy(float)
        tau45 = sub["tau_C3_C4_C5_C6_deg"].to_numpy(float)
        rows.append(
            {
                "subset": name,
                "n_frames": len(sub),
                "fraction_O1_H02_same_face": float((sub["O1_H02_face_relation"] == "same").mean()),
                "fraction_O1_H02_opposite_face": float((sub["O1_H02_face_relation"] == "opposite").mean()),
                "C4_O1_median_A": float(sub["C4_O1_A"].median()),
                "tau_C2_C3_C4_C5_median_deg": float(np.median(tau)),
                "tau_C2_C3_C4_C5_p05_deg": float(np.percentile(tau, 5)),
                "tau_C2_C3_C4_C5_p95_deg": float(np.percentile(tau, 95)),
                "tau_C2_C3_C4_C5_circular_mean_deg": circular_mean(tau),
                "tau_C2_C3_C4_C5_circular_R": circular_r(tau),
                "tau_C3_C4_C5_C6_median_deg": float(np.median(tau45)),
                "tau_C3_C4_C5_C6_p05_deg": float(np.percentile(tau45, 5)),
                "tau_C3_C4_C5_C6_p95_deg": float(np.percentile(tau45, 95)),
                "tau_C3_C4_C5_C6_circular_mean_deg": circular_mean(tau45),
                "tau_C3_C4_C5_C6_circular_R": circular_r(tau45),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    u = mda.Universe(str(TOP), str(TRAJ))
    atoms = {
        name: one(u, f"resname UL1 and name {name}")
        for name in ["C1", "O1", "C2", "C3", "C4", "C5", "C6", "C10", "H02"]
    }
    rows = []
    for ts in u.trajectory:
        c3, c4, c5 = atoms["C3"], atoms["C4"], atoms["C5"]
        h_o1 = signed_height_to_c3_c4_c5_plane(c3, c4, c5, atoms["O1"])
        h_h02 = signed_height_to_c3_c4_c5_plane(c3, c4, c5, atoms["H02"])
        product = h_o1 * h_h02
        rows.append(
            {
                "frame": int(ts.frame),
                "time_ps": float(ts.time),
                "time_ns": float(ts.time / 1000.0),
                "C4_O1_A": float(np.linalg.norm(atoms["C4"].position - atoms["O1"].position)),
                "tau_C2_C3_C4_C5_deg": dihedral_deg(atoms["C2"], atoms["C3"], atoms["C4"], atoms["C5"]),
                "tau_H02_C4_C3_C2_deg": dihedral_deg(atoms["H02"], atoms["C4"], atoms["C3"], atoms["C2"]),
                "tau_C3_C4_C5_C6_deg": dihedral_deg(atoms["C3"], atoms["C4"], atoms["C5"], atoms["C6"]),
                "tau_C3_C4_C5_C10_deg": dihedral_deg(atoms["C3"], atoms["C4"], atoms["C5"], atoms["C10"]),
                "phi_C1_O1_C4_C5_deg": dihedral_deg(atoms["C1"], atoms["O1"], atoms["C4"], atoms["C5"]),
                "h02_C1_O1_C4_H02_deg": dihedral_deg(atoms["C1"], atoms["O1"], atoms["C4"], atoms["H02"]),
                "O1_height_to_C3C4C5_plane_A": h_o1,
                "H02_height_to_C3C4C5_plane_A": h_h02,
                "O1_H02_face_product_A2": product,
                "O1_H02_face_relation": "same" if product > 0 else "opposite",
            }
        )
    df = pd.DataFrame(rows)
    df["tau_C2_C3_C4_C5_rotamer_bin"] = df["tau_C2_C3_C4_C5_deg"].map(rotamer_bin)
    df.to_csv(OUTDIR / "ul1_c3c4_c4c5_face_timeseries.csv", index=False)

    summary = summarize(df)
    summary.to_csv(OUTDIR / "ul1_rotamer_face_summary.csv", index=False)

    occ = (
        df.assign(near3=df["C4_O1_A"] <= 3.0, near35=df["C4_O1_A"] <= 3.5)
        .groupby(["near3", "near35", "tau_C2_C3_C4_C5_rotamer_bin", "O1_H02_face_relation"])
        .size()
        .reset_index(name="n_frames")
    )
    occ.to_csv(OUTDIR / "ul1_rotamer_face_occupancy.csv", index=False)

    # Select the shortest C4-O1 examples from each face/rotamer bin plus global shortest frames.
    candidates = []
    near = df[df["C4_O1_A"] <= 3.0].copy()
    for _, sub in near.groupby(["O1_H02_face_relation", "tau_C2_C3_C4_C5_rotamer_bin"]):
        candidates.append(sub.nsmallest(2, "C4_O1_A"))
    candidates.append(df.nsmallest(6, "C4_O1_A"))
    selected = pd.concat(candidates, ignore_index=True).drop_duplicates("frame").nsmallest(14, "C4_O1_A")
    selected.to_csv(OUTDIR / "selected_ul1_rotamer_face_frames.csv", index=False)
    write_selected_pdbs(u, selected, OUTDIR)

    print("Wrote", OUTDIR)
    print(summary.to_string(index=False))
    print("\\nNear <=3A face counts:")
    print(near["O1_H02_face_relation"].value_counts().to_string())
    print("\\nNear <=3A C3-C4 rotamer counts:")
    print(near["tau_C2_C3_C4_C5_rotamer_bin"].value_counts().to_string())


if __name__ == "__main__":
    main()
