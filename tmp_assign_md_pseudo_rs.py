#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd


ROOT = Path("/mnt/e/TJ/260706_mcpb_m2_200ns")
TOP = ROOT / "md_200ns/pbc/md_200ns_prot_mcpb_center_start.gro"
TRAJ = ROOT / "md_200ns/pbc/md_200ns_prot_mcpb_fit.xtc"
OUTDIR = ROOT / "md_200ns/analysis_reactive_C4O1_20260708/ul1_rotamer_face_analysis/product_calibrated_pseudo_RS"


def one(u: mda.Universe, selection: str):
    ag = u.select_atoms(selection)
    if len(ag) != 1:
        raise SystemExit(f"{selection!r} matched {len(ag)} atoms, expected 1")
    return ag[0]


def signed_triple(center: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    return float(np.linalg.det(np.vstack([a - center, b - center, c - center])))


def signed_tetra(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> float:
    return float(np.linalg.det(np.vstack([a - d, b - d, c - d])))


def assignment(value: float) -> str:
    # Calibrated from RDKit product references:
    # R product: signed_triple(O1, C5_phenyl, C3_chain about C4) < 0
    # S product: signed_triple(O1, C5_phenyl, C3_chain about C4) > 0
    if value < 0:
        return "R_like_by_product_reference"
    if value > 0:
        return "S_like_by_product_reference"
    return "undefined"


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cutoff in [3.50, 3.25, 3.00, 2.90, 2.80, 2.75, 2.70]:
        sub = df[df["C4_O1_A"] <= cutoff]
        rows.append(row_summary(sub, f"C4_O1_le_{cutoff:.2f}A"))
    sorted_df = df.sort_values("C4_O1_A")
    for n in [5, 10, 20, 50, 100, 500, 1000]:
        rows.append(row_summary(sorted_df.head(n), f"top_{n}_shortest_C4_O1"))
    return pd.DataFrame(rows)


def write_representative_pdbs_and_pml(u: mda.Universe, df: pd.DataFrame) -> None:
    reps = []
    near = df[df["C4_O1_A"] <= 3.0].copy()
    for label in ["R_like_by_product_reference", "S_like_by_product_reference"]:
        sub = near[near["pseudo_RS_by_product_ref"] == label].sort_values("C4_O1_A")
        if len(sub):
            reps.append(sub.iloc[0])
    if not reps:
        return

    rep_df = pd.DataFrame(reps)
    rep_df.to_csv(OUTDIR / "closest_R_like_vs_S_like_frames.csv", index=False)

    all_atoms = u.select_atoms("all")
    site = u.select_atoms("resname UL1 FE1 HD1 HD2 GU1 AT1 HH1 or resid 322 333 255 256")
    pdb_rows = []
    for row in rep_df.itertuples(index=False):
        u.trajectory[int(row.frame)]
        short = "R_like" if str(row.pseudo_RS_by_product_ref).startswith("R") else "S_like"
        tag = f"frame{int(row.frame):05d}_time{float(row.time_ns):07.3f}ns_C4O1_{float(row.C4_O1_A):.2f}_{short}"
        full_pdb = OUTDIR / f"{tag}_protein_mcpb.pdb"
        site_pdb = OUTDIR / f"{tag}_active_site.pdb"
        all_atoms.write(str(full_pdb))
        site.write(str(site_pdb))
        pdb_rows.append((short, tag, full_pdb, row))

    pml = OUTDIR / "view_closest_R_like_vs_S_like_near_attack.pml"
    with pml.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("reinitialize\n")
        handle.write("bg_color white\n")
        handle.write("set stick_radius, 0.15\n")
        handle.write("set sphere_scale, 0.22\n")
        handle.write("set dash_radius, 0.05\n")
        handle.write("set label_size, 16\n")
        handle.write("set label_color, black\n")
        handle.write("set cartoon_transparency, 0.72\n")
        for short, tag, pdb, row in pdb_rows:
            obj = f"{short}_frame{int(row.frame):05d}"
            carbon_color = "gray70" if short == "R_like" else "gray90"
            ul1_color = "cyan" if short == "R_like" else "salmon"
            handle.write(f"load {pdb.as_posix()}, {obj}\n")
            handle.write(f"dss {obj}\n")
            handle.write(f"hide everything, {obj}\n")
            handle.write(f"show cartoon, {obj} and polymer\n")
            handle.write(f"color {carbon_color}, {obj} and elem C\n")
            handle.write(f"show sticks, {obj} and (resn UL1+FE1+HD1+HD2+GU1+AT1+HH1 or resi 1+173+256+335+364+365+366)\n")
            handle.write(f"show spheres, {obj} and resn FE1\n")
            handle.write(f"color orange, {obj} and resn FE1\n")
            handle.write(f"color {ul1_color}, {obj} and resn UL1 and elem C\n")
            handle.write(f"color red, {obj} and resn UL1 and name O1\n")
            handle.write(f"color yellow, {obj} and resn UL1 and name C4\n")
            handle.write(f"color gray40, {obj} and resn UL1 and name H02\n")
            handle.write(f"distance {short}_C4_O1, {obj} and resn UL1 and name C4, {obj} and resn UL1 and name O1\n")
            handle.write(f"label {obj} and resn UL1 and name C4, \"{short}; frame {int(row.frame)}; C4-O1={float(row.C4_O1_A):.2f} A\"\n")
        handle.write("set dash_color, black\n")
        handle.write("zoom (resn UL1+FE1+HD1+HD2+GU1+AT1+HH1 or resi 1+173+256+335+364+365+366), 7\n")
        handle.write("orient (resn UL1+FE1+HD1+HD2+GU1+AT1+HH1 or resi 1+173+256+335+364+365+366)\n")


def row_summary(sub: pd.DataFrame, name: str) -> dict[str, object]:
    total = len(sub)
    counts = sub["pseudo_RS_by_product_ref"].value_counts()
    return {
        "subset": name,
        "n_frames": total,
        "R_like_n": int(counts.get("R_like_by_product_reference", 0)),
        "R_like_fraction": float(counts.get("R_like_by_product_reference", 0) / total) if total else np.nan,
        "S_like_n": int(counts.get("S_like_by_product_reference", 0)),
        "S_like_fraction": float(counts.get("S_like_by_product_reference", 0) / total) if total else np.nan,
        "median_C4_O1_A": float(sub["C4_O1_A"].median()) if total else np.nan,
        "median_signed_triple": float(sub["signed_triple_O1_C5phenyl_C3chain_about_C4"].median()) if total else np.nan,
    }


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    u = mda.Universe(str(TOP), str(TRAJ))
    atoms = {name: one(u, f"resname UL1 and name {name}") for name in ["O1", "C3", "C4", "C5", "H02"]}

    rows = []
    for ts in u.trajectory:
        c4 = atoms["C4"].position.copy()
        o1 = atoms["O1"].position.copy()
        c5 = atoms["C5"].position.copy()
        c3 = atoms["C3"].position.copy()
        h02 = atoms["H02"].position.copy()
        triple = signed_triple(c4, o1, c5, c3)
        tetra = signed_tetra(o1, c5, c3, h02)
        rows.append(
            {
                "frame": int(ts.frame),
                "time_ps": float(ts.time),
                "time_ns": float(ts.time / 1000.0),
                "C4_O1_A": float(np.linalg.norm(c4 - o1)),
                "signed_triple_O1_C5phenyl_C3chain_about_C4": triple,
                "signed_tetra_O1_C5phenyl_C3chain_using_H02": tetra,
                "pseudo_RS_by_product_ref": assignment(triple),
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(OUTDIR / "md_product_calibrated_pseudo_RS_timeseries.csv", index=False)
    summarize(df).to_csv(OUTDIR / "md_product_calibrated_pseudo_RS_summary.csv", index=False)

    selected = df[df["frame"].isin([22506, 33055, 19517, 4483, 9122, 10080])]
    selected = selected.sort_values("C4_O1_A")
    selected.to_csv(OUTDIR / "selected_frame_pseudo_RS_assignment.csv", index=False)
    write_representative_pdbs_and_pml(u, df)

    print("Wrote", OUTDIR)
    print(selected.to_string(index=False))
    print()
    print(pd.read_csv(OUTDIR / "md_product_calibrated_pseudo_RS_summary.csv").to_string(index=False))


if __name__ == "__main__":
    main()
