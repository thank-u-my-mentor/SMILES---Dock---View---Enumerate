#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import shutil
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


ROOT = Path("/mnt/e/TJ/260706_mcpb_m2_200ns")
PBC_TOP = ROOT / "md_200ns/pbc/md_200ns_prot_mcpb_center_start.gro"
PBC_XTC = ROOT / "md_200ns/pbc/md_200ns_prot_mcpb_fit.xtc"
ANALYSIS = ROOT / "md_200ns/analysis_reactive_C4O1_20260708"
FACE_DIR = ANALYSIS / "ul1_rotamer_face_analysis"
TIMESERIES = FACE_DIR / "ul1_c3c4_c4c5_face_timeseries.csv"
OUTDIR = FACE_DIR / "same_vs_opposite_near_attack_view"
PRODUCT_DIR = ANALYSIS / "product_RS_reference"


PRODUCTS = [
    ("R_product_user_C@@H", "[H]/N=C(O1)/CC[C@@H]1C2=CC=CC=C2"),
    ("S_product_user_C@H", "[H]/N=C(O1)/CC[C@H]1C2=CC=CC=C2"),
]


def norm(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n == 0.0:
        return v
    return v / n


def signed_height(c3: np.ndarray, c4: np.ndarray, c5: np.ndarray, point: np.ndarray) -> float:
    # Same orientation as tmp_analyze_ul1_rotamer_face.py:
    # normal = cross(C3-C4, C5-C4), origin C4.
    normal = norm(np.cross(c3 - c4, c5 - c4))
    return float(np.dot(point - c4, normal))


def signed_triple(center: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    # Oriented triple product for three ordered substituent vectors around the chiral center.
    return float(np.linalg.det(np.vstack([a - center, b - center, c - center])))


def signed_tetra(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> float:
    # Same ordered substituents, but using H as the fourth reference point.
    return float(np.linalg.det(np.vstack([a - d, b - d, c - d])))


def circular_delta_deg(a: float, b: float) -> float:
    return ((a - b + 180.0) % 360.0) - 180.0


def write_threshold_stats(df: pd.DataFrame) -> None:
    thresholds = [3.50, 3.25, 3.00, 2.90, 2.80, 2.75, 2.70, 2.65]
    rows = []
    for cutoff in thresholds:
        sub = df[df["C4_O1_A"] <= cutoff]
        rows.append(
            {
                "subset": f"C4_O1_le_{cutoff:.2f}A",
                "n_frames": len(sub),
                "same_n": int((sub["O1_H02_face_relation"] == "same").sum()),
                "same_fraction": float((sub["O1_H02_face_relation"] == "same").mean()) if len(sub) else math.nan,
                "opposite_n": int((sub["O1_H02_face_relation"] == "opposite").sum()),
                "opposite_fraction": float((sub["O1_H02_face_relation"] == "opposite").mean()) if len(sub) else math.nan,
                "median_C4_O1_A": float(sub["C4_O1_A"].median()) if len(sub) else math.nan,
            }
        )
    sorted_df = df.sort_values("C4_O1_A")
    for n in [5, 10, 20, 50, 100, 500, 1000]:
        sub = sorted_df.head(n)
        rows.append(
            {
                "subset": f"top_{n}_shortest_C4_O1",
                "n_frames": len(sub),
                "same_n": int((sub["O1_H02_face_relation"] == "same").sum()),
                "same_fraction": float((sub["O1_H02_face_relation"] == "same").mean()),
                "opposite_n": int((sub["O1_H02_face_relation"] == "opposite").sum()),
                "opposite_fraction": float((sub["O1_H02_face_relation"] == "opposite").mean()),
                "median_C4_O1_A": float(sub["C4_O1_A"].median()),
            }
        )
    pd.DataFrame(rows).to_csv(OUTDIR / "near_attack_face_threshold_stats.csv", index=False)


def choose_representatives(df: pd.DataFrame) -> pd.DataFrame:
    near = df[df["C4_O1_A"] <= 3.0].copy()
    rows = []
    for relation in ["opposite", "same"]:
        sub = near[near["O1_H02_face_relation"] == relation]
        if sub.empty:
            continue
        row = sub.sort_values("C4_O1_A").iloc[0].copy()
        row["selection_reason"] = f"shortest_C4_O1_in_{relation}_face_near3A"
        rows.append(row)
    reps = pd.DataFrame(rows)
    reps.to_csv(OUTDIR / "same_vs_opposite_representative_frames.csv", index=False)
    return reps


def write_full_frame_pdbs(reps: pd.DataFrame) -> None:
    u = mda.Universe(str(PBC_TOP), str(PBC_XTC))
    all_atoms = u.select_atoms("all")
    for row in reps.itertuples(index=False):
        u.trajectory[int(row.frame)]
        tag = (
            f"frame{int(row.frame):05d}_time{float(row.time_ns):07.3f}ns_"
            f"C4O1_{float(row.C4_O1_A):.2f}_{row.O1_H02_face_relation}"
        )
        out = OUTDIR / f"{tag}_protein_mcpb.pdb"
        all_atoms.write(str(out))

        active_src = FACE_DIR / "selected_active_site_pdbs"
        active_matches = list(active_src.glob(f"frame{int(row.frame):05d}_*_active_site.pdb"))
        if active_matches:
            shutil.copy2(active_matches[0], OUTDIR / f"{tag}_active_site.pdb")


def write_view_pml(reps: pd.DataFrame) -> None:
    objects = []
    for row in reps.itertuples(index=False):
        tag = (
            f"frame{int(row.frame):05d}_time{float(row.time_ns):07.3f}ns_"
            f"C4O1_{float(row.C4_O1_A):.2f}_{row.O1_H02_face_relation}"
        )
        objects.append((row.O1_H02_face_relation, tag, OUTDIR / f"{tag}_protein_mcpb.pdb", row))

    pml = OUTDIR / "view_same_vs_opposite_near_attack.pml"
    with pml.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("reinitialize\n")
        handle.write("bg_color white\n")
        handle.write("set stick_radius, 0.15\n")
        handle.write("set sphere_scale, 0.22\n")
        handle.write("set dash_radius, 0.05\n")
        handle.write("set label_size, 16\n")
        handle.write("set label_color, black\n")
        handle.write("set cartoon_transparency, 0.72\n")
        for relation, tag, pdb, row in objects:
            obj = f"{relation}_frame{int(row.frame):05d}"
            color = "gray85" if relation == "opposite" else "gray55"
            handle.write(f"load {pdb.as_posix()}, {obj}\n")
            handle.write(f"dss {obj}\n")
            handle.write(f"hide everything, {obj}\n")
            handle.write(f"show cartoon, {obj} and polymer\n")
            handle.write(f"color {color}, {obj} and elem C\n")
            handle.write(f"show sticks, {obj} and (resn UL1+FE1+HD1+HD2+GU1+AT1+HH1 or resi 1+173+256+335+364+365+366)\n")
            handle.write(f"show spheres, {obj} and resn FE1\n")
            handle.write(f"color orange, {obj} and resn FE1\n")
            handle.write(f"color cyan, {obj} and resn UL1 and elem C\n")
            handle.write(f"color gray60, {obj} and resn UL1 and name H02\n")
            handle.write(f"color red, {obj} and resn UL1 and name O1\n")
            handle.write(f"color yellow, {obj} and resn UL1 and name C4\n")
            handle.write(f"distance {relation}_C4_O1, {obj} and resn UL1 and name C4, {obj} and resn UL1 and name O1\n")
            handle.write(f"distance {relation}_O1_H02, {obj} and resn UL1 and name O1, {obj} and resn UL1 and name H02\n")
            handle.write(f"label {obj} and resn UL1 and name C4, \"{relation} frame {int(row.frame)}; C4-O1={float(row.C4_O1_A):.2f} A\"\n")
        handle.write("hide labels, *_O1_H02\n")
        handle.write("set dash_color, black\n")
        handle.write("zoom (resn UL1+FE1+HD1+HD2+GU1+AT1+HH1 or resi 1+173+256+335+364+365+366), 7\n")
        handle.write("orient (resn UL1+FE1+HD1+HD2+GU1+AT1+HH1 or resi 1+173+256+335+364+365+366)\n")
        handle.write("set_view (\\\n")
        handle.write("     0.680000,    0.170000,   -0.713000,\\\n")
        handle.write("    -0.520000,    0.790000,   -0.325000,\\\n")
        handle.write("     0.517000,    0.589000,    0.621000,\\\n")
        handle.write("     0.000000,    0.000000, -110.000000,\\\n")
        handle.write("    43.800000,   37.300000,   44.600000,\\\n")
        handle.write("    80.000000,  140.000000,  -20.000000 )\n")


def embed_product(name: str, smiles: str) -> dict[str, object]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse {smiles}")
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 20260710
    params.useRandomCoords = True
    if AllChem.EmbedMolecule(mol_h, params) != 0:
        raise RuntimeError(f"Embedding failed for {name}")
    if AllChem.MMFFHasAllMoleculeParams(mol_h):
        AllChem.MMFFOptimizeMolecule(mol_h, maxIters=1000)
    else:
        AllChem.UFFOptimizeMolecule(mol_h, maxIters=1000)

    Chem.MolToMolFile(mol_h, str(PRODUCT_DIR / f"{name}.sdf"))
    Chem.MolToPDBFile(mol_h, str(PRODUCT_DIR / f"{name}.pdb"))

    chiral_idx, cip = centers[0]
    center_h = mol_h.GetAtomWithIdx(chiral_idx)
    conf = mol_h.GetConformer()
    pos = lambda idx: np.array(conf.GetAtomPosition(idx), dtype=float)

    o_idx = None
    phenyl_idx = None
    aliphatic_idx = None
    h_idx = None
    for nbr in center_h.GetNeighbors():
        idx = nbr.GetIdx()
        atomic = nbr.GetAtomicNum()
        if atomic == 1:
            h_idx = idx
        elif atomic == 8:
            o_idx = idx
        elif atomic == 6 and nbr.GetIsAromatic():
            phenyl_idx = idx
        elif atomic == 6:
            aliphatic_idx = idx

    missing = [
        label
        for label, idx in [("O", o_idx), ("phenyl_ipso", phenyl_idx), ("aliphatic_C", aliphatic_idx), ("H", h_idx)]
        if idx is None
    ]
    if missing:
        raise RuntimeError(f"{name}: missing substituent mapping {missing}")

    center = pos(chiral_idx)
    h_o = signed_height(pos(aliphatic_idx), center, pos(phenyl_idx), pos(o_idx))
    h_h = signed_height(pos(aliphatic_idx), center, pos(phenyl_idx), pos(h_idx))
    relation = "same" if h_o * h_h > 0 else "opposite"
    triple = signed_triple(center, pos(o_idx), pos(phenyl_idx), pos(aliphatic_idx))
    tetra = signed_tetra(pos(o_idx), pos(phenyl_idx), pos(aliphatic_idx), pos(h_idx))

    return {
        "name": name,
        "input_smiles": smiles,
        "rdkit_isomeric_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
        "rdkit_cip": cip,
        "chiral_center_atom_index_0based": chiral_idx,
        "mapped_O_atom_index_0based": o_idx,
        "mapped_phenyl_ipso_atom_index_0based": phenyl_idx,
        "mapped_aliphatic_C_atom_index_0based": aliphatic_idx,
        "mapped_H_atom_index_0based": h_idx,
        "O_height_to_aliphatic_center_phenyl_plane_A": h_o,
        "H_height_to_aliphatic_center_phenyl_plane_A": h_h,
        "O_H_face_relation": relation,
        "signed_triple_O_phenyl_aliphatic_about_center": triple,
        "signed_tetra_O_phenyl_aliphatic_using_H": tetra,
    }


def write_product_refs() -> None:
    PRODUCT_DIR.mkdir(parents=True, exist_ok=True)
    rows = [embed_product(name, smiles) for name, smiles in PRODUCTS]
    with (PRODUCT_DIR / "product_RS_reference_summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    with (PRODUCT_DIR / "view_product_RS_reference.pml").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("reinitialize\n")
        handle.write("bg_color white\n")
        handle.write("set stick_radius, 0.15\n")
        for row in rows:
            obj = row["name"]
            handle.write(f"load {(PRODUCT_DIR / (obj + '.pdb')).as_posix()}, {obj}\n")
            handle.write(f"hide everything, {obj}\nshow sticks, {obj}\n")
            handle.write(f"color gray80, {obj} and elem C\n")
            handle.write(f"label {obj} and index {int(row['chiral_center_atom_index_0based']) + 1}, \"{row['rdkit_cip']}\"\n")
        handle.write("zoom all, 4\n")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(TIMESERIES)
    write_threshold_stats(df)
    reps = choose_representatives(df)
    write_full_frame_pdbs(reps)
    write_view_pml(reps)
    write_product_refs()

    print(f"Wrote {OUTDIR}")
    print(f"Wrote {PRODUCT_DIR}")
    print("\nRepresentative frames:")
    print(reps[["frame", "time_ns", "C4_O1_A", "O1_H02_face_relation", "tau_C2_C3_C4_C5_deg", "phi_C1_O1_C4_C5_deg"]].to_string(index=False))
    print("\nThreshold stats:")
    print(pd.read_csv(OUTDIR / "near_attack_face_threshold_stats.csv").to_string(index=False))
    print("\nProduct reference:")
    print(pd.read_csv(PRODUCT_DIR / "product_RS_reference_summary.csv").to_string(index=False))


if __name__ == "__main__":
    main()
