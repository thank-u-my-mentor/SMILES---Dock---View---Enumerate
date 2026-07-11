#!/usr/bin/env python3
"""Cluster MD frames by mechanistic active-site geometry features.

This is meant for reactive-geometry classification, not global stability.
It clusters scalar features such as C4-O1 distance, Fe-donor distances,
approach angles, and key-residue contacts.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import calc_angles, calc_dihedrals
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler


ATOM_SELECTIONS = {
    "UL1_C4": "resname UL1 and name C4",
    "UL1_O1": "resname UL1 and name O1",
    "UL1_C1": "resname UL1 and name C1",
    "UL1_C3": "resname UL1 and name C3",
    "UL1_C5": "resname UL1 and name C5",
    "FE": "resname FE1 and name FE",
    "UL1_N1": "resname UL1 and name N1",
    "HD1_NE2": "resname HD1 and name NE2",
    "HD2_NE2": "resname HD2 and name NE2",
    "GU1_OE1": "resname GU1 and name OE1",
    "AT1_O2": "resname AT1 and name O2",
    "HH1_O": "resname HH1 and name O",
}

KEY_RESIDUES = {
    "PHE333": "resid 333 and not name H*",
    "HD2_256": "resid 256 and not name H*",
    "PHE322": "resid 322 and not name H*",
    "GLN255": "resid 255 and not name H*",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--top", required=True)
    parser.add_argument("--traj", required=True)
    parser.add_argument("--distance-csv", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--stride", type=int, default=10, help="Use every Nth frame for the all-frame set")
    parser.add_argument("--short-cutoff", type=float, default=3.0, help="Angstrom cutoff for the reactive subset")
    parser.add_argument("--k-min", type=int, default=2)
    parser.add_argument("--k-max", type=int, default=8)
    return parser.parse_args()


def one_atom(u: mda.Universe, selection: str):
    ag = u.select_atoms(selection)
    if len(ag) != 1:
        raise SystemExit(f"Selection {selection!r} matched {len(ag)} atoms, expected 1")
    return ag[0]


def dist(a, b) -> float:
    return float(np.linalg.norm(a.position - b.position))


def angle_deg(a, b, c) -> float:
    value = calc_angles(a.position[None, :], b.position[None, :], c.position[None, :])[0]
    return float(np.degrees(value))


def dihedral_deg(a, b, c, d) -> float:
    value = calc_dihedrals(a.position[None, :], b.position[None, :], c.position[None, :], d.position[None, :])[0]
    return float(np.degrees(value))


def min_residue_distance(residue_atoms, atoms) -> float:
    coords = residue_atoms.positions
    targets = np.array([atom.position for atom in atoms])
    dmin = float("inf")
    for target in targets:
        dmin = min(dmin, float(np.linalg.norm(coords - target, axis=1).min()))
    return dmin


def build_features(u: mda.Universe, frame_rows: pd.DataFrame) -> pd.DataFrame:
    atoms = {name: one_atom(u, sel) for name, sel in ATOM_SELECTIONS.items()}
    residues = {name: u.select_atoms(sel) for name, sel in KEY_RESIDUES.items()}
    for name, ag in residues.items():
        if len(ag) == 0:
            raise SystemExit(f"Residue selection matched no atoms: {name} {KEY_RESIDUES[name]}")

    rows = []
    for row in frame_rows.itertuples(index=False):
        u.trajectory[int(row.frame)]
        c4 = atoms["UL1_C4"]
        o1 = atoms["UL1_O1"]
        fe = atoms["FE"]
        record = {
            "frame": int(row.frame),
            "time_ps": float(row.time_ps),
            "time_ns": float(row.time_ns),
            "C4_O1_distance_A": float(row.C4_O1_distance_A),
            "Fe_C4_A": dist(fe, c4),
            "Fe_O1_A": dist(fe, o1),
            "Fe_UL1_N1_A": dist(fe, atoms["UL1_N1"]),
            "Fe_HD1_NE2_A": dist(fe, atoms["HD1_NE2"]),
            "Fe_HD2_NE2_A": dist(fe, atoms["HD2_NE2"]),
            "Fe_GU1_OE1_A": dist(fe, atoms["GU1_OE1"]),
            "Fe_AT1_O2_A": dist(fe, atoms["AT1_O2"]),
            "Fe_HH1_O_A": dist(fe, atoms["HH1_O"]),
            "angle_C3_C4_O1_deg": angle_deg(atoms["UL1_C3"], c4, o1),
            "angle_C5_C4_O1_deg": angle_deg(atoms["UL1_C5"], c4, o1),
            "angle_C4_O1_C1_deg": angle_deg(c4, o1, atoms["UL1_C1"]),
            "dihedral_C3_C4_O1_C1_deg": dihedral_deg(atoms["UL1_C3"], c4, o1, atoms["UL1_C1"]),
        }
        for resname, ag in residues.items():
            record[f"{resname}_min_to_C4O1_A"] = min_residue_distance(ag, [c4, o1])
        rows.append(record)
    return pd.DataFrame(rows)


def choose_k_and_cluster(features: pd.DataFrame, k_min: int, k_max: int):
    meta_cols = {"frame", "time_ps", "time_ns"}
    feature_cols = [c for c in features.columns if c not in meta_cols]
    X = features[feature_cols].to_numpy(dtype=float)
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)
    rows = []
    labels_by_k = {}
    centers_by_k = {}
    for k in range(k_min, min(k_max, len(features) - 1) + 1):
        model = KMeans(n_clusters=k, n_init=50, random_state=20260708)
        labels = model.fit_predict(Xs)
        sil = silhouette_score(Xs, labels) if len(set(labels)) > 1 else float("nan")
        counts = pd.Series(labels).value_counts().sort_values(ascending=False)
        rows.append(
            {
                "k": k,
                "silhouette": float(sil),
                "largest_cluster_size": int(counts.iloc[0]),
                "largest_cluster_fraction": float(counts.iloc[0] / len(features)),
                "smallest_cluster_size": int(counts.iloc[-1]),
            }
        )
        labels_by_k[k] = labels
        centers_by_k[k] = model.cluster_centers_
    scan = pd.DataFrame(rows).sort_values(["silhouette", "k"], ascending=[False, True])
    best_k = int(scan.iloc[0]["k"])
    return feature_cols, Xs, scan, best_k, labels_by_k[best_k], centers_by_k[best_k]


def medoid_table(features: pd.DataFrame, feature_cols, Xs: np.ndarray, labels: np.ndarray) -> pd.DataFrame:
    rows = []
    total = len(features)
    for cid in sorted(set(labels), key=lambda x: (-np.sum(labels == x), x)):
        members = np.where(labels == cid)[0]
        center = Xs[members].mean(axis=0)
        d = np.linalg.norm(Xs[members] - center, axis=1)
        medoid_i = int(members[int(np.argmin(d))])
        subset = features.iloc[members]
        medoid = features.iloc[medoid_i]
        row = {
            "cluster_id": int(cid),
            "n_frames": int(len(members)),
            "population_fraction": float(len(members) / total),
            "medoid_sample_index": medoid_i,
            "medoid_original_frame": int(medoid["frame"]),
            "medoid_time_ns": float(medoid["time_ns"]),
            "medoid_C4_O1_A": float(medoid["C4_O1_distance_A"]),
            "min_C4_O1_A": float(subset["C4_O1_distance_A"].min()),
            "median_C4_O1_A": float(subset["C4_O1_distance_A"].median()),
            "max_C4_O1_A": float(subset["C4_O1_distance_A"].max()),
        }
        for c in feature_cols:
            row[f"mean_{c}"] = float(subset[c].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def safe_num(x: float) -> str:
    return f"{x:.3f}".replace(".", "p")


def export_medoids(u: mda.Universe, table: pd.DataFrame, outdir: Path, prefix: str) -> None:
    medoid_dir = outdir / f"{prefix}_medoid_pdbs"
    medoid_dir.mkdir(parents=True, exist_ok=True)
    for rank, row in enumerate(table.itertuples(index=False), start=1):
        u.trajectory[int(row.medoid_original_frame)]
        filename = (
            medoid_dir
            / f"rank{rank:02d}_cluster{int(row.cluster_id):02d}_pop{float(row.population_fraction):.3f}"
            f"_frame{int(row.medoid_original_frame):05d}_time_{safe_num(float(row.medoid_time_ns))}ns"
            f"_C4O1_{safe_num(float(row.medoid_C4_O1_A))}A.pdb"
        )
        u.atoms.write(str(filename))


def export_pml(outdir: Path, prefix: str, table: pd.DataFrame) -> None:
    medoid_dir = outdir / f"{prefix}_medoid_pdbs"
    pml = outdir / f"view_{prefix}_feature_cluster_medoids.pml"
    lines = ["reinitialize", "bg_color white"]
    for rank, _ in enumerate(table.itertuples(index=False), start=1):
        files = sorted(medoid_dir.glob(f"rank{rank:02d}_cluster*.pdb"))
        if files:
            lines.append(f"load {files[0].as_posix().replace('/mnt/e/', 'E:/')}, {prefix}_cluster{rank:02d}")
    lines += [
        "hide everything",
        "select active_site, (resn UL1 or resn FE1 or resn HD1 or resn HD2 or resn GU1 or resn AT1 or resn HH1)",
        "show cartoon, polymer",
        "show sticks, active_site",
        "show spheres, active_site and resn FE1",
        "set stick_radius, 0.15",
        "set sphere_scale, 0.35, active_site and resn FE1",
        "color gray75, polymer",
        "color orange, active_site and resn FE1",
        "color gray90, active_site and resn UL1",
        "distance C4_O1, active_site and resn UL1 and name C4, active_site and resn UL1 and name O1",
        "zoom active_site, 8",
    ]
    pml.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_subset(u: mda.Universe, rows: pd.DataFrame, outdir: Path, prefix: str, k_min: int, k_max: int) -> None:
    features = build_features(u, rows)
    features.to_csv(outdir / f"{prefix}_features.csv", index=False)
    feature_cols, Xs, scan, best_k, labels, _ = choose_k_and_cluster(features, k_min, k_max)
    scan.to_csv(outdir / f"{prefix}_k_scan.csv", index=False)
    assigned = features.copy()
    assigned["cluster_id"] = labels
    assigned.to_csv(outdir / f"{prefix}_frame_assignments_k{best_k}.csv", index=False)
    table = medoid_table(features, feature_cols, Xs, labels)
    table.to_csv(outdir / f"{prefix}_cluster_table_k{best_k}.csv", index=False)
    export_medoids(u, table, outdir, f"{prefix}_k{best_k}")
    export_pml(outdir, f"{prefix}_k{best_k}", table)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    dist = pd.read_csv(args.distance_csv)
    all_rows = dist.iloc[:: args.stride].reset_index(drop=True)
    short_rows = dist[dist["C4_O1_distance_A"] <= args.short_cutoff].reset_index(drop=True)
    u = mda.Universe(args.top, args.traj)

    run_subset(u, all_rows, outdir, "all_stride", args.k_min, args.k_max)
    run_subset(u, short_rows, outdir, f"short_C4O1_le_{str(args.short_cutoff).replace('.', 'p')}A", args.k_min, args.k_max)

    readme = [
        "# Reactive Geometry Feature Clustering",
        "",
        "This clusters mechanistic scalar features rather than whole-protein CA RMSD.",
        "",
        "Features include C4-O1 distance, Fe-donor distances, Fe-C4/Fe-O1 distances, approach angles, a C3-C4-O1-C1 dihedral, and min distances from PHE333, HD2/His256, PHE322, and GLN255 to C4/O1.",
        "",
        f"All-frame set: every {args.stride} frames = {len(all_rows)} frames.",
        f"Reactive short-distance set: C4-O1 <= {args.short_cutoff:.3f} A = {len(short_rows)} frames.",
        "",
        "Medoid filenames use `p` as decimal point, e.g. `time_112p530ns` means 112.530 ns.",
    ]
    (outdir / "README_feature_clustering.md").write_text("\n".join(readme) + "\n", encoding="utf-8")
    print(outdir)


if __name__ == "__main__":
    main()
