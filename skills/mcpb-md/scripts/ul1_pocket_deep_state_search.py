#!/usr/bin/env python3
"""Search richer UL1-pocket feature spaces for R/S-enriched pocket states.

This script intentionally goes beyond residue min-distance/contact flags. It
streams a trajectory, builds several mutation-compatible pocket representations,
clusters each representation with PCA/UMAP + HDBSCAN parameter sweeps, and only
then overlays product-calibrated R-like/S-like labels as posterior labels.

Feature families:

1. `local_coords`: residue centroids and closest atoms in a UL1-fixed local 3D
   coordinate frame.
2. `ul1_distance_matrix`: distances from each residue to every UL1 heavy atom.
3. `pair_network`: pairwise residue-centroid distances and local pocket moments.
4. `shape_grid`: coarse 3D atom-density grid around UL1 in the same local frame.

The output is an evaluation table, not a claim that any one method is correct.
"""

from __future__ import annotations

import argparse
import itertools
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from MDAnalysis import Universe
from MDAnalysis.lib.distances import distance_array
from sklearn.cluster import HDBSCAN
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

try:
    import umap
except Exception:  # pragma: no cover - optional
    umap = None


SOLVENT_RESNAMES = {
    "WAT",
    "HOH",
    "SOL",
    "TIP3",
    "TIP3P",
    "NA",
    "Na+",
    "NAA",
    "CL",
    "Cl-",
    "CLA",
    "K",
    "K+",
}
BACKBONE_HEAVY = {"N", "CA", "C", "O", "OXT"}
RESIDUE_CLASSES = {
    "aromatic": {"PHE", "TYR", "TRP", "HIS", "HID", "HIE", "HIP", "HD1", "HD2"},
    "cationic": {"ARG", "LYS", "LYN", "HIP"},
    "anionic": {"ASP", "GLU", "GU1"},
    "polar": {"SER", "THR", "ASN", "GLN", "CYS", "HIS", "HID", "HIE", "HIP", "HD1", "HD2", "TYR", "TRP"},
    "hydrophobic": {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "PRO", "TRP"},
    "core": {"FE1", "AT1", "HH1", "GU1"},
}
AROMATIC_RING_ATOMS = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "ND1", "CE1", "NE2", "CD2"},
    "HID": {"CG", "ND1", "CE1", "NE2", "CD2"},
    "HIE": {"CG", "ND1", "CE1", "NE2", "CD2"},
    "HIP": {"CG", "ND1", "CE1", "NE2", "CD2"},
    "HD1": {"CG", "ND1", "CE1", "NE2", "CD2"},
    "HD2": {"CG", "ND1", "CE1", "NE2", "CD2"},
}
UL1_KEY_ATOMS = ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "N1", "O1"]
GRID_CLASSES = ["all", "aromatic", "polar_charged", "hydrophobic", "backbone"]


@dataclass
class ResidueSpec:
    resid: int
    resname: str
    reskey: str
    heavy_indices: np.ndarray
    side_indices: np.ndarray
    backbone_indices: np.ndarray
    ring_indices: np.ndarray
    class_label: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Deep UL1-pocket state search with multiple feature loops.")
    parser.add_argument("--engine", required=True)
    parser.add_argument("--top", required=True)
    parser.add_argument("--traj", required=True)
    parser.add_argument("--rs-csv", required=True)
    parser.add_argument("--residue-summary", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--target-resname", default="UL1")
    parser.add_argument("--candidate-count", type=int, default=35)
    parser.add_argument("--pair-count", type=int, default=24)
    parser.add_argument("--grid-bins", type=int, default=6)
    parser.add_argument("--grid-radius", type=float, default=6.5)
    parser.add_argument("--frame-stride", type=int, default=1)
    parser.add_argument("--max-frames", type=int, default=0)
    parser.add_argument("--save-matrices", action="store_true")
    parser.add_argument("--skip-umap", action="store_true")
    return parser.parse_args()


def is_hydrogen_name(name: str) -> bool:
    return name.strip().upper().startswith("H")


def normalize_label(value: str) -> str:
    text = str(value).strip()
    upper = text.upper()
    if upper.startswith("R"):
        return "R-like"
    if upper.startswith("S"):
        return "S-like"
    return text


def residue_class(resname: str) -> str:
    if resname in RESIDUE_CLASSES["core"]:
        return "core"
    if resname in RESIDUE_CLASSES["aromatic"]:
        return "aromatic"
    if resname in RESIDUE_CLASSES["cationic"]:
        return "cationic"
    if resname in RESIDUE_CLASSES["anionic"]:
        return "anionic"
    if resname in RESIDUE_CLASSES["polar"]:
        return "polar"
    if resname in RESIDUE_CLASSES["hydrophobic"]:
        return "hydrophobic"
    return "other"


def atom_grid_class(resname: str, atom_name: str) -> list[int]:
    classes = [0]
    name = atom_name.upper()
    if resname in RESIDUE_CLASSES["aromatic"]:
        classes.append(1)
    if resname in RESIDUE_CLASSES["polar"] or resname in RESIDUE_CLASSES["cationic"] or resname in RESIDUE_CLASSES["anionic"]:
        classes.append(2)
    if resname in RESIDUE_CLASSES["hydrophobic"]:
        classes.append(3)
    if name in BACKBONE_HEAVY:
        classes.append(4)
    return classes


def safe_box(ts) -> np.ndarray | None:
    dims = getattr(ts, "dimensions", None)
    if dims is None:
        return None
    arr = np.asarray(dims, dtype=float)
    if arr.shape[0] < 6 or np.any(~np.isfinite(arr[:3])) or np.any(arr[:3] <= 0):
        return None
    return arr


def min_image_delta(points: np.ndarray, origin: np.ndarray, box: np.ndarray | None) -> np.ndarray:
    delta = points - origin.reshape(1, 3)
    if box is None:
        return delta
    lengths = box[:3]
    angles = box[3:6]
    if np.allclose(angles, [90.0, 90.0, 90.0], atol=1e-2):
        delta = delta - lengths.reshape(1, 3) * np.round(delta / lengths.reshape(1, 3))
    return delta


def normalize_vec(vec: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vec))
    if not np.isfinite(norm) or norm < 1.0e-8:
        return np.array([1.0, 0.0, 0.0], dtype=float)
    return vec / norm


def local_frame(c4: np.ndarray, o1: np.ndarray, phenyl_centroid: np.ndarray) -> np.ndarray:
    e1 = normalize_vec(o1 - c4)
    v2 = phenyl_centroid - c4
    e3 = np.cross(e1, v2)
    if np.linalg.norm(e3) < 1.0e-8:
        e3 = np.array([0.0, 0.0, 1.0], dtype=float)
    e3 = normalize_vec(e3)
    e2 = normalize_vec(np.cross(e3, e1))
    return np.vstack([e1, e2, e3])


def to_local(points: np.ndarray, origin: np.ndarray, axes: np.ndarray, box: np.ndarray | None) -> np.ndarray:
    delta = min_image_delta(points, origin, box)
    return delta @ axes.T


def centroid(pos: np.ndarray) -> np.ndarray:
    if pos.size == 0:
        return np.full(3, np.nan)
    return np.mean(pos, axis=0)


def plane_normal(pos: np.ndarray) -> np.ndarray:
    if pos.shape[0] < 3:
        return np.zeros(3, dtype=float)
    centered = pos - pos.mean(axis=0)
    try:
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
    except np.linalg.LinAlgError:
        return np.zeros(3, dtype=float)
    n = vh[-1]
    if np.linalg.norm(n) < 1.0e-8:
        return np.zeros(3, dtype=float)
    return n / np.linalg.norm(n)


def load_rs(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "frame" not in df.columns:
        raise ValueError(f"{path} lacks frame column")
    if "pseudo_RS_by_product_ref" in df.columns:
        raw = df["pseudo_RS_by_product_ref"].astype(str)
    elif "RS_label" in df.columns:
        raw = df["RS_label"].astype(str)
    else:
        raise ValueError(f"{path} lacks R/S label column")
    df["RS_label"] = raw.map(normalize_label)
    keep = [col for col in ["frame", "time_ps", "time_ns", "C4_O1_A", "RS_label"] if col in df.columns]
    return df[keep].copy()


def load_candidates(summary_path: Path, engine: str, count: int) -> pd.DataFrame:
    df = pd.read_csv(summary_path)
    if "engine" in df.columns and engine in set(df["engine"].astype(str)):
        df = df[df["engine"].astype(str) == engine].copy()
    required = {"resid", "resname", "reskey"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{summary_path} lacks {sorted(missing)}")
    sort_cols = [c for c in ["seen_fraction", "occupancy_6A", "median_min_to_UL1_A"] if c in df.columns]
    if "seen_fraction" in sort_cols:
        df = df.sort_values(["seen_fraction", "median_min_to_UL1_A"], ascending=[False, True])
    else:
        df = df.sort_values(sort_cols, ascending=False)
    return df.head(count)[["resid", "resname", "reskey"]].drop_duplicates().reset_index(drop=True)


def build_residue_specs(u: Universe, candidates: pd.DataFrame) -> list[ResidueSpec]:
    specs = []
    for row in candidates.itertuples(index=False):
        matches = [res for res in u.residues if int(res.resid) == int(row.resid) and str(res.resname) == str(row.resname)]
        if not matches:
            matches = [res for res in u.residues if int(res.resid) == int(row.resid)]
        if not matches:
            continue
        res = matches[0]
        heavy = np.array([atom.ix for atom in res.atoms if not is_hydrogen_name(atom.name)], dtype=int)
        if heavy.size == 0:
            continue
        side = np.array(
            [
                atom.ix
                for atom in res.atoms
                if not is_hydrogen_name(atom.name) and atom.name.upper() not in BACKBONE_HEAVY
            ],
            dtype=int,
        )
        if side.size == 0:
            side = heavy.copy()
        backbone = np.array(
            [
                atom.ix
                for atom in res.atoms
                if not is_hydrogen_name(atom.name) and atom.name.upper() in BACKBONE_HEAVY
            ],
            dtype=int,
        )
        if backbone.size == 0:
            backbone = heavy.copy()
        ring_names = AROMATIC_RING_ATOMS.get(str(res.resname), set())
        ring = np.array(
            [
                atom.ix
                for atom in res.atoms
                if not is_hydrogen_name(atom.name) and atom.name.upper() in ring_names
            ],
            dtype=int,
        )
        specs.append(
            ResidueSpec(
                resid=int(res.resid),
                resname=str(res.resname),
                reskey=str(row.reskey),
                heavy_indices=heavy,
                side_indices=side,
                backbone_indices=backbone,
                ring_indices=ring,
                class_label=residue_class(str(res.resname)),
            )
        )
    return specs


def build_pool(u: Universe, target_resname: str):
    indices = []
    class_channels = []
    names = []
    reskeys = []
    for atom in u.atoms:
        if atom.resname == target_resname or atom.resname in SOLVENT_RESNAMES or is_hydrogen_name(atom.name):
            continue
        indices.append(atom.ix)
        class_channels.append(atom_grid_class(str(atom.resname), str(atom.name)))
        names.append(str(atom.name))
        reskeys.append(f"{atom.resname}{int(atom.resid)}")
    return np.array(indices, dtype=int), class_channels, names, reskeys


def feature_names(specs: list[ResidueSpec], pair_count: int, grid_bins: int) -> dict[str, list[str]]:
    local_cols = []
    for spec in specs:
        for prefix in ["sc", "bb", "closest"]:
            local_cols.extend([f"{spec.reskey}__{prefix}_{axis}" for axis in ["x", "y", "z", "r"]])
        local_cols.extend([f"{spec.reskey}__ring_centroid_{axis}" for axis in ["x", "y", "z", "r"]])
        local_cols.extend([f"{spec.reskey}__ring_normal_{axis}" for axis in ["x", "y", "z"]])

    dist_cols = []
    for spec in specs:
        for atom in UL1_KEY_ATOMS:
            dist_cols.append(f"{spec.reskey}__sc_to_UL1_{atom}_A")
        for atom in UL1_KEY_ATOMS:
            dist_cols.append(f"{spec.reskey}__min_heavy_to_UL1_{atom}_A")
        for atom in UL1_KEY_ATOMS:
            dist_cols.append(f"{spec.reskey}__mean_heavy_to_UL1_{atom}_A")

    pair_specs = specs[:pair_count]
    pair_cols = []
    for a, b in itertools.combinations(pair_specs, 2):
        pair_cols.append(f"{a.reskey}__to__{b.reskey}__sc_centroid_A")
    for spec in specs:
        for target in ["C4", "O1", "phenyl"]:
            pair_cols.append(f"{spec.reskey}__sc_to_{target}_A")
    for cls in GRID_CLASSES:
        for prefix in ["count", "mean_x", "mean_y", "mean_z", "std_x", "std_y", "std_z"]:
            pair_cols.append(f"pocket_{cls}__{prefix}")

    grid_cols = []
    for cls in GRID_CLASSES:
        for ix in range(grid_bins):
            for iy in range(grid_bins):
                for iz in range(grid_bins):
                    grid_cols.append(f"grid_{cls}_{ix}_{iy}_{iz}")

    return {
        "local_coords": local_cols,
        "ul1_distance_matrix": dist_cols,
        "pair_network": pair_cols,
        "shape_grid": grid_cols,
    }


def allocate_matrices(n_frames: int, names: dict[str, list[str]]) -> dict[str, np.ndarray]:
    return {name: np.zeros((n_frames, len(cols)), dtype=np.float32) for name, cols in names.items()}


def dist_matrix(a: np.ndarray, b: np.ndarray, box: np.ndarray | None) -> np.ndarray:
    if a.size == 0 or b.size == 0:
        return np.empty((0, 0), dtype=float)
    return distance_array(a, b, box=box)


def fill_features_for_frame(
    u: Universe,
    specs: list[ResidueSpec],
    pool_indices: np.ndarray,
    pool_channels: list[list[int]],
    target_atoms,
    matrices: dict[str, np.ndarray],
    row_i: int,
    args: argparse.Namespace,
) -> dict[str, float]:
    ts = u.trajectory.ts
    box = safe_box(ts)
    target = target_atoms["target"]
    c4 = target_atoms["c4"].positions.astype(float)[0]
    o1 = target_atoms["o1"].positions.astype(float)[0]
    phenyl_pos = target_atoms["phenyl"].positions.astype(float)
    phenyl_cent = centroid(phenyl_pos)
    axes = local_frame(c4, o1, phenyl_cent)
    ul1_key_pos = np.vstack([target.select_atoms(f"name {atom}").positions.astype(float)[0] for atom in UL1_KEY_ATOMS])

    local_values = []
    dist_values = []
    side_centroids = []
    side_valid = []

    c4_o1_pos = np.vstack([c4, o1])
    for spec in specs:
        heavy_pos = u.atoms[spec.heavy_indices].positions.astype(float)
        side_pos = u.atoms[spec.side_indices].positions.astype(float)
        bb_pos = u.atoms[spec.backbone_indices].positions.astype(float)
        sc_cent = centroid(side_pos)
        bb_cent = centroid(bb_pos)
        side_centroids.append(sc_cent)
        side_valid.append(np.all(np.isfinite(sc_cent)))

        heavy_to_c4o1 = dist_matrix(heavy_pos, c4_o1_pos, box)
        if heavy_to_c4o1.size:
            closest_idx = int(np.unravel_index(np.argmin(heavy_to_c4o1), heavy_to_c4o1.shape)[0])
            closest_pos = heavy_pos[closest_idx]
            closest_dist = float(np.min(heavy_to_c4o1))
        else:
            closest_pos = np.full(3, np.nan)
            closest_dist = args.grid_radius + 3.0

        for point in [sc_cent, bb_cent, closest_pos]:
            loc = to_local(point.reshape(1, 3), c4, axes, box)[0] if np.all(np.isfinite(point)) else np.zeros(3)
            local_values.extend([loc[0], loc[1], loc[2], float(np.linalg.norm(loc))])
        local_values[-1] = closest_dist

        if spec.ring_indices.size:
            ring_pos = u.atoms[spec.ring_indices].positions.astype(float)
            ring_cent = centroid(ring_pos)
            ring_loc = to_local(ring_cent.reshape(1, 3), c4, axes, box)[0]
            n_global = plane_normal(ring_pos)
            n_local = n_global @ axes.T if np.linalg.norm(n_global) > 0 else np.zeros(3)
            local_values.extend([ring_loc[0], ring_loc[1], ring_loc[2], float(np.linalg.norm(ring_loc))])
            local_values.extend([n_local[0], n_local[1], n_local[2]])
        else:
            local_values.extend([0.0] * 7)

        sc_to_ul1 = dist_matrix(sc_cent.reshape(1, 3), ul1_key_pos, box)[0]
        heavy_to_ul1 = dist_matrix(heavy_pos, ul1_key_pos, box)
        dist_values.extend(sc_to_ul1.tolist())
        dist_values.extend(np.min(heavy_to_ul1, axis=0).tolist())
        dist_values.extend(np.mean(heavy_to_ul1, axis=0).tolist())

    side_centroids_arr = np.asarray(side_centroids, dtype=float)
    pair_values = []
    pair_n = min(args.pair_count, len(specs))
    if pair_n > 1:
        pair_d = dist_matrix(side_centroids_arr[:pair_n], side_centroids_arr[:pair_n], box)
        for i, j in itertools.combinations(range(pair_n), 2):
            pair_values.append(float(pair_d[i, j]))
    for sc_cent in side_centroids_arr:
        d = dist_matrix(sc_cent.reshape(1, 3), np.vstack([c4, o1, phenyl_cent]), box)[0]
        pair_values.extend(d.tolist())

    # Pocket shape moments and grid occupancy in the UL1-fixed local frame.
    pool_pos = u.atoms[pool_indices].positions.astype(float)
    pool_local = to_local(pool_pos, c4, axes, box)
    radii = np.linalg.norm(pool_local, axis=1)
    inside = radii <= args.grid_radius
    selected_local = pool_local[inside]
    selected_channels = [pool_channels[i] for i in np.flatnonzero(inside)]

    class_points: dict[int, list[np.ndarray]] = {i: [] for i in range(len(GRID_CLASSES))}
    for loc, channels in zip(selected_local, selected_channels):
        for channel in channels:
            class_points[channel].append(loc)

    for channel in range(len(GRID_CLASSES)):
        pts = np.asarray(class_points[channel], dtype=float)
        if pts.size == 0:
            pair_values.extend([0.0] * 7)
        else:
            pair_values.append(float(len(pts)))
            pair_values.extend(np.mean(pts, axis=0).tolist())
            pair_values.extend(np.std(pts, axis=0).tolist())

    grid = np.zeros((len(GRID_CLASSES), args.grid_bins, args.grid_bins, args.grid_bins), dtype=np.float32)
    span = 2.0 * args.grid_radius
    for loc, channels in zip(selected_local, selected_channels):
        shifted = (loc + args.grid_radius) / span
        if np.any(shifted < 0.0) or np.any(shifted >= 1.0):
            continue
        ix, iy, iz = np.floor(shifted * args.grid_bins).astype(int)
        ix = int(np.clip(ix, 0, args.grid_bins - 1))
        iy = int(np.clip(iy, 0, args.grid_bins - 1))
        iz = int(np.clip(iz, 0, args.grid_bins - 1))
        for channel in channels:
            grid[channel, ix, iy, iz] += 1.0
    if selected_local.shape[0] > 0:
        grid /= float(selected_local.shape[0])

    matrices["local_coords"][row_i, :] = np.asarray(local_values, dtype=np.float32)
    matrices["ul1_distance_matrix"][row_i, :] = np.asarray(dist_values, dtype=np.float32)
    matrices["pair_network"][row_i, :] = np.asarray(pair_values, dtype=np.float32)
    matrices["shape_grid"][row_i, :] = grid.reshape(-1)

    return {
        "mda_frame": int(ts.frame),
        "mda_time_ps": float(getattr(ts, "time", np.nan)),
        "pocket_atom_count_6p5A": int(selected_local.shape[0]),
    }


def cluster_score(labels: np.ndarray, y: np.ndarray) -> dict[str, float]:
    labels = np.asarray(labels)
    y = np.asarray(y).astype(int)
    valid = labels != -1
    n_total = len(labels)
    n_valid = int(valid.sum())
    global_r = float(y.mean()) if n_total else np.nan
    unique = sorted(set(labels[valid].tolist()))
    if n_valid == 0 or len(unique) == 0:
        return {
            "n_clusters": 0,
            "noise_fraction": 1.0,
            "coverage": 0.0,
            "weighted_abs_R_enrichment": 0.0,
            "max_abs_R_enrichment": 0.0,
            "purity": np.nan,
            "NMI": 0.0,
            "ARI": 0.0,
            "score": 0.0,
        }
    weighted_abs = 0.0
    max_abs = 0.0
    majority = 0
    for label in unique:
        mask = labels == label
        n = int(mask.sum())
        r_count = int(y[mask].sum())
        s_count = n - r_count
        r_frac = r_count / n
        delta = abs(r_frac - global_r)
        weighted_abs += (n / n_valid) * delta
        max_abs = max(max_abs, delta)
        majority += max(r_count, s_count)
    purity = majority / n_valid
    nmi = normalized_mutual_info_score(y[valid], labels[valid]) if len(unique) > 1 else 0.0
    ari = adjusted_rand_score(y[valid], labels[valid]) if len(unique) > 1 else 0.0
    coverage = n_valid / n_total
    # Score favors broad coverage and state-specific R/S enrichment, but still records all raw metrics.
    score = weighted_abs * coverage * math.log1p(len(unique))
    return {
        "n_clusters": int(len(unique)),
        "noise_fraction": 1.0 - coverage,
        "coverage": coverage,
        "weighted_abs_R_enrichment": weighted_abs,
        "max_abs_R_enrichment": max_abs,
        "purity": purity,
        "NMI": nmi,
        "ARI": ari,
        "score": score,
    }


def cluster_details(labels: np.ndarray, y: np.ndarray, c4o1: np.ndarray) -> pd.DataFrame:
    rows = []
    labels = np.asarray(labels)
    y = np.asarray(y).astype(int)
    for label in sorted(set(labels.tolist())):
        mask = labels == label
        n = int(mask.sum())
        if n == 0:
            continue
        r_count = int(y[mask].sum())
        rows.append(
            {
                "state": int(label),
                "frames": n,
                "R_like_frames": r_count,
                "S_like_frames": n - r_count,
                "R_like_fraction": r_count / n,
                "median_C4_O1_A": float(np.nanmedian(c4o1[mask])),
                "C4_O1_le_3A_fraction": float(np.mean(c4o1[mask] <= 3.0)),
            }
        )
    return pd.DataFrame(rows).sort_values("frames", ascending=False)


def block_logistic_auc(X: np.ndarray, y: np.ndarray, n_blocks: int = 5) -> float:
    if len(np.unique(y)) < 2:
        return np.nan
    blocks = np.array_split(np.arange(len(y)), n_blocks)
    aucs = []
    for block in blocks:
        train = np.ones(len(y), dtype=bool)
        train[block] = False
        test = ~train
        if len(np.unique(y[train])) < 2 or len(np.unique(y[test])) < 2:
            continue
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X[train])
        X_test = scaler.transform(X[test])
        clf = LogisticRegression(max_iter=1000, class_weight="balanced", solver="liblinear")
        clf.fit(X_train, y[train])
        prob = clf.predict_proba(X_test)[:, 1]
        aucs.append(roc_auc_score(y[test], prob))
    return float(np.mean(aucs)) if aucs else np.nan


def evaluate_method(
    engine: str,
    method: str,
    X: np.ndarray,
    meta: pd.DataFrame,
    outdir: Path,
    skip_umap: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    y = (meta["RS_label"].astype(str) == "R-like").to_numpy(dtype=int)
    c4o1 = meta["C4_O1_A"].to_numpy(dtype=float)
    X = np.nan_to_num(X.astype(np.float32), nan=0.0, posinf=0.0, neginf=0.0)
    # Remove zero-variance columns before scaling.
    variances = np.var(X, axis=0)
    keep = variances > 1.0e-8
    X = X[:, keep]
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)
    n_pca = min(20, Xs.shape[0], Xs.shape[1])
    pca = PCA(n_components=n_pca, random_state=42)
    pcs = pca.fit_transform(Xs)

    representations = {
        "PCA2": pcs[:, :2],
        "PCA10": pcs[:, : min(10, pcs.shape[1])],
        "PCA20": pcs[:, : min(20, pcs.shape[1])],
    }
    umap_xy = None
    if not skip_umap and umap is not None and Xs.shape[0] >= 50:
        reducer = umap.UMAP(n_neighbors=60, min_dist=0.05, metric="euclidean", random_state=42)
        umap_xy = reducer.fit_transform(pcs[:, : min(20, pcs.shape[1])])
        representations["UMAP2"] = umap_xy

    rows = []
    best = None
    min_cluster_sizes = [50, 100, 300, 800, 1500]
    min_samples_list = [10, 30, 80]
    for rep_name, Z in representations.items():
        for min_cluster_size in min_cluster_sizes:
            if min_cluster_size >= len(y):
                continue
            for min_samples in min_samples_list:
                try:
                    labels = HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples).fit_predict(Z)
                except Exception:
                    continue
                score = cluster_score(labels, y)
                row = {
                    "engine": engine,
                    "method": method,
                    "representation": rep_name,
                    "min_cluster_size": min_cluster_size,
                    "min_samples": min_samples,
                    **score,
                }
                rows.append(row)
                if best is None or row["score"] > best["row"]["score"]:
                    best = {"row": row, "labels": labels, "representation": rep_name}

    score_df = pd.DataFrame(rows)
    if best is None:
        return score_df, pd.DataFrame(), pd.DataFrame()

    best_labels = best["labels"]
    details = cluster_details(best_labels, y, c4o1)
    for key, value in best["row"].items():
        details[key] = value

    embed = meta.copy()
    embed["PC1"] = pcs[:, 0]
    embed["PC2"] = pcs[:, 1] if pcs.shape[1] > 1 else 0.0
    if umap_xy is not None:
        embed["UMAP1"] = umap_xy[:, 0]
        embed["UMAP2"] = umap_xy[:, 1]
    embed["best_state"] = best_labels.astype(int)
    embed["method"] = method
    embed["best_representation"] = best["representation"]

    auc = block_logistic_auc(pcs[:, : min(20, pcs.shape[1])], y, n_blocks=5)
    score_df["blocked_logistic_auc_on_PCA20"] = auc
    details["blocked_logistic_auc_on_PCA20"] = auc

    return score_df, details, embed


def plot_scores(scores: pd.DataFrame, out: Path) -> None:
    if scores.empty:
        return
    best = scores.sort_values("score", ascending=False).groupby("method").head(1).copy()
    best = best.sort_values("score", ascending=True)
    fig, ax = plt.subplots(figsize=(8.0, 4.8), dpi=220)
    ax.barh(best["method"], best["score"], color="#4C78A8")
    ax.set_xlabel("Best HDBSCAN posterior R/S enrichment score", fontsize=12)
    ax.set_ylabel("Feature loop", fontsize=12)
    ax.set_title("Which pocket representation separates R/S better?", fontsize=14)
    ax.tick_params(labelsize=10)
    for i, row in enumerate(best.itertuples(index=False)):
        ax.text(
            row.score,
            i,
            f"  purity={row.purity:.2f}, NMI={row.NMI:.3f}",
            va="center",
            fontsize=9,
        )
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_embedding(embed: pd.DataFrame, out: Path, coord_prefix: str, title: str) -> None:
    c1 = f"{coord_prefix}1"
    c2 = f"{coord_prefix}2"
    if embed.empty or c1 not in embed.columns or c2 not in embed.columns:
        return
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 5.0), dpi=220)
    colors = {"R-like": "#D55E00", "S-like": "#0072B2"}
    ax = axes[0]
    for label, sub in embed.groupby("RS_label"):
        ax.scatter(sub[c1], sub[c2], s=4, alpha=0.35, color=colors.get(label, "#666"), label=label, rasterized=True)
    ax.set_title("Posterior R/S", fontsize=12)
    ax.set_xlabel(c1)
    ax.set_ylabel(c2)
    ax.legend(markerscale=4, frameon=False, fontsize=9)
    ax = axes[1]
    for state, sub in embed.groupby("best_state"):
        color = "#999999" if int(state) < 0 else None
        ax.scatter(sub[c1], sub[c2], s=4, alpha=0.35, color=color, label=f"State {int(state)}", rasterized=True)
    ax.set_title("Best HDBSCAN states", fontsize=12)
    ax.set_xlabel(c1)
    ax.set_ylabel(c2)
    ax.legend(markerscale=4, frameon=False, fontsize=8, loc="best")
    fig.suptitle(title, fontsize=14)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.engine}: {args.top} + {args.traj}", flush=True)
    u = Universe(args.top, args.traj)
    rs = load_rs(Path(args.rs_csv))
    candidates = load_candidates(Path(args.residue_summary), args.engine, args.candidate_count)
    specs = build_residue_specs(u, candidates)
    if len(specs) < 3:
        raise ValueError("Too few candidate residues were found in the topology")
    candidates.to_csv(outdir / f"{args.engine}_deep_state_candidate_residues.csv", index=False)

    target = u.select_atoms(f"resname {args.target_resname}")
    if len(target) == 0:
        raise ValueError(f"No target residue named {args.target_resname}")
    target_atoms = {
        "target": target,
        "c4": target.select_atoms("name C4"),
        "o1": target.select_atoms("name O1"),
        "phenyl": target.select_atoms("name C5 C6 C7 C8 C9 C10"),
    }
    if len(target_atoms["c4"]) != 1 or len(target_atoms["o1"]) != 1 or len(target_atoms["phenyl"]) < 3:
        raise ValueError("UL1 atom selection failed for C4/O1/phenyl")

    pool_indices, pool_channels, _, _ = build_pool(u, args.target_resname)
    total = len(u.trajectory[:: args.frame_stride])
    n_frames = min(args.max_frames, total) if args.max_frames and args.max_frames > 0 else total
    names = feature_names(specs, min(args.pair_count, len(specs)), args.grid_bins)
    matrices = allocate_matrices(n_frames, names)
    meta_rows = []

    for row_i, ts in enumerate(u.trajectory[:: args.frame_stride]):
        if row_i >= n_frames:
            break
        if row_i % 1000 == 0:
            print(f"[{args.engine}] feature frame {int(ts.frame)} ({row_i + 1}/{n_frames})", flush=True)
        meta_rows.append(
            fill_features_for_frame(
                u,
                specs,
                pool_indices,
                pool_channels,
                target_atoms,
                matrices,
                row_i,
                args,
            )
        )

    meta = pd.DataFrame(meta_rows)
    meta = meta.merge(rs, left_on="mda_frame", right_on="frame", how="left")
    meta.insert(0, "engine", args.engine)
    meta.to_csv(outdir / f"{args.engine}_deep_state_frame_metadata.csv", index=False)

    if args.save_matrices:
        for method, X in matrices.items():
            np.savez_compressed(outdir / f"{args.engine}_{method}_matrix.npz", X=X, columns=np.array(names[method]))

    all_scores = []
    all_details = []
    best_embeds = []
    for method, X in matrices.items():
        print(f"[{args.engine}] evaluating {method}: {X.shape}", flush=True)
        scores, details, embed = evaluate_method(args.engine, method, X, meta, outdir, args.skip_umap)
        scores.to_csv(outdir / f"{args.engine}_{method}_hdbscan_sweep_scores.csv", index=False)
        details.to_csv(outdir / f"{args.engine}_{method}_best_state_RS_details.csv", index=False)
        if not embed.empty:
            embed.to_csv(outdir / f"{args.engine}_{method}_best_state_embeddings.csv.gz", index=False, compression="gzip")
            plot_embedding(embed, outdir / f"{args.engine}_{method}_best_PCA_embedding.png", "PC", f"{args.engine} {method}")
            if "UMAP1" in embed.columns:
                plot_embedding(embed, outdir / f"{args.engine}_{method}_best_UMAP_embedding.png", "UMAP", f"{args.engine} {method}")
        all_scores.append(scores)
        all_details.append(details)
        if not embed.empty:
            best_embeds.append(embed)

    score_df = pd.concat(all_scores, ignore_index=True) if all_scores else pd.DataFrame()
    detail_df = pd.concat(all_details, ignore_index=True) if all_details else pd.DataFrame()
    score_df.to_csv(outdir / f"{args.engine}_all_deep_feature_hdbscan_scores.csv", index=False)
    detail_df.to_csv(outdir / f"{args.engine}_all_deep_feature_best_state_details.csv", index=False)
    plot_scores(score_df, outdir / f"{args.engine}_deep_feature_loop_score_comparison.png")

    readme = [
        f"# {args.engine} deep UL1-pocket state search",
        "",
        "This run evaluates four feature loops:",
        "",
        "1. `local_coords`: residue centroids and closest atoms in a UL1-fixed local 3D frame.",
        "2. `ul1_distance_matrix`: residue-to-every-UL1-heavy-atom distance fingerprints.",
        "3. `pair_network`: pairwise residue-centroid distances and pocket atom-cloud moments.",
        "4. `shape_grid`: coarse 3D atom-density grid around UL1.",
        "",
        "R/S labels are merged only after feature extraction and are used for posterior enrichment scoring.",
        "The best-state score is not a mechanism proof; it ranks which unsupervised state definition most enriches R-like vs S-like frames.",
        "Blocked logistic AUC on PCA20 is included only as a diagnostic for whether the feature space contains R/S information under time-block validation.",
        "",
        "Main files:",
        f"- `{args.engine}_all_deep_feature_hdbscan_scores.csv`",
        f"- `{args.engine}_all_deep_feature_best_state_details.csv`",
        f"- `{args.engine}_deep_feature_loop_score_comparison.png`",
    ]
    (outdir / f"README_{args.engine}_deep_state_search.md").write_text("\n".join(readme) + "\n", encoding="utf-8")

    print(f"Done: {outdir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
