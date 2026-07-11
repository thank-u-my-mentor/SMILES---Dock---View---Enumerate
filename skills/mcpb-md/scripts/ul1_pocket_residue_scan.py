#!/usr/bin/env python3
"""UL1-centered pocket residue feature scan for MCPB/MD trajectories.

The script streams an MD trajectory frame by frame, finds all non-solvent
residues within a UL1-centered cutoff, extracts residue-level geometric
features, and then builds unsupervised pocket-state embeddings. R/S labels are
joined only after feature extraction, so they are used as posterior labels
rather than as clustering inputs.
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from MDAnalysis import Universe
from MDAnalysis.lib.distances import capped_distance, distance_array
from sklearn.cluster import HDBSCAN, OPTICS
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

try:
    import umap
except Exception:  # pragma: no cover - optional dependency
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
    "cationic": {"ARG", "LYS", "HIP"},
    "anionic": {"ASP", "GLU", "GU1"},
    "polar": {"SER", "THR", "ASN", "GLN", "CYS", "HIS", "HID", "HIE", "HIP", "HD1", "HD2", "TYR", "TRP"},
    "hydrophobic": {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "PRO", "TRP"},
    "metal": {"FE1"},
    "cofactor": {"AT1", "HH1"},
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


@dataclass(frozen=True)
class ResidueMeta:
    resindex: int
    resid: int
    resname: str
    reskey: str
    class_label: str
    is_aromatic: int
    is_cationic: int
    is_anionic: int
    is_polar: int
    is_hydrophobic: int
    is_metal: int
    is_cofactor: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scan all trajectory frames for UL1-centered pocket residue features."
    )
    parser.add_argument("--engine", required=True, help="Label such as gromacs or amber_pmemd.")
    parser.add_argument("--top", required=True, help="Topology/structure file readable by MDAnalysis.")
    parser.add_argument("--traj", required=True, help="Trajectory file readable by MDAnalysis.")
    parser.add_argument("--rs-csv", required=True, help="Frame-wise C4-O1/R-like/S-like CSV.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--target-resname", default="UL1")
    parser.add_argument("--max-cutoff", type=float, default=6.0)
    parser.add_argument("--contact-cutoffs", type=float, nargs="+", default=[4.0, 5.0, 6.0])
    parser.add_argument("--c4o1-cutoffs", type=float, nargs="+", default=[3.5, 4.0, 4.5, 5.0])
    parser.add_argument("--frame-stride", type=int, default=1, help="1 means every saved frame.")
    parser.add_argument("--max-frames", type=int, default=0, help="0 means all frames.")
    parser.add_argument("--min-occupancy", type=float, default=0.01)
    parser.add_argument("--max-residues-for-matrix", type=int, default=80)
    parser.add_argument("--hdbscan-min-cluster-size", type=int, default=300)
    parser.add_argument("--hdbscan-min-samples", type=int, default=30)
    parser.add_argument("--skip-umap", action="store_true")
    return parser.parse_args()


def is_hydrogen_name(name: str) -> bool:
    stripped = name.strip().upper()
    return stripped.startswith("H")


def safe_box(ts) -> np.ndarray | None:
    dims = getattr(ts, "dimensions", None)
    if dims is None:
        return None
    arr = np.asarray(dims, dtype=float)
    if arr.shape[0] < 6 or np.any(~np.isfinite(arr[:3])) or np.any(arr[:3] <= 0):
        return None
    return arr


def residue_class(resname: str) -> tuple[str, dict[str, int]]:
    flags = {key: int(resname in values) for key, values in RESIDUE_CLASSES.items()}
    if flags["metal"]:
        label = "metal"
    elif flags["cofactor"]:
        label = "cofactor"
    elif flags["aromatic"]:
        label = "aromatic"
    elif flags["cationic"]:
        label = "cationic"
    elif flags["anionic"]:
        label = "anionic"
    elif flags["polar"]:
        label = "polar"
    elif flags["hydrophobic"]:
        label = "hydrophobic"
    else:
        label = "other"
    return label, flags


def make_residue_meta(residue) -> ResidueMeta:
    label, flags = residue_class(residue.resname)
    return ResidueMeta(
        resindex=int(residue.resindex),
        resid=int(residue.resid),
        resname=str(residue.resname),
        reskey=f"{residue.resname}{int(residue.resid)}",
        class_label=label,
        is_aromatic=flags["aromatic"],
        is_cationic=flags["cationic"],
        is_anionic=flags["anionic"],
        is_polar=flags["polar"],
        is_hydrophobic=flags["hydrophobic"],
        is_metal=flags["metal"],
        is_cofactor=flags["cofactor"],
    )


def centroid(positions: np.ndarray) -> np.ndarray:
    if positions.size == 0:
        return np.full(3, np.nan)
    return np.mean(positions, axis=0)


def plane_normal(positions: np.ndarray) -> np.ndarray | None:
    if positions.shape[0] < 3:
        return None
    centered = positions - positions.mean(axis=0)
    try:
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
    except np.linalg.LinAlgError:
        return None
    normal = vh[-1]
    norm = np.linalg.norm(normal)
    if norm == 0 or not np.isfinite(norm):
        return None
    return normal / norm


def normal_angle_deg(n1: np.ndarray | None, n2: np.ndarray | None) -> float:
    if n1 is None or n2 is None:
        return np.nan
    cosv = float(np.clip(abs(np.dot(n1, n2)), 0.0, 1.0))
    return math.degrees(math.acos(cosv))


def point_distance(a: np.ndarray, b: np.ndarray, box: np.ndarray | None) -> float:
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        return np.nan
    return float(distance_array(a.reshape(1, 3), b.reshape(1, 3), box=box)[0, 0])


def min_distance_info(
    residue_positions: np.ndarray,
    residue_names: np.ndarray,
    target_positions: np.ndarray,
    target_names: np.ndarray,
    box: np.ndarray | None,
) -> tuple[float, str, str]:
    if residue_positions.size == 0 or target_positions.size == 0:
        return np.nan, "", ""
    dist = distance_array(residue_positions, target_positions, box=box)
    flat_idx = int(np.argmin(dist))
    i, j = np.unravel_index(flat_idx, dist.shape)
    return float(dist[i, j]), str(residue_names[i]), str(target_names[j])


def load_rs_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "frame" not in df.columns:
        raise ValueError(f"{path} lacks a frame column")
    if "pseudo_RS_by_product_ref" in df.columns:
        raw = df["pseudo_RS_by_product_ref"].astype(str)
    elif "RS_label" not in df.columns:
        raise ValueError(f"{path} lacks pseudo_RS_by_product_ref or RS_label")
    else:
        raw = df["RS_label"].astype(str)
    df["RS_label_raw"] = raw
    df["RS_label"] = raw.map(normalize_rs_label)
    keep = [col for col in ["frame", "time_ps", "time_ns", "C4_O1_A", "RS_label"] if col in df.columns]
    keep.append("RS_label_raw")
    if "signed_triple_O1_C5phenyl_C3chain_about_C4" in df.columns:
        keep.append("signed_triple_O1_C5phenyl_C3chain_about_C4")
    return df[keep].copy()


def normalize_rs_label(value: str) -> str:
    text = str(value).strip()
    upper = text.upper()
    if upper.startswith("R"):
        return "R-like"
    if upper.startswith("S"):
        return "S-like"
    return text


def write_readme(path: Path, lines: Iterable[str]) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_static_groups(u: Universe, target_resname: str):
    target = u.select_atoms(f"resname {target_resname}")
    if len(target) == 0:
        raise ValueError(f"No atoms found for target resname {target_resname}")
    target_heavy = target.select_atoms("not name H*")
    if len(target_heavy) == 0:
        raise ValueError(f"No heavy atoms found for target resname {target_resname}")
    c4 = target.select_atoms("name C4")
    o1 = target.select_atoms("name O1")
    phenyl = target.select_atoms("name C5 C6 C7 C8 C9 C10")
    if len(c4) != 1 or len(o1) != 1 or len(phenyl) < 3:
        raise ValueError(
            f"Expected UL1 C4/O1/phenyl atoms, got C4={len(c4)}, O1={len(o1)}, phenyl={len(phenyl)}"
        )

    pool_atoms = []
    for atom in u.atoms:
        if atom.resname == target_resname:
            continue
        if atom.resname in SOLVENT_RESNAMES:
            continue
        if is_hydrogen_name(atom.name):
            continue
        pool_atoms.append(atom)
    pool = u.atoms[[atom.ix for atom in pool_atoms]]
    if len(pool) == 0:
        raise ValueError("No non-solvent heavy atoms available for pocket scan")

    res_to_pool_indices: dict[int, np.ndarray] = {}
    res_to_side_indices: dict[int, np.ndarray] = {}
    res_to_ring_indices: dict[int, np.ndarray] = {}
    meta: dict[int, ResidueMeta] = {}
    pool_ix_to_resindex = np.array([atom.resindex for atom in pool], dtype=int)

    for residue in u.residues:
        if residue.resname == target_resname or residue.resname in SOLVENT_RESNAMES:
            continue
        heavy_atoms = [atom for atom in residue.atoms if not is_hydrogen_name(atom.name)]
        heavy_pool_ix = np.array([i for i, atom in enumerate(pool) if atom.resindex == residue.resindex], dtype=int)
        if heavy_pool_ix.size == 0 or len(heavy_atoms) == 0:
            continue
        res_to_pool_indices[int(residue.resindex)] = heavy_pool_ix
        side_names = {atom.name for atom in heavy_atoms if atom.name.upper() not in BACKBONE_HEAVY}
        side_pool_ix = np.array(
            [
                i
                for i, atom in enumerate(pool)
                if atom.resindex == residue.resindex and atom.name.upper() in side_names
            ],
            dtype=int,
        )
        if side_pool_ix.size == 0:
            side_pool_ix = heavy_pool_ix
        res_to_side_indices[int(residue.resindex)] = side_pool_ix

        ring_names = AROMATIC_RING_ATOMS.get(residue.resname, set())
        ring_pool_ix = np.array(
            [
                i
                for i, atom in enumerate(pool)
                if atom.resindex == residue.resindex and atom.name.upper() in ring_names
            ],
            dtype=int,
        )
        res_to_ring_indices[int(residue.resindex)] = ring_pool_ix
        meta[int(residue.resindex)] = make_residue_meta(residue)

    return {
        "target": target,
        "target_heavy": target_heavy,
        "c4": c4,
        "o1": o1,
        "phenyl": phenyl,
        "pool": pool,
        "pool_ix_to_resindex": pool_ix_to_resindex,
        "res_to_pool_indices": res_to_pool_indices,
        "res_to_side_indices": res_to_side_indices,
        "res_to_ring_indices": res_to_ring_indices,
        "meta": meta,
    }


def scan_frames(args: argparse.Namespace, u: Universe, groups: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict] = []
    frame_rows: list[dict] = []
    total_frames = len(u.trajectory)
    max_frames = args.max_frames if args.max_frames and args.max_frames > 0 else total_frames

    for analyzed_count, ts in enumerate(u.trajectory[:: args.frame_stride]):
        if analyzed_count >= max_frames:
            break
        if analyzed_count % 1000 == 0:
            print(
                f"[{args.engine}] frame {int(ts.frame)} ({analyzed_count + 1}/{max_frames if max_frames else total_frames})",
                flush=True,
            )

        box = safe_box(ts)
        target_heavy_pos = groups["target_heavy"].positions.astype(np.float64)
        target_heavy_names = np.asarray(groups["target_heavy"].names, dtype=object)
        pool_pos = groups["pool"].positions.astype(np.float64)
        pairs, pair_dist = capped_distance(
            target_heavy_pos,
            pool_pos,
            max_cutoff=float(args.max_cutoff),
            box=box,
            return_distances=True,
        )
        if len(pairs) == 0:
            frame_rows.append(
                {
                    "engine": args.engine,
                    "frame": int(ts.frame),
                    "time_ps_mda": float(getattr(ts, "time", np.nan)),
                    "n_res_4A": 0,
                    "n_res_5A": 0,
                    "n_res_6A": 0,
                }
            )
            continue

        c4_pos = groups["c4"].positions.astype(np.float64)[0]
        o1_pos = groups["o1"].positions.astype(np.float64)[0]
        c4o1_pos = np.vstack([c4_pos, o1_pos])
        c4o1_names = np.asarray(["C4", "O1"], dtype=object)
        phenyl_pos = groups["phenyl"].positions.astype(np.float64)
        phenyl_cent = centroid(phenyl_pos)
        phenyl_normal = plane_normal(phenyl_pos)

        contacted_res = groups["pool_ix_to_resindex"][pairs[:, 1]]
        unique_resindices = np.unique(contacted_res)
        min_by_res = {}
        for resindex in unique_resindices:
            mask = contacted_res == resindex
            min_by_res[int(resindex)] = float(np.min(pair_dist[mask]))

        n_by_cutoff = {
            f"n_res_{str(cut).replace('.', 'p')}A": int(sum(dist <= cut for dist in min_by_res.values()))
            for cut in args.contact_cutoffs
        }
        frame_rows.append(
            {
                "engine": args.engine,
                "frame": int(ts.frame),
                "time_ps_mda": float(getattr(ts, "time", np.nan)),
                **n_by_cutoff,
            }
        )

        for resindex in unique_resindices:
            resindex = int(resindex)
            meta = groups["meta"].get(resindex)
            if meta is None:
                continue
            pool_idx = groups["res_to_pool_indices"][resindex]
            side_idx = groups["res_to_side_indices"][resindex]
            ring_idx = groups["res_to_ring_indices"].get(resindex, np.array([], dtype=int))

            res_pos = pool_pos[pool_idx]
            res_names = np.asarray(groups["pool"].names[pool_idx], dtype=object)
            side_cent = centroid(pool_pos[side_idx])
            ring_pos = pool_pos[ring_idx] if ring_idx.size else np.empty((0, 3))
            ring_normal = plane_normal(ring_pos)

            min_ul1, closest_res_atom, closest_ul1_atom = min_distance_info(
                res_pos, res_names, target_heavy_pos, target_heavy_names, box
            )
            min_c4o1, closest_c4o1_res_atom, closest_c4o1_atom = min_distance_info(
                res_pos, res_names, c4o1_pos, c4o1_names, box
            )

            row = {
                "engine": args.engine,
                "frame": int(ts.frame),
                "time_ps_mda": float(getattr(ts, "time", np.nan)),
                "resindex": meta.resindex,
                "resid": meta.resid,
                "resname": meta.resname,
                "reskey": meta.reskey,
                "class_label": meta.class_label,
                "is_aromatic": meta.is_aromatic,
                "is_cationic": meta.is_cationic,
                "is_anionic": meta.is_anionic,
                "is_polar": meta.is_polar,
                "is_hydrophobic": meta.is_hydrophobic,
                "is_metal": meta.is_metal,
                "is_cofactor": meta.is_cofactor,
                "min_to_UL1_heavy_A": min_ul1,
                "min_to_C4O1_A": min_c4o1,
                "sidechain_centroid_to_C4_A": point_distance(side_cent, c4_pos, box),
                "sidechain_centroid_to_O1_A": point_distance(side_cent, o1_pos, box),
                "sidechain_centroid_to_UL1_phenyl_centroid_A": point_distance(side_cent, phenyl_cent, box),
                "closest_atom_to_UL1": closest_res_atom,
                "closest_UL1_atom": closest_ul1_atom,
                "closest_atom_to_C4O1": closest_c4o1_res_atom,
                "closest_C4O1_atom": closest_c4o1_atom,
                "ring_centroid_to_C4_A": point_distance(centroid(ring_pos), c4_pos, box) if ring_idx.size else np.nan,
                "ring_centroid_to_UL1_phenyl_centroid_A": point_distance(centroid(ring_pos), phenyl_cent, box)
                if ring_idx.size
                else np.nan,
                "ring_plane_to_UL1_phenyl_plane_deg": normal_angle_deg(ring_normal, phenyl_normal),
            }
            for cut in args.contact_cutoffs:
                suffix = str(cut).replace(".", "p")
                row[f"contact_UL1_{suffix}A"] = int(min_ul1 <= cut)
            for cut in args.c4o1_cutoffs:
                suffix = str(cut).replace(".", "p")
                row[f"contact_C4O1_{suffix}A"] = int(min_c4o1 <= cut)
            rows.append(row)

    return pd.DataFrame(rows), pd.DataFrame(frame_rows)


def add_rs_labels(feature_df: pd.DataFrame, frame_df: pd.DataFrame, rs_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    labels = rs_df.rename(columns={"frame": "frame"}).copy()
    feature_df = feature_df.merge(labels, on="frame", how="left")
    frame_df = frame_df.merge(labels, on="frame", how="left")
    return feature_df, frame_df


def summarize_residues(feature_df: pd.DataFrame, total_frames: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    if feature_df.empty:
        return pd.DataFrame(), pd.DataFrame()
    agg = feature_df.groupby(["engine", "resid", "resname", "reskey", "class_label"], dropna=False).agg(
        frames_seen=("frame", "nunique"),
        occupancy_6A=("contact_UL1_6p0A", "mean"),
        occupancy_5A=("contact_UL1_5p0A", "mean"),
        occupancy_4A=("contact_UL1_4p0A", "mean"),
        median_min_to_UL1_A=("min_to_UL1_heavy_A", "median"),
        mean_min_to_UL1_A=("min_to_UL1_heavy_A", "mean"),
        median_min_to_C4O1_A=("min_to_C4O1_A", "median"),
        mean_min_to_C4O1_A=("min_to_C4O1_A", "mean"),
        median_sidechain_centroid_to_C4_A=("sidechain_centroid_to_C4_A", "median"),
        median_sidechain_centroid_to_O1_A=("sidechain_centroid_to_O1_A", "median"),
        median_sidechain_centroid_to_phenyl_A=("sidechain_centroid_to_UL1_phenyl_centroid_A", "median"),
        median_ring_plane_angle_deg=("ring_plane_to_UL1_phenyl_plane_deg", "median"),
        closest_atom_mode=("closest_atom_to_UL1", lambda s: s.mode().iloc[0] if not s.mode().empty else ""),
        closest_UL1_atom_mode=("closest_UL1_atom", lambda s: s.mode().iloc[0] if not s.mode().empty else ""),
    )
    agg = agg.reset_index()
    agg["frames_total"] = int(total_frames)
    agg["seen_fraction"] = agg["frames_seen"] / float(total_frames)
    agg = agg.sort_values(["seen_fraction", "median_min_to_UL1_A"], ascending=[False, True])

    rs_agg = (
        feature_df.groupby(["engine", "RS_label", "resid", "resname", "reskey", "class_label"], dropna=False)
        .agg(
            frames_seen=("frame", "nunique"),
            occupancy_6A=("contact_UL1_6p0A", "mean"),
            occupancy_5A=("contact_UL1_5p0A", "mean"),
            occupancy_4A=("contact_UL1_4p0A", "mean"),
            median_min_to_UL1_A=("min_to_UL1_heavy_A", "median"),
            mean_min_to_UL1_A=("min_to_UL1_heavy_A", "mean"),
            median_min_to_C4O1_A=("min_to_C4O1_A", "median"),
            mean_min_to_C4O1_A=("min_to_C4O1_A", "mean"),
            median_sidechain_centroid_to_C4_A=("sidechain_centroid_to_C4_A", "median"),
            median_sidechain_centroid_to_O1_A=("sidechain_centroid_to_O1_A", "median"),
            median_sidechain_centroid_to_phenyl_A=("sidechain_centroid_to_UL1_phenyl_centroid_A", "median"),
            median_ring_plane_angle_deg=("ring_plane_to_UL1_phenyl_plane_deg", "median"),
        )
        .reset_index()
        .sort_values(["resid", "RS_label"])
    )
    return agg, rs_agg


def build_frame_matrix(
    feature_df: pd.DataFrame,
    frame_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    min_occupancy: float,
    max_residues: int,
    fill_distance: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if feature_df.empty or summary_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    candidates = summary_df[summary_df["seen_fraction"] >= min_occupancy].copy()
    if len(candidates) > max_residues:
        candidates = candidates.head(max_residues)
    candidate_keys = candidates["reskey"].tolist()

    matrix = frame_df[["engine", "frame", "time_ps_mda", "RS_label", "C4_O1_A"]].drop_duplicates("frame").copy()
    numeric_features = [
        "min_to_UL1_heavy_A",
        "min_to_C4O1_A",
        "sidechain_centroid_to_C4_A",
        "sidechain_centroid_to_O1_A",
        "sidechain_centroid_to_UL1_phenyl_centroid_A",
        "ring_plane_to_UL1_phenyl_plane_deg",
        "contact_UL1_4p0A",
        "contact_UL1_5p0A",
        "contact_UL1_6p0A",
    ]

    for reskey in candidate_keys:
        sub = feature_df[feature_df["reskey"] == reskey]
        for feat in numeric_features:
            if feat not in sub.columns:
                continue
            col = f"{reskey}__{feat}"
            tmp = sub[["frame", feat]].drop_duplicates("frame").rename(columns={feat: col})
            matrix = matrix.merge(tmp, on="frame", how="left")

    distance_cols = [c for c in matrix.columns if c.endswith("_A") or c.endswith("_deg")]
    contact_cols = [c for c in matrix.columns if "__contact_" in c]
    for col in distance_cols:
        if col == "C4_O1_A":
            continue
        if col.endswith("_deg"):
            matrix[col] = matrix[col].fillna(90.0)
        else:
            matrix[col] = matrix[col].fillna(fill_distance)
    for col in contact_cols:
        matrix[col] = matrix[col].fillna(0.0)
    return matrix, candidates


def embed_and_cluster(matrix: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    if matrix.empty:
        return matrix
    feat_cols = [c for c in matrix.columns if "__" in c]
    X = matrix[feat_cols].replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy(dtype=float)
    if X.shape[1] == 0 or X.shape[0] < 5:
        return matrix
    Xs = StandardScaler().fit_transform(X)
    pca_components = min(10, Xs.shape[0], Xs.shape[1])
    pca = PCA(n_components=pca_components, random_state=42)
    pcs = pca.fit_transform(Xs)
    matrix = matrix.copy()
    matrix["PC1"] = pcs[:, 0]
    matrix["PC2"] = pcs[:, 1] if pcs.shape[1] > 1 else 0.0
    for i, var in enumerate(pca.explained_variance_ratio_[: min(5, pca_components)], start=1):
        matrix[f"PC{i}_explained_var"] = var

    min_cluster_size = min(args.hdbscan_min_cluster_size, max(5, Xs.shape[0] // 10))
    min_samples = min(args.hdbscan_min_samples, max(2, min_cluster_size // 4))
    try:
        clusterer = HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples)
        labels = clusterer.fit_predict(pcs[:, : min(5, pcs.shape[1])])
        method = "sklearn_HDBSCAN_on_PCA"
    except Exception:
        clusterer = OPTICS(min_samples=min_samples, xi=0.05, min_cluster_size=min_cluster_size)
        labels = clusterer.fit_predict(pcs[:, : min(5, pcs.shape[1])])
        method = "OPTICS_on_PCA_fallback"
    matrix["pocket_cluster"] = labels.astype(int)
    matrix["pocket_cluster_method"] = method

    if not args.skip_umap and umap is not None and Xs.shape[0] >= 10:
        reducer = umap.UMAP(n_neighbors=40, min_dist=0.15, metric="euclidean", random_state=42)
        uv = reducer.fit_transform(Xs)
        matrix["UMAP1"] = uv[:, 0]
        matrix["UMAP2"] = uv[:, 1]
    return matrix


def plot_bar(summary: pd.DataFrame, out: Path, title: str) -> None:
    if summary.empty:
        return
    top = summary.head(30).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.0, 8.0), dpi=220)
    ax.barh(top["reskey"], top["seen_fraction"] * 100.0, color="#4C78A8")
    ax.set_xlabel("Frames within 6 Å of UL1 (%)", fontsize=13)
    ax.set_ylabel("Residue", fontsize=13)
    ax.set_title(title, fontsize=14)
    ax.tick_params(labelsize=10)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_embedding(matrix: pd.DataFrame, out: Path, coord_a: str, coord_b: str, title: str) -> None:
    if matrix.empty or coord_a not in matrix.columns or coord_b not in matrix.columns:
        return
    rs_colors = {"R-like": "#D55E00", "S-like": "#0072B2"}
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 5.2), dpi=220)
    ax = axes[0]
    for label, sub in matrix.groupby("RS_label", dropna=False):
        ax.scatter(
            sub[coord_a],
            sub[coord_b],
            s=4,
            alpha=0.35,
            color=rs_colors.get(str(label), "#666666"),
            label=str(label),
            rasterized=True,
        )
    ax.set_title("Posterior R/S labels", fontsize=13)
    ax.set_xlabel(coord_a, fontsize=12)
    ax.set_ylabel(coord_b, fontsize=12)
    ax.legend(markerscale=4, fontsize=10, frameon=False, loc="best")
    ax.tick_params(labelsize=10)

    ax = axes[1]
    labels = sorted(matrix["pocket_cluster"].dropna().unique()) if "pocket_cluster" in matrix.columns else []
    cmap = plt.get_cmap("tab10")
    for i, label in enumerate(labels):
        sub = matrix[matrix["pocket_cluster"] == label]
        color = "#999999" if int(label) < 0 else cmap(i % 10)
        ax.scatter(sub[coord_a], sub[coord_b], s=4, alpha=0.35, color=color, label=f"State {int(label)}", rasterized=True)
    ax.set_title("Unsupervised pocket states", fontsize=13)
    ax.set_xlabel(coord_a, fontsize=12)
    ax.set_ylabel(coord_b, fontsize=12)
    if labels:
        ax.legend(markerscale=4, fontsize=9, frameon=False, loc="best")
    ax.tick_params(labelsize=10)
    fig.suptitle(title, fontsize=14)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def cluster_rs_summary(matrix: pd.DataFrame) -> pd.DataFrame:
    if matrix.empty or "pocket_cluster" not in matrix.columns:
        return pd.DataFrame()
    rows = []
    for cluster, sub in matrix.groupby("pocket_cluster"):
        total = len(sub)
        r_like = int((sub["RS_label"] == "R-like").sum())
        s_like = int((sub["RS_label"] == "S-like").sum())
        near3 = int((sub["C4_O1_A"] <= 3.0).sum()) if "C4_O1_A" in sub.columns else 0
        rows.append(
            {
                "pocket_cluster": int(cluster),
                "frames": total,
                "R_like_frames": r_like,
                "S_like_frames": s_like,
                "R_like_fraction": r_like / total if total else np.nan,
                "C4_O1_le_3A_frames": near3,
                "C4_O1_le_3A_fraction": near3 / total if total else np.nan,
                "median_C4_O1_A": sub["C4_O1_A"].median() if "C4_O1_A" in sub.columns else np.nan,
            }
        )
    return pd.DataFrame(rows).sort_values("frames", ascending=False)


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.engine}: {args.top} + {args.traj}", flush=True)
    u = Universe(args.top, args.traj)
    groups = build_static_groups(u, args.target_resname)
    rs_df = load_rs_table(Path(args.rs_csv))

    feature_df, frame_df = scan_frames(args, u, groups)
    feature_df, frame_df = add_rs_labels(feature_df, frame_df, rs_df)
    total_analyzed_frames = int(frame_df["frame"].nunique()) if not frame_df.empty else 0

    long_path = outdir / f"{args.engine}_ul1_pocket_frame_residue_features.csv.gz"
    frame_path = outdir / f"{args.engine}_ul1_pocket_frame_counts.csv"
    feature_df.to_csv(long_path, index=False, compression="gzip")
    frame_df.to_csv(frame_path, index=False)

    summary_df, rs_summary_df = summarize_residues(feature_df, total_analyzed_frames)
    summary_path = outdir / f"{args.engine}_ul1_pocket_residue_summary.csv"
    rs_summary_path = outdir / f"{args.engine}_ul1_pocket_residue_summary_by_RS.csv"
    summary_df.to_csv(summary_path, index=False)
    rs_summary_df.to_csv(rs_summary_path, index=False)

    matrix, candidate_df = build_frame_matrix(
        feature_df,
        frame_df,
        summary_df,
        min_occupancy=args.min_occupancy,
        max_residues=args.max_residues_for_matrix,
        fill_distance=args.max_cutoff + 2.0,
    )
    candidate_df.to_csv(outdir / f"{args.engine}_ul1_pocket_matrix_residues.csv", index=False)
    matrix = embed_and_cluster(matrix, args)
    matrix_path = outdir / f"{args.engine}_ul1_pocket_frame_matrix_embeddings.csv.gz"
    matrix.to_csv(matrix_path, index=False, compression="gzip")
    cluster_df = cluster_rs_summary(matrix)
    cluster_df.to_csv(outdir / f"{args.engine}_ul1_pocket_cluster_RS_summary.csv", index=False)

    plot_bar(summary_df, outdir / f"{args.engine}_ul1_pocket_top_residue_occupancy.png", f"{args.engine}: UL1-centered pocket residues")
    plot_embedding(matrix, outdir / f"{args.engine}_ul1_pocket_PCA_states_by_RS.png", "PC1", "PC2", f"{args.engine}: PCA pocket states")
    plot_embedding(matrix, outdir / f"{args.engine}_ul1_pocket_UMAP_states_by_RS.png", "UMAP1", "UMAP2", f"{args.engine}: UMAP pocket states")

    readme_lines = [
        f"# {args.engine} UL1-centered pocket residue scan",
        "",
        f"- Frames analyzed: {total_analyzed_frames}",
        f"- Frame stride: {args.frame_stride} (`1` means every saved frame; no postprocessing stride).",
        f"- UL1 around cutoff: {args.max_cutoff:.1f} A.",
        "- R/S labels were merged after feature extraction. They were not used to define PCA/UMAP/HDBSCAN states.",
        "- `min_to_UL1_heavy_A` and `min_to_C4O1_A` are per-frame closest-pair distances, not trajectory-wide minima.",
        "- Side-chain centroid distances, ring-plane angles, residue class flags, closest atom identities, and contact occupancies are included to avoid relying on a single hand-picked feature.",
        "- Pocket clusters are unsupervised density states from the residue-feature matrix. Use `*_cluster_RS_summary.csv` to check posterior R/S enrichment.",
        "",
        "## Key outputs",
        f"- `{long_path.name}`: long table; one row per contacted residue per frame.",
        f"- `{summary_path.name}`: residue occupancy and median/mean distance summary.",
        f"- `{rs_summary_path.name}`: same summary split by R-like/S-like.",
        f"- `{matrix_path.name}`: fixed-width frame matrix with PCA/UMAP coordinates and pocket cluster labels.",
    ]
    write_readme(outdir / f"README_{args.engine}_ul1_pocket_scan.md", readme_lines)

    print(f"Done: {outdir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
