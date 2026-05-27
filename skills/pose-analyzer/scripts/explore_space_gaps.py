#!/usr/bin/env python
"""Identify under-explored regions in score spaces and recommend frontier seeds.

This script reads the feature matrix from score_space_model.py, maps molecules onto
Pure-SMILES-space and structural-interaction-space, detects low-density gaps,
and ranks "frontier seeds" -- molecules that sit on the boundary between known
and unknown regions. These seeds are ideal starting points for
druglike-pocket-refiner to generate gap-filling candidates.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def safe_float(v: str) -> float:
    try:
        if v in ("", None):
            return math.nan
        return float(v)
    except Exception:
        return math.nan


def grid_density(points: np.ndarray, n_bins: int = 20) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return integer density grid and bin edges."""
    min_x, max_x = float(points[:, 0].min()), float(points[:, 0].max())
    min_y, max_y = float(points[:, 1].min()), float(points[:, 1].max())
    pad_x = (max_x - min_x) * 0.05 + 1e-6
    pad_y = (max_y - min_y) * 0.05 + 1e-6
    x_edges = np.linspace(min_x - pad_x, max_x + pad_x, n_bins + 1)
    y_edges = np.linspace(min_y - pad_y, max_y + pad_y, n_bins + 1)
    density = np.zeros((n_bins, n_bins), dtype=int)
    for x, y in points:
        ix = min(n_bins - 1, max(0, np.digitize(x, x_edges) - 1))
        iy = min(n_bins - 1, max(0, np.digitize(y, y_edges) - 1))
        density[iy, ix] += 1
    return density, x_edges, y_edges


def knn_mean_distance(points: np.ndarray, k: int) -> np.ndarray:
    """Mean distance to k nearest neighbors for each point (O(n^2), ok for n<3000)."""
    n = len(points)
    if n <= 1:
        return np.zeros(n)
    diff = points[:, np.newaxis, :] - points[np.newaxis, :, :]
    dists_sq = (diff ** 2).sum(axis=2)
    dists = np.sqrt(dists_sq)
    sorted_dists = np.sort(dists, axis=1)
    kk = min(k, n - 1)
    return sorted_dists[:, 1 : kk + 1].mean(axis=1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-matrix", type=Path, required=True,
                        help="Path to binding_score_feature_matrix.csv from score_space_model.py")
    parser.add_argument("--outdir", type=Path, default=Path("~/vina_task2/space_exploration"))
    parser.add_argument("--grid-bins", type=int, default=20,
                        help="Grid resolution for density maps")
    parser.add_argument("--density-threshold", type=int, default=1,
                        help="Grid cells with <= this many molecules are considered gaps")
    parser.add_argument("--frontier-k", type=int, default=5,
                        help="K nearest neighbors for frontier scoring")
    parser.add_argument("--top-frontier", type=int, default=20,
                        help="Top frontier seeds to recommend for refinement")
    parser.add_argument("--top-gap-candidates", type=int, default=100,
                        help="Top gap regions to report")
    parser.add_argument("--space-priority", choices=["pure", "structural", "both"], default="structural",
                        help="Which space to prioritize for gap detection. "
                             "'structural' is recommended because it encodes both chemistry and pocket interaction.")
    args = parser.parse_args()

    rows = read_csv(args.feature_matrix.expanduser())

    pure_coords = []
    struct_coords = []
    valid_rows = []
    for row in rows:
        p1 = safe_float(row.get("pure_smiles_pc1"))
        p2 = safe_float(row.get("pure_smiles_pc2"))
        s1 = safe_float(row.get("structural_interaction_pc1"))
        s2 = safe_float(row.get("structural_interaction_pc2"))
        if not any(math.isnan(v) for v in (p1, p2, s1, s2)):
            pure_coords.append([p1, p2])
            struct_coords.append([s1, s2])
            valid_rows.append(row)

    if len(valid_rows) < 10:
        raise ValueError(f"too few valid points: {len(valid_rows)}")

    pure_arr = np.array(pure_coords)
    struct_arr = np.array(struct_coords)

    # Compute densities in both spaces
    pure_density, pure_xe, pure_ye = grid_density(pure_arr, args.grid_bins)
    struct_density, struct_xe, struct_ye = grid_density(struct_arr, args.grid_bins)

    # Frontier scoring: mean distance to k nearest neighbors in each space
    pure_gap = knn_mean_distance(pure_arr, args.frontier_k)
    struct_gap = knn_mean_distance(struct_arr, args.frontier_k)

    # Combine according to priority
    if args.space_priority == "pure":
        combined_gap = pure_gap
    elif args.space_priority == "structural":
        combined_gap = struct_gap
    else:
        combined_gap = (pure_gap + struct_gap) / 2.0

    frontier_scores: list[dict[str, object]] = []
    for i, row in enumerate(valid_rows):
        frontier_scores.append({
            **row,
            "frontier_score": round(float(combined_gap[i]), 6),
            "pure_gap": round(float(pure_gap[i]), 6),
            "struct_gap": round(float(struct_gap[i]), 6),
        })

    frontier_scores.sort(key=lambda r: float(r["frontier_score"]), reverse=True)

    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Write frontier seeds
    frontier_fields = [
        "seq_id", "canonical_smiles", "frontier_score", "pure_gap", "struct_gap",
        "official_binding_score", "predicted_binding_score", "qed",
        "chembl_scaffold_similarity", "affinity_kcal_mol",
    ]
    with (outdir / "frontier_seeds.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=frontier_fields)
        writer.writeheader()
        for row in frontier_scores[:args.top_frontier]:
            writer.writerow({k: row.get(k, "") for k in frontier_fields})

    with (outdir / "frontier_seeds.smi").open("w", encoding="utf-8") as f:
        for row in frontier_scores[:args.top_frontier]:
            smi = row.get("canonical_smiles", "")
            sid = row.get("seq_id", "unknown")
            if smi:
                f.write(f"{smi} {sid}_frontier\n")

    # Identify gap regions (grid cells below density threshold)
    gap_records: list[dict[str, object]] = []
    for space_name, density, xe, ye, arr in [
        ("pure", pure_density, pure_xe, pure_ye, pure_arr),
        ("structural", struct_density, struct_xe, struct_ye, struct_arr),
    ]:
        for iy in range(density.shape[0]):
            for ix in range(density.shape[1]):
                if density[iy, ix] <= args.density_threshold:
                    cx = float((xe[ix] + xe[ix + 1]) / 2)
                    cy = float((ye[iy] + ye[iy + 1]) / 2)
                    # nearest known point
                    dists = np.sqrt(((arr - np.array([cx, cy])) ** 2).sum(axis=1))
                    nearest_idx = int(np.argmin(dists))
                    nearest_dist = float(dists[nearest_idx])
                    gap_records.append({
                        "space": space_name,
                        "grid_x": ix,
                        "grid_y": iy,
                        "center_x": round(cx, 6),
                        "center_y": round(cy, 6),
                        "nearest_seq_id": valid_rows[nearest_idx].get("seq_id", ""),
                        "nearest_smiles": valid_rows[nearest_idx].get("canonical_smiles", ""),
                        "nearest_official_score": valid_rows[nearest_idx].get("official_binding_score", ""),
                        "distance_to_center": round(nearest_dist, 6),
                        "known_density": int(density[iy, ix]),
                    })

    gap_records.sort(key=lambda r: r["distance_to_center"], reverse=True)

    gap_fields = [
        "space", "grid_x", "grid_y", "center_x", "center_y",
        "nearest_seq_id", "nearest_smiles", "nearest_official_score",
        "distance_to_center", "known_density",
    ]
    with (outdir / "gap_regions.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=gap_fields)
        writer.writeheader()
        for row in gap_records[:args.top_gap_candidates]:
            writer.writerow(row)

    # Summary stats
    pure_gaps = sum(1 for r in gap_records if r["space"] == "pure")
    struct_gaps = sum(1 for r in gap_records if r["space"] == "structural")
    print(
        f"points={len(valid_rows)} pure_gaps={pure_gaps} structural_gaps={struct_gaps} "
        f"frontier_seeds={args.top_frontier} wrote={outdir}",
        flush=True,
    )


if __name__ == "__main__":
    main()
