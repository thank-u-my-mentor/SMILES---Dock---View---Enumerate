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
    # Compute full pairwise squared distances
    diff = points[:, np.newaxis, :] - points[np.newaxis, :, :]
    dists_sq = (diff ** 2).sum(axis=2)
    dists = np.sqrt(dists_sq)
    # Sort each row; column 0 is self (distance 0)
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
    parser.add_argument(--continue)
