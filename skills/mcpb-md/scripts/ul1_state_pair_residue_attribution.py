#!/usr/bin/env python3
"""Attribute a pair of UL1-pocket states to residue-level geometry changes.

This script is the interpretability step after `ul1_pocket_deep_state_search.py`.
It reuses an already chosen state assignment, recomputes mutation-compatible
`local_coords` features only for the selected state-pair frames, then asks:

1. Which residues move or reorient between the two states?
2. Which residues are most important for separating the two state labels under
   a simple time-blocked ablation test?

The ablation is exploratory, not causal. It is designed for presentation-grade
mechanistic triage before considering heavier ML/GNN models.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from MDAnalysis import Universe
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from ul1_pocket_deep_state_search import (
    allocate_matrices,
    build_pool,
    build_residue_specs,
    feature_names,
    fill_features_for_frame,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Residue attribution for a UL1-pocket state pair.")
    parser.add_argument("--engine", default="gromacs_full")
    parser.add_argument("--top", required=True)
    parser.add_argument("--traj", required=True)
    parser.add_argument("--state-embedding", required=True)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--state-a", type=int, default=8, help="Positive/R-enriched state label.")
    parser.add_argument("--state-b", type=int, default=1, help="Negative/S-enriched state label.")
    parser.add_argument("--state-a-name", default="R_enriched_state8")
    parser.add_argument("--state-b-name", default="S_enriched_state1")
    parser.add_argument("--target-resname", default="UL1")
    parser.add_argument("--pair-count", type=int, default=24)
    parser.add_argument("--grid-bins", type=int, default=6)
    parser.add_argument("--grid-radius", type=float, default=6.5)
    parser.add_argument("--c4o1-min", type=float, default=float("nan"))
    parser.add_argument("--c4o1-max", type=float, default=float("nan"))
    parser.add_argument("--n-blocks", type=int, default=5)
    parser.add_argument("--top-n", type=int, default=20)
    return parser.parse_args()


def clean_label(text: str) -> str:
    return str(text).replace(" ", "_").replace("/", "_")


def state_filter(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    required = {"mda_frame", "best_state", "RS_label", "C4_O1_A"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"State embedding lacks required columns: {sorted(missing)}")
    keep = df[df["best_state"].astype(int).isin([args.state_a, args.state_b])].copy()
    if math.isfinite(args.c4o1_min):
        keep = keep[keep["C4_O1_A"].astype(float) >= args.c4o1_min].copy()
    if math.isfinite(args.c4o1_max):
        keep = keep[keep["C4_O1_A"].astype(float) <= args.c4o1_max].copy()
    keep["contrast_label"] = np.where(
        keep["best_state"].astype(int) == args.state_a,
        args.state_a_name,
        args.state_b_name,
    )
    keep = keep.sort_values("mda_frame").reset_index(drop=True)
    if keep.empty:
        raise ValueError("No frames remain after state/C4-O1 filtering")
    return keep


def pooled_sd(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    va = float(np.nanvar(a, ddof=1))
    vb = float(np.nanvar(b, ddof=1))
    n1 = len(a)
    n2 = len(b)
    denom = n1 + n2 - 2
    if denom <= 0:
        return float("nan")
    value = math.sqrt(((n1 - 1) * va + (n2 - 1) * vb) / denom)
    return value if value > 1.0e-12 else float("nan")


def cohen_d(a: np.ndarray, b: np.ndarray) -> float:
    psd = pooled_sd(a, b)
    if not np.isfinite(psd):
        return 0.0
    return float((np.nanmean(a) - np.nanmean(b)) / psd)


def norm_delta(mean_a: np.ndarray, mean_b: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(mean_a, dtype=float) - np.asarray(mean_b, dtype=float)))


def vector_effect(df: pd.DataFrame, cols: list[str], mask_a: np.ndarray, mask_b: np.ndarray) -> tuple[float, float, float, float]:
    if not all(col in df.columns for col in cols):
        return (0.0, 0.0, 0.0, 0.0)
    a = df.loc[mask_a, cols].to_numpy(dtype=float)
    b = df.loc[mask_b, cols].to_numpy(dtype=float)
    mean_a = np.nanmean(a, axis=0)
    mean_b = np.nanmean(b, axis=0)
    deltas = mean_a - mean_b
    effects = np.array([cohen_d(a[:, i], b[:, i]) for i in range(len(cols))], dtype=float)
    return float(deltas[0]), float(deltas[1]), float(deltas[2]), float(np.linalg.norm(effects))


def angle_between_vectors(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na < 1.0e-8 or nb < 1.0e-8:
        return float("nan")
    # Plane normals are sign-ambiguous, so use abs(dot).
    dot = abs(float(np.dot(a / na, b / nb)))
    dot = max(-1.0, min(1.0, dot))
    return float(np.degrees(np.arccos(dot)))


def blocked_auc(X: np.ndarray, y: np.ndarray, frames: np.ndarray, n_blocks: int) -> float:
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=int)
    frames = np.asarray(frames, dtype=int)
    if len(np.unique(y)) < 2 or X.shape[1] == 0:
        return float("nan")
    order = np.argsort(frames)
    blocks = np.array_split(order, max(2, n_blocks))
    aucs = []
    for test_idx in blocks:
        train_mask = np.ones(len(y), dtype=bool)
        train_mask[test_idx] = False
        test_mask = ~train_mask
        if len(np.unique(y[train_mask])) < 2 or len(np.unique(y[test_mask])) < 2:
            continue
        Xt = np.nan_to_num(X[train_mask], nan=0.0, posinf=0.0, neginf=0.0)
        Xv = np.nan_to_num(X[test_mask], nan=0.0, posinf=0.0, neginf=0.0)
        var = np.var(Xt, axis=0)
        keep = var > 1.0e-8
        if not np.any(keep):
            continue
        scaler = StandardScaler()
        Xt = scaler.fit_transform(Xt[:, keep])
        Xv = scaler.transform(Xv[:, keep])
        clf = LogisticRegression(max_iter=1000, class_weight="balanced", solver="liblinear")
        clf.fit(Xt, y[train_mask])
        prob = clf.predict_proba(Xv)[:, 1]
        aucs.append(roc_auc_score(y[test_mask], prob))
    return float(np.mean(aucs)) if aucs else float("nan")


def build_local_feature_table(
    u: Universe,
    selected: pd.DataFrame,
    candidates: pd.DataFrame,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    specs = build_residue_specs(u, candidates)
    if len(specs) < 3:
        raise ValueError("Too few candidate residues were found in the topology")

    target = u.select_atoms(f"resname {args.target_resname}")
    target_atoms = {
        "target": target,
        "c4": target.select_atoms("name C4"),
        "o1": target.select_atoms("name O1"),
        "phenyl": target.select_atoms("name C5 C6 C7 C8 C9 C10"),
    }
    if len(target_atoms["c4"]) != 1 or len(target_atoms["o1"]) != 1 or len(target_atoms["phenyl"]) < 3:
        raise ValueError("UL1 atom selection failed for C4/O1/phenyl")

    pool_indices, pool_channels, _, _ = build_pool(u, args.target_resname)
    names = feature_names(specs, min(args.pair_count, len(specs)), args.grid_bins)
    matrices = allocate_matrices(len(selected), names)
    helper_args = SimpleNamespace(
        grid_radius=args.grid_radius,
        grid_bins=args.grid_bins,
        pair_count=min(args.pair_count, len(specs)),
    )

    wanted = set(selected["mda_frame"].astype(int).tolist())
    row_for_frame = {int(frame): i for i, frame in enumerate(selected["mda_frame"].astype(int).tolist())}
    meta_rows = [None] * len(selected)
    for ts in u.trajectory:
        frame = int(ts.frame)
        if frame not in wanted:
            continue
        row_i = row_for_frame[frame]
        if row_i % 1000 == 0:
            print(f"[state attribution] feature frame {frame} ({row_i + 1}/{len(selected)})", flush=True)
        meta_rows[row_i] = fill_features_for_frame(
            u,
            specs,
            pool_indices,
            pool_channels,
            target_atoms,
            matrices,
            row_i,
            helper_args,
        )
    missing = [frame for frame, row in zip(selected["mda_frame"].astype(int), meta_rows) if row is None]
    if missing:
        raise ValueError(f"Trajectory did not provide {len(missing)} requested frames, e.g. {missing[:5]}")

    local_cols = names["local_coords"]
    local_df = pd.DataFrame(matrices["local_coords"], columns=local_cols)
    selected_meta = selected.reset_index(drop=True).copy()
    recalculated_meta = pd.DataFrame(meta_rows).add_prefix("recalc_")
    feature_df = pd.concat([selected_meta, recalculated_meta, local_df], axis=1)
    return feature_df, matrices["local_coords"], local_cols


def residue_attribution(feature_df: pd.DataFrame, candidates: pd.DataFrame, cols: list[str], args: argparse.Namespace) -> pd.DataFrame:
    mask_a = feature_df["contrast_label"].astype(str).eq(args.state_a_name).to_numpy()
    mask_b = feature_df["contrast_label"].astype(str).eq(args.state_b_name).to_numpy()
    rows = []
    for row in candidates.itertuples(index=False):
        reskey = str(row.reskey)
        prefix = f"{reskey}__"
        res_cols = [col for col in cols if col.startswith(prefix)]
        if not res_cols:
            continue
        record = {
            "resid": int(row.resid),
            "resname": str(row.resname),
            "reskey": reskey,
            "n_features": len(res_cols),
        }
        for label in [args.state_a_name, args.state_b_name]:
            state_mask = feature_df["contrast_label"].astype(str).eq(label)
            record[f"{label}_frames"] = int(state_mask.sum())
            record[f"{label}_R_like_fraction"] = float((feature_df.loc[state_mask, "RS_label"].astype(str) == "R-like").mean())
            record[f"{label}_median_C4_O1_A"] = float(feature_df.loc[state_mask, "C4_O1_A"].median())

        for prefix_name in ["sc", "bb", "closest", "ring_centroid"]:
            xyz = [f"{reskey}__{prefix_name}_{axis}" for axis in ["x", "y", "z"]]
            dx, dy, dz, eff = vector_effect(feature_df, xyz, mask_a, mask_b)
            record[f"{prefix_name}_dx_A"] = dx
            record[f"{prefix_name}_dy_A"] = dy
            record[f"{prefix_name}_dz_A"] = dz
            record[f"{prefix_name}_xyz_delta_A"] = math.sqrt(dx * dx + dy * dy + dz * dz)
            record[f"{prefix_name}_xyz_effect_norm"] = eff

        for scalar in ["sc_r", "bb_r", "closest_r", "ring_centroid_r"]:
            col = f"{reskey}__{scalar}"
            if col in feature_df.columns:
                a = feature_df.loc[mask_a, col].to_numpy(dtype=float)
                b = feature_df.loc[mask_b, col].to_numpy(dtype=float)
                record[f"{scalar}_mean_delta_A"] = float(np.nanmean(a) - np.nanmean(b))
                record[f"{scalar}_cohen_d"] = cohen_d(a, b)

        ncols = [f"{reskey}__ring_normal_{axis}" for axis in ["x", "y", "z"]]
        if all(col in feature_df.columns for col in ncols):
            na = feature_df.loc[mask_a, ncols].to_numpy(dtype=float)
            nb = feature_df.loc[mask_b, ncols].to_numpy(dtype=float)
            record["ring_normal_angle_deg"] = angle_between_vectors(np.nanmean(na, axis=0), np.nanmean(nb, axis=0))

        record["structural_score"] = (
            record.get("sc_xyz_effect_norm", 0.0)
            + 0.50 * record.get("bb_xyz_effect_norm", 0.0)
            + 0.75 * abs(record.get("closest_r_cohen_d", 0.0))
            + 0.75 * record.get("ring_centroid_xyz_effect_norm", 0.0)
        )
        rows.append(record)
    return pd.DataFrame(rows).sort_values("structural_score", ascending=False)


def ablation(feature_matrix: np.ndarray, cols: list[str], candidates: pd.DataFrame, feature_df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    y = (feature_df["contrast_label"].astype(str) == args.state_a_name).to_numpy(dtype=int)
    frames = feature_df["mda_frame"].to_numpy(dtype=int)
    full_auc = blocked_auc(feature_matrix, y, frames, args.n_blocks)
    rows = []
    for row in candidates.itertuples(index=False):
        prefix = f"{row.reskey}__"
        remove = np.array([col.startswith(prefix) for col in cols], dtype=bool)
        if not np.any(remove):
            continue
        reduced_auc = blocked_auc(feature_matrix[:, ~remove], y, frames, args.n_blocks)
        rows.append(
            {
                "resid": int(row.resid),
                "resname": str(row.resname),
                "reskey": str(row.reskey),
                "removed_features": int(remove.sum()),
                "full_blocked_auc": full_auc,
                "without_residue_blocked_auc": reduced_auc,
                "ablation_auc_drop": float(full_auc - reduced_auc) if np.isfinite(full_auc) and np.isfinite(reduced_auc) else float("nan"),
            }
        )
    return pd.DataFrame(rows).sort_values("ablation_auc_drop", ascending=False)


def setup_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 13,
            "axes.labelsize": 15,
            "axes.titlesize": 16,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "legend.fontsize": 12,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def plot_bar(df: pd.DataFrame, value_col: str, title: str, xlabel: str, out: Path, top_n: int, color: str) -> None:
    data = df.head(top_n).iloc[::-1].copy()
    fig, ax = plt.subplots(figsize=(7.2, 7.2), dpi=240)
    labels = data["reskey"].astype(str).tolist()
    ax.barh(labels, data[value_col].astype(float), color=color, edgecolor="black", linewidth=0.4)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_heatmap(df: pd.DataFrame, out: Path, top_n: int) -> None:
    metrics = [
        "sc_xyz_delta_A",
        "sc_xyz_effect_norm",
        "closest_r_mean_delta_A",
        "closest_r_cohen_d",
        "ring_centroid_xyz_delta_A",
        "ring_centroid_xyz_effect_norm",
        "ablation_auc_drop",
    ]
    available = [col for col in metrics if col in df.columns]
    data = df.head(top_n).set_index("reskey")[available].copy()
    values = data.to_numpy(dtype=float)
    scaled = values.copy()
    for j in range(scaled.shape[1]):
        max_abs = np.nanmax(np.abs(scaled[:, j]))
        if not np.isfinite(max_abs) or max_abs == 0:
            max_abs = 1.0
        scaled[:, j] = scaled[:, j] / max_abs
    fig, ax = plt.subplots(figsize=(8.2, 7.2), dpi=240)
    im = ax.imshow(scaled, aspect="auto", cmap="viridis", vmin=-1.0, vmax=1.0)
    ax.set_xticks(range(len(available)))
    ax.set_xticklabels(available, rotation=45, ha="right")
    ax.set_yticks(range(len(data.index)))
    ax.set_yticklabels(data.index)
    ax.set_title("State-pair residue contrast: normalized columns")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("Column-scaled signed magnitude")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            val = values[i, j]
            if np.isfinite(val):
                color = "white" if abs(scaled[i, j]) > 0.55 else "black"
                ax.text(j, i, f"{val:.2g}", ha="center", va="center", fontsize=8.5, color=color)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def write_readme(out: Path, args: argparse.Namespace, selected: pd.DataFrame, attr: pd.DataFrame, abl: pd.DataFrame) -> None:
    top_attr = attr.head(8)[["reskey", "structural_score", "sc_xyz_delta_A", "closest_r_mean_delta_A"]].to_string(index=False)
    top_abl = abl.head(8)[["reskey", "ablation_auc_drop", "without_residue_blocked_auc"]].to_string(index=False)
    text = f"""# UL1 State-Pair Residue Attribution

Engine: `{args.engine}`

State A: `{args.state_a}` / `{args.state_a_name}`

State B: `{args.state_b}` / `{args.state_b_name}`

Frames used: {len(selected)}

C4-O1 filter: min={args.c4o1_min}, max={args.c4o1_max}

## Meaning

This is a contrast between two already discovered unsupervised pocket states.
It is not a direct causal model. R/S labels are posterior labels from the
product-calibrated scalar triple product.

`structural_score` combines side-chain centroid effect, backbone centroid
effect, closest-heavy-atom distance effect, and aromatic ring-centroid effect.
`ablation_auc_drop` asks how much a simple time-blocked linear classifier loses
when all local-coordinate columns for one residue are removed.

## Top Structural Contrasts

```text
{top_attr}
```

## Top Ablation Drops

```text
{top_abl}
```

## Main Files

- `state_pair_residue_attribution.csv`
- `state_pair_residue_ablation.csv`
- `state_pair_selected_frame_metadata.csv`
- `state_pair_feature_values.csv.gz`
- `state_pair_top_residue_structural_score.png`
- `state_pair_ablation_auc_drop.png`
- `state_pair_feature_delta_heatmap.png`
"""
    out.write_text(text, encoding="utf-8")


def main() -> int:
    args = parse_args()
    setup_matplotlib()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    embedding = pd.read_csv(args.state_embedding)
    selected = state_filter(embedding, args)
    candidates = pd.read_csv(args.candidates)
    required = {"resid", "resname", "reskey"}
    missing = required - set(candidates.columns)
    if missing:
        raise ValueError(f"Candidate table lacks {sorted(missing)}")

    print(f"Loading trajectory: {args.top} + {args.traj}", flush=True)
    print(f"Selected {len(selected)} frames from states {args.state_a}/{args.state_b}", flush=True)
    u = Universe(args.top, args.traj)
    feature_df, X_local, local_cols = build_local_feature_table(u, selected, candidates, args)

    attr = residue_attribution(feature_df, candidates, local_cols, args)
    abl = ablation(X_local, local_cols, candidates, feature_df, args)
    merged = attr.merge(
        abl[["reskey", "full_blocked_auc", "without_residue_blocked_auc", "ablation_auc_drop"]],
        on="reskey",
        how="left",
    )
    merged["combined_interpretability_score"] = merged["structural_score"] + 10.0 * merged["ablation_auc_drop"].fillna(0.0)
    merged = merged.sort_values("combined_interpretability_score", ascending=False)

    selected.to_csv(outdir / "state_pair_selected_frame_metadata.csv", index=False)
    feature_df.to_csv(outdir / "state_pair_feature_values.csv.gz", index=False, compression="gzip")
    attr.to_csv(outdir / "state_pair_residue_attribution.csv", index=False)
    abl.to_csv(outdir / "state_pair_residue_ablation.csv", index=False)
    merged.to_csv(outdir / "state_pair_residue_attribution_with_ablation.csv", index=False)

    plot_bar(
        merged.sort_values("structural_score", ascending=False),
        "structural_score",
        "Residue Geometry Contrast",
        "Structural contrast score",
        outdir / "state_pair_top_residue_structural_score.png",
        args.top_n,
        "#4C78A8",
    )
    plot_bar(
        merged.sort_values("ablation_auc_drop", ascending=False),
        "ablation_auc_drop",
        "Residue Ablation For State Separation",
        "Drop in blocked AUC after removing residue",
        outdir / "state_pair_ablation_auc_drop.png",
        args.top_n,
        "#F58518",
    )
    plot_heatmap(merged.sort_values("combined_interpretability_score", ascending=False), outdir / "state_pair_feature_delta_heatmap.png", args.top_n)
    write_readme(outdir / "README_state_pair_residue_attribution.md", args, selected, attr, abl)

    print(f"Wrote {outdir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
