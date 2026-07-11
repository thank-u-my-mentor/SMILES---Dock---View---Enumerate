#!/usr/bin/env python3
"""Combine GROMACS and Amber UL1-pocket scan outputs.

This script consumes the outputs from `ul1_pocket_residue_scan.py` and builds
engine-comparison tables/plots. It keeps the unsupervised pocket states separate
from posterior R/S labels.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ENGINE_DISPLAY = {
    "gromacs_full": "GROMACS",
    "amber_pmemd_full": "Amber/pmemd",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize full-frame UL1 pocket residue scans.")
    parser.add_argument("--gromacs-dir", required=True)
    parser.add_argument("--amber-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--near-cutoff", type=float, default=3.0)
    parser.add_argument("--top-n", type=int, default=35)
    return parser.parse_args()


def find_one(folder: Path, suffix: str) -> Path:
    hits = sorted(folder.glob(f"*{suffix}"))
    if len(hits) != 1:
        raise FileNotFoundError(f"Expected exactly one *{suffix} in {folder}, found {len(hits)}")
    return hits[0]


def load_engine(folder: Path) -> dict[str, pd.DataFrame]:
    tables = {
        "summary": pd.read_csv(find_one(folder, "_ul1_pocket_residue_summary.csv")),
        "by_rs": pd.read_csv(find_one(folder, "_ul1_pocket_residue_summary_by_RS.csv")),
        "frames": pd.read_csv(find_one(folder, "_ul1_pocket_frame_counts.csv")),
        "features": pd.read_csv(find_one(folder, "_ul1_pocket_frame_residue_features.csv.gz")),
        "clusters": pd.read_csv(find_one(folder, "_ul1_pocket_cluster_RS_summary.csv")),
    }
    engine = str(tables["summary"]["engine"].dropna().iloc[0])
    for table in tables.values():
        if "engine" not in table.columns:
            table.insert(0, "engine", engine)
    return tables


def display_engine(value: str) -> str:
    return ENGINE_DISPLAY.get(value, value)


def add_display_engine(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["engine_display"] = df["engine"].map(display_engine)
    return df


def near_summary(features: pd.DataFrame, frames: pd.DataFrame, near_cutoff: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    near_frames = frames[frames["C4_O1_A"] <= near_cutoff].copy()
    near_keys = near_frames[["engine", "frame"]].drop_duplicates().copy()
    sub = features.merge(near_keys, on=["engine", "frame"], how="inner").copy()

    total_by_engine = near_frames.groupby("engine")["frame"].nunique().to_dict()
    total_by_engine_rs = near_frames.groupby(["engine", "RS_label"])["frame"].nunique().to_dict()

    group_cols = ["engine", "resid", "resname", "reskey", "class_label"]
    near_all = (
        sub.groupby(group_cols, dropna=False)
        .agg(
            near_frames_seen=("frame", "nunique"),
            median_min_to_UL1_A=("min_to_UL1_heavy_A", "median"),
            median_min_to_C4O1_A=("min_to_C4O1_A", "median"),
            median_sidechain_centroid_to_C4_A=("sidechain_centroid_to_C4_A", "median"),
            median_sidechain_centroid_to_O1_A=("sidechain_centroid_to_O1_A", "median"),
            median_sidechain_centroid_to_phenyl_A=("sidechain_centroid_to_UL1_phenyl_centroid_A", "median"),
            median_ring_plane_angle_deg=("ring_plane_to_UL1_phenyl_plane_deg", "median"),
        )
        .reset_index()
    )
    near_all["near_total_frames"] = near_all["engine"].map(total_by_engine)
    near_all["near_seen_fraction"] = near_all["near_frames_seen"] / near_all["near_total_frames"]

    group_rs_cols = ["engine", "RS_label", "resid", "resname", "reskey", "class_label"]
    near_rs = (
        sub.groupby(group_rs_cols, dropna=False)
        .agg(
            near_frames_seen=("frame", "nunique"),
            median_min_to_UL1_A=("min_to_UL1_heavy_A", "median"),
            median_min_to_C4O1_A=("min_to_C4O1_A", "median"),
            median_sidechain_centroid_to_C4_A=("sidechain_centroid_to_C4_A", "median"),
            median_sidechain_centroid_to_O1_A=("sidechain_centroid_to_O1_A", "median"),
            median_sidechain_centroid_to_phenyl_A=("sidechain_centroid_to_UL1_phenyl_centroid_A", "median"),
            median_ring_plane_angle_deg=("ring_plane_to_UL1_phenyl_plane_deg", "median"),
        )
        .reset_index()
    )
    near_rs["near_RS_total_frames"] = near_rs.apply(
        lambda r: total_by_engine_rs.get((r["engine"], r["RS_label"]), np.nan), axis=1
    )
    near_rs["near_RS_seen_fraction"] = near_rs["near_frames_seen"] / near_rs["near_RS_total_frames"]
    return near_all, near_rs


def delta_table(by_rs: pd.DataFrame, prefix: str, metrics: list[str]) -> pd.DataFrame:
    rows = []
    key_cols = ["engine", "resid", "resname", "reskey", "class_label"]
    for metric in metrics:
        piv = by_rs.pivot_table(index=key_cols, columns="RS_label", values=metric, aggfunc="first")
        piv = piv.reset_index()
        if "R-like" not in piv.columns or "S-like" not in piv.columns:
            continue
        piv[f"{prefix}_R_like"] = piv["R-like"]
        piv[f"{prefix}_S_like"] = piv["S-like"]
        piv[f"{prefix}_RminusS"] = piv["R-like"] - piv["S-like"]
        piv["metric"] = metric
        rows.append(piv[key_cols + ["metric", f"{prefix}_R_like", f"{prefix}_S_like", f"{prefix}_RminusS"]])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def make_presence_plot(all_summary: pd.DataFrame, near_all: pd.DataFrame, out: Path, top_n: int) -> None:
    all_summary = all_summary.copy()
    near_all = near_all.copy()
    all_summary["context"] = "all frames"
    all_summary["value"] = all_summary["seen_fraction"]
    near_all["context"] = "C4-O1 <= 3.0 Å"
    near_all["value"] = near_all["near_seen_fraction"]
    cols = ["engine", "reskey", "class_label", "context", "value"]
    plot_df = pd.concat([all_summary[cols], near_all[cols]], ignore_index=True)
    rank = plot_df.groupby("reskey")["value"].max().sort_values(ascending=False).head(top_n).index
    plot_df = plot_df[plot_df["reskey"].isin(rank)].copy()
    plot_df["column"] = plot_df["engine"].map(display_engine) + "\n" + plot_df["context"]
    table = plot_df.pivot_table(index="reskey", columns="column", values="value", aggfunc="max").fillna(0.0)
    table = table.loc[rank]

    fig, ax = plt.subplots(figsize=(8.5, 8.5), dpi=240)
    im = ax.imshow(table.values * 100.0, aspect="auto", cmap="viridis", vmin=0, vmax=100)
    ax.set_xticks(range(table.shape[1]))
    ax.set_xticklabels(table.columns, rotation=35, ha="right", fontsize=10)
    ax.set_yticks(range(table.shape[0]))
    ax.set_yticklabels(table.index, fontsize=9)
    ax.set_title("UL1-centered residue presence", fontsize=15)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("Frames where residue is within 6 Å of UL1 (%)", fontsize=11)
    cbar.ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def make_delta_heatmap(delta: pd.DataFrame, metric: str, out: Path, top_reskeys: list[str], title: str) -> None:
    sub = delta[delta["metric"] == metric].copy()
    if sub.empty:
        return
    sub["column"] = sub["engine"].map(display_engine) + "\n" + sub["context"]
    table = sub.pivot_table(index="reskey", columns="column", values="RminusS", aggfunc="first")
    table = table.reindex(top_reskeys).dropna(how="all")
    if table.empty:
        return
    vmax = np.nanpercentile(np.abs(table.values), 95)
    vmax = float(vmax if np.isfinite(vmax) and vmax > 0 else 1.0)

    fig, ax = plt.subplots(figsize=(8.0, 8.5), dpi=240)
    im = ax.imshow(table.values, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(table.shape[1]))
    ax.set_xticklabels(table.columns, rotation=35, ha="right", fontsize=10)
    ax.set_yticks(range(table.shape[0]))
    ax.set_yticklabels(table.index, fontsize=9)
    ax.set_title(title, fontsize=15)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("R-like median minus S-like median (Å or degree)", fontsize=10)
    cbar.ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def make_cluster_plot(clusters: pd.DataFrame, out: Path) -> None:
    if clusters.empty:
        return
    clusters = clusters.copy()
    clusters["engine_display"] = clusters["engine"].map(display_engine)
    clusters["label"] = clusters["engine_display"] + " state " + clusters["pocket_cluster"].astype(str)
    clusters = clusters.sort_values(["engine_display", "frames"], ascending=[True, False])
    x = np.arange(len(clusters))

    fig, ax1 = plt.subplots(figsize=(9.0, 5.2), dpi=240)
    ax1.bar(x - 0.18, clusters["R_like_fraction"] * 100.0, width=0.36, color="#4C78A8", label="R-like")
    ax1.bar(
        x + 0.18,
        clusters["C4_O1_le_3A_fraction"] * 100.0,
        width=0.36,
        color="#F58518",
        label="C4-O1 <= 3 Å",
    )
    ax1.set_ylabel("Fraction within pocket state (%)", fontsize=12)
    ax1.set_xticks(x)
    ax1.set_xticklabels(clusters["label"], rotation=35, ha="right", fontsize=9)
    ax1.tick_params(labelsize=10)
    ax1.legend(frameon=False, fontsize=10, loc="upper right")
    ax1.set_title("Posterior R/S and near-attack enrichment by unsupervised pocket state", fontsize=13)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gmx = load_engine(Path(args.gromacs_dir))
    amber = load_engine(Path(args.amber_dir))

    all_summary = pd.concat([gmx["summary"], amber["summary"]], ignore_index=True)
    all_by_rs = pd.concat([gmx["by_rs"], amber["by_rs"]], ignore_index=True)
    frames = pd.concat([gmx["frames"], amber["frames"]], ignore_index=True)
    clusters = pd.concat([gmx["clusters"], amber["clusters"]], ignore_index=True)
    features = pd.concat([gmx["features"], amber["features"]], ignore_index=True)

    all_summary = add_display_engine(all_summary)
    all_by_rs = add_display_engine(all_by_rs)
    frames = add_display_engine(frames)
    clusters = add_display_engine(clusters)
    features = add_display_engine(features)

    near_all, near_rs = near_summary(features, frames, args.near_cutoff)
    near_all = add_display_engine(near_all)
    near_rs = add_display_engine(near_rs)

    all_summary.to_csv(outdir / "combined_ul1_pocket_residue_summary_all_frames.csv", index=False)
    all_by_rs.to_csv(outdir / "combined_ul1_pocket_residue_summary_by_RS_all_frames.csv", index=False)
    near_all.to_csv(outdir / "combined_ul1_pocket_residue_summary_near3A.csv", index=False)
    near_rs.to_csv(outdir / "combined_ul1_pocket_residue_summary_by_RS_near3A.csv", index=False)
    frames.to_csv(outdir / "combined_ul1_pocket_frame_counts.csv", index=False)
    clusters.to_csv(outdir / "combined_ul1_pocket_cluster_RS_summary.csv", index=False)

    metrics = [
        "median_min_to_UL1_A",
        "median_min_to_C4O1_A",
        "median_sidechain_centroid_to_C4_A",
        "median_sidechain_centroid_to_O1_A",
        "median_sidechain_centroid_to_phenyl_A",
        "median_ring_plane_angle_deg",
    ]
    all_delta = delta_table(all_by_rs, "all", metrics)
    if not all_delta.empty:
        all_delta = all_delta.rename(
            columns={"all_R_like": "R_like", "all_S_like": "S_like", "all_RminusS": "RminusS"}
        )
        all_delta["context"] = "all frames"
    near_delta = delta_table(near_rs, "near3A", metrics)
    if not near_delta.empty:
        near_delta = near_delta.rename(
            columns={"near3A_R_like": "R_like", "near3A_S_like": "S_like", "near3A_RminusS": "RminusS"}
        )
        near_delta["context"] = "C4-O1 <= 3.0 Å"
    deltas = pd.concat([all_delta, near_delta], ignore_index=True)
    deltas["engine_display"] = deltas["engine"].map(display_engine)
    deltas.to_csv(outdir / "combined_ul1_pocket_RS_median_feature_deltas.csv", index=False)

    rank = all_summary.groupby("reskey")["seen_fraction"].max().sort_values(ascending=False).head(args.top_n).index.tolist()
    make_presence_plot(all_summary, near_all, outdir / "combined_ul1_pocket_residue_presence_heatmap.png", args.top_n)
    make_delta_heatmap(
        deltas,
        "median_sidechain_centroid_to_C4_A",
        outdir / "combined_RS_delta_sidechain_centroid_to_C4_heatmap.png",
        rank,
        "R-like vs S-like pocket shifts: side-chain centroid to UL1 C4",
    )
    make_delta_heatmap(
        deltas,
        "median_sidechain_centroid_to_phenyl_A",
        outdir / "combined_RS_delta_sidechain_centroid_to_phenyl_heatmap.png",
        rank,
        "R-like vs S-like pocket shifts: side-chain centroid to UL1 phenyl",
    )
    make_delta_heatmap(
        deltas,
        "median_ring_plane_angle_deg",
        outdir / "combined_RS_delta_aromatic_ring_plane_angle_heatmap.png",
        rank,
        "R-like vs S-like pocket shifts: aromatic ring-plane angle",
    )
    make_cluster_plot(clusters, outdir / "combined_pocket_state_posterior_RS_near_attack.png")

    readme = [
        "# Full-frame UL1-centered pocket residue scan",
        "",
        "This folder combines GROMACS and Amber/pmemd outputs from `ul1_pocket_residue_scan.py`.",
        "",
        "Important interpretation points:",
        "",
        "- Every saved frame was used; no stride was applied.",
        "- Residues are selected by trajectory coordinates, analogous to a PyMOL `around 6 Å` selection around UL1, but computed for every frame.",
        "- R/S labels are posterior labels from the product-calibrated scalar triple product. They were not used to define PCA/UMAP/HDBSCAN pocket states.",
        "- `seen_fraction` means the fraction of frames where that residue entered the 6 Å UL1-centered pocket.",
        "- `near_seen_fraction` means the same quantity, but only among C4-O1 <= 3.0 Å frames.",
        "- `RminusS` in the delta table means R-like median minus S-like median for a feature. Negative distance deltas mean that feature is closer in R-like frames.",
        "- `min_to_*` values are per-frame closest-pair descriptors. They are included, but the workflow also includes side-chain centroids, phenyl centroid distances, residue classes, closest atom identities, and aromatic plane angles.",
        "",
        "Key files:",
        "",
        "- `combined_ul1_pocket_residue_summary_all_frames.csv`",
        "- `combined_ul1_pocket_residue_summary_by_RS_all_frames.csv`",
        "- `combined_ul1_pocket_residue_summary_near3A.csv`",
        "- `combined_ul1_pocket_residue_summary_by_RS_near3A.csv`",
        "- `combined_ul1_pocket_RS_median_feature_deltas.csv`",
        "- `combined_ul1_pocket_cluster_RS_summary.csv`",
    ]
    (outdir / "README_fullframe_UL1_pocket_scan.md").write_text("\n".join(readme) + "\n", encoding="utf-8")
    print(f"Done: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
