#!/usr/bin/env python3
"""Summarize deep UL1-pocket state-search runs across engines."""

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
    parser = argparse.ArgumentParser(description="Combine deep UL1-pocket state-search outputs.")
    parser.add_argument("--gromacs-dir", required=True)
    parser.add_argument("--amber-dir", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def find_one(folder: Path, suffix: str) -> Path:
    hits = sorted(folder.glob(f"*{suffix}"))
    if len(hits) != 1:
        raise FileNotFoundError(f"Expected one *{suffix} in {folder}, found {len(hits)}")
    return hits[0]


def load_run(folder: Path) -> dict[str, pd.DataFrame]:
    scores = pd.read_csv(find_one(folder, "_all_deep_feature_hdbscan_scores.csv"))
    details = pd.read_csv(find_one(folder, "_all_deep_feature_best_state_details.csv"))
    meta = pd.read_csv(find_one(folder, "_deep_state_frame_metadata.csv"))
    engine = str(scores["engine"].iloc[0])
    embeds = {}
    for method in sorted(scores["method"].unique()):
        path = folder / f"{engine}_{method}_best_state_embeddings.csv.gz"
        if path.exists():
            embeds[method] = pd.read_csv(path)
    return {"scores": scores, "details": details, "meta": meta, "embeds": embeds}


def display_engine(engine: str) -> str:
    return ENGINE_DISPLAY.get(engine, engine)


def best_by_method(scores: pd.DataFrame) -> pd.DataFrame:
    return scores.sort_values("score", ascending=False).groupby(["engine", "method"], as_index=False).head(1)


def enriched_states(details: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    engine = str(details["engine"].iloc[0])
    global_r = float((meta["RS_label"] == "R-like").mean())
    rows = []
    for method, sub in details.groupby("method"):
        # Use only the states associated with that method's best parameter set.
        sub = sub.copy()
        sub["global_R_like_fraction"] = global_r
        sub["abs_R_enrichment"] = (sub["R_like_fraction"] - global_r).abs()
        sub["R_enrichment"] = sub["R_like_fraction"] - global_r
        sub["engine_display"] = display_engine(engine)
        rows.append(sub.sort_values("abs_R_enrichment", ascending=False).head(12))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def time_block_stability(embeds: dict[str, pd.DataFrame], details: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    engine = str(details["engine"].iloc[0])
    global_r = float((meta["RS_label"] == "R-like").mean())
    for method, embed in embeds.items():
        if embed.empty or "best_state" not in embed.columns:
            continue
        det = details[details["method"] == method].copy()
        det["abs_R_enrichment"] = (det["R_like_fraction"] - global_r).abs()
        top_states = det.sort_values("abs_R_enrichment", ascending=False).head(8)["state"].tolist()
        embed = embed.copy()
        embed["time_block"] = pd.qcut(np.arange(len(embed)), q=5, labels=False, duplicates="drop")
        for state in top_states:
            sub_state = embed[embed["best_state"] == state]
            for block, sub in sub_state.groupby("time_block"):
                frames = len(sub)
                rows.append(
                    {
                        "engine": engine,
                        "engine_display": display_engine(engine),
                        "method": method,
                        "state": int(state),
                        "time_block": int(block),
                        "frames": frames,
                        "R_like_fraction": float((sub["RS_label"] == "R-like").mean()) if frames else np.nan,
                        "global_R_like_fraction": global_r,
                        "median_C4_O1_A": float(sub["C4_O1_A"].median()) if frames else np.nan,
                    }
                )
    return pd.DataFrame(rows)


def plot_best_scores(best: pd.DataFrame, out: Path) -> None:
    if best.empty:
        return
    methods = ["local_coords", "ul1_distance_matrix", "shape_grid", "pair_network"]
    engines = [eng for eng in ["gromacs_full", "amber_pmemd_full"] if eng in set(best["engine"])]
    x = np.arange(len(methods))
    width = 0.38
    fig, ax = plt.subplots(figsize=(9.2, 5.4), dpi=240)
    for i, eng in enumerate(engines):
        sub = best[best["engine"] == eng].set_index("method")
        values = [sub.loc[m, "score"] if m in sub.index else 0.0 for m in methods]
        ax.bar(x + (i - 0.5) * width, values, width=width, label=display_engine(eng))
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=25, ha="right", fontsize=10)
    ax.set_ylabel("Best posterior R/S enrichment score", fontsize=12)
    ax.set_title("Deep pocket-state feature loops", fontsize=15)
    ax.legend(frameon=False, fontsize=10)
    ax.tick_params(labelsize=10)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_metric_bars(best: pd.DataFrame, out: Path) -> None:
    if best.empty:
        return
    methods = ["local_coords", "ul1_distance_matrix", "shape_grid", "pair_network"]
    engines = [eng for eng in ["gromacs_full", "amber_pmemd_full"] if eng in set(best["engine"])]
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.8), dpi=240, sharex=True)
    metrics = [("purity", "Purity"), ("NMI", "NMI"), ("coverage", "Coverage")]
    x = np.arange(len(methods))
    width = 0.38
    for ax, (metric, label) in zip(axes, metrics):
        for i, eng in enumerate(engines):
            sub = best[best["engine"] == eng].set_index("method")
            values = [sub.loc[m, metric] if m in sub.index else 0.0 for m in methods]
            ax.bar(x + (i - 0.5) * width, values, width=width, label=display_engine(eng))
        ax.set_title(label, fontsize=12)
        ax.set_xticks(x)
        ax.set_xticklabels(methods, rotation=30, ha="right", fontsize=9)
        ax.tick_params(labelsize=9)
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle("Best HDBSCAN state quality by feature loop", fontsize=15)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_enriched_states(states: pd.DataFrame, out: Path) -> None:
    if states.empty:
        return
    best_methods = (
        states.sort_values("score", ascending=False)
        .groupby("engine", as_index=False)
        .head(1)[["engine", "method"]]
        .drop_duplicates()
    )
    rows = []
    for row in best_methods.itertuples(index=False):
        sub = states[(states["engine"] == row.engine) & (states["method"] == row.method)]
        rows.append(sub.sort_values("abs_R_enrichment", ascending=False).head(12))
    plot_df = pd.concat(rows, ignore_index=True)
    plot_df["label"] = (
        plot_df["engine_display"]
        + "\n"
        + plot_df["method"]
        + " state "
        + plot_df["state"].astype(str)
    )
    plot_df = plot_df.sort_values(["engine_display", "R_like_fraction"], ascending=[True, True])
    fig, ax = plt.subplots(figsize=(9.4, 8.5), dpi=240)
    colors = np.where(plot_df["R_like_fraction"] >= plot_df["global_R_like_fraction"], "#D55E00", "#0072B2")
    ax.barh(plot_df["label"], plot_df["R_like_fraction"] * 100.0, color=colors)
    for i, row in enumerate(plot_df.itertuples(index=False)):
        ax.text(row.R_like_fraction * 100.0 + 1.0, i, f"n={row.frames}", va="center", fontsize=8)
    ax.set_xlim(0, 112)
    ax.set_xlabel("R-like fraction in state (%)", fontsize=12)
    ax.set_title("Most R/S-enriched states from best feature loops", fontsize=14)
    ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    runs = [load_run(Path(args.gromacs_dir)), load_run(Path(args.amber_dir))]
    scores = pd.concat([run["scores"] for run in runs], ignore_index=True)
    details = pd.concat([run["details"] for run in runs], ignore_index=True)
    best = best_by_method(scores).copy()
    best["engine_display"] = best["engine"].map(display_engine)
    best = best.sort_values("score", ascending=False)

    states = []
    block_rows = []
    for run in runs:
        states.append(enriched_states(run["details"], run["meta"]))
        block_rows.append(time_block_stability(run["embeds"], run["details"], run["meta"]))
    state_df = pd.concat(states, ignore_index=True) if states else pd.DataFrame()
    block_df = pd.concat(block_rows, ignore_index=True) if block_rows else pd.DataFrame()

    scores.to_csv(outdir / "combined_deep_state_all_hdbscan_scores.csv", index=False)
    details.to_csv(outdir / "combined_deep_state_all_best_state_details.csv", index=False)
    best.to_csv(outdir / "combined_deep_state_best_by_method.csv", index=False)
    state_df.to_csv(outdir / "combined_deep_state_top_RS_enriched_states.csv", index=False)
    block_df.to_csv(outdir / "combined_deep_state_time_block_stability.csv", index=False)

    plot_best_scores(best, outdir / "combined_deep_state_best_scores.png")
    plot_metric_bars(best, outdir / "combined_deep_state_best_metrics.png")
    plot_enriched_states(state_df, outdir / "combined_deep_state_top_enriched_states.png")

    readme = [
        "# Combined deep UL1-pocket state search",
        "",
        "This summary compares four feature loops:",
        "",
        "- `local_coords`: residue centroids/closest atoms/ring normals in a UL1-fixed local 3D frame.",
        "- `ul1_distance_matrix`: residue-to-every-UL1-heavy-atom distance fingerprints.",
        "- `pair_network`: pairwise residue-centroid distances plus atom-cloud moments.",
        "- `shape_grid`: coarse 3D atom-density grids around UL1.",
        "",
        "The clusters are unsupervised HDBSCAN states from PCA/UMAP representations. R/S labels are overlaid only afterward.",
        "",
        "Important files:",
        "",
        "- `combined_deep_state_best_by_method.csv`: best parameter set per engine and feature loop.",
        "- `combined_deep_state_top_RS_enriched_states.csv`: states most enriched for R-like or S-like frames.",
        "- `combined_deep_state_time_block_stability.csv`: whether enriched states persist across five time blocks.",
        "- `combined_deep_state_best_scores.png`: score comparison.",
    ]
    (outdir / "README_combined_deep_state_search.md").write_text("\n".join(readme) + "\n", encoding="utf-8")
    print(f"Done: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
