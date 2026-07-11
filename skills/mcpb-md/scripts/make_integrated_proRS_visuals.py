#!/usr/bin/env python3
"""Build integrated pro-R/pro-S visual summaries across GROMACS and Amber."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ANGSTROM = "\u00c5"
LE = "\u2264"
SQUARE_LARGE = (5.4, 5.4)
DATASET_COLORS = {"GROMACS": "#2b8cbe", "Amber pmemd": "#d95f02", "Combined": "#636363"}
PRO_COLORS = {"R-like": "#4c78a8", "S-like": "#e45756"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary-dir", required=True)
    parser.add_argument("--comparison-dir", required=True)
    parser.add_argument("--arial-font", default="/mnt/c/Windows/Fonts/arial.ttf")
    parser.add_argument("--near-cutoff", type=float, default=3.0)
    return parser.parse_args()


def setup_style(arial_font: str) -> str:
    font_name = "DejaVu Sans"
    font_path = Path(arial_font)
    if font_path.exists():
        fm.fontManager.addfont(str(font_path))
        font_name = fm.FontProperties(fname=str(font_path)).get_name()
    plt.rcParams.update(
        {
            "font.family": font_name,
            "font.size": 10.8,
            "axes.labelsize": 12.2,
            "axes.titlesize": 11.8,
            "xtick.labelsize": 10.8,
            "ytick.labelsize": 10.8,
            "legend.fontsize": 10.8,
            "axes.linewidth": 1.05,
            "xtick.major.width": 1.05,
            "ytick.major.width": 1.05,
            "xtick.major.size": 4.2,
            "ytick.major.size": 4.2,
            "savefig.facecolor": "white",
        }
    )
    return font_name


def savefig(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=360)
    plt.close(fig)


def top_legend(fig: plt.Figure, handles, labels, ncol: int = 2) -> None:
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=ncol,
        frameon=False,
        handlelength=1.6,
        handletextpad=0.45,
        columnspacing=1.35,
        markerscale=3.0,
        scatterpoints=1,
    )


def dataset_label(raw: str) -> str:
    return {"gromacs": "GROMACS", "amber_pmemd": "Amber pmemd"}.get(raw, raw)


def short_feature_label(feature: str) -> str:
    return (
        feature.replace("PHE336 ring-C4", "F336 ring-C4")
        .replace("PHE336/F336", "F336 min")
        .replace("PHE347/F347", "F347 min")
        .replace("Q269/GLN255", "Q269")
        .replace("H270/HD2", "H270")
    )


def load_rs_summary(comparison_dir: Path, near_cutoff: float) -> pd.DataFrame:
    rs = pd.read_csv(comparison_dir.parent.parent.parent.parent / "amber_200ns_pmemd_compare" / "analysis_vs_gromacs_manual_20260710" / "amber_vs_gromacs_pseudo_RS_summary.csv")
    # Fallback to a sibling path is intentionally avoided here; the manual summary is the authoritative lightweight source.
    wanted = {"all_frames", f"C4_O1_le_{near_cutoff:.2f}A", f"C4_O1_le_{near_cutoff:.1f}A"}
    sub = rs[rs["subset"].isin(wanted)].copy()
    sub["dataset_label"] = sub["dataset"].map(dataset_label)
    sub["subset_label"] = np.where(sub["subset"] == "all_frames", "All frames", f"C4-O1 {LE} {near_cutoff:.1f} {ANGSTROM}")
    sub.to_csv(comparison_dir / "integrated_proRS_fraction_summary.csv", index=False)
    return sub


def load_rs_summary_robust(summary_dir: Path, comparison_dir: Path, near_cutoff: float) -> pd.DataFrame:
    manual = comparison_dir.parents[3] / "amber_200ns_pmemd_compare" / "analysis_vs_gromacs_manual_20260710" / "amber_vs_gromacs_pseudo_RS_summary.csv"
    if not manual.exists():
        manual = Path(str(summary_dir).replace("/md_200ns/analysis_reactive_C4O1_20260708/publication_style_summary", "/amber_200ns_pmemd_compare/analysis_vs_gromacs_manual_20260710/amber_vs_gromacs_pseudo_RS_summary.csv"))
    rs = pd.read_csv(manual)
    wanted = {"all_frames", f"C4_O1_le_{near_cutoff:.2f}A", f"C4_O1_le_{near_cutoff:.1f}A"}
    sub = rs[rs["subset"].isin(wanted)].copy()
    sub["dataset_label"] = sub["dataset"].map(dataset_label)
    sub["subset_label"] = np.where(sub["subset"] == "all_frames", "All frames", f"C4-O1 {LE} {near_cutoff:.1f} {ANGSTROM}")
    sub.to_csv(comparison_dir / "integrated_proRS_fraction_summary.csv", index=False)
    return sub


def plot_fraction_and_model_performance(rs: pd.DataFrame, comparison_dir: Path, near_cutoff: float) -> None:
    group_metrics = pd.read_csv(comparison_dir / "gromacs_active_site_group_vs_proRS_metrics.csv")
    clf = pd.read_csv(comparison_dir / "pocket_feature_RS_classifier_metrics.csv")

    fig, axes = plt.subplots(2, 1, figsize=(5.4, 5.4), gridspec_kw={"height_ratios": [1.05, 1.0]})
    order = [
        ("GROMACS", "All frames"),
        ("GROMACS", f"C4-O1 {LE} {near_cutoff:.1f} {ANGSTROM}"),
        ("Amber pmemd", "All frames"),
        ("Amber pmemd", f"C4-O1 {LE} {near_cutoff:.1f} {ANGSTROM}"),
    ]
    labels = []
    r_vals = []
    s_vals = []
    ns = []
    for dataset, subset in order:
        row = rs[(rs["dataset_label"] == dataset) & (rs["subset_label"] == subset)]
        if row.empty:
            continue
        one = row.iloc[0]
        labels.append(f"{dataset}\n{subset}")
        r_vals.append(float(one["R_like_fraction"]) * 100.0)
        s_vals.append(float(one["S_like_fraction"]) * 100.0)
        ns.append(int(one["n_frames"]))
    x = np.arange(len(labels))
    axes[0].bar(x, r_vals, color=PRO_COLORS["R-like"], alpha=0.75, label="R-like")
    axes[0].bar(x, s_vals, bottom=r_vals, color=PRO_COLORS["S-like"], alpha=0.75, label="S-like")
    for i, n in enumerate(ns):
        axes[0].text(i, 103, f"n={n}", ha="center", va="bottom", fontsize=5.5)
    axes[0].set_ylim(0, 112)
    axes[0].set_ylabel("Frame composition (%)")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=0, fontsize=6.8)
    axes[0].tick_params(axis="y", labelsize=7.4)
    axes[0].yaxis.label.set_size(8.2)
    axes[0].text(0.02, 0.90, "R/S composition", transform=axes[0].transAxes, ha="left", va="top", fontsize=7.2)

    near_metric = group_metrics[group_metrics["subset"].str.contains(f"{near_cutoff:.1f}")].iloc[0]
    perf_labels = ["R-majority\nbaseline", "GROMACS\ngroups", "GROMACS\ncontact RF", "Amber\ncontact RF", "Combined\ncontact RF"]
    perf_vals = [
        float(near_metric["majority_baseline_accuracy"]) * 100.0,
        float(near_metric["group_majority_accuracy"]) * 100.0,
        float(clf.loc[clf["dataset"] == "GROMACS", "balanced_accuracy"].iloc[0]) * 100.0,
        float(clf.loc[clf["dataset"] == "Amber pmemd", "balanced_accuracy"].iloc[0]) * 100.0,
        float(clf.loc[clf["dataset"] == "Combined", "balanced_accuracy"].iloc[0]) * 100.0,
    ]
    colors = ["0.72", "#9ecae1", DATASET_COLORS["GROMACS"], DATASET_COLORS["Amber pmemd"], DATASET_COLORS["Combined"]]
    axes[1].bar(np.arange(len(perf_vals)), perf_vals, color=colors, alpha=0.84)
    axes[1].set_ylim(45, 85)
    axes[1].set_ylabel("Balanced/majority accuracy (%)")
    axes[1].set_xticks(np.arange(len(perf_vals)))
    axes[1].set_xticklabels(perf_labels, rotation=25, ha="right", fontsize=6.8)
    axes[1].tick_params(axis="y", labelsize=7.4)
    axes[1].yaxis.label.set_size(8.2)
    axes[1].text(0.02, 0.92, "Can features recover R/S?", transform=axes[1].transAxes, ha="left", va="top", fontsize=7.2)
    handles, leg_labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        leg_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=2,
        frameon=False,
        fontsize=7.4,
        handlelength=1.4,
        handletextpad=0.35,
        columnspacing=1.0,
    )
    fig.subplots_adjust(top=0.88, bottom=0.18, left=0.16, right=0.97, hspace=0.52)
    savefig(fig, comparison_dir / "pub_integrated_proRS_fraction_and_model_performance_Arial.png")


def plot_pocket_heatmap(comparison_dir: Path, near_cutoff: float, statistic: str = "median") -> None:
    dist = pd.read_csv(comparison_dir / "amber_vs_gromacs_proRS_pocket_distance_summary.csv")
    if statistic not in {"median", "mean"}:
        raise ValueError("statistic must be 'median' or 'mean'")
    value_col = f"{statistic}_distance_A"
    features = list(dict.fromkeys(dist["pocket_feature"]))
    columns = [
        ("GROMACS", "R-like"),
        ("GROMACS", "S-like"),
        ("Amber pmemd", "R-like"),
        ("Amber pmemd", "S-like"),
    ]
    matrix = np.full((len(features), len(columns)), np.nan)
    for i, feature in enumerate(features):
        for j, (dataset, label) in enumerate(columns):
            row = dist[(dist["dataset"] == dataset) & (dist["pocket_feature"] == feature) & (dist["pro_rs_label"] == label)]
            if not row.empty:
                matrix[i, j] = float(row[value_col].iloc[0])

    fig, ax = plt.subplots(figsize=SQUARE_LARGE)
    im = ax.imshow(matrix, cmap="viridis", aspect="auto")
    ax.set_xticks(np.arange(len(columns)))
    ax.set_xticklabels(["GMX\nR-like", "GMX\nS-like", "Amber\nR-like", "Amber\nS-like"])
    ax.set_yticks(np.arange(len(features)))
    ax.set_yticklabels([short_feature_label(f) for f in features])
    ax.tick_params(axis="x", pad=2)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if np.isfinite(matrix[i, j]):
                rgba = im.cmap(im.norm(matrix[i, j]))
                luminance = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
                text_color = "white" if luminance < 0.47 else "black"
                ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=10.0, color=text_color)
    cbar = fig.colorbar(im, ax=ax)
    label = "Median" if statistic == "median" else "Mean"
    cbar.set_label(f"{label} distance ({ANGSTROM})")
    ax.set_title(f"{label} pocket distances in near-attack frames, C4-O1 {LE} {near_cutoff:.1f} {ANGSTROM}", pad=7)
    fig.subplots_adjust(top=0.90, bottom=0.16, left=0.24, right=0.88)
    suffix = "" if statistic == "median" else "_mean"
    savefig(fig, comparison_dir / f"pub_integrated_proRS_pocket_distance_heatmap{suffix}_Arial.png")


def plot_delta_importance(comparison_dir: Path, near_cutoff: float) -> pd.DataFrame:
    delta = pd.read_csv(comparison_dir / "amber_vs_gromacs_proRS_pocket_delta.csv")
    importance = pd.read_csv(comparison_dir / "pocket_feature_RS_classifier_importance.csv")
    merged = delta.merge(importance[["dataset", "pocket_feature", "importance"]], on=["dataset", "pocket_feature"], how="inner")
    merged.to_csv(comparison_dir / "integrated_proRS_delta_importance_table.csv", index=False)

    fig, axes = plt.subplots(
        2,
        1,
        figsize=SQUARE_LARGE,
        gridspec_kw={"height_ratios": [1.15, 0.95]},
    )
    ax = axes[0]
    ax_zoom = axes[1]

    for dataset, sub in merged.groupby("dataset", sort=False):
        ax.scatter(
            sub["R_minus_S_median_A"],
            sub["importance"],
            s=80,
            color=DATASET_COLORS.get(dataset, "0.5"),
            alpha=0.78,
            edgecolors="white",
            linewidths=0.6,
            label=dataset,
        )
        for row in sub.itertuples():
            if row.pocket_feature != "PHE336 ring-C4":
                continue
            engine = "GMX" if row.dataset == "GROMACS" else "Amber"
            label = f"{engine} {short_feature_label(row.pocket_feature)}"
            dx = 0.015
            dy = -0.015 if row.dataset == "GROMACS" else 0.012
            ax.text(row.R_minus_S_median_A + dx, row.importance + dy, label, fontsize=6.9, va="center")
    ax.axvline(0.0, color="0.35", lw=0.85)
    ax.set_xlabel("")
    ax.set_ylabel("Random-forest feature importance")
    ax.set_title(f"Near attack C4-O1 {LE} {near_cutoff:.1f} {ANGSTROM}; negative x = closer in R-like", pad=7)
    ax.set_xlim(-0.56, 0.07)
    ax.set_ylim(0.0, 0.78)

    low = merged[merged["importance"] < 0.16].copy()
    label_offsets = {
        ("GROMACS", "PHE336/F336"): (8, 19),
        ("GROMACS", "PHE347/F347"): (-31, -16),
        ("GROMACS", "H270/HD2"): (8, 10),
        ("GROMACS", "Q269/GLN255"): (-44, 2),
        ("Amber pmemd", "PHE336/F336"): (8, 16),
        ("Amber pmemd", "PHE347/F347"): (-58, 6),
        ("Amber pmemd", "H270/HD2"): (20, 19),
        ("Amber pmemd", "Q269/GLN255"): (16, -17),
    }
    for dataset, sub in low.groupby("dataset", sort=False):
        ax_zoom.scatter(
            sub["R_minus_S_median_A"],
            sub["importance"],
            s=62,
            color=DATASET_COLORS.get(dataset, "0.5"),
            alpha=0.78,
            edgecolors="white",
            linewidths=0.55,
        )
        for row in sub.itertuples():
            if row.pocket_feature not in {"H270/HD2", "Q269/GLN255"}:
                continue
            engine = "GMX" if row.dataset == "GROMACS" else "Amber"
            label = f"{engine} {short_feature_label(row.pocket_feature)}"
            xytext = label_offsets.get((row.dataset, row.pocket_feature), (8, 8))
            ax_zoom.annotate(
                label,
                xy=(row.R_minus_S_median_A, row.importance),
                xytext=xytext,
                textcoords="offset points",
                ha="left" if xytext[0] >= 0 else "right",
                va="center",
                fontsize=6.1,
                arrowprops={"arrowstyle": "-", "lw": 0.45, "color": "0.35", "shrinkA": 0, "shrinkB": 4},
            )
    ax_zoom.axvline(0.0, color="0.35", lw=0.8)
    ax_zoom.set_xlim(-0.09, 0.035)
    ax_zoom.set_ylim(0.045, 0.135)
    ax_zoom.set_ylabel("Importance")
    ax_zoom.set_xlabel(f"Median distance difference, R-like - S-like ({ANGSTROM})")
    ax_zoom.set_title("Zoom: lower-importance contact features", pad=4)
    handles, labels = ax.get_legend_handles_labels()
    top_legend(fig, handles, labels, ncol=2)
    fig.subplots_adjust(top=0.84, bottom=0.14, left=0.17, right=0.97, hspace=0.58)
    savefig(fig, comparison_dir / "pub_integrated_proRS_delta_vs_importance_Arial.png")
    return merged


def write_readme(comparison_dir: Path, font_used: str, merged: pd.DataFrame) -> None:
    top = merged.sort_values("importance", ascending=False).head(6)
    lines = [
        "# Integrated Pro-R/Pro-S Visuals",
        "",
        f"Font requested/used by Matplotlib: `{font_used}`.",
        "",
        "These integrated plots combine the GROMACS R/S split from `publication_style_summary` with the Amber-vs-GROMACS comparison data. They are meant to answer three linked questions:",
        "",
        "1. Are R-like and S-like frames populated similarly in GROMACS and Amber?",
        "2. Do the old unsupervised Group 1/2/3 labels actually recover R/S?",
        "3. Which pocket contact features both differ between R/S and help classify R/S?",
        "",
        "## Figures",
        "",
        "- `pub_integrated_proRS_fraction_and_model_performance_Arial.png`: stacked R/S composition plus a direct comparison of group-based classification versus contact-feature random forest.",
        "- `pub_integrated_proRS_pocket_distance_heatmap_Arial.png`: R-like/S-like median pocket distances for GROMACS and Amber on the same heatmap.",
        "- `pub_integrated_proRS_pocket_distance_heatmap_mean_Arial.png`: the same view using mean distances.",
        "- `pub_integrated_proRS_delta_vs_importance_Arial.png`: combines physical direction (`R-like - S-like` distance) with statistical classifier importance.",
        "",
        "## Interpretation",
        "",
        "The key distinction is that Group 1/2/3 are unsupervised geometry groups, while R-like/S-like labels come from the product-calibrated scalar triple product. The groups partially enrich R/S, but they are not a clean R/S classifier.",
        "",
        "The pocket heatmaps use per-frame contact features. For labels such as `F336 min`, the value in each frame is the minimum distance from the selected residue atoms to the reactive C4/O1 atom pair, not the minimum over the full trajectory. The heatmap statistic is then computed across frames. The median version is robust to occasional excursions; the mean version is a conventional average and is more sensitive to tails.",
        "",
        "The most useful integrated view is the delta-vs-importance plot. A residue/contact feature is strongest when it has both a sizable R-S distance difference and high classifier importance. In the current WT-like system, the F336 ring-C4 feature dominates by this criterion.",
        "",
        "For future mutant series, the F336 ring-centroid feature should be replaced by mutation-compatible pocket-position descriptors, because F336R no longer has an aromatic ring. A robust cross-mutant model should center each frame on UL1/UNL1 and represent each pocket position by distances/contact occupancies to C4, O1, and the substrate phenyl centroid, plus residue class descriptors such as charge, polarity, aromaticity, and side-chain centroid/guanidinium contact where appropriate.",
        "",
        "Top integrated delta/importance rows:",
        "",
        top.to_string(index=False),
        "",
    ]
    (comparison_dir / "README_integrated_proRS_visuals.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    summary_dir = Path(args.summary_dir)
    comparison_dir = Path(args.comparison_dir)
    comparison_dir.mkdir(parents=True, exist_ok=True)
    font_used = setup_style(args.arial_font)
    rs = load_rs_summary_robust(summary_dir, comparison_dir, args.near_cutoff)
    plot_fraction_and_model_performance(rs, comparison_dir, args.near_cutoff)
    plot_pocket_heatmap(comparison_dir, args.near_cutoff, statistic="median")
    plot_pocket_heatmap(comparison_dir, args.near_cutoff, statistic="mean")
    merged = plot_delta_importance(comparison_dir, args.near_cutoff)
    write_readme(comparison_dir, font_used, merged)
    print(comparison_dir)


if __name__ == "__main__":
    main()
