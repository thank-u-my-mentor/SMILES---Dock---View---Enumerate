#!/usr/bin/env python3
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
import pandas as pd


OUT = Path("/mnt/e/TJ/260706_mcpb_m2_200ns/amber_200ns_pmemd_compare/analysis_vs_gromacs_manual_20260710")


def setup_style() -> None:
    arial = Path("/mnt/c/Windows/Fonts/arial.ttf")
    if arial.exists():
        font_manager.fontManager.addfont(str(arial))
        mpl.rcParams["font.family"] = "Arial"
    else:
        mpl.rcParams["font.family"] = "DejaVu Sans"
    mpl.rcParams["axes.unicode_minus"] = False
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42


def plot_summary() -> None:
    df = pd.read_csv(OUT / "amber_vs_gromacs_pseudo_RS_summary.csv")
    keep = ["all_frames", "C4_O1_le_3.00A", "C4_O1_le_2.80A", "top_100_shortest_C4_O1"]
    labels = {
        "all_frames": "All frames",
        "C4_O1_le_3.00A": "C4-O1 <= 3.0 A",
        "C4_O1_le_2.80A": "C4-O1 <= 2.8 A",
        "top_100_shortest_C4_O1": "Top 100 shortest",
    }
    sub = df[df["subset"].isin(keep)].copy()
    sub["subset_label"] = sub["subset"].map(labels)
    sub["dataset_label"] = sub["dataset"].map({"gromacs": "GROMACS", "amber_pmemd": "Amber pmemd"})

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), dpi=300)
    colors = {"GROMACS": "#4C78A8", "Amber pmemd": "#F58518"}
    x = range(len(keep))
    width = 0.36
    for offset, dataset in [(-width / 2, "GROMACS"), (width / 2, "Amber pmemd")]:
        d = sub[sub["dataset_label"] == dataset].set_index("subset").loc[keep]
        axes[0].bar([i + offset for i in x], d["R_like_fraction"] * 100, width=width, label=dataset, color=colors[dataset])
        axes[1].bar([i + offset for i in x], d["C4_O1_median_A"], width=width, label=dataset, color=colors[dataset])

    axes[0].set_ylabel("R-like fraction (%)")
    axes[0].set_ylim(0, 100)
    axes[1].set_ylabel("Median C4-O1 distance (A)")
    axes[1].set_ylim(2.5, 4.6)
    for ax in axes:
        ax.set_xticks(list(x))
        ax.set_xticklabels([labels[k] for k in keep], rotation=25, ha="right")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "amber_vs_gromacs_near_attack_Rlike_summary.png")
    fig.savefig(OUT / "amber_vs_gromacs_near_attack_Rlike_summary.pdf")
    plt.close(fig)


def plot_crosstab() -> None:
    df = pd.read_csv(OUT / "face_vs_pseudo_RS_crosstab.csv")
    sub = df[df["subset"] == "top_100_shortest"].copy()
    sub["dataset_label"] = sub["dataset"].map({"gromacs": "GROMACS", "amber_pmemd": "Amber pmemd"})
    sub["class_label"] = sub["face"] + " / " + sub["pseudo_RS"].str.replace("_by_product_reference", "", regex=False).str.replace("_like", "-like", regex=False)
    order = [
        "same / R-like",
        "same / S-like",
        "opposite / R-like",
        "opposite / S-like",
    ]
    colors = ["#9ecae9", "#fdd0a2", "#3182bd", "#e6550d"]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), dpi=300, sharey=True)
    for ax, dataset in zip(axes, ["GROMACS", "Amber pmemd"]):
        d = sub[sub["dataset_label"] == dataset].set_index("class_label").reindex(order).fillna(0)
        ax.bar(order, d["fraction_of_subset"] * 100, color=colors)
        ax.set_title(dataset)
        ax.set_ylim(0, 100)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=30, ha="right")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("Fraction in top 100 shortest C4-O1 frames (%)")
    fig.tight_layout()
    fig.savefig(OUT / "face_vs_pseudo_RS_top100_crosstab.png")
    fig.savefig(OUT / "face_vs_pseudo_RS_top100_crosstab.pdf")
    plt.close(fig)


def main() -> None:
    setup_style()
    plot_summary()
    plot_crosstab()
    print(f"Wrote plots to {OUT}")


if __name__ == "__main__":
    main()
