#!/usr/bin/env python3
"""Compare Amber pmemd and GROMACS MD summaries for the TJ MCPB project."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, adjusted_mutual_info_score, balanced_accuracy_score, normalized_mutual_info_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict


ANGSTROM = "\u00c5"
LE = "\u2264"
SQUARE = (4.4, 4.4)
SQUARE_READABLE = (5.4, 5.4)
SQUARE_HIST = (6.8, 6.8)
DATASET_COLORS = {"GROMACS": "#2b8cbe", "Amber pmemd": "#d95f02"}
PRO_COLORS = {"R-like": "#4c78a8", "S-like": "#e45756"}
TRIPLE_COL = "signed_triple_O1_C5phenyl_C3chain_about_C4"
PRO_RS_COL = "pseudo_RS_by_product_ref"
PRO_RS_LABELS = {
    "R_like_by_product_reference": "R-like",
    "S_like_by_product_reference": "S-like",
}
POCKET_FEATURES = [
    ("PHE322_min_to_C4O1_A", "PHE336/F336"),
    ("PHE333_min_to_C4O1_A", "PHE347/F347"),
    ("HD2_256_min_to_C4O1_A", "H270/HD2"),
    ("GLN255_min_to_C4O1_A", "Q269/GLN255"),
    ("PHE322_ring_centroid_to_C4_A", "PHE336 ring-C4"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmx-summary-dir", required=True)
    parser.add_argument("--amber-dir", required=True)
    parser.add_argument("--outdir", required=True)
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
            "axes.labelsize": 12.5,
            "axes.titlesize": 12.0,
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


def safe_dataset_name(dataset: str) -> str:
    return dataset.lower().replace(" ", "_")


def top_legend(fig: plt.Figure, handles, labels, ncol: int = 2) -> None:
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=ncol,
        frameon=False,
        handlelength=1.8,
        handletextpad=0.45,
        columnspacing=1.6,
        markerscale=3.2,
        scatterpoints=1,
    )


def make_histogram_axis_readable(ax: plt.Axes) -> None:
    ax.tick_params(axis="both", labelsize=17.0, width=1.35, length=5.8)
    ax.xaxis.label.set_size(20.0)
    ax.yaxis.label.set_size(20.0)
    ax.title.set_size(20.0)
    for spine in ax.spines.values():
        spine.set_linewidth(1.35)


def read_cpptraj_dat(path: Path, value_name: str, frame_dt_ps: float = 5.0) -> pd.DataFrame:
    rows = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("@"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            frame = int(float(parts[0]))
            rows.append((frame, (frame - 1) * frame_dt_ps / 1000.0, float(parts[1])))
    return pd.DataFrame(rows, columns=["frame", "time_ns", value_name])


def load_rmsd(gmx_dir: Path, amber_dir: Path) -> pd.DataFrame:
    gmx = read_cpptraj_dat(gmx_dir / "rmsd_source" / "gmx_CA_RMSD_pbcfit.dat", "rmsd_A")
    gmx["dataset"] = "GROMACS"
    amber = read_cpptraj_dat(amber_dir / "amber_CA_RMSD_A.dat", "rmsd_A")
    amber["dataset"] = "Amber pmemd"
    return pd.concat([gmx, amber], ignore_index=True)


def load_c4o1(gmx_dir: Path, amber_dir: Path) -> pd.DataFrame:
    cv_path = gmx_dir / "cv_source" / "reactive_geometry_CVs.csv"
    if cv_path.exists():
        gmx_source = pd.read_csv(cv_path)
    else:
        gmx_source = pd.read_csv(
            gmx_dir.parent
            / "ul1_rotamer_face_analysis"
            / "product_calibrated_pseudo_RS"
            / "md_product_calibrated_pseudo_RS_timeseries.csv"
        )
    gmx = gmx_source[["frame", "time_ns", "C4_O1_A"]].copy()
    gmx["dataset"] = "GROMACS"
    amber = read_cpptraj_dat(amber_dir / "amber_C4_O1_distance_A.dat", "C4_O1_A")
    amber["dataset"] = "Amber pmemd"
    return pd.concat([gmx, amber], ignore_index=True)


def load_pro_rs(gmx_dir: Path, amber_dir: Path) -> pd.DataFrame:
    gmx = pd.read_csv(gmx_dir.parent / "ul1_rotamer_face_analysis" / "product_calibrated_pseudo_RS" / "md_product_calibrated_pseudo_RS_timeseries.csv")
    gmx["dataset"] = "GROMACS"
    amber = pd.read_csv(amber_dir / "amber_product_calibrated_pseudo_RS_timeseries.csv")
    amber["dataset"] = "Amber pmemd"
    return pd.concat([gmx, amber], ignore_index=True)


def read_rmsf(path: Path, dataset: str) -> pd.DataFrame:
    rows = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("@"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            rows.append((int(float(parts[0])), float(parts[1]), dataset))
    return pd.DataFrame(rows, columns=["residue_index", "rmsf_A", "dataset"])


def load_rmsf(outdir: Path) -> pd.DataFrame | None:
    source = outdir / "rmsf_source"
    gmx_path = source / "gromacs_CA_RMSF_byres.dat"
    amber_path = source / "amber_CA_RMSF_byres.dat"
    if not gmx_path.exists() or not amber_path.exists():
        return None
    return pd.concat(
        [
            read_rmsf(gmx_path, "GROMACS"),
            read_rmsf(amber_path, "Amber pmemd"),
        ],
        ignore_index=True,
    )


def read_cpptraj_series(path: Path, value_name: str, frame_offset: int = 0) -> pd.DataFrame:
    rows = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("@"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            rows.append((int(float(parts[0])) + frame_offset, float(parts[1])))
    return pd.DataFrame(rows, columns=["frame", value_name])


def load_amber_pocket_features(outdir: Path, amber_dir: Path) -> pd.DataFrame | None:
    source = outdir / "pocket_source"
    feature_files = {
        "PHE322_min_to_C4O1_A": source / "amber_PHE322_min_to_C4O1.dat",
        "PHE333_min_to_C4O1_A": source / "amber_PHE333_min_to_C4O1.dat",
        "HD2_256_min_to_C4O1_A": source / "amber_HD2_256_min_to_C4O1.dat",
        "GLN255_min_to_C4O1_A": source / "amber_GLN255_min_to_C4O1.dat",
        "PHE322_ring_centroid_to_C4_A": source / "amber_PHE322_ring_centroid_to_C4.dat",
    }
    if any(not path.exists() for path in feature_files.values()):
        return None

    merged = None
    for col, path in feature_files.items():
        # cpptraj writes frames as 1..N, while the calibrated Amber R/S table is 0..N-1.
        one = read_cpptraj_series(path, col, frame_offset=-1)
        merged = one if merged is None else merged.merge(one, on="frame", how="inner")
    pro = pd.read_csv(amber_dir / "amber_product_calibrated_pseudo_RS_timeseries.csv")
    pro["pro_rs_label"] = pro["pseudo_RS_by_product_ref"].map(
        {
            "R_like_by_product_reference": "R-like",
            "S_like_by_product_reference": "S-like",
        }
    )
    merged = merged.merge(pro[["frame", "time_ns", "C4_O1_A", "pro_rs_label"]], on="frame", how="inner")
    merged["dataset"] = "Amber pmemd"
    merged["sampling"] = "all frames"
    return merged


def load_gromacs_pocket_features(gmx_dir: Path) -> pd.DataFrame:
    group_frames = pd.read_csv(gmx_dir / "active_site_group_frame_assignments.csv")
    pro = pd.read_csv(gmx_dir.parent / "ul1_rotamer_face_analysis" / "product_calibrated_pseudo_RS" / "md_product_calibrated_pseudo_RS_timeseries.csv")
    pro["pro_rs_label"] = pro["pseudo_RS_by_product_ref"].map(
        {
            "R_like_by_product_reference": "R-like",
            "S_like_by_product_reference": "S-like",
        }
    )
    # The active-site feature clustering used every 10th GROMACS frame, so pocket features here are strided.
    merged = group_frames.merge(pro[["frame", "time_ns", "C4_O1_A", "pro_rs_label"]], on="frame", how="inner", suffixes=("", "_proRS"))
    merged["PHE322_ring_centroid_to_C4_A"] = merged["PHE322_ring_centroid_to_C4_A"].astype(float)
    merged["dataset"] = "GROMACS"
    merged["sampling"] = "stride 10"
    return merged


def running_median(y: np.ndarray, window: int = 401) -> np.ndarray:
    return pd.Series(y).rolling(window, center=True, min_periods=1).median().to_numpy()


def plot_rmsd(rmsd: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    fig, ax = plt.subplots(figsize=SQUARE)
    handles = []
    labels = []
    rows = []
    for dataset, sub in rmsd.groupby("dataset", sort=False):
        sub = sub.sort_values("time_ns")
        color = DATASET_COLORS[dataset]
        ax.plot(sub["time_ns"], sub["rmsd_A"], color=color, lw=0.18, alpha=0.20)
        line = ax.plot(sub["time_ns"], running_median(sub["rmsd_A"].to_numpy(float)), color=color, lw=1.25, label=dataset)[0]
        handles.append(line)
        labels.append(dataset)
        last50 = sub[sub["time_ns"] >= sub["time_ns"].max() - 50.0]
        rows.append(
            {
                "dataset": dataset,
                "n_frames": int(len(sub)),
                "initial_RMSD_A": float(sub["rmsd_A"].iloc[0]),
                "RMSD_1ns_A": float(sub.iloc[(sub["time_ns"] - 1.0).abs().argmin()]["rmsd_A"]),
                "median_RMSD_A": float(sub["rmsd_A"].median()),
                "last50ns_median_RMSD_A": float(last50["rmsd_A"].median()),
                "last_RMSD_A": float(sub["rmsd_A"].iloc[-1]),
            }
        )
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(f"Protein C-alpha RMSD ({ANGSTROM})")
    ax.set_xlim(0, max(rmsd["time_ns"]))
    top_legend(fig, handles, labels, ncol=2)
    fig.subplots_adjust(top=0.84, bottom=0.16, left=0.18, right=0.96)
    savefig(fig, outdir / "pub_amber_vs_gromacs_CA_RMSD_Arial.png")
    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "amber_vs_gromacs_RMSD_summary.csv", index=False)
    return summary


def plot_c4o1(c4o1: pd.DataFrame, outdir: Path, near_cutoff: float) -> pd.DataFrame:
    fig, ax = plt.subplots(figsize=SQUARE_HIST)
    rows = []
    bins = np.linspace(2.5, 7.0, 91)
    for dataset, sub in c4o1.groupby("dataset", sort=False):
        y = sub["C4_O1_A"].to_numpy(float)
        weights = np.ones_like(y) * 100.0 / len(y)
        ax.hist(y, bins=bins, weights=weights, histtype="step", lw=2.1, color=DATASET_COLORS[dataset], label=dataset)
        rows.append(
            {
                "dataset": dataset,
                "n_frames": int(len(y)),
                "median_C4_O1_A": float(np.median(y)),
                "mean_C4_O1_A": float(np.mean(y)),
                "min_C4_O1_A": float(np.min(y)),
                f"fraction_C4O1_le_{near_cutoff:.1f}A": float(np.mean(y <= near_cutoff)),
            }
        )
    ax.axvline(near_cutoff, color="0.35", ls=":", lw=1.25)
    ax.set_xlabel(f"UL1 C4-O1 distance ({ANGSTROM})")
    ax.set_ylabel("Frequency (%)")
    make_histogram_axis_readable(ax)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=2,
        frameon=False,
        fontsize=17.0,
        handlelength=2.2,
        handletextpad=0.55,
        columnspacing=1.7,
    )
    fig.subplots_adjust(top=0.84, bottom=0.17, left=0.17, right=0.96)
    savefig(fig, outdir / "pub_amber_vs_gromacs_C4O1_histogram_Arial.png")
    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "amber_vs_gromacs_C4O1_summary.csv", index=False)
    return summary


def plot_single_c4o1_histograms(c4o1: pd.DataFrame, outdir: Path, near_cutoff: float) -> None:
    for dataset, sub in c4o1.groupby("dataset", sort=False):
        y = sub["C4_O1_A"].to_numpy(float)
        fig, ax = plt.subplots(figsize=SQUARE_HIST)
        weights = np.ones_like(y) * 100.0 / len(y)
        ax.hist(y, bins=85, weights=weights, color="#8ecae6", edgecolor="#4a90b5", linewidth=0.75)
        ax.set_xlabel(f"UL1 C4-O1 distance ({ANGSTROM})")
        ax.set_ylabel("Frequency (%)")
        ax.set_title(dataset, pad=4)
        make_histogram_axis_readable(ax)
        fig.subplots_adjust(top=0.88, bottom=0.17, left=0.17, right=0.96)
        savefig(fig, outdir / f"pub_{safe_dataset_name(dataset)}_C4O1_histogram_Arial.png")


def add_pro_rs_label(data: pd.DataFrame) -> pd.DataFrame:
    out = data.copy()
    out["pro_rs_label"] = out[PRO_RS_COL].map(PRO_RS_LABELS).fillna(out[PRO_RS_COL])
    return out


def draw_triple_scatter(ax: plt.Axes, data: pd.DataFrame, near_cutoff: float, title: str | None = None) -> tuple[list, list]:
    handles = []
    labels = []
    for label in ["R-like", "S-like"]:
        sub = data[data["pro_rs_label"] == label]
        sc = ax.scatter(
            sub[TRIPLE_COL],
            sub["C4_O1_A"],
            s=5.2,
            alpha=0.15,
            linewidths=0,
            color=PRO_COLORS[label],
            label=label,
        )
        handles.append(sc)
        labels.append(label)
    ax.axvline(0.0, color="0.35", lw=1.05, ls="--")
    ax.axhline(near_cutoff, color="0.35", lw=1.05, ls=":")
    ax.set_xlim(-8.6, 8.6)
    ax.set_ylim(2.5, 5.15)
    if title:
        ax.set_title(title, pad=5)
    return handles, labels


def plot_pro_rs_scatter_by_dataset(pro_rs: pd.DataFrame, outdir: Path, near_cutoff: float) -> None:
    data = add_pro_rs_label(pro_rs)
    for dataset, sub in data.groupby("dataset", sort=False):
        fig, ax = plt.subplots(figsize=SQUARE_READABLE)
        handles, labels = draw_triple_scatter(ax, sub, near_cutoff, dataset)
        ax.set_xlabel(r"Triple product V around C4 ($\AA^3$)")
        ax.set_ylabel(f"C4-O1 distance ({ANGSTROM})")
        top_legend(fig, handles, labels, ncol=2)
        fig.subplots_adjust(top=0.84, bottom=0.16, left=0.18, right=0.96)
        savefig(fig, outdir / f"pub_{safe_dataset_name(dataset)}_proRS_triple_product_vs_C4O1_Arial.png")

    datasets = [x for x in ["GROMACS", "Amber pmemd"] if x in set(data["dataset"])]
    fig, axes = plt.subplots(len(datasets), 1, figsize=(5.4, 6.4), sharex=True, sharey=True)
    if len(datasets) == 1:
        axes = [axes]
    handles = labels = None
    for ax, dataset in zip(axes, datasets):
        sub = data[data["dataset"] == dataset]
        handles, labels = draw_triple_scatter(ax, sub, near_cutoff, dataset)
        ax.set_ylabel(f"C4-O1 ({ANGSTROM})")
    axes[-1].set_xlabel(r"Triple product V around C4 ($\AA^3$)")
    if handles is not None and labels is not None:
        top_legend(fig, handles, labels, ncol=2)
    fig.subplots_adjust(top=0.88, bottom=0.12, left=0.16, right=0.96, hspace=0.24)
    savefig(fig, outdir / "pub_amber_vs_gromacs_proRS_triple_product_vs_C4O1_Arial.png")


def plot_pro_rs_thresholds(amber_dir: Path, outdir: Path) -> pd.DataFrame:
    summary = pd.read_csv(amber_dir / "amber_vs_gromacs_pseudo_RS_summary.csv")
    cut = summary[summary["subset"].str.match(r"C4_O1_le_[0-9.]+A")].copy()
    cut["cutoff_A"] = cut["subset"].str.extract(r"C4_O1_le_([0-9.]+)A").astype(float)
    cut["dataset_label"] = cut["dataset"].map({"gromacs": "GROMACS", "amber_pmemd": "Amber pmemd"}).fillna(cut["dataset"])
    cut["low_count"] = cut["n_frames"] < 100

    fig, ax = plt.subplots(figsize=SQUARE)
    handles = []
    labels = []
    ax.axvspan(2.65, 2.805, color="0.92", zorder=0)
    ax.text(2.725, 8.0, "low-n tail", ha="center", va="bottom", fontsize=7, color="0.35")
    for dataset, sub in cut.groupby("dataset_label", sort=False):
        sub = sub.sort_values("cutoff_A")
        line = ax.plot(sub["cutoff_A"], sub["R_like_fraction"] * 100.0, lw=1.15, color=DATASET_COLORS[dataset], label=dataset)[0]
        confident = sub[~sub["low_count"]]
        low = sub[sub["low_count"]]
        ax.scatter(confident["cutoff_A"], confident["R_like_fraction"] * 100.0, s=18, color=DATASET_COLORS[dataset], zorder=3)
        ax.scatter(
            low["cutoff_A"],
            low["R_like_fraction"] * 100.0,
            s=22,
            facecolors="white",
            edgecolors=DATASET_COLORS[dataset],
            linewidths=0.9,
            zorder=4,
        )
        for row in sub.itertuples():
            if row.n_frames < 100:
                ax.text(row.cutoff_A, row.R_like_fraction * 100.0 + 3.0, f"n={row.n_frames}", ha="center", va="bottom", fontsize=6.8, color=DATASET_COLORS[dataset])
        handles.append(line)
        labels.append(dataset)
    ax.set_xlabel(f"Near-attack cutoff C4-O1 ({ANGSTROM})")
    ax.set_ylabel("R-like fraction (%)")
    ax.set_ylim(0, 105)
    ax.invert_xaxis()
    top_legend(fig, handles, labels, ncol=2)
    fig.subplots_adjust(top=0.84, bottom=0.16, left=0.17, right=0.96)
    savefig(fig, outdir / "pub_amber_vs_gromacs_Rlike_fraction_by_cutoff_Arial.png")
    cut.to_csv(outdir / "amber_vs_gromacs_Rlike_fraction_by_cutoff.csv", index=False)
    return cut


def plot_rmsf(rmsf: pd.DataFrame | None, outdir: Path) -> pd.DataFrame:
    if rmsf is None or rmsf.empty:
        return pd.DataFrame()

    fig, ax = plt.subplots(figsize=SQUARE)
    rows = []
    handles = []
    labels = []
    for dataset, sub in rmsf.groupby("dataset", sort=False):
        sub = sub.sort_values("residue_index")
        color = DATASET_COLORS[dataset]
        line = ax.plot(sub["residue_index"], sub["rmsf_A"], color=color, lw=0.8, alpha=0.85, label=dataset)[0]
        handles.append(line)
        labels.append(dataset)
        max_row = sub.loc[sub["rmsf_A"].idxmax()]
        rows.append(
            {
                "dataset": dataset,
                "n_CA_residues": int(len(sub)),
                "median_RMSF_A": float(sub["rmsf_A"].median()),
                "mean_RMSF_A": float(sub["rmsf_A"].mean()),
                "max_RMSF_A": float(sub["rmsf_A"].max()),
                "max_RMSF_residue_index": int(max_row["residue_index"]),
            }
        )
    ax.set_xlabel("Residue index")
    ax.set_ylabel(f"C-alpha RMSF ({ANGSTROM})")
    ax.set_xlim(rmsf["residue_index"].min(), rmsf["residue_index"].max())
    top_legend(fig, handles, labels, ncol=2)
    fig.subplots_adjust(top=0.84, bottom=0.16, left=0.17, right=0.96)
    savefig(fig, outdir / "pub_amber_vs_gromacs_CA_RMSF_Arial.png")

    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "amber_vs_gromacs_CA_RMSF_summary.csv", index=False)
    return summary


def plot_single_rmsf(rmsf: pd.DataFrame | None, outdir: Path) -> None:
    if rmsf is None or rmsf.empty:
        return
    for dataset, sub in rmsf.groupby("dataset", sort=False):
        sub = sub.sort_values("residue_index")
        color = DATASET_COLORS[dataset]
        fig, ax = plt.subplots(figsize=SQUARE)
        ax.plot(sub["residue_index"], sub["rmsf_A"], color=color, lw=0.8, alpha=0.88)
        ax.set_xlabel("Residue index")
        ax.set_ylabel(f"C-alpha RMSF ({ANGSTROM})")
        ax.set_title(dataset, pad=4)
        fig.subplots_adjust(top=0.88, bottom=0.16, left=0.17, right=0.96)
        savefig(fig, outdir / f"pub_{safe_dataset_name(dataset)}_CA_RMSF_Arial.png")


def plot_pocket_pro_rs_delta(gmx_dir: Path, amber_dir: Path, outdir: Path, near_cutoff: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    gmx = load_gromacs_pocket_features(gmx_dir)
    amber = load_amber_pocket_features(outdir, amber_dir)
    if amber is None or amber.empty:
        return pd.DataFrame(), pd.DataFrame()
    data = pd.concat([gmx, amber], ignore_index=True, sort=False)
    near = data[data["C4_O1_A"] <= near_cutoff].copy()
    if near.empty:
        return pd.DataFrame(), pd.DataFrame()

    rows = []
    for dataset, dsub in near.groupby("dataset", sort=False):
        sampling = str(dsub["sampling"].iloc[0])
        for pro_label, psub in dsub.groupby("pro_rs_label"):
            for col, label in POCKET_FEATURES:
                rows.append(
                    {
                        "dataset": dataset,
                        "sampling": sampling,
                        "subset": f"C4O1_le_{near_cutoff:.1f}A",
                        "pro_rs_label": pro_label,
                        "pocket_feature": label,
                        "n_frames": int(len(psub)),
                        "median_distance_A": float(psub[col].median()),
                        "mean_distance_A": float(psub[col].mean()),
                        "q25_distance_A": float(psub[col].quantile(0.25)),
                        "q75_distance_A": float(psub[col].quantile(0.75)),
                    }
                )
    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "amber_vs_gromacs_proRS_pocket_distance_summary.csv", index=False)

    delta_rows = []
    for dataset, dsub in summary.groupby("dataset", sort=False):
        for feature, fsub in dsub.groupby("pocket_feature", sort=False):
            med = fsub.set_index("pro_rs_label")["median_distance_A"]
            n = fsub.set_index("pro_rs_label")["n_frames"]
            if "R-like" not in med.index or "S-like" not in med.index:
                continue
            delta_rows.append(
                {
                    "dataset": dataset,
                    "sampling": str(fsub["sampling"].iloc[0]),
                    "pocket_feature": feature,
                    "R_like_median_A": float(med["R-like"]),
                    "S_like_median_A": float(med["S-like"]),
                    "R_minus_S_median_A": float(med["R-like"] - med["S-like"]),
                    "R_like_n": int(n["R-like"]),
                    "S_like_n": int(n["S-like"]),
                }
            )
    delta = pd.DataFrame(delta_rows)
    delta.to_csv(outdir / "amber_vs_gromacs_proRS_pocket_delta.csv", index=False)

    features = [label for _, label in POCKET_FEATURES]
    x = np.arange(len(features))
    width = 0.34
    fig, ax = plt.subplots(figsize=(5.4, 5.4))
    offsets = {"GROMACS": -width / 2, "Amber pmemd": width / 2}
    handles = []
    labels = []
    for dataset in ["GROMACS", "Amber pmemd"]:
        sub = delta[delta["dataset"] == dataset].set_index("pocket_feature").reindex(features)
        vals = sub["R_minus_S_median_A"].to_numpy(float)
        bars = ax.bar(x + offsets[dataset], vals, width=width, color=DATASET_COLORS[dataset], alpha=0.78, label=dataset)
        handles.append(bars[0])
        labels.append(dataset)
    ax.axhline(0.0, color="0.25", lw=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(features, rotation=30, ha="right")
    ax.set_ylabel(f"Median distance difference, R-like - S-like ({ANGSTROM})")
    ax.text(
        0.02,
        0.98,
        f"C4-O1 {LE} {near_cutoff:.1f} {ANGSTROM}; negative means closer in R-like",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    top_legend(fig, handles, labels, ncol=2)
    fig.subplots_adjust(top=0.84, bottom=0.26, left=0.18, right=0.96)
    savefig(fig, outdir / "pub_amber_vs_gromacs_proRS_pocket_delta_Arial.png")
    return summary, delta


def pocket_feature_dataset(gmx_dir: Path, amber_dir: Path, outdir: Path, near_cutoff: float) -> pd.DataFrame:
    gmx = load_gromacs_pocket_features(gmx_dir)
    amber = load_amber_pocket_features(outdir, amber_dir)
    data = gmx if amber is None or amber.empty else pd.concat([gmx, amber], ignore_index=True, sort=False)
    data = data[data["C4_O1_A"] <= near_cutoff].copy()
    data = data.dropna(subset=[col for col, _ in POCKET_FEATURES] + ["pro_rs_label"])
    return data


def analyze_pocket_feature_classifier(gmx_dir: Path, amber_dir: Path, outdir: Path, near_cutoff: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    data = pocket_feature_dataset(gmx_dir, amber_dir, outdir, near_cutoff)
    if data.empty:
        return pd.DataFrame(), pd.DataFrame()

    feature_cols = [col for col, _ in POCKET_FEATURES]
    feature_names = [label for _, label in POCKET_FEATURES]
    metric_rows = []
    importance_rows = []
    for dataset_name, dsub in list(data.groupby("dataset", sort=False)) + [("Combined", data)]:
        dsub = dsub.dropna(subset=feature_cols + ["pro_rs_label"]).copy()
        y = (dsub["pro_rs_label"] == "R-like").astype(int).to_numpy()
        if len(np.unique(y)) < 2 or min(np.bincount(y)) < 5:
            continue
        x = dsub[feature_cols].to_numpy(float)
        n_splits = min(5, int(min(np.bincount(y))))
        model = RandomForestClassifier(
            n_estimators=500,
            random_state=20260710,
            class_weight="balanced",
            max_depth=4,
            min_samples_leaf=5,
        )
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=20260710)
        prob = cross_val_predict(model, x, y, cv=cv, method="predict_proba")[:, 1]
        pred = (prob >= 0.5).astype(int)
        metric_rows.append(
            {
                "dataset": dataset_name,
                "subset": f"C4O1_le_{near_cutoff:.1f}A",
                "n_frames": int(len(dsub)),
                "R_like_n": int(y.sum()),
                "S_like_n": int((1 - y).sum()),
                "roc_auc": float(roc_auc_score(y, prob)),
                "accuracy": float(accuracy_score(y, pred)),
                "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
                "model": "RandomForest on pocket contact features",
            }
        )
        model.fit(x, y)
        for name, importance in zip(feature_names, model.feature_importances_):
            importance_rows.append(
                {
                    "dataset": dataset_name,
                    "subset": f"C4O1_le_{near_cutoff:.1f}A",
                    "pocket_feature": name,
                    "importance": float(importance),
                }
            )

    metrics = pd.DataFrame(metric_rows)
    importances = pd.DataFrame(importance_rows)
    metrics.to_csv(outdir / "pocket_feature_RS_classifier_metrics.csv", index=False)
    importances.to_csv(outdir / "pocket_feature_RS_classifier_importance.csv", index=False)

    if not importances.empty:
        # Plot only real engines. The pooled Combined rows remain in CSV/README as a robustness check,
        # but should not look like a third physical trajectory in presentation figures.
        datasets = [x for x in ["GROMACS", "Amber pmemd"] if x in set(importances["dataset"])]
        x = np.arange(len(feature_names))
        width = 0.30
        offsets = np.linspace(-width, width, len(datasets))
        fig, ax = plt.subplots(figsize=(5.4, 5.4))
        colors = {"GROMACS": DATASET_COLORS["GROMACS"], "Amber pmemd": DATASET_COLORS["Amber pmemd"], "Combined": "#636363"}
        handles = []
        labels = []
        for dataset_name, offset in zip(datasets, offsets):
            sub = importances[importances["dataset"] == dataset_name].set_index("pocket_feature").reindex(feature_names)
            bars = ax.bar(x + offset, sub["importance"].to_numpy(float), width=width, color=colors[dataset_name], alpha=0.78, label=dataset_name)
            handles.append(bars[0])
            labels.append(dataset_name)
        ax.set_xticks(x)
        ax.set_xticklabels(feature_names, rotation=30, ha="right")
        ax.set_ylabel("Random-forest feature importance")
        ax.text(
            0.02,
            0.98,
            f"Separate-engine models, C4-O1 {LE} {near_cutoff:.1f} {ANGSTROM}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=10.0,
        )
        top_legend(fig, handles, labels, ncol=min(3, len(labels)))
        fig.subplots_adjust(top=0.84, bottom=0.26, left=0.16, right=0.96)
        savefig(fig, outdir / "pub_pocket_feature_RS_classifier_importance_Arial.png")
    return metrics, importances


def plot_triple_hist(pro_rs: pd.DataFrame, outdir: Path, near_cutoff: float) -> pd.DataFrame:
    rows = []
    fig, axes = plt.subplots(2, 1, figsize=(5.4, 5.4), sharex=True, sharey=True)
    bins = np.linspace(-8.5, 8.5, 86)
    for ax, (dataset, sub) in zip(axes, pro_rs.groupby("dataset", sort=False)):
        near = sub[sub["C4_O1_A"] <= near_cutoff]
        for label, color in [("R_like_by_product_reference", PRO_COLORS["R-like"]), ("S_like_by_product_reference", PRO_COLORS["S-like"])]:
            vals = near.loc[near["pseudo_RS_by_product_ref"] == label, TRIPLE_COL].to_numpy(float)
            if len(vals):
                weights = np.ones_like(vals) * 100.0 / len(near)
                ax.hist(vals, bins=bins, weights=weights, histtype="stepfilled", alpha=0.32, color=color, label=label.replace("_by_product_reference", "").replace("_", "-"))
        ax.axvline(0.0, color="0.35", ls="--", lw=0.9)
        ax.text(0.02, 0.88, dataset, transform=ax.transAxes, ha="left", va="top")
        total = len(near)
        for pro_label, pro_sub in near.groupby("pseudo_RS_by_product_ref"):
            rows.append(
                {
                    "dataset": dataset,
                    "subset": f"C4O1_le_{near_cutoff:.1f}A",
                    "pro_rs_label": pro_label.replace("_by_product_reference", "").replace("_", "-"),
                    "n_frames": int(len(pro_sub)),
                    "fraction": float(len(pro_sub) / total) if total else np.nan,
                    "median_triple_A3": float(pro_sub[TRIPLE_COL].median()),
                }
            )
    axes[-1].set_xlabel(r"Triple product V around C4 ($\AA^3$)")
    for ax in axes:
        ax.set_ylabel("Frequency (%)")
    handles, labels = axes[0].get_legend_handles_labels()
    top_legend(fig, handles, ["R-like", "S-like"], ncol=2)
    fig.subplots_adjust(top=0.88, bottom=0.13, left=0.16, right=0.96, hspace=0.18)
    savefig(fig, outdir / "pub_amber_vs_gromacs_triple_product_near_attack_Arial.png")
    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "amber_vs_gromacs_triple_product_near_attack_summary.csv", index=False)
    return summary


def load_group_pro_rs(gmx_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    group_frames = pd.read_csv(gmx_dir / "active_site_group_frame_assignments.csv")
    pro = pd.read_csv(gmx_dir.parent / "ul1_rotamer_face_analysis" / "product_calibrated_pseudo_RS" / "md_product_calibrated_pseudo_RS_timeseries.csv")
    pro["pro_rs_label"] = pro["pseudo_RS_by_product_ref"].map(
        {
            "R_like_by_product_reference": "R-like",
            "S_like_by_product_reference": "S-like",
        }
    )
    merged = group_frames.merge(pro[["frame", "C4_O1_A", TRIPLE_COL, "pro_rs_label"]], on="frame", suffixes=("", "_proRS"))
    return group_frames, merged


def plot_group_and_pro_rs_same_axes(gmx_dir: Path, outdir: Path, near_cutoff: float) -> None:
    _, merged = load_group_pro_rs(gmx_dir)
    if TRIPLE_COL not in merged.columns:
        return
    merged = merged.dropna(subset=[TRIPLE_COL, "C4_O1_A_proRS", "pro_rs_label", "group"])
    if merged.empty:
        return

    group_palette = {"Group 1": "#1b9e77", "Group 2": "#4c78a8", "Group 3": "#e45756"}
    fig, axes = plt.subplots(2, 1, figsize=(5.4, 6.4), sharex=True, sharey=True)

    handles = []
    labels = []
    for label in ["R-like", "S-like"]:
        sub = merged[merged["pro_rs_label"] == label]
        sc = axes[0].scatter(
            sub[TRIPLE_COL],
            sub["C4_O1_A_proRS"],
            s=8.0,
            alpha=0.24,
            linewidths=0,
            color=PRO_COLORS[label],
            label=label,
        )
        handles.append(sc)
        labels.append(label)
    axes[0].set_title("Same strided GROMACS frames colored by R/S", pad=5)
    axes[0].set_ylabel(f"C4-O1 ({ANGSTROM})")

    group_handles = []
    group_labels = []
    for group, sub in merged.sort_values("group_order").groupby("group", sort=False):
        sc = axes[1].scatter(
            sub[TRIPLE_COL],
            sub["C4_O1_A_proRS"],
            s=8.0,
            alpha=0.34,
            linewidths=0,
            color=group_palette.get(group, "0.5"),
            label=group,
        )
        group_handles.append(sc)
        group_labels.append(group)
    axes[1].set_title("Same frames colored by active-site Group 1/2/3", pad=5)
    axes[1].set_xlabel(r"Triple product V around C4 ($\AA^3$)")
    axes[1].set_ylabel(f"C4-O1 ({ANGSTROM})")

    for ax in axes:
        ax.axvline(0.0, color="0.35", lw=1.0, ls="--")
        ax.axhline(near_cutoff, color="0.35", lw=1.0, ls=":")
        ax.set_xlim(-8.6, 8.6)
        ax.set_ylim(2.5, 5.15)
    axes[0].legend(handles, labels, loc="upper right", frameon=False, markerscale=3.0, handletextpad=0.2)
    axes[1].legend(group_handles, group_labels, loc="upper right", frameon=False, markerscale=2.6, handletextpad=0.2)
    fig.subplots_adjust(top=0.92, bottom=0.11, left=0.16, right=0.96, hspace=0.28)
    savefig(fig, outdir / "pub_gromacs_group_vs_proRS_same_axes_Arial.png")


def summarize_groups(gmx_dir: Path, outdir: Path, near_cutoff: float) -> pd.DataFrame:
    _, merged = load_group_pro_rs(gmx_dir)
    rows = []
    metric_rows = []
    for subset_name, sub in [("all_strided_frames", merged), (f"C4O1_le_{near_cutoff:.1f}A_strided", merged[merged["C4_O1_A_proRS"] <= near_cutoff])]:
        if not sub.empty and sub["group"].nunique() > 1 and sub["pro_rs_label"].nunique() > 1:
            counts_by_group = pd.crosstab(sub["group"], sub["pro_rs_label"])
            majority_correct = int(counts_by_group.max(axis=1).sum())
            baseline_correct = int(sub["pro_rs_label"].value_counts().max())
            metric_rows.append(
                {
                    "subset": subset_name,
                    "n_frames": int(len(sub)),
                    "n_groups_present": int(sub["group"].nunique()),
                    "majority_baseline_accuracy": float(baseline_correct / len(sub)),
                    "group_majority_accuracy": float(majority_correct / len(sub)),
                    "normalized_mutual_information": float(normalized_mutual_info_score(sub["pro_rs_label"], sub["group"])),
                    "adjusted_mutual_information": float(adjusted_mutual_info_score(sub["pro_rs_label"], sub["group"])),
                }
            )
        for group, grp in sub.groupby("group"):
            counts = grp["pro_rs_label"].value_counts()
            rows.append(
                {
                    "subset": subset_name,
                    "group": group,
                    "n_frames": int(len(grp)),
                    "R_like_n": int(counts.get("R-like", 0)),
                    "S_like_n": int(counts.get("S-like", 0)),
                    "R_like_fraction": float(counts.get("R-like", 0) / len(grp)) if len(grp) else np.nan,
                    "median_C4_O1_A": float(grp["C4_O1_A_proRS"].median()),
                    "median_PHE336_C4_A": float(grp["PHE322_ring_centroid_to_C4_A"].median()),
                }
            )
    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "gromacs_active_site_group_vs_proRS_summary.csv", index=False)
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(outdir / "gromacs_active_site_group_vs_proRS_metrics.csv", index=False)

    fig, axes = plt.subplots(2, 1, figsize=(5.4, 5.4), sharex=True, sharey=True)
    groups = sorted(summary.loc[summary["subset"] == "all_strided_frames", "group"].unique())
    for ax, subset_name, title in [
        (axes[0], "all_strided_frames", "All stride frames"),
        (axes[1], f"C4O1_le_{near_cutoff:.1f}A_strided", f"C4-O1 {LE} {near_cutoff:.1f} {ANGSTROM} stride frames"),
    ]:
        sub = summary[summary["subset"] == subset_name].set_index("group")
        heights = [float(sub.loc[g, "R_like_fraction"] * 100.0) if g in sub.index else 0.0 for g in groups]
        ns = [int(sub.loc[g, "n_frames"]) if g in sub.index else 0 for g in groups]
        colors = ["#4c78a8" if n else "0.85" for n in ns]
        ax.bar(groups, heights, color=colors, alpha=0.78)
        for i, (height, n) in enumerate(zip(heights, ns)):
            ax.text(i, min(height + 4.0, 103), f"n={n}", ha="center", va="bottom", fontsize=8)
        ax.text(0.02, 0.88, title, transform=ax.transAxes, ha="left", va="top")
        ax.axhline(50, color="0.45", lw=0.8, ls=":")
        ax.set_ylim(0, 110)
        ax.set_ylabel("R-like fraction (%)")
    axes[1].set_xlabel("GROMACS active-site group")
    fig.subplots_adjust(top=0.96, bottom=0.13, left=0.17, right=0.96, hspace=0.18)
    savefig(fig, outdir / "pub_gromacs_active_site_group_Rlike_fraction_Arial.png")

    fig, axes = plt.subplots(2, 1, figsize=(5.4, 5.4), sharex=True, sharey=True)
    for ax, subset_name, title in [
        (axes[0], "all_strided_frames", "All stride frames"),
        (axes[1], f"C4O1_le_{near_cutoff:.1f}A_strided", f"C4-O1 {LE} {near_cutoff:.1f} {ANGSTROM} stride frames"),
    ]:
        sub = summary[summary["subset"] == subset_name].set_index("group").reindex(groups)
        r_vals = sub["R_like_n"].fillna(0).to_numpy(float)
        s_vals = sub["S_like_n"].fillna(0).to_numpy(float)
        totals = r_vals + s_vals
        r_pct = np.divide(r_vals, totals, out=np.zeros_like(r_vals), where=totals > 0) * 100.0
        s_pct = np.divide(s_vals, totals, out=np.zeros_like(s_vals), where=totals > 0) * 100.0
        ax.bar(groups, r_pct, color=PRO_COLORS["R-like"], alpha=0.72, label="R-like")
        ax.bar(groups, s_pct, bottom=r_pct, color=PRO_COLORS["S-like"], alpha=0.72, label="S-like")
        for i, n in enumerate(totals.astype(int)):
            ax.text(i, 103, f"n={n}", ha="center", va="bottom", fontsize=7.5)
        ax.set_ylim(0, 112)
        ax.set_ylabel("Group composition (%)")
        ax.text(0.02, 0.88, title, transform=ax.transAxes, ha="left", va="top")
    axes[1].set_xlabel("GROMACS active-site group")
    handles, labels = axes[0].get_legend_handles_labels()
    top_legend(fig, handles, labels, ncol=2)
    fig.subplots_adjust(top=0.86, bottom=0.13, left=0.17, right=0.96, hspace=0.18)
    savefig(fig, outdir / "pub_gromacs_group_vs_proRS_composition_Arial.png")
    return summary


def write_readme(
    outdir: Path,
    font_used: str,
    rmsd_summary: pd.DataFrame,
    rmsf_summary: pd.DataFrame,
    c4o1_summary: pd.DataFrame,
    threshold_summary: pd.DataFrame,
    triple_summary: pd.DataFrame,
    group_summary: pd.DataFrame,
    pocket_summary: pd.DataFrame,
    pocket_delta: pd.DataFrame,
    group_metrics: pd.DataFrame,
    classifier_metrics: pd.DataFrame,
    classifier_importance: pd.DataFrame,
) -> None:
    lines = [
        "# Amber vs GROMACS MD Comparison",
        "",
        f"Font requested/used by Matplotlib: `{font_used}`.",
        "",
        "## What This Folder Adds",
        "",
        "These figures place the completed GROMACS 200 ns trajectory and the Amber pmemd trajectory on the same axes. "
        "The goal is to show that both MD engines produced stable trajectories and to compare whether the product-calibrated R-like/S-like tendencies agree.",
        "",
        "## Figure Guide",
        "",
        "- `pub_amber_vs_gromacs_CA_RMSD_Arial.png`: global production stability; both engines should plateau rather than drift.",
        "- `pub_amber_vs_gromacs_CA_RMSF_Arial.png`: per-residue flexibility comparison after C-alpha fitting.",
        "- `pub_gromacs_CA_RMSF_Arial.png` and `pub_amber_pmemd_CA_RMSF_Arial.png`: single-engine per-residue RMSF views without smoothing.",
        "- `pub_amber_vs_gromacs_C4O1_histogram_Arial.png`: direct comparison of the reactive C4-O1 approach distribution.",
        "- `pub_gromacs_C4O1_histogram_Arial.png` and `pub_amber_pmemd_C4O1_histogram_Arial.png`: single-engine C4-O1 histograms without percentile/median/cutoff guide lines; numerical summaries remain in CSV/README.",
        "- `pub_amber_vs_gromacs_Rlike_fraction_by_cutoff_Arial.png`: R-like fraction as the C4-O1 cutoff is tightened; low-n tail points are marked separately.",
        "- `pub_gromacs_proRS_triple_product_vs_C4O1_Arial.png`, `pub_amber_pmemd_proRS_triple_product_vs_C4O1_Arial.png`, and `pub_amber_vs_gromacs_proRS_triple_product_vs_C4O1_Arial.png`: product-calibrated R/S scalar-triple-product scatter plots on the same axes.",
        "- `pub_gromacs_group_vs_proRS_same_axes_Arial.png`: the same GROMACS strided frames colored once by R/S and once by Group 1/2/3, showing whether unsupervised clusters align with reaction-facing labels.",
        "- `pub_gromacs_group_vs_proRS_composition_Arial.png`: whether unsupervised active-site groups cleanly separate R-like/S-like frames.",
        "- `pub_amber_vs_gromacs_proRS_pocket_delta_Arial.png`: which pocket contacts differ between R-like and S-like near-attack frames.",
        "- `pub_pocket_feature_RS_classifier_importance_Arial.png`: separate GROMACS and Amber contact-feature classifier importances; residues are treated as interpretable contact nodes rather than raw xyz coordinates. Pooled `Combined` rows remain in CSV as a robustness check, but are not plotted as a third trajectory.",
        "",
        "## Stability",
        "",
        rmsd_summary.to_string(index=False),
        "",
        "RMSD is shown as raw traces plus a running median. The running median is a sliding-window median used to suppress fast thermal noise and reveal plateau/drift behavior.",
        "",
        "C-alpha RMSF was calculated after C-alpha fitting using cpptraj `atomicfluct @CA byres`. It compares per-residue flexibility, not global drift.",
        "The RMSF plots show direct per-residue values without a running median.",
        "",
        rmsf_summary.to_string(index=False) if rmsf_summary is not None and not rmsf_summary.empty else "RMSF summary was not generated.",
        "",
        "## C4-O1 Approach",
        "",
        c4o1_summary.to_string(index=False),
        "",
        "## R-like/S-like Tendencies",
        "",
        "R-like/S-like labels come from the product-calibrated scalar triple product around C4. They are not literal products formed in classical MD.",
        "",
        threshold_summary[["dataset_label", "subset", "n_frames", "R_like_fraction", "S_like_fraction", "median_signed_triple"]].to_string(index=False),
        "",
        "Cutoffs below about 2.8 A have small frame counts, so the apparent endpoint swings should be treated as descriptive rather than statistically strong. "
        "In particular, the GROMACS 2.70 A point contains only 3 frames, so its lower R-like fraction is not evidence that the trajectory truly becomes nonselective at the shortest distances.",
        "",
        "Near-attack triple-product split:",
        "",
        triple_summary.to_string(index=False),
        "",
        "## Do The GROMACS Groups Explain R/S?",
        "",
        "The original active-site groups are useful geometry clusters, but they are not identical to R-like/S-like labels. "
        "They were generated by KMeans clustering on mechanistic active-site features rather than whole-protein RMSD: C4-O1 distance, Fe-C4/Fe-O1, Fe-donor distances, C4/O1 approach angles, one C3-C4-O1-C1 torsion, and minimum distances from PHE333/F347, HD2/H270, PHE322/F336, and GLN255/Q269 to C4/O1. "
        "Raw clusters were relabeled so the cluster with the shortest median C4-O1 distance is Group 1; the remaining groups are ordered by population. "
        "The table below checks whether each group is enriched in R-like or S-like frames.",
        "",
        "`Stride 10` means the active-site feature clustering used every 10th saved GROMACS trajectory frame. "
        "The production trajectory was saved every 5 ps/frame, so stride 10 corresponds to one analyzed frame every 50 ps. "
        "This is a postprocessing subsampling choice to reduce autocorrelation and file size. It is not the MD integration timestep and not a separate simulation.",
        "",
        group_summary.to_string(index=False),
        "",
        "Interpretation: Group 1/2/3 are not reliable pro-R/pro-S predictors. "
        "The same-axis scatter plot shows that Group labels are mainly controlled by active-site geometry, especially C4-O1 distance and pocket basin. "
        "Group 1 is the common near-attack basin and is R-like enriched, but that does not mean Group 1 equals pro-R. "
        "Group 3 is S-like in near-attack frames but has only 8 stride frames, so it should be treated as a rare-state hint rather than a stable stereochemical class. "
        "For stereoselectivity analysis, the product-calibrated R/S label should remain the primary axis; Group 1/2/3 should be treated as auxiliary pocket-state context.",
        "",
        "Quantitative group-vs-R/S diagnostics:",
        "",
        group_metrics.to_string(index=False) if group_metrics is not None and not group_metrics.empty else "Group-vs-R/S metrics were not generated.",
        "",
        "## Pocket Residue Differences Between R-like And S-like Frames",
        "",
        "The table below compares median pocket distances in the near-attack subset. "
        "Negative `R_minus_S_median_A` means that feature is closer in R-like frames; positive means it is closer in S-like frames. "
        "GROMACS pocket features come from the stride-10 active-site feature table; Amber pocket features were extracted from all 40000 frames with cpptraj `mindist`/`distance geom`.",
        "",
        pocket_delta.to_string(index=False) if pocket_delta is not None and not pocket_delta.empty else "Pocket delta summary was not generated.",
        "",
        "## GNN-like Pocket Feature Classifier",
        "",
        "A true deep GNN is possible, but for a single 200 ns trajectory it would be easy to overfit and hard to interpret. "
        "This first-pass `GNN-like` analysis uses residue/contact features as graph edges and trains an interpretable random-forest classifier to predict R-like versus S-like labels. "
        "This asks whether pocket contacts contain enough information to recover the R/S-like split.",
        "",
        "This random-forest analysis is exploratory. Each frame is one row; the label is R-like versus S-like; and the inputs are a small hand-selected set of pocket distances. "
        "Because MD frames are time-correlated and the current features are not a complete pocket graph, the classifier should not be used as proof of causality or as a final predictive model. "
        "It is most useful as a ranking/audit tool: which selected distances are repeatedly useful for separating the already-defined R/S labels?",
        "",
        "`PHE336 ring-C4` importance means the random forest often used that feature to split frames into purer R-like/S-like subsets, measured as mean decrease in Gini impurity across trees. "
        "It is a statistical classifier importance, not an interaction energy and not proof that PHE336 alone causes selectivity. "
        "It also depends on the current residue identity: if F336 is mutated to Arg, a literal aromatic-ring centroid feature no longer exists, so cross-mutant analysis should use mutation-compatible features.",
        "",
        "For mutant comparisons, the next more general pocket model should center every frame on UL1/UNL1 and describe the surrounding pocket with residue-type-aware but mutation-compatible descriptors: residue center/side-chain heavy-atom minimum distance to UL1 C4/O1/phenyl centroid, residue charge/polarity/aromatic class, contact occupancy, and local pocket point-cloud/PCA or UMAP coordinates. "
        "That design allows F336->Arg to remain comparable as a pocket position/contact node even though the aromatic-ring feature is replaced by side-chain centroid/guanidinium-contact features.",
        "",
        classifier_metrics.to_string(index=False) if classifier_metrics is not None and not classifier_metrics.empty else "Classifier metrics were not generated.",
        "",
        "Random-forest feature importances:",
        "",
        classifier_importance.to_string(index=False) if classifier_importance is not None and not classifier_importance.empty else "Classifier importances were not generated.",
        "",
    ]
    (outdir / "README_amber_vs_gromacs_comparison.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    gmx_dir = Path(args.gmx_summary_dir)
    amber_dir = Path(args.amber_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    font_used = setup_style(args.arial_font)

    rmsd_summary = plot_rmsd(load_rmsd(gmx_dir, amber_dir), outdir)
    rmsf = load_rmsf(outdir)
    rmsf_summary = plot_rmsf(rmsf, outdir)
    plot_single_rmsf(rmsf, outdir)
    c4o1 = load_c4o1(gmx_dir, amber_dir)
    c4o1_summary = plot_c4o1(c4o1, outdir, args.near_cutoff)
    plot_single_c4o1_histograms(c4o1, outdir, args.near_cutoff)
    threshold_summary = plot_pro_rs_thresholds(amber_dir, outdir)
    pro_rs = load_pro_rs(gmx_dir, amber_dir)
    plot_pro_rs_scatter_by_dataset(pro_rs, outdir, args.near_cutoff)
    triple_summary = plot_triple_hist(pro_rs, outdir, args.near_cutoff)
    group_summary = summarize_groups(gmx_dir, outdir, args.near_cutoff)
    plot_group_and_pro_rs_same_axes(gmx_dir, outdir, args.near_cutoff)
    pocket_summary, pocket_delta = plot_pocket_pro_rs_delta(gmx_dir, amber_dir, outdir, args.near_cutoff)
    group_metrics = pd.read_csv(outdir / "gromacs_active_site_group_vs_proRS_metrics.csv")
    classifier_metrics, classifier_importance = analyze_pocket_feature_classifier(gmx_dir, amber_dir, outdir, args.near_cutoff)
    write_readme(
        outdir,
        font_used,
        rmsd_summary,
        rmsf_summary,
        c4o1_summary,
        threshold_summary,
        triple_summary,
        group_summary,
        pocket_summary,
        pocket_delta,
        group_metrics,
        classifier_metrics,
        classifier_importance,
    )
    print(outdir)


if __name__ == "__main__":
    main()
