#!/usr/bin/env python3
"""Create publication-style MD summary plots for the TJ MCPB project."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


ANGSTROM = "\u00c5"
ALPHA = "\u03b1"
PHI = "\u03c6"
PI = "\u03c0"
DEG = "\u00b0"
LE = "\u2264"
GE = "\u2265"
PLUS_MINUS = "\u00b1"
SQUARE = (4.4, 4.4)
SQUARE_LARGE = (5.4, 5.4)
NEUTRAL_GRAY = "0.72"

PHE336_DIST_COL = "PHE322_ring_centroid_to_C4_A"
PHE336_ANGLE_COL = "C4_radical_plane_normal_vs_PHE322_ring_normal_angle_deg"

POCKET_COLS = [
    "PHE322_min_to_C4O1_A",
    "PHE333_min_to_C4O1_A",
    "HD2_256_min_to_C4O1_A",
    "GLN255_min_to_C4O1_A",
]
POCKET_LABELS = [
    "PHE336\n(GMX 322)",
    "PHE347\n(GMX 333)",
    "H270/HD2\n(GMX 256)",
    "Q269\n(GMX 255)",
]

PRO_RS_TS = Path("ul1_rotamer_face_analysis/product_calibrated_pseudo_RS/md_product_calibrated_pseudo_RS_timeseries.csv")
TRIPLE_COL = "signed_triple_O1_C5phenyl_C3chain_about_C4"
PRO_RS_COL = "pseudo_RS_by_product_ref"
PRO_RS_LABELS = {
    "R_like_by_product_reference": "R-like",
    "S_like_by_product_reference": "S-like",
}

DEPRECATED_OUTPUTS = [
    "pub_residue_contact_frequency_short_vs_long_Arial.png",
    "pub_residue_mean_min_distance_short_vs_long_Arial.png",
    "pub_PHE322_pi_geometry_distance_vs_angle_Arial.png",
    "pub_PHE322_pi_radical_distance_hist_by_phi_basin_Arial.png",
    "PHE322_pi_radical_summary_by_phi_basin.csv",
    "pub_cluster_phi360_vs_C4O1_Arial.png",
    "pub_cluster_pocket_distance_boxplots_Arial.png",
    "cluster_pocket_interaction_summary.csv",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-dir", required=True)
    parser.add_argument("--cv-csv", required=True)
    parser.add_argument("--cluster-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--short-cutoff", type=float, default=3.0)
    parser.add_argument("--pi-distance-cutoff", type=float, default=5.1)
    parser.add_argument("--arial-font", default="/mnt/c/Windows/Fonts/arial.ttf")
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


def savefig(fig: plt.Figure, outpath: Path) -> None:
    fig.savefig(outpath, dpi=360)
    plt.close(fig)


def add_top_legend(fig: plt.Figure, handles, labels, ncol: int = 2) -> None:
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


def archive_deprecated_outputs(outdir: Path) -> None:
    archive = outdir / "deprecated_old_short_long_and_C_labels"
    for name in DEPRECATED_OUTPUTS:
        for path in [outdir / name, (outdir / name).with_suffix(".pdf")]:
            if not path.exists():
                continue
            archive.mkdir(exist_ok=True)
            target = archive / path.name
            if target.exists():
                target.unlink()
            path.replace(target)


def phi360(series: pd.Series) -> pd.Series:
    return (pd.to_numeric(series, errors="coerce") + 360.0) % 360.0


def plot_rmsd(rmsd: pd.DataFrame, outdir: Path) -> dict[str, float]:
    t = rmsd["time_ps"].to_numpy(float) / 1000.0
    y = rmsd["rmsd_A"].to_numpy(float)
    window_frames = 401
    window_ns = window_frames * np.median(np.diff(t))
    roll = pd.Series(y).rolling(window_frames, center=True, min_periods=1).median().to_numpy()
    fig, ax = plt.subplots(figsize=SQUARE)
    raw_line = ax.plot(t, y, lw=0.22, color=NEUTRAL_GRAY, alpha=0.65, label="Raw")[0]
    med_line = ax.plot(t, roll, lw=1.15, color="#2b8cbe", label=f"{window_ns:.1f} ns running median")[0]
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(f"Protein C{ALPHA} RMSD ({ANGSTROM})")
    ax.set_xlim(t.min(), t.max())
    add_top_legend(fig, [raw_line, med_line], [raw_line.get_label(), med_line.get_label()], ncol=2)
    fig.subplots_adjust(top=0.84, bottom=0.16, left=0.17, right=0.96)
    savefig(fig, outdir / "pub_RMSD_CA_Arial.png")
    return {
        "rmsd_initial_A": float(y[0]),
        "rmsd_1ns_A": float(y[np.argmin(abs(t - 1.0))]),
        "rmsd_median_A": float(np.median(y)),
        "rmsd_last50ns_median_A": float(np.median(y[t >= (t.max() - 50.0)])),
        "rmsd_max_A": float(np.max(y)),
    }


def load_rmsd_or_previous_summary(analysis_dir: Path, outdir: Path) -> tuple[pd.DataFrame | None, dict[str, float]]:
    rmsd_csv = analysis_dir / "protein_CA_RMSD.csv"
    if rmsd_csv.exists():
        return pd.read_csv(rmsd_csv), {"rmsd_source_available": 1.0}

    rmsd_dat = outdir / "rmsd_source" / "gmx_CA_RMSD_pbcfit.dat"
    if rmsd_dat.exists():
        return load_cpptraj_rmsd_dat(rmsd_dat), {"rmsd_source_available": 1.0}

    rmsd_dat = outdir / "rmsd_source" / "gmx_CA_RMSD.dat"
    if rmsd_dat.exists():
        return load_cpptraj_rmsd_dat(rmsd_dat), {"rmsd_source_available": 1.0}

    previous = outdir / "README_publication_style_summary.md"
    stats = {
        "rmsd_source_available": 0.0,
        "rmsd_initial_A": float("nan"),
        "rmsd_1ns_A": float("nan"),
        "rmsd_median_A": float("nan"),
        "rmsd_last50ns_median_A": float("nan"),
        "rmsd_max_A": float("nan"),
    }
    if not previous.exists():
        return None, stats

    text = previous.read_text(encoding="utf-8", errors="ignore")
    patterns = {
        "rmsd_initial_A": r"Initial.*?RMSD:\s*([0-9.]+)",
        "rmsd_1ns_A": r"RMSD at 1 ns:\s*([0-9.]+)",
        "rmsd_median_A": r"Median.*?RMSD:\s*([0-9.]+)",
        "rmsd_last50ns_median_A": r"Last 50 ns median.*?RMSD:\s*([0-9.]+)",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, text)
        if match:
            stats[key] = float(match.group(1))
    return None, stats


def load_cpptraj_rmsd_dat(path: Path) -> pd.DataFrame:
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
            rows.append((frame, (frame - 1) * 5.0, float(parts[1])))
    return pd.DataFrame(rows, columns=["frame", "time_ps", "rmsd_A"])


def load_c4o1_distance(analysis_dir: Path, cv: pd.DataFrame) -> pd.DataFrame:
    dist_csv = analysis_dir / "UL1_C4_O1_distance_timeseries.csv"
    if dist_csv.exists():
        return pd.read_csv(dist_csv)
    return pd.DataFrame(
        {
            "time_ps": cv["time_ps"],
            "C4_O1_distance_A": cv["C4_O1_A"],
        }
    )


def load_pro_rs(analysis_dir: Path) -> pd.DataFrame | None:
    path = analysis_dir / PRO_RS_TS
    if not path.exists():
        return None
    df = pd.read_csv(path)
    df["pro_rs_label"] = df[PRO_RS_COL].map(PRO_RS_LABELS).fillna(df[PRO_RS_COL])
    return df


def plot_c4o1_hist(dist: pd.DataFrame, outdir: Path, short_cutoff: float) -> dict[str, float]:
    y = dist["C4_O1_distance_A"].to_numpy(float)
    q10 = float(np.quantile(y, 0.10))
    med = float(np.median(y))
    fig, ax = plt.subplots(figsize=SQUARE)
    weights = np.ones_like(y) * 100.0 / len(y)
    ax.hist(y, bins=85, weights=weights, color="#8ecae6", edgecolor="#4a90b5", linewidth=0.35)
    q_line = ax.axvline(q10, color="#d7301f", ls="--", lw=1.0, label=f"10th percentile {q10:.2f} {ANGSTROM}")
    m_line = ax.axvline(med, color="0.20", ls="-", lw=1.0, label=f"Median {med:.2f} {ANGSTROM}")
    c_line = ax.axvline(short_cutoff, color="#2b8cbe", ls=":", lw=1.1, label=f"Near-attack {short_cutoff:.1f} {ANGSTROM}")
    ax.set_xlabel(f"UL1 C4-O1 distance ({ANGSTROM})")
    ax.set_ylabel("Frequency (%)")
    add_top_legend(fig, [q_line, m_line, c_line], [q_line.get_label(), m_line.get_label(), c_line.get_label()], ncol=1)
    fig.subplots_adjust(top=0.70, bottom=0.16, left=0.17, right=0.96)
    savefig(fig, outdir / "pub_C4_O1_distance_histogram_Arial.png")
    return {
        "C4O1_q10_A": q10,
        "C4O1_median_A": med,
        "C4O1_fraction_le_short_cutoff": float((y <= short_cutoff).mean()),
    }


def plot_phi360_scatter(cv: pd.DataFrame, selected_csv: Path, outdir: Path, short_cutoff: float) -> None:
    df = cv.copy()
    df["phi360"] = phi360(df["dihedral_C1_O1_C4_C5_deg"])
    df["absH"] = df["dihedral_C1_O1_C4_H02_deg"].abs()
    short = df[df["C4_O1_A"] <= short_cutoff]
    fig, ax = plt.subplots(figsize=SQUARE)
    all_sc = ax.scatter(df["phi360"], df["absH"], s=2.1, c="0.25", alpha=0.035, linewidths=0, marker="o", label="All frames")
    short_sc = ax.scatter(
        short["phi360"],
        short["absH"],
        s=3.7,
        c="#1b9e77",
        alpha=0.24,
        linewidths=0,
        marker="o",
        label=f"C4-O1 {LE} {short_cutoff:.1f} {ANGSTROM}",
    )
    handles = [all_sc, short_sc]
    labels = ["All frames", f"C4-O1 {LE} {short_cutoff:.1f} {ANGSTROM}"]
    if selected_csv.exists():
        selected = pd.read_csv(selected_csv)
        rep_sc = ax.scatter(
            selected["phi_C1O1C4C5_360_deg"],
            selected["abs_C1O1C4H02_deg"],
            s=24,
            c="#d95f02",
            edgecolors="white",
            linewidths=0.45,
            zorder=5,
            label="Representative frames",
        )
        handles.append(rep_sc)
        labels.append("Representative frames")
    ax.axvline(180.0, color="0.65", lw=0.75, ls="--")
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 180)
    ax.set_xticks(np.arange(0, 361, 60))
    ax.set_yticks(np.arange(0, 181, 30))
    ax.set_xlabel(f"Ordered torsion C1-O1-C4-C5(Ph), {PHI}0-360 ({DEG})")
    ax.set_ylabel(f"|Ordered torsion C1-O1-C4-H02| ({DEG})")
    add_top_legend(fig, handles, labels, ncol=1)
    fig.subplots_adjust(top=0.73, bottom=0.16, left=0.16, right=0.96)
    savefig(fig, outdir / "pub_phi360_prochirality_scatter_Arial.png")


def plot_phe336_pi_hist(cv: pd.DataFrame, outdir: Path, short_cutoff: float, pi_cutoff: float) -> pd.DataFrame:
    df = cv.copy()
    df["phi360"] = phi360(df["dihedral_C1_O1_C4_C5_deg"])
    short = df[df["C4_O1_A"] <= short_cutoff].copy()
    basin_low = "phi0-360 < 180 deg"
    basin_high = "phi0-360 >= 180 deg"
    short["basin"] = np.where(short["phi360"] < 180.0, basin_low, basin_high)
    colors = {basin_low: "#4c78a8", basin_high: "#f58518"}
    fig, axes = plt.subplots(2, 1, figsize=SQUARE_LARGE, sharex=True)
    summary_rows = []
    for ax, basin in zip(axes, [basin_low, basin_high]):
        sub = short[short["basin"] == basin]
        d = sub[PHE336_DIST_COL].to_numpy(float)
        weights = np.ones_like(d) * 100.0 / len(d) if len(d) else np.array([])
        ax.hist(d, bins=np.linspace(3.0, 8.0, 51), weights=weights, color="#9ecae1", edgecolor=colors[basin], linewidth=0.45)
        ax.axvline(pi_cutoff, color="0.20", lw=0.9, ls="--")
        pct = float((d <= pi_cutoff).mean() * 100.0) if len(d) else float("nan")
        mean = float(np.mean(d)) if len(d) else float("nan")
        sd = float(np.std(d, ddof=1)) if len(d) > 1 else float("nan")
        ax.text(
            0.97,
            0.88,
            f"{pct:.1f}% within\n{pi_cutoff:.1f} {ANGSTROM}\n"
            f"d = {mean:.2f} {PLUS_MINUS} {sd:.2f} {ANGSTROM}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8,
        )
        ax.set_title(basin.replace("deg", DEG))
        ax.set_xlabel(f"{PI}(PHE336)-C4 distance ({ANGSTROM})")
        summary_rows.append(
            {
                "basin": basin,
                "n_frames": len(d),
                "mean_distance_A": mean,
                "sd_distance_A": sd,
                "percent_within_pi_cutoff": pct,
                "median_angle_deg": float(sub[PHE336_ANGLE_COL].median()) if len(sub) else float("nan"),
            }
        )
    axes[0].set_ylabel("Frequency (%)")
    axes[1].set_ylabel("Frequency (%)")
    fig.subplots_adjust(top=0.96, bottom=0.13, left=0.16, right=0.97, hspace=0.28)
    savefig(fig, outdir / "pub_PHE336_pi_radical_distance_hist_by_phi_basin_Arial.png")
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(outdir / "PHE336_pi_radical_summary_by_phi_basin.csv", index=False)
    return summary


def plot_phe336_geometry_scatter(cv: pd.DataFrame, outdir: Path, short_cutoff: float, pi_cutoff: float) -> None:
    short = cv[cv["C4_O1_A"] <= short_cutoff].copy()
    fig, ax = plt.subplots(figsize=SQUARE)
    sc = ax.scatter(
        short[PHE336_DIST_COL],
        short[PHE336_ANGLE_COL],
        c=short["C4_O1_A"],
        s=5.0,
        cmap="viridis_r",
        alpha=0.55,
        linewidths=0,
    )
    ax.axvline(pi_cutoff, color="0.30", lw=0.85, ls="--")
    ax.axhline(45.0, color="0.50", lw=0.75, ls=":")
    ax.set_xlabel(f"{PI}(PHE336)-C4 distance ({ANGSTROM})")
    ax.set_ylabel(f"C4 plane / PHE336 ring angle ({DEG})")
    ax.set_ylim(0, 90)
    ax.text(0.02, 0.98, f"Dashed line: {PI} cutoff {pi_cutoff:.1f} {ANGSTROM}", transform=ax.transAxes, ha="left", va="top", fontsize=8)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label(f"C4-O1 distance ({ANGSTROM})")
    fig.subplots_adjust(top=0.96, bottom=0.15, left=0.16, right=0.84)
    savefig(fig, outdir / "pub_PHE336_pi_geometry_distance_vs_angle_Arial.png")


def merge_assignments_with_cv(assign: pd.DataFrame, cv: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "frame",
        "C4_O1_A",
        "dihedral_C1_O1_C4_C5_deg",
        "dihedral_C1_O1_C4_H02_deg",
        "pseudo_CIP_O1_C5_C3_H02_chiral_volume_A3",
        PHE336_DIST_COL,
        PHE336_ANGLE_COL,
    ]
    merged = assign.merge(cv[keep], on="frame", how="left", suffixes=("", "_cv"))
    merged["phi360"] = phi360(merged["dihedral_C1_O1_C4_C5_deg"])
    merged["absH"] = merged["dihedral_C1_O1_C4_H02_deg"].abs()
    return merged


def assign_group_labels(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for cid, sub in data.groupby("cluster_id"):
        rows.append(
            {
                "raw_cluster_id": int(cid),
                "n_frames": int(len(sub)),
                "population_fraction": float(len(sub) / len(data)),
                "median_C4_O1_A": float(sub["C4_O1_A"].median()),
                "median_phi360_deg": float(sub["phi360"].median()),
                "median_PHE336_centroid_C4_A": float(sub[PHE336_DIST_COL].median()),
                "median_PHE336_angle_deg": float(sub[PHE336_ANGLE_COL].median()),
                "median_PHE336_min_to_C4O1_A": float(sub["PHE322_min_to_C4O1_A"].median()),
                "median_PHE347_min_to_C4O1_A": float(sub["PHE333_min_to_C4O1_A"].median()),
                "median_H270_HD2_min_to_C4O1_A": float(sub["HD2_256_min_to_C4O1_A"].median()),
                "median_Q269_min_to_C4O1_A": float(sub["GLN255_min_to_C4O1_A"].median()),
            }
        )
    summary = pd.DataFrame(rows)
    nearest = int(summary.sort_values("median_C4_O1_A").iloc[0]["raw_cluster_id"])
    remaining = summary[summary["raw_cluster_id"] != nearest].sort_values("n_frames", ascending=False)
    ordered = [nearest] + [int(x) for x in remaining["raw_cluster_id"].tolist()]
    mapping = {cid: f"Group {idx + 1}" for idx, cid in enumerate(ordered)}
    labeled = data.copy()
    labeled["group"] = labeled["cluster_id"].map(mapping)
    labeled["group_order"] = labeled["cluster_id"].map({cid: i for i, cid in enumerate(ordered)})
    summary["group"] = summary["raw_cluster_id"].map(mapping)
    summary["group_order"] = summary["raw_cluster_id"].map({cid: i for i, cid in enumerate(ordered)})
    summary = summary.sort_values("group_order").drop(columns=["group_order"])
    return labeled, summary


def plot_pro_rs_triple(pro_rs: pd.DataFrame, outdir: Path, short_cutoff: float) -> pd.DataFrame:
    fig, ax = plt.subplots(figsize=SQUARE)
    colors = {"R-like": "#4c78a8", "S-like": "#e45756"}
    handles = []
    labels = []
    for label, sub in pro_rs.groupby("pro_rs_label", sort=False):
        sc = ax.scatter(
            sub[TRIPLE_COL],
            sub["C4_O1_A"],
            s=2.2,
            alpha=0.14,
            linewidths=0,
            color=colors.get(label, "0.45"),
            label=label,
        )
        handles.append(sc)
        labels.append(label)
    ax.axvline(0.0, color="0.35", lw=0.85, ls="--")
    ax.axhline(short_cutoff, color="0.35", lw=0.85, ls=":")
    ax.set_xlabel(r"Triple product V around C4 ($\AA^3$)")
    ax.set_ylabel(f"C4-O1 distance ({ANGSTROM})")
    add_top_legend(fig, handles, labels, ncol=2)
    fig.subplots_adjust(top=0.84, bottom=0.16, left=0.17, right=0.96)
    savefig(fig, outdir / "pub_proRS_triple_product_vs_C4O1_Arial.png")

    rows = []
    for subset_name, sub in [
        ("all_frames", pro_rs),
        (f"C4O1_le_{short_cutoff:.1f}A", pro_rs[pro_rs["C4_O1_A"] <= short_cutoff]),
    ]:
        for label, grp in sub.groupby("pro_rs_label"):
            rows.append(
                {
                    "subset": subset_name,
                    "pro_rs_label": label,
                    "n_frames": int(len(grp)),
                    "fraction_in_subset": float(len(grp) / len(sub)) if len(sub) else float("nan"),
                    "median_C4_O1_A": float(grp["C4_O1_A"].median()),
                    "median_triple_A3": float(grp[TRIPLE_COL].median()),
                }
            )
    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "proRS_triple_product_summary.csv", index=False)
    return summary


def merge_group_and_pro_rs(group_frames: pd.DataFrame, pro_rs: pd.DataFrame) -> pd.DataFrame:
    keep = ["frame", "time_ns", "C4_O1_A", TRIPLE_COL, "pro_rs_label"]
    return group_frames.merge(pro_rs[keep], on="frame", how="inner", suffixes=("", "_proRS"))


def plot_pro_rs_pocket(group_frames: pd.DataFrame, pro_rs: pd.DataFrame, outdir: Path, short_cutoff: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = merge_group_and_pro_rs(group_frames, pro_rs)
    near = merged[merged["C4_O1_A_proRS"] <= short_cutoff].copy()
    if near.empty:
        return pd.DataFrame(), pd.DataFrame()

    colors = {"R-like": "#4c78a8", "S-like": "#e45756"}
    fig, ax = plt.subplots(figsize=SQUARE_LARGE)
    positions = []
    values = []
    box_colors = []
    tick_positions = []
    pos = 1.0
    for label, col in zip(POCKET_LABELS, POCKET_COLS):
        group_positions = []
        for pro_label in ["R-like", "S-like"]:
            vals = near.loc[near["pro_rs_label"] == pro_label, col].to_numpy(float)
            values.append(vals)
            positions.append(pos)
            box_colors.append(colors[pro_label])
            group_positions.append(pos)
            pos += 0.92
        tick_positions.append(float(np.mean(group_positions)))
        pos += 0.55
    bp = ax.boxplot(values, positions=positions, widths=0.48, patch_artist=True, showfliers=False)
    for patch, color in zip(bp["boxes"], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.48)
        patch.set_edgecolor("0.25")
    for element in ["whiskers", "caps", "medians"]:
        for artist in bp[element]:
            artist.set_color("0.25")
            artist.set_linewidth(0.8)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(POCKET_LABELS, rotation=0)
    ax.set_ylabel(f"Minimum distance to C4/O1 ({ANGSTROM})")
    legend_handles = [plt.Line2D([0], [0], color=colors[x], lw=6, alpha=0.55, label=x) for x in ["R-like", "S-like"]]
    add_top_legend(fig, legend_handles, ["R-like", "S-like"], ncol=2)
    fig.subplots_adjust(top=0.86, bottom=0.19, left=0.15, right=0.97)
    savefig(fig, outdir / "pub_proRS_pocket_distance_boxplots_Arial.png")

    rows = []
    for pro_label, sub in near.groupby("pro_rs_label"):
        for label, col in zip(POCKET_LABELS, POCKET_COLS):
            rows.append(
                {
                    "subset": f"C4O1_le_{short_cutoff:.1f}A",
                    "pro_rs_label": pro_label,
                    "pocket_feature": label.replace("\n", " "),
                    "n_frames": int(len(sub)),
                    "median_distance_A": float(sub[col].median()),
                    "mean_distance_A": float(sub[col].mean()),
                    "q25_distance_A": float(sub[col].quantile(0.25)),
                    "q75_distance_A": float(sub[col].quantile(0.75)),
                }
            )
    pro_summary = pd.DataFrame(rows)
    pro_summary.to_csv(outdir / "proRS_pocket_distance_summary.csv", index=False)

    crosstab = pd.crosstab(near["group"], near["pro_rs_label"], normalize="index")
    counts = pd.crosstab(near["group"], near["pro_rs_label"])
    combined = crosstab.add_suffix("_fraction").join(counts.add_suffix("_n"))
    combined.insert(0, "n_frames", counts.sum(axis=1))
    combined.to_csv(outdir / "proRS_by_active_site_group_summary.csv")
    return pro_summary, combined.reset_index()


def plot_pocket_state_enrichment(group_frames: pd.DataFrame, pro_rs: pd.DataFrame, outdir: Path, short_cutoff: float) -> pd.DataFrame:
    merged = merge_group_and_pro_rs(group_frames, pro_rs)
    near = merged[merged["C4_O1_A_proRS"] <= short_cutoff].copy()
    features = POCKET_COLS + [PHE336_DIST_COL, PHE336_ANGLE_COL]
    near = near.dropna(subset=features + ["pro_rs_label"])
    if len(near) < 9:
        return pd.DataFrame()

    x = StandardScaler().fit_transform(near[features].to_numpy(float))
    n_clusters = 3
    raw_state = KMeans(n_clusters=n_clusters, random_state=7, n_init=50).fit_predict(x)
    near["raw_pocket_state"] = raw_state

    state_rows = []
    for raw, sub in near.groupby("raw_pocket_state"):
        r_frac = float((sub["pro_rs_label"] == "R-like").mean())
        state_rows.append(
            {
                "raw_pocket_state": int(raw),
                "n_frames": int(len(sub)),
                "R_like_fraction": r_frac,
                "median_PHE336_min_to_C4O1_A": float(sub["PHE322_min_to_C4O1_A"].median()),
                "median_PHE347_min_to_C4O1_A": float(sub["PHE333_min_to_C4O1_A"].median()),
                "median_H270_HD2_min_to_C4O1_A": float(sub["HD2_256_min_to_C4O1_A"].median()),
                "median_Q269_min_to_C4O1_A": float(sub["GLN255_min_to_C4O1_A"].median()),
                "median_PHE336_centroid_C4_A": float(sub[PHE336_DIST_COL].median()),
            }
        )
    state_summary = pd.DataFrame(state_rows).sort_values(["R_like_fraction", "n_frames"], ascending=[False, False])
    mapping = {raw: f"Pocket state {idx + 1}" for idx, raw in enumerate(state_summary["raw_pocket_state"].tolist())}
    near["pocket_state"] = near["raw_pocket_state"].map(mapping)
    state_summary["pocket_state"] = state_summary["raw_pocket_state"].map(mapping)

    counts = pd.crosstab(near["pocket_state"], near["pro_rs_label"])
    fractions = counts.div(counts.sum(axis=1), axis=0)
    states = state_summary["pocket_state"].tolist()
    fig, ax = plt.subplots(figsize=SQUARE)
    bottom = np.zeros(len(states))
    colors = {"R-like": "#4c78a8", "S-like": "#e45756"}
    handles = []
    labels = []
    for label in ["R-like", "S-like"]:
        vals = fractions.reindex(states).get(label, pd.Series(0.0, index=states)).to_numpy(float) * 100.0
        bars = ax.bar(states, vals, bottom=bottom, color=colors[label], alpha=0.72, label=label)
        bottom += vals
        handles.append(bars[0])
        labels.append(label)
    for idx, state in enumerate(states):
        n = int(counts.reindex(states).sum(axis=1).loc[state])
        ax.text(idx, 102, f"n={n}", ha="center", va="bottom", fontsize=7.5)
    ax.set_ylim(0, 112)
    ax.set_ylabel("Fraction within pocket state (%)")
    add_top_legend(fig, handles, labels, ncol=2)
    fig.subplots_adjust(top=0.84, bottom=0.18, left=0.16, right=0.96)
    savefig(fig, outdir / "pub_pocket_state_proRS_enrichment_Arial.png")

    near.to_csv(outdir / "pocket_state_near_attack_frame_assignments.csv", index=False)
    state_summary.to_csv(outdir / "pocket_state_proRS_summary.csv", index=False)
    return state_summary


def plot_cluster_geometry(cluster_dir: Path, cv: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    all_assign = pd.read_csv(cluster_dir / "all_stride_frame_assignments_k3.csv")
    all_m = merge_assignments_with_cv(all_assign, cv)
    all_m, summary = assign_group_labels(all_m)

    palette = {"Group 1": "#1b9e77", "Group 2": "#4c78a8", "Group 3": "#e45756"}
    fig, axes = plt.subplots(2, 1, figsize=SQUARE_LARGE)
    handles = []
    labels = []
    for group, sub in all_m.sort_values("group_order").groupby("group", sort=False):
        sc = axes[0].scatter(sub["phi360"], sub["C4_O1_A"], s=4.5, alpha=0.42, linewidths=0, color=palette.get(group), label=group)
        handles.append(sc)
        labels.append(group)
    axes[0].axvline(180.0, color="0.68", lw=0.75, ls="--")
    axes[0].axhline(3.0, color="0.30", lw=0.75, ls=":")
    axes[0].set_xlim(0, 360)
    axes[0].set_xlabel(f"Ordered torsion {PHI}0-360 ({DEG})")
    axes[0].set_ylabel(f"C4-O1 distance ({ANGSTROM})")
    axes[0].set_title("Active-site geometry groups", pad=3)

    for group, sub in all_m.sort_values("group_order").groupby("group", sort=False):
        axes[1].scatter(sub["C4_O1_A"], sub[PHE336_DIST_COL], s=4.5, alpha=0.42, linewidths=0, color=palette.get(group), label=group)
    axes[1].axvline(3.0, color="0.30", lw=0.75, ls=":")
    axes[1].axhline(5.1, color="0.30", lw=0.75, ls="--")
    axes[1].set_xlabel(f"C4-O1 distance ({ANGSTROM})")
    axes[1].set_ylabel(f"{PI}(PHE336)-C4 distance ({ANGSTROM})")
    axes[1].set_title("PHE336 contact by group")
    add_top_legend(fig, handles, labels, ncol=3)
    fig.subplots_adjust(top=0.90, bottom=0.10, left=0.15, right=0.96, hspace=0.33)
    savefig(fig, outdir / "pub_active_site_group_geometry_map_Arial.png")

    fig, ax = plt.subplots(figsize=SQUARE)
    groups = list(summary["group"])
    x = np.arange(len(groups))
    width = 0.36
    ax.bar(x - width / 2, summary["population_fraction"] * 100.0, width, color="#4c78a8", label="Population")
    ax2 = ax.twinx()
    ax2.plot(x + width / 2, summary["median_C4_O1_A"], color="#d7301f", marker="o", lw=1.2, label="Median C4-O1")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{g}\n(raw {cid})" for g, cid in zip(groups, summary["raw_cluster_id"])])
    ax.set_ylabel("Population (%)")
    ax2.set_ylabel(f"Median C4-O1 ({ANGSTROM})")
    ax.set_title("Active-site group summary")
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    add_top_legend(fig, lines + lines2, labels + labels2, ncol=1)
    fig.subplots_adjust(top=0.82, bottom=0.17, left=0.16, right=0.84)
    savefig(fig, outdir / "pub_active_site_group_population_summary_Arial.png")

    fig, ax = plt.subplots(figsize=SQUARE_LARGE)
    positions = []
    values = []
    colors = []
    tick_positions = []
    pos = 1.0
    for label, col in zip(POCKET_LABELS, POCKET_COLS):
        group_positions = []
        for group in groups:
            values.append(all_m.loc[all_m["group"] == group, col].to_numpy(float))
            positions.append(pos)
            colors.append(palette.get(group, "0.5"))
            group_positions.append(pos)
            pos += 1.0
        tick_positions.append(np.mean(group_positions))
        pos += 0.70
    bp = ax.boxplot(values, positions=positions, widths=0.55, patch_artist=True, showfliers=False)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.45)
        patch.set_edgecolor("0.25")
    for element in ["whiskers", "caps", "medians"]:
        for artist in bp[element]:
            artist.set_color("0.25")
            artist.set_linewidth(0.8)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(POCKET_LABELS, rotation=0)
    ax.set_ylabel(f"Minimum distance to C4/O1 ({ANGSTROM})")
    ax.set_title("Pocket residue distances by active-site group", pad=6)
    legend_handles = [plt.Line2D([0], [0], color=palette[g], lw=6, alpha=0.55, label=g) for g in groups]
    add_top_legend(fig, legend_handles, [h.get_label() for h in legend_handles], ncol=len(groups))
    fig.subplots_adjust(top=0.85, bottom=0.20, left=0.15, right=0.96)
    savefig(fig, outdir / "pub_active_site_group_pocket_distance_boxplots_Arial.png")

    summary.to_csv(outdir / "active_site_group_summary.csv", index=False)
    all_m.to_csv(outdir / "active_site_group_frame_assignments.csv", index=False)
    return summary


def write_readme(
    outdir: Path,
    font_used: str,
    summaries: dict[str, float],
    phe_summary: pd.DataFrame,
    group_summary: pd.DataFrame,
    pro_rs_summary: pd.DataFrame | None,
    pocket_state_summary: pd.DataFrame | None,
    short_cutoff: float,
) -> None:
    rmsd_source_note = (
        "RMSD source: regenerated with cpptraj from the PBC-fitted protein/MCPB GROMACS trajectory using the C-alpha atom mask."
        if summaries.get("rmsd_source_available", 0.0) == 1.0
        else "Note: the RMSD source CSV was not present in this folder, so the existing RMSD figure was retained and the numeric values were recovered from the previous README."
    )
    lines = [
        "# Publication-Style MD Summary",
        "",
        f"Font requested/used by Matplotlib: `{font_used}`.",
        "",
        "## Stability",
        "",
        f"- Initial C-alpha RMSD: {summaries['rmsd_initial_A']:.3f} A",
        f"- C-alpha RMSD at 1 ns: {summaries['rmsd_1ns_A']:.3f} A",
        f"- Median C-alpha RMSD: {summaries['rmsd_median_A']:.3f} A",
        f"- Last 50 ns median C-alpha RMSD: {summaries['rmsd_last50ns_median_A']:.3f} A",
        "",
        "The rapid rise from 0 to about 1 A is expected when production is compared with the first production frame. "
        "The key question is whether RMSD plateaus rather than continuously drifting.",
        "",
        "Running median explanation: the blue RMSD curve is a sliding-window median, not a separate trajectory. "
        "For each time point, the median RMSD of the neighboring frames is plotted to suppress fast thermal noise "
        "and emphasize slow drift or plateau behavior. It is more resistant to isolated spikes than a running mean.",
        "",
        rmsd_source_note,
        "",
        "## C4-O1 Distance",
        "",
        f"- 10th percentile: {summaries['C4O1_q10_A']:.3f} A",
        f"- Median: {summaries['C4O1_median_A']:.3f} A",
        f"- Fraction with C4-O1 <= {short_cutoff:.1f} A: {summaries['C4O1_fraction_le_short_cutoff']*100:.2f}%",
        "",
        "The red dashed line in the histogram is the 10th percentile, not the median.",
        "",
        "## Dihedral Direction",
        "",
        "An unordered plane-plane angle is unsigned and ranges from 0 to 180 degrees. "
        "A molecular torsion/dihedral is signed because the four ordered atoms define an orientation. "
        "Here phi0-360 is computed as `(signed_dihedral_C1_O1_C4_C5 + 360) % 360`. "
        "It is a basin descriptor, not an absolute R/S assignment by itself.",
        "",
        "## Triple-Product Pro-R/Pro-S Descriptor",
        "",
        "The scalar triple product used here is `V = det[(O1-C4), (C5-C4), (C3-C4)]`. "
        "Its sign captures which side of the local C4-centered geometry the substrate occupies. "
        "The sign has been used only after product-reference calibration, so the output labels are written as "
        "`R-like` and `S-like` instead of claiming that an unreacted MD frame has already formed an R or S product.",
        "",
        "Triple-product summary:",
        "",
        pro_rs_summary.to_string(index=False) if pro_rs_summary is not None and not pro_rs_summary.empty else "Not generated.",
        "",
        "## PHE336 pi-radical Contact",
        "",
        phe_summary.to_string(index=False),
        "",
        "Residue numbering note: GROMACS/Amber residue PHE322 corresponds to original-structure F336. "
        "The same +14 offset is used for PHE333/F347, HD2 256/H270, and GLN255/Q269 in the plot labels.",
        "",
        "## Active-Site Groups",
        "",
        "Old labels C0/C1/C2 have been replaced by Group 1/2/3 to avoid confusion with atom names C1/C2. "
        "Groups are based on the existing active-site feature clustering, not global protein C-alpha RMSD.",
        "",
        group_summary.to_string(index=False),
        "",
        "## Pocket-State Enrichment",
        "",
        "Pocket states were clustered from pocket contact features rather than raw xyz coordinates. "
        "The first-pass features are distances involving F336/PHE336, F347/PHE347, H270/HD2, Q269, "
        "plus PHE336 ring/C4 geometry. This is a lightweight GNN-like featurization: residues become interpretable "
        "contact features, and the resulting states are checked for R-like/S-like enrichment.",
        "",
        pocket_state_summary.to_string(index=False) if pocket_state_summary is not None and not pocket_state_summary.empty else "Not generated.",
        "",
        "Deprecated exploratory short-vs-long plots were moved to `deprecated_old_short_long_and_C_labels/`.",
        "",
    ]
    (outdir / "README_publication_style_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    archive_deprecated_outputs(outdir)

    font_used = setup_style(args.arial_font)
    analysis_dir = Path(args.analysis_dir)
    cv = pd.read_csv(args.cv_csv)
    rmsd, rmsd_fallback = load_rmsd_or_previous_summary(analysis_dir, outdir)
    dist = load_c4o1_distance(analysis_dir, cv)

    summaries = {}
    if rmsd is not None:
        summaries.update(rmsd_fallback)
        summaries.update(plot_rmsd(rmsd, outdir))
    else:
        summaries.update(rmsd_fallback)
    summaries.update(plot_c4o1_hist(dist, outdir, args.short_cutoff))
    selected_csv = analysis_dir / "prochirality_phi360_representatives" / "selected_frames_phi360.csv"
    plot_phi360_scatter(cv, selected_csv, outdir, args.short_cutoff)
    phe_summary = plot_phe336_pi_hist(cv, outdir, args.short_cutoff, args.pi_distance_cutoff)
    plot_phe336_geometry_scatter(cv, outdir, args.short_cutoff, args.pi_distance_cutoff)
    group_summary = plot_cluster_geometry(Path(args.cluster_dir), cv, outdir)
    pro_rs = load_pro_rs(analysis_dir)
    pro_rs_summary = None
    pro_rs_pocket_summary = None
    pocket_state_summary = None
    if pro_rs is not None:
        all_assign = pd.read_csv(Path(args.cluster_dir) / "all_stride_frame_assignments_k3.csv")
        all_m = merge_assignments_with_cv(all_assign, cv)
        all_m, _ = assign_group_labels(all_m)
        pro_rs_summary = plot_pro_rs_triple(pro_rs, outdir, args.short_cutoff)
        pro_rs_pocket_summary, _ = plot_pro_rs_pocket(all_m, pro_rs, outdir, args.short_cutoff)
        pocket_state_summary = plot_pocket_state_enrichment(all_m, pro_rs, outdir, args.short_cutoff)
    write_readme(outdir, font_used, summaries, phe_summary, group_summary, pro_rs_summary, pocket_state_summary, args.short_cutoff)
    print(outdir)


if __name__ == "__main__":
    main()
