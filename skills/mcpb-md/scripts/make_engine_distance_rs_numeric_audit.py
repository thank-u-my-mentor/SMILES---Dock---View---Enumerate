#!/usr/bin/env python3
"""Write exact C4-O1/R-S numeric audits for GROMACS vs Amber trajectories."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


RS_MAP = {
    "R_like_by_product_reference": "R-like",
    "S_like_by_product_reference": "S-like",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gromacs-rs-csv", required=True)
    parser.add_argument("--amber-rs-csv", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def load_one(path: Path, dataset: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df[["frame", "time_ns", "C4_O1_A", "signed_triple_O1_C5phenyl_C3chain_about_C4", "pseudo_RS_by_product_ref"]].copy()
    df["dataset"] = dataset
    df["pro_rs_label"] = df["pseudo_RS_by_product_ref"].map(RS_MAP)
    return df


def count_rs(sub: pd.DataFrame) -> dict[str, float | int]:
    counts = sub["pro_rs_label"].value_counts()
    n = int(len(sub))
    r = int(counts.get("R-like", 0))
    s = int(counts.get("S-like", 0))
    return {
        "n_frames": n,
        "R_like_n": r,
        "S_like_n": s,
        "R_like_fraction": float(r / n) if n else np.nan,
        "S_like_fraction": float(s / n) if n else np.nan,
    }


def summarize_engine(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    ranges = [
        ("all_frames", np.full(len(df), True)),
        ("C4O1_le_3.0A", df["C4_O1_A"] <= 3.0),
        ("C4O1_lt_3.0A", df["C4_O1_A"] < 3.0),
        ("C4O1_2.5_to_4.0A", (df["C4_O1_A"] >= 2.5) & (df["C4_O1_A"] < 4.0)),
        ("C4O1_3.0_to_4.0A", (df["C4_O1_A"] >= 3.0) & (df["C4_O1_A"] < 4.0)),
        ("C4O1_ge_4.0A", df["C4_O1_A"] >= 4.0),
    ]
    for subset, mask in ranges:
        sub = df.loc[mask].copy()
        row = {"dataset": str(df["dataset"].iloc[0]), "subset": subset}
        row.update(count_rs(sub))
        row.update(
            {
                "fraction_of_dataset": float(len(sub) / len(df)) if len(df) else np.nan,
                "median_C4_O1_A": float(sub["C4_O1_A"].median()) if len(sub) else np.nan,
                "mean_C4_O1_A": float(sub["C4_O1_A"].mean()) if len(sub) else np.nan,
                "min_C4_O1_A": float(sub["C4_O1_A"].min()) if len(sub) else np.nan,
                "median_triple_A3": float(sub["signed_triple_O1_C5phenyl_C3chain_about_C4"].median()) if len(sub) else np.nan,
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def bin_counts(df: pd.DataFrame) -> pd.DataFrame:
    edges = [2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, np.inf]
    labels = ["2.50-2.75", "2.75-3.00", "3.00-3.25", "3.25-3.50", "3.50-3.75", "3.75-4.00", "4.00-4.25", "4.25-4.50", "4.50-4.75", "4.75-5.00", ">=5.00"]
    rows = []
    binned = pd.cut(df["C4_O1_A"], bins=edges, labels=labels, right=False, include_lowest=True)
    for label in labels:
        sub = df.loc[binned == label]
        row = {"dataset": str(df["dataset"].iloc[0]), "C4_O1_bin_A": label}
        row.update(count_rs(sub))
        row["fraction_of_dataset"] = float(len(sub) / len(df)) if len(df) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def threshold_counts(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cutoff in [2.70, 2.75, 2.80, 2.90, 3.00, 3.25, 3.50, 3.75, 4.00]:
        sub = df[df["C4_O1_A"] <= cutoff]
        row = {"dataset": str(df["dataset"].iloc[0]), "cutoff_A": cutoff}
        row.update(count_rs(sub))
        row["fraction_of_dataset"] = float(len(sub) / len(df)) if len(df) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def compare_bins(bins: pd.DataFrame) -> pd.DataFrame:
    pivot = bins.pivot(index="C4_O1_bin_A", columns="dataset", values=["n_frames", "fraction_of_dataset", "R_like_fraction"])
    pivot.columns = [f"{a}_{b}".replace(" ", "_") for a, b in pivot.columns]
    pivot = pivot.reset_index()
    if {"fraction_of_dataset_GROMACS", "fraction_of_dataset_Amber_pmemd"}.issubset(pivot.columns):
        pivot["Amber_minus_GROMACS_fraction"] = pivot["fraction_of_dataset_Amber_pmemd"] - pivot["fraction_of_dataset_GROMACS"]
    if {"R_like_fraction_GROMACS", "R_like_fraction_Amber_pmemd"}.issubset(pivot.columns):
        pivot["Amber_minus_GROMACS_R_like_fraction"] = pivot["R_like_fraction_Amber_pmemd"] - pivot["R_like_fraction_GROMACS"]
    return pivot


def write_markdown(outdir: Path, engine_summary: pd.DataFrame, bins: pd.DataFrame, thresholds: pd.DataFrame, comparison: pd.DataFrame) -> None:
    lines = [
        "# Exact C4-O1 / R-S Numeric Audit",
        "",
        "This audit uses all saved frames, not stride-10 subsampling.",
        "",
        "## Engine/Subsets",
        "",
        engine_summary.to_string(index=False),
        "",
        "## Distance Bin Counts",
        "",
        bins.to_string(index=False),
        "",
        "## Threshold Counts",
        "",
        thresholds.to_string(index=False),
        "",
        "## Amber Minus GROMACS Bin Differences",
        "",
        comparison.to_string(index=False),
        "",
    ]
    (outdir / "README_exact_C4O1_RS_numeric_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gmx = load_one(Path(args.gromacs_rs_csv), "GROMACS")
    amber = load_one(Path(args.amber_rs_csv), "Amber pmemd")
    data = pd.concat([gmx, amber], ignore_index=True)

    engine_summary = pd.concat([summarize_engine(gmx), summarize_engine(amber)], ignore_index=True)
    bins = pd.concat([bin_counts(gmx), bin_counts(amber)], ignore_index=True)
    thresholds = pd.concat([threshold_counts(gmx), threshold_counts(amber)], ignore_index=True)
    comparison = compare_bins(bins)

    engine_summary.to_csv(outdir / "exact_C4O1_RS_engine_subset_counts.csv", index=False)
    bins.to_csv(outdir / "exact_C4O1_RS_distance_bin_counts.csv", index=False)
    thresholds.to_csv(outdir / "exact_C4O1_RS_threshold_counts.csv", index=False)
    comparison.to_csv(outdir / "exact_C4O1_RS_bin_engine_differences.csv", index=False)
    write_markdown(outdir, engine_summary, bins, thresholds, comparison)
    print(outdir)


if __name__ == "__main__":
    main()
