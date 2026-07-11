#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd


def one(u: mda.Universe, selection: str):
    ag = u.select_atoms(selection)
    if len(ag) != 1:
        raise SystemExit(f"{selection!r} matched {len(ag)} atoms, expected 1")
    return ag[0]


def signed_triple(center: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    return float(np.linalg.det(np.vstack([a - center, b - center, c - center])))


def assignment(value: float) -> str:
    # Calibrated from product references:
    # R product: signed_triple(O1, C5 phenyl, C3 chain about C4) < 0
    # S product: signed_triple(O1, C5 phenyl, C3 chain about C4) > 0
    return "R_like_by_product_reference" if value < 0 else "S_like_by_product_reference"


def summarize(df: pd.DataFrame, label: str) -> pd.DataFrame:
    rows = []
    for cutoff in [3.50, 3.25, 3.00, 2.90, 2.80, 2.75, 2.70]:
        sub = df[df["C4_O1_A"] <= cutoff]
        rows.append(row_summary(sub, label, f"C4_O1_le_{cutoff:.2f}A"))
    sorted_df = df.sort_values("C4_O1_A")
    for n in [5, 10, 20, 50, 100, 500, 1000]:
        rows.append(row_summary(sorted_df.head(n), label, f"top_{n}_shortest_C4_O1"))
    rows.append(row_summary(df, label, "all_frames"))
    return pd.DataFrame(rows)


def row_summary(sub: pd.DataFrame, dataset: str, subset: str) -> dict[str, object]:
    total = len(sub)
    counts = sub["pseudo_RS_by_product_ref"].value_counts()
    return {
        "dataset": dataset,
        "subset": subset,
        "n_frames": total,
        "C4_O1_median_A": float(sub["C4_O1_A"].median()) if total else np.nan,
        "C4_O1_mean_A": float(sub["C4_O1_A"].mean()) if total else np.nan,
        "C4_O1_min_A": float(sub["C4_O1_A"].min()) if total else np.nan,
        "R_like_n": int(counts.get("R_like_by_product_reference", 0)),
        "R_like_fraction": float(counts.get("R_like_by_product_reference", 0) / total) if total else np.nan,
        "S_like_n": int(counts.get("S_like_by_product_reference", 0)),
        "S_like_fraction": float(counts.get("S_like_by_product_reference", 0) / total) if total else np.nan,
        "median_signed_triple": float(sub["signed_triple_O1_C5phenyl_C3chain_about_C4"].median()) if total else np.nan,
    }


def analyze_amber_ul1(prmtop: Path, nc: Path) -> pd.DataFrame:
    u = mda.Universe(str(prmtop), str(nc))
    atoms = {name: one(u, f"name {name}") for name in ["O1", "C3", "C4", "C5", "H02"]}
    rows = []
    for ts in u.trajectory:
        c4 = atoms["C4"].position.copy()
        o1 = atoms["O1"].position.copy()
        c5 = atoms["C5"].position.copy()
        c3 = atoms["C3"].position.copy()
        triple = signed_triple(c4, o1, c5, c3)
        rows.append(
            {
                "frame": int(ts.frame),
                "time_ps": float(ts.time),
                "time_ns": float(ts.time / 1000.0),
                "C4_O1_A": float(np.linalg.norm(c4 - o1)),
                "signed_triple_O1_C5phenyl_C3chain_about_C4": triple,
                "pseudo_RS_by_product_ref": assignment(triple),
            }
        )
    return pd.DataFrame(rows)


def read_cpptraj_rmsd(path: Path) -> dict[str, float | int]:
    if not path.exists():
        return {}
    df = pd.read_csv(path, comment="#", sep=r"\s+", header=None)
    if df.shape[1] < 2:
        return {}
    values = df.iloc[:, 1].astype(float)
    return {
        "amber_CA_RMSD_n_frames": int(len(values)),
        "amber_CA_RMSD_median_A": float(values.median()),
        "amber_CA_RMSD_mean_A": float(values.mean()),
        "amber_CA_RMSD_last_A": float(values.iloc[-1]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--amber-ul1-prmtop", required=True, type=Path)
    parser.add_argument("--amber-ul1-nc", required=True, type=Path)
    parser.add_argument("--gromacs-timeseries", required=True, type=Path)
    parser.add_argument("--amber-ca-rmsd", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    amber = analyze_amber_ul1(args.amber_ul1_prmtop, args.amber_ul1_nc)
    amber.to_csv(args.outdir / "amber_product_calibrated_pseudo_RS_timeseries.csv", index=False)
    amber_summary = summarize(amber, "amber_pmemd")
    amber_summary.to_csv(args.outdir / "amber_product_calibrated_pseudo_RS_summary.csv", index=False)

    gmx = pd.read_csv(args.gromacs_timeseries)
    gmx = gmx[["frame", "time_ps", "time_ns", "C4_O1_A", "signed_triple_O1_C5phenyl_C3chain_about_C4", "pseudo_RS_by_product_ref"]].copy()
    gmx_summary = summarize(gmx, "gromacs")

    combined = pd.concat([gmx_summary, amber_summary], ignore_index=True)
    combined.to_csv(args.outdir / "amber_vs_gromacs_pseudo_RS_summary.csv", index=False)

    rmsd_stats = read_cpptraj_rmsd(args.amber_ca_rmsd)
    pd.DataFrame([rmsd_stats]).to_csv(args.outdir / "amber_CA_RMSD_summary.csv", index=False)

    print("Wrote", args.outdir)
    print(combined[combined["subset"].isin(["all_frames", "C4_O1_le_3.00A", "C4_O1_le_2.80A", "top_100_shortest_C4_O1"])].to_string(index=False))
    if rmsd_stats:
        print(rmsd_stats)


if __name__ == "__main__":
    main()
