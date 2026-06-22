#!/usr/bin/env python
"""Plot residue-ligand distance-fluctuation velocity distributions."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute V = delta(distance)/delta(time) distributions from "
            "residue_ligand_distance_timeseries.csv. Negative V means the residue "
            "approaches the ligand; positive V means it moves away."
        )
    )
    parser.add_argument("--timeseries", required=True, type=Path)
    parser.add_argument("--rdc-csv", type=Path, help="Optional RDC ranking CSV used to choose default residues")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument(
        "--residues",
        nargs="*",
        help="Residue labels exactly as in the timeseries header, e.g. GLU:SYSTEM:319",
    )
    parser.add_argument("--top-n", type=int, default=6, help="Use top N residues from --rdc-csv when --residues is omitted")
    parser.add_argument("--max-x", type=float, default=20.0, help="Plot x range in V x 1e-2 A/ps units")
    parser.add_argument(
        "--low-dynamic-threshold",
        type=float,
        default=5.0,
        help="Threshold line in V x 1e-2 A/ps units; default +/-5 equals +/-0.05 A/ps",
    )
    parser.add_argument("--bins", type=int, default=44)
    parser.add_argument("--title", default="Distribution density of distance fluctuations between residues and substrate")
    return parser.parse_args()


def choose_residues(args: argparse.Namespace, columns: list[str]) -> list[str]:
    residue_columns = [column for column in columns if column != "time_ns"]
    if args.residues:
        missing = [residue for residue in args.residues if residue not in residue_columns]
        if missing:
            raise SystemExit(f"Residues not found in timeseries header: {missing}")
        return args.residues
    if args.rdc_csv and args.rdc_csv.exists():
        with args.rdc_csv.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        ranked = [row["residue"] for row in rows if row.get("residue") in residue_columns]
        if ranked:
            return ranked[: args.top_n]
    return residue_columns[: args.top_n]


def velocity_distribution(times_ns: np.ndarray, distances: np.ndarray) -> np.ndarray:
    times_ps = times_ns * 1000.0
    dt = np.diff(times_ps)
    dd = np.diff(distances)
    mask = np.isfinite(dt) & np.isfinite(dd) & (dt > 0)
    if not np.any(mask):
        return np.asarray([], dtype=float)
    return dd[mask] / dt[mask]


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    table = pd.read_csv(args.timeseries)
    if "time_ns" not in table.columns:
        raise SystemExit("timeseries CSV must contain a time_ns column")
    residues = choose_residues(args, list(table.columns))
    times_ns = table["time_ns"].to_numpy(dtype=float)

    velocity_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    distributions: dict[str, np.ndarray] = {}
    for residue in residues:
        distances = table[residue].to_numpy(dtype=float)
        velocities = velocity_distribution(times_ns, distances)
        distributions[residue] = velocities
        for idx, value in enumerate(velocities):
            velocity_rows.append(
                {
                    "residue": residue,
                    "interval_index": idx,
                    "time_mid_ns": round(float((times_ns[idx] + times_ns[idx + 1]) / 2.0), 5),
                    "velocity_A_per_ps": round(float(value), 8),
                    "velocity_x1e_minus2_A_per_ps": round(float(value / 0.01), 5),
                }
            )
        if velocities.size:
            summary_rows.append(
                {
                    "residue": residue,
                    "n_intervals": int(velocities.size),
                    "mean_velocity_A_per_ps": round(float(np.mean(velocities)), 8),
                    "median_velocity_A_per_ps": round(float(np.median(velocities)), 8),
                    "std_velocity_A_per_ps": round(float(np.std(velocities)), 8),
                    "mean_velocity_x1e_minus2_A_per_ps": round(float(np.mean(velocities) / 0.01), 5),
                    "std_velocity_x1e_minus2_A_per_ps": round(float(np.std(velocities) / 0.01), 5),
                    "fraction_approaching": round(float(np.mean(velocities < 0)), 5),
                    "fraction_departing": round(float(np.mean(velocities > 0)), 5),
                    "p05_x1e_minus2": round(float(np.percentile(velocities / 0.01, 5)), 5),
                    "p95_x1e_minus2": round(float(np.percentile(velocities / 0.01, 95)), 5),
                }
            )

    pd.DataFrame(velocity_rows).to_csv(args.outdir / "residue_ligand_velocity_distribution.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(args.outdir / "residue_ligand_velocity_distribution_summary.csv", index=False)

    bins = np.linspace(-args.max_x, args.max_x, args.bins + 1)
    colors = [
        "#6c757d",
        "#e6a57e",
        "#e9d985",
        "#8ab17d",
        "#80b1d3",
        "#b39ddb",
        "#f28482",
        "#84a59d",
    ]

    plt.figure(figsize=(10, 5.8))
    for index, residue in enumerate(residues):
        velocities = distributions.get(residue, np.asarray([]))
        if velocities.size == 0:
            continue
        values = velocities / 0.01
        plt.hist(
            values,
            bins=bins,
            density=True,
            alpha=0.34,
            color=colors[index % len(colors)],
            label=residue,
            edgecolor=colors[index % len(colors)],
            linewidth=0.4,
        )

    threshold = args.low_dynamic_threshold
    plt.axvline(-threshold, color="#666666", linestyle="--", linewidth=1.0)
    plt.axvline(threshold, color="#666666", linestyle="--", linewidth=1.0)
    plt.axvline(0.0, color="#222222", linewidth=0.8)
    ymax = plt.ylim()[1]
    plt.text(-args.max_x * 0.86, ymax * 0.55, "High dynamic change", ha="left", va="center")
    plt.text(0.0, ymax * 0.55, "Low dynamic change", ha="center", va="center")
    plt.text(args.max_x * 0.55, ymax * 0.55, "High dynamic change", ha="left", va="center")
    plt.xlabel(r"V x $10^{-2}$ (A/ps)")
    plt.ylabel("Density")
    plt.title(args.title, color="#c1121f", fontsize=12)
    plt.legend(ncol=min(3, max(1, len(residues))), fontsize=8, frameon=False, loc="upper center")
    plt.tight_layout()
    plt.savefig(args.outdir / "velocity_distribution_density.png", dpi=220)
    plt.close()

    print(f"residues={','.join(residues)}")
    print(f"plot={args.outdir / 'velocity_distribution_density.png'}")
    print(f"summary={args.outdir / 'residue_ligand_velocity_distribution_summary.csv'}")


if __name__ == "__main__":
    main()
