#!/usr/bin/env python
"""Measure Fe-core coordination distances from a GROMACS trajectory."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import MDAnalysis as mda
except ImportError as exc:  # pragma: no cover
    raise SystemExit("MDAnalysis is required") from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--structure", required=True, type=Path, help="Structure file readable by MDAnalysis, e.g. .gro")
    parser.add_argument("--trajectory", required=True, type=Path, help="Trajectory file, e.g. .xtc")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument(
        "--pair",
        action="append",
        required=True,
        help="Pair as label,atom_number_1,atom_number_2 using 1-based GROMACS atom numbers",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    universe = mda.Universe(str(args.structure), str(args.trajectory))

    pairs: list[tuple[str, int, int]] = []
    for raw in args.pair:
        label, first, second = raw.split(",", 2)
        pairs.append((label, int(first), int(second)))

    rows: list[dict[str, float | str]] = []
    for ts in universe.trajectory:
        row: dict[str, float | str] = {"time_ns": round(float(ts.time) / 1000.0, 5)}
        for label, first, second in pairs:
            atom1 = universe.atoms[first - 1]
            atom2 = universe.atoms[second - 1]
            distance_a = float(np.linalg.norm(atom1.position - atom2.position))
            row[label] = round(distance_a, 5)
        rows.append(row)

    table = pd.DataFrame(rows)
    table.to_csv(args.outdir / "fe_core_distance_timeseries.csv", index=False)

    summary_rows = []
    for label, _, _ in pairs:
        values = table[label].to_numpy(dtype=float)
        summary_rows.append(
            {
                "pair": label,
                "mean_A": round(float(np.mean(values)), 5),
                "std_A": round(float(np.std(values)), 5),
                "min_A": round(float(np.min(values)), 5),
                "max_A": round(float(np.max(values)), 5),
                "n_frames": int(values.size),
            }
        )
    pd.DataFrame(summary_rows).to_csv(args.outdir / "fe_core_distance_summary.csv", index=False)
    print(args.outdir / "fe_core_distance_summary.csv")


if __name__ == "__main__":
    main()
