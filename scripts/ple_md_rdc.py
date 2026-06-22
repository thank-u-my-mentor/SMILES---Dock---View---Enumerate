#!/usr/bin/env python
"""Compute residue-to-ligand distance-change metrics from a GROMACS trajectory."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

try:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance, distance_array
except ImportError as exc:  # pragma: no cover - runtime dependency check
    raise SystemExit("MDAnalysis is required: conda install -c conda-forge mdanalysis") from exc


@dataclass
class ResidueTrace:
    key: str
    segid: str
    resid: int
    resname: str
    distances: list[float]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-residue ligand distance-change metrics. Positive slope/RDC "
            "means the residue moves away from the ligand; negative means it approaches."
        )
    )
    parser.add_argument("--topology", required=True, type=Path, help="Topology/structure file, e.g. .tpr or .gro")
    parser.add_argument("--trajectory", required=True, type=Path, help="Trajectory file, e.g. .xtc")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory")
    parser.add_argument("--ligand-selection", default="resname LIG", help="MDAnalysis selection for substrate/ligand")
    parser.add_argument(
        "--residue-selection",
        default="protein",
        help="MDAnalysis selection for residues to score; default protein",
    )
    parser.add_argument("--distance-cutoff", type=float, default=10.0, help="Only rank residues with min distance <= cutoff A")
    parser.add_argument("--stride", type=int, default=1, help="Use every Nth trajectory frame")
    parser.add_argument("--start-ns", type=float, help="Ignore frames before this time in ns")
    parser.add_argument("--end-ns", type=float, help="Ignore frames after this time in ns")
    parser.add_argument("--top-n", type=int, default=30, help="Number of residues in summary plots")
    parser.add_argument("--plot-traces", action="store_true", help="Write per-residue distance traces for top residues")
    parser.add_argument(
        "--reference",
        choices=["first", "mean-first-10pct"],
        default="first",
        help="Reference distance used for RDC normalization",
    )
    return parser.parse_args()


def heavy(atomgroup):
    heavy_group = atomgroup.select_atoms("not name H* and not type H")
    return heavy_group if len(heavy_group) else atomgroup


def residue_key(residue) -> str:
    segid = residue.segid.strip() or "-"
    return f"{residue.resname}:{segid}:{residue.resid}"


def min_distance(residue_atoms, ligand_atoms, box) -> float:
    if len(residue_atoms) == 0 or len(ligand_atoms) == 0:
        return math.nan
    # Small residue/ligand groups are cheap; distance_array keeps periodic boxes if present.
    distances = distance_array(residue_atoms.positions, ligand_atoms.positions, box=box)
    if distances.size == 0:
        return math.nan
    return float(np.nanmin(distances))


def linear_fit(times_ns: np.ndarray, distances: np.ndarray) -> tuple[float, float, float]:
    mask = np.isfinite(times_ns) & np.isfinite(distances)
    if mask.sum() < 2:
        return math.nan, math.nan, math.nan
    x = times_ns[mask]
    y = distances[mask]
    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else math.nan
    return float(slope), float(intercept), float(r2)


def safe_reference(distances: np.ndarray, mode: str) -> float:
    finite = distances[np.isfinite(distances)]
    if finite.size == 0:
        return math.nan
    if mode == "mean-first-10pct":
        n = max(1, int(math.ceil(0.1 * finite.size)))
        return float(np.mean(finite[:n]))
    return float(finite[0])


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_bar(rows: list[dict[str, object]], field: str, outpath: Path, title: str, top_n: int) -> None:
    ranked = sorted(rows, key=lambda row: abs(float(row[field])), reverse=True)[:top_n]
    if not ranked:
        return
    labels = [str(row["residue"]) for row in ranked][::-1]
    values = [float(row[field]) for row in ranked][::-1]
    colors = ["#2a9d8f" if value < 0 else "#b23a48" for value in values]
    height = max(5, 0.28 * len(labels) + 1.5)
    plt.figure(figsize=(9, height))
    plt.barh(labels, values, color=colors)
    plt.axvline(0.0, color="black", linewidth=0.8)
    plt.xlabel(field)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def plot_distance_traces(
    traces: dict[str, np.ndarray],
    times_ns: np.ndarray,
    rows: list[dict[str, object]],
    outpath: Path,
    top_n: int,
) -> None:
    ranked = sorted(rows, key=lambda row: abs(float(row["rdc_per_ns"])), reverse=True)[:top_n]
    if not ranked:
        return
    plt.figure(figsize=(10, 6))
    for row in ranked:
        key = str(row["residue"])
        if key not in traces:
            continue
        plt.plot(times_ns, traces[key], linewidth=1.2, label=key)
    plt.xlabel("time (ns)")
    plt.ylabel("min heavy-atom distance to ligand (A)")
    plt.title(f"Top {len(ranked)} residue-ligand distance traces")
    plt.legend(fontsize=7, ncol=2)
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    universe = mda.Universe(str(args.topology), str(args.trajectory))
    ligand = heavy(universe.select_atoms(args.ligand_selection))
    if len(ligand) == 0:
        raise SystemExit(f"No ligand atoms matched selection: {args.ligand_selection!r}")

    selected = universe.select_atoms(args.residue_selection)
    if len(selected) == 0:
        raise SystemExit(f"No residue atoms matched selection: {args.residue_selection!r}")

    residues = []
    for residue in selected.residues:
        atoms = heavy(residue.atoms)
        if len(atoms):
            residues.append((residue_key(residue), residue, atoms))
    if not residues:
        raise SystemExit("No residues with atoms found after heavy-atom filtering")

    times_ns: list[float] = []
    traces: dict[str, list[float]] = {key: [] for key, _, _ in residues}

    for frame_index, ts in enumerate(universe.trajectory):
        if frame_index % max(args.stride, 1) != 0:
            continue
        time_ns = float(ts.time) / 1000.0
        if args.start_ns is not None and time_ns < args.start_ns:
            continue
        if args.end_ns is not None and time_ns > args.end_ns:
            continue
        times_ns.append(time_ns)
        box = ts.dimensions if ts.dimensions is not None and np.any(ts.dimensions[:3]) else None
        for key, _, atoms in residues:
            traces[key].append(min_distance(atoms, ligand, box))

    if len(times_ns) < 2:
        raise SystemExit("Need at least two trajectory frames after filtering")

    times = np.asarray(times_ns, dtype=float)
    rows: list[dict[str, object]] = []
    trace_arrays: dict[str, np.ndarray] = {}
    for key, residue, _ in residues:
        distances = np.asarray(traces[key], dtype=float)
        trace_arrays[key] = distances
        finite = distances[np.isfinite(distances)]
        if finite.size < 2:
            continue
        slope, intercept, r2 = linear_fit(times, distances)
        ref = safe_reference(distances, args.reference)
        final_distance = float(finite[-1])
        initial_distance = float(finite[0])
        rdc_per_ns = slope / ref if np.isfinite(ref) and ref != 0 else math.nan
        rows.append(
            {
                "residue": key,
                "resname": residue.resname,
                "segid": residue.segid.strip() or "-",
                "resid": int(residue.resid),
                "initial_distance_A": round(initial_distance, 4),
                "final_distance_A": round(final_distance, 4),
                "reference_distance_A": round(ref, 4) if np.isfinite(ref) else "",
                "mean_distance_A": round(float(np.mean(finite)), 4),
                "min_distance_A": round(float(np.min(finite)), 4),
                "max_distance_A": round(float(np.max(finite)), 4),
                "delta_final_minus_initial_A": round(final_distance - initial_distance, 4),
                "slope_A_per_ns": round(slope, 6) if np.isfinite(slope) else "",
                "rdc_per_ns": round(rdc_per_ns, 8) if np.isfinite(rdc_per_ns) else "",
                "fit_r2": round(r2, 5) if np.isfinite(r2) else "",
                "n_frames": int(finite.size),
                "within_cutoff": bool(float(np.min(finite)) <= args.distance_cutoff),
            }
        )

    all_csv = args.outdir / "residue_ligand_rdc_all.csv"
    cutoff_csv = args.outdir / "residue_ligand_rdc_within_cutoff.csv"
    write_csv(all_csv, sorted(rows, key=lambda row: (float(row["min_distance_A"]), int(row["resid"]))))

    cutoff_rows = [row for row in rows if row["within_cutoff"]]
    cutoff_rows.sort(key=lambda row: abs(float(row["rdc_per_ns"])), reverse=True)
    write_csv(cutoff_csv, cutoff_rows)

    # Full distance time series for interactive dashboards.
    all_keys = [str(row["residue"]) for row in sorted(rows, key=lambda row: int(row["resid"]))]
    with (args.outdir / "residue_ligand_distance_timeseries_all.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_ns", *all_keys])
        for i, time_ns in enumerate(times):
            writer.writerow([round(float(time_ns), 5), *[round(float(trace_arrays[key][i]), 4) for key in all_keys]])

    # Distance time series for nearby residues; useful for manual plotting and sanity checks.
    nearby_keys = [str(row["residue"]) for row in cutoff_rows]
    with (args.outdir / "residue_ligand_distance_timeseries.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_ns", *nearby_keys])
        for i, time_ns in enumerate(times):
            writer.writerow([round(float(time_ns), 5), *[round(float(trace_arrays[key][i]), 4) for key in nearby_keys]])

    plot_bar(cutoff_rows, "rdc_per_ns", args.outdir / "top_abs_rdc_per_ns.png", "Top absolute residue-ligand RDC", args.top_n)
    plot_bar(
        cutoff_rows,
        "slope_A_per_ns",
        args.outdir / "top_abs_slope_A_per_ns.png",
        "Top absolute residue-ligand distance slope",
        args.top_n,
    )
    if args.plot_traces:
        plot_distance_traces(
            trace_arrays,
            times,
            cutoff_rows,
            args.outdir / "top_rdc_distance_traces.png",
            min(args.top_n, 12),
        )

    metadata = {
        "topology": str(args.topology),
        "trajectory": str(args.trajectory),
        "ligand_selection": args.ligand_selection,
        "residue_selection": args.residue_selection,
        "distance_cutoff_A": args.distance_cutoff,
        "stride": args.stride,
        "start_ns": args.start_ns,
        "end_ns": args.end_ns,
        "frames_used": len(times),
        "time_min_ns": float(times.min()),
        "time_max_ns": float(times.max()),
        "rdc_definition": "rdc_per_ns = linear_fit_slope_A_per_ns / reference_distance_A",
        "sign": "positive means moving away from ligand; negative means approaching ligand",
    }
    (args.outdir / "rdc_metadata.json").write_text(
        __import__("json").dumps(metadata, indent=2),
        encoding="utf-8",
    )

    print(f"frames_used={len(times)} time_ns={times.min():.3f}-{times.max():.3f}")
    print(f"ligand_atoms={len(ligand)} residues_scored={len(rows)} nearby_residues={len(cutoff_rows)}")
    print(f"all_csv={all_csv}")
    print(f"cutoff_csv={cutoff_csv}")
    if cutoff_rows:
        print("top_abs_rdc:")
        for row in cutoff_rows[:10]:
            print(
                f"  {row['residue']} rdc_per_ns={row['rdc_per_ns']} "
                f"slope_A_per_ns={row['slope_A_per_ns']} min_A={row['min_distance_A']}"
            )


if __name__ == "__main__":
    main()
