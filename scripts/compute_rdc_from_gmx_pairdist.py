#!/usr/bin/env python3
"""Compute residue-ligand RDC metrics from a GROMACS pairdist XVG file."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path


def parse_xvg(path: Path) -> tuple[list[str], list[float], list[list[float]]]:
    labels: list[str] = []
    times: list[float] = []
    values: list[list[float]] = []
    series_re = re.compile(r'@\s+s(\d+)\s+legend\s+"([^"]+)"')
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("@"):
            match = series_re.match(line.strip())
            if match:
                idx = int(match.group(1))
                while len(labels) <= idx:
                    labels.append("")
                labels[idx] = match.group(2)
            continue
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        times.append(float(parts[0]))
        values.append([float(v) for v in parts[1:]])
    if not values:
        raise ValueError(f"No numeric data rows found in {path}")
    n_columns = len(values[0])
    if len(labels) != n_columns or any(not label for label in labels):
        labels = [f"series_{i}" for i in range(n_columns)]
    return labels, times, values


def normalize_label(label: str) -> tuple[str, str, str, int]:
    # GROMACS pairdist legends commonly look like:
    # "min dist. Protein_1-LIG" or "Protein_1-LIG".
    cleaned = label.replace("min dist.", "").strip()
    cleaned = cleaned.split("-", 1)[0].strip()
    # Try to preserve resname/resid from selection labels such as ALA-12.
    match = re.search(r"([A-Za-z]{3,4})[-_:]?([A-Za-z]?:)?(\d+)", cleaned)
    if match:
        resname = match.group(1).upper()
        resid = int(match.group(3))
        segid = "A"
        residue = f"{resname}:{segid}:{resid}"
        return residue, resname, segid, resid
    # Fall back to GROMACS residue-group index.
    match = re.search(r"(\d+)", cleaned)
    resid = int(match.group(1)) if match else 0
    residue = f"RES:A:{resid}"
    return residue, "RES", "A", resid


def linear_fit(x: list[float], y: list[float]) -> tuple[float, float, float]:
    pairs = [(a, b) for a, b in zip(x, y) if math.isfinite(a) and math.isfinite(b)]
    if len(pairs) < 2:
        return math.nan, math.nan, math.nan
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]
    n = len(xs)
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    var_x = sum((v - mean_x) ** 2 for v in xs)
    if var_x == 0:
        return math.nan, math.nan, math.nan
    cov = sum((a - mean_x) * (b - mean_y) for a, b in zip(xs, ys))
    slope = cov / var_x
    intercept = mean_y - slope * mean_x
    pred = [slope * a + intercept for a in xs]
    ss_res = sum((b - p) ** 2 for b, p in zip(ys, pred))
    ss_tot = sum((b - mean_y) ** 2 for b in ys)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else math.nan
    return slope, intercept, r2


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xvg", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--cutoff", type=float, default=10.0)
    parser.add_argument("--time-unit", choices=["ps", "ns"], default="ps")
    parser.add_argument("--label-csv", type=Path, help="Optional CSV with columns series,residue,resname,segid,resid")
    parser.add_argument("--metadata", type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    labels, raw_times, matrix = parse_xvg(args.xvg)
    times_ns = [t / 1000.0 if args.time_unit == "ps" else t for t in raw_times]
    label_override: dict[int, tuple[str, str, str, int]] = {}
    if args.label_csv and args.label_csv.exists():
        with args.label_csv.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                idx = int(row["series"])
                label_override[idx] = (
                    row["residue"],
                    row["resname"],
                    row.get("segid") or "A",
                    int(row["resid"]),
                )

    rows: list[dict[str, object]] = []
    traces_by_residue: dict[str, list[float]] = {}
    for col, label in enumerate(labels):
        distances_nm = [row[col] for row in matrix]
        distances = [v * 10.0 for v in distances_nm]  # GROMACS pairdist writes nm.
        residue, resname, segid, resid = label_override.get(col, normalize_label(label))
        finite = [v for v in distances if math.isfinite(v)]
        if len(finite) < 2:
            continue
        slope, _, r2 = linear_fit(times_ns, distances)
        initial = finite[0]
        final = finite[-1]
        ref = initial if initial else math.nan
        rdc = slope / ref if math.isfinite(ref) and ref != 0 else math.nan
        traces_by_residue[residue] = distances
        rows.append(
            {
                "residue": residue,
                "resname": resname,
                "segid": segid,
                "resid": resid,
                "initial_distance_A": round(initial, 4),
                "final_distance_A": round(final, 4),
                "reference_distance_A": round(ref, 4) if math.isfinite(ref) else "",
                "mean_distance_A": round(sum(finite) / len(finite), 4),
                "min_distance_A": round(min(finite), 4),
                "max_distance_A": round(max(finite), 4),
                "delta_final_minus_initial_A": round(final - initial, 4),
                "slope_A_per_ns": round(slope, 6) if math.isfinite(slope) else "",
                "rdc_per_ns": round(rdc, 8) if math.isfinite(rdc) else "",
                "fit_r2": round(r2, 5) if math.isfinite(r2) else "",
                "n_frames": len(finite),
                "within_cutoff": min(finite) <= args.cutoff,
            }
        )

    rows_by_min = sorted(rows, key=lambda row: (float(row["min_distance_A"]), int(row["resid"])))
    rows_by_rdc = sorted(rows, key=lambda row: abs(float(row["rdc_per_ns"] or 0.0)), reverse=True)
    write_csv(args.outdir / "residue_ligand_rdc_all.csv", rows_by_min)
    write_csv(args.outdir / "residue_ligand_rdc_within_cutoff.csv", [r for r in rows_by_rdc if r["within_cutoff"]])

    sorted_keys = [str(row["residue"]) for row in sorted(rows, key=lambda row: int(row["resid"]))]
    with (args.outdir / "residue_ligand_distance_timeseries_all.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_ns", *sorted_keys])
        for i, time_ns in enumerate(times_ns):
            writer.writerow([round(time_ns, 5), *[round(traces_by_residue[key][i], 4) for key in sorted_keys]])

    if args.metadata:
        metadata = {
            "source_xvg": str(args.xvg),
            "n_frames": len(times_ns),
            "n_residues": len(rows),
            "cutoff_A": args.cutoff,
            "method": "gmx pairdist minimum heavy-atom residue-to-LIG distance; RDC=slope/reference_distance",
        }
        args.metadata.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"rows={len(rows)} frames={len(times_ns)} out={args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
