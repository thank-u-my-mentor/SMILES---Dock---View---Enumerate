#!/usr/bin/env python3
"""Stdlib-only 2D density grid scan for reactive CV plots.

This fallback version avoids pandas/matplotlib so it can run from a plain
Windows Python when WSL or scientific Python packages are unavailable.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import struct
import zlib
from pathlib import Path


XCOL = "dihedral_C1_O1_C4_C5_deg"
YCOL = "dihedral_C1_O1_C4_H02_deg"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--short-cutoff", type=float, default=3.0)
    parser.add_argument("--all-frames", action="store_true", help="Use all frames instead of filtering by C4-O1 distance.")
    parser.add_argument("--bins", default="60,80,100,120,160,200,240,300,360,480")
    parser.add_argument("--size", type=int, default=900)
    parser.add_argument("--sigma", type=float, default=1.0)
    return parser.parse_args()


def read_points(path: Path, short_cutoff: float, all_frames: bool = False) -> list[tuple[float, float]]:
    points: list[tuple[float, float]] = []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                d = float(row["C4_O1_A"])
                x = float(row[XCOL])
                y = abs(float(row[YCOL]))
            except (KeyError, TypeError, ValueError):
                continue
            if (all_frames or d <= short_cutoff) and -180.0 <= x <= 180.0 and 0.0 <= y <= 180.0:
                points.append((x, y))
    return points


def histogram(points: list[tuple[float, float]], bins: int) -> list[list[int]]:
    grid = [[0 for _ in range(bins)] for _ in range(bins)]
    for x, y in points:
        ix = int((x + 180.0) / 360.0 * bins)
        iy = int(y / 180.0 * bins)
        ix = max(0, min(bins - 1, ix))
        iy = max(0, min(bins - 1, iy))
        grid[iy][ix] += 1
    return grid


def gaussian_kernel1d(sigma: float) -> list[float]:
    if sigma <= 0:
        return [1.0]
    radius = max(2, int(round(4.0 * sigma)))
    xs = list(range(-radius, radius + 1))
    vals = [math.exp(-(x * x) / (2.0 * sigma * sigma)) for x in xs]
    total = sum(vals)
    return [v / total for v in vals]


def smooth_grid(grid: list[list[int]], sigma: float) -> list[list[float]]:
    kernel = gaussian_kernel1d(sigma)
    radius = len(kernel) // 2
    h = len(grid)
    w = len(grid[0])
    temp = [[0.0 for _ in range(w)] for _ in range(h)]
    out = [[0.0 for _ in range(w)] for _ in range(h)]
    for y in range(h):
        for x in range(w):
            acc = 0.0
            for k, weight in enumerate(kernel):
                xx = min(w - 1, max(0, x + k - radius))
                acc += grid[y][xx] * weight
            temp[y][x] = acc
    for y in range(h):
        for x in range(w):
            acc = 0.0
            for k, weight in enumerate(kernel):
                yy = min(h - 1, max(0, y + k - radius))
                acc += temp[yy][x] * weight
            out[y][x] = acc
    return out


def metrics(grid: list[list[int]]) -> dict[str, float]:
    vals = [v for row in grid for v in row if v > 0]
    cells = len(grid) * len(grid[0])
    if not vals:
        return {
            "occupied_cells": 0,
            "occupied_fraction": 0.0,
            "mean_count_occupied": 0.0,
            "single_count_fraction_occupied": 0.0,
            "max_count": 0,
        }
    return {
        "occupied_cells": len(vals),
        "occupied_fraction": len(vals) / cells,
        "mean_count_occupied": sum(vals) / len(vals),
        "single_count_fraction_occupied": sum(1 for v in vals if v == 1) / len(vals),
        "max_count": max(vals),
    }


def percentile(vals: list[float], pct: float) -> float:
    if not vals:
        return 1.0
    vals = sorted(vals)
    idx = int(round((len(vals) - 1) * pct / 100.0))
    return vals[max(0, min(len(vals) - 1, idx))]


def png_chunk(tag: bytes, data: bytes) -> bytes:
    return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)


def write_png_rgb(path: Path, pixels: list[list[tuple[int, int, int]]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    height = len(pixels)
    width = len(pixels[0])
    raw = bytearray()
    for row in pixels:
        raw.append(0)
        for r, g, b in row:
            raw.extend((r, g, b))
    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    data = b"\x89PNG\r\n\x1a\n"
    data += png_chunk(b"IHDR", ihdr)
    data += png_chunk(b"IDAT", zlib.compress(bytes(raw), 9))
    data += png_chunk(b"IEND", b"")
    path.write_bytes(data)


def render_grid(grid_values: list[list[float]], out_png: Path, width: int, height: int, add_zero_line: bool = True) -> None:
    logs = [math.log10(v + 1.0) for row in grid_values for v in row]
    vmax = max(0.2, percentile(logs, 99.5))
    bins_y = len(grid_values)
    bins_x = len(grid_values[0])
    pixels: list[list[tuple[int, int, int]]] = []
    for py in range(height):
        row: list[tuple[int, int, int]] = []
        # Put y=180 at image top.
        iy = bins_y - 1 - int(py / height * bins_y)
        iy = max(0, min(bins_y - 1, iy))
        for px in range(width):
            ix = int(px / width * bins_x)
            ix = max(0, min(bins_x - 1, ix))
            value = math.log10(grid_values[iy][ix] + 1.0) / vmax
            value = max(0.0, min(1.0, value))
            gray = int(round(255 - 245 * value))
            if add_zero_line and abs(px - width // 2) <= 1:
                gray = min(gray, 150)
            row.append((gray, gray, gray))
        pixels.append(row)
    # Black border.
    for x in range(width):
        pixels[0][x] = (0, 0, 0)
        pixels[-1][x] = (0, 0, 0)
    for y in range(height):
        pixels[y][0] = (0, 0, 0)
        pixels[y][-1] = (0, 0, 0)
    write_png_rgb(out_png, pixels)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    raw_dir = outdir / "raw_counts"
    soft_dir = outdir / "soft_sigma1"
    raw_dir.mkdir(parents=True, exist_ok=True)
    soft_dir.mkdir(parents=True, exist_ok=True)

    points = read_points(Path(args.csv), args.short_cutoff, all_frames=args.all_frames)
    bins_list = [int(v.strip()) for v in args.bins.split(",") if v.strip()]
    width = args.size
    height = max(1, args.size // 2)
    rows = []
    for bins in bins_list:
        grid = histogram(points, bins)
        met = metrics(grid)
        rows.append({"bins": bins, **met})
        render_grid([[float(v) for v in row] for row in grid], raw_dir / f"raw_bins_{bins:03d}.png", width, height)
        smoothed = smooth_grid(grid, args.sigma)
        render_grid(smoothed, soft_dir / f"soft_sigma{args.sigma:g}_bins_{bins:03d}.png", width, height)

    metrics_csv = outdir / "grid_resolution_scan_metrics.csv"
    with metrics_csv.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["bins", "occupied_cells", "occupied_fraction", "mean_count_occupied", "single_count_fraction_occupied", "max_count"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    useful = [row for row in rows if row["mean_count_occupied"] >= 2.0 and row["single_count_fraction_occupied"] <= 0.55]
    recommended = useful[-1]["bins"] if useful else rows[0]["bins"]
    subset_label = "All frames" if args.all_frames else f"C4-O1 <= {args.short_cutoff:.1f} A"
    readme = [
        "# Reactive CV Grid Resolution Scan",
        "",
        f"Subset: {subset_label}",
        f"Frames used: {len(points)}",
        "",
        "Raw images show true occupied-bin sampling. Soft images apply only sigma=1 bin smoothing.",
        "The vertical gray line marks 0 degrees for C1-O1-C4-C5(Ph).",
        "",
        f"Practical upper limit by mean occupied count >= 2 and singleton fraction <= 55%: about {recommended}x{recommended} bins.",
        "",
        "Interpretation: once `single_count_fraction_occupied` becomes high, finer bins mainly resolve individual saved frames rather than a stable density.",
        "",
        "bins, occupied_fraction, mean_count_occupied, single_count_fraction_occupied, max_count",
    ]
    for row in rows:
        readme.append(
            f"{row['bins']}, {row['occupied_fraction']:.6f}, {row['mean_count_occupied']:.3f}, "
            f"{row['single_count_fraction_occupied']:.3f}, {int(row['max_count'])}"
        )
    (outdir / "README_grid_resolution_scan.md").write_text("\n".join(readme), encoding="utf-8")
    print(os.fspath(outdir))


if __name__ == "__main__":
    main()
