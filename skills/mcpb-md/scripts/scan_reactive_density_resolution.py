#!/usr/bin/env python3
"""Scan 2D density grid resolution for reactive MD CV plots."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--short-cutoff", type=float, default=3.0)
    parser.add_argument("--bins", default="60,80,100,120,160,200,240,300,360,480")
    parser.add_argument("--sigma", type=float, default=1.0)
    return parser.parse_args()


def gaussian_kernel1d(sigma: float) -> np.ndarray:
    if sigma <= 0:
        return np.array([1.0])
    radius = max(2, int(round(4.0 * sigma)))
    x = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-(x * x) / (2.0 * sigma * sigma))
    return kernel / kernel.sum()


def smooth2d(values: np.ndarray, sigma: float) -> np.ndarray:
    kernel = gaussian_kernel1d(sigma)
    tmp = np.apply_along_axis(lambda col: np.convolve(col, kernel, mode="same"), 0, values)
    return np.apply_along_axis(lambda row: np.convolve(row, kernel, mode="same"), 1, tmp)


def hist_metrics(counts: np.ndarray) -> dict[str, float]:
    occupied = counts[counts > 0]
    cells = counts.size
    if len(occupied) == 0:
        return {
            "occupied_cells": 0,
            "occupied_fraction": 0.0,
            "mean_count_occupied": 0.0,
            "single_count_fraction_occupied": 0.0,
            "max_count": 0.0,
        }
    return {
        "occupied_cells": int(len(occupied)),
        "occupied_fraction": float(len(occupied) / cells),
        "mean_count_occupied": float(occupied.mean()),
        "single_count_fraction_occupied": float((occupied == 1).mean()),
        "max_count": float(occupied.max()),
    }


def plot_scan(
    x: np.ndarray,
    y: np.ndarray,
    bins_list: list[int],
    out_png: Path,
    mode: str,
    sigma: float,
    xlim=(-180.0, 180.0),
    ylim=(0.0, 180.0),
) -> pd.DataFrame:
    ncols = 5
    nrows = math.ceil(len(bins_list) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.0 * ncols, 2.7 * nrows), sharex=True, sharey=True)
    axes = np.asarray(axes).reshape(-1)
    rows = []
    for ax, bins in zip(axes, bins_list):
        counts, xedges, yedges = np.histogram2d(x, y, bins=bins, range=[xlim, ylim])
        metrics = hist_metrics(counts)
        rows.append({"bins": bins, "mode": mode, **metrics})
        grid = counts.T
        if mode == "soft":
            grid = smooth2d(grid, sigma)
        image = np.log10(grid + 1.0)
        ax.imshow(
            image,
            origin="lower",
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
            cmap="Greys",
            aspect="auto",
            interpolation="nearest" if mode == "raw" else "bicubic",
            vmin=0,
            vmax=max(0.5, np.nanpercentile(image, 99.5)),
        )
        ax.axvline(0.0, color="0.72", lw=0.5, ls="--")
        ax.set_title(
            f"{bins}x{bins}\nocc {metrics['occupied_fraction']*100:.1f}%, "
            f"mean {metrics['mean_count_occupied']:.1f}, "
            f"1-count {metrics['single_count_fraction_occupied']*100:.0f}%",
            fontsize=8,
        )
        ax.tick_params(labelsize=7)
    for ax in axes[len(bins_list) :]:
        ax.axis("off")
    fig.supxlabel("Dihedral C1-O1-C4-C5(Ph) (deg)", fontsize=10)
    fig.supylabel("|Dihedral C1-O1-C4-H02| (deg)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=260)
    plt.close(fig)
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.csv)
    short = df[df["C4_O1_A"] <= args.short_cutoff].copy()
    short["abs_dihedral_C1_O1_C4_H02_deg"] = short["dihedral_C1_O1_C4_H02_deg"].abs()
    data = short[["dihedral_C1_O1_C4_C5_deg", "abs_dihedral_C1_O1_C4_H02_deg"]].replace([np.inf, -np.inf], np.nan).dropna()
    x = data["dihedral_C1_O1_C4_C5_deg"].to_numpy(float)
    y = data["abs_dihedral_C1_O1_C4_H02_deg"].to_numpy(float)

    bins_list = [int(v.strip()) for v in args.bins.split(",") if v.strip()]
    raw = plot_scan(
        x,
        y,
        bins_list,
        outdir / "grid_resolution_scan_raw_counts_C4O1_le_3A.png",
        mode="raw",
        sigma=args.sigma,
    )
    soft = plot_scan(
        x,
        y,
        bins_list,
        outdir / "grid_resolution_scan_soft_sigma1_C4O1_le_3A.png",
        mode="soft",
        sigma=args.sigma,
    )
    table = pd.concat([raw, soft], ignore_index=True)
    table.to_csv(outdir / "grid_resolution_scan_metrics.csv", index=False)

    useful = raw[(raw["mean_count_occupied"] >= 2.0) & (raw["single_count_fraction_occupied"] <= 0.55)]
    recommended = int(useful["bins"].max()) if len(useful) else int(raw.iloc[0]["bins"])
    lines = [
        "# Reactive CV Grid Resolution Scan",
        "",
        f"Subset: C4-O1 <= {args.short_cutoff:.1f} A",
        f"Frames used: {len(data)}",
        "",
        "The raw-count panel exposes the true sampling limit. Once most occupied cells contain only one frame, increasing bins mostly adds speckle rather than information.",
        "",
        f"Practical upper limit by a conservative criterion mean occupied count >= 2 and singleton fraction <= 55%: about {recommended}x{recommended} bins.",
        "",
        "Columns in `grid_resolution_scan_metrics.csv`:",
        "- `occupied_fraction`: fraction of all 2D grid cells containing at least one frame.",
        "- `mean_count_occupied`: mean frame count among non-empty cells.",
        "- `single_count_fraction_occupied`: fraction of non-empty cells with exactly one frame.",
        "",
        raw.to_string(index=False),
        "",
    ]
    (outdir / "README_grid_resolution_scan.md").write_text("\n".join(lines), encoding="utf-8")
    print(outdir)


if __name__ == "__main__":
    main()
