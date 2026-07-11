#!/usr/bin/env python3
"""Make publication-style soft gray density plots for reactive CVs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec


plt.rcParams.update(
    {
        "font.size": 8,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "axes.linewidth": 0.8,
        "savefig.facecolor": "white",
    }
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--short-cutoff", type=float, default=3.0)
    parser.add_argument("--bins", type=int, default=260)
    parser.add_argument("--sigma", type=float, default=2.2)
    return parser.parse_args()


def gaussian_kernel1d(sigma: float, radius: int | None = None) -> np.ndarray:
    if sigma <= 0:
        return np.array([1.0])
    if radius is None:
        radius = max(3, int(round(4.0 * sigma)))
    x = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-(x * x) / (2.0 * sigma * sigma))
    return kernel / kernel.sum()


def smooth1d(values: np.ndarray, sigma: float) -> np.ndarray:
    kernel = gaussian_kernel1d(sigma)
    return np.convolve(values, kernel, mode="same")


def smooth2d(values: np.ndarray, sigma: float) -> np.ndarray:
    kernel = gaussian_kernel1d(sigma)
    tmp = np.apply_along_axis(lambda row: np.convolve(row, kernel, mode="same"), 0, values)
    return np.apply_along_axis(lambda row: np.convolve(row, kernel, mode="same"), 1, tmp)


def clean_xy(df: pd.DataFrame, xcol: str, ycol: str) -> pd.DataFrame:
    out = df[[xcol, ycol, "C4_O1_A"]].copy()
    out[xcol] = pd.to_numeric(out[xcol], errors="coerce")
    out[ycol] = pd.to_numeric(out[ycol], errors="coerce")
    out["C4_O1_A"] = pd.to_numeric(out["C4_O1_A"], errors="coerce")
    return out.replace([np.inf, -np.inf], np.nan).dropna()


def joint_density_plot(
    df: pd.DataFrame,
    xcol: str,
    ycol: str,
    out_png: Path,
    title: str,
    xlabel: str,
    ylabel: str,
    bins: int,
    sigma: float,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    short_cutoff: float | None = None,
) -> None:
    data = clean_xy(df, xcol, ycol)
    if short_cutoff is not None:
        data = data[data["C4_O1_A"] <= short_cutoff]
    if xlim is not None:
        data = data[(data[xcol] >= xlim[0]) & (data[xcol] <= xlim[1])]
    if ylim is not None:
        data = data[(data[ycol] >= ylim[0]) & (data[ycol] <= ylim[1])]

    x = data[xcol].to_numpy(float)
    y = data[ycol].to_numpy(float)
    if len(x) < 10:
        raise SystemExit(f"Not enough points for {out_png.name}: {len(x)}")

    hist, xedges, yedges = np.histogram2d(x, y, bins=bins, range=[xlim, ylim])
    smooth = smooth2d(hist.T, sigma)
    density = np.log10(smooth + 1.0)

    fig = plt.figure(figsize=(4.7, 4.55))
    grid = GridSpec(
        2,
        2,
        width_ratios=[4.0, 0.58],
        height_ratios=[0.58, 4.0],
        hspace=0.03,
        wspace=0.03,
    )
    ax_top = fig.add_subplot(grid[0, 0])
    ax = fig.add_subplot(grid[1, 0], sharex=ax_top)
    ax_right = fig.add_subplot(grid[1, 1], sharey=ax)

    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    ax.imshow(
        density,
        origin="lower",
        extent=extent,
        aspect="auto",
        cmap="Greys",
        interpolation="bicubic",
        vmin=0,
        vmax=np.nanpercentile(density, 99.5),
    )
    ax.scatter(x, y, s=0.15, c="black", alpha=0.010, linewidths=0)

    xhist, xhist_edges = np.histogram(x, bins=bins, range=xlim)
    yhist, yhist_edges = np.histogram(y, bins=bins, range=ylim)
    xcenters = 0.5 * (xhist_edges[:-1] + xhist_edges[1:])
    ycenters = 0.5 * (yhist_edges[:-1] + yhist_edges[1:])
    ax_top.plot(xcenters, smooth1d(xhist.astype(float), sigma), color="0.45", lw=1.25)
    ax_right.plot(smooth1d(yhist.astype(float), sigma), ycenters, color="0.45", lw=1.25)

    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.axvline(0.0, color="0.72", lw=0.55, ls="--")
    ax_top.axvline(0.0, color="0.78", lw=0.5, ls="--")
    ax_top.axis("off")
    ax_right.axis("off")
    fig.savefig(out_png, dpi=420, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.csv)

    df["abs_dihedral_C1_O1_C4_H02_deg"] = df["dihedral_C1_O1_C4_H02_deg"].abs()
    df["abs_dihedral_FE_N1_C1_C2_deg"] = df["dihedral_FE_N1_C1_C2_deg"].abs()
    df["abs_dihedral_C3_C4_O1_C1_deg"] = df["dihedral_C3_C4_O1_C1_deg"].abs()

    specs = [
        (
            "dihedral_C1_O1_C4_C5_deg",
            "abs_dihedral_C1_O1_C4_H02_deg",
            "soft_density_product_like_prochirality_all_frames.png",
            None,
            "All frames: product-like prochirality torsions",
            "Dihedral C1-O1-C4-C5(Ph) (deg)",
            "|Dihedral C1-O1-C4-H02| (deg)",
            (-180.0, 180.0),
            (0.0, 180.0),
        ),
        (
            "dihedral_C1_O1_C4_C5_deg",
            "abs_dihedral_C1_O1_C4_H02_deg",
            "soft_density_product_like_prochirality_C4O1_le_3A.png",
            args.short_cutoff,
            f"C4-O1 <= {args.short_cutoff:.1f} A: product-like prochirality torsions",
            "Dihedral C1-O1-C4-C5(Ph) (deg)",
            "|Dihedral C1-O1-C4-H02| (deg)",
            (-180.0, 180.0),
            (0.0, 180.0),
        ),
        (
            "dihedral_N1_C1_C2_C3_deg",
            "abs_dihedral_FE_N1_C1_C2_deg",
            "soft_density_literature_like_FeN_axis_all_frames.png",
            None,
            "All frames: literature-like Fe-N1/C1-C2 orientation",
            "Dihedral N1-C1-C2-C3 (deg)",
            "|Dihedral Fe-N1-C1-C2| (deg)",
            (-180.0, 180.0),
            (0.0, 180.0),
        ),
        (
            "dihedral_C1_O1_C4_C5_deg",
            "pseudo_CIP_O1_C5_C3_H02_chiral_volume_A3",
            "soft_density_C1O1C4C5_vs_pseudo_CIP_volume_all_frames.png",
            None,
            "All frames: signed torsion vs pseudo-CIP volume",
            "Dihedral C1-O1-C4-C5(Ph) (deg)",
            "Pseudo-CIP chiral volume (A^3)",
            (-180.0, 180.0),
            None,
        ),
        (
            "PHE322_ring_centroid_to_C4_A",
            "C4_radical_plane_normal_vs_PHE322_ring_normal_angle_deg",
            "soft_density_PHE322_pi_radical_geometry_all_frames.png",
            None,
            "All frames: PHE322 pi/radical geometry",
            "PHE322 centroid-C4 distance (A)",
            "C4 radical plane/PHE322 ring normal angle (deg)",
            None,
            (0.0, 90.0),
        ),
    ]
    for xcol, ycol, png, cutoff, title, xlabel, ylabel, xlim, ylim in specs:
        joint_density_plot(
            df=df,
            xcol=xcol,
            ycol=ycol,
            out_png=outdir / png,
            title=title,
            xlabel=xlabel,
            ylabel=ylabel,
            bins=args.bins,
            sigma=args.sigma,
            xlim=xlim,
            ylim=ylim,
            short_cutoff=cutoff,
        )

    short = df[df["C4_O1_A"] <= args.short_cutoff]
    lines = [
        "# Soft Gray Reactive CV Density Plots",
        "",
        "These plots use smoothed log-density plus very faint black points to mimic the soft gray cloud style in literature figures.",
        "They are density/occupancy views, not free-energy maps.",
        "",
        "Main product-like prochirality proposal:",
        "- x = `dihedral_C1_O1_C4_C5_deg`: signed phenyl orientation around the putative O1-C4 forming bond.",
        "- y = `abs(dihedral_C1_O1_C4_H02_deg)`: absolute H02 orientation around the same product-like axis.",
        "- The sign of x can be used as a pro-face descriptor only after calibration against an explicit product stereochemical model.",
        "",
        f"All frames: {len(df)}",
        f"C4-O1 <= {args.short_cutoff:.1f} A frames: {len(short)}",
        "",
        "Short-frame sign summary:",
        f"- median C1-O1-C4-C5 dihedral: {short['dihedral_C1_O1_C4_C5_deg'].median():.2f} deg",
        f"- fraction with C1-O1-C4-C5 > 0: {(short['dihedral_C1_O1_C4_C5_deg'] > 0).mean():.3f}",
        f"- median pseudo-CIP volume: {short['pseudo_CIP_O1_C5_C3_H02_chiral_volume_A3'].median():.3f} A^3",
        f"- fraction with pseudo-CIP volume > 0: {(short['pseudo_CIP_O1_C5_C3_H02_chiral_volume_A3'] > 0).mean():.3f}",
        "",
    ]
    (outdir / "README_soft_density.md").write_text("\n".join(lines), encoding="utf-8")
    print(outdir)


if __name__ == "__main__":
    main()
