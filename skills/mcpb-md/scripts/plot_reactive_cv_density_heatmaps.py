#!/usr/bin/env python3
"""Plot raw-density heatmaps for reactive CVs without free-energy conversion."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--short-cutoff", type=float, default=3.0)
    return parser.parse_args()


def plot_density(df: pd.DataFrame, xcol: str, ycol: str, out_png: Path, title: str, short_cutoff: float, bins=90):
    data = pd.DataFrame(
        {
            xcol: df[xcol],
            ycol: df[ycol],
            "_short_marker_C4_O1_A": df["C4_O1_A"],
        }
    ).replace([np.inf, -np.inf], np.nan).dropna()
    x = data[xcol].to_numpy(float)
    y = data[ycol].to_numpy(float)
    counts, xedges, yedges = np.histogram2d(x, y, bins=bins)
    log_counts = np.log10(counts.T + 1.0)

    fig, ax = plt.subplots(figsize=(6.7, 5.4))
    im = ax.pcolormesh(xedges, yedges, log_counts, cmap="magma", shading="auto")
    short = data[data["_short_marker_C4_O1_A"] <= short_cutoff]
    if len(short):
        ax.scatter(short[xcol], short[ycol], s=2.0, c="#33ccff", alpha=0.20, linewidths=0, label=f"C4-O1 <= {short_cutoff:.1f} A")
        ax.legend(frameon=False, loc="best", fontsize=8)
    ax.set_xlabel(xcol)
    ax.set_ylabel(ycol)
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("log10(frame count + 1)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=260)
    plt.close(fig)


def plot_1d(df: pd.DataFrame, col: str, out_png: Path, title: str, short_cutoff: float):
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    ax.hist(df[col].dropna(), bins=90, color="#777777", alpha=0.72, label="all frames")
    short = df[df["C4_O1_A"] <= short_cutoff]
    ax.hist(short[col].dropna(), bins=70, color="#33ccff", alpha=0.55, label=f"C4-O1 <= {short_cutoff:.1f} A")
    ax.set_xlabel(col)
    ax.set_ylabel("Frame count")
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=240)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.csv)
    df["abs_dihedral_FE_N1_C1_C2_deg"] = df["dihedral_FE_N1_C1_C2_deg"].abs()

    specs = [
        (
            "dihedral_N1_C1_C2_C3_deg",
            "abs_dihedral_FE_N1_C1_C2_deg",
            "density_abs_FeN1C1C2_vs_N1C1C2C3.png",
            "Density: |Fe-N1-C1-C2| vs N1-C1-C2-C3",
        ),
        (
            "C4_O1_A",
            "C4_face_signed_height_A",
            "density_C4O1_vs_O1_face_signed_height.png",
            "Density: C4-O1 distance vs O1 attack face",
        ),
        (
            "C4_O1_A",
            "H02_face_signed_height_A",
            "density_C4O1_vs_H02_face_signed_height.png",
            "Density: C4-O1 distance vs H02 face",
        ),
        (
            "C4_O1_A",
            "O1_H02_same_face_product_A2",
            "density_C4O1_vs_O1_H02_same_face_product.png",
            "Density: C4-O1 distance vs O1/H02 same-face product",
        ),
        (
            "C4_O1_A",
            "pseudo_CIP_O1_C5_C3_H02_chiral_volume_A3",
            "density_C4O1_vs_pseudo_CIP_chiral_volume.png",
            "Density: C4-O1 distance vs product-like pseudo-CIP chiral volume",
        ),
        (
            "C4_O1_A",
            "lactone_plane_H02_C5_same_side_product_A2",
            "density_C4O1_vs_lactone_plane_H02_C5_product.png",
            "Density: C4-O1 distance vs H/Ph sides of C4-O1-C1 plane",
        ),
        (
            "dihedral_C1_O1_C4_C5_deg",
            "dihedral_C1_O1_C4_H02_deg",
            "density_C1O1C4C5_vs_C1O1C4H02.png",
            "Density: product-like O1-C4 dihedrals for Ph and H02",
        ),
        (
            "dihedral_C1_O1_C4_C5_deg",
            "pseudo_CIP_O1_C5_C3_H02_chiral_volume_A3",
            "density_C1O1C4C5_vs_pseudo_CIP_volume.png",
            "Density: Ph-side O1-C4 dihedral vs pseudo-CIP volume",
        ),
        (
            "C4_O1_A",
            "PHE322_ring_centroid_to_C4_A",
            "density_C4O1_vs_PHE322_ring_centroid_C4.png",
            "Density: C4-O1 distance vs PHE322 centroid-C4",
        ),
        (
            "PHE322_C4_abs_height_to_ring_A",
            "PHE322_C4_lateral_offset_A",
            "density_PHE322_height_vs_lateral_offset.png",
            "Density: C4 position relative to PHE322 ring plane",
        ),
        (
            "PHE322_ring_centroid_to_C4_A",
            "C4_radical_plane_normal_vs_PHE322_ring_normal_angle_deg",
            "density_PHE322_centroid_C4_vs_radical_ring_normal_angle.png",
            "Density: PHE322-C4 distance vs radical/ring normal alignment",
        ),
    ]
    for xcol, ycol, png, title in specs:
        plot_density(df, xcol, ycol, outdir / png, title, args.short_cutoff)

    plot_1d(df, "O1_H02_same_face_product_A2", outdir / "hist_O1_H02_same_face_product.png", "O1/H02 face relation", args.short_cutoff)
    plot_1d(df, "pseudo_CIP_O1_C5_C3_H02_chiral_volume_A3", outdir / "hist_pseudo_CIP_chiral_volume.png", "Product-like pseudo-CIP chiral volume", args.short_cutoff)
    plot_1d(df, "lactone_plane_H02_C5_same_side_product_A2", outdir / "hist_lactone_plane_H02_C5_product.png", "H02/Ph side relation to C4-O1-C1 plane", args.short_cutoff)
    plot_1d(df, "C4_face_signed_height_A", outdir / "hist_O1_attack_face_signed_height.png", "O1 attack face signed height", args.short_cutoff)

    short = df[df["C4_O1_A"] <= args.short_cutoff]
    rows = [
        "# Reactive CV Density Heatmaps",
        "",
        "These plots show raw trajectory occupancy, not -kT ln(P) free energy.",
        "Color scale is `log10(frame count + 1)`. Cyan points mark short C4-O1 frames.",
        "",
        f"Short-frame cutoff: C4-O1 <= {args.short_cutoff:.3f} A",
        f"All frames: {len(df)}",
        f"Short frames: {len(short)}",
        "",
        "Useful interpretation:",
        "- `O1_H02_same_face_product_A2 > 0`: O1 and H02 lie on the same signed side of the C3-C4-C5 local plane.",
        "- `O1_H02_same_face_product_A2 < 0`: O1 and H02 lie on opposite signed sides.",
        "- `pseudo_CIP_O1_C5_C3_H02_chiral_volume_A3`: signed product-like C4 chirality volume using tentative priority O1 > C5(Ph) > C3(chain) > H02.",
        "- `lactone_plane_H02_C5_same_side_product_A2`: H02 and C5(Ph) side relation relative to the C4-O1-C1 plane proposed for product-like ring closure.",
        "- The sign convention is geometric; assign pro-R/pro-S only after calibration against a product stereochemical model.",
        "",
        "Short-frame summary:",
        short[
            [
                "C4_O1_A",
                "C4_face_signed_height_A",
                "H02_face_signed_height_A",
                "O1_H02_same_face_product_A2",
                "pseudo_CIP_O1_C5_C3_H02_chiral_volume_A3",
                "lactone_plane_H02_C5_same_side_product_A2",
                "dihedral_C1_O1_C4_C5_deg",
                "dihedral_C1_O1_C4_H02_deg",
                "PHE322_ring_centroid_to_C4_A",
                "PHE322_C4_abs_height_to_ring_A",
                "PHE322_C4_lateral_offset_A",
                "C4_radical_plane_normal_vs_PHE322_ring_normal_angle_deg",
            ]
        ]
        .describe(percentiles=[0.1, 0.5, 0.9])
        .to_string(),
        "",
    ]
    (outdir / "README_density_heatmaps.md").write_text("\n".join(rows), encoding="utf-8")
    print(outdir)


if __name__ == "__main__":
    main()
