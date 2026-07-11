#!/usr/bin/env python3
"""Plot outputs from postprocess_reactive_md.py."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-dir", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    d = Path(args.analysis_dir)
    dist = pd.read_csv(d / "UL1_C4_O1_distance_timeseries.csv")
    rmsd = pd.read_csv(d / "protein_CA_RMSD.csv")
    enrich = pd.read_csv(d / "residue_C4O1_short_distance_enrichment.csv")

    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(dist["time_ns"], dist["C4_O1_distance_A"], lw=0.55, color="#1f77b4")
    ax.axhline(dist["C4_O1_distance_A"].quantile(0.10), color="#d62728", ls="--", lw=1.0, label="10% quantile")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("UL1 C4-O1 distance (A)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(d / "plot_C4_O1_distance_timeseries.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(rmsd["time_ps"] / 1000.0, rmsd["rmsd_A"], lw=0.6, color="#2ca02c")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Protein CA RMSD (A)")
    fig.tight_layout()
    fig.savefig(d / "plot_protein_CA_RMSD.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.hist(dist["C4_O1_distance_A"], bins=80, color="#4c78a8", alpha=0.85)
    ax.axvline(dist["C4_O1_distance_A"].quantile(0.10), color="#d62728", ls="--", lw=1.0)
    ax.set_xlabel("UL1 C4-O1 distance (A)")
    ax.set_ylabel("Frame count")
    fig.tight_layout()
    fig.savefig(d / "plot_C4_O1_distance_histogram.png", dpi=220)
    plt.close(fig)

    top = enrich.head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.barh(top["residue"], top["contact_enrichment_short_minus_long"], color="#9467bd")
    ax.set_xlabel("Contact enrichment in short-distance frames")
    ax.set_ylabel("Residue")
    fig.tight_layout()
    fig.savefig(d / "plot_top_residue_contact_enrichment.png", dpi=220)
    plt.close(fig)

    print(d)


if __name__ == "__main__":
    main()
