#!/usr/bin/env python3
"""Refresh a product-calibrated pro-R/pro-S triple-product scatter plot."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import make_md_publication_style_summary as summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pro-rs-csv", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--near-cutoff", type=float, default=3.0)
    parser.add_argument("--arial-font", default="/mnt/c/Windows/Fonts/arial.ttf")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    pro = pd.read_csv(args.pro_rs_csv)
    pro["pro_rs_label"] = pro[summary.PRO_RS_COL].map(summary.PRO_RS_LABELS)
    summary.setup_style(args.arial_font)
    summary.plot_pro_rs_triple(pro, Path(args.outdir), args.near_cutoff)
    print(args.outdir)


if __name__ == "__main__":
    main()
