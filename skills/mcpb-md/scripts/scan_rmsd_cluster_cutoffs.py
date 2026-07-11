#!/usr/bin/env python3
"""Scan RMSD-metric agglomerative clustering cutoffs for sampled MD frames."""

from __future__ import annotations

import argparse
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--top", required=True)
    parser.add_argument("--traj", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--stride", type=int, default=100)
    parser.add_argument("--cutoffs", default="0.8,1.0,1.2,1.5,2.0", help="Angstrom, comma-separated")
    parser.add_argument("--select", default="protein and name CA")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cutoffs = [float(x) for x in args.cutoffs.split(",") if x.strip()]

    t = md.load(args.traj, top=args.top, stride=args.stride)
    atom_indices = t.topology.select(args.select)
    if len(atom_indices) == 0:
        raise SystemExit(f"Selection matched no atoms: {args.select}")

    n = t.n_frames
    dist_nm = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        dist_nm[i, :] = md.rmsd(t, t, frame=i, atom_indices=atom_indices)
    dist_A = dist_nm * 10.0

    rows = []
    for cutoff in cutoffs:
        model = AgglomerativeClustering(
            n_clusters=None,
            metric="precomputed",
            linkage="average",
            distance_threshold=cutoff,
        )
        labels = model.fit_predict(dist_A)
        counts = np.array([np.sum(labels == lab) for lab in sorted(set(labels))], dtype=int)
        largest_label = sorted(set(labels))[int(np.argmax(counts))]
        members = np.where(labels == largest_label)[0]
        sub = dist_A[np.ix_(members, members)]
        medoid = int(members[int(np.argmin(sub.mean(axis=1)))])
        rows.append(
            {
                "cutoff_A": cutoff,
                "n_clusters": int(len(counts)),
                "largest_cluster_size": int(counts.max()),
                "largest_cluster_fraction": float(counts.max() / n),
                "largest_cluster_medoid_sample_index": medoid,
                "largest_cluster_medoid_original_frame": int(medoid * args.stride),
                "largest_cluster_medoid_time_ps": float(t.time[medoid]),
                "largest_cluster_medoid_time_ns": float(t.time[medoid] / 1000.0),
            }
        )
    df = pd.DataFrame(rows)
    df.to_csv(outdir / "rmsd_cluster_cutoff_scan.csv", index=False)
    print(outdir / "rmsd_cluster_cutoff_scan.csv")


if __name__ == "__main__":
    main()
