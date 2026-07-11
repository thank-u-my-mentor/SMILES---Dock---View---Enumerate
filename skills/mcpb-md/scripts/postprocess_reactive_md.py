#!/usr/bin/env python3
"""Post-process MCPB/GROMACS MD for reactive-distance analysis.

Outputs:
- protein CA RMSD over time
- reactive atom-pair distance over time
- closest-frame PDB snapshots
- residue contact enrichment against short reactive-distance frames
- RMSD-metric clustering on sampled frames
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path

import MDAnalysis as mda
import mdtraj as md
import numpy as np
import pandas as pd
from MDAnalysis.analysis import rms
from scipy.stats import spearmanr
from sklearn.cluster import AgglomerativeClustering


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--top", required=True, help="Topology/structure file, e.g. .gro")
    parser.add_argument("--traj", required=True, help="Trajectory file, e.g. fitted .xtc")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--atom-a", default="resname UL1 and name C4")
    parser.add_argument("--atom-b", default="resname UL1 and name O1")
    parser.add_argument("--rmsd-select", default="protein and name CA")
    parser.add_argument("--contact-cutoff", type=float, default=4.0, help="Angstrom")
    parser.add_argument("--contact-stride", type=int, default=10, help="Trajectory frame stride for residue enrichment")
    parser.add_argument("--short-quantile", type=float, default=0.10)
    parser.add_argument("--long-quantile", type=float, default=0.50)
    parser.add_argument("--n-closest", type=int, default=20)
    parser.add_argument("--cluster-stride", type=int, default=100, help="Trajectory frame stride")
    parser.add_argument("--cluster-cutoff", type=float, default=2.0, help="CA RMSD cutoff in Angstrom")
    return parser.parse_args()


def require_one(selection, label: str):
    if len(selection) != 1:
        raise SystemExit(f"Selection {label!r} matched {len(selection)} atoms, expected 1")
    return selection[0]


def residue_label(residue) -> str:
    return f"{residue.resname}{residue.resid}"


def compute_rmsd(u: mda.Universe, select: str, outdir: Path) -> pd.DataFrame:
    if len(u.select_atoms(select)) == 0:
        raise SystemExit(f"RMSD selection matched no atoms: {select}")
    r = rms.RMSD(u, u, select=select, ref_frame=0)
    r.run()
    df = pd.DataFrame(
        {
            "frame": r.results.rmsd[:, 1].astype(int),
            "time_ps": r.results.rmsd[:, 1] * np.nan,
            "rmsd_A": r.results.rmsd[:, 2],
        }
    )
    # MDAnalysis RMSD stores [frame, time(ps), rmsd(A)].
    df["frame"] = r.results.rmsd[:, 0].astype(int)
    df["time_ps"] = r.results.rmsd[:, 1]
    df.to_csv(outdir / "protein_CA_RMSD.csv", index=False)
    return df


def compute_pair_distance(u: mda.Universe, atom_a, atom_b, outdir: Path) -> pd.DataFrame:
    rows = []
    for ts in u.trajectory:
        dist = float(np.linalg.norm(atom_a.position - atom_b.position))
        rows.append((int(ts.frame), float(ts.time), dist))
    df = pd.DataFrame(rows, columns=["frame", "time_ps", "C4_O1_distance_A"])
    df["time_ns"] = df["time_ps"] / 1000.0
    df.to_csv(outdir / "UL1_C4_O1_distance_timeseries.csv", index=False)
    return df


def write_closest_snapshots(u: mda.Universe, dist_df: pd.DataFrame, outdir: Path, n: int) -> pd.DataFrame:
    snapdir = outdir / "closest_C4_O1_frames"
    snapdir.mkdir(parents=True, exist_ok=True)
    closest = dist_df.nsmallest(n, "C4_O1_distance_A").copy()
    records = []
    for rank, row in enumerate(closest.itertuples(index=False), start=1):
        u.trajectory[int(row.frame)]
        filename = snapdir / f"rank{rank:02d}_frame{int(row.frame):05d}_t{row.time_ns:08.3f}ns_d{row.C4_O1_distance_A:05.2f}A.pdb"
        u.atoms.write(str(filename))
        records.append(
            {
                "rank": rank,
                "frame": int(row.frame),
                "time_ps": float(row.time_ps),
                "time_ns": float(row.time_ns),
                "C4_O1_distance_A": float(row.C4_O1_distance_A),
                "pdb": str(filename),
            }
        )
    out = pd.DataFrame(records)
    out.to_csv(outdir / "closest_C4_O1_frames.csv", index=False)
    return out


def residue_contact_enrichment(
    u: mda.Universe,
    atom_a,
    atom_b,
    dist_df: pd.DataFrame,
    outdir: Path,
    contact_cutoff: float,
    contact_stride: int,
    short_quantile: float,
    long_quantile: float,
) -> pd.DataFrame:
    short_thr = float(dist_df["C4_O1_distance_A"].quantile(short_quantile))
    long_thr = float(dist_df["C4_O1_distance_A"].quantile(long_quantile))
    dist_sample = dist_df.iloc[::contact_stride].reset_index(drop=True)
    rxn = dist_sample["C4_O1_distance_A"].to_numpy(dtype=float)
    short_mask = rxn <= short_thr
    long_mask = rxn >= long_thr

    residues = []
    for residue in u.residues:
        if residue.resname == "UL1":
            continue
        atoms = residue.atoms.select_atoms("not name H*")
        if len(atoms) == 0:
            continue
        residues.append((residue, atoms.indices))

    n_contact_frames = len(dist_sample)
    min_d = np.empty((len(residues), n_contact_frames), dtype=np.float32)
    pair_indices = np.array([atom_a.index, atom_b.index], dtype=int)

    for iframe, ts in enumerate(u.trajectory[::contact_stride]):
        pair_pos = u.atoms[pair_indices].positions
        for ires, (_, atom_indices) in enumerate(residues):
            pos = u.atoms[atom_indices].positions
            d1 = np.linalg.norm(pos - pair_pos[0], axis=1).min()
            d2 = np.linalg.norm(pos - pair_pos[1], axis=1).min()
            min_d[ires, iframe] = min(d1, d2)

    rows = []
    for ires, (residue, atom_indices) in enumerate(residues):
        vals = min_d[ires].astype(float)
        corr = spearmanr(vals, rxn, nan_policy="omit").correlation
        if corr is None or not math.isfinite(float(corr)):
            corr = 0.0
        contact = vals <= contact_cutoff
        freq_short = float(contact[short_mask].mean()) if short_mask.any() else 0.0
        freq_long = float(contact[long_mask].mean()) if long_mask.any() else 0.0
        mean_short = float(vals[short_mask].mean()) if short_mask.any() else float("nan")
        mean_long = float(vals[long_mask].mean()) if long_mask.any() else float("nan")
        enrichment = freq_short - freq_long
        score = enrichment * max(0.0, float(corr))
        rows.append(
            {
                "residue": residue_label(residue),
                "resname": residue.resname,
                "resid": residue.resid,
                "n_heavy_atoms": len(atom_indices),
                "contact_cutoff_A": contact_cutoff,
                "contact_stride_frames": contact_stride,
                "short_distance_threshold_A": short_thr,
                "long_distance_threshold_A": long_thr,
                "contact_freq_short": freq_short,
                "contact_freq_long": freq_long,
                "contact_enrichment_short_minus_long": enrichment,
                "mean_min_dist_short_A": mean_short,
                "mean_min_dist_long_A": mean_long,
                "spearman_residue_min_dist_vs_C4O1": float(corr),
                "proximity_score": score,
            }
        )
    df = pd.DataFrame(rows).sort_values(
        ["proximity_score", "contact_enrichment_short_minus_long"],
        ascending=[False, False],
    )
    df.to_csv(outdir / "residue_C4O1_short_distance_enrichment.csv", index=False)
    return df


def cluster_sampled_frames(top: str, traj: str, outdir: Path, stride: int, cutoff_A: float) -> pd.DataFrame:
    t = md.load(traj, top=top, stride=stride)
    atom_indices = t.topology.select("protein and name CA")
    if len(atom_indices) == 0:
        raise SystemExit("mdtraj clustering selection matched no CA atoms")
    n = t.n_frames
    dist_nm = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        dist_nm[i, :] = md.rmsd(t, t, frame=i, atom_indices=atom_indices)
    dist_A = dist_nm * 10.0

    model = AgglomerativeClustering(
        n_clusters=None,
        metric="precomputed",
        linkage="average",
        distance_threshold=cutoff_A,
    )
    labels = model.fit_predict(dist_A)
    rows = []
    for label in sorted(set(labels)):
        members = np.where(labels == label)[0]
        sub = dist_A[np.ix_(members, members)]
        medoid_local = int(np.argmin(sub.mean(axis=1)))
        medoid_sample = int(members[medoid_local])
        orig_frame = medoid_sample * stride
        rows.append(
            {
                "cluster_id": int(label),
                "n_sampled_frames": int(len(members)),
                "population_fraction": float(len(members) / n),
                "medoid_sample_index": medoid_sample,
                "medoid_original_frame": orig_frame,
                "medoid_time_ps": float(t.time[medoid_sample]),
                "medoid_time_ns": float(t.time[medoid_sample] / 1000.0),
                "mean_intra_cluster_RMSD_A": float(sub.mean()) if len(members) > 1 else 0.0,
            }
        )
    df = pd.DataFrame(rows).sort_values("n_sampled_frames", ascending=False)
    df.to_csv(outdir / "rmsd_metric_clusters_CA_stride.csv", index=False)
    if not df.empty:
        snapdir = outdir / "cluster_medoid_frames"
        snapdir.mkdir(parents=True, exist_ok=True)
        full = mda.Universe(top, traj)
        for row in df.head(10).itertuples(index=False):
            full.trajectory[int(row.medoid_original_frame)]
            fn = snapdir / f"cluster{int(row.cluster_id):03d}_pop{float(row.population_fraction):.3f}_frame{int(row.medoid_original_frame):05d}_t{float(row.medoid_time_ns):08.3f}ns.pdb"
            full.atoms.write(str(fn))
    return df


def write_summary(outdir: Path, dist_df: pd.DataFrame, rmsd_df: pd.DataFrame, closest: pd.DataFrame, enrich: pd.DataFrame, clusters: pd.DataFrame) -> None:
    q = dist_df["C4_O1_distance_A"].quantile([0, 0.01, 0.05, 0.10, 0.50, 0.90, 0.95, 0.99, 1.0])
    lines = []
    lines.append("# Reactive MD Postprocessing Summary\n")
    lines.append("## C4-O1 Distance Quantiles\n")
    lines.append(q.to_string())
    lines.append("\n\n## RMSD Summary\n")
    lines.append(rmsd_df["rmsd_A"].describe().to_string())
    lines.append("\n\n## Closest C4-O1 Frames\n")
    lines.append(closest.head(20).to_string(index=False))
    lines.append("\n\n## Top Residue Short-Distance Enrichment\n")
    cols = [
        "residue",
        "contact_freq_short",
        "contact_freq_long",
        "contact_enrichment_short_minus_long",
        "spearman_residue_min_dist_vs_C4O1",
        "proximity_score",
        "mean_min_dist_short_A",
        "mean_min_dist_long_A",
    ]
    lines.append(enrich[cols].head(30).to_string(index=False))
    lines.append("\n\n## RMSD-Metric Clusters\n")
    lines.append(clusters.head(20).to_string(index=False))
    lines.append("\n")
    (outdir / "SUMMARY.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    u = mda.Universe(args.top, args.traj)
    atom_a = require_one(u.select_atoms(args.atom_a), args.atom_a)
    atom_b = require_one(u.select_atoms(args.atom_b), args.atom_b)

    rmsd_df = compute_rmsd(u, args.rmsd_select, outdir)
    u.trajectory[0]
    dist_df = compute_pair_distance(u, atom_a, atom_b, outdir)
    u.trajectory[0]
    closest = write_closest_snapshots(u, dist_df, outdir, args.n_closest)
    u.trajectory[0]
    enrich = residue_contact_enrichment(
        u,
        atom_a,
        atom_b,
        dist_df,
        outdir,
        args.contact_cutoff,
        args.contact_stride,
        args.short_quantile,
        args.long_quantile,
    )
    clusters = cluster_sampled_frames(args.top, args.traj, outdir, args.cluster_stride, args.cluster_cutoff)
    write_summary(outdir, dist_df, rmsd_df, closest, enrich, clusters)
    print(outdir / "SUMMARY.md")


if __name__ == "__main__":
    main()
