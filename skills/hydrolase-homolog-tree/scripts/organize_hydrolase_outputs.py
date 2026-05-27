#!/usr/bin/env python3
"""Move hydrolase workflow intermediates into work/ folders without deleting them."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from typing import Iterable


def ensure_inside(parent: Path, child: Path) -> None:
    parent = parent.resolve()
    child = child.resolve()
    if parent != child and parent not in child.parents:
        raise ValueError(f"Refusing to move path outside {parent}: {child}")


def move_existing(paths: Iterable[Path], dest: Path, root: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    ensure_inside(root, dest)
    for src in paths:
        if not src.exists():
            continue
        ensure_inside(root, src)
        target = dest / src.name
        if target.exists():
            continue
        print(f"Move {src} -> {target}")
        shutil.move(str(src), str(target))


def main() -> None:
    parser = argparse.ArgumentParser(description="Organize SSN/tree intermediates while keeping upload files visible.")
    parser.add_argument("ssn_dir", help="SSN output directory, e.g. hydrolase_repro_ssn.")
    args = parser.parse_args()

    ssn_dir = Path(args.ssn_dir).resolve()
    tree_dir = ssn_dir / "tree_analysis"
    if not tree_dir.exists():
        raise FileNotFoundError(tree_dir)

    move_existing(
        [
            ssn_dir / "filtered_for_cluster.fasta",
            ssn_dir / "mmseqs_avA.m8",
            ssn_dir / "mmseqs_cluster_all_seqs.fasta",
            ssn_dir / "mmseqs_cluster_cluster.tsv",
            ssn_dir / "mmseqs_cluster_rep_seq.fasta",
            ssn_dir / "representatives_aln.fasta",
            ssn_dir / "mmseqs_tmp",
            ssn_dir / "mmseqs_tmp_avA",
        ],
        ssn_dir / "work" / "ssn_intermediates",
        ssn_dir,
    )

    move_existing(
        [tree_dir / "alignment.fasta", tree_dir / "input_sequences.fasta"],
        tree_dir / "work" / "msa",
        ssn_dir,
    )

    move_existing(
        list(tree_dir.glob("iqtree_run.*")),
        tree_dir / "work" / "iqtree",
        ssn_dir,
    )

    print("Organized hydrolase output folders.")


if __name__ == "__main__":
    main()
