#!/usr/bin/env python3
"""Build the optional iTOL annotation layer for an existing tree_analysis run."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional


def run(cmd: List[str]) -> None:
    print("\n$ " + " ".join(cmd))
    result = subprocess.run(cmd, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(cmd)}")


def copy_if_exists(src: Path, dst: Path) -> Optional[Path]:
    if not src.exists():
        return None
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    print(f"Copied {src} -> {dst}")
    return dst


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Build consolidated hydrolase tree iTOL annotations.")
    parser.add_argument("--tree-dir", required=True, help="tree_analysis directory containing phylogenetic_tree.newick.")
    parser.add_argument("--nodes", required=True, help="ssn_nodes.csv from the SSN output directory.")
    parser.add_argument("--base-xlsx", default=None, help="Optional curated base.xlsx for candidate/core-label annotations.")
    parser.add_argument("--soluprot-csv", default=None, help="Optional user-downloaded SoluProt/NetSolP CSV, usually input.csv.")
    parser.add_argument("--id-mapping", default=None, help="Optional soluprot_id_mapping.tsv from prepare_soluprot_tree_inputs.py.")
    parser.add_argument("--low-threshold", type=float, default=0.2)
    parser.add_argument("--tool-label", default="SoluProt")
    args = parser.parse_args()

    tree_dir = Path(args.tree_dir).resolve()
    nodes = Path(args.nodes).resolve()
    annotations_dir = tree_dir / "annotations"
    annotations_dir.mkdir(parents=True, exist_ok=True)

    top_level_files: List[Path] = []
    base_out = annotations_dir / "base_xlsx"
    if args.base_xlsx:
        run([
            sys.executable,
            str(script_dir / "annotate_base_xlsx_itol.py"),
            str(Path(args.base_xlsx).resolve()),
            "--nodes", str(nodes),
            "--outdir", str(base_out),
        ])
        copied = copy_if_exists(base_out / "itol_base_candidate_highlight.txt", tree_dir / "itol_base_candidate_highlight.txt")
        if copied:
            top_level_files.append(copied)
        copied = copy_if_exists(base_out / "itol_core_short_name_text.txt", tree_dir / "itol_core_short_name_text.txt")
        if copied:
            top_level_files.append(copied)

    sol_out = annotations_dir / "soluprot"
    if args.soluprot_csv:
        cmd = [
            sys.executable,
            str(script_dir / "annotate_solubility_itol.py"),
            str(Path(args.soluprot_csv).resolve()),
            "--nodes", str(nodes),
            "--outdir", str(sol_out),
            "--tool-label", args.tool_label,
            "--low-threshold", str(args.low_threshold),
            "--symbols-only",
            "--no-color-strip",
        ]
        if args.id_mapping:
            cmd.extend(["--id-mapping", str(Path(args.id_mapping).resolve())])
        run(cmd)
        copied = copy_if_exists(sol_out / "itol_soluprot_gradient_symbols.txt", tree_dir / "itol_soluprot_gradient_symbols.txt")
        if copied:
            top_level_files.append(copied)

    print("\nAnnotation output:")
    print(f"  audit/intermediate directory: {annotations_dir}")
    for path in top_level_files:
        print(f"  iTOL upload file: {path}")


if __name__ == "__main__":
    main()
