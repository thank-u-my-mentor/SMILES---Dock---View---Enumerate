#!/usr/bin/env python3
"""
End-to-end wrapper:
seed FASTA -> HMMER homologs -> SSN/MMseqs representatives -> IQ-TREE/iTOL files.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import List, Optional


def run(cmd: List[str], cwd: Optional[Path] = None) -> None:
    print("\n$ " + " ".join(cmd))
    result = subprocess.run(cmd, cwd=str(cwd) if cwd else None, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(cmd)}")


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Run the hydrolase homolog-tree workflow.")
    parser.add_argument("seed_fasta", help="Curated characterized hydrolase seed FASTA.")
    parser.add_argument("--output", default="hydrolase_homolog_tree_run")
    parser.add_argument("--hmmer-mode", choices=["auto", "web-phmmer", "local-hmmsearch"], default="auto")
    parser.add_argument("--db-fasta", default=None, help="Optional local protein FASTA database for profile HMM search.")
    parser.add_argument("--max-hits", type=int, default=2500)
    parser.add_argument("--max-total", type=int, default=3000)
    parser.add_argument("--max-rep", type=int, default=2500)
    parser.add_argument("--cdhit-identity", type=float, default=0.35)
    parser.add_argument("--cdhit-coverage", type=float, default=0.75)
    parser.add_argument("--max-recon", type=int, default=3000)
    parser.add_argument("--ssn-script", default=str(script_dir / "ssn_pipeline_mmseq2.py"))
    parser.add_argument("--tree-script", default=str(script_dir / "tree_pipeline.py"))
    parser.add_argument("--solubility-csv", default=None, help="Optional NetSolP/SoluProt CSV to convert into iTOL annotation.")
    args = parser.parse_args()

    output = Path(args.output).resolve()
    hmmer_dir = output / "01_hmmer"
    ssn_dir = output / "02_ssn"
    tree_dir = ssn_dir / "tree_analysis"
    output.mkdir(parents=True, exist_ok=True)

    fetch_cmd = [
        sys.executable,
        str(script_dir / "fetch_hmmer_homologs.py"),
        str(Path(args.seed_fasta).resolve()),
        "--outdir", str(hmmer_dir),
        "--mode", args.hmmer_mode,
        "--max-hits", str(args.max_hits),
        "--max-total", str(args.max_total),
    ]
    if args.db_fasta:
        fetch_cmd.extend(["--db-fasta", str(Path(args.db_fasta).resolve())])
    run(fetch_cmd)

    combined_fasta = hmmer_dir / "homologs_plus_core.fasta"
    core_fasta = hmmer_dir / "core.fasta"
    ssn_cmd = [
        sys.executable,
        str(Path(args.ssn_script).resolve()),
        str(combined_fasta),
        "--core-fasta", str(core_fasta),
        "--output", str(ssn_dir),
        "--max-rep", str(args.max_rep),
        "--cdhit-identity", str(args.cdhit_identity),
        "--cdhit-coverage", str(args.cdhit_coverage),
        "--max-recon", str(args.max_recon),
    ]
    run(ssn_cmd)

    representatives = ssn_dir / "representatives.fasta"
    if not representatives.exists():
        raise FileNotFoundError(f"Missing representatives FASTA: {representatives}")

    tree_cmd = [
        sys.executable,
        str(Path(args.tree_script).resolve()),
        str(representatives),
        "--output", str(tree_dir),
    ]
    # tree_pipeline.py intentionally reads ssn_nodes.csv from current working directory.
    run(tree_cmd, cwd=ssn_dir)

    if args.solubility_csv:
        sol_cmd = [
            sys.executable,
            str(script_dir / "annotate_solubility_itol.py"),
            str(Path(args.solubility_csv).resolve()),
            "--nodes", str(ssn_dir / "ssn_nodes.csv"),
            "--outdir", str(tree_dir),
            "--tool-label", "Solubility",
        ]
        run(sol_cmd)

    print("\nDone.")
    print(f"Combined homolog FASTA: {combined_fasta}")
    print(f"SSN nodes: {ssn_dir / 'ssn_nodes.csv'}")
    print(f"SSN XGMML: {ssn_dir / 'ssn_network.xgmml'}")
    print(f"Tree: {tree_dir / 'phylogenetic_tree.newick'}")
    print(f"iTOL files: {tree_dir}")


if __name__ == "__main__":
    main()
