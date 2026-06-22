#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional


SCRIPT_DIR = Path(__file__).resolve().parent
HYDROLASE_SCRIPTS = Path("/mnt/e/Codex/skills/hydrolase-homolog-tree/scripts")
DEFAULT_CONDA_BIN = Path("/mnt/l/WSL/softwares/conda_envs/md/bin")


def env_with_tools() -> Dict[str, str]:
    env = dict(os.environ)
    parts = [str(DEFAULT_CONDA_BIN), env.get("PATH", "")]
    env["PATH"] = ":".join(p for p in parts if p)
    env.pop("KIMI_API_KEY", None)
    return env


def run(cmd: List[str], cwd: Optional[Path] = None, env: Optional[Dict[str, str]] = None) -> None:
    print("\n$ " + " ".join(cmd))
    result = subprocess.run(cmd, cwd=str(cwd) if cwd else None, env=env, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(cmd)}")


def require_tool(name: str, env: Dict[str, str]) -> Optional[str]:
    found = shutil.which(name, path=env.get("PATH"))
    return found


def write_simple_tree_from_core(core_fasta: Path, outdir: Path, env: Dict[str, str], threads: int) -> None:
    tree_dir = outdir / "tree_analysis"
    tree_dir.mkdir(parents=True, exist_ok=True)
    aln = tree_dir / "alignment.fasta"
    tree = tree_dir / "phylogenetic_tree.newick"
    clustalo = require_tool("clustalo", env)
    iqtree = require_tool("iqtree", env) or require_tool("iqtree2", env)
    if not clustalo:
        raise RuntimeError("clustalo was not found. Install it or add it to PATH.")
    if not iqtree:
        raise RuntimeError("iqtree/iqtree2 was not found. Install it or add it to PATH.")
    run([clustalo, "-i", str(core_fasta), "-o", str(aln), "--outfmt", "fasta", "--force", "--threads", str(threads)], env=env)
    run([
        iqtree,
        "-s", str(aln),
        "-m", "LG+F+R4",
        "-B", "1000",
        "-T", str(threads),
        "-redo",
        "--prefix", str(tree_dir / "iqtree_run"),
    ], env=env)
    produced = tree_dir / "iqtree_run.treefile"
    if produced.exists():
        shutil.copyfile(produced, tree)
    else:
        raise FileNotFoundError(f"Missing IQ-TREE output: {produced}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Fast Metal-F project to phylogenetic tree workflow.")
    parser.add_argument("--metal-project", default="/home/qin/Metal-F_project")
    parser.add_argument("--outdir", default="/mnt/e/Tree-Metal-F")
    parser.add_argument("--min-length", type=int, default=120)
    parser.add_argument("--exclude-pdb", nargs="*", default=[])
    parser.add_argument("--mode", choices=["core-tree", "hmmer-tree", "existing-all-tree"], default="core-tree")
    parser.add_argument("--all-fasta", default="", help="Existing homolog/all FASTA for --mode existing-all-tree.")
    parser.add_argument("--hmmer-mode", choices=["ebi-hmmsearch", "web-phmmer", "local-hmmsearch"], default="ebi-hmmsearch")
    parser.add_argument("--hmmer-database", default="uniprot")
    parser.add_argument("--hmmer-evalue", default="1e-10")
    parser.add_argument("--db-fasta", default="")
    parser.add_argument("--max-hits", type=int, default=800)
    parser.add_argument("--max-total", type=int, default=1200)
    parser.add_argument("--max-rep", type=int, default=200)
    parser.add_argument("--cdhit-identity", type=float, default=0.35)
    parser.add_argument("--cdhit-coverage", type=float, default=0.75)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    env = env_with_tools()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    run([
        sys.executable,
        str(SCRIPT_DIR / "build_metal_f_tree_core.py"),
        "--metal-project", args.metal_project,
        "--outdir", str(outdir),
        "--min-length", str(args.min_length),
        "--exclude-pdb", *args.exclude_pdb,
    ], env=env)
    core_fasta = outdir / "tree_core.fasta"

    run([
        sys.executable,
        str(HYDROLASE_SCRIPTS / "construct_sequence_qc.py"),
        str(core_fasta),
        "--sequence-type", "protein",
        "--output", str(outdir / "tree_core.construct_qc.csv"),
    ], env=env)

    if args.mode == "core-tree":
        ssn_dir = outdir / "ssn"
        ssn_dir.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(core_fasta, ssn_dir / "representatives.fasta")
        write_simple_tree_from_core(ssn_dir / "representatives.fasta", outdir, env, args.threads)
    else:
        hmmer_dir = outdir / "hmmer"
        if args.mode == "existing-all-tree":
            if not args.all_fasta:
                raise ValueError("--all-fasta is required with --mode existing-all-tree.")
            hmmer_dir.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(Path(args.all_fasta), hmmer_dir / "homologs_plus_core.fasta")
            shutil.copyfile(core_fasta, hmmer_dir / "core.fasta")
        else:
            if args.hmmer_mode == "ebi-hmmsearch":
                run([
                    sys.executable,
                    str(SCRIPT_DIR / "fetch_hmmer_hmmsearch.py"),
                    str(core_fasta),
                    "--outdir", str(hmmer_dir),
                    "--database", args.hmmer_database,
                    "--evalue", args.hmmer_evalue,
                    "--max-hits", str(args.max_hits),
                    "--max-total", str(args.max_total),
                    "--threads", str(args.threads),
                ], env=env)
            else:
                fetch_cmd = [
                    sys.executable,
                    str(HYDROLASE_SCRIPTS / "fetch_hmmer_homologs.py"),
                    str(core_fasta),
                    "--outdir", str(hmmer_dir),
                    "--mode", args.hmmer_mode,
                    "--max-hits", str(args.max_hits),
                    "--max-total", str(args.max_total),
                ]
                if args.db_fasta:
                    fetch_cmd.extend(["--db-fasta", args.db_fasta])
                run(fetch_cmd, env=env)

        ssn_dir = outdir / "ssn"
        run([
            sys.executable,
            str(SCRIPT_DIR / "select_representatives_mmseqs.py"),
            str(hmmer_dir / "homologs_plus_core.fasta"),
            "--core-fasta", str(hmmer_dir / "core.fasta"),
            "--output", str(ssn_dir),
            "--metadata", str(hmmer_dir / "homolog_metadata.csv"), str(outdir / "tree_core_metadata.csv"),
            "--max-rep", str(args.max_rep),
            "--min-seq-id", str(args.cdhit_identity),
            "--coverage", str(args.cdhit_coverage),
            "--threads", str(args.threads),
        ], env=env)
        write_simple_tree_from_core(ssn_dir / "representatives.fasta", outdir, env, args.threads)

    run([
        sys.executable,
        str(SCRIPT_DIR / "offline_kingdom_itol.py"),
        "--nodes", str(outdir / "ssn" / "ssn_nodes.csv"),
        "--fasta", str(outdir / "ssn" / "representatives.fasta"),
        "--seed-metadata", str(outdir / "tree_core_metadata.csv"),
        "--outdir", str(outdir / "tree_analysis"),
    ], env=env)
    run([
        sys.executable,
        str(SCRIPT_DIR / "prepare_soluprot_tree_inputs.py"),
        str(outdir / "ssn" / "representatives.fasta"),
        "--outdir", str(outdir / "soluprot_inputs"),
    ], env=env)

    (outdir / "README_tree_workflow_outputs.txt").write_text(
        "Metal-F tree workflow outputs:\n"
        "- tree_core.fasta: broad local Metal-F seed FASTA, not geometry-filtered down to 7 records.\n"
        "- tree_core_metadata.csv: source PDB/UniProt/organism/kingdom metadata.\n"
        "- ssn/representatives.fasta: sequences used for the tree.\n"
        "- tree_analysis/alignment.fasta: ClustalO alignment.\n"
        "- tree_analysis/phylogenetic_tree.newick: tree for iTOL.\n"
        "- tree_analysis/itol_kingdom_color_strip_offline.txt: offline kingdom color strip.\n"
        "- soluprot_inputs/representatives_for_soluprot.fasta: manual SoluProt input; run SoluProt later.\n",
        encoding="utf-8",
    )

    print("\nDone.")
    print(f"Core FASTA: {core_fasta}")
    print(f"Tree FASTA: {outdir / 'ssn' / 'representatives.fasta'}")
    print(f"Tree: {outdir / 'tree_analysis' / 'phylogenetic_tree.newick'}")
    print(f"iTOL kingdom: {outdir / 'tree_analysis' / 'itol_kingdom_color_strip_offline.txt'}")
    print(f"SoluProt input for later: {outdir / 'soluprot_inputs' / 'representatives_for_soluprot.fasta'}")


if __name__ == "__main__":
    main()
