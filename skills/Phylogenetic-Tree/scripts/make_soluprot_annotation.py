#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


HYDROLASE_SCRIPT = Path("/mnt/e/Codex/skills/hydrolase-homolog-tree/scripts/annotate_solubility_itol.py")


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert manual SoluProt/NetSolP CSV to iTOL annotation files.")
    parser.add_argument("--input", default="/mnt/e/Tree-Metal-F/soluprot_inputs/input.csv")
    parser.add_argument("--outdir", default="/mnt/e/Tree-Metal-F/tree_analysis")
    parser.add_argument("--nodes", default="")
    parser.add_argument("--id-mapping", default="/mnt/e/Tree-Metal-F/soluprot_inputs/soluprot_id_mapping.tsv")
    parser.add_argument("--tool-label", default="SoluProt")
    parser.add_argument("--low-threshold", default="0.2")
    args = parser.parse_args()

    cmd = [
        sys.executable,
        str(HYDROLASE_SCRIPT),
        args.input,
        "--outdir",
        args.outdir,
        "--tool-label",
        args.tool_label,
        "--low-threshold",
        str(args.low_threshold),
    ]
    if args.nodes:
        cmd.extend(["--nodes", args.nodes])
    if args.id_mapping:
        cmd.extend(["--id-mapping", args.id_mapping])
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
