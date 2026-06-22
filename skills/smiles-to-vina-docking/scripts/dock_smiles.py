#!/usr/bin/env python
"""Dock one SMILES with Vina and append it to dock_history.csv.

This is the smallest daily-use entry point:

    python dock_smiles.py --smiles "<SMILES>" --history-dir ~/vina_task2/dock_history

On the first run, also provide --receptor, --config, and --vina. They are saved
in history_config.json, so later runs can be short.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import dock_utils as dl


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smiles", required=True, help="Input SMILES to dock")
    parser.add_argument("--nickname", default="", help="Optional user-facing label to store in the ledger")
    parser.add_argument("--history-dir", "--ledger-dir", dest="ledger_dir", metavar="HISTORY_DIR", type=Path, default=dl.DEFAULT_LEDGER_DIR, help="Directory for dock_history.csv and pose/log files; defaults to ~/vina_task2/dock_history")
    parser.add_argument("--receptor", type=Path, help="Receptor PDBQT; saved after first run")
    parser.add_argument("--config", type=Path, help="Vina config txt; saved after first run")
    parser.add_argument("--vina", help="Vina executable; saved after first run")
    parser.add_argument("--meeko", help="Meeko ligand-preparation executable; saved after first run")
    parser.add_argument("--obabel", help="OpenBabel executable; saved after first run")
    parser.add_argument("--engine-arg", dest="engine_args", action="append", default=[], help="Extra docking-engine argument; repeat for multiple tokens, e.g. --engine-arg=--cnn_scoring --engine-arg=refinement")
    parser.add_argument("--max-rounds", type=int, default=0, help="Analog generation rounds; 0 means dock only the input SMILES")
    parser.add_argument("--batch-size", type=int, default=24, help="Maximum new analogs to dock per generation round")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--cpu", type=int, help="Vina worker count; defaults to about physical-core count, capped at 16")
    parser.add_argument("--exhaustiveness", type=int)
    parser.add_argument("--num-modes", type=int, default=9)
    parser.add_argument("--energy-range", type=int)
    parser.add_argument("--internal-cluster-rmsd-cutoff", type=float, default=3.0)
    parser.add_argument("--edit-mode", default="default", help="default/add/delete/shrink/replace/drastic or comma-separated modes")
    parser.add_argument("--add", dest="add_only", action="store_true")
    parser.add_argument("--delete", dest="delete_only", action="store_true")
    parser.add_argument("--shrink", dest="shrink_only", action="store_true")
    parser.add_argument("--replace", dest="replace_only", action="store_true")
    parser.add_argument("--drastic", dest="drastic_only", action="store_true")
    parser.add_argument("--deterministic-batch", action="store_true")
    parser.add_argument("--redock-existing", action="store_true", help="Dock again even if an equivalent molecule already exists")
    parser.add_argument("--min-analog-qed", type=float, default=0.20, help="Pre-docking filter for generated analogs only")
    parser.add_argument("--max-analog-mw", type=float, default=650.0, help="Pre-docking MW filter for generated analogs only")
    parser.add_argument("--max-analog-rot-bonds", type=int, default=14, help="Pre-docking rotatable-bond filter for generated analogs only")
    parser.add_argument("--max-analog-tpsa", type=float, default=180.0, help="Pre-docking TPSA filter for generated analogs only")
    args = parser.parse_args()
    dl.collect(args)


if __name__ == "__main__":
    main()
