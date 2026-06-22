#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/e/Pose_create/product_R_harmonic20ns
HISTORY="$PROJECT/dock_history"
mkdir -p "$HISTORY"

/mnt/l/WSL/softwares/conda_envs/md/bin/python \
  /mnt/e/Codex/skills/smiles-to-vina-docking/scripts/dock_smiles.py \
  --smiles 'O=C1O[C@H](CC1)C2=CC=CC=C2' \
  --nickname product_R \
  --history-dir "$HISTORY" \
  --receptor /mnt/e/Pose_create/ple_project_AFe_noHHH_degree1/receptor/pdb2r5v_noHHH_keepFe.pdb \
  --config /mnt/e/Pose_create/ple_project_AFe_noHHH_degree1/configs/A_chain_Fe_center_21A.txt \
  --vina /home/qin/gnina/build/bin/gnina \
  --max-rounds 0 \
  --batch-size 0 \
  --exhaustiveness 8 \
  --num-modes 9 \
  --engine-arg=--cnn_scoring \
  --engine-arg=refinement
