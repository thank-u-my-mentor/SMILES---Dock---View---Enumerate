#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/e/Pose_create/product_R_harmonic20ns
STATE=ple_project
OUTDIR="$PROJECT/md_handoff/R_product_mode1_hppAFe_harmonic"

mkdir -p "$PROJECT/$STATE"
if [ ! -d "$PROJECT/$STATE/dock_history" ]; then
  cp -a "$PROJECT/dock_history" "$PROJECT/$STATE/dock_history"
fi

/mnt/l/WSL/softwares/conda_envs/md/bin/python -m pocket_ligand_explorer.cli md-handoff \
  --project-dir "$PROJECT" \
  --state-dir "$STATE" \
  --pdb /mnt/e/Pose_create/pdb2r5v_chainA_GLU108_Hpp_heavy_plus_AFe_noHHH.pdb \
  --seq-id S000001 \
  --pose-mode 1 \
  --outdir "$OUTDIR" \
  --ligand-charge 0 \
  --receptor-md-mode protein-only \
  --protonation-map /mnt/e/Pose_create/pdb2r5v_Hpp_pH7p4_first_pass_amber_protonation_map.csv \
  --keep-hetatm-resname FE \
  --keep-hetatm-element FE \
  --ambertools-bin /mnt/l/WSL/conda_envs/AmberTools25/bin \
  --gmx gmx_mpi \
  --time 20 \
  --dt-ps 0.002 \
  --maxwarn 2 \
  --mdrun-args="-nb gpu -pme gpu" \
  --extra-ld-library-path /usr/local/cuda/lib64
