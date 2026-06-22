#!/usr/bin/env bash
set -euo pipefail

export LD_LIBRARY_PATH=/home/qin/softwares/cp2k-2026.1/build_minlib/src:/usr/local/cuda/lib64:${LD_LIBRARY_PATH:-}
GMX=/home/qin/softwares/gromacs-2026.0-gpu/bin/gmx_mpi
PY=/mnt/l/WSL/softwares/conda_envs/md/bin/python
GMXDIR=/mnt/e/Pose_create/product_R_harmonic20ns/md_handoff/R_product_mode1_hppAFe_harmonic/gromacs
VIS=$GMXDIR/visual

mkdir -p "$VIS"
cd "$GMXDIR"

"$PY" /mnt/e/Codex/scripts/write_gmx_visual_index.py \
  --gro product_R_harmonic_20ns.gro \
  --out "$VIS/product_R_visual.ndx"

"$GMX" trjconv \
  -s product_R_harmonic_20ns.tpr \
  -f product_R_harmonic_20ns.xtc \
  -n "$VIS/product_R_visual.ndx" \
  -o "$VIS/product_R_centered_compact.xtc" \
  -pbc mol \
  -ur compact \
  -center <<'EOF'
Protein_LIG_FE
Protein_LIG_FE
EOF

"$GMX" trjconv \
  -s product_R_harmonic_20ns.tpr \
  -f product_R_harmonic_20ns.gro \
  -n "$VIS/product_R_visual.ndx" \
  -o "$VIS/product_R_centered_compact.gro" \
  -pbc mol \
  -ur compact \
  -center <<'EOF'
Protein_LIG_FE
Protein_LIG_FE
EOF

ls -lh "$VIS"
