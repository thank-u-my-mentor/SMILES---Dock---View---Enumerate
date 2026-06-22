#!/usr/bin/env bash
set -euo pipefail

GMXDIR=/mnt/e/Pose_create/product_R_harmonic20ns/md_handoff/R_product_mode1_hppAFe_harmonic/gromacs
cd "$GMXDIR"

export LD_LIBRARY_PATH=/home/qin/softwares/cp2k-2026.1/build_minlib/src:/usr/local/cuda/lib64:${LD_LIBRARY_PATH:-}
GMX_BIN=/home/qin/softwares/gromacs-2026.0-gpu/bin/gmx_mpi

for ext in log xtc edr cpt; do
  if [ ! -e "product_R_harmonic_20ns.$ext" ] && [ -e "product_R_harmonic_20ns_failed_no_bond_insertion.$ext" ]; then
    mv "product_R_harmonic_20ns_failed_no_bond_insertion.$ext" "product_R_harmonic_20ns.$ext"
  fi
done

if pgrep -f "product_R_harmonic_20ns.*mdrun" >/dev/null 2>&1; then
  echo "product_R_harmonic_20ns mdrun already appears to be running" >&2
  pgrep -af "product_R_harmonic_20ns.*mdrun" >&2
  exit 3
fi

nohup "$GMX_BIN" mdrun \
  -deffnm product_R_harmonic_20ns \
  -s product_R_harmonic_20ns.tpr \
  -cpi product_R_harmonic_20ns.cpt \
  -append \
  -nb gpu \
  -pme gpu \
  -cpt 5 \
  -v \
  > product_R_harmonic_20ns.stdout.log 2>&1 &

echo $! > product_R_harmonic_20ns.pid
echo "pid=$(cat product_R_harmonic_20ns.pid)"
echo "log=$GMXDIR/product_R_harmonic_20ns.log"
