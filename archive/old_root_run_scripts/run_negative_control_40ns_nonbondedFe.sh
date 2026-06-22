#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/mnt/e/Pose_create/ple_project_AFe_noHHH_degree1/md_handoff/S000003_mode8_hppAFe_terfix/gromacs"
GMX="/home/qin/softwares/gromacs-2026.0-gpu/bin/gmx_mpi"
CP2K_LIB="/home/qin/softwares/cp2k-2026.1/build_minlib/src"

cd "$WORKDIR"
export LD_LIBRARY_PATH="${CP2K_LIB}:${LD_LIBRARY_PATH:-}"

exec "$GMX" mdrun \
  -deffnm negative_control_40ns_nonbondedFe \
  -s negative_control_40ns_nonbondedFe.tpr \
  -ntomp 16 \
  -pin on
