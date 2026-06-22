#!/usr/bin/env bash
set -euo pipefail

ROOT=/mnt/e/Pose_create/product_R_harmonic20ns/md_handoff/R_product_mode1_hppAFe_harmonic
GMXDIR="$ROOT/gromacs"
cd "$GMXDIR"

/mnt/l/WSL/softwares/conda_envs/md/bin/python /mnt/e/Codex/scripts/make_fe_core_harmonic_top.py \
  --input-top complex_GMX.top \
  --output-top complex_GMX_fe3core_harmonic.top \
  --bond 5142,2302,0.2106,100000,HID160_NE2_FE_positive_control \
  --bond 5142,3486,0.2055,100000,HID240_NE2_FE_positive_control \
  --bond 5142,4719,0.2084,100000,GLU319_OE1_FE_positive_control

cp md_20ns.mdp product_R_fe3core_harmonic_20ns_10ps.mdp

GMX_BIN="${GMX_BIN:-}"
if [ -z "$GMX_BIN" ]; then
  for candidate in \
    /usr/local/gromacs/bin/gmx_mpi \
    /usr/local/gromacs/bin/gmx \
    /usr/local/bin/gmx_mpi \
    /usr/local/bin/gmx \
    /mnt/l/WSL/softwares/gromacs/bin/gmx_mpi \
    /mnt/l/WSL/softwares/gromacs/bin/gmx \
    /mnt/l/WSL/softwares/gromacs_cp2k/bin/gmx_mpi \
    /mnt/l/WSL/softwares/gromacs_cp2k/bin/gmx \
    /mnt/l/WSL/softwares/conda_envs/md/bin/gmx_mpi \
    /mnt/l/WSL/softwares/conda_envs/md/bin/gmx \
    gmx_mpi \
    gmx
  do
    if command -v "$candidate" >/dev/null 2>&1; then
      GMX_BIN="$(command -v "$candidate")"
      break
    fi
    if [ -x "$candidate" ]; then
      GMX_BIN="$candidate"
      break
    fi
  done
fi

if [ -z "$GMX_BIN" ]; then
  echo "Could not find gmx/gmx_mpi. Source GMXRC or set GMX_BIN=/path/to/gmx before running." >&2
  exit 2
fi

echo "gmx_bin=$GMX_BIN"
"$GMX_BIN" --version | sed -n '1,12p'

MDRUN_ARGS=(-nb gpu -pme gpu)

"$GMX_BIN" grompp -f minim.mdp -c complex_GMX.gro -p complex_GMX_fe3core_harmonic.top -o product_R_harmonic_minim.tpr -maxwarn 2 -v
"$GMX_BIN" mdrun -deffnm product_R_harmonic_minim -v

"$GMX_BIN" grompp -f nvt.mdp -c product_R_harmonic_minim.gro -r product_R_harmonic_minim.gro -p complex_GMX_fe3core_harmonic.top -o product_R_harmonic_nvt.tpr -maxwarn 2
"$GMX_BIN" mdrun -deffnm product_R_harmonic_nvt -v "${MDRUN_ARGS[@]}"

"$GMX_BIN" grompp -f npt.mdp -c product_R_harmonic_nvt.gro -r product_R_harmonic_nvt.gro -t product_R_harmonic_nvt.cpt -p complex_GMX_fe3core_harmonic.top -o product_R_harmonic_npt.tpr -maxwarn 2
"$GMX_BIN" mdrun -deffnm product_R_harmonic_npt -v "${MDRUN_ARGS[@]}"

"$GMX_BIN" grompp -f product_R_fe3core_harmonic_20ns_10ps.mdp -c product_R_harmonic_npt.gro -t product_R_harmonic_npt.cpt -p complex_GMX_fe3core_harmonic.top -o product_R_harmonic_20ns.tpr -maxwarn 2
"$GMX_BIN" mdrun -deffnm product_R_harmonic_20ns -v "${MDRUN_ARGS[@]}"
