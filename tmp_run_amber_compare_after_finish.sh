#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/e/TJ/260706_mcpb_m2_200ns"
AMBER_DIR="${PROJECT}/amber_200ns_pmemd_compare"
PID_FILE="${AMBER_DIR}/logs/md_200ns_pmemd_cudaSPFP.pid"
OUT_LOG="${AMBER_DIR}/md_200ns_pmemd_cudaSPFP.out"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUTDIR="${AMBER_DIR}/analysis_vs_gromacs_${STAMP}"
mkdir -p "${OUTDIR}"

echo "[watch] started at $(date)" | tee "${OUTDIR}/watch_and_compare.log"

if [[ -f "${PID_FILE}" ]]; then
  PID="$(cat "${PID_FILE}")"
  while ps -p "${PID}" >/dev/null 2>&1; do
    echo "[watch] pmemd pid ${PID} still running at $(date)" | tee -a "${OUTDIR}/watch_and_compare.log"
    sleep 300
  done
else
  echo "[watch] no pid file found; checking output directly" | tee -a "${OUTDIR}/watch_and_compare.log"
fi

echo "[watch] pmemd no longer running at $(date)" | tee -a "${OUTDIR}/watch_and_compare.log"

if ! grep -q "Total wall time" "${OUT_LOG}" && ! grep -q "FINAL PERFORMANCE" "${OUT_LOG}" && ! grep -q "TIMINGS" "${OUT_LOG}"; then
  echo "[watch] WARNING: normal final timing marker not found yet. Last 80 lines:" | tee -a "${OUTDIR}/watch_and_compare.log"
  tail -80 "${OUT_LOG}" | tee -a "${OUTDIR}/watch_and_compare.log"
  exit 1
fi

source /mnt/l/WSL/conda_envs/AmberTools25/amber.sh

cat > "${OUTDIR}/extract_amber_for_compare.cpptraj" <<EOF
parm ${AMBER_DIR}/solv.prmtop
trajin ${AMBER_DIR}/md_200ns_pmemd_cudaSPFP.nc
autoimage
rms first @CA out ${OUTDIR}/amber_CA_RMSD_A.dat
distance C4_O1 :1@C4 :1@O1 out ${OUTDIR}/amber_C4_O1_distance_A.dat
strip !(:1) parmout ${OUTDIR}/amber_UL1.prmtop
trajout ${OUTDIR}/amber_UL1.nc netcdf
run
EOF

echo "[cpptraj] extracting Amber comparison data" | tee -a "${OUTDIR}/watch_and_compare.log"
cpptraj -i "${OUTDIR}/extract_amber_for_compare.cpptraj" > "${OUTDIR}/cpptraj_extract.log" 2>&1

echo "[python] comparing Amber vs GROMACS" | tee -a "${OUTDIR}/watch_and_compare.log"
/mnt/l/WSL/softwares/conda_envs/md/bin/python /mnt/e/Codex/tmp_compare_amber_gromacs.py \
  --amber-ul1-prmtop "${OUTDIR}/amber_UL1.prmtop" \
  --amber-ul1-nc "${OUTDIR}/amber_UL1.nc" \
  --gromacs-timeseries "${PROJECT}/md_200ns/analysis_reactive_C4O1_20260708/ul1_rotamer_face_analysis/product_calibrated_pseudo_RS/md_product_calibrated_pseudo_RS_timeseries.csv" \
  --amber-ca-rmsd "${OUTDIR}/amber_CA_RMSD_A.dat" \
  --outdir "${OUTDIR}" \
  > "${OUTDIR}/python_compare.log" 2>&1

echo "[done] analysis complete at $(date)" | tee -a "${OUTDIR}/watch_and_compare.log"
