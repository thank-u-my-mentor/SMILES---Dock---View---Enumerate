#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/mnt/e/Pose_create/ple_project_AFe_noHHH_degree1/md_handoff/S000003_mode8_hppAFe_terfix/gromacs"
RUN="negative_control_40ns_nonbondedFe"

cd "$WORKDIR"
echo "== process =="
if [[ -f "${RUN}.pid" ]]; then
  ps -p "$(cat "${RUN}.pid")" -o pid,ppid,stat,etime,cmd || true
else
  echo "no ${RUN}.pid"
fi

echo
echo "== gpu =="
if [[ -x /usr/lib/wsl/lib/nvidia-smi ]]; then
  /usr/lib/wsl/lib/nvidia-smi
else
  nvidia-smi || true
fi

echo
echo "== last progress =="
if [[ -f "${RUN}.log" ]]; then
  grep -E "Step|Performance|Time:|ns/day|Writing final coordinates|Finished mdrun" "${RUN}.log" | tail -n 30
  echo
  tail -n 25 "${RUN}.log"
else
  echo "no ${RUN}.log yet"
fi

echo
echo "== outputs =="
ls -lh "${RUN}"* 2>/dev/null || true
