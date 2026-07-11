#!/usr/bin/env bash
set -eo pipefail

base=/mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_chgmod_tests_20260707_v2
src=/mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_cutoff23
origlog=/mnt/e/TJ/tj2_unl_NH_C4radical_M6_20260629_v1/mcpb/gamess_logs/TJ2unlNH_large_mk_uhf_direct_fullgrid_8core.log

mkdir -p "$base"

for mode in 3b 3c 3d; do
  d="$base/$mode"
  mkdir -p "$d"
  cp -a "$src"/. "$d"/
  cp "$origlog" "$d/large_mk_M6.log"
  cd "$d"
  set +u
  source /mnt/l/WSL/conda_envs/AmberTools25/amber.sh >/dev/null 2>&1
  set -u
  if ! MCPB.py -i mcpb.in -s "$mode" --logf large_mk_M6.log > "mcpb_${mode}_M6log_constraint_test.log" 2>&1; then
    echo "MODE_${mode}_FAILED"
    tail -80 "mcpb_${mode}_M6log_constraint_test.log"
    exit 1
  fi
  echo "MODE_${mode}_OK"
  grep -aE "Check the large model|Check the standard model|Good\\. There are|Generating the .*mol2|RESP Charge fitting" \
    "mcpb_${mode}_M6log_constraint_test.log" | tail -40 || true
done

echo "DONE $base"
