#!/usr/bin/env bash
set -euo pipefail

GMXDIR=/mnt/e/Pose_create/product_R_harmonic20ns/md_handoff/R_product_mode1_hppAFe_harmonic/gromacs
PIDFILE="$GMXDIR/product_R_harmonic_20ns.pid"

if [ -f "$PIDFILE" ]; then
  pid=$(cat "$PIDFILE")
  ps -p "$pid" -o pid,etime,pcpu,pmem,cmd || true
else
  echo "pidfile missing: $PIDFILE"
fi

echo
echo "Latest MD log:"
tail -n 36 "$GMXDIR/product_R_harmonic_20ns.log" 2>/dev/null || true

echo
echo "Trajectory/checkpoint:"
ls -lh "$GMXDIR"/product_R_harmonic_20ns.{xtc,cpt,gro,edr,log} 2>/dev/null || true

echo
echo "GPU:"
nvidia-smi --query-gpu=utilization.gpu,memory.used,power.draw --format=csv,noheader 2>/dev/null || true
