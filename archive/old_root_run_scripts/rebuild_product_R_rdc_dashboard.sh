#!/usr/bin/env bash
set -euo pipefail

GMXDIR=/mnt/e/Pose_create/product_R_harmonic20ns/md_handoff/R_product_mode1_hppAFe_harmonic/gromacs
OUT=/mnt/e/Pose_create/product_R_harmonic20ns/md_analysis/product_R_harmonic_20ns_rdc
REFPDB=/mnt/e/Pose_create/product_R_harmonic20ns/md_handoff/R_product_mode1_hppAFe_harmonic/input/receptor_md.pdb
PY=/mnt/l/WSL/softwares/conda_envs/md/bin/python

"$PY" /mnt/e/Codex/scripts/ple_md_rdc.py \
  --topology "$GMXDIR/product_R_harmonic_20ns.gro" \
  --trajectory "$GMXDIR/product_R_harmonic_20ns.xtc" \
  --outdir "$OUT" \
  --ligand-selection "resname LIG" \
  --distance-cutoff 10 \
  --stride 1 \
  --top-n 30 \
  --plot-traces

"$PY" /mnt/e/Codex/scripts/fe_core_distance_qc.py \
  --structure "$GMXDIR/product_R_harmonic_20ns.gro" \
  --trajectory "$GMXDIR/product_R_harmonic_20ns.xtc" \
  --outdir "$OUT" \
  --pair Fe_HID160_NE2,5142,2302 \
  --pair Fe_HID240_NE2,5142,3486 \
  --pair Fe_GLU319_OE1,5142,4719

"$PY" /mnt/e/Codex/scripts/map_rdc_labels_to_reference_pdb.py \
  --rdc-dir "$OUT" \
  --reference-pdb "$REFPDB" \
  --gro "$GMXDIR/product_R_harmonic_20ns.gro"

"$PY" /mnt/e/Codex/scripts/build_rdc_dashboard.py \
  --rdc-dir "$OUT" \
  --title "氨基酸-底物 RDC 动态分析 Dashboard" \
  --subtitle "R 构型产物的 harmonic Fe-core control MD：20 ns，10 ps 帧间隔。该页面用于探索性动态分析，不作为催化效率的硬分类器。" \
  --velocity-threshold 3.0

echo "dashboard=$OUT/index.html"
