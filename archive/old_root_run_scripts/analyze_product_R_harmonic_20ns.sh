#!/usr/bin/env bash
set -euo pipefail

ROOT=/mnt/e/Pose_create/product_R_harmonic20ns/md_handoff/R_product_mode1_hppAFe_harmonic
GMXDIR="$ROOT/gromacs"
OUT=/mnt/e/Pose_create/product_R_harmonic20ns/md_analysis/product_R_harmonic_20ns_rdc
VEL=/mnt/e/Pose_create/product_R_harmonic20ns/md_analysis/product_R_harmonic_20ns_rdc_velocity
PY=/mnt/l/WSL/softwares/conda_envs/md/bin/python

if [ ! -s "$GMXDIR/product_R_harmonic_20ns.gro" ]; then
  echo "Production .gro not found yet. Wait for mdrun to finish." >&2
  exit 2
fi
if [ ! -s "$GMXDIR/product_R_harmonic_20ns.xtc" ]; then
  echo "Production .xtc not found yet." >&2
  exit 2
fi

"$PY" /mnt/e/Codex/scripts/ple_md_rdc.py \
  --topology "$GMXDIR/product_R_harmonic_20ns.gro" \
  --trajectory "$GMXDIR/product_R_harmonic_20ns.xtc" \
  --outdir "$OUT" \
  --ligand-selection "resname LIG" \
  --distance-cutoff 10 \
  --stride 1 \
  --top-n 30 \
  --plot-traces

"$PY" /mnt/e/Codex/scripts/plot_rdc_velocity_distribution.py \
  --timeseries "$OUT/residue_ligand_distance_timeseries.csv" \
  --rdc-csv "$OUT/residue_ligand_rdc_within_cutoff.csv" \
  --outdir "$VEL" \
  --top-n 6 \
  --max-x 20 \
  --bins 44 \
  --title "Distribution density of distance fluctuations between residues and R-product"

"$PY" /mnt/e/Codex/scripts/build_velocity_distribution_report.py \
  --rdc-dir "$OUT" \
  --velocity-dir "$VEL" \
  --title "Distribution Density of Distance Fluctuations Between Residues and R-Product" \
  --subtitle "R-product harmonic Fe-core positive-control MD, 20 ns, 10 ps frame spacing; ligand selection: <code>resname LIG</code>" \
  --top-n 12

"$PY" /mnt/e/Codex/scripts/fe_core_distance_qc.py \
  --structure "$GMXDIR/product_R_harmonic_20ns.gro" \
  --trajectory "$GMXDIR/product_R_harmonic_20ns.xtc" \
  --outdir "$OUT" \
  --pair Fe_HID160_NE2,5142,2302 \
  --pair Fe_HID240_NE2,5142,3486 \
  --pair Fe_GLU319_OE1,5142,4719

echo "report=$OUT/index.html"
