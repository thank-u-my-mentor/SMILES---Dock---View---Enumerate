#!/usr/bin/env bash
set -euo pipefail

src="/mnt/e/TJ/tj2_unl_NH_C4radical_M6_20260629_v1/mcpb/TJ2unlNH_M2_hessian_from_punchhess_opt12_diis_shift.inp"
dst="/mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_small_mk_M2"
stem="TJ2unlNH_M2_small_mk_from_opt12_UHF_B3LYP_fullgrid"

mkdir -p "$dst/gamess_logs" "$dst/notes"

write_header() {
  local exetyp="$1"
  {
    printf ' $SYSTEM MEMDDI=400 MWORDS=200 $END\n'
    printf ' $CONTRL\n'
    printf '   SCFTYP=UHF DFTTYP=B3LYP RUNTYP=ENERGY EXETYP=%s\n' "$exetyp"
    printf '   ICHARG=1 MULT=2 MAXIT=200 COORD=UNIQUE UNITS=ANGS\n'
    printf ' $END\n'
    printf ' $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END\n'
    printf ' $SCF\n'
    printf '   DIRSCF=.T. DIIS=.T. SOSCF=.F. DAMP=.T. SHIFT=.T.\n'
    printf '   ETHRSH=10.0 MAXDII=30\n'
    printf ' $END\n'
    printf ' $DFT\n'
    printf '   NRAD=96 NLEB=302 NRAD0=96 NLEB0=302 SWOFF=0.0 SWITCH=0.0\n'
    printf ' $END\n'
    printf ' $ELPOT IEPOT=1 WHERE=PDC $END\n'
    printf ' $PDC PTSEL=CONNOLLY CONSTR=NONE $END\n'
    printf ' $GUESS GUESS=MOREAD NORB=672 $END\n'
  }
}

write_data_vec() {
  awk '
    /^ \$DATA/ {p=1}
    p {print}
    p && /^ \$END/ {exit}
  ' "$src"
  awk '
    /^ \$VEC/ {p=1}
    p {print}
    p && /^ \$END/ {exit}
  ' "$src"
}

{
  write_header RUN
  write_data_vec
} > "$dst/${stem}.inp"

{
  write_header CHECK
  write_data_vec
} > "$dst/${stem}_check.inp"

atom_count=$(awk '
  /^ \$DATA/ {p=1; skip=2; next}
  p && /^ \$END/ {print n; exit}
  p {
    if (skip>0) {skip--; next}
    if ($1 ~ /^[A-Z][a-z]?$/ && NF >= 5) n++
  }
' "$dst/${stem}.inp")

{
  echo "Small_model M2 MK test built from successful Hessian input:"
  echo "$src"
  echo "Run input: $dst/${stem}.inp"
  echo "Check input: $dst/${stem}_check.inp"
  echo "DATA atom count: $atom_count"
  echo "UHF/B3LYP, charge +1, multiplicity 2, full DFT grid, MOREAD NORB=672."
  echo "Purpose: test whether spin-consistent small_model M2 ESP/MK can converge cleanly."
} > "$dst/notes/small_mk_summary.txt"

echo "$dst/${stem}.inp"
echo "$dst/${stem}_check.inp"
echo "atoms=$atom_count"
