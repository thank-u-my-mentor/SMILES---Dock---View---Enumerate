# MCPB / Fe-Radical Bug Avoidance Guide

This note is intentionally AI-readable and project-practical. Use it before
running an Fe-containing MCPB.py workflow from PLE, especially for non-heme iron
enzymes, radical ligands, H++-processed PDB files, and hand-built catalytic
intermediates.

## 1. Electron Count And Multiplicity

Always check the actual model that GAMESS will run, not the earlier hand-built
PDB or spin-screen model.

Formula:

```text
total_electrons = sum(atomic_numbers) - total_charge
```

Parity rule:

```text
even electron count -> odd multiplicity is allowed: 1, 3, 5, 7, ...
odd electron count  -> even multiplicity is allowed: 2, 4, 6, ...
```

Why:

```text
MULT = 2S + 1
```

For `MULT=7`, `S=3`, which corresponds to an even-electron high-spin state in a
simple spin-counting picture. Therefore a GAMESS job with `MULT=7` must have an
even total electron count. If the model has 289 electrons, `MULT=7` is
inconsistent; changing the charge to `+1` gives 288 electrons and makes `MULT=7`
electron-parity consistent.

For the TJ MCPB-generated models:

```text
small model: 62 atoms, charge=+1, MULT=7, 288 electrons
large model: 92 atoms, charge=+1, MULT=7, 456 electrons
```

Important trap: a custom spin-screen model and the MCPB-generated small/large
models can have different atom counts, caps, hydrogens, waters, and total
charges. Do not blindly reuse the spin-screen charge in MCPB. Recompute electron
counts after `MCPB.py -s 1` has generated `*_small.pdb` and `*_large.pdb`.

## 2. Radical Ligand H-Deletion

If a ligand is intended to be a carbon radical, the hydrogen removed from the
radical carbon must stay removed throughout the workflow.

For the TJ intermediate, the intended radical correction was:

```text
UNL A1 C4 loses one H, e.g. H4A
UNL remains odd-electron before the full metal-site charge is assigned
```

This radical is exactly why charge/multiplicity errors are easy to trigger. The
UNL fragment changes electron parity, but the final MCPB model parity also
depends on Fe oxidation state, ACT charge, Glu/His capping, the coordinating
water, and MCPB's own truncation/capping atoms.

Do not let RDKit/OpenBabel silently regenerate the missing radical H when
preparing the final MCPB ligand coordinates.

## 3. Chain Selection

For homodimeric enzymes, use a chain-specific MCPB model unless the metal site
really spans chains.

TJ lesson:

```text
Use chain A only for the Fe site.
Do not feed both chain A and chain B to MCPB if they share residue numbers.
```

Reason: `pymsmt`/MCPB may index residues internally by residue number and can
merge or confuse chain A/B residues with identical residue IDs. This corrupts
side-chain truncation and ligand mapping.

## 4. H++ And Protonation

H++ is used to decide protein protonation, but the MCPB PDB must still be Amber
and MCPB friendly.

Checklist:

```text
1. Remove irrelevant chains/ligands before H++ if needed.
2. Keep the catalytic Fe and required ligands only if H++ can process them.
3. Preserve TER records between chains.
4. Convert H++ histidine decisions to Amber names: HID/HIE/HIP.
5. For Fe-coordinating histidines donating through NE2, use HID unless the
   local chemistry clearly says otherwise.
6. Keep a protonation_map.csv and the original pkout as provenance.
```

Do not assume GLU/ASP are protonated only because they are near Fe. Use H++,
visual inspection, and the catalytic mechanism. For Fe ligands, the residue
name and protonation state must be chemically intentional.

## 5. Six-Coordinate Fe Model

For non-heme Fe catalytic intermediates, decide the coordination sphere before
MCPB. Do not let MD discover Fe coordination from a nonbonded ion model.

TJ M7 model:

```text
Fe(III), high-spin M7 model
HIS187 NE2
HIS270 NE2
GLU349 OE1
ACT501 O2
UNL N1
HOH875 O
```

If the literature intermediate contains a ligand, radical, azide, nitrene,
water, dioxygen, or superoxo in the first coordination sphere, include it in the
MCPB model. The Fe bonded model is only as chemically meaningful as the
coordination sphere supplied to MCPB.

## 6. Names And Files MCPB Expects

Keep names consistent.

```text
PDB Fe residue name: FE
Fe mol2 residue name: FE
coordinating water: HOH with O/H1/H2 atoms if treated as ligand
non-standard ligands: ACT, UNL, HOH in naa_mol2files
ligand frcmod files: ACT.frcmod, UNL.frcmod
```

Example:

```text
naa_mol2files ACT.mol2 UNL.mol2 HOH.mol2
frcmod_files ACT.frcmod UNL.frcmod
```

If a residue is heavy-only and MCPB side-chain truncation needs a cap atom such
as `HA`, add the minimum helper atoms and document that they are modeling helper
atoms, not experimental coordinates.

## 7. GAMESS / MCPB Runtime Checks

A GAMESS job is usable only if the log ends normally:

```text
EXECUTION OF GAMESS TERMINATED NORMALLY
```

Common benign warnings:

```text
IEEE_UNDERFLOW_FLAG
IEEE_DENORMAL
```

These are usually floating-point noise, not a failed calculation.

A geometry optimization that stops because it hit `NSTEP` is not automatically
usable. Prefer a log that prints:

```text
......END OF GEOMETRY SEARCH......
```

and has small gradients. If not converged, continue from the last reliable
computed coordinate block, then export the converged structure.

MCPB-generated `*_small_fc.inp` can have an empty `$DATA` block. Correct flow:

```text
1. Run *_small_opt.inp.
2. Extract optimized coordinates.
3. Fill *_small_fc.inp with those coordinates.
4. Run Hessian / force constants.
5. Run large-model MK/RESP job.
6. Run MCPB.py steps 2, 3, 4.
7. Only then build Amber/GROMACS MD topology.
```

For GAMESS parallelism, more cores are not always faster. Start with 8-16 cores
on a strong desktop CPU; if DDI errors or poor scaling appear, reduce cores.

When launching GAMESS from tmux, keep the tmux command simple. Do not inline a
long quoted `rungms` pipeline inside `tmux new-session`; nested quoting can split
`ddikick.x`, the job name, or shell variables and create misleading errors such
as:

```text
ddikick: command not found
ZFY_M7_FE_small_opt: command not found
```

Preferred pattern:

```text
run_03_small_opt_foreground.sh  # contains the real rungms command
run_03_small_opt_tmux.sh        # only runs: bash run_03_small_opt_foreground.sh
```

Verify a real job with `pgrep -af ZFY_M7_FE_small_opt` or
`pgrep -af ddikick`. A healthy run shows one `ddikick.x` parent and one
`gamess.00.x` worker per requested core.

For SCF failures in a Fe-radical small optimization, do not immediately proceed
to Hessian. A practical M7-only retry pattern is:

```text
NCORES=8
$CONTRL ... MAXIT=200
$SCF DIRSCF=.T. DIIS=.T. DAMP=.T. SHIFT=.T. ETHRSH=10.0 MAXDII=30 $END
$STATPT NSTEP=300 OPTTOL=0.0002 $END
```

GAMESS input-card pitfalls observed in the ZFY M7 retry:

```text
MAXIT must be <= 200.
Do not put very long $SCF cards on one line; old GAMESS card parsing can lose $END.
Do not enable DIIS and SOSCF together; GAMESS reports them as mutually exclusive.
```

The helper `scripts/prepare_gamess_stable_retry.py` writes the stable retry
input and tmux runner. For ZFY, the M7-only stable retry job name is:

```text
ZFY_M7_FE_small_opt_stable
```

If B3LYP/6-31G(d) M7 fails in the first SCF, try a deliberately cheap
pre-optimization ladder before changing the chemistry:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/prepare_gamess_preopt_retry.py \
  --mcpb-dir /mnt/e/ZFY/mcpb_m7_small_opt/mcpb \
  --source-job ZFY_M7_FE_small_opt \
  --preopt-job ZFY_M7_FE_small_preopt_rohf_sto3g_soscf \
  --charge 1 \
  --mult 7 \
  --method rohf \
  --basis sto3g \
  --runtype optimize \
  --scf-mode soscf \
  --cores 8
```

Do not mistake a shell `exit_code=0` or `TERMINATED NORMALLY` for a usable QM
state. GAMESS can terminate normally after an unconverged single-point run. The
log must not contain either of these:

```text
SCF IS UNCONVERGED
SCF HAS NOT CONVERGED
```

The ZFY 2R5V M7 Fe-radical small model showed this pattern on 2026-06-22:

```text
UHF/STO-3G optimize:      SCF unconverged, <S^2> = 17.082
SVWN/STO-3G optimize:     SCF unconverged, <S^2> = 12.002
ROHF/STO-3G optimize:     SCF unconverged, <S^2> = 12.000
ROHF/STO-3G/SOSCF energy: normal termination but SCF unconverged, <S^2> = 12.000
ROHF/STO-3G/SOSCF opt:    SCF unconverged before any geometry step
```

Interpretation: the requested `MULT=7` spin expectation is plausible, but the
current electronic model is not self-consistent enough for Hessian or MCPB force
constants. Stop blind parameter hopping and audit the chemistry:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/audit_mcpb_project.py \
  --mcpb-dir /mnt/e/ZFY/mcpb_m7_small_opt/mcpb \
  --group ZFY_M7_FE \
  --cutoff 3.0
```

For the current ZFY model, the audit shows the six Fe donors are geometrically
reasonable but chemically demanding:

```text
UNK O2  1.979 A
HID241 NE2  2.060 A
GLU320 OE1  2.086 A
HID161 NE2  2.101 A
ACT O1  2.142 A
UNK N1  2.234 A
```

Before retrying long GAMESS jobs, verify whether `UNK` should really bind Fe
through both `N1` and `O2`, whether ACT should bind via `O1` or `O2`, and whether
the total model charge `+1` is still correct for Fe(III) + ACT(-1) + UNK radical
plus the chosen GLU/HID truncation.

The revised ZFY `2.pdb` scaffold added a coordinating water and moved the UNK
geometry so the first-shell model became chemically cleaner:

```text
Project: /mnt/e/ZFY/mcpb_m7_2pdb_stable
Input:   /mnt/e/ZFY/zfy_mcpb_scaffold/2.pdb
Fe:      FE A4113
Charge:  +1
MULT:    7

UNK A501 N1      2.047 A
HID A241 NE2     2.060 A
HOH A500 O       2.063 A
ACT A502 O1      2.071 A
GLU A320 OE1     2.086 A
HID A161 NE2     2.101 A
```

MCPB.py step 1 passed for this revised scaffold:

```text
small model: 56 atoms, 284 electrons
large model: 86 atoms, 452 electrons
```

When a coordinating `HOH` appears within the Fe cutoff, include `HOH.mol2` in
`naa_mol2files` but do not add a water `frcmod` to `frcmod_files`:

```text
naa_mol2files UNK.mol2 ACT.mol2 HOH.mol2
frcmod_files UNK.frcmod ACT.frcmod
```

The helper `prepare_fe_radical_mcpb_opt_project.py` now writes this HOH mol2
automatically. Without this, MCPB.py step 1 fails with:

```text
HOH is required in naa_mol2files but not provided.
```

Fast single-point smoke tests on the revised `2.pdb` model still did not prove a
usable QM state:

```text
ROHF/STO-3G/SOSCF energy: normal termination but SCF unconverged, energy=0
UHF/STO-3G/DIIS energy:   normal termination but SCF unconverged, <S^2>=19.825
SVWN/STO-3G/DIIS energy:  strong SCF oscillation; stop instead of burning CPU
```

Therefore the active path is a stable B3LYP/6-31G(d) M7 small optimization:

```bash
cd /mnt/e/ZFY/mcpb_m7_2pdb_stable/mcpb
env NCORES=8 SESSION=zfy2_m7_small_opt_stable \
  bash run_ZFY2_M7_FE_small_opt_stable_tmux.sh
bash monitor_ZFY2_M7_FE_small_opt_stable.sh
tail -f gamess_logs/ZFY2_M7_FE_small_opt_stable.log
```

At launch on 2026-06-22, this run reached `U-B3LYP SCF CALCULATION` and the DIIS
error dropped to about `8.5e-3` by SCF iteration 16, but it ultimately failed at
the first geometry point:

```text
SCF IS UNCONVERGED, TOO MANY ITERATIONS
S-SQUARED = 12.038
FAILURE TO LOCATE STATIONARY POINT, SCF HAS NOT CONVERGED
EXECUTION OF GAMESS TERMINATED -ABNORMALLY-
```

Do not proceed to Hessian from this log. A cheaper `B3LYP/3-21G` preoptimization
also looked promising through about SCF iteration 50, but after DFT switched on
it diverged badly with huge energy/density jumps. Stop that run rather than
letting it burn CPU.

The current fallback is an even more conservative `UHF/3-21G` M7 preoptimization
without DFT switching:

```bash
cd /mnt/e/ZFY/mcpb_m7_2pdb_stable/mcpb
env NCORES=8 SESSION=zfy2_m7_preopt_uhf321g \
  bash run_ZFY2_M7_FE_small_preopt_uhf_321g_tmux.sh
bash monitor_ZFY2_M7_FE_small_preopt_uhf_321g.sh
tail -f gamess_logs/ZFY2_M7_FE_small_preopt_uhf_321g.log
```

Only continue toward Hessian after a job reaches a real geometry-search result
with no `SCF IS UNCONVERGED` and no `SCF HAS NOT CONVERGED`.

The `UHF/3-21G` M7 fallback also failed for the revised ZFY `2.pdb` model:

```text
SCF IS UNCONVERGED, TOO MANY ITERATIONS
S-SQUARED = 21.885
FAILURE TO LOCATE STATIONARY POINT, SCF HAS NOT CONVERGED
EXECUTION OF GAMESS TERMINATED -ABNORMALLY-
```

At this point do not keep forcing `MULT=7`. For Fe(III) high spin plus an
organic radical, `MULT=5` is a chemically plausible antiferromagnetically
coupled state and should be tested as a separate clean project while preserving
the same donor set and total charge:

```text
Project: /mnt/e/ZFY/mcpb_m5_2pdb_stable
Input:   /mnt/e/ZFY/zfy_mcpb_scaffold/2.pdb
Charge:  +1
MULT:    5

small model: 56 atoms, 284 electrons
large model: 86 atoms, 452 electrons

UNK A501 N1      2.047 A
HID A241 NE2     2.060 A
HOH A500 O       2.063 A
ACT A502 O1      2.071 A
GLU A320 OE1     2.086 A
HID A161 NE2     2.101 A
```

The first M5 stable B3LYP/6-31G(d) small optimization was started as:

```bash
cd /mnt/e/ZFY/mcpb_m5_2pdb_stable/mcpb
env NCORES=8 SESSION=zfy2_m5_small_opt_stable \
  bash run_ZFY2_M5_FE_small_opt_stable_tmux.sh
tail -f gamess_logs/ZFY2_M5_FE_small_opt_stable.log
```

Early SCF behavior looked more promising than the failed M7 retries: by
iteration 15, DIIS error dropped from about `0.69` to about `0.0077`, and the
alpha/beta occupations were `144/140`, matching `MULT=5`. This is not yet a
success; it only means the M5 route is now the correct active test. Proceed to
Hessian only after the M5 log reaches a real optimized geometry with no SCF
failure markers.

The first full-level M5 `B3LYP/6-31G(d)` small optimization later diverged after
DFT was switched on. Around SCF iteration 163-181 the energy jumped from roughly
`-3024` to positive values and back, density changes reached hundreds, and DIIS
was repeatedly disabled by GAMESS. This run was stopped manually and saved as:

```text
/mnt/e/ZFY/mcpb_m5_2pdb_stable/mcpb/gamess_logs/
  ZFY2_M5_FE_small_opt_stable_diverged_manualstop.log
```

Do not use this log for Hessian or MCPB force constants.

For the revised ZFY `5_for_tleap.pdb` M5 model, the useful rescue pattern was:

```text
1. Keep the chemistry fixed: charge=+1, MULT=5, same six Fe donors.
2. Do not continue directly to geometry optimization.
3. First run fixed-geometry ENERGY with DIIS disabled and SOSCF enabled.
4. If the first SOSCF ENERGY reaches MAXIT but writes a large .dat file, use
   its $VEC block with GUESS=MOREAD for a second fixed-geometry ENERGY run.
5. Only after the MOREAD ENERGY says DENSITY CONVERGED should optimization be
   retried from the converged vector.
```

Observed ZFY M5 results on 2026-06-22:

```text
DIIS ENERGY:
  reached near -3024.182 at iteration 42
  GAMESS switched DIIS off after a tiny energy rise
  SCF then exploded; stop this route

SOSCF ENERGY from Huckel:
  stable, did not explode
  reached MAXIT=200 but still printed SCF IS UNCONVERGED
  wrote /home/qin/softwares/gamess/restart/ZFY5_M5_FE_scf_energy_soscf.dat

SOSCF ENERGY with GUESS=MOREAD from that .dat:
  DENSITY CONVERGED
  FINAL U-B3LYP ENERGY = -3024.1031719382 after 67 iterations
  S-SQUARED = 12.156
```

Key lesson: in this Fe-radical case, "GAMESS terminated normally" was not enough
for the first SOSCF run because it still contained `SCF IS UNCONVERGED`. The
usable state was the second `GUESS=MOREAD` run, which explicitly printed
`DENSITY CONVERGED`. Use that converged vector as the starting point for the next
geometry optimization attempt.

Older low-level fallback tests were less useful than the converged B3LYP MOREAD
vector:

```text
M5/UHF/STO-3G energy-only:
  normal shell/GAMESS termination, but SCF was still unconverged
  FINAL UHF ENERGY = 0.0000000000
  <S^2> = 16.932

M5/ROHF/STO-3G/SOSCF energy-only:
  worse behavior; energy rose and orbital gradient increased
  manually stopped
```

Interpretation: after `ZFY5_M5_FE_scf_energy_soscf_moread1` prints
`DENSITY CONVERGED`, do not go back to the failed STO-3G ladder. The active path
is a short conservative `OPTIMIZE` from the converged B3LYP vector:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/prepare_gamess_moread_optimize.py \
  --mcpb-dir /mnt/e/ZFY/mcpb_5for_tleap_M5_try/mcpb \
  --source-job ZFY5_M5_FE_small_opt \
  --out-job ZFY5_M5_FE_small_opt_moread1_nstep60 \
  --dat /home/qin/softwares/gamess/restart/ZFY5_M5_FE_scf_energy_soscf_moread1.dat \
  --charge 1 \
  --mult 5 \
  --norb 630 \
  --nstep 60 \
  --opttol 0.0002 \
  --cores 8 \
  --session zfy5_m5_opt_moread1_n60

cd /mnt/e/ZFY/mcpb_5for_tleap_M5_try/mcpb
env NCORES=8 SESSION=zfy5_m5_opt_moread1_n60 \
  bash run_ZFY5_M5_FE_small_opt_moread1_nstep60_tmux.sh
bash monitor_ZFY5_M5_FE_small_opt_moread1_nstep60.sh
tail -f gamess_logs/ZFY5_M5_FE_small_opt_moread1_nstep60.log
```

This retry intentionally uses:

```text
RUNTYP=OPTIMIZE
GUESS=MOREAD
DIIS=.F.
SOSCF=.T.
DAMP=.T.
SHIFT=.T.
NSTEP=60
```

If it reaches `END OF GEOMETRY SEARCH`, export the optimized coordinates and then
prepare the Hessian input. If it reaches `NSTEP` without full convergence but
prints evaluated `NSERCH` coordinate blocks with stable Fe-donor distances, use a
continuation from the last evaluated geometry plus the newest `.dat` vector. Do
not use a predicted next geometry whose energy/gradient was never evaluated.
If this MOREAD optimization again explodes or prints `SCF HAS NOT CONVERGED`,
then stop and re-audit the chemical model: UNK donor choice, ACT donor choice,
water orientation, total charge, and radical hydrogen deletion.

For PDB-vs-mol2 ligand connectivity confusion, audit the topology source rather
than trusting PyMOL sticks. In the ZFY `5_for_tleap.pdb` model, the PDB `CONECT`
records contained a bogus `UNK C1-C3` bond, and PyMOL displayed it. The actual
`UNK.mol2` template did not contain `C1-C3`; it had `C1-C2`, `C1-O1`, `C1-H02`,
`C1-H03`, and `C3-C2/C4/C8`. For Amber/tleap/MCPB ligand templates, a bad PDB
`CONECT` is mostly a visualization problem if the mol2 template and atom names
are correct. A bad mol2 bond is serious and must be fixed before MD.

## 8. Do Not Start Production MD Until These Pass

Before any 50 ns, 100 ns, or 400 ns production run:

```text
1. MCPB step 4 produced final frcmod/mol2/lib files without fatal errors.
2. tleap builds solvated topology without missing parameters.
3. The Fe-donor bonded pairs exist in the Amber topology.
4. Energy minimization does not explode.
5. NVT and NPT finish.
6. PBC-fixed trajectory is generated for PyMOL/VMD analysis.
7. Fe-donor distances remain chemically reasonable in a short test.
```

For visual analysis, always keep both:

```text
raw xtc/trr
pbc-fixed xtc/trr
```

Never judge PyMOL trajectory shape from an unfixed PBC trajectory.

## 9. TJ Directory Cleanup Rule

The clean TJ template should be:

```text
E:/TJ/
  1T47.pdb
  1T47_chainA_Hpp_try_keepFe_NTD_original_names.pkout.txt
  mcpb_m7_400ns/
  archive_legacy_YYYYMMDD/
  TJ_WORKFLOW_README.md
```

Current canonical project:

```text
E:/TJ/mcpb_m7_400ns/
  01_inputs/
  02_qm_hessian/
  03_mcpb/
  04_amber/
  05_gromacs/
  06_analysis/
  logs/
  README.md
```

Likely archive candidates after the current MCPB/GAMESS jobs are secure:

```text
616_ligand_inspect/
fe_ligand_building_blocks/
md_analysis/
md_handoff/
prepared/
pymol/
qm_ligand_b3lyp_631gdp/
qm_ligand_b3lyp_631gdp_ACT/
acetate.pdb
nbutamide.pdb
dock_result.pdb
616.pse
616.png
616_2.png
```

Keep as provenance until final MD has started successfully:

```text
feiii_radical_preflight/
hplusplus_inputs/
0616_initial_guess*.pdb
protonation_map*.csv
```

Do not delete old directories while a GAMESS/MCPB job is running. Move them into
an archive directory first, then delete only after the final clean workflow is
reproducible from `mcpb_m7_400ns/`.
