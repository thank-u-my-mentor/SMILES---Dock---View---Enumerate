# Fe-Radical MCPB Workflow Notes

Use these notes when PLE is preparing a non-heme Fe active-site intermediate in
which substrate or product fragments are intentionally coordinated to Fe and one
fragment may be a radical.

## Core Order

1. Start from a curated active-site PDB, not the raw docking receptor.
2. Remove old occupant ligands that block the target site, then add back the
   catalytic Fe and intentional coordinating ligands.
3. Run H++ or manual protonation first, then preserve the resulting Amber names
   for Fe-contact histidines, aspartates, and glutamates.
4. Do a cheap ligand-format preflight with AmberTools before expensive QM:
   antechamber, parmchk2, and a ligand-only `tleap check`.
5. Build the MCPB.py model only after ligand names, atom counts, charges, and
   radical multiplicity are commandable.
6. Screen plausible spin multiplicities with short GAMESS optimizations before
   paying for the final force-constant calculation.

The AmberTools preflight only checks formatting and ordinary GAFF-style ligand
parameters. It does not validate Fe coordination chemistry.

## Six-Coordinate Non-Heme Fe

Prefer a six-coordinate model when the literature, crystal structure, or
mechanism implies an octahedral non-heme Fe center. A typical reactive
intermediate may include:

- two histidine imidazole donors;
- one glutamate or aspartate carboxylate donor;
- one substrate/product N or O donor;
- one co-substrate carboxylate donor;
- one water, dioxygen, superoxo, azide, nitrene, or other mechanism-specific
  sixth ligand.

For Fe(III) high-spin plus an organic radical, do not assume one multiplicity is
automatically correct. Practical first candidates are:

- `MULT=5`: antiferromagnetic coupling between high-spin Fe(III) and the radical;
- `MULT=7`: ferromagnetic coupling between high-spin Fe(III) and the radical.

These multiplicities require an even total electron count. If GAMESS prints
`CHECK YOUR INPUT CHARGE AND MULTIPLICITY`, revise the model charge,
protonation, ligand set, or spin-state candidates before running a long job. In
the TJ 0616 six-coordinate model, `ICHARG=1` gave 271 electrons and was rejected
for `MULT=5/7`; `ICHARG=0` gave 272 electrons and entered the geometry search.

Compare final geometry, Fe-donor distances, spin density, `<S^2>`, and relative
energy. If one multiplicity breaks a coordination bond during optimization, do
not use that model for MCPB parameters.

After both spin-screen logs finish, compare them with:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/compare_gamess_spin_screen.py \
  --workdir /mnt/e/TJ/feiii_radical_preflight/gamess_spin_screen
```

Interpretation order:

1. Both logs must contain `EXECUTION OF GAMESS TERMINATED NORMALLY`.
2. If the optimizer reports `FAILURE TO LOCATE STATIONARY POINT`, the job did
   not crash, but the geometry search hit `NSTEP`; treat the result as a screen,
   not a final Hessian-ready structure.
3. Compare relative energy only between runs with the same model, charge, basis,
   functional, and starting geometry.
4. Reject a formally lower-energy model if Fe-donor distances show a broken
   coordination bond or an implausible ligand rearrangement.
5. Check spin contamination. For multiplicity `M`, expected
   `<S^2> = S(S+1)` where `S=(M-1)/2`; therefore `M5` expects about `6.0` and
   `M7` expects about `12.0`.

For the TJ 0616 `NCORES=16` screen, both `M5` and `M7` finished normally but
both hit the 80-step geometry-search limit. `M5` was lower by only about
`0.12 kcal/mol`, while `M7` had `<S^2>` very close to the ideal septet value.
Both retained six Fe-donor distances around 1.94-2.28 A except the coordinating
water at about 2.12 A. This is a chemically useful screen, but it is not yet the
final MCPB force-constant geometry.

When exporting the chosen spin state, do not blindly use the last coordinate
block in a GAMESS log. If the optimizer hits `NSTEP`, GAMESS prints a "next
predicted" coordinate set whose energy and gradient are unknown. Export the last
evaluated `NSERCH` geometry instead:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/export_gamess_optimized_model.py \
  --workdir /mnt/e/TJ/feiii_radical_preflight/gamess_spin_screen \
  --log logs/FEIII_radical_M7_fastopt.log \
  --full-pdb /mnt/e/TJ/feiii_radical_preflight/complex_AFeIII_UNL_C4radical_dropH4A.pdb \
  --outdir /mnt/e/TJ/feiii_radical_preflight/m7_export
```

This writes:

- `M7_optimized_small_model_gamess_frame.pdb`: evaluated QM geometry in the
  GAMESS coordinate frame;
- `M7_optimized_small_model_aligned_to_initial.pdb`: the same geometry aligned
  back to the initial PDB frame;
- `complex_AFeIII_UNL_C4radical_M7_metalcenter_patch.pdb`: full protein with
  safely mapped Fe-center atoms patched to the M7 aligned coordinates;
- `M7_initial_vs_optimized_compare.tsv` and `M7_export_report.md`: initial-vs-M7
  atom movements and Fe-donor distance changes.

The full-protein patch is a geometry handoff, not a replacement for MCPB.py. Use
it to inspect the M7 active-site model and as the coordinate source for the next
MCPB Hessian/parameterization stage.

If the selected `M7` optimization reaches `NSTEP` before convergence but keeps
reasonable Fe-donor geometry, continue the same spin state from the last
evaluated coordinate rather than restarting from the hand-built PDB:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/prepare_gamess_continuation.py \
  --workdir /mnt/e/TJ/feiii_radical_preflight/gamess_spin_screen \
  --source-log logs/FEIII_radical_M7_fastopt.log \
  --outdir /mnt/e/TJ/feiii_radical_preflight/gamess_m7_cont60 \
  --nstep 60 \
  --charge 0 \
  --mult 7 \
  --job-name FEIII_radical_M7_cont60
cd /mnt/e/TJ/feiii_radical_preflight/gamess_m7_cont60
env NCORES=16 SESSION=feiii_m7_cont60 bash start_continuation_tmux.sh
```

The continuation input is generated from the last evaluated `NSERCH` geometry,
not from GAMESS's unknown "next predicted" coordinates after an NSTEP-limit
warning. Monitor with:

```bash
cd /mnt/e/TJ/feiii_radical_preflight/gamess_m7_cont60
bash monitor_continuation.sh
tmux attach -t feiii_m7_cont60
```

## Radical Ligands

Make the radical explicit in the structure before preflight. If the desired
radical is on a benzyl carbon, remove one H from that carbon and run the ligand
or Fe model as an odd-electron system, usually doublet for ligand-only checks.
Do not model an N-centered radical if the same N is intended to be a neutral
amide-like Fe donor unless the mechanism specifically requires that chemistry.

For a long ligand where radical delocalization matters, keep the full radical
ligand in the QM model. Truncate protein side chains instead:

- His -> methyl-imidazole;
- Glu/Asp -> acetate;
- acetate/co-substrate -> acetate;
- radical ligand -> full ligand when aryl/benzyl spin delocalization matters;
- sixth water -> water.

Free ligand B3LYP optimization is only a geometry hygiene check. It gives the
ligand's relaxed internal conformation without protein, Fe, ligand field, or
pocket steric constraints. It is not a realistic active-site intermediate by
itself.

## Naming And Charges

Avoid `ACE` as a free acetate ligand in Amber input because Amber uses `ACE` for
an N-terminal acetyl cap. Prefer the official PDB CCD residue name `ACT` for
acetate, or another non-reserved three-letter ligand name if the chemistry is not
ordinary acetate.

Use charge and multiplicity deliberately:

- acetate `ACT`: usually charge `-1`, multiplicity `1`;
- neutral amide radical ligand after deleting one C-H: often charge `0`,
  multiplicity `2` for ligand-only preflight;
- full Fe(III) radical model: screen `MULT=5` and `MULT=7`;
- Fe(II) non-radical model: different spin states may apply and should not be
  copied from the Fe(III)-radical case.

## GAMESS Details

GAMESS can replace Gaussian for the QM parts of MCPB.py when Gaussian is not
available. Use the GAMESS syntax for B3LYP/6-31G(d,p):

```text
$CONTRL SCFTYP=UHF DFTTYP=B3LYP RUNTYP=OPTIMIZE ICHARG=... MULT=... $END
$BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=1 $END
```

Do not let `$CONTRL` lines become too long; older GAMESS-style input is sensitive
to card formatting. A completed optimization should contain both
`END OF GEOMETRY SEARCH` and `TERMINATED NORMALLY`.

For MCPB.py force constants, the Hessian/force-constant job is the expensive
step. Run it only after a smaller optimization and spin-state screen look
chemically plausible.

## GAMESS Runtime And CPU Use

For WSL-local GAMESS jobs, prefer a persistent `tmux` session over one-shot
`wsl -e ... nohup ...` launches. One-shot Windows-to-WSL launches can make it
look as if the job stopped when the terminal command returns, while the real
GAMESS state is only visible from the log and `ddikick.x` process. A running
optimization is healthy when the log advances from one-electron integrals to
two-electron integrals and then to the SCF block:

```text
END OF ONE-ELECTRON INTEGRALS
END OF TWO-ELECTRON INTEGRALS
U-B3LYP SCF CALCULATION
ITER EX      TOTAL ENERGY ...
```

Do not call such a job finished until the log contains
`EXECUTION OF GAMESS TERMINATED NORMALLY` or
`EXECUTION OF GAMESS TERMINATED -ABNORMALLY-`.

When running the paired `MULT=5`/`MULT=7` screen, a finished `M5` log does not
mean the whole screen is finished. Check both `.status` files:

```bash
cd /mnt/e/TJ/feiii_radical_preflight/gamess_spin_screen
cat logs/FEIII_radical_M5_fastopt.status
cat logs/FEIII_radical_M7_fastopt.status
tmux capture-pane -t feiii_spin_screen -p -S -80
```

For the TJ 0616 `NCORES=16` run, `M5` finished normally after about 11 h 29 min
48 s wall time and the runner immediately started `M7`. The overall spin screen
is still running while `logs/FEIII_radical_M7_fastopt.status` only has a start
line and `ddikick.x ... FEIII_radical_M7_fastopt ... cpus=16` is present.

CPU use is limited first by `NCORES`, not by the number of logical CPUs in the
machine. On a Ryzen 9950X with 32 logical threads, `NCORES=4` can look like only
about 12-15% total CPU even if those four worker processes are busy. To approach
roughly 50-70% total CPU, start the spin screen with about `NCORES=16`:

```bash
cd /mnt/e/TJ/feiii_radical_preflight/gamess_spin_screen
tmux kill-session -t feiii_spin_screen 2>/dev/null || true
env NCORES=16 SESSION=feiii_spin_screen bash start_fast_spin_screen_tmux.sh
```

Use this only after the input has passed charge/multiplicity and early SCF
checks. More cores do not guarantee linear speedup for DFT, and some GAMESS
stages remain partly serial or I/O-bound. If the job becomes unstable or slower,
drop to `NCORES=8` or `NCORES=12`. Keep `MWORDS` moderate; requesting excessive
replicated memory per process can trigger silent WSL or DDI instability.

## TJ 0616 Preflight Pattern

For the TJ hand-built intermediate, the intended six-coordinate Fe site was:

- `ACT A501 O2`
- `GLU A349 OE1`
- `HOH B875 O`
- `UNL A1 N1`
- `HIS A187 NE2`
- `HIS A270 NE2`

The current radical correction is to remove one hydrogen from `UNL A1 C4`, e.g.
`H4A`, then run ACT as charge `-1` singlet and UNL as charge `0` doublet for
ligand-only AmberTools preflight. A warning from antechamber that radical bond
types may be wrong is expected and must be manually reviewed before final MCPB.

The bundled PLE helper `scripts/prepare_fe_radical_preflight.py` generated:

```text
/mnt/e/TJ/feiii_radical_preflight/
  complex_AFeIII_UNL_C4radical_dropH4A.pdb
  ligands/ACT_A501.pdb
  ligands/UNL_A1_original.pdb
  ligands/UNL_A1_C4radical.pdb
  metal_site/chainA_fe_donors.tsv
  amber_preflight/run_amber_ligand_preflight.sh
```

The current helper is still a command-line script, not a top-level `ple`
subcommand. Promote it to the CLI only after the radical atom, ligand charges,
and donor list are expressed as stable general parameters.

After preflight, use `scripts/prepare_fe_radical_gamess_screen.py` to generate a
truncated Fe(III)-radical screening model and two GAMESS inputs:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/prepare_fe_radical_gamess_screen.py \
  --pdb /mnt/e/TJ/feiii_radical_preflight/complex_AFeIII_UNL_C4radical_dropH4A.pdb \
  --outdir /mnt/e/TJ/feiii_radical_preflight/gamess_spin_screen \
  --charge 0 \
  --nstep 80
```

It writes:

```text
gamess_spin_screen/
  feiii_radical_truncated_model.pdb
  feiii_radical_truncated_model.xyz
  FEIII_radical_M5_fastopt.inp
  FEIII_radical_M7_fastopt.inp
  run_fast_spin_screen.sh
  monitor_fast_spin_screen.sh
  WORKFLOW_FE_RADICAL_MCPB.md
```

The generated workflow markdown contains a Mermaid flowchart and an explanation
of how MCPB.py turns QM-derived metal-site geometry/force constants into Amber
bonded Fe parameters for later GROMACS MD.

## TJ M7 MCPB Clean Project Lessons

For the TJ 0616 Fe(III)-radical intermediate, the clean project directory is:

```text
/mnt/e/TJ/mcpb_m7_400ns/
```

The key working command is:

```bash
cd /mnt/e/TJ/mcpb_m7_400ns/03_mcpb
bash start_mcpb_qm_tmux.sh
bash monitor_mcpb_gamess_qm.sh
tmux attach -t tj_mcpb_m7_qm
```

Important fixes that must be preserved in future generalized PLE code:

1. Use a chain-A-only MCPB PDB for this homodimeric enzyme unless the metal site
   explicitly spans chains. `pymsmt`/MCPB indexes residues internally by residue
   id and can merge chain A and chain B residues with the same number, which
   corrupts side-chain truncation.
2. Keep the coordinating water as the sixth ligand. In this project the crystal
   water was `HOH B875`; for the chain-A-only MCPB PDB it is remapped to
   `HOH A875`, and H1/H2 coordinates are copied from the M7 small model.
3. Fe-coordinating histidines that donate through `NE2` are written as `HID` in
   the Amber/MCPB PDB. Without explicit HID/HIE/HIP names, MCPB may not find
   charges in the Amber residue library.
4. If the H++/curated PDB is heavy-only, MCPB side-chain truncation may require
   an approximate `HA` on coordinating residues. Add only the minimum atoms
   needed for MCPB capping, and document that these are modeling helper atoms.
5. The Fe residue name must match the ion mol2 residue name. Here the PDB Fe is
   normalized from `FE2` to residue `FE`, matching `FE.mol2`.
6. Include all non-standard coordinating residues in `naa_mol2files`:

```text
naa_mol2files ACT.mol2 UNL.mol2 HOH.mol2
frcmod_files ACT.frcmod UNL.frcmod
```

7. The custom 70-atom M7 screen and the MCPB-generated small model are not the
   same model. The screen identifies a chemically plausible spin state and
   geometry, but the final MCPB force constants should be computed on MCPB's
   own generated small model. For the current MCPB small/large models,
   `charge=+1`, `MULT=7` gives even electron counts:

```text
small model: 62 atoms, 288 electrons
large model: 92 atoms, 456 electrons
```

8. MCPB-generated `*_small_fc.inp` can have an empty `$DATA` block. The normal
   workflow is: run `*_small_opt.inp`, extract optimized coordinates, fill
   `*_small_fc.inp`, then run Hessian. The helper
   `scripts/write_mcpb_gamess_pipeline.py` automates this.

## ZFY 2R5V Intermediate Small-Optimization Pattern

For a hand-built 2R5V-derived ZFY Fe-radical intermediate, use the clean MCPB
scaffold rather than the raw PyMOL-edited PDB:

```text
/mnt/e/ZFY/zfy_mcpb_scaffold/ZFY_2R5V_MCPB_start_modify_NH_clean_for_MCPB.pdb
```

The intended first-pass six-coordinate model is:

```text
HID A161 NE2
HID A241 NE2
GLU A320 OE1
UNK A501 N1
UNK A501 O2
ACT A502 O1
```

Generate the compact MCPB small-optimization project with:

```bash
/mnt/l/WSL/softwares/conda_envs/md/bin/python \
  /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/prepare_fe_radical_mcpb_opt_project.py \
  --pdb /mnt/e/ZFY/zfy_mcpb_scaffold/ZFY_2R5V_MCPB_start_modify_NH_clean_for_MCPB.pdb \
  --outdir /mnt/e/ZFY/mcpb_m7_small_opt \
  --group ZFY_M7_FE \
  --metal-chain A \
  --metal-resseq 4113 \
  --metal-atom FE \
  --cutoff 2.8 \
  --charge 1 \
  --mult 7 \
  --cores 16 \
  --session zfy_mcpb_m7_opt \
  --ligand UNK /mnt/e/ZFY/zfy_mcpb_scaffold/UNK_NH_clean_for_MCPB.mol2 /mnt/e/ZFY/zfy_mcpb_scaffold/UNK_NH_clean_for_MCPB.frcmod \
  --ligand ACT /mnt/e/ZFY/zfy_mcpb_scaffold/ACT_for_MCPB.mol2 /mnt/e/ZFY/zfy_mcpb_scaffold/ACT_for_MCPB.frcmod
```

Then run only the Hessian-before geometry optimization:

```bash
cd /mnt/e/ZFY/mcpb_m7_small_opt/mcpb
bash run_01_mcpb_step1.sh
bash run_02_prepare_gamess_inputs.sh
env NCORES=16 bash run_03_small_opt_tmux.sh
tail -f gamess_logs/ZFY_M7_FE_small_opt.log
```

Do not launch `*_small_fc` Hessian until `ZFY_M7_FE_small_opt.log` terminates
normally and the Fe-donor distances are inspected. In the first ZFY run, MCPB
step 1 found a 55-atom small model with 276 electrons and an 85-atom large model
with 444 electrons, so `charge=+1, MULT=7` was electron-parity consistent.
