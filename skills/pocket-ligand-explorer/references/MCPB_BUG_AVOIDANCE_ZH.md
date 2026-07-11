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

Hard rule for MCPB: H++ must not override the first Fe coordination sphere.
Before `MCPB.py -s 1`, re-audit every residue with a donor atom within the metal
cutoff and force the Amber protonation/name from the metal geometry and
mechanism:

```text
Fe-bound His via NE2 -> neutral HID: ND1-H present, NE2-H absent
Fe-bound His via ND1 -> neutral HIE: NE2-H present, ND1-H absent
HIP near Fe          -> only keep if a positively charged, non-donor His is
                        chemically intended; do not let HIP protonate the donor N
Fe-bound Cys         -> often CYM in MCPB examples, not neutral CYS, when it is
                        deprotonated/thiolate-bound
Fe-bound Asp/Glu     -> usually deprotonated coordinating carboxylate unless the
                        mechanism explicitly needs ASH/GLH
```

This is why PLE runs `fix_fe_histidine_protonation.py` and
`audit_mcpb_histidine_protonation.py` before MCPB. H++ is still useful for the
ordinary protein background, but metal-site residues are curated by the
coordination model. A H++-produced `HIP` or a histidine with no imidazole N-H in
the Fe first shell is a red flag, not an instruction to trust blindly.

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

GAMESS input cards are sensitive to long lines. Do not put a crowded `$CONTRL`
or `$SCF` group on one long line when adding `RUNTYP`, `ICHARG`, `MULT`,
`MAXIT`, `COORD`, and `UNITS`. If the log says:

```text
**** ERROR READING INPUT GROUP $CONTRL *****
THE PROBLEM IS WITH THIS INPUT LINE, NEAR THE X MARKER
$SCF ...
```

then `$CONTRL` was likely truncated before `$END`, so the next `$SCF` line was
mistakenly read as part of `$CONTRL`. Use multi-line groups:

```text
$CONTRL
 SCFTYP=UHF DFTTYP=B3LYP RUNTYP=ENERGY
 ICHARG=1 MULT=7 MAXIT=200 COORD=UNIQUE UNITS=ANGS
$END
$SCF
 DIRSCF=.T. DIIS=.T. DAMP=.T. SOSCF=.T.
 ETHRSH=2.0 MAXDII=20
$END
```

In this local GAMESS build, `$CONTRL MAXIT` must be between 0 and 200. Also
avoid `COORD=CART` for continuation-style metal-site work: GAMESS warns that it
will rotate coordinates to principal axes. Use `COORD=UNIQUE UNITS=ANGS` when
you want the input coordinates to remain the coordinate frame used by the later
PDB/MCPB audits.

## 7.1 Fe-Radical GAMESS Grid Law

For non-heme Fe radical/intermediate GAMESS jobs in PLE, use the fine DFT grid
from the first SCF iteration. Treat this as the default rule, not as a tunable
performance option:

```text
$DFT
 NRAD=96 NLEB=302 NRAD0=96 NLEB0=302 SWOFF=0.0
$END
```

Reason: several ZFY Fe-radical runs looked stable on GAMESS' cheaper early
grid, then collapsed when the calculation switched back to the fine grid. The
failure signature was:

```text
DFT CODE IS SWITCHING BACK TO THE FINE GRID
energy jumps by hundreds or thousands of Hartree
DENSITY CHANGE jumps to tens or hundreds
SCF oscillates or never reaches DENSITY CONVERGED
```

The correct interpretation is not "the chemistry almost converged". It means
the electronic state was only stable on the temporary coarse grid. Do not start
`GRADIENT`, `OPTIMIZE`, Hessian, or MCPB force-constant extraction from such a
vector.

Durable PLE rule:

```text
1. Run fixed-geometry ENERGY/SCF.
2. Reuse the converged $VEC with GUESS=MOREAD.
3. Use the actual NORB from the current GAMESS log; never reuse an old hardcoded
   value such as 630 from another model.
4. Run GRADIENT and OPTIMIZE with the fine grid enabled from the first step.
5. Accept a wavefunction only after DENSITY CONVERGED appears and no
   SCF IS UNCONVERGED / SCF HAS NOT CONVERGED marker appears.
6. Reject early electronic cliffs. If the first 10-25 SCF iterations show the
   total energy becoming tens to hundreds of Hartree less stable while density
   change stays above about 0.02 and orbital gradient stays high, stop the job.
   Do not wait for a miracle convergence, and never pass that $VEC to OPT.
```

For a brand-new Fe-radical spin branch, the first fresh SCF should use the
slowest stable mode:

```text
prepare_gamess_fixed_energy.py --scf-mode damp

$SCF
 DIRSCF=.T. DIIS=.F. SOSCF=.F. DAMP=.T. SHIFT=.T.
$END
```

If the ordinary `damp` branch still lets the Fe-radical wavefunction jump around
or repeatedly falls into a DIIS-like recovery path, use the stricter emergency
mode:

```text
prepare_gamess_fixed_energy.py --scf-mode very-damp

$SCF
 DIRSCF=.T. DIIS=.F. SOSCF=.F. DAMP=.T. SHIFT=.T.
 ETHRSH=0.0 MAXDII=1
$END
```

This is intentionally slower and may need many iterations. Its purpose is not
speed; it is to get a fixed-geometry `DENSITY CONVERGED` wavefunction without
coarse-grid shortcuts or aggressive extrapolation. Use it before trying OPT on a
new Fe-radical/nitrene branch, especially after changing ligand hydrogens,
histidine protonation, or the formal spin multiplicity.

Only consider SOSCF or DIIS after one of these is true:

```text
1. A matching model already has a clean converged $VEC and the correct NORB.
2. The damp-only branch shows a sane, monotonic-enough trend but is too slow.
3. You are deliberately reproducing an old benchmark log and record that choice.
```

Do not use a formally parity-allowed but collapsing spin branch just because it
ran for many iterations. The TJ repaired-H model made this explicit:
`charge=+1, MULT=2` was legal for 303 electrons, but chemically/electronically
unusable.

Do not overinterpret a single GAMESS banner such as:

```text
* * * INITIATING DIIS PROCEDURE * * *
```

In this local build, that phrase can still appear even when the input contains
`DIIS=.F.`. Treat the input card as the requested strategy, but judge the run by
the iteration table. A valid branch must not show early energy/density cliffs.
Use the guard script for fresh Fe-radical SCF attempts:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/guard_gamess_scf.py \
  --log gamess_logs/<job>.log \
  --kill-pattern <job> \
  --kill
```

The helper generated by `prepare_gamess_fixed_energy.py` writes
`run_<job>_guarded.sh`; prefer that for long fresh-SCF attempts so bad
wavefunctions are stopped automatically.

The guard must not be too impatient. For fresh HUCKEL starts on Fe-radical
models, the first 5-10 SCF iterations can swing violently before the density
settles. The default guard now waits for at least 25 parsed SCF rows unless
there is catastrophic divergence. A branch that needs 100-300 iterations to
converge is acceptable if the density/metric trend eventually improves and no
`SCF IS UNCONVERGED` marker appears.

This is slower per SCF iteration, but it avoids wasting hours on a route that
only fails at the coarse-to-fine switch. In plain language: start with the
"high-resolution map" even though each step is more expensive, because the route
is stable.

The 2026-06-29 TJ M2 verification is the canonical negative example. After
protein-H repair changed the electron parity, `charge=+1, MULT=2` was formally
allowed, but its fresh fine-grid HUCKEL/SOSCF run immediately drifted from about
`-2995` Hartree to about `-2808` Hartree by iteration 24. Density and orbital
gradient remained large and damping exploded. This is not "slow progress"; it is
an electronic-state cliff. Such a run must be stopped and recorded as a failed
spin/electronic-state branch.

The PLE MOREAD helpers now enforce this by default for new Fe-radical work:

```text
prepare_gamess_fixed_energy.py       -> fresh HUCKEL fixed-geometry SCF, fine grid from start
prepare_gamess_moread_energy.py       -> fine grid by default
prepare_gamess_moread_gradient.py     -> fine grid always
prepare_gamess_moread_optimize.py     -> fine grid always
prepare_gamess_moread_harmonic_optimize.py -> fine grid always
```

If the previous SCF did not reach `DENSITY CONVERGED`, do not use its `$VEC`
for `MOREAD`. Start a clean fixed-geometry SCF instead:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/prepare_gamess_fixed_energy.py \
  --mcpb-dir /path/to/mcpb \
  --source-job TJ_HID_M7_FE_small_opt \
  --out-job TJ_HID_M7_FE_scf_energy_fresh_finegrid_damp \
  --charge 1 --mult 7 --scf-mode damp --cores 8
```

This keeps the MCPB small-model coordinates, changes `RUNTYP=ENERGY`, removes
`$STATPT`, `$VEC`, `$ELPOT`, and `$PDC`, uses `GUESS=HUCKEL`, and forces
`NRAD0=96 NLEB0=302 SWOFF=0.0`. It is slower than a good `MOREAD`, but much
safer than inheriting an unconverged or coarse-grid-contaminated wavefunction.

Only use a coarse-grid path to reproduce an old log. If a script exposes an
escape hatch such as `--allow-coarse-grid`, treat it as a debugging override,
not a production workflow.

## 7.2 Guarded Geometry Optimization For Radical Ligands

Do not launch a free `NSTEP=1000` OPT immediately after the first fixed-geometry
SCF converges. A converged SCF only says the electronic state is self-consistent
at the current coordinates; it does not guarantee that an unconstrained QM
geometry optimization will preserve the intended covalent connectivity of a
hand-built radical intermediate.

ZFY lesson:

```text
The bidentate UNK model reached a formal OPT endpoint, but the optimized small
model broke the UNK O1-C9 covalent bond. GAMESS did not make a software mistake:
the QM Hamiltonian is not constrained by the mol2 bond table, so it can optimize
toward bond cleavage or a different reaction-like structure if the supplied
intermediate is chemically unstable.
```

Safer first OPT protocol:

```text
1. Run fixed-geometry SCF first.
2. Start OPT from the SCF $VEC via GUESS=MOREAD.
3. Use short OPT segments first: NSTEP=10-20, not NSTEP=1000.
4. Keep harmonic restraints on all Fe first-shell donor distances.
5. Also protect Fe first-shell protein scaffold bonds that must not drift:
   His `CA-CB`, `CB-CG`, and imidazole ring bonds; Glu/Asp `CA-CB`,
   `CB-CG`, `CG-CD`, and carboxylate `CD-OE/OD` bonds. This avoids a repeat of
   the TJ-style failure where a histidine-like custom residue kept a formal Fe
   contact but lost chemically valid backbone-sidechain geometry.
6. Also protect nonreactive ligand internal bonds that must not break, especially
   carbonyl/amide/ester bonds near Fe and radical-bearing scaffold bonds.
7. After every short segment, export an aligned PDB and audit ligand geometry
   against the mol2 bond table.
8. Hard-stop before Hessian if any mol2 bond is too long, too short, missing, or
   changed by more than the chosen tolerance.
```

Use `prepare_gamess_moread_harmonic_optimize.py` for this guarded first segment.
It can write both metal-donor and ligand-internal harmonic distance restraints:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/prepare_gamess_moread_harmonic_optimize.py \
  --mcpb-dir /path/to/mcpb \
  --source-job ZFY2_HOH_M7_FE_small_opt \
  --out-job ZFY2_HOH_M7_FE_small_opt_guarded10 \
  --dat /home/qin/softwares/gamess/restart/ZFY2_HOH_M7_FE_small_scf_energy_diis.dat \
  --small-pdb small_model_named_for_constraints.pdb \
  --metal FE:A:4113:FE \
  --donor HID:A:161:NE2 \
  --donor HID:A:241:NE2 \
  --donor GLU:A:320:OE1 \
  --donor UNK:A:501:N1 \
  --donor ACT:A:502:O2 \
  --donor HOH:A:503:O \
  --protect-bond UNK:A:501:O1,UNK:A:501:C9 \
  --protect-bond UNK:A:501:N1,UNK:A:501:C9 \
  --protect-bond UNK:A:501:C9,UNK:A:501:O2 \
  --charge 1 --mult 7 --nstep 10 --scf-mode diis --cores 8
```

For the first ZFY monodentate-UNK/HOH test, the practical starting force constant
was `500` in the GAMESS harmonic-restraint input for Fe-donor, His/Glu scaffold,
and key ligand bonds. Treat this as a short-segment guardrail, not as a final MD
force-field parameter. After the first 10-20 steps, reduce restraints gradually
or rerun a less constrained continuation only if the bond audit stays clean.

This is not meant to fake the final force field. It is a chemical guardrail for
the first few optimization steps. If the guarded OPT can only stay intact because
very strong ligand/protein-internal constraints are present, that is a warning
that the hand-built intermediate may not be a stable QM structure. In that case
do not proceed blindly to Hessian; revise the coordination model, protonation,
charge, or radical placement.

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

Before Hessian, also audit every coordinating ligand against its `mol2` bond
table. GAMESS does not know or preserve the `mol2` bond list; in a strained
Fe-radical model it can converge by breaking an internal ligand covalent bond.
This is not a file-format error, but it makes the geometry unusable for MCPB
force constants unless the goal is to study that bond-breaking reaction itself.

Also audit Fe-bound protein side chains inside the MCPB small model. A
coordinating Glu/Asp can be truncated to a capped side-chain model such as
`CH3-CB-CG-CD-OE*`. If only the ligand or Fe-ligand distance is protected, the
unprotected carboxylate side chain can collapse into an unphysical almost-linear
geometry. A TJ Fe-nitrene retry showed this hard-stop pattern:

Important distinction: this does not mean "every missing hydrogen is wrong".
MCPB.py intentionally truncates protein residues in the small model and caps
the cut site. A side chain extracted from the protein may become a methyl-capped
fragment rather than a full amino acid with backbone atoms. That is normal for
MCPB. The hard rule is narrower:

```text
Positions that should have H in the capped small model must have H.
Positions that should be unprotonated, such as Fe-donor His NE2/ND1 or
deprotonated carboxylate OE/OD atoms, must not be counted as missing-H errors.
```

Therefore after `MCPB.py -s 1`, run:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/audit_mcpb_small_model_hydrogens.py \
  --small-pdb TJ2nit_small.pdb \
  --mol2-dir /path/to/mcpb \
  --out logs/mcpb_small_model_hydrogen_audit.tsv \
  --max-warnings 2
```

The default policy is:

```text
warning_count <= 2  -> continue only after review
warning_count >  2  -> hard error; do not prepare OPT or Hessian
```

This protects against the exact failure where a GLU small-model fragment has
missing hydrogens on `CH3`, `CB`, or `CG`, while avoiding false positives for
legitimately unprotonated Fe donor atoms or deprotonated carboxylate oxygens.
`tleap` can add hydrogens later for a final Amber topology, but it cannot repair
a GAMESS small model that was already optimized or Hessian-calculated with the
wrong valence.

自动修复也必须有边界。PLE 的修复策略是：

```text
MCPB step1 前：
  可以自动修标准蛋白残基 ATOM，例如 Fe-bound His 的 HID/HIE 互变异构、
  GLU/ASP 的 CA/CB/CG 脂肪族氢。

永远不要自动修：
  UNL/UNK/ACT/nitrene/radical/substrate 等 HETATM 或非标准底物。
  这些分子的 H 数、自由基位置、Fe=N 或 Fe-N 配位形式必须由人工和 mol2
  模板共同定义。

MCPB step1 后：
  先生成 `<group>_small_visual_check.pdb` 给 PyMOL 人眼检查；
  再运行 small-model hydrogen/electron audit。
```

`build_mcpb_visual_check_pdb.py` 会把 MCPB small model、ligand mol2 键表、
以及 Fe-donor 距离连接写进一个只用于可视化的 PDB。这个文件不是力场文件，
也不会改变 MCPB/GAMESS/Amber 使用的化学参数；它的意义是让人类在 OPT 或
Hessian 前最后看一眼小模型是不是少氢、断键、配位数错、底物放错。

氢原子数量会通过电子数影响 multiplicity：

```text
total_electrons = sum(atomic_numbers) - model_charge
even electrons -> allowed MULT is odd: 1, 3, 5, 7, ...
odd electrons  -> allowed MULT is even: 2, 4, 6, ...
```

所以“补一个本该存在的蛋白氢”可能会让以前看似能跑的 `MULT=7` 立刻变成
奇偶不匹配。这不是坏事，而是说明旧模型的电子数可能是被缺氢错误骗出来的。

```text
initial small GLU349:
  CH3-CB-CG angle = 115.60 deg
  CB-CG-CD angle  = 109.76 deg
  CB-CG           = 1.6655 A

optimized small GLU349:
  CH3-CB-CG angle = 179.42 deg
  CB-CG-CD angle  = 174.36 deg
  CB-CG           = 1.2354 A
```

This is not a valid glutamate-like geometry and must not be used for Hessian.
It can look like an allene/double-bond artifact in PyMOL, but the practical
interpretation is simpler: the small-model optimization has left the intended
protein-side-chain chemistry. Stop the Hessian, discard that OPT geometry, and
rerun with explicit protection of the Fe-bound acidic side-chain internal bonds:

```text
Glu/Asp scaffold guards:
  cap/CA-CB
  CB-CG (or CB-CG/CD depending on residue)
  CG-CD
  CD-OE1 / CD-OE2
  Fe-O donor distance
```

Do not judge this only after patching back into the full protein. The optimized
small model itself must pass these side-chain geometry audits before its Hessian
is useful.

Example audit:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/audit_ligand_mol2_bonds_in_pdb.py \
  --mol2 /mnt/e/ZFY/mcpb_project/input/UNK.mol2 \
  --pdb /mnt/e/ZFY/mcpb_project/opt_compare/optimized_small_model_aligned.pdb \
  --reference-pdb /mnt/e/ZFY/mcpb_project/opt_compare/initial_small_model.pdb \
  --resname UNK \
  --out /mnt/e/ZFY/mcpb_project/opt_compare/UNK_mol2_bond_audit.tsv
```

Hard-stop criteria:

```text
BROKEN_LONG_BOND  -> reject geometry for Hessian
TOO_SHORT         -> reject geometry for Hessian
large ligand rearrangement around a nonreactive covalent bond -> review manually
```

ZFY lesson: the bidentate UNK model converged electronically and geometrically,
but `UNK O1-C9` stretched from about `1.52 A` to `3.41 A`. That means the
optimization found a bond-broken ligand state. The next safer model removes UNK
O2 as a required donor and restores six coordination with a water ligand:

```text
HID A161 NE2
HID A241 NE2
GLU A320 OE1
UNK A501 N1
ACT A502 O2
HOH A503 O
```

For GAMESS parallelism, more cores are not always faster. Start with 8-16 cores
on a strong desktop CPU; if DDI errors or poor scaling appear, reduce cores.
For a stable Fe-radical OPT that is already reducing gradients smoothly, do not
kill the current segment just to change the core count. Let the current short
`NSTEP` segment finish, then prepare a continuation from the last evaluated
`NSERCH` geometry plus the newest GAMESS `.dat` `$VEC` block. A practical
desktop pattern is:

```text
1. Run a cautious first MOREAD OPT segment, e.g. NSTEP=20, NCORES=8.
2. If gradients steadily decrease and Fe-donor distances remain intact, prepare
   a faster continuation, e.g. NSTEP=60-80, NCORES=12-16.
3. Before Hessian, require a real optimized geometry or at least a documented
   small-gradient plateau with stable Fe-donor distances.
```

Do not make a "continuation" that only reuses `$VEC` while returning to the old
starting coordinates. The continuation must use the last evaluated coordinate
block from the GAMESS log, not the unknown "next predicted" coordinates after an
`NSTEP` limit, and not the original MCPB small-opt geometry.

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

When an unconstrained MOREAD optimization still fails but the fixed-geometry
SCF is already converged, use a chemically constrained retry before changing the
model again. For the ZFY `5_for_tleap.pdb` M5 case, the active constrained route
is:

```bash
/mnt/l/WSL/softwares/conda_envs/md/bin/python \
  /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/prepare_gamess_moread_harmonic_optimize.py \
  --mcpb-dir /mnt/e/ZFY/mcpb_5for_tleap_M5_try/mcpb \
  --source-job ZFY5_M5_FE_small_opt \
  --out-job ZFY5_M5_FE_small_opt_moread1_harm6_nstep80 \
  --dat /home/qin/softwares/gamess/restart/ZFY5_M5_FE_scf_energy_soscf_moread1.dat \
  --small-pdb /mnt/e/ZFY/mcpb_5for_tleap_M5_try/mcpb/ZFY5_M5_FE_small.pdb \
  --metal FE:A:4113:FE \
  --donor HID:A:161:NE2 \
  --donor HID:A:241:NE2 \
  --donor GLU:A:320:OE1 \
  --donor HOH:A:500:O \
  --donor UNK:A:501:N1 \
  --donor ACT:A:502:O1 \
  --charge 1 \
  --mult 5 \
  --norb 630 \
  --nstep 80 \
  --opttol 0.0002 \
  --force 500 \
  --cores 4 \
  --session zfy5_m5_harm6_n80

cd /mnt/e/ZFY/mcpb_5for_tleap_M5_try/mcpb
env NCORES=4 SESSION=zfy5_m5_harm6_n80 \
  bash run_ZFY5_M5_FE_small_opt_moread1_harm6_nstep80_tmux.sh
tail -f gamess_logs/ZFY5_M5_FE_small_opt_moread1_harm6_nstep80.log
```

This path uses the converged `GUESS=MOREAD` vector, keeps `DIIS=.F.` and
`SOSCF=.T.`, and adds harmonic bond restraints between Fe and all six first-shell
donors. It does not freeze the entire His/Glu scaffold; instead it prevents
Fe-His/Glu/UNK/ACT/HOH coordination from flying apart while allowing the small
model to relax. The helper reads the atom order from `ZFY5_M5_FE_small.pdb` and
writes a TSV with the exact GAMESS atom indices. In the current ZFY project the
restraints are:

```text
Fe#56-HID161 NE2#9   target 2.1006 A
Fe#56-HID241 NE2#18  target 2.0597 A
Fe#56-GLU320 OE1#25  target 2.0860 A
Fe#56-HOH500 O#27    target 2.0634 A
Fe#56-UNK501 N1#29   target 2.0435 A
Fe#56-ACT502 O1#50   target 2.0941 A
```

Use `NCORES=4` for this retry when a separate TJ GROMACS production is running.
Four GAMESS workers should look like roughly 12-15% total CPU on a 32-thread
Ryzen desktop; this is intentional resource sharing, not a failure. The log must
show both of these early markers before the run is considered launched:

```text
HARMONIC RESTRAINTS: BONDS
GUESS =MOREAD
```

This constrained retry was tested on 2026-06-22/23 and did not produce a usable
optimized geometry:

```text
Job: ZFY5_M5_FE_small_opt_moread1_harm6_nstep80
Wall time: about 2 h 29 min 53 s
NCORES: 4
Result: EXECUTION OF GAMESS TERMINATED -ABNORMALLY-
Failure marker: SCF IS UNCONVERGED, TOO MANY ITERATIONS
Failure marker: FAILURE TO LOCATE STATIONARY POINT, SCF HAS NOT CONVERGED
```

Important interpretation: the harmonic restraints themselves worked. At the
failure point the restraint report showed essentially zero displacement from
all six target Fe-donor distances:

```text
Fe-HID161 NE2   delta -0.00001 A
Fe-HID241 NE2   delta -0.00003 A
Fe-GLU320 OE1   delta -0.00004 A
Fe-HOH500 O     delta -0.00004 A
Fe-UNK501 N1    delta -0.00002 A
Fe-ACT502 O1    delta -0.00004 A
```

Therefore this failure should not be interpreted as "the ligands flew away".
It means the current Fe-radical electronic model still cannot obtain a stable
SCF solution at the first geometry point under this B3LYP/SOSCF/MOREAD setting.

To separate "optimizer problem" from "electronic model problem", a fixed-geometry
gradient diagnostic was generated with:

```bash
/mnt/l/WSL/softwares/conda_envs/md/bin/python \
  /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/prepare_gamess_moread_gradient.py \
  --mcpb-dir /mnt/e/ZFY/mcpb_5for_tleap_M5_try/mcpb \
  --source-job ZFY5_M5_FE_small_opt \
  --out-job ZFY5_M5_FE_moread1_gradient \
  --dat /home/qin/softwares/gamess/restart/ZFY5_M5_FE_scf_energy_soscf_moread1.dat \
  --charge 1 \
  --mult 5 \
  --norb 630 \
  --cores 4 \
  --session zfy5_m5_moread1_gradient
```

This job used the same coordinates and the same converged MOREAD vector but
`RUNTYP=GRADIENT`, so it did not move any atoms. It also diverged during the SCF
before producing a usable gradient: by iterations 8-13, the energy oscillated
from about `-945` to `+2769` a.u. and the density change rose to roughly
`1333`. The diagnostic was stopped manually. This strongly suggests the current
ZFY Fe-radical electronic model is not robust enough for Hessian/MCPB force
constants, even though one previous fixed-geometry ENERGY job found a formally
converged density. Do not use the `ZFY5_M5_FE_moread1_gradient` log for Hessian,
MCPB force constants, or MD parameters.

The next ZFY decision should be chemical/electronic-model triage, not another
blind geometry optimization with the same model:

```text
1. Re-check total charge and radical/proton assignment.
2. Re-check whether UNK should be monodentate or bidentate.
3. Re-check whether ACT O1 or O2 is the intended Fe donor.
4. Consider a smaller first-shell QM model before adding the full UNK aryl tail.
5. Consider ORCA/Gaussian-style broken-symmetry or a different functional only
   after the chemical model is explicit.
```

On 2026-06-23, a revised ZFY `5_for_tleap.pdb` model removed the first-shell
water and changed the intended Fe first shell to UNK bidentate coordination:

```text
Input:
  /mnt/e/ZFY/zfy_mcpb_scaffold/tleap_preflight_5_pdbcoords_mol2chem/5_for_tleap.pdb

Project:
  /mnt/e/ZFY/mcpb_5for_tleap_bidentate_M5_scf

Fe first-shell donors within 2.8 A:
  ACT A502 O2      1.834 A
  UNK A501 O2      1.914 A
  HID A241 NE2     2.042 A
  UNK A501 N1      2.118 A
  GLU A320 OE1     2.167 A
  HID A161 NE2     2.171 A
```

MCPB.py step 1 reached the small-model generation before failing in the
standard/large model with a ligand-library naming mismatch:

```text
KeyError: 'UNK-H01'
```

Interpretation: for the immediate "can this electronic model even self-consist?"
question, the generated small model is still useful. Do not proceed to full
MCPB step 2/3/4 until the `UNK` PDB-vs-mol2 hydrogen names are made consistent
for the standard/large model. The generated small model had:

```text
small model: 55 atoms, 276 electrons
charge=+1, MULT=5
NORB=615 inferred from SOSCF rotation counts
```

The active SCF smoke test is fixed-geometry B3LYP/6-31G(d), not optimization:

```bash
cd /mnt/e/ZFY/mcpb_5for_tleap_bidentate_M5_scf/mcpb
tail -f gamess_logs/ZFY5B_M5_FE_scf_energy_soscf.log
cat gamess_logs/ZFY5B_M5_FE_scf_energy_soscf_to_ZFY5B_M5_FE_scf_energy_soscf_moread1.chain_status
```

Because the first SOSCF ENERGY run may only generate a restart vector, a watcher
was launched to automatically run a second fixed-geometry MOREAD ENERGY if the
first run does not print `DENSITY CONVERGED`:

```text
first job:  ZFY5B_M5_FE_scf_energy_soscf
second job: ZFY5B_M5_FE_scf_energy_soscf_moread1
watcher tmux session: zfy5b_m5_scf_chain
```

This uses the new helper:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/run_gamess_moread_energy_chain.py ...
```

Only if one of these fixed-geometry ENERGY jobs reaches `DENSITY CONVERGED`
without `SCF IS UNCONVERGED` should the next diagnostic be `RUNTYP=GRADIENT`.
If the fixed-geometry gradient is also stable, then attempt a tightly controlled
optimization. Do not jump directly from this smoke test to Hessian.

To fix the `UNK-H01` class of errors, do not hand-edit random hydrogen names in
PyMOL. Treat the PDB as coordinates and the mol2 as the naming/chemistry
template. Ignore misleading PDB `CONECT` records and rename ligand hydrogens by
their nearest heavy-atom parent so PDB names match the mol2 library:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/rename_pdb_ligand_atoms_from_mol2.py \
  --pdb /mnt/e/ZFY/zfy_mcpb_scaffold/tleap_preflight_5_pdbcoords_mol2chem/5_for_tleap.pdb \
  --mol2 /mnt/e/ZFY/zfy_mcpb_scaffold/tleap_preflight_5_pdbcoords_mol2chem/UNK.mol2 \
  --resname UNK \
  --parent-mode distance \
  --out /mnt/e/ZFY/zfy_mcpb_scaffold/tleap_preflight_5_pdbcoords_mol2chem/5_for_tleap_UNKmol2names.pdb
```

For the 2026-06-23 ZFY `UNK` template, this maps:

```text
C1: H01 -> H02
C1: H02 -> H03
C2: H05 -> H08
C4: H07 -> H09
C5: H08 -> H10
C6: H09 -> H11
C7: H10 -> H12
C8: H11 -> H13
```

Using the fixed PDB, MCPB.py step 1 successfully generates all three models:

```text
small model:    55 atoms, 276 electrons
standard model: 58 atoms
large model:    85 atoms, 444 electrons
```

After ACT was moved farther from Fe, the current ZFY bidentate smoke-test path
is:

```text
Input:
  /mnt/e/ZFY/zfy_mcpb_scaffold/5_for_tleap.pdb

Fixed-name PDB:
  /mnt/e/ZFY/zfy_mcpb_scaffold/5_for_tleap_UNKmol2names_20260623_140354.pdb

Project:
  /mnt/e/ZFY/mcpb_5for_tleap_bidentate_M5_scf_20260623_140354

Charge/spin:
  charge=+1, MULT=5

Fe first-shell donors within 2.8 A:
  UNK A501 O2      1.914 A
  HID A241 NE2     2.042 A
  ACT A502 O2      2.117 A
  UNK A501 N1      2.118 A
  GLU A320 OE1     2.167 A
  HID A161 NE2     2.171 A

MCPB.py step 1:
  small model: 276 electrons
  large model: 444 electrons
  ORBITAL PRINTING OPTION NORB=615
```

The active job is fixed-geometry `RUNTYP=ENERGY`, not geometry optimization:

```bash
cd /mnt/e/ZFY/mcpb_5for_tleap_bidentate_M5_scf_20260623_140354/mcpb
tail -f gamess_logs/ZFY5B_M5_FE_scf_energy_soscf.log
cat gamess_logs/ZFY5B_M5_FE_scf_energy_soscf_to_ZFY5B_M5_FE_scf_energy_soscf_moread1.chain_status
tail -f gamess_logs/ZFY5B_M5_FE_scf_energy_soscf_moread1.log
```

This run uses `DIIS=.F.`, `SOSCF=.T.`, `DAMP=.T.`, and `SHIFT=.T.` for the
first fixed-geometry SCF. A watcher is running in tmux; if the first job writes
a `$VEC` restart but does not reach `DENSITY CONVERGED`, it automatically
launches `ZFY5B_M5_FE_scf_energy_soscf_moread1` with `GUESS=MOREAD`. Evaluate
success only from the log markers, not from shell exit code:

```text
success marker: DENSITY CONVERGED
failure markers: SCF IS UNCONVERGED / SCF HAS NOT CONVERGED
```

Result on 2026-06-23: this bidentate M5 fixed-geometry chain did not produce a
usable converged state. The first `ZFY5B_M5_FE_scf_energy_soscf` job terminated
normally at the shell/GAMESS level but printed:

```text
SCF IS UNCONVERGED, TOO MANY ITERATIONS
FINAL U-B3LYP ENERGY IS 0.0000000000 AFTER 200 ITERATIONS
S-SQUARED = 10.535
```

The automatic `GUESS=MOREAD` retry,
`ZFY5B_M5_FE_scf_energy_soscf_moread1`, initially looked encouraging on the
coarse DFT grid: the density change dropped to about `2.3e-4` by SCF iteration
37. However, when GAMESS switched back to the fine grid, the SCF became
nonphysical and entered a large two-cycle oscillation:

```text
DFT CODE IS SWITCHING BACK TO THE FINE GRID
...
energy jumps from about -2949 a.u. to +2695 / -333 a.u.
density change remains hundreds
```

This job was manually stopped and must not be used for `GRADIENT`, `OPTIMIZE`,
Hessian, MCPB force constants, or MD parameters:

```text
/mnt/e/ZFY/mcpb_5for_tleap_bidentate_M5_scf_20260623_140354/mcpb/
  gamess_logs/ZFY5B_M5_FE_scf_energy_soscf_moread1_manualstop_*.note
```

Important lesson: "the log looked good before the fine-grid switch" is not a
success. For these Fe-radical DFT jobs, only treat the SCF as usable after the
final fine-grid stage prints `DENSITY CONVERGED` and the log has no
`SCF IS UNCONVERGED` marker. Do not start a faster/more aggressive
optimization from a vector that diverged after grid restoration.

A later fine-grid-only MOREAD fixed-geometry run did succeed for the same ZFY
bidentate M5 project:

```text
Project:
  /mnt/e/ZFY/mcpb_5for_tleap_bidentate_M5_scf_20260623_140354/mcpb

Job:
  ZFY5B_M5_FE_scf_energy_soscf_finegrid_moread1

DFT grid:
  NRAD=96 NLEB=302 NRAD0=96 NLEB0=302 SWOFF=0.0

Result:
  DENSITY CONVERGED
  FINAL U-B3LYP ENERGY = -2949.2280370328
  58 SCF iterations
  normal GAMESS termination
```

This is a real fixed-geometry electronic state, but it still does not make
geometry optimization cheap. A bounded "fast OPT" retry was launched from the
converged fine-grid vector:

```text
Job:
  ZFY5B_M5_FE_opt_finegrid_moread1_nstep40_fast

Settings:
  RUNTYP=OPTIMIZE
  GUESS=MOREAD NORB=615
  NSTEP=40
  OPTTOL=0.0005
  DIIS=.F.
  SOSCF=.T.
  fine-grid DFT from the first iteration
```

The OPT did not crash immediately and initially behaved much better than the
coarse-grid attempt: the first geometry point slowly reduced the density change
from about `2e-3` to near the `2e-5` threshold. However it repeatedly slipped
back above the threshold, restarted SOSCF, and by SCF iteration 185 jumped to:

```text
DENSITY CHANGE = 0.050283398
ORB. GRAD      = 0.003465836
```

The run was manually stopped on 2026-06-23 after about 187 first-point SCF
iterations to avoid an unbounded optimization. Interpretation:

```text
fine-grid ENERGY success: yes
fast/cheap OPT from that vector: not yet practical
Hessian/MCPB force constants: do not run from this state
```

Do not read the manual `SIGINT`/`ddikick fatal` tail as the primary chemical
failure; it was an intentional stop after the first OPT SCF point became too
sticky. The actionable conclusion is that this ZFY bidentate Fe-radical model
needs a different optimization strategy before Hessian, such as a fixed-geometry
gradient diagnostic, a smaller first-shell QM model, explicit harmonic
Fe-donor constraints plus a cleaner SCF mode, or a chemistry/model revision.

If the neutral `UNK` radical bidentate model (`UNK` mol2 charge near 0,
MCPB charge `+1`, `MULT=5`) fails fixed-geometry GRADIENT/OPT, test the
chemically distinct anionic-UNK branch instead of only changing SCF knobs. The
diagnostic branch generated on 2026-06-23 is:

```text
Project:
  /mnt/e/ZFY/mcpb_5for_tleap_bidentate_UNKneg1_charge0_M6_scf_20260623_v1

Source:
  /mnt/e/ZFY/mcpb_5for_tleap_bidentate_M5_scf_20260623_140354

Change:
  UNK.mol2 total charge: -1
  MCPB small/large charge: 0
  multiplicity: 6
```

This is deliberately a fixed-geometry SCF smoke test, not a final ligand
parameterization. It keeps atom names, atom types, bonds, and coordinates
unchanged, then redistributes the extra -1 charge over `UNK` donor atoms
`N1,O2,C9` for a quick electronic-stability test. Because MCPB step 1 reports:

```text
small model: 55 atoms, 277 electrons
large model: 85 atoms, 445 electrons
```

`MULT=5` is not valid for this charge-zero model. Odd-electron models need an
even multiplicity, so the first test uses `MULT=6`.

Use the durable helper, not a root-level temporary script:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/prepare_mcpb_charge_variant_scf.py \
  --source-dir /mnt/e/ZFY/mcpb_5for_tleap_bidentate_M5_scf_20260623_140354 \
  --outdir /mnt/e/ZFY/mcpb_5for_tleap_bidentate_UNKneg1_charge0_M6_scf_20260623_v1 \
  --group ZFY5B_UNKneg1_M6_FE \
  --ligand UNK \
  --ligand-charge -1 \
  --charge-atoms N1,O2,C9 \
  --model-charge 0 \
  --mult 6 \
  --cores 4 \
  --run-step1

cd /mnt/e/ZFY/mcpb_5for_tleap_bidentate_UNKneg1_charge0_M6_scf_20260623_v1/mcpb
env NCORES=4 SESSION=zfy_unkneg1_m6_scf \
  bash run_ZFY5B_UNKneg1_M6_FE_scf_energy_soscf_finegrid_tmux.sh
tail -f gamess_logs/ZFY5B_UNKneg1_M6_FE_scf_energy_soscf_finegrid.log
```

Bug fixed in the helper: one-line GAMESS cards such as
`$STATPT NSTEP=1000 $END` must be removed as a single line. Otherwise the next
card, e.g. `$BASIS GBASIS=N31 ...`, can be accidentally deleted and GAMESS will
fail immediately with:

```text
ILLEGAL BASIS FUNCTION TYPE=C NGAUSS=6
```

That error is an input-generation bug, not evidence against the chemical model.

Do not assume AmberTools `antechamber -c bcc` can regenerate a radical-anion
ligand template. Testing `UNK=-1` with `-m 2` on this ZFY branch showed:

```text
antechamber reports: total electrons 87, net charge -1
sqm.in contains: spin=2 qmcharge=-1
sqm.out reports: QMMM spin multiplicity has legal range 1..1
```

Interpretation: this local SQM/AM1-BCC path is effectively closed-shell only
for this use case. It is fine for ordinary ACT(-1) and neutral closed-shell
ligands, but it cannot automatically produce a production-quality BCC mol2 for
the `UNK(-1)` doublet radical anion. For that branch, either use a deliberate
manual/diagnostic mol2 charge model, or fit charges through a QM/RESP workflow
outside the closed-shell AM1-BCC shortcut.

Important correction from the TJ audit on 2026-06-24:

The successful TJ MCPB handoff does not prove that the radical ligand mol2 must
be net `-1`. In `/mnt/e/TJ/260622_mcpb_m7_200ns/03_mcpb`, the actual values are:

```text
ACT.mol2 charge_sum = -1.000001
UNL.mol2 charge_sum = -0.000000
mcpb.in smmodel_chg/lgmodel_chg = +1
mcpb.in smmodel_spin/lgmodel_spin = 7
small model: 62 atoms, 288 electrons
large model: 92 atoms, 456 electrons
```

So the TJ pattern is: anionic `ACT`, neutral radical `UNL`, full MCPB QM cluster
`charge=+1, MULT=7`. The radical is represented by deleting one hydrogen from
the radical carbon, not by making `UNL.mol2` net `-1`. Do not automatically
convert a ZFY radical ligand into `UNK=-1` or the whole MCPB model into
`charge=0, MULT=6` just because the ligand is radical.

TJ recheck on 2026-06-25 for the HID-fixed restart
`/mnt/e/TJ/240624_tj_m7_hidfix_scf_from_initial/mcpb`:

```text
Fe-His:
  HID A187 NE2-bound neutral HID -> ok_ne2_bound_neutral_hid
  HID A270 NE2-bound neutral HID -> ok_ne2_bound_neutral_hid

mol2 charge sums:
  ACT.mol2 = -1.000001
  UNL.mol2 = -0.000000
  FE.mol2  = +3.000000
  HOH.mol2 =  0.000000

MCPB model:
  smmodel_chg/lgmodel_chg = +1
  small model = 296 electrons
  large model = 466 electrons
```

Therefore both `MULT=7` and `MULT=5` are parity-allowed for this even-electron
model. In that restart, the fresh fine-grid `MULT=7` SOSCF single point reached
normal GAMESS termination but remained `SCF IS UNCONVERGED` at `MAXIT=200`; it
is not usable for `MOREAD`, `GRADIENT`, `OPT`, or Hessian. The next fallback is
fresh fine-grid `MULT=5` SOSCF with the same charge/model, not a chemically
different UNL charge branch.

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

## 8.1 TJ MCPB Handoff Lessons: TER, Coordinating Water, And Frame Saving

The TJ Fe-radical MCPB handoff exposed two Amber/GROMACS pitfalls that should be
treated as standard checks before long production MD.

First, MCPB.py can write many internal `TER` records in `*_mcpbpy.pdb`. Amber
can interpret those internal breaks as real chain termini and silently add fake
terminal atoms such as `OXT`, `H1`, `H2`, or `H3` inside the protein chain. In
TJ this produced an unphysical overlap near `TRP 71 OXT` / `GLY 72 N`, and
GROMACS minimization failed with an infinite force. If a MCPB-generated PDB has
many internal protein `TER` records, clean it before `tleap`:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/clean_mcpb_internal_ter.py \
  --in-pdb TJ_M7_FE_mcpbpy.pdb \
  --out-pdb TJ_M7_FE_mcpbpy_no_internal_ter.pdb
```

Then make the `tleap` input load the cleaned PDB. In TJ, this reduced the
solvated `tleap` handoff to `Errors = 0` and only 9 warnings, after which
ACPYPE, EM, NVT, and NPT all passed.

Second, a MCPB coordinating water can appear as a residue such as `HH1` with
oxygen type `Y5` and hydrogen type `h`. AmberTools may have no local parameters
for `Y5-h` and `h-Y5-h`. Patch the MCPB `tleap` input with the PLE helper:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/patch_mcpb_tleap_water_h.py \
  --amber-dir 04_amber \
  --tleap-in leap_tj_mcpb_fixed.in \
  --tleap-out leap_tj_mcpb_no_internal_ter_hvdw.in \
  --loadpdb TJ_M7_FE_mcpbpy_no_internal_ter.pdb \
  --output-prefix TJ_M7_FE_noTER
```

This writes a small local frcmod with TIP3P-like `Y5-h` bond and `h-Y5-h`
angle terms plus zero Lennard-Jones for the water hydrogens. It is a practical
handoff patch for the coordinated water atom types generated by MCPB.

For long production trajectories, use literature-style frame intervals unless
there is a specific fast-motion question. A good default is:

```text
2 ps: high-resolution analysis or short/key-window runs; large files
5 ps: practical default for 100-400 ns production
10 ps: coarse screening when disk and analysis cost matter more than detail
```

At `dt = 0.002 ps`, a 5 ps compressed trajectory interval means:

```text
nstxout-compressed = 2500
nstenergy          = 2500
nstlog             = 2500
```

For 200 ns this gives about 40,000 frames. In the TJ MCPB production handoff,
GROMACS estimated about 9.4 GB output at this 5 ps setting. A 2 ps setting would
give about 100,000 frames and substantially larger trajectories.

The current TJ production directory after the no-TER repair is:

```text
/mnt/e/TJ/mcpb_m7_400ns/05_gromacs/TJ_M7_FE_noTER_solv.amb2gmx
```

Monitor the run without entering tmux:

```bash
cd /mnt/e/TJ/mcpb_m7_400ns/05_gromacs/TJ_M7_FE_noTER_solv.amb2gmx
tail -f logs/mdrun_production_200ns.log
```

Do not confuse the successful TJ handoff with the current ZFY work. On
2026-06-22 the ZFY `5_for_tleap.pdb` M5 optimization restarted from a converged
MOREAD vector, but the geometry optimization still failed with:

```text
FAILURE TO LOCATE STATIONARY POINT, SCF HAS NOT CONVERGED
EXECUTION OF GAMESS TERMINATED -ABNORMALLY-
```

That ZFY log is not valid for Hessian, MCPB force constants, or production MD.
The next ZFY step is chemical-model/preoptimization review, not Hessian.

Third, inspect MCPB custom residue bonds before trusting production MD. In the
TJ 200 ns MCPB run, PyMOL showed an odd-looking `HD1 173` residue. This was not
merely a PyMOL atom-color or stick-display bug: `HD1` is the MCPB/Amber custom
residue name for one Fe-coordinating histidine, and the GROMACS topology was
missing the `CA-CB` bond for that residue:

```text
resnr 173 HD1:
  N-CA    present
  CA-C    present
  C-O     present
  CA-CB   missing
  CB-CG   present
  NE2-FE  present

resnr 256 HD2:
  CA-CB   present
```

The missing `CA-CB` bond allowed the HD1 side chain to drift away from its
backbone during MD:

```text
HD1 173 CA-CB distance:
  initial GROMACS gro: 1.897 A
  EM:                 3.034 A
  NVT:                3.776 A
  NPT:                3.316 A
  raw production end: 6.329 A
```

Therefore this trajectory must be treated as compromised for conclusions that
depend on the local Fe/HIS187 environment. PBC-fixing or PyMOL style changes
cannot repair a missing covalent bond in the topology. Before rerunning TJ
production, patch the topology or regenerate the MCPB/Amber handoff so `HD1`
has a normal `CA-CB` bond/angle/dihedral connectivity like `HD2`, then rerun
EM/NVT/NPT and only then production. Add an automated pre-production check:

```text
For each custom histidine-like residue HD*/HID-derived:
  CA-CB topology bond must exist
  CA-CB distance after EM must stay near 1.5 A
  Fe-donor bonded distances must remain near the MCPB target
```

Root cause tracing on 2026-06-23 showed the failure happened before `tleap` and
before ACPYPE/GROMACS conversion:

```text
Original/H++ His187 CA-CB:                  1.534 A
M7 metal-center patch full complex CA-CB:   1.892 A
mcpb_original.pdb His187 CA-CB:             1.892 A
MCPB-generated HD1.mol2 CA-CB:              1.892 A, bond missing

Original/H++ His270 CA-CB:                  1.538 A
M7 metal-center patch full complex CA-CB:   1.612 A
MCPB-generated HD2.mol2 CA-CB:              1.612 A, bond present
```

`HD1.mol2` itself lacked the `CA-CB` bond, while `HD2.mol2` contained it.
Therefore `tleap` did not report an error because it loaded a syntactically
valid user-supplied nonstandard residue template; a missing internal covalent
bond in that template is chemically wrong, but not a `tleap` syntax failure.
The likely upstream cause was the "patch optimized small-model coordinates back
onto the full protein" step: it moved His187 side-chain atoms enough to stretch
the CA-CB distance beyond MCPB/Amber mol2 bond perception, so MCPB generated an
internally disconnected custom histidine residue.

Future rule: when patching QM-optimized metal-center coordinates back into a
full protein, do not blindly overwrite backbone-adjacent atoms such as `CA` or
`CB` unless the whole residue is moved by a rigid transform that preserves
internal residue geometry. For side-chain metal ligands, prefer patching only
the metal-donor side-chain atoms required by the coordination model, then run an
explicit covalent-connectivity audit before MCPB step 4 and again after ACPYPE:

```text
1. Compare original vs patched CA-CB for every coordinating HIS/GLU/ASP.
2. Reject any custom amino-acid mol2 missing backbone/side-chain bonds:
   HIS-like: N-CA, CA-C, C-O, CA-CB, CB-CG, ring bonds.
   GLU/ASP-like: CA-CB, CB-CG, CG-CD/OE bonds as appropriate.
3. Reject any post-tleap topology where those bonds are absent.
4. Do not rely on `tleap` warnings alone for this class of error.
```

The practical TJ repair on 2026-06-23 was:

```text
Keep:
  TJ_M7_FE_mcpbpy.frcmod
  MCPB Fe-donor force constants from the GAMESS Hessian/Seminario workflow

Reject as direct MD coordinates:
  TJ_M7_FE_mcpbpy.pdb with HD1 A187 CA-CB = 1.892 A

Repair:
  1. Patch HD1.mol2 so the CA-CB bond exists.
  2. Generate TJ_M7_FE_mcpbpy_coordfix.pdb by preserving the M7 side-chain
     direction but shortening HD1 CA-CB back to the reference 1.534 A.
  3. Remove internal protein TER records:
     TJ_M7_FE_mcpbpy_coordfix_no_internal_ter.pdb
  4. Rebuild tleap with source leaprc.gaff2 plus the local HH1 water-hydrogen
     patch.
  5. Replace old MCPB bond commands such as mol.187.NE2/mol.431.FE with the
     actual tleap-loaded unit indices from desc mol:
       AT1=501, HD1=673, HD2=756, GU1=835, FE1=864, HH1=1004, UL1=1005.
```

Use the packaged helper instead of rewriting temporary scripts:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/tleap_index_tools.py probe \
  --base-leap leap_tj_mcpb_fixed.in \
  --out leap_probe.in \
  --targets 187 270 349 431 501 875 \
  --atoms FE N1 O O2 NE2 OE1

python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/tleap_index_tools.py scan \
  --base-leap leap_tj_mcpb_fixed.in \
  --out leap_scan.in \
  --start 1 \
  --end 1005

python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/tleap_index_tools.py rewrite-bonds \
  --pdb TJ_M7_FE_mcpbpy_coordfix_no_internal_ter.pdb \
  --leap-in leap_tj_mcpb_fixed.in \
  --out leap_tj_mcpb_fixed_unit_indices.in \
  --chain A
```

The `rewrite-bonds` mode only rewrites simple `bond mol.<number>.<atom>`
selectors where the number is uniquely mapped from the loaded PDB residue order.
If it prints `UNRESOLVED`, inspect with `probe`/`scan` and patch manually.

The corrected handoff passed:

```text
MCPB mol2 audit:
  HD1 CA-CB present

MCPB PDB audit:
  HD1 A187 CA-CB = 1.534 A

GROMACS initial topology/gro audit:
  HD1 173 CA-CB topology bond = true
  HD1 173 CA-CB = 1.534 A

After NPT:
  HD1 173 CA-CB topology bond = true
  HD1 173 CA-CB = 1.590 A
  HD2 and GU1 also passed the same backbone/side-chain checks.
```

This is the model to follow: Hessian-derived force constants are reusable, but
bad full-protein patched coordinates are not sacred. Treat the parameters and
coordinates as two different products of the workflow.

## 9. Fe-Coordinating Histidine Protonation Audit

Fe-bound histidines must be audited as chemical residues, not only as atoms that
appear near the metal. In an ordinary neutral histidine that coordinates Fe
through `NE2`, the expected Amber-like state is usually `HID`: `NE2` donates to
Fe and has no H, while the opposite imidazole nitrogen `ND1` carries one H.
If coordination is through `ND1`, the analogous neutral state is usually `HIE`,
with `NE2-H`.

Do not silently accept a custom MCPB histidine residue where both imidazole
nitrogens lack H unless the intended chemistry is explicitly an imidazolate-like
deprotonated ligand. This changes the local charge model, Fe-His bond strength,
hydrogen-bond network, and mechanistic interpretation.

TJ audit on 2026-06-24 found a second serious template issue after the `CA-CB`
repair was complete:

```text
03_mcpb/HD1.mol2 atoms:
  N CA C O CB CG CD2 ND1 CE1 NE2 HA
  no ND1-H, no NE2-H, no ring C-H atoms

03_mcpb/HD2.mol2 atoms:
  N CA C O CB CG CD2 ND1 CE1 NE2 HA
  no ND1-H, no NE2-H, no ring C-H atoms

05_gromacs/.../production_200ns_centered_protein_compact.gro:
  173HD1 contains only heavy atoms plus HA
  256HD2 contains only heavy atoms plus HA
```

Therefore the completed TJ 200 ns production can still be postprocessed as the
"current MCPB-bound Fe control" trajectory, but its Fe-His chemistry is not a
standard neutral His model. Any RDC, distance, or pose interpretation from this
run must carry this caveat.

Future pre-production audit:

```text
For each Fe-coordinating custom histidine-like residue HD*/HX*:
  1. CA-CB and ring covalent bonds must exist.
  2. If the intended state is neutral HID/HIE, exactly one imidazole N-H must exist.
     NE2-bound Fe-His: ND1-H present, NE2-H absent.
     ND1-bound Fe-His: NE2-H present, ND1-H absent.
  3. Ring carbon hydrogens should not disappear unless the MCPB residue template
     and charges intentionally use a special united/capped model.
  4. Compare atom count and N-H placement against the H++/Amber input residue.
  5. Reject the handoff before production if the histidine protonation state is
     not chemically intentional.
```

Use the packaged audit helper before accepting a MCPB handoff:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/audit_mcpb_histidine_protonation.py \
  --pdb /path/to/mcpb/ZFY5B_M7_FE_small.pdb \
  --cutoff 2.8

python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/audit_mcpb_histidine_protonation.py \
  --mol2 /path/to/03_mcpb/HD1.mol2 /path/to/03_mcpb/HD2.mol2 \
  --mol2-donor NE2
```

The ZFY bidentate M7 branch on 2026-06-24 showed the exact failure mode this
helper is meant to catch:

```text
HID A161: Fe-NE2 = 2.171 A, Fe-ND1 = 4.070 A, expected neutral state = HID,
          but no ND1-H or NE2-H was present.
HID A241: Fe-NE2 = 2.042 A, Fe-ND1 = 4.068 A, expected neutral state = HID,
          but no ND1-H or NE2-H was present.
```

The same day, `RUNTYP=GRADIENT` for that branch failed before any usable
gradient was produced:

```text
SCF IS UNCONVERGED, TOO MANY ITERATIONS
NO GRADIENT, SCF DID NOT CONVERGE
```

Interpretation: do not proceed to OPT/Hessian by calling this a gradient
failure. Rebuild or patch the Fe-His protonation/electronic model first, then
rerun fixed-geometry SCF and only then GRADIENT/OPT.

The MCPB project generators now apply this rule before step 1 by default:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/fix_fe_histidine_protonation.py \
  --pdb hand_built_or_hpp_complex.pdb \
  --out mcpb_original_fe_his_fixed.pdb \
  --cutoff 2.8
```

Default rule:

```text
Fe-NE2 histidine -> HID, add ND1-H, remove any NE2-H.
Fe-ND1 histidine -> HIE, add NE2-H, remove any ND1-H.
```

This is now called automatically by
`prepare_fe_radical_mcpb_opt_project.py` unless `--no-fix-fe-his` is passed, and
by the TJ canonical project generator before writing `03_mcpb/mcpb_original.pdb`.
After MCPB step 4, the generated `HD*/HE*` mol2 files must still be audited,
because residue renaming alone is not enough. The mol2 template must contain the
actual N-H atom and bond, not only the residue name `HID`.

## 10. TJ Directory Cleanup Rule

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

## 11. Reusing A Previous MCPB PDB As A New Input

Do not feed a previous MCPB/Amber custom-residue PDB directly into a fresh MCPB
run. Names such as `HD1`, `HD2`, `GU1`, `UL1`, `AT1`, `HH1`, and `FE1` are
downstream custom names, not a clean human-readable starting model. Normalize
them first:

```text
HD1/HD2 -> HID
GU1     -> GLU
UL1     -> UNL
AT1     -> ACT
HH1     -> HOH
FE1     -> FE
```

For TJ-style rechecks, use:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/standardize_tj_mcpb_pdb.py \
  --pdb /mnt/e/TJ/2.pdb \
  --out /mnt/e/TJ/<project>/input/2_standard_names.pdb

python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/fix_fe_histidine_protonation.py \
  --pdb /mnt/e/TJ/<project>/input/2_standard_names.pdb \
  --out /mnt/e/TJ/<project>/input/2_standard_names_HIDfixed.pdb \
  --cutoff 2.8
```

The expected Fe-bound neutral histidine state is still:

```text
Fe-NE2 donor: HID, ND1-H present, NE2-H absent.
```

In the 2026-06-26 TJ `2.pdb` recheck, the clean donor list after normalization
and HID repair was:

```text
UNL A1 N1      1.910 A
GLU A349 OE1   1.945 A
ACT A501 O2    2.025 A
HOH A875 O     2.027 A
HID A270 NE2   2.168 A
HID A187 NE2   2.477 A
```

MCPB step 1 then produced a small model with 72 atoms and 298 electrons for
`charge=+1, MULT=7`.

## 12. Nitrene-like UNL/UNK Templates: Mol2 Is Topology, Not Spin Truth

For hand-built Fe-radical or Fe-nitrene intermediates, do not expect a mol2
file to faithfully encode where the unpaired electron lives. In practice, mol2
is used by Amber/MCPB as a residue topology and partial-charge template: atom
names, atom types, bond graph, and charges. The open-shell electronic state is
defined in the QM job through `ICHARG`, `MULT`, and the UHF/DFT wavefunction.

If the user moves the radical character from benzylic carbon to a metal-bound
nitrogen and removes the N-H bond, the ligand mol2 must be made consistent with
the curated PDB before MCPB.py step 1. A common failure mode is:

```text
PDB:  UNL has no HN1; a hydrogen named H02 is bonded/placed near C4
mol2: UNL still contains HN1 and a N1-HN1 bond
```

This is not a harmless visualization difference. MCPB.py/tleap will treat the
mol2 graph as the ligand template, so the final model silently contradicts the
intended Fe=N/nitrene-like chemistry. Patch the mol2 graph first, then audit:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/patch_ligand_mol2_atom_from_pdb.py \
  --mol2 old/UNL.mol2 \
  --pdb /mnt/e/TJ/2.pdb \
  --out input/UNL.mol2 \
  --resname UNL \
  --rename HN1 H02 \
  --replace-bond C4 H02 \
  --new-type hc \
  --set-type C4=c3 \
  --set-type H4B=hc

python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/audit_ligand_mol2_bonds_in_pdb.py \
  --mol2 input/UNL.mol2 \
  --pdb /mnt/e/TJ/2.pdb \
  --resname UNL \
  --out logs/UNL_pdb_vs_patched_mol2.tsv
```

After patching, rerun MCPB.py step 1 and check electron parity. In the TJ
2026-06-26 nitrene-like retry, the patched neutral `UNL.mol2` retained 23 atoms
and total charge near zero. MCPB step 1 reported:

```text
small model: 72 atoms, 298 electrons
large model: 104 atoms, 468 electrons
```

Thus `charge=+1, MULT=7` and `charge=+1, MULT=5` remain parity-consistent. If
an H atom is truly deleted instead of moved, the electron count changes; do not
reuse the old charge/multiplicity blindly.

For this local GAMESS build, `$CONTRL MAXIT` must be between 0 and 200. Setting
`MAXIT=220` causes a parse-time abnormal termination before any SCF iteration:

```text
ERROR: MAXIT MUST BE BETWEEN 0 AND 200
```

Treat that as a script/input bug, not electronic-model failure. Use
`MAXIT=200` for long SCF attempts and judge the run by `DENSITY CONVERGED`,
`SCF IS UNCONVERGED`, and the SCF iteration trend.

## 12. SCF Monitoring Should Be Scripted

Avoid judging GAMESS SCF health by raw `tail` alone. The column labels scroll
away, and shell quoting regularly hides the useful lines. Use:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/summarize_gamess_scf.py \
  gamess_logs/<job>.log --last 30
```

For the same 2026-06-26 TJ recheck:

- `M7 fresh fine-grid SOSCF` was stopped early because orbital gradient and
  damping worsened quickly.
- `M7 fresh fine-grid DIIS` reached density values near `0.005` but DIIS error
  stalled around `0.3`, so it was stopped as a poor trend.
- `M5 fresh fine-grid DIIS` temporarily approached low density changes
  (`1e-5` to `1e-4`) but then jumped to an electronic-collapse state. Its log
  was preserved as
  `TJ2_HID_M5_FE_scf_energy_fresh_finegrid_diis_stopped_after_electronic_collapse.log`.

This means HID protonation repair was necessary but not sufficient. When both
M7 and M5 fresh guesses fail or collapse, do not brute-force the same input.
Next options are: preserve/reuse a near-converged `$VEC` if the punch file is
complete, tighten SCF with a smaller-step strategy, or re-check the input
geometry/electronic model before starting OPT/Hessian.
