---
name: pocket-ligand-explorer
description: Combine SMILES-to-Vina docking and pose analysis into a no-official-score ligand exploration workflow. Use when the user wants a project-directory CLI that auto-detects a unique protein PDB/PDBQT, accepts manual SMILES or ligand files, creates an initial large Vina box, docks and mutates ligands, clusters pose landing regions, shrinks the box from pose evidence, continues analog exploration, and writes pose-analysis HTML/PyMOL outputs without relying on official binding scores.
---

# Pocket Ligand Explorer

Use this skill for arbitrary protein + arbitrary ligand exploration when there is
no external activity table and no official binding score. The workflow is:

1. Treat one directory as the project root.
2. Auto-detect the unique `*.pdb` or `*.pdbqt` receptor in that directory.
3. Accept starting ligands from `--smiles`, `--smiles-file`, `--sdf`, `--mol`, or
   a CSV column.
4. Create a large whole-protein Vina config if no box is supplied.
5. Dock the starting ligands and a small conservative analog batch.
6. Read pose coordinates, find the common pose landing region, and write a
   focused Vina config.
7. Run another analog/docking round in the focused box.
8. Run `pose-analyzer` with `--score-column affinity_kcal_mol` and build an
   offline HTML report with PyMOL `.pml` links and optional `.pse` sessions.
9. Optionally build the pose-analyzer molecule-space dashboard from the whole
   docking history for clustering, dimensionality reduction, affinity/CNN-pose
   triage, and contact inspection.

Do not introduce any external or official score column into this skill. Use
Vina/GNINA affinity only as an internal triage signal and always describe it as
docking-based, not experimental truth.

## Main Command

The bundled package lives in `package/`. After local editable install:

```bash
python -m pip install -e /path/to/pocket-ligand-explorer/package
```

The CLI reuses the sibling skills `smiles-to-vina-docking` and `pose-analyzer`.
When installed outside the shared `skills/` folder, set one of:

```bash
ple run --skills-dir /path/to/skills --smiles "CCO"
```

or:

```bash
export CODEX_SKILLS_DIR=/path/to/skills
```

Recommended environment:

```bash
conda create -n ple -c conda-forge python=3.11 rdkit meeko openbabel numpy scikit-learn pymol-open-source vina
conda activate ple
python -m pip install -e /path/to/pocket-ligand-explorer/package
```

The package also declares the PyPI `vina` dependency, so `pip install` can provide
Vina on systems where that wheel works. The package exposes `ple-vina`, a small
Vina-compatible CLI wrapper around the PyPI Vina Python API. The CLI still
auto-searches existing Vina or GNINA executables first.

Default docking-engine discovery checks:

- `PATH`
- `PLE_VINA`, `VINA_EXE`, `VINA`, `GNINA`, `GLINA`
- project root, `project/vina_bin`, `project/bin`
- `ple_project/vina_bin`
- active conda `bin/` or `Scripts/`
- `~/vina_task2/vina_bin`
- bundled `ple-vina` wrapper from the PyPI `vina` package

When running inside a conda environment named `gnina` or `glina`, `ple` prefers
`gnina`/`glina` from that environment or PATH before falling back to the older
`~/vina_task2/vina_bin/vina` location. You can also force this behavior with
`PLE_ENGINE=gnina`.

If Vina/GNINA is installed elsewhere, use `--engine-exe /path/to/vina`.
For GNINA CNN refinement-style scoring, pass:

```bash
ple --smiles "CCO" --engine-exe gnina --gnina-cnn-scoring refinement
```

When GNINA/GLINA is used with `--gnina-cnn-scoring refinement`, PLE defaults to a
two-stage strategy: the initial whole-protein scout uses cheaper
`--cnn_scoring rescore` and a small wide batch, then the pose-derived focused box
uses `--cnn_scoring refinement`. This avoids spending most of the run inside a
large 80 A search box.

If the receptor is a raw PDB, `ple run` prepares
`ple_project/receptor/<name>.pdbqt` automatically with OpenBabel by default.

Users can run from inside a project directory:

```bash
ple run --smiles "CCO"
```

Short form also works:

```bash
ple --smiles "CCO"
```

or explicitly:

```bash
ple run --project-dir /path/to/project --smiles "CCO"
```

The default exploration strength is `--degree 1`, intended to keep the current
small-run scale: wide-box docking, focused-box docking, pose analysis, HTML
report, and PyMOL `.pml` views. Increase the degree instead of memorizing many
low-level flags:

```bash
ple --smiles "CCO" --degree 2
ple --smiles "CCO" --degree 3
```

Degree presets:

- `--degree 1`: about 50 attempted rows, conservative edits, pose report and
  PyMOL `.pml` outputs.
- `--degree 2`: about 100 attempted rows, still conservative.
- `--degree 3`: adds drastic edits and selects a few good existing history
  parents using affinity plus `inner_rmsd` and `cnn_pose_score` triage.
- `--degree 4`: broader degree 3-style exploration, up to about 200 extra
  attempts depending on duplicates and property filters.
- `--degree 5`: strongest preset; more history-guided parents and drastic
  exploration.

All degree presets include `replace` edits by default. This keeps exploration
from becoming add-only: analog generation should interleave additions with
amide/linker bioisosteres, aromatic C-to-N swaps, and phenyl-to-small
heteroaryl replacements when chemically possible.

Manual overrides still work, for example:

```bash
ple --smiles "CCO" --degree 3 --focus-batch-size 80
```

For metalloenzymes or crystal structures that contain an old occupant ligand,
define the docking site explicitly instead of starting from the whole-protein box.
Use `--drop-receptor-resname` to remove the blocking ligand from the receptor copy,
and `--site-resname`/`--site-chain` to center the first box on the catalytic atom or
residue. Example for a 2R5V-derived receptor where `HHH` occupies the Fe pocket:

```bash
ple --smiles "O=C(NOC(C)=O)CCCC1=CC=CC=C1" \
  --drop-receptor-resname HHH \
  --site-resname FE \
  --site-chain A \
  --site-box-size 21 \
  --engine-exe gnina \
  --gnina-cnn-scoring refinement \
  --no-gnina-wide-scout \
  --receptor-prep none
```

The generated receptor copy is written under `ple_project/receptor/`, the config is
centered on the matched site atom/residue, and the original PDB is left untouched.
If multiple residues match the site selector, add `--site-resseq` or
`--site-atom-name`.

Automatic focused boxes are inferred from the common landing region of the best
successful poses, not from every scattered pose in the whole-protein scout. The
default focused box target is `21 x 21 x 21 A`, with a maximum automatic side
length of `21 A`. Rebuild only the focused config from an existing history:

```bash
ple focus-box --state-dir ple_project
```

Build the global molecule-space dashboard when you want the advanced
pose-analyzer visualization:

```bash
ple --smiles "CCO" --degree 3 --build-space-dashboard
```

This writes `ple_project/score_space_dashboard/index.html`. In PLE mode this is
an unsupervised clustering/dimensionality dashboard, not a supervised binding
score model. Use `affinity_kcal_mol`, `cnn_pose_score`, `inner_rmsd`, contacts,
and descriptors together for triage.

Rebuild analysis outputs from an existing `dock_history.csv` without docking:

```bash
ple --analysis-only --build-space-dashboard
```

Use this after a long run when you only want to regenerate the pose report,
PyMOL files, or global dashboard.

Prepare a selected docking result for the later MD workflow:

```bash
ple md-handoff --state-dir ple_project --pdb target.pdb
```

This creates `ple_project/md_handoff/<seq_id>/` with copied receptor/ligand
inputs, AmberTools `tleap` input, ACPYPE conversion scripts, GROMACS minim/NVT/NPT
and production `.mdp` files, and a small RMSD/Rg FEL plotting script. It does
not run the long production MD by itself. Run the generated scripts in order
inside the conda MD environment:

```bash
cd ple_project/md_handoff/<seq_id>
./00_check_env.sh
./01_prepare_amber_to_gmx.sh
./02_run_gromacs.sh
./03_analyze_basic.sh
```

For real project work, prefer the higher-level flat-layout command:

```bash
ple md-run \
  --project-dir /mnt/e/TJ \
  --state-dir ple_project_product_R \
  --pdb receptor_hpp_curated_keepFe.pdb \
  --seq-id S000001 \
  --pose-mode 6 \
  --protonation-map protonation_map.csv \
  --receptor-md-mode full \
  --keep-hetatm-element Fe \
  --ambertools-bin /mnt/l/WSL/conda_envs/AmberTools25/bin \
  --gmx /home/qin/softwares/gromacs-2026.0-gpu/bin/gmx_mpi \
  --gmxrc /home/qin/softwares/gromacs-2026.0-gpu/bin/GMXRC \
  --time 50 \
  --mdrun-arg=-nb --mdrun-arg=gpu --mdrun-arg=-pme --mdrun-arg=gpu
```

`ple md-run` writes one coherent handoff in `<project-dir>/md_handoff/` instead
of nesting under `ple_project_product_R/md_handoff/S000001_mode6_...`. It then
runs `00_check_env.sh`, `01_prepare_amber_to_gmx.sh`, `02_run_gromacs.sh`, and,
for foreground runs, the generated `03_postprocess.sh`. Background runs chain
`02_run_gromacs.sh` and `03_postprocess.sh` through
`04_run_gromacs_then_postprocess.sh` unless `--skip-postprocess` is set. The
postprocess script calls `ple md-fix-pbc`, `ple pymol-metal`, and `ple md-rdc`,
producing:

- `<project-dir>/md_handoff/gromacs/md_nojump.xtc`
- `<project-dir>/md_handoff/gromacs/md_centered_compact.xtc`
- `<project-dir>/md_handoff/gromacs/md_centered_compact.gro`
- `<project-dir>/pymol/show_metal_centered.pml`
- `<project-dir>/md_analysis/rdc/index.html`

Use `--prepare-only` to write the full commandable workflow without running MD.
Use `--background` for long production jobs; this writes `md_handoff/logs/md.pid`
and logs the chained MD/PBC/RDC/PyMOL run to
`md_handoff/logs/run_gromacs_background.log`.

Use `--pose-mode` when the pose-analysis/PyMOL report points to a specific
GNINA/Vina model such as `seqS000003_mode8`:

```bash
ple md-handoff \
  --state-dir ple_project \
  --pdb target.pdb \
  --seq-id S000003 \
  --pose-mode 8
```

For local GROMACS builds that need an initialization script or extra shared
library path, write those into the handoff once. If AmberTools/ACPYPE lives in a
separate environment, point to that `bin/` directory with `--ambertools-bin`; the
generated scripts keep the current MD Python first and append the AmberTools
commands to `PATH`.

```bash
ple md-handoff \
  --state-dir ple_project \
  --pdb target.pdb \
  --ambertools-bin /mnt/l/WSL/conda_envs/AmberTools25/bin \
  --gmx /home/qin/softwares/gromacs-2026.0-gpu/bin/gmx_mpi \
  --gmxrc /home/qin/softwares/gromacs-2026.0-gpu/bin/GMXRC \
  --time 10 \
  --extra-ld-library-path /home/qin/softwares/cp2k-2026.1/build_minlib/src \
  --mdrun-arg=-nb --mdrun-arg=gpu --mdrun-arg=-pme --mdrun-arg=gpu
```

`--time 10` means a conventional 10 ns production MD stage; `--time 20` means
20 ns. It is an alias for the older `--production-ns` option and controls the
generated production `.mdp` file.

Use repeated `--mdrun-arg` for tokens starting with `-`; it avoids shell/argparse
quoting problems.

The default MD receptor mode is `--receptor-md-mode protein-only`. It preserves
the original PDB as `input/receptor.pdb`, writes the actual MD input as
`input/receptor_md.pdb`, and records removed waters/cofactors/phosphate in
`input/receptor_md_report.json`. Fe/FE records are kept by default because a
catalytic iron center is usually part of the active system, not disposable
solvent. Use `--drop-default-md-metals` only for a temporary protein-only smoke
test.

Manual H++ or literature-curated protonation can be applied directly to the
Amber input PDB:

```bash
ple md-handoff \
  --state-dir ple_project \
  --pdb target.pdb \
  --protonate A:123=HIE \
  --protonate A:45=ASH
```

or from a text/CSV file:

```bash
ple md-handoff --protonation-map hpp_amber_names.csv
```

Accepted one-line rules include `A:123=HID`, `A:45=ASH`, and `128=GLH`. CSV maps
can use columns such as `chain,resseq,resname`. Use Amber residue names directly:
`ASH`, `GLH`, `HID`, `HIE`, `HIP`, `CYM`, `CYX`, etc.

If H++ produced a `.pkout` or copied `.pkout.txt`, convert it into a reusable
Amber protonation map first:

```bash
ple hpp-map \
  --project-dir /path/to/project \
  --pkout hplusplus.pkout.txt \
  --pdb prepared/receptor_no_old_ligand_keepFe.pdb \
  --out protonation_map.csv \
  --chain A \
  --ph 7.4
```

This writes `protonation_map.csv` plus
`protonation_map_metal_contacts.csv`. Explicit H++ names such as `HID`, `HIE`,
and `HIP` are preserved; remaining `HIS` entries are converted with conservative
pH-based heuristics and tagged in the `reason` column. Fe contacts within the
metal cutoff are audited so catalytic histidines can be reviewed before MD.

Extra metal/cofactor parameters can be copied into the handoff and loaded by
`tleap`:

```bash
ple md-handoff \
  --state-dir ple_project \
  --pdb target.pdb \
  --keep-hetatm-resname FE \
  --amber-frcmod iron_center.frcmod \
  --amber-lib iron_center.lib \
  --amber-prep iron_center.prepi
```

The generated `leap.in` loads `frcmod.ions234lm_126_tip3p` by default, matching
the earlier working AmberTools notebook style. This helps ordinary ions, but it
is not a complete force-field model for a catalytic iron center with a defined
oxidation state, spin state, ligand field, and coordinating residues. For serious
metalloenzyme MD, feed in literature parameters or MCPB.py-derived
`frcmod/lib/prep/mol2` files, then rerun `./01_prepare_amber_to_gmx.sh`.

Prepare a first MCPB.py metal-center draft after a successful MD handoff:

```bash
ple metal-model \
  --pdb ple_project/md_handoff/S000003/input/receptor_md.pdb \
  --ligand-pdb ple_project/md_handoff/S000003/amber/ligand_pose.pdb \
  --ligand-mol2 ple_project/md_handoff/S000003/amber/ligand.mol2 \
  --outdir ple_project/metal_model/S000003_FeIII \
  --metal-element FE \
  --metal-chain A \
  --metal-resseq 4113 \
  --ligand-role nonbonded \
  --smmodel-chg 2 \
  --smmodel-spin 6 \
  --lgmodel-chg 2 \
  --lgmodel-spin 6 \
  --software-version g16 \
  --cutoff 2.8 \
  --ambertools-bin /mnt/l/WSL/conda_envs/AmberTools25/bin
```

`ple metal-model` defaults to `--metal-charge 3.0`, i.e. Fe(III), and writes
`FE.mol2`, `mcpb_original.pdb`, `metal_site.pdb`,
`metal_coordination_candidates.csv`, `mcpb.in`, `run_mcpb_step1.sh`, and
`run_mcpb_workflow_skeleton.sh`. For ferric non-heme iron, use an explicit
high-spin assumption such as `--smmodel-spin 6 --lgmodel-spin 6`; do not leave
MCPB.py to auto-pick singlet spin for a Fe(III) catalytic center.

By default, docked ligands are treated as `--ligand-role nonbonded`: PLE still
reports nearby ligand N/O/S atoms in `metal_coordination_candidates.csv`, but it
does not include the docked ligand in the MCPB bonded metal-site fitting. This
keeps exploratory substrates free during MD. Only pass
`--ligand-role coordinating` when the ligand is intentionally part of the metal
coordination sphere; in that case also pass its antechamber `ligand.mol2` so PLE
can copy it as `LIG.mol2` and rewrite ligand PDB atom names to match the mol2
atom names.

For the current 2R5V-derived S000003 mode 8 Fe(III) protein-only draft, MCPB.py
step 1 successfully identified `HID A161 NE2`, `HID A241 NE2`, and
`GLU A320 OE1` within 2.8 A of `FE A4113`, wrote high-spin `2 6` Gaussian input
files, and left the docked ligand as nonbonded. The remaining
production-quality work is MCPB.py step 2/3 with completed QM output files, then
step 4 to generate bonded-model files for the final Amber/GROMACS handoff.

For Fe-bound radical intermediates where the substrate/product fragments are
part of the coordination sphere, read the skill-local
`references/FE_RADICAL_MCPB_ZH.md` before preparing MCPB input. The short rule
is: first fix ligand names and radical atom counts, then run a ligand-only
AmberTools preflight, then build a six-coordinate Fe model and screen plausible
spin multiplicities before the expensive Hessian. For Fe(III) high-spin plus an
organic radical, test both `MULT=5` and `MULT=7` unless the literature or a
validated electronic-structure result already fixes the coupling. Do not treat a
free ligand B3LYP optimization as the final active-site geometry; it is only a
geometry hygiene check outside the Fe/protein ligand field.

For metalloenzyme production handoff, prefer a curated receptor PDB, not the raw
docking receptor. The preferred order is: repair pocket-relevant missing loops,
run or manually curate H++ protonation, remove the old crystal occupant ligand,
then add back the catalytic metal/cofactor records that belong to the simulated
system. For the 2R5V-derived Fe/HHH case, the practical first-pass receptor is a
chain-A repaired H++ heavy-only PDB with `GLU A108` present, `HHH` removed, chain B
removed, and chain A `FE A4113` re-added before MD handoff:

```bash
ple md-handoff \
  --state-dir ple_project_AFe_noHHH_degree1 \
  --pdb /mnt/e/Pose_create/pdb2r5v_chainA_GLU108_Hpp_heavy_plus_AFe_noHHH.pdb \
  --seq-id S000003 \
  --pose-mode 8 \
  --ambertools-bin /mnt/l/WSL/conda_envs/AmberTools25/bin \
  --protonation-map /mnt/e/Pose_create/pdb2r5v_Hpp_pH7p4_first_pass_amber_protonation_map.csv \
  --keep-hetatm-element Fe \
  --time 10
```

When keeping metals as HETATM records, the MD receptor writer inserts `TER`
between protein `ATOM` records and metal/cofactor `HETATM` records. This prevents
Amber from treating the metal as part of the covalent protein chain and helps C
terminal residues such as `OXT` map to the correct terminal residue template.

The ligand preparation step prefers the selected docked PDBQT pose coordinates
while using the original SDF chemistry/connectivity. GNINA/Vina/OpenBabel PDBQT
atom order can differ from RDKit SDF atom order. Therefore the handoff copies the
prepared ligand PDBQT (`input/ligand_prepared.pdbqt`) and maps its serials back to
the SDF heavy atoms by matching the original coordinates and elements; `REMARK
SMILES IDX` is used as a fallback when available. The script writes
`amber/ligand_pose_mapping.tsv` and prints the maximum heavy-atom bond length. A
normal mapped ligand should not show multi-Angstrom covalent bonds such as 4-12 A.

For ferric non-heme exploratory MD before MCPB.py is available, a pragmatic
positive-control handoff is:

1. Normalize single-atom iron records to Amber ion names before `tleap`: `FE/FE`
   means Fe(III), while `FE2/FE2` means Fe(II). Crystal residue names such as
   `FE2` can mean "iron component" in the PDB rather than a desired Fe(II) MD
   oxidation state, so choose explicitly.
2. Keep H++-curated catalytic histidines such as `HID`/`HIE` from
   `protonation_map.csv`; audit `protonation_map_metal_contacts.csv`.
3. After ACPYPE writes `complex_GMX.gro/top`, add three strong harmonic Fe-core
   bonds to the known Fe ligands, e.g. two histidine `NE2` atoms and one
   glutamate/aspartate oxygen. This is not MCPB.py and must be reported as a
   restraint/positive-control model, but it prevents the known Fe-His/Glu core
   from drifting like a naked ion.

The old Codex helper scripts used for the TJ Fe(III) handoff are transitional
prototypes. Prefer PLE commands first: `ple md-run`, `ple md-fix-pbc`,
`ple md-rdc`, and `ple pymol-metal`. Keep one-off project scripts only as
debugging artifacts or archive material.

Before launching long MCPB-derived production MD, audit both the MCPB custom
residue templates and the GROMACS topology. `tleap` can accept a custom residue
mol2 that is syntactically valid but chemically broken, such as a histidine-like
`HD1` residue missing its `CA-CB` bond. Use:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/audit_mcpb_custom_residue_bonds.py \
  --mcpb-dir /path/to/03_mcpb \
  --pdb /path/to/03_mcpb/TJ_M7_FE_mcpbpy_coordfix_no_internal_ter.pdb

python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/audit_gromacs_custom_residue_bonds.py \
  --gro /path/to/gromacs/npt.gro \
  --top /path/to/gromacs/TJ_M7_FE_coordfix_noTER_solv_GMX.top
```

If a QM small-model patch distorted one full-protein residue but the Hessian and
MCPB force constants are still valid, keep the `frcmod`/force constants and fix
the bad full-protein coordinates separately. For the TJ HD1 problem, the
accepted repair was to patch `HD1.mol2`, clamp only `A:187` `CA-CB` back to the
reference length, clean internal MCPB `TER` records, then rebuild Amber/GROMACS:

```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/repair_mcpb_pdb_residue_coords.py \
  --target-pdb TJ_M7_FE_mcpbpy.pdb \
  --reference-pdb initial_radical_complex.pdb \
  --out-pdb TJ_M7_FE_mcpbpy_coordfix.pdb \
  --replace-residue A:187 \
  --clamp-ca-cb
```

Treat MCPB force constants and full-protein coordinates as separate products:
the former can remain useful even when a careless coordinate patch must be
repaired before MD.

Docketed poses can have perfectly valid heavy-atom geometry but impossible
hydrogen contacts after Amber/OpenBabel adds hydrogens. If minimization starts
with astronomical `LJ (SR)` such as `1e16`, inspect the closest ligand-protein
contacts in `complex_GMX.gro`. In the TJ 1T47 handoff, the selected pose had a
ligand hydrogen only 0.010 nm from a Val side-chain hydrogen. Moving only the
offending ligand hydrogen while preserving its parent X-H bond length allowed
steepest descent to converge without moving the docked heavy-atom pose.

The generated GROMACS script runs energy minimization without the user-supplied
GPU PME flags because `-pme gpu` is not valid for the steepest-descent minimizer.
NVT/NPT/production then use the repeated `--mdrun-arg` tokens. The preparation
script also writes and includes `posre_complex.itp` so `-DPOSRES` is a real
protein/ligand position restraint, not an unused macro.

Prefer flat project directories for real use. A project root such as `/mnt/e/TJ`
should expose the important outputs directly:

```text
/mnt/e/TJ/
  dock_history/
  pose_analysis/
  pose_report/
  md_handoff/
    input/
    amber/
    gromacs/
  md_analysis/
    rdc/
  pymol/
```

Avoid burying production artifacts under nested names such as
`ple_project_product_R/md_handoff/S000001_mode6_hppFe_harmonic50ns_v2` unless
running multiple experimental handoffs side by side. When a specific selected
pose is known, pass `--outdir /mnt/e/TJ/md_handoff` to `ple md-handoff` and
write analysis to `/mnt/e/TJ/md_analysis/<analysis_name>`.

After production MD, fix PBC before PyMOL/VMD viewing:

```bash
ple md-fix-pbc \
  --gmxdir /mnt/e/TJ/md_handoff/gromacs \
  --gmx /home/qin/softwares/gromacs-2026.0-gpu/bin/gmx_mpi \
  --gmxrc /home/qin/softwares/gromacs-2026.0-gpu/bin/GMXRC \
  --extra-ld-library-path /home/qin/softwares/cp2k-2026.1/build_minlib/src
```

This writes `md_nojump.xtc`, `md_centered_compact.xtc`, and
`md_centered_compact.gro`. The default groups are `Protein` for centering and
`System` for output/nojump, avoiding fragile numeric index choices. Use the
centered/compact pair for visualization; keep the raw `md.xtc` for provenance.
For MCPB/ACPYPE direct GROMACS runners generated by
`scripts/write_gromacs_mcpb_run.py`, the same post-production convention is:
`production_*_nojump.xtc`,
`production_*_centered_protein_compact.xtc`, and
`production_*_centered_protein_compact.gro`.

Build an RDC dashboard from a production trajectory:

```bash
ple md-rdc \
  --gmxdir /mnt/e/TJ/md_handoff/gromacs \
  --outdir /mnt/e/TJ/md_analysis/rdc \
  --topology /mnt/e/TJ/md_handoff/gromacs/md.gro \
  --trajectory /mnt/e/TJ/md_handoff/gromacs/md.xtc \
  --gro /mnt/e/TJ/md_handoff/gromacs/md.gro \
  --reference-pdb /mnt/e/TJ/md_handoff/input/receptor_md.pdb
```

For PyMOL metal-site viewing:

```bash
ple pymol-metal \
  --structure /mnt/e/TJ/md_handoff/gromacs/md_centered_compact.gro \
  --trajectory /mnt/e/TJ/md_handoff/gromacs/md_centered_compact.xtc \
  --out /mnt/e/TJ/pymol/show_metal_centered.pml
```

The generated PML defines `show_metal`. It shows protein cartoon, ligand sticks,
Fe as a small sphere (`sphere_scale 0.1`), side chains within 4 A as sticks,
carbon in PyMOL green, and keeps N/O/S heteroatom colors close to PyMOL defaults
with black distance dashes and larger black labels.
It intentionally avoids recoloring the whole protein so a user's `pymolrc` style
is less likely to be overwritten.

The handoff intentionally stops at real tool errors. The common first failures
are missing AmberTools/ACPYPE/ParmEd, ligand charge mistakes, nonstandard
residue names, metals, heme, calcium/iron coordination, or missing custom
parameters. Fix the first failing log under `logs/` before running production.

For H++/protonation preprocessing, sequence gaps in crystal structures are
common and often occur because a loop has no clear electron density. A docking
engine may tolerate the missing loop, but H++ can stop because pKa estimates near
chain breaks become unreliable. Do not treat this as a docking bug. First decide
whether the gap is near the Fe center, ligand pose, or catalytic residues:

- If the gap is far from the pocket, make an H++-only cleaned copy or manually
  curate only pocket/catalytic protonation states, then pass those names with
  `--protonate` or `--protonation-map`.
- If the gap is near Fe, the ligand, or catalytic acid/base residues, repair or
  model the missing residues before serious MD.
- For production MD, avoid leaving broken covalent chain geometry inside the
  simulated protein unless the missing segment is intentionally truncated with
  proper terminal treatment.

PLE can prepare isolated one-residue internal gap repairs for Rosetta:

```bash
ple repair-loops --pdb target.pdb --chain A --resseq 108 --run-rosetta
```

The command compares `SEQRES` with observed `ATOM` residues, writes
`missing_residues.csv`, and automatically targets only isolated single-residue
internal gaps. For each target it builds a chain-only PDB with the missing
residue inserted from the SEQRES identity using Biopython internal coordinates,
fits that tripeptide template onto the observed left/right flanks, writes a
Rosetta `loopmodel` KIC/refine script, and writes `audit_loopmodel.py`.

For `pdb2r5v.pdb`, H++ reports a chain A break between residues 107 and 109. The
SEQRES-derived repair target is `A:108 GLU`. Direct Rosetta Remodel
length-change blueprints were tested but were less reliable on this structure:
they either failed closure or hit a Remodel internal crash. The robust local
route was to insert GLU108 from a fitted `SER-GLU-ALA` template, then run
Rosetta `loopmodel` over residues 104-112. A successful generated run produced
347 chain-A residues with local sequence `...GQHSEAAVTT...`, peptide C-N
distances near 1.2-1.34 A around GLU108, and Rosetta `chainbreak 0.0178211`.

Do not expect this command to auto-fix long missing loops. In `pdb2r5v.pdb`,
chain B contains multiple internal gaps; those are reported but not automatically
modelled by this one-residue repair mode. For H++ on this structure, use the
repaired chain A copy or deliberately decide whether chain B should be removed,
modelled with heavier sampling, or treated as a separate curated system.

Initialize a project skeleton:

```bash
ple init --project-dir /path/to/project
```

Prepare only the receptor:

```bash
ple prepare-receptor
```

If the project directory has exactly one receptor-like file, it is used
automatically. If it has multiple `*.pdb`/`*.pdbqt` files, require `--pdb`.

## Expected Project Layout

```text
my_project/
  target.pdb or target.pdbqt
  ligands.sdf                 optional
  ple_project/
    configs/
      initial_large_box.txt
      focused_box_round1.txt
    dock_history/
    pose_analysis/
    pose_report/index.html
    score_space_model/          optional with --build-space-dashboard
    score_space_dashboard/index.html
    md_handoff/                 optional from ple md-handoff
```

The project state is written under `ple_project/` so the receptor/ligand inputs
stay easy to see.

## Practical Defaults

- Initial whole-protein box padding: `8 A`
- Initial max box side: `80 A`
- Focused box target side: `21 A`
- Focused box max side: `21 A`
- Focused box padding around the selected pose cluster: `6 A`
- Focused box cluster radius: `18 A`
- Exploration degree: `1`
- First analog batch: `24`
- Focused analog batch: `48`
- Exhaustiveness: `8`
- Num modes: `9`
- Docking engine command: `auto`, override with `--engine-exe`
- Receptor PDB preparation: `--receptor-prep auto`, using OpenBabel rigid PDBQT
  fallback (`obabel receptor.pdb -O receptor.pdbqt -xr -xc -xn`)

Use `--dry-run` to print planned commands without docking.
`--dry-run` also skips the docking-engine existence check.

## Input Guidance

Preferred input is manual SMILES:

```bash
ple run --smiles "CCO" --nickname ethanol
```

For files:

- `--smiles-file ligands.smi`: one SMILES per line, optional second token as
  nickname.
- `--csv ligands.csv --smiles-column smiles`: CSV input.
- `--sdf ligands.sdf` or `--mol ligand.mol`: use RDKit to extract SMILES.

PDF is not a reliable ligand input format. If a ligand is only in a PDF, extract
or redraw it into SMILES/SDF first; do not ask this skill to infer chemistry from
PDF images.

## Receptor Handling

Auto-detection:

- one `*.pdbqt`: use it directly for Vina.
- one `*.pdb`: prepare `ple_project/receptor/<stem>.pdbqt` with OpenBabel by
  default, then use that PDBQT for Vina.
- multiple receptor files: fail and ask for `--pdb`.

This receptor preparation is a practical fallback, not expert curation.
Protonation, charges, waters, cofactors, metals, and chains are project-specific.
For routine or publication-quality docking, prefer storing a curated
`target.pdbqt` in the project directory or passing it with `--pdb`.

## Box Logic

When no `--config` is supplied, compute the receptor atom bounding box from
PDB/PDBQT coordinates:

- center = receptor coordinate midpoint
- size = coordinate span + `2 * --initial-padding`
- clamp each side to `--initial-max-size`

After initial docking, compute a focused box from ligand pose atom coordinates in
successful pose PDBQT files:

- center = pose atom bounding midpoint
- size = pose atom span + `2 * --focus-padding`
- enforce `--focus-min-size`

The focused box is evidence-driven. If poses scatter across the protein, keep the
large box and warn that the target pocket is not yet stable.

## Related Skills

Use scripts from:

- `smiles-to-vina-docking/scripts/dock_smiles.py`
- `pose-analyzer/scripts/infer_binding_mode_families.py`
- `pose-analyzer/scripts/build_pose_report.py`

The package locates these sibling skills automatically when installed inside the
same `skills/` directory. If needed, pass `--skills-dir /path/to/skills`.

## Metal / MD References

For Fe-containing catalytic intermediates, H++ handoff, radical ligands,
GAMESS-backed MCPB.py, and TJ-style directory cleanup, read:

- `references/FE_RADICAL_MCPB_ZH.md`
- `references/MCPB_BUG_AVOIDANCE_ZH.md`
- `references/TJ_CANONICAL_LAYOUT_ZH.md`

These notes are the durable memory for repeated pitfalls: chain A/B confusion,
H++ residue names, radical H deletion, electron-count/multiplicity parity,
MCPB-generated small/large model charge changes, empty `*_small_fc.inp` files,
PBC-fixed trajectory requirements before PyMOL/RDC analysis, and keeping TJ-style
project directories small enough for other people to reproduce.

For MCPB small models, do not use a blunt "missing H" rule. After `MCPB.py -s 1`,
audit only positions that should chemically have hydrogens in the capped small
model, such as His ring/cap carbons and Glu/Asp CB/CG/cap carbons; do not flag
Fe-donor nitrogens or deprotonated carboxylate oxygens that should be
unprotonated. Use `scripts/audit_mcpb_small_model_hydrogens.py`; if
`warning_count > 2`, hard-stop before GAMESS OPT/Hessian.

Before `MCPB.py -s 1`, repair hydrogens only on standard protein residues with
`scripts/repair_mcpb_protein_hydrogens.py`. The repair scope is deliberately
limited to protein `ATOM` records such as Fe-bound His and nearby Glu/Asp
aliphatic hydrogens. Never auto-add or auto-delete hydrogens on nonstandard
ligands/radicals such as `UNL`, `UNK`, `ACT`, nitrene-like fragments, or
substrates; those are mol2/manual-curation territory. After step 1, PLE should
write a human visual-check file with
`scripts/build_mcpb_visual_check_pdb.py`, e.g. `<group>_small_visual_check.pdb`,
before enforcing the hard audit. This gives a PyMOL-inspectable active-site QM
model even when the later electron/multiplicity or hydrogen audit blocks the
run.

Hydrogen count affects spin multiplicity indirectly through the electron count:
`total_electrons = sum(atomic_numbers) - model_charge`. Even electron counts
require odd multiplicities such as 1/3/5/7; odd electron counts require even
multiplicities such as 2/4/6. Therefore adding or removing one H can make a
previous M5/M7 setting invalid.

Formal charge/multiplicity parity is necessary but not sufficient. A fresh
Fe-radical GAMESS SCF branch must start with the slow stable mode
`prepare_gamess_fixed_energy.py --scf-mode damp` unless a matching clean
converged `$VEC` already exists. It must also show a sane early trend. If the
first 10-25 iterations show an electronic-state cliff, for example total energy
becoming tens to hundreds of Hartree less stable while density/orbital-gradient
remain large, stop that spin branch and never pass its `$VEC` into OPT/Hessian.
For long fresh-SCF attempts, run the generated `run_<job>_guarded.sh` so PLE's
`guard_gamess_scf.py` can automatically stop early electronic cliffs.
The guard should be patient for fresh HUCKEL starts: do not kill merely because
the first 5-10 rows swing; wait for at least about 25 rows unless the log shows
catastrophic divergence.
