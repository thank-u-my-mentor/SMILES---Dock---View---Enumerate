---
name: pose-analyzer
description: Downstream structural analysis for docked Vina poses. Analyze pose mode families,
  tune the Vina grid box iteratively, visualize pocket contacts in PyMOL, and export pocket
  atom/residue datasets for geometric analysis or later ML. Designed to work with
  smiles-to-vina-docking and druglike-pocket-refiner as the final structural triage step.
---

# Pose Analyzer

Use this skill **after** docking and optional drug-like refinement, when you need to:

- Understand which Vina pose modes are structurally consistent across candidates.
- Iteratively optimize the Vina grid box based on structural coverage signals.
- Visualize ligand-receptor contacts and predicted binding surfaces in PyMOL.
- Export clean pocket_graph.json records for later geometric or GNN analysis.

This skill does **not** predict binding scores de novo. It consumes scores or rankings that
already exist in your workflow (official assay scores, docking contact efficiencies, or
druglike-refinement scores) and uses them to guide structural interpretation and box tuning.

## Prerequisites

A docking history from `smiles-to-vina-docking`, usually containing:

- `dock_history.csv`
- `binding_score_analysis.csv` (or any scored-ligand table with canonical SMILES)
- ligand pose PDBQT files under `poses/`
- receptor PDBQT or PDB file
- Vina config file with `center_*` and `size_*`

## Workflow Position

```
smiles-to-vina-docking  ->  druglike-pocket-refiner  ->  pose-analyzer
        |                          |                        |
   generate & dock          trim / cyclize /          analyze poses,
   keep history              bioisosteric edits         optimize box,
   (properties gated)        (QED + synth aware)        visualize, export
```

You can also run this skill **before** refinement to obtain a tuned box, then feed that box
back into a new docking round.

## Core Idea

1. Read the scored-ligand table (`binding_score_analysis.csv` or equivalent).
2. Pick the current best-scored ligand as the structural reference standard.
3. Parse all Vina modes from that ligand's pose file.
4. For every ligand in the scored set, parse every mode and compute receptor-contact fingerprints.
5. Compare all scored modes to each reference-standard mode to find structurally similar families.
6. Identify which reference mode family is most common among high-scoring ligands.
7. Export pocket atom/residue datasets around each selected common mode.
8. Optionally tune the Vina grid box to maximize correlation between score and predicted
   binding-surface coverage.

## Main Script

Basic structural analysis:

```bash
python /path/to/pose-analyzer/scripts/infer_binding_mode_families.py \
  --analysis-csv ~/vina_task2/dock_history/binding_score_analysis.csv \
  --history-csv ~/vina_task2/dock_history/dock_history.csv \
  --receptor ~/vina_task2/vina_bin/target.pdbqt \
  --config ~/vina_task2/vina_bin/idz6F_config.txt \
  --outdir ~/vina_task2/pose_analysis
```

With iterative grid-box tuning:

```bash
python /path/to/pose-analyzer/scripts/infer_binding_mode_families.py \
  --analysis-csv ~/vina_task2/dock_history/binding_score_analysis.csv \
  --history-csv ~/vina_task2/dock_history/dock_history.csv \
  --receptor ~/vina_task2/vina_bin/target.pdbqt \
  --config ~/vina_task2/vina_bin/idz6F_config.txt \
  --outdir ~/vina_task2/pose_analysis \
  --tune-grid-box \
  --grid-target-pearson 0.7 \
  --grid-tune-max-steps 10000 \
  --grid-tune-metric residue
```

Every continued shell line must end with `\`. If the line before `--tune-grid-box` does not
end with `\`, bash will run `--tune-grid-box` as a separate command and print `command not found`.
`--config` is required for actual grid-box tuning; without it the script can only analyze the
full receptor surface and prints a warning.

Outputs:

- `mode_family_summary.csv`: reference mode families and their support among scored ligands.
- `scored_mode_matches.csv`: every scored-ligand mode matched to the closest reference mode.
- `consensus_residue_contacts.csv`: residues repeatedly touched by high-score modes.
- `pocket_residues.csv`: residue-level pocket membership for each reference mode.
- `receptor_surface_residues.csv`: approximate surface-like residue table for the full receptor.
- `predicted_binding_surface_residues.csv`: residues in the Vina grid box plus a small pocket-shell padding.
- `predicted_binding_surface.pdb`: receptor atoms for the predicted binding surface, loaded into PyMOL as `predicted_binding_surface`.
- `pockets/standard_mode_*/pocket_atoms.pdb`: receptor atoms near the standard mode.
- `pockets/standard_mode_*/ligand_mode.pdb`: ligand atoms for that mode.
- `pockets/standard_mode_*/pocket_graph.json`: ligand atoms, receptor atoms, distances, contacts, and features for later ML.
- `pockets/standard_mode_*/view.pml`: PyMOL view for that mode.
- `visualizations/rank*_surface*_seq*_mode*.pml`: ranked PyMOL views ordered by surface-contact fraction, then support.

Add `--save-pse` to also execute the generated PyMOL scripts and save `.pse` sessions.
Without `--save-pse`, the skill writes portable `.pml` files only.

## Grid-Box Tuning

`--config` is strongly recommended. The script reads Vina `center_*` and `size_*` values, then
defines a predicted binding surface from residues inside that grid plus `--grid-surface-padding` A,
default 4 A. This deliberately includes likely pocket-forming residues near the box edge, not
only the residues touched by one ligand pose. Use `--surface-only-grid-residues` only when you
want a stricter SASA-exposed subset.

Add `--tune-grid-box` when the current Vina box is only a rough guess. The script performs a
local search without editing the original config file:

- center moves by `--grid-center-step` A, default 0.5 A
- size changes by `--grid-size-step` A, default 1.0 A
- each size dimension is clamped by `--grid-min-size`, default 18 A
- the chosen box maximizes Pearson correlation between the supplied score column and predicted
  binding-surface coverage
- search uses deterministic multi-start itineraries via `--grid-tune-restarts`
- `--grid-target-pearson` defaults to 0.7; `--grid-tune-max-steps` defaults to 1000
- `--grid-tune-metric atom` optimizes coverage by non-hydrogen receptor atoms;
  `--grid-tune-metric residue` optimizes by residues. PML visualization still shows residues.

`grid_box_tuning.csv` reports `pearson_r`, `objective_score`, and `loss`. The terminal prints
`grid_tuning_result=...` with Pearson, Spearman, loss, steps, center, and size.

Interpretation warning: if 1000-step tuning cannot approach the target Pearson, do not overfit
the box. Treat that as evidence that one box-coverage feature is not enough; combine it with
interaction features (H-bond-like contacts, hydrophobic contacts, aromatic/CH-pi contacts,
ligand burial, mode-family support) in downstream analysis.

The tuning audit is written to `grid_box_tuning.csv`, and the selected box is written to
`selected_grid_box.csv`.

## Practical Interpretation

Treat this as **structural triage**, not final mechanistic proof.

Useful signals:

- High support count for one reference mode family.
- Repeated contact residues among high-scoring ligands.
- Pocket residue coverage: how much of the receptor residue set is touched within 4-5 A.
- Predicted binding-surface coverage: how much of the Vina-grid binding surface is touched.
- Surface overlap: whether the ligand contacts surface-like residues or is buried.
- More ligand atoms within 4 A of receptor atoms.
- More H-bond-like, hydrophobic, aromatic, and CH-pi-like contacts.
- Low minimum ligand-receptor distance without obvious steric crash.

If an external score (e.g. official assay score) and Vina affinity disagree, prefer the
structural signals from this skill over raw Vina affinity.

`receptor_surface_residues.csv` uses PyMOL solvent-accessible surface area when PyMOL is available.
The fallback is residue-neighbor density and is marked as `surface_method=neighbor_density`.
Prefer PyMOL SASA or FreeSASA/MSMS when surface exposure matters.

The PyMOL views use `predicted_binding_surface` in lime and `predicted_surface_contact_4a` /
contact residues in `tv_orange`. Ranked PML filenames include both binding-surface coverage and
surface-contact fraction so the source pose remains traceable.

## Future ML Direction

The first reliable ML target is not raw ligand SMILES. Use exported `pocket_graph.json` records:

- ligand atoms as one node type,
- receptor pocket atoms as another,
- coordinates preserved,
- edges by distance bins or RBF distance encoding,
- labels from whatever score column you supplied,
- optional pose-family labels from `mode_family_summary.csv`.

The user's older `GNN_eigen_extract_failed.py` is useful conceptually for geometry embeddings
and EGNN-style models, but the first version of this skill should generate clean pocket data
before training.
