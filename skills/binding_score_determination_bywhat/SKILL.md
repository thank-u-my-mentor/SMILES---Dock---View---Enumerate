---
name: binding_score_determination_bywhat
description: Analyze official binding scores against Vina docked poses to infer which ligand modes, receptor contacts, and pocket environments are favored; export pocket atom datasets for later geometric or SE(3)-equivariant GNN analysis.
---

# Binding Score Determination By What

Use this skill when the user wants to understand why official Task2 binding scores prefer some docked ligands or poses over others.

This skill assumes a docking history already exists, usually from `smiles-to-vina-docking`, with:

- `dock_history.csv`
- `binding_score_analysis.csv`
- ligand pose PDBQT files under `poses/`
- receptor PDBQT or PDB file

## Core Idea
1. Read `binding_score_analysis.csv`.
2. Pick the row with the best `official_binding_score` as the current best standard.
3. Parse all modes from that ligand's `pose_path`.
4. For every scored ligand, parse every Vina mode and compute receptor-contact fingerprints.
5. Compare all scored modes to each best-standard mode.
6. Find which best-standard mode is most common among high-scoring ligands.
7. Export pocket atom/residue datasets around each selected common mode.

The output is intended to answer:
- Which Vina mode family appears most compatible with the official binding score?
- Which residues repeatedly contact high-score ligands?
- Are high-score ligands deeply buried, surface-exposed, or contacting a consistent pocket shell?
- Which pocket JSON/PDB exports should a future SE(3)-equivariant GNN consume?

## Main Script

Use:
```bash
python /path/to/binding_score_determination_bywhat/scripts/infer_binding_mode_families.py \
  --analysis-csv ~/vina_task2/dock_history/binding_score_analysis.csv \
  --history-csv ~/vina_task2/dock_history/dock_history.csv \
  --receptor ~/vina_task2/vina_bin/target.pdbqt \
  --config ~/vina_task2/vina_bin/idz6F_config.txt \
  --outdir ~/vina_task2/binding_score_mode_analysis
```

For a tuned binding-surface box:
```bash
python /path/to/binding_score_determination_bywhat/scripts/infer_binding_mode_families.py \
  --analysis-csv ~/vina_task2/dock_history/binding_score_analysis.csv \
  --history-csv ~/vina_task2/dock_history/dock_history.csv \
  --receptor ~/vina_task2/vina_bin/target.pdbqt \
  --config ~/vina_task2/vina_bin/idz6F_config.txt \
  --outdir ~/vina_task2/binding_score_mode_analysis \
  --tune-grid-box \
  --grid-target-pearson 0.7 \
  --grid-tune-max-steps 10000 \
  --grid-tune-metric residue
```

Every continued shell line must end with `\`. If the line before `--tune-grid-box` does not end with `\`, bash will run `--tune-grid-box` as a separate command and print `command not found`. `--config` is required for actual grid-box tuning; without it the script can only analyze the full receptor surface and prints a warning.
Outputs:

- `mode_family_summary.csv`: best-standard mode families and their support among scored ligands.
- `scored_mode_matches.csv`: every scored ligand mode matched to the closest best-standard mode.
- `consensus_residue_contacts.csv`: residues repeatedly touched by high-score modes.
- `pocket_residues.csv`: residue-level pocket membership for each best-standard mode.
- `receptor_surface_residues.csv`: approximate surface-like residue table for the full receptor.
- `predicted_binding_surface_residues.csv`: residues in the Vina grid box plus a small pocket-shell padding; this is the main "likely binding surface" table.
- `predicted_binding_surface.pdb`: receptor atoms for the predicted binding surface, loaded into PyMOL views as `predicted_binding_surface`.
- `pockets/standard_mode_*/pocket_atoms.pdb`: receptor atoms near the standard mode.
- `pockets/standard_mode_*/ligand_mode.pdb`: ligand atoms for that mode.
- `pockets/standard_mode_*/pocket_graph.json`: ligand, receptor atoms, distances, contact residues, and features for later ML.
- `pockets/standard_mode_*/view.pml`: PyMOL view for that mode.
- `visualizations/rank*_surface*_seq*_mode*.pml`: ranked PyMOL views ordered by surface-contact fraction, then support.
Add `--save-pse` to also execute the generated PyMOL scripts and save `.pse` sessions. Without `--save-pse`, the skill writes portable `.pml` files only.

`--config` is strongly recommended. The script reads Vina `center_*` and `size_*` values, then defines a predicted binding surface from residues inside that grid plus `--grid-surface-padding` A, default 4 A. This deliberately includes likely pocket-forming residues near the box edge, not only the residues touched by one ligand pose. Use `--surface-only-grid-residues` only when you want a stricter SASA-exposed subset.

Add `--tune-grid-box` when the current Vina box is only a rough guess. The script performs a local search without editing the original config file:
- center moves by `--grid-center-step` A, default 0.5 A
- size changes by `--grid-size-step` A, default 1.0 A
- each size dimension is clamped by `--grid-min-size`, default 18 A
- the chosen box maximizes Pearson correlation between official binding score and predicted binding-surface coverage
- the search uses deterministic multi-start itineraries via `--grid-tune-restarts`, so it is not limited to one greedy path from the original box
- `--grid-target-pearson` defaults to 0.7, and `--grid-tune-max-steps` defaults to 1000 candidate boxes
- `--grid-tune-metric atom` optimizes coverage by non-hydrogen receptor atoms; `--grid-tune-metric residue` optimizes by residues. PML visualization still shows residues either way.
- `grid_box_tuning.csv` reports `pearson_r`, `objective_score`, and `loss`; currently `objective_score = pearson_r` and `loss = max(0, target_pearson - pearson_r)`
- the terminal prints `grid_tuning_result=...` with Pearson, Spearman, loss, steps, center, and size so the run quality is visible immediately.

Interpretation warning: if 1000-step tuning cannot approach the target Pearson, do not overfit the box. Treat that as evidence that one box-coverage feature is not enough and combine it with interaction features such as H-bond-like contacts, hydrophobic contacts, aromatic/CH-pi contacts, ligand burial, and mode-family support.

The tuning audit is written to `grid_box_tuning.csv`, and the selected box is written to `selected_grid_box.csv`.

## Practical Interpretation

Treat this as structural triage, not final mechanistic proof.

Useful signals:

- High support count for one standard mode family.
- Repeated contact residues among high official scores.
- Pocket residue coverage: how much of the receptor residue set is touched by the ligand within 4-5 A.
- Predicted binding-surface coverage: how much of the Vina-grid binding surface is touched within 4-5 A.
- Surface overlap: whether the ligand contacts surface-like residues or is buried in an interior pocket.
- More ligand atoms within 4 A of receptor atoms.
- More H-bond-like, hydrophobic, aromatic, and CH-pi-like contacts.
- Low minimum ligand-receptor distance without obvious steric crash.

If official binding score and Vina affinity disagree, prefer this structural analysis over raw Vina affinity.

`receptor_surface_residues.csv` uses PyMOL solvent-accessible surface area when PyMOL is available. The fallback is residue-neighbor density and is marked as `surface_method=neighbor_density`. Prefer PyMOL SASA or FreeSASA/MSMS when surface exposure matters.

The PyMOL views use `predicted_binding_surface` in lime and `predicted_surface_contact_4a` / contact residues in `tv_orange`. Ranked PML filenames include both binding-surface coverage and surface-contact fraction so the source pose remains traceable.

## Future GNN Direction

The first reliable ML target is not raw ligand SMILES. Use exported `pocket_graph.json` records:
- ligand atoms as one node type,
- receptor pocket atoms as another,
- coordinates preserved,
- edges by distance bins or RBF distance encoding,
- labels from `official_binding_score`,
- optional pose-family labels from `mode_family_summary.csv`.

The user's older `GNN_eigen_extract_failed.py` is useful conceptually for geometry embeddings and EGNN-style models, but the first version of this skill should generate clean pocket data before training.
