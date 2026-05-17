---
name: pose-analyzer
description: |
  Downstream structural analysis for docked Vina poses. Analyze pose mode families,
  tune the Vina grid box iteratively, visualize pocket contacts in PyMOL, and export
  pocket atom/residue datasets for geometric analysis or later ML. It also builds a
  Random-Forest-based binding-score exploration model from scored history rows.
  Designed to work with smiles-to-vina-docking and druglike-pocket-refiner as the
  final structural triage step.
---

# Pose Analyzer

Use this skill **after** docking and optional drug-like refinement, when you need to:

- Understand which Vina pose modes are structurally consistent across candidates.
- Iteratively optimize the Vina grid box based on structural coverage signals.
- Visualize ligand-receptor contacts and predicted binding surfaces in PyMOL.
- Explore the scored chemical space with a lightweight supervised model.
- Export clean `pocket_graph.json` records for later geometric or GNN analysis.

This skill does **not** predict binding scores de novo. It consumes scores or rankings
that already exist in your workflow (official assay scores, docking contact
efficiencies, or druglike-refinement scores) and uses them to guide structural
interpretation, box tuning, and exploratory modeling.

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
   (properties gated)        (QED + synth aware)        visualize, model, export
```

You can also run this skill **before** refinement to obtain a tuned box, then feed
that box back into a new docking round.

## Core Idea

1. Read the scored-ligand table (`binding_score_analysis.csv` or equivalent).
2. Pick the current best-scored ligand as the structural reference standard.
3. Parse all Vina modes from that ligand's pose file.
4. For every ligand in the scored set, parse every mode and compute
   receptor-contact fingerprints.
5. Compare all scored modes to each reference-standard mode to find structurally
   similar families.
6. Identify which reference mode family is most common among high-scoring ligands.
7. Export pocket atom/residue datasets around each selected common mode.
8. Optionally tune the Vina grid box to maximize correlation between score and
   predicted binding-surface coverage.
9. Optionally train a lightweight Random Forest regressor on scored history rows
   to explore which structural features correlate with the teacher signal.

## Main Script: Structural Analysis

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

Every continued shell line must end with `\`. If the line before `--tune-grid-box`
does not end with `\`, bash will run `--tune-grid-box` as a separate command and
print `command not found`. `--config` is required for actual grid-box tuning;
without it the script can only analyze the full receptor surface and prints a
warning.

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

Add `--save-pse` to also execute the generated PyMOL scripts and save `.pse`
sessions. Without `--save-pse`, the skill writes portable `.pml` files only.

## Score Space Model (Exploratory Supervised Analysis)

`score_space_model.py` trains a **Random Forest Regressor** on scored historical
molecules to explore which features correlate with the teacher signal
(`official_binding_score`). It is **diagnostic**, not a production QSAR model.

### What model is used?

- **Regressor**: scikit-learn `RandomForestRegressor` (`n_estimators=400` trees,
  `min_samples_leaf=5`, `random_state` fixed for reproducibility).
- **Preprocessing**: median imputation for missing values; no standardization
  because tree models are scale-invariant.
- **Validation**: 5-fold cross-validation on the training set (reported as
  `cv_r2_mean±std` and `cv_mae_mean±std`), plus a held-out test set.
- **Features**: a curated set of ~20 structural and physicochemical descriptors
  (see below). **PCA coordinates are used only for visualization**, not as model
  inputs.

### Two visualization spaces (PCA, not model features)

Both spaces are 2-D projections for human inspection; the Random Forest consumes
raw features, not these compressed coordinates.

| Space | Input to PCA | What proximity means |
|:---|:---|:---|
| **Pure-SMILES-space** | 2048-bit Morgan fingerprint (radius 2) | Similar chemical substructure / scaffold |
| **structural-interaction-space** | Docking contacts, physicochemical descriptors, ChEMBL scaffold similarity | Similar binding behavior and drug-like profile |

### Feature list

The model draws from three sources:

1. **Docking pose quality**: `affinity_kcal_mol`, `inner_rmsd`, `whole_rmsd`
2. **Pocket contacts**: `hbond_count`, `hydrophobic_count`, `pi_contact_count`,
   `ch_pi_contact_count`, plus 3-D pose-context features when available
   (`contact_residue_count_4a`, `surface_contact_fraction_4a`,
   `receptor_residue_coverage_4a`, `ligand_atom_contact_fraction_4a`,
   `min_ligand_receptor_distance`)
3. **Molecular properties & ChEMBL**: `logp`, `hbd`, `hba`, `tpsa`, `rot_bonds`,
   `qed`, `aromatic_rings`, `chembl_scaffold_similarity`

Redundant features (`heavy_atoms`, `mode_count`, `formal_charge`,
`druglike_refinement_score`, `vdw_contact_count`, etc.) are excluded to reduce
multicollinearity.

### Running the model

```bash
python /path/to/pose-analyzer/scripts/score_space_model.py \
  --history-csv ~/vina_task2/dock_history/dock_history.csv \
  --analysis-csv ~/vina_task2/dock_history/binding_score_analysis.csv \
  --reference-smiles-csv ~/vina_task2/references/chembl_kinase_smiles.csv \
  --outdir ~/vina_task2/score_space_model \
  --train-size 60 \
  --test-size 20
```

Outputs:

- `binding_score_feature_matrix.csv`: one row per historical molecule with teacher
  score, predicted score, Pure-SMILES-space coordinates,
  structural-interaction-space coordinates, descriptors, docking contacts, and
  ChEMBL scaffold similarity.
- `activity_cliff_report.csv`: pairs of molecules that are chemically or
  space-near but have large official binding-score gaps. The report includes
  Morgan Tanimoto, MCS fraction, Pure-SMILES/structural-space distances, Vina
  direction, and rule-based penalty hypotheses such as possible Vina false
  positive, lost H-bond network, less consistent pose family, overpacking,
  lipophilicity/desolvation shift, or conformational entropy penalty.
- `binding_score_model_metrics.json`: model metrics plus `activity_cliff_rows`
  and a reminder that activity-cliff explanations are hypotheses for pose/SAR
  inspection, not proof of physical mechanism.
- `high_low_group_feature_comparison.csv`: descriptive comparison of high-score
  (`official_binding_score > 0.4` by default) versus low-score
  (`0 <= official_binding_score <= 0.25` by default) molecules across docking,
  contact, descriptor, ChEMBL, and embedding-space features.
- `high_low_group_members.csv`: the molecules assigned to the high/low groups,
  with SMILES and space coordinates for manual inspection.
- `high_low_group_summary.json`: group counts, thresholds, embedding-center
  distance, and per-group summary statistics.

Activity-cliff defaults are intentionally conservative enough to catch useful
SAR contradictions:

```bash
--activity-cliff-min-score-delta 0.15 \
--activity-cliff-min-tanimoto 0.55 \
--activity-cliff-min-mcs-fraction 0.80 \
--activity-cliff-max-pure-distance 0.45 \
--activity-cliff-max-structural-distance 0.50
```

Add `--activity-cliff-require-vina-signal` when you only want cliffs where Vina
is similar or points in the wrong direction relative to the official score.

High/low group thresholds can be adjusted:

```bash
--high-score-threshold 0.4 --low-score-max 0.25
```

This is the first lightweight step toward S1-style pocket embedding analysis:
it compares existing pose/contact/descriptors and embedding coordinates before
training a GNN. Later extensions should add explicit pocket geometry vectors,
residue physicochemical composition, water/nearby atom counts, and conservation
or entropy features from pocket-graph JSON.
- `binding_score_model_metrics.json`: train/test counts, feature list, test MAE/R2,
  **5-fold CV R2/MAE**, adaptive score bins, feature importance, and a
  plain-language `model_interpretation`.
- `pure_smiles_space.html`: chemical-space map colored by adaptive binding-score bins.
- `structural_interaction_space.html`: docking/contact/descriptor-space map colored
  the same way.

### Interpreting model quality

| Signal | Interpretation |
|:---|:---|
| `cv_r2_mean < 0.3` and `test_r2 < 0.3` | Model is unreliable. Use the dashboard to inspect why current features fail to explain the score (wrong pose? noisy assay? missing pocket context?). |
| `cv_r2_mean 0.3–0.6` or `test_r2 0.3–0.6` | Weak but usable hypothesis generator. Look at feature importance and spatial clusters, but do not bet on predicted scores alone. |
| `cv_r2_mean > 0.6` and `test_r2 > 0.6` | Usable screening signal, but still needs external validation on a fresh batch. |

**Important:** predicted scores should never override experimental data or
structural intuition. They are a navigation aid, not a oracle.

### Dashboard

Build an interactive localhost dashboard from the feature matrix:

```bash
python /path/to/pose-analyzer/scripts/build_score_space_dashboard.py \
  --feature-matrix ~/vina_task2/score_space_model/binding_score_feature_matrix.csv \
  --outdir ~/vina_task2/score_space_dashboard
cd ~/vina_task2/score_space_dashboard
python -m http.server 8765 --bind 127.0.0.1
```

Open `http://127.0.0.1:8765`. If the server prints `Serving HTTP on 0.0.0.0`,
do **not** open `http://0.0.0.0:8765` in the browser. `0.0.0.0` is a server-side
listening address meaning "accept connections on all local interfaces"; it is not
a real destination address. Use `127.0.0.1`, `localhost`, or the WSL IP from
`hostname -I`.

The dashboard can switch between `Pure-SMILES-space` and
`structural-interaction-space`, recolor points by official score, ChEMBL scaffold
similarity, QED, affinity, MW, or experimental predicted score, and inspect each
molecule without treating the weak prediction model as authoritative.

### Dashboard color logic

- **Gray** points: the selected color field is missing.
- **Large outlined** points: have an official teacher binding score.
- **Small** points: unscored history rows.
- **`official binding score`**: adaptive quantile-based bins (20/40/60/80th
  percentile of scored data). Colors: red → amber → muted green → blue → green.
- **`predicted score`**: data-driven bins from the predicted distribution itself,
  independent of official score bins.
- **`Vina affinity`**: Vina-style directionality (`<=-12` dark green, `-12 to -10`
  green, `-10 to -8` amber, `>-8` orange).
- **`MW`**: diverging window scale. Red for `<250` or `>600`, green for the
  optimal oral-drug window `300–500`, yellow for transition zones.
- **`ChEMBL scaffold similarity`** and **`QED`**: continuous low-to-high scale.

## Grid-Box Tuning

`--config` is strongly recommended. The script reads Vina `center_*` and `size_*`
values, then defines a predicted binding surface from residues inside that grid
plus `--grid-surface-padding` A, default 4 A. This deliberately includes likely
pocket-forming residues near the box edge, not only the residues touched by one
ligand pose. Use `--surface-only-grid-residues` only when you want a stricter
SASA-exposed subset.

Add `--tune-grid-box` when the current Vina box is only a rough guess. The script
performs a local search without editing the original config file:

- center moves by `--grid-center-step` A, default 0.5 A
- size changes by `--grid-size-step` A, default 1.0 A
- each size dimension is clamped by `--grid-min-size`, default 18 A
- the chosen box maximizes Pearson correlation between the supplied score column
  and predicted binding-surface coverage
- search uses deterministic multi-start itineraries via `--grid-tune-restarts`
- `--grid-target-pearson` defaults to 0.7; `--grid-tune-max-steps` defaults to 1000
- `--grid-tune-metric atom` optimizes coverage by non-hydrogen receptor atoms;
  `--grid-tune-metric residue` optimizes by residues. PML visualization still shows
  residues.

`grid_box_tuning.csv` reports `pearson_r`, `objective_score`, and `loss`. The
terminal prints `grid_tuning_result=...` with Pearson, Spearman, loss, steps,
center, and size.

Interpretation warning: if 1000-step tuning cannot approach the target Pearson, do
not overfit the box. Treat that as evidence that one box-coverage feature is not
enough; combine it with interaction features (H-bond-like contacts, hydrophobic
contacts, aromatic/CH-pi contacts, ligand burial, mode-family support) in
downstream analysis.

The tuning audit is written to `grid_box_tuning.csv`, and the selected box is
written to `selected_grid_box.csv`.

## Practical Interpretation

Treat this as **structural triage**, not final mechanistic proof.

Useful signals:

- High support count for one reference mode family.
- Repeated contact residues among high-scoring ligands.
- Pocket residue coverage: how much of the receptor residue set is touched within
  4–5 A.
- Predicted binding-surface coverage: how much of the Vina-grid binding surface is
  touched.
- Surface overlap: whether the ligand contacts surface-like residues or is buried.
- More ligand atoms within 4 A of receptor atoms.
- More H-bond-like, hydrophobic, aromatic, and CH-pi-like contacts.
- Low minimum ligand-receptor distance without obvious steric crash.

If an external score (e.g. official assay score) and Vina affinity disagree,
prefer the structural signals from this skill over raw Vina affinity.

`receptor_surface_residues.csv` uses PyMOL solvent-accessible surface area when
PyMOL is available. The fallback is residue-neighbor density and is marked as
`surface_method=neighbor_density`. Prefer PyMOL SASA or FreeSASA/MSMS when surface
exposure matters.

The PyMOL views use `predicted_binding_surface` in lime and
`predicted_surface_contact_4a` / contact residues in `tv_orange`. Ranked PML
filenames include both binding-surface coverage and surface-contact fraction so
the source pose remains traceable.

## Future ML Direction

The first reliable ML target is not raw ligand SMILES. Use exported
`pocket_graph.json` records:

- ligand atoms as one node type,
- receptor pocket atoms as another,
- coordinates preserved,
- edges by distance bins or RBF distance encoding,
- labels from whatever score column you supplied,
- optional pose-family labels from `mode_family_summary.csv`.

The user's older `GNN_eigen_extract_failed.py` is useful conceptually for geometry
embeddings and EGNN-style models, but the first version of this skill should
generate clean pocket data before training.
